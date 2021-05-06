//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Matei David (matei@cs.toronto.edu)
//---------------------------------------------------------

// Reference:
// http://stackoverflow.com/questions/14086417/how-to-write-custom-input-stream-in-c

/* 
   Certain parts of the code, specifically for BGZF support, were taken out of bgzf.h and bgzf.c
   from htslib https://github.com/samtools/htslib and so must include the MIT License below
*/

/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013-2020 Genome Research Ltd

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#pragma once

#include <cassert>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <strict_fstream.hpp>
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>

#if __cplusplus == 201103L
#include <zstr_make_unique_polyfill.h>
#endif

#if (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) \
  || defined(__LITTLE_ENDIAN__) \
  || defined(HTS_x86) \
  || defined(__ARMEL__) || defined(__THUMBEL__) || defined(__AARCH64EL__) \
  || defined(_MIPSEL) || defined(__MIPSEL) || defined(__MIPSEL__)
#define ZSTR_IS_LITTLE_ENDIAN
#endif

namespace zstr
{

static const std::size_t default_buff_size = (std::size_t)1 << 20;

static const int BGZF_BLOCK_SIZE = 0xff00; // make sure compressBound(BGZF_BLOCK_SIZE) < BGZF_MAX_BLOCK_SIZE
static const int BGZF_MAX_BLOCK_SIZE = 0x10000;
static const int BGZF_BLOCK_HEADER_LENGTH = 18;
static const int BGZF_BLOCK_FOOTER_LENGTH = 8;

static const std::size_t bgzf_default_buff_size = BGZF_MAX_BLOCK_SIZE;

/* BGZF/GZIP header (specialized from RFC 1952; little endian):
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 | 31|139|  8|  4|              0|  0|255|      6| 66| 67|      2|BLK_LEN|
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
  BGZF extension:
                ^                              ^   ^   ^
                |                              |   |   |
               FLG.EXTRA                     XLEN  B   C

  BGZF format is compatible with GZIP. It limits the size of each compressed
  block to 2^16 bytes and adds and an extra "BC" field in the gzip header which
  records the size.

*/
static const uint8_t BGZF_MAGIC_HEADER[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

typedef struct
{
    uint64_t uaddr;  // offset w.r.t. uncompressed data
    uint64_t caddr;  // offset w.r.t. compressed data
}
bgzidx1_t;

class packer {
public:

    static inline void packInt16(uint8_t *buffer, uint16_t value)
    {
    #ifdef ZSTR_IS_LITTLE_ENDIAN
        *((uint16_t *)buffer) = value;
    #else
        buffer[0] = value;
        buffer[1] = value >> 8;
    #endif
    }

    static inline uint16_t unpackInt16(const uint8_t *buffer)
    {
    #ifdef ZSTR_IS_LITTLE_ENDIAN
        return *((uint16_t *)buffer);
    #else
        return buffer[0] | buffer[1] << 8;
    #endif
    }

    static inline void packInt32(uint8_t *buffer, uint32_t value)
    {
    #ifdef ZSTR_IS_LITTLE_ENDIAN
        *((uint32_t*)buffer) = value;
    #else
        buffer[0] = value;
        buffer[1] = value >> 8;
        buffer[2] = value >> 16;
        buffer[3] = value >> 24;
    #endif
    }

    static inline uint32_t unpackInt32(const uint8_t *buffer) {
    #ifdef ZSTR_IS_LITTLE_ENDIAN
        return *((uint32_t*)buffer);
    #else
        return buffer[0] | buffer[1] << 8 | buffer[2] << 16 | buffer[3] << 24;
    #endif
    }

};


/// Exception class thrown by failed zlib operations.
class Exception
    : public std::ios_base::failure
{
public:
    static std::string error_to_message(z_stream * zstrm_p, int ret, const char *more_msg = "")
    {
        std::string msg = "zstr zlib: ";
        switch (ret)
        {
        case Z_STREAM_ERROR:
            msg += "Z_STREAM_ERROR: ";
            break;
        case Z_DATA_ERROR:
            msg += "Z_DATA_ERROR: ";
            break;
        case Z_MEM_ERROR:
            msg += "Z_MEM_ERROR: ";
            break;
        case Z_VERSION_ERROR:
            msg += "Z_VERSION_ERROR: ";
            break;
        case Z_BUF_ERROR:
            msg += "Z_BUF_ERROR: ";
            break;
        default:
            std::ostringstream oss;
            oss << ret;
            msg += "[" + oss.str() + "]: ";
            break;
        }
        if (zstrm_p->msg) {
            msg += zstrm_p->msg;
        }
        msg += " ("
                "next_in: " +
                std::to_string(uintptr_t(zstrm_p->next_in)) +
                ", avail_in: " +
                std::to_string(uintptr_t(zstrm_p->avail_in)) +
                ", next_out: " +
                std::to_string(uintptr_t(zstrm_p->next_out)) +
                ", avail_out: " +
                std::to_string(uintptr_t(zstrm_p->avail_out)) +
                ")" + more_msg;
        return msg;
    }

    Exception(z_stream * zstrm_p, int ret, const char *msg = "")
        : std::ios_base::failure(error_to_message(zstrm_p, ret, msg))
    {
    }
}; // class Exception

namespace detail
{

class z_stream_wrapper
    : public z_stream
{
public:
    z_stream_wrapper(bool _is_input, int _level, int _window_bits)
        : is_input(_is_input)
    {
        this->zalloc = Z_NULL;
        this->zfree = Z_NULL;
        this->opaque = Z_NULL;
        int ret;
        if (is_input)
        {
            this->avail_in = 0;
            this->next_in = Z_NULL;
            ret = inflateInit2(this, _window_bits ? _window_bits : 15+32);
        }
        else
        {
            ret = deflateInit2(this, _level, Z_DEFLATED, _window_bits ? _window_bits : 15+16, 8, Z_DEFAULT_STRATEGY);
        }
        if (ret != Z_OK) throw Exception(this, ret);
    }
    ~z_stream_wrapper()
    {
        if (is_input)
        {
            inflateEnd(this);
        }
        else
        {
            deflateEnd(this);
        }
    }
private:
    bool is_input;
}; // class z_stream_wrapper

} // namespace detail

class bgzf_virtual_file_pointer 
{
    // a virtual file pointer represents a byte in the uncompressed stream.
    // and can be used to efficiently seek in the bgzf_ostream to that byte
    // but does not know the precise uncompressed offset from the beginning of the stream.
    // physical_pos must always be the position in the compressed stream where
    //     a new BGZF block header starts. 
    // block_pos is the position in the uncompressed block where the desired 
    //     uncompressed offset resides 
    // compressed files in excess of 256 TB are not supported by this virtual file pointer
    uint64_t physical_pos : 48;
    uint64_t block_pos    : 16;

public:
    bgzf_virtual_file_pointer() : physical_pos{}, block_pos{} {}
    bgzf_virtual_file_pointer(size_t file_offset, uint16_t block_offset)
    : physical_pos(file_offset)
    , block_pos(block_offset) {}

    bool operator==(const bgzf_virtual_file_pointer o) const {
        return physical_pos == o.physical_pos && block_pos == o.block_pos;
    }
    bool operator!=(const bgzf_virtual_file_pointer o) const {
        return ! (*this == o);
    }
    bool operator<(const bgzf_virtual_file_pointer o) const {
        if (physical_pos < o.physical_pos) return true;
        else if (physical_pos == o.physical_pos) return block_pos < o.block_pos;
        else return false;
    }
    bool operator<=(const bgzf_virtual_file_pointer o) const {
        return (*this == o) | (*this < o);
    }
    size_t get_file_offset() const { return physical_pos; }
    uint16_t get_block_offset() const { return block_pos; }
    static const bgzf_virtual_file_pointer get_invalid() {
        return bgzf_virtual_file_pointer(0xffffff, (uint16_t)0xff);
    }
};

class bgzf_index 
{
public:

    bgzf_index() {}
    bgzf_index(std::istream &is) {
        read_index(is);
    }
    bgzf_index(const bgzf_index &copy) = default;
    bgzf_index(bgzf_index &&move) = default;
    bgzf_index &operator=(const bgzf_index &copy) = default;
    bgzf_index &operator=(bgzf_index &&move) = default;

    void rebaseline(uint64_t uncompressed_offset, uint64_t compressed_offset) {
        assert(is_valid());
        for(auto &entry : raw_index) {
            entry.uaddr += uncompressed_offset;
            entry.caddr += compressed_offset;
        }
        assert(is_valid());
    }

    std::ostream& write_index(std::ostream &os) {
        for(const auto &entry : raw_index) {
            os << entry.uaddr;
            os << entry.caddr;
        }
        return os;
    }
    void read_index(std::istream &is) {
        while (!is.eof()) {
            bgzidx1_t entry = read_next_from_index(is);
            append(entry);
        }
    }
    static bgzidx1_t read_next_from_index(std::istream &is) {
        bgzidx1_t entry;
        is >> entry.uaddr;
        is >> entry.caddr;
        return entry;
    }
    void append_incremental_block(int uncompressed_len, int compressed_len) {
        assert(compressed_len >= BGZF_BLOCK_HEADER_LENGTH + BGZF_BLOCK_FOOTER_LENGTH);
        assert(uncompressed_len >= 0);
        assert(compressed_len <= BGZF_MAX_BLOCK_SIZE);
        assert(uncompressed_len <= BGZF_BLOCK_SIZE);
        bgzidx1_t last_block{};
        if (!raw_index.empty()) {
            last_block = raw_index.back();
        }
        last_block.uaddr += uncompressed_len;
        last_block.caddr += compressed_len;
        append(last_block);
    }
    void append_absolute_block(int uncompressed_len, int compressed_len) {
        bgzidx1_t blk;
        blk.uaddr = uncompressed_len;
        blk.caddr = compressed_len;
        append(blk);
    }
    void append(bgzidx1_t entry) {
        if (raw_index.empty()) {
            raw_index.push_back({0,0}); // always start with 0,0 index
        }
        raw_index.push_back(entry);
        auto & last = raw_index[raw_index.size()-2];
        assert(last.caddr <= raw_index.back().caddr);
        assert(last.uaddr <= raw_index.back().uaddr);
    }
    
    typedef std::vector<bgzidx1_t>::const_iterator idx_iter;
    bgzf_virtual_file_pointer find_uncompressed_pointer(uint64_t uncompressed_offset) const {
        if (raw_index.empty()) return bgzf_virtual_file_pointer::get_invalid();
        assert(is_valid());
        bgzidx1_t test{};
        test.uaddr = uncompressed_offset;
        idx_iter iter = raw_index.begin();
        iter = std::lower_bound(raw_index.begin(), raw_index.end(), test, cmp_uncompressed);
        // get the one just before
        if (iter == raw_index.cend()) {
            // not found 
            test = raw_index.back();
        } else if (iter->uaddr == uncompressed_offset) {
            // exact match
        } else {
            test = *(--iter);
        }
        if (test.uaddr + BGZF_BLOCK_SIZE < uncompressed_offset || test.uaddr > uncompressed_offset) {
            // not in this block
            return bgzf_virtual_file_pointer::get_invalid();
        }
        assert(test.uaddr <= uncompressed_offset);
        assert(test.uaddr + BGZF_MAX_BLOCK_SIZE >= uncompressed_offset);
        uint16_t block_offset = uncompressed_offset - test.uaddr;
        return bgzf_virtual_file_pointer(test.caddr, (uint16_t)block_offset);
    }

    static bool cmp_compressed(const bgzidx1_t a, const bgzidx1_t b) {
        return a.caddr < b.caddr;
    }

    static bool cmp_uncompressed(const bgzidx1_t a, const bgzidx1_t b) {
        return a.uaddr < b.uaddr;
    }

    // verifies the index is monotonic in both compressed and uncompressed
    bool is_valid() const {
        bool is_monotonic = true;
        idx_iter iter, last = raw_index.begin();
        iter=last;
        while (iter != raw_index.end()) {
            if (iter != last) {
                is_monotonic &= !cmp_compressed(*iter, *last);
                is_monotonic &= !cmp_uncompressed(*iter, *last);
            }
            last = iter;
            iter++;
        }
        return is_monotonic;
    }

    size_t size() const { return raw_index.size(); }
    idx_iter begin() const { return raw_index.begin(); }
    idx_iter end() const { return raw_index.end(); }

protected:
    std::vector<bgzidx1_t> raw_index;
};

class _istreambuf
    : public std::streambuf
{
public:
    _istreambuf(std::streambuf * _sbuf_p, std::size_t _buff_size, bool _auto_detect, int _window_bits, bool _is_bgzf = false)
        : sbuf_p(_sbuf_p),
          zstrm_p(nullptr),
          buff_size(_buff_size),
          auto_detect(_auto_detect),
          auto_detect_run(false),
          is_text(false),
          window_bits(_window_bits),
          is_bgzf(_is_bgzf)
    {
        assert(sbuf_p);
        in_buff = std::make_unique<char[]>(buff_size);
        out_buff = std::make_unique<char[]>(buff_size);
        reset();
    }

    _istreambuf(const _istreambuf &) = delete;
    _istreambuf & operator = (const _istreambuf &) = delete;

    void reset() {
        zstrm_p.reset();
        in_buff_start = in_buff.get();
        in_buff_end = in_buff.get();
        setg(out_buff.get(), out_buff.get(), out_buff.get());
    }

/*
protected:
    virtual pos_type seekoff(off_type off, std::ios_base::seekdir dir,
                     std::ios_base::openmode which) override
    {
        std::cerr << "seekoff(" << off << ")" << std::endl;
        if (off != 0 || dir != std::ios_base::cur) {
            return std::streambuf::seekoff(off, dir, which);
        }

        if (!zstrm_p) {
            return 0;
        }

        return zstrm_p->total_out - in_avail();
    }
    virtual pos_type seekpos(pos_type pos,
                     std::ios_base::openmode which = std::ios_base::in) override
    {
        std::cerr << "seekpos(" << pos << ")" << std::endl;
        return std::streambuf::seekpos(pos, which);
        
    }
    */
    
public:
    std::streambuf::int_type underflow() override
    {
        std::cerr << "underflow" << std::endl;
        if (this->gptr() == this->egptr())
        {
            std::cerr << "underflow2" << std::endl;
            if (zstrm_p) std::cerr << "zstrm_p:" << (void*) zstrm_p.get() << " avail_in=" << zstrm_p->avail_in << " avail_out=" << zstrm_p->avail_out << std::endl;
            // pointers for free region in output buffer
            char * out_buff_free_start = out_buff.get();
            int tries = 0;
            do
            {
                if (++tries > 1000) {
                    throw std::ios_base::failure("Failed to fill buffer after 1000 tries");
                }

                if (!zstrm_p) {
                    in_buff_start = in_buff.get();
                    in_buff_end = in_buff.get();
                }

                // read more input if none available
                if (in_buff_start == in_buff_end)
                {
                    // empty input buffer: refill from the start
                    in_buff_start = in_buff.get();
                    std::streamsize sz = sbuf_p->sgetn(in_buff.get(), buff_size);
                    in_buff_end = in_buff_start + sz;
                    std::cerr << "read " << sz << " bytes from raw stream is_bgzf=" << is_bgzf << std::endl;
                    if (in_buff_end == in_buff_start) break; // end of input
                }
                // auto detect if the stream contains text or deflate data
                if (auto_detect && ! auto_detect_run)
                {
                    auto_detect_run = true;
                    unsigned char b0 = *reinterpret_cast< unsigned char * >(in_buff_start);
                    unsigned char b1 = *reinterpret_cast< unsigned char * >(in_buff_start + 1);
                    // Ref:
                    // http://en.wikipedia.org/wiki/Gzip
                    // http://stackoverflow.com/questions/9050260/what-does-a-zlib-header-look-like
                    is_text = ! (in_buff_start + 2 <= in_buff_end
                                 && ((b0 == 0x1F && b1 == 0x8B)         // gzip header
                                     || (b0 == 0x78 && (b1 == 0x01      // zlib header
                                                        || b1 == 0x9C
                                                        || b1 == 0xDA))));
                    if (is_bgzf && memcmp(in_buff_start, BGZF_MAGIC_HEADER, BGZF_BLOCK_HEADER_LENGTH - 2) != 0) {
                        std::cerr << "WARNING: Expecting a BGZF but did not discover the BGZF magic header. Treating this as normal zip file" << std::endl;
                        is_bgzf = false;
                    }
                }
                if (is_text)
                {
                    // simply swap in_buff and out_buff, and adjust pointers
                    std::cerr << "Discovered text" << std::endl;
                    assert(in_buff_start == in_buff.get());
                    std::swap(in_buff, out_buff);
                    out_buff_free_start = in_buff_end;
                    in_buff_start = in_buff.get();
                    in_buff_end = in_buff.get();
                }
                else
                {
                    // run inflate() on input
                    if (! zstrm_p) zstrm_p = std::make_unique<detail::z_stream_wrapper>(true, Z_DEFAULT_COMPRESSION, window_bits);
                    zstrm_p->next_in = reinterpret_cast< decltype(zstrm_p->next_in) >(in_buff_start);
                    zstrm_p->avail_in = uint32_t(in_buff_end - in_buff_start);
                    zstrm_p->next_out = reinterpret_cast< decltype(zstrm_p->next_out) >(out_buff_free_start);
                    zstrm_p->avail_out = uint32_t((out_buff.get() + buff_size) - out_buff_free_start);
                    std::cerr << "Inflating avail_in=" << zstrm_p->avail_in << " to avail_out=" << zstrm_p->avail_out << " buf->in_avail()=" << this->sbuf_p->in_avail() << std::endl;
                    std::cerr << "First 24 bytes: " << (void*) *((int64_t*) zstrm_p->next_in) << " " << (void*) *((int64_t*) (zstrm_p->next_in+8)) << " " << (void*) *((int64_t*) (zstrm_p->next_in+16)) << std::endl;
                    if (is_bgzf && zstrm_p->total_in == 0 && zstrm_p->avail_in >= BGZF_BLOCK_HEADER_LENGTH-2 && memcmp(zstrm_p->next_in, BGZF_MAGIC_HEADER, BGZF_BLOCK_HEADER_LENGTH - 2) != 0) {
                        std::cerr << "WARNING: Autodetect failed... did not discover the BGZF magic header. Treating this as normal zip file..." << std::endl;
                        is_bgzf = false;
                    }
                    int ret = inflate(zstrm_p.get(), is_bgzf ? Z_SYNC_FLUSH : Z_NO_FLUSH);
                    // process return code
                    if (is_bgzf && ret == Z_STREAM_END ) {
                        // we finished a GZIP member
                        // scratch for peeking to see if the file is over
                        bool has_more = zstrm_p->avail_in > 0;
                        if (!has_more) {
                            int c = this->sbuf_p->snextc();
                            if (c != EOF) {
                                has_more = true;
                                this->sbuf_p->sputbackc((char) c);
                            }
                        }
                        if (has_more) {
                            //std::cerr << "input has more" << std::endl;
                            int reset_ret = inflateReset(zstrm_p.get());
                            if (reset_ret != Z_OK) {
                                throw Exception(zstrm_p.get(), ret);
                            }

                        } else {
                            // we consumed all the input data and hit Z_STREAM_END
                            // so stop looping, even if we never fill the output buffer
                            std::cerr << "input has NO more avail_in=" << zstrm_p->avail_in << " to avail_out=" << zstrm_p->avail_out << " buf->in_avail()=" << this->sbuf_p->in_avail() <<std::endl;
                            // stream ended, deallocate inflator
                            zstrm_p.reset();
                            break;
                        }
                    }
                    if (ret != Z_OK && ret != Z_STREAM_END) throw Exception(zstrm_p.get(), ret);
                    // update in&out pointers following inflate()
                    std::cerr << "After Inflating avail_in=" << zstrm_p->avail_in << " avail_out=" << zstrm_p->avail_out << std::endl;
                    in_buff_start = reinterpret_cast< decltype(in_buff_start) >(zstrm_p->next_in);
                    in_buff_end = in_buff_start + zstrm_p->avail_in;
                    out_buff_free_start = reinterpret_cast< decltype(out_buff_free_start) >(zstrm_p->next_out);
                    assert(out_buff_free_start + zstrm_p->avail_out == out_buff.get() + buff_size);

                    if (!is_bgzf && ret == Z_STREAM_END) {
                        // if stream ended, deallocate inflator
                        zstrm_p.reset();
                    }
                }
            } while (out_buff_free_start == out_buff.get());
            // 2 exit conditions:
            // - end of input: there might or might not be output available
            // - out_buff_free_start != out_buff: output available
            this->setg(out_buff.get(), out_buff.get(), out_buff_free_start);
        }
        return this->gptr() == this->egptr()
            ? traits_type::eof()
            : traits_type::to_int_type(*this->gptr());
    }

protected:
    std::streambuf * sbuf_p;
    std::unique_ptr<char[]> in_buff;
    char * in_buff_start;
    char * in_buff_end;
    std::unique_ptr<char[]> out_buff;
    std::unique_ptr<detail::z_stream_wrapper> zstrm_p;
    std::size_t buff_size;
    bool auto_detect;
    bool auto_detect_run;
    bool is_text;
    int window_bits;
    bool is_bgzf; // for bgzf format support

}; // class _istreambuf

class istreambuf : public _istreambuf 
{
public:
    istreambuf(std::streambuf * _sbuf_p,
               std::size_t _buff_size = default_buff_size, bool _auto_detect = true, int _window_bits = 0)
               : _istreambuf(_sbuf_p, _buff_size, _auto_detect, _window_bits) {
    } 
}; // class istreambuf

class bgzf_istreambuf : public _istreambuf
{
public:
    // use the buff size and flags for bgzf support
    bgzf_istreambuf(std::streambuf * _sbuf_p,
               std::size_t = 0, bool = true, int = 0)
               : _istreambuf(_sbuf_p, bgzf_default_buff_size, true, 15+16, true) {
    }

    void seek_to_bgzf_pointer(bgzf_virtual_file_pointer vfp) {
        assert(!this->zstrm_p || this->zstrm_p->avail_out == this->buff_size);
        std::cerr << "seek_to_bgzf_pointer gptr=" << (void*) this->gptr() << " egptr=" << (void*) this->egptr() << std::endl;
        auto pos = this->sbuf_p->pubseekpos(vfp.get_file_offset(), std::ios_base::in);
        if (pos != (std::streampos) vfp.get_file_offset()) {
            throw Exception(this->zstrm_p.get(), Z_ERRNO);
        }
        std::cerr << "seek_to_bgzf_pointer gptr=" << (void*) this->gptr() << " egptr=" << (void*) this->egptr() << std::endl;
        std::cerr << "Seeked in the compressed stream pos=" << pos << std::endl;
        this->zstrm_p.reset();
        auto uf = this->underflow();
        if (uf == traits_type::eof()) {
            throw Exception(zstrm_p.get(), Z_ERRNO, "Seek attempted past EOF");
        }
        
        auto block_offset = vfp.get_block_offset();
        if (block_offset) {
            if (!this->zstrm_p || block_offset > this->zstrm_p->avail_out) {
                throw Exception(zstrm_p.get(), Z_DATA_ERROR, " block_offset is larger than avail_out");
            }
            // start at the block offset
            this->zstrm_p->next_out += block_offset;
            this->zstrm_p->avail_out -= block_offset;
        }
        std::cerr << "successful underflow in the compressed stream pos=" << pos << " starting at block offset=" << block_offset << " avail_out=" << (this->zstrm_p ? this->zstrm_p->avail_out : -1) << std::endl;
    }

    bgzf_virtual_file_pointer find_next_bgzf_block(size_t compressed_offset) {
        std::cerr << "find_next_bgzf_block gptr=" << (void*) this->gptr() << " egptr=" << (void*) this->egptr() << std::endl;
        auto pos = this->sbuf_p->pubseekpos(compressed_offset, std::ios_base::in);
        std::cerr << "find_next_bgzf_block gptr=" << (void*) this->gptr() << " egptr=" << (void*) this->egptr() << std::endl;
        if (pos != (std::streampos) compressed_offset) {
            throw Exception(this->zstrm_p.get(), Z_ERRNO);
        }
        auto offset = this->find_next_bgzf_block2();
        if (offset >= 0) {
            return bgzf_virtual_file_pointer(compressed_offset + offset, 0);
        } else {
            return bgzf_virtual_file_pointer::get_invalid();
        }
    }
    std::streamsize find_next_bgzf_block2() {
        // read 2x the max buffer size as a block must be fully contained within that range
        auto test_buf = std::make_unique<char[]>(bgzf_default_buff_size * 2);
        std::streamsize sz = sbuf_p->sgetn(test_buf.get(), bgzf_default_buff_size * 2);
        std::cerr << "bgzf_istreambuf::find_next_bgzf_block2 sgetn got " << sz << std::endl;
        std::streamsize offset = 0;
        for( ; offset < (std::streamsize) bgzf_default_buff_size + 1; offset++) {
            if (offset > sz - BGZF_BLOCK_HEADER_LENGTH + BGZF_BLOCK_FOOTER_LENGTH) {
                // end of file
                break;
            }
            //std::cerr << "find_next_bgzf_block testing offset=" << offset << std::endl;
            // fail quickly looking for the next magic header
            if (memcmp(test_buf.get() + offset, BGZF_MAGIC_HEADER, BGZF_BLOCK_HEADER_LENGTH - 2) == 0) {
                
                // test for a sane compressed length within the data already read (last two bytes of BGZF_MAGIC_HEADER)
                uint16_t clen = packer::unpackInt16((uint8_t*)test_buf.get() + offset + BGZF_BLOCK_HEADER_LENGTH - 2);
                if (clen > sz - offset) continue; // erroneous length exceeds the EOF

                // possible hit, try to decompress without an error
                try {
                    this->zstrm_p = std::make_unique<detail::z_stream_wrapper>(true, Z_DEFAULT_COMPRESSION, window_bits);
                    this->zstrm_p->next_in = reinterpret_cast< decltype(this->zstrm_p->next_in) >(test_buf.get() + offset);
                    this->zstrm_p->avail_in = sz - offset;
                    this->zstrm_p->next_out = reinterpret_cast< decltype(this->zstrm_p->next_out) >(out_buff.get());
                    this->zstrm_p->avail_out = buff_size;
                    this->in_buff_start = reinterpret_cast< decltype(this->in_buff_start) >(this->zstrm_p->next_in);
                    this->in_buff_end = in_buff_start + this->zstrm_p->avail_in;
                    this->underflow();
                    // success! cleanup & return
                    // uncompressed output is lost, and one must call seek_to_bgzf_pointer to actually reposition the pointer.
                    std::cerr << "found header at +" << offset << std::endl;
                    return offset;

                } catch(...) {
                    // okay to rarely fail
                }
            }
        }
        this->reset();
        if (sz < (std::streamsize) bgzf_default_buff_size * 2) {
            // not found but at eof
            std::cerr << "found EOF at +" << sz << std::endl;
            return sz;
        } else {
            // not found, so the bgzf stream is in error
            throw Exception(this->zstrm_p.get(), Z_DATA_ERROR, " No BGZF block found");
        }
    }


}; // bgzf_istreambuf

class _ostreambuf
    : public std::streambuf
{
public:
    _ostreambuf(std::streambuf * _sbuf_p,
               std::size_t _buff_size, int _level, int _window_bits, bool _is_bgzf = false)
        : sbuf_p(_sbuf_p),
          zstrm_p(std::make_unique<detail::z_stream_wrapper>(false, _level, _window_bits)),
          buff_size(_buff_size),
          level(_level),
          is_bgzf(_is_bgzf)
    {
        assert(sbuf_p);
        //std::cerr << "ostreambuf constuctor buff_size=" << buff_size << " window_bits=" << _window_bits << " is_bgzf=" << _is_bgzf << std::endl;
        in_buff = std::make_unique<char[]>(buff_size);
        out_buff = std::make_unique<char[]>(buff_size);
        setp(in_buff.get(), in_buff.get() + (is_bgzf ? BGZF_BLOCK_SIZE : buff_size));
    }

    _ostreambuf(const _ostreambuf &) = delete;
    _ostreambuf & operator = (const _ostreambuf &) = delete;

    std::streamsize write_sink() {
        std::streamsize write_size = (zstrm_p && zstrm_p->next_out) ? reinterpret_cast< decltype(out_buff.get()) >(zstrm_p->next_out) - out_buff.get() : 0;
        if (write_size > 0) {
            assert(zstrm_p && zstrm_p->next_out);
            assert(out_buff.get() + write_size == reinterpret_cast< decltype(out_buff.get()) >(zstrm_p->next_out));
            //std::cerr << "Writing " << (buff_size - zstrm_p->avail_out) << " == " << write_size << std::endl;
            std::streamsize sz = sbuf_p->sputn(out_buff.get(), write_size);
            if (sz != write_size)
            {
                throw Exception(zstrm_p.get(), Z_ERRNO);
            }
        } else {
            //std::cerr << "Nothing to write zstrm_p=" << zstrm_p.get() << " out_buff=" << (void*) out_buff.get() << " next_out=" << (void*) (zstrm_p ? zstrm_p->next_out : 0) << " avail_out=" << (zstrm_p ? buff_size - zstrm_p->avail_out : -1)<< std::endl;
        }
        
        zstrm_p->next_out = reinterpret_cast< decltype(zstrm_p->next_out) >(out_buff.get());
        zstrm_p->avail_out = buff_size;
        return write_size;
    }

    int deflate_loop(int flush)
    {
        while (true)
        {
            std::streamsize sz = write_sink();
            //std::cerr << "Looping flush=" << flush << " with sz=" << sz << " avail_in=" << zstrm_p->avail_in << " flush=" << flush << " next_out=" << (size_t) (reinterpret_cast< decltype(out_buff.get()) >(zstrm_p->next_out) - out_buff.get()) << std::endl;
            zstrm_p->next_out = reinterpret_cast< decltype(zstrm_p->next_out) >(out_buff.get());
            zstrm_p->avail_out = uint32_t(buff_size);
            int ret;
            if (is_bgzf) {
                ret = bgzf_deflate(zstrm_p.get());
            } else {
                ret = deflate(zstrm_p.get(), flush);
            }
            //std::cerr << "avail_in=" << zstrm_p->avail_in << " next_out=" << buff_size - zstrm_p->avail_out << std::endl;
            if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR) {
                failed = true;
                throw Exception(zstrm_p.get(), ret);
            }
            sz = write_sink();
            if (ret == Z_STREAM_END || ret == Z_BUF_ERROR || sz == 0)
            {
                break;
            }
        }
        return 0;
    }

    virtual ~_ostreambuf()
    {
        // flush the zlib stream
        //
        // NOTE: Errors here (sync() return value not 0) are ignored, because we
        // cannot throw in a destructor. This mirrors the behaviour of
        // std::basic_filebuf::~basic_filebuf(). To see an exception on error,
        // close the ofstream with an explicit call to close(), and do not rely
        // on the implicit call in the destructor.
        //
        if (!failed) try {
            sync();
            if (is_bgzf) {
                // write one last empty block
                //std::cerr << "Writing last empty block for BGZF" << std::endl;
                bgzf_deflate(zstrm_p.get());
                std::streamsize sz = write_sink();
                assert(sz == 28 && "last empty block is 28 bytes");
            }
        } catch (...) {}
    }
    std::streambuf::int_type overflow(std::streambuf::int_type c = traits_type::eof()) override
    {
        //std::cerr << "overflow" << std::endl;
        zstrm_p->next_in = reinterpret_cast< decltype(zstrm_p->next_in) >(pbase());
        zstrm_p->avail_in = uint32_t(pptr() - pbase());
        while (zstrm_p->avail_in > 0)
        {
            //std::cerr << "entering deflate_loop avail_in=" << zstrm_p->avail_in << std::endl;
            int r = deflate_loop(Z_NO_FLUSH);
            if (r != 0)
            {
                setp(nullptr, nullptr);
                return traits_type::eof();
            }
        }
        setp(in_buff.get(), in_buff.get() + (is_bgzf ? BGZF_BLOCK_SIZE : buff_size));
        return traits_type::eq_int_type(c, traits_type::eof()) ? traits_type::eof() : sputc(char_type(c));
    }
    int sync() override
    {
        //std::cerr << "sync() caled" << std::endl;
        // first, call overflow to clear in_buff
        overflow();
        if (! pptr()) return -1;
        // then, call deflate asking to finish the zlib stream
        //zstrm_p->next_in = nullptr;
        //std::cerr << "sync, finishing loop again" << std::endl;
        assert(zstrm_p->avail_in == 0);
        //zstrm_p->avail_in = 0;
        if (deflate_loop(Z_FINISH) != 0) return -1;
        deflateReset(zstrm_p.get());
        write_sink();
        return 0;
    }

    int bgzf_deflate(z_stream *zs) {
        // every call results in a full block with header and footer
        
        zs->zalloc = NULL;
        zs->zfree = NULL;
        zs->msg = NULL;
    
        // cap input to BGZF_BLOCK_SIZE
        auto slen = zs->avail_in;
        if (slen > BGZF_BLOCK_SIZE) {
            slen = BGZF_BLOCK_SIZE; // maximum length allowed here
            zs->avail_in = slen;
        }
        if (zs->avail_in == 0 and zs->next_in == nullptr) {
            zs->next_in = reinterpret_cast< decltype(this->zstrm_p->next_in) >(in_buff.get());
        }
        auto src = zs->next_in;

        uint8_t*dst = reinterpret_cast< decltype(zs->next_out) >(zs->next_out);
        auto orig_avail_out = zs->avail_out;
        //std::cerr << "orig_avail_out=" << orig_avail_out << std::endl;
        assert(orig_avail_out >= BGZF_MAX_BLOCK_SIZE);
        // offset the actual writing by the header & avail_out by the header & footer lengths
        zs->next_out = dst + BGZF_BLOCK_HEADER_LENGTH;
        uInt start_out = BGZF_MAX_BLOCK_SIZE - BGZF_BLOCK_HEADER_LENGTH - BGZF_BLOCK_FOOTER_LENGTH;
        zs->avail_out = start_out;

        // initialize the new block
        int ret = deflateInit2(zs, this->level, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); // -15 to disable zlib header/footer manged here
        if (ret!=Z_OK) {
            throw Exception(zstrm_p.get(), ret);
        }
        // compress the body
        if ((ret = deflate(zs, Z_FINISH)) != Z_STREAM_END) {
            throw Exception(zstrm_p.get(), ret);
        }
        assert(zs->avail_in == 0 && "The entire available input was compressed");
        assert(zs->next_in == src + slen && "next_in is at the end");
        auto write_len = start_out - zs->avail_out;
        assert(zs->total_out == write_len);
        uInt dlen = write_len + BGZF_BLOCK_HEADER_LENGTH + BGZF_BLOCK_FOOTER_LENGTH;
        assert(dlen <= 1<<16);
        assert(zs->next_out == dst + BGZF_BLOCK_HEADER_LENGTH + write_len);

        // finish the block
        if ((ret = deflateEnd(zs)) != Z_OK) {
            throw Exception(zstrm_p.get(), ret);
        }
        

        // write the header
        memcpy(dst, BGZF_MAGIC_HEADER, BGZF_BLOCK_HEADER_LENGTH); // the last two bytes are a place holder for the length of the block
        packer::packInt16(&dst[16], dlen - 1); // write the compressed length; -1 to fit 2 bytes
        
        // write the footer with checksum and uncompressed length
        uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)src, slen);
        packer::packInt32((uint8_t*)&dst[dlen - 8], crc);
        packer::packInt32((uint8_t*)&dst[dlen - 4], slen);

        // cleanup zs
        //std::cerr << "Deflated " << slen << " bytes, wrote " << write_len << " plus header/footer to " << dlen << std::endl;
        zs->next_out += BGZF_BLOCK_FOOTER_LENGTH;
        zs->avail_out = orig_avail_out - dlen;
        zs->total_out += BGZF_BLOCK_HEADER_LENGTH + BGZF_BLOCK_FOOTER_LENGTH;
        assert(zs->next_out == dst + dlen);
        index.append_incremental_block( slen, dlen );
        std::cerr << "Block uncomp:" << slen << " comp:" << dlen << std::endl;
        
        return Z_STREAM_END;

    }

protected:
    std::streambuf * sbuf_p = nullptr;
    std::unique_ptr<char[]> in_buff;
    std::unique_ptr<char[]> out_buff;
    std::unique_ptr<detail::z_stream_wrapper> zstrm_p;
    std::size_t buff_size;
    int level;
    bool failed = false;
    bool is_bgzf; // for bgzf format
    bgzf_index index;

}; // class _ostreambuf

class ostreambuf : public _ostreambuf 
{
public:
    ostreambuf(std::streambuf * _sbuf_p,
               std::size_t _buff_size = default_buff_size, int _level = Z_DEFAULT_COMPRESSION, int _window_bits = 0)
    : _ostreambuf(_sbuf_p, _buff_size, _level, _window_bits) {
        assert(!is_bgzf);
    }
}; // class ostreambuf

class bgzf_ostreambuf : public _ostreambuf
{
public:
    bgzf_ostreambuf(std::streambuf * _sbuf_p,
               std::size_t = 0, int _level = Z_DEFAULT_COMPRESSION, int = 0) 
    : _ostreambuf(_sbuf_p, bgzf_default_buff_size, _level, 15+16, true) {
        assert(this->is_bgzf);
    }

    void output_index(std::ostream &os, pos_type uncompressed_offset = 0, pos_type compressed_offset = 0) {
        this->index.rebaseline(uncompressed_offset, compressed_offset);
        this->index.write_index(os);
    }

}; // class bgzf_ostreambuf

template<typename _istreambuf>
class _istream : public std::istream
{
public:
    _istream(std::istream & is, std::size_t _buff_size, bool _auto_detect, int _window_bits)
        : std::istream(new _istreambuf(is.rdbuf(), _buff_size, _auto_detect, _window_bits))
    {
        exceptions(std::ios_base::badbit);
    }
    explicit _istream(std::streambuf * sbuf_p)
        : std::istream(new _istreambuf(sbuf_p))
    {
        exceptions(std::ios_base::badbit);
    }
    virtual ~_istream()
    {
        delete rdbuf();
    }
}; // class _istream

class istream : public _istream<istreambuf> 
{
public:
    istream(std::istream & is,
            std::size_t _buff_size = default_buff_size, bool _auto_detect = true, int _window_bits = 0)
        : _istream(is, _buff_size, _auto_detect, _window_bits) {}
}; // class istream

class bgzf_istream : public _istream<bgzf_istreambuf> 
{
public:
    bgzf_istream(std::istream & is,
            std::size_t = 0, bool _auto_detect = true, int = 0)
        : _istream(is, bgzf_default_buff_size, _auto_detect, 15+16) {}

    static bool is_bgzf(std::istream &is) {
        // peeks at the stream, reading the first few bytes for the magic bgzf header
        char buf[19];
        if (is.fail() || is.eof()) return false;
        is.read(buf,16);
        if (is.fail() || is.eof()) return false;
        bool ret = memcmp(buf, BGZF_MAGIC_HEADER, 16) == 0;
        // put back those 16 bytes
        for (int i = 0; i < 16; i++)
            is.putback(buf[16-1-i]);
        return ret;
    }

}; // class bgzf_istream

template<typename _ostreambuf>
class _ostream : public std::ostream
{
public:
    _ostream(std::ostream & os,
            std::size_t _buff_size = default_buff_size, int _level = Z_DEFAULT_COMPRESSION, int _window_bits = 0)
        : std::ostream(new _ostreambuf(os.rdbuf(), _buff_size, _level, _window_bits))
    {
        exceptions(std::ios_base::badbit);
    }
    explicit _ostream(std::streambuf * sbuf_p)
        : std::ostream(new _ostreambuf(sbuf_p))
    {
        exceptions(std::ios_base::badbit);
    }
    virtual ~_ostream()
    {
        delete rdbuf();
    }
}; // class _ostream

class ostream : public _ostream<ostreambuf> 
{
public:
    ostream(std::ostream & os,
            std::size_t _buff_size = default_buff_size, int _level = Z_DEFAULT_COMPRESSION, int _window_bits = 0)
        : _ostream(os, _buff_size, _level, _window_bits) {}
}; // class ostream

class bgzf_ostream : public _ostream<bgzf_ostreambuf> {
public:
    bgzf_ostream(std::ostream & os,
            std::size_t = 0, int _level = Z_DEFAULT_COMPRESSION, int = 0)
        : _ostream(os, bgzf_default_buff_size, _level, 15+16) {}
}; // class bgzf_ostream


namespace detail
{

template < typename FStream_Type >
struct strict_fstream_holder
{
    strict_fstream_holder(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in)
        : _fs(filename, mode)
    {}
    FStream_Type _fs;
}; // class strict_fstream_holder

} // namespace detail

template<typename _istreambuf>
class _ifstream
    : private detail::strict_fstream_holder< strict_fstream::ifstream >,
      public std::istream
{
public:
    explicit _ifstream(const std::string filename, std::ios_base::openmode mode = std::ios_base::in, size_t buff_size = _istreambuf::default_buff_size)
        : detail::strict_fstream_holder< strict_fstream::ifstream >(filename, mode),
          std::istream(new _istreambuf(_fs.rdbuf(), buff_size))
    {
        exceptions(std::ios_base::badbit);
    }
    virtual ~_ifstream()
    {
        if (_fs.is_open()) close();
        if (rdbuf()) delete rdbuf();
    }
    virtual void close()
    {
        _fs.close();
    }
    // Return the position within the compressed file
    std::streampos tellg()
    {
        return _fs.tellg();
    }
    std::streampos compressed_tellg()
    {
        return tellg();
    }
    _ifstream & seekg(std::streampos pos)
    {
        this->clear();
        _fs.clear();
        _fs.seekg(pos);
        return *this;
    }
}; // class ifstream


class ifstream : public _ifstream<istreambuf> 
{
public:
    explicit ifstream(const std::string filename, std::ios_base::openmode mode = std::ios_base::in, size_t buff_size = default_buff_size)
        : _ifstream(filename, mode, buff_size)
    {
        exceptions(std::ios_base::badbit);
    }
};

class bgzf_ifstream : public _ifstream<bgzf_istreambuf> 
{
public:
    explicit bgzf_ifstream(const std::string filename, std::ios_base::openmode mode = std::ios_base::in, size_t = 0)
        : _ifstream(filename, mode, bgzf_default_buff_size)
    {
        exceptions(std::ios_base::badbit);
    }

    void seek_to_bgzf_pointer(bgzf_virtual_file_pointer vfp) {
        std::cerr << "bgzf_ifstream::seek_to_bgzf_pointer Seeking to " << vfp.get_file_offset() << " at " << vfp.get_block_offset() << std::endl;
        bgzf_istreambuf * bgzf_isb = (bgzf_istreambuf *) rdbuf();
        this->clear();
        this->seekg(vfp.get_file_offset());
        bgzf_isb->reset();
        bgzf_isb->underflow();
        auto offset = vfp.get_block_offset();
        if (offset > 0) {
            auto buf = std::make_unique<char[]>(offset+1);
            this->read(buf.get(), offset);
            if (offset != this->gcount()) {
                throw Exception(nullptr, Z_DATA_ERROR, "Could not read to block offset");
            }
        }
        std::cerr << "seek_to_bgzf_pointer at now" << std::endl;
        //bgzf_isb->seek_to_bgzf_pointer(vfp);
    }

    bgzf_virtual_file_pointer find_next_bgzf_block2(size_t compressed_offset) {
        bgzf_istreambuf * bgzf_isb = (bgzf_istreambuf *) rdbuf();
        return bgzf_isb->find_next_bgzf_block(compressed_offset);
    }

    bgzf_virtual_file_pointer find_next_bgzf_block(size_t compressed_offset) {

        std::cerr << "bgzf_ifstream::find_next_bgzf_block(" << compressed_offset << ")" << std::endl;
        this->clear();
        this->seekg(compressed_offset);

        bgzf_istreambuf * bgzf_isb = (bgzf_istreambuf *) rdbuf();
        
        auto offset = bgzf_isb->find_next_bgzf_block2();
        if (offset >= 0) {
            std::cerr << "Found valid block at " << compressed_offset + offset << std::endl;
            return bgzf_virtual_file_pointer(compressed_offset + offset, 0);
        } else {
            std::cerr << "DID NOT FIND valid block after " << compressed_offset << std::endl;
            return bgzf_virtual_file_pointer::get_invalid();
        }
        
    }

};

template<typename base_ostreambuf>
class _ofstream
    : private detail::strict_fstream_holder< strict_fstream::ofstream >,
      public std::ostream
{
public:
    explicit _ofstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out,
                      int level = Z_DEFAULT_COMPRESSION, size_t buff_size = default_buff_size)
        : detail::strict_fstream_holder< strict_fstream::ofstream >(filename, mode | std::ios_base::binary),
          std::ostream(new base_ostreambuf(_fs.rdbuf(), buff_size, level))
    {
        exceptions(std::ios_base::badbit);
    }
    _ofstream& flush() {
        std::ostream::flush();
        _fs.flush();
        return *this;
    }
    virtual ~_ofstream()
    {
        if (_fs.is_open()) close();
        if (rdbuf()) delete rdbuf();
    }
    virtual void close()
    {
        std::ostream::flush();
        _fs.close();
    }
    // Return the position within the compressed file
    std::streampos compressed_tellp()
    {
        return _fs.tellp();
    }

}; // class ofstream

class ofstream : public _ofstream<ostreambuf> 
{
public:
    explicit ofstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out,
                      int level = Z_DEFAULT_COMPRESSION, size_t buff_size = default_buff_size)
        : _ofstream(filename, mode, level, buff_size) {

        }
};

class bgzf_ofstream : public _ofstream<bgzf_ostreambuf> 
{
public:
        explicit bgzf_ofstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out,
                      int level = Z_DEFAULT_COMPRESSION)
        : _ofstream(filename, mode, level, bgzf_default_buff_size) {

        }
        void output_index(std::ostream &os, pos_type uncompressed_offset = 0, pos_type compressed_offset = 0) {
            bgzf_ostreambuf *bgzf_osbuf = (bgzf_ostreambuf *) this->rdbuf();
            bgzf_osbuf->output_index(os, uncompressed_offset, compressed_offset);
        }
    
};

} // namespace zstr

