#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include "zstr.hpp"

void usage(std::ostream& os, const std::string& prog_name)
{
    os << "Use: " << prog_name << " bgzf.gz ..." << std::endl
       << "Synposis:" << std::endl
       << "  Reads the bgzf file and verifies the find_next_bgzf_block and seek_to_bgzf_pointer methods work" << std::endl;
}

void cat_stream(std::istream& is, std::ostream& os)
{
    const std::streamsize buff_size = 1 << 16;
    char * buff = new char [buff_size];
    while (true)
    {
        is.read(buff, buff_size);
        std::streamsize cnt = is.gcount();
        if (cnt == 0) break;
        os.write(buff, cnt);
    }
    delete [] buff;
} // cat_stream

void decompress_files(const std::vector< std::string >& file_v, const std::string& output_file)
{
    //
    // Set up sink ostream
    //
    std::unique_ptr< std::ofstream > ofs_p;
    std::ostream * os_p = &std::cout;
    if (not output_file.empty())
    {
        ofs_p = std::unique_ptr< std::ofstream >(new strict_fstream::ofstream(output_file));
        os_p = ofs_p.get();
    }
    //
    // Process files
    //
    for (const auto& f : file_v)
    {
        //
        // If `f` is a file, create a zstr::ifstream, else (it is stdin) create a zstr::istream wrapper
        //
        std::unique_ptr< std::istream > is_p =
            (f != "-"
             ? std::unique_ptr< std::istream >(new zstr::bgzf_ifstream(f))
             : std::unique_ptr< std::istream >(new zstr::bgzf_istream(std::cin)));
        //
        // Cat stream
        //
        cat_stream(*is_p, *os_p);
    }
} // decompress_files

bool verify_random_access(const std::string &f) {
    // first find all the magic blocks in the file and index it
    std::ifstream raw(f);
    zstr::bgzf_index idx, idx2;
    char buf[zstr::BGZF_BLOCK_HEADER_LENGTH + zstr::BGZF_BLOCK_FOOTER_LENGTH + 1];
    bool is_verified = true;
    uint64_t uncmp_pos = 0;
    uint64_t cmp_pos = 0;
    std::cerr << "Scanning " << f << " for BGZF blocks" << std::endl;
    raw.read(buf,zstr::BGZF_BLOCK_HEADER_LENGTH);
    auto bytes = raw.gcount();
    while(is_verified) {
        
        if (bytes == 0) {
            // try reading more
            raw.read(buf,1);
            bytes = raw.gcount();
            if (raw.eof()) { std::cerr << "Found natural end at " << raw.tellg() << std::endl; break; }
            if (bytes > 0 ) { is_verified = false; std::cerr << "Found unnatural end at " << raw.tellg() << std::endl; break;}
            break;
        }
        if (raw.fail()) { is_verified = false; std::cerr << "Got failure reading BGZF header" << std::endl; break; }
        if (bytes != 18) { is_verified = false; std::cerr << "Did not read full BGZF header" << std::endl; break; }
        if (memcmp(buf, zstr::BGZF_MAGIC_HEADER, zstr::BGZF_BLOCK_HEADER_LENGTH-2) != 0) { is_verified = false; std::cerr << "Block header mismatch at " << raw.tellg() << std::endl; break; }
        auto cmp_len = zstr::packer::unpackInt16((uint8_t*)buf + zstr::BGZF_BLOCK_HEADER_LENGTH - 2) + 1;
        auto seek_pos = cmp_pos + cmp_len - zstr::BGZF_BLOCK_FOOTER_LENGTH;
        //std::cerr << "block scanning Seeking to " << seek_pos << " cmp_len=" << cmp_len << std::endl;
        raw.seekg(seek_pos);
        assert((size_t) raw.tellg() == seek_pos);
        raw.read(buf,zstr::BGZF_BLOCK_FOOTER_LENGTH + zstr::BGZF_BLOCK_HEADER_LENGTH);
        bytes = raw.gcount();
        
        if (bytes < zstr::BGZF_BLOCK_FOOTER_LENGTH) { is_verified = false; std::cerr << "Did not read full BGZF footer" << std::endl; break; }
        auto uncmp_len = zstr::packer::unpackInt32((uint8_t*)buf + zstr::BGZF_BLOCK_FOOTER_LENGTH - 4);
        //std::cerr << "at " << cmp_pos << " uncomp=" << uncmp_len << std::endl;
        cmp_pos += cmp_len;
        uncmp_pos += uncmp_len;
        idx.append_absolute_block(uncmp_pos, cmp_pos);
        idx2.append_incremental_block(uncmp_len, cmp_len);
        //std::cerr << "at " << cmp_pos << " uncomp=" << uncmp_len << std::endl;
        if (bytes == zstr::BGZF_BLOCK_FOOTER_LENGTH + zstr::BGZF_BLOCK_HEADER_LENGTH) {
            // shift
            memmove(buf, buf + zstr::BGZF_BLOCK_FOOTER_LENGTH, zstr::BGZF_BLOCK_HEADER_LENGTH);
            bytes = zstr::BGZF_BLOCK_HEADER_LENGTH;
        } else if (bytes == zstr::BGZF_BLOCK_FOOTER_LENGTH) {
            // natural end
            break;
        }
        if (raw.fail()) { is_verified = false; std::cerr << "Got failure reading BGZF footer" << std::endl; break; }
    }
    if (!is_verified) {
        std::cerr << f << " is not a valid BGZF file!" << std::endl;
        return is_verified;
    }
    std::cout << "Found " << idx.size() << " blocks with total uncompressed size of " << uncmp_pos << std::endl;
    raw.close();
    assert(idx.size() == idx2.size());
    auto it1 = idx.begin();
    auto it2 = idx2.begin();
    for(size_t i = 0; i < idx.size(); i++) {
        assert(it1->caddr == it2->caddr);
        assert(it1->uaddr == it2->uaddr);
        it1++; it2++;
    }

    std::ostringstream oss;

    if (!zstr::bgzf_ifstream::is_bgzf(f)) {
        std::cerr << f << " is NOT a BGZF file" << std::endl;
    }

    zstr::bgzf_ifstream in(f);
    assert((size_t)oss.tellp() == 0);
    //std::cerr << "try @0 in.tellg()=" << in.tellg() << " uncomp.tellp()=" << oss.tellp()<< std::endl;
    assert((size_t)in.tellg() == 0);
    cat_stream(in, oss);
    //std::cerr << "in.tellg()=" << in.tellg() << " uncomp.tellp()=" << oss.tellp()<< std::endl;
    assert((size_t)oss.tellp() == uncmp_pos);
    assert((size_t)in.tellg() == cmp_pos);
    
    in.seekg(0);
    std::ostringstream().swap(oss);
    //std::cerr << "try @0 again in.tellg()=" << in.tellg() << " uncomp.tellp()=" << oss.tellp() << std::endl;
    assert((size_t)oss.tellp() == 0);
    assert((size_t)in.tellg() == 0);
    cat_stream(in, oss);
    //std::cerr << "in.tellg()=" << in.tellg() << " uncomp.tellp()=" << oss.tellp() << std::endl;
    assert((size_t)oss.tellp() == uncmp_pos);
    assert((size_t)in.tellg() == cmp_pos);
    

    //in.seekg(0);
    //std::ostringstream().swap(oss);
    //std::cerr << "try partial in.tellg()=" << in.tellg() << " uncomp.tellp()=" << oss.tellp()<< std::endl;

    // now verify the blocks can be found
    auto last = idx.end();
    int blk_i = 0;
    for(auto iter = idx.begin(); iter != idx.end(); iter++) {
        const auto &blk = *iter;
        assert(blk.uaddr <= uncmp_pos);
        // test exact
        //std::cerr << "Looking at blk.uaddr=" << blk.uaddr << " blk.caddr=" << blk.caddr << std::endl;
        auto vfp = in.find_next_bgzf_block(blk.caddr);
        //auto vfp = zstr::bgzf_virtual_file_pointer(blk.caddr, 0);
        //assert(vfp == vfp2);
        //std::cerr << "Found vfp=" << vfp.get_file_offset() << " at " << vfp.get_block_offset() << " after pos=" << blk.caddr << std::endl;
        if (vfp.get_file_offset() != blk.caddr) { is_verified = false; std::cerr << "Expected first byte of block did not match!" << std::endl; break; }
        if (vfp.get_block_offset() != 0) { is_verified = false; std::cerr << "Expected first byte of block did not match block offset!" << std::endl; break; }
    
        
        if (last != idx.end() && iter->caddr > last->caddr) {
            // test 1 over previous block
            std::vector<size_t> attempts;
            attempts.push_back(last->caddr+1);
            attempts.push_back((last->caddr + iter->caddr) / 2);
            attempts.push_back(iter->caddr - 1);
            for (auto attempt : attempts) {
                auto test_vfp = in.find_next_bgzf_block(attempt);
                if (test_vfp != vfp) { is_verified = false; std::cerr << "Expected find within previous block to find this one!" << std::endl; break;}
            }
        }
        
        //assert((size_t)in.tellg() == vfp.get_file_offset());
        std::ostringstream().swap(oss);
        assert((uint64_t) oss.tellp() == 0);
        //std::cerr << "Seeking to " << vfp.get_file_offset() << " at " << vfp.get_block_offset() << std::endl;
        in.seek_to_bgzf_pointer(vfp);
        //std::cerr << "Checking block " << blk_i << " at cmp_offset " << vfp.get_file_offset() << " at " << vfp.get_block_offset() << std::endl;
        cat_stream(in, oss);
        //std::cerr << "in.tellg()=" << in.tellg() << " uncomp.tellp()=" << oss.tellp() << std::endl;
        if ((uint64_t) oss.tellp() != uncmp_pos - blk.uaddr) { is_verified = false; std::cerr << "Got wrong remainder of uncompressed stream. expected " << uncmp_pos - blk.uaddr << " got " << oss.tellp() << " blk.uaddr=" << blk.uaddr << std::endl; break; }

        if (oss.tellp() > 0) {
            vfp = zstr::bgzf_virtual_file_pointer(vfp.get_file_offset(), 1);
            //std::cerr << "Seeking to " << vfp.get_file_offset() << " at " << vfp.get_block_offset() << std::endl;
            in.seek_to_bgzf_pointer(vfp);
            //std::cerr << "Checking block " << blk_i << " at cmp_offset " << vfp.get_file_offset() << " at " << vfp.get_block_offset() << std::endl;
            std::ostringstream().swap(oss);
            assert((uint64_t) oss.tellp() == 0);
            cat_stream(in, oss);
            if ((uint64_t) oss.tellp() != uncmp_pos - blk.uaddr - 1) { is_verified = false; std::cerr << "2Got wrong remainder of uncompressed stream. expected " << uncmp_pos - blk.uaddr << " got " << oss.tellp()  << std::endl; break; }

        }

        last = iter;
        blk_i++;
    }
    return is_verified;
}

void compress_files(const std::vector< std::string >& file_v, const std::string& output_file)
{
    //
    // Set up compression sink ostream
    //
    std::unique_ptr< std::ostream > os_p =
        (not output_file.empty()
         ? std::unique_ptr< std::ostream >(new zstr::bgzf_ofstream(output_file))
         : std::unique_ptr< std::ostream >(new zstr::bgzf_ostream(std::cout)));
    //
    // Process files
    //
    for (const auto& f : file_v)
    {
        //
        // If `f` is a file, create an ifstream, else read stdin
        //
        std::unique_ptr< std::ifstream > ifs_p;
        std::istream * is_p = &std::cin;
        if (f != "-")
        {
            ifs_p = std::unique_ptr< std::ifstream >(new strict_fstream::ifstream(f));
            is_p = ifs_p.get();
        }
        //
        // Cat stream
        //
        cat_stream(*is_p, *os_p);
    }
} // compress_files

int main(int argc, char * argv[])
{
    std::string output_file;
    int c;
    while ((c = getopt(argc, argv, "h?")) != -1)
    {
        switch (c)
        {
        case '?':
        case 'h':
            usage(std::cout, argv[0]);
            std::exit(EXIT_SUCCESS);
            break;
        default:
            usage(std::cerr, argv[0]);
            std::exit(EXIT_FAILURE);
        }
    }
    //
    // Gather files to process
    //
    std::vector< std::string > file_v(&argv[optind], &argv[argc]);
    //
    // With no other arguments, process stdin
    //
    if (file_v.empty()) file_v.push_back("-");
    //
    // Perform compression/decompression
    //
    decompress_files(file_v, "/dev/null");

    int ret = 0;
    for(auto &f : file_v)
        ret |= verify_random_access(f) ? 0 : 1;
    return ret;
}
