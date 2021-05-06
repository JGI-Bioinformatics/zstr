.SUFFIXES:
MAKEFLAGS += -r
SHELL := /bin/bash
.DELETE_ON_ERROR:
.PHONY: all test clean

all: test-strict_fstream ztxtpipe zpipe zc

%: %.cpp
	g++ -std=c++11 -O0 -g3 -ggdb -fno-eliminate-unused-debug-types -Wall -Wextra -pedantic -I../src -o $@ $^ -lz

test: ztxtpipe zpipe zc bgzc bgzf_random_access
	diff -q ztxtpipe.cpp <(./ztxtpipe <ztxtpipe.cpp)
	diff -q ztxtpipe.cpp <(gzip <ztxtpipe.cpp | ./ztxtpipe)
	diff -q ztxtpipe.cpp <(gzip -9 <ztxtpipe.cpp | ./ztxtpipe)
	diff -q <(cat ztxtpipe.cpp ztxtpipe.cpp) <(cat ztxtpipe.cpp ztxtpipe.cpp | gzip | ./ztxtpipe)
	diff -q <(cat ztxtpipe.cpp ztxtpipe.cpp) <({ gzip <ztxtpipe.cpp; gzip <ztxtpipe.cpp; } | ./ztxtpipe)
	diff -q zpipe.cpp <(./zpipe <zpipe.cpp)
	diff -q zpipe.cpp <(gzip <zpipe.cpp | ./zpipe)
	diff -q zpipe.cpp <(gzip -9 <zpipe.cpp | ./zpipe)
	diff -q <(cat zpipe.cpp zpipe.cpp) <(cat zpipe.cpp zpipe.cpp | gzip | ./zpipe)
	diff -q <(cat zpipe.cpp zpipe.cpp) <({ gzip <zpipe.cpp; gzip <zpipe.cpp; } | ./zpipe)
	diff -q <(<zpipe.cpp | gzip) <(<zpipe.cpp | gzip | gzip | ./zpipe)
	diff -q zc.cpp <(./zc <zc.cpp)
	diff -q zc.cpp <(./zc - <zc.cpp)
	diff -q zc.cpp <(./zc - - <zc.cpp)
	diff -q zc.cpp <(./zc zc.cpp)
	diff -q zc.cpp <(<zc.cpp gzip | ./zc)
	diff -q zc.cpp <(<zc.cpp gzip | ./zc -)
	diff -q zc.cpp <(<zc.cpp gzip | ./zc - -)
	diff -q zc.cpp <(./zc <(<zc.cpp gzip))
	diff -q zc.cpp <(<zc.cpp ./zc -c | zcat)
	diff -q zc.cpp <(<zc.cpp ./zc -c - | zcat)
	diff -q zc.cpp <(<zc.cpp ./zc -c - - | zcat)
	diff -q zc.cpp <(./zc -c zc.cpp | zcat)
	diff -q <(cat zc.cpp zc.cpp) <(./zc -c zc.cpp - <zc.cpp | zcat)
	diff -q <(cat zc.cpp zc.cpp) <(./zc -c - zc.cpp <zc.cpp | zcat)
	diff -q <(cat zc.cpp zc.cpp) <(./zc -c zc.cpp zc.cpp | zcat)
	diff -q <(cat zc.cpp zc.cpp) <({ ./zc -c zc.cpp; ./zc -c zc.cpp; } | zcat)
	diff -q <(cat zc.cpp zc.cpp) <({ gzip <zc.cpp; ./zc -c zc.cpp; } | zcat)
	diff -q <(cat zc.cpp zc.cpp) <({ ./zc -c zc.cpp; gzip <zc.cpp; } | zcat)
	tar -cf .tmp ztxtpipe zpipe zc bgzc
	bgzip -f -c .tmp >.tmp.gz
	gzip -t .tmp.gz
	diff -q .tmp <(./zc <.tmp.gz)
	diff -q .tmp <(./bgzc <.tmp.gz)
	diff -q .tmp <(./bgzc .tmp.gz)
	./bgzc -c .tmp > .tmp2.gz
	gzip -t .tmp2.gz
	diff -q .tmp <(./zc <.tmp2.gz)
	diff -q .tmp <(./bgzc <.tmp2.gz)
	diff -q .tmp <(./bgzc .tmp2.gz)
	gzip -c .tmp > .tmp3.gz
	gzip -t .tmp3.gz
	diff -q .tmp <(./zc <.tmp3.gz)
	diff -q .tmp <(./bgzc <.tmp3.gz) && echo "Expected a BGZF magic header warning."
	diff -q .tmp <(./bgzc .tmp3.gz) && echo "Expected a BGZF magic header warning."
	./bgzf_random_access .tmp.gz
	./bgzf_random_access .tmp2.gz
	./bgzf_random_access .tmp3.gz || echo "Expected that to fail"
	rm -f .tmp*
	@echo "all passed"

clean:
	rm -rf test-strict_fstream ztxtpipe zpipe zc bgzc bgzf_random_access
