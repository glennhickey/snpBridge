# this Makefile is largely derived from https://github.com/adamnovak/corg/blob/master/Makefile
.PHONY: all clean

CXX=g++
INCLUDES=-Ivg -Ivg/gssw/src -Ivg/protobuf/build/include -Ivg/gcsa2 -Ivg/cpp -Ivg/sdsl-lite/install/include -Ivg/vcflib/src -Ivg/vcflib -Ivg/vcflib/tabixpp/htslib -Ivg/progress_bar -Ivg/sparsehash/build/include -Ivg/lru_cache -Ivg/fastahack -Ivg/xg -Ivg/xg/sdsl-lite/build/include -Ivg/rocksdb/include
CXXFLAGS=-O3 -std=c++11 -fopenmp -g $(INCLUDES)
#CXXFLAGS=-O0 -std=c++11 -fopenmp -g $(INCLUDES) -DDEBUG
LDSEARCH=-Lvg -Lvg/xg -Lvg/xg/sdsl-lite/build/lib -Lvg/xg/sdsl-lite/build/external/libdivsufsort/lib
LDFLAGS=-lm -lpthread -lz -lbz2 -lsnappy -ldivsufsort -ljansson $(LDSEARCH)
LIBVG=vg/libvg.a
LIBXG=vg/xg/libxg.a
LIBPROTOBUF=vg/protobuf/libprotobuf.a
LIBSDSL=vg/sdsl-lite/install/lib/libsdsl.a
LIBGSSW=vg/gssw/src/libgssw.a
LIBSNAPPY=vg/snappy/libsnappy.a
LIBROCKSDB=vg/rocksdb/librocksdb.a
LIBHTS=vg/htslib/libhts.a
LIBGCSA2=vg/gcsa2/libgcsa2.a
LIBVCFLIB=vg/vcflib/libvcflib.a
VGLIBS=$(LIBVG) $(LIBXG) $(LIBVCFLIB) $(LIBGSSW) $(LIBSNAPPY) $(LIBROCKSDB) $(LIBHTS) $(LIBGCSA2) $(LIBSDSL) $(LIBPROTOBUF)

#Some little adjustments to build on OSX
#(tested with gcc4.9 and jansson installed from MacPorts)
SYS=$(shell uname -s)
ifeq (${SYS},Darwin)
	LDFLAGS:=$(LDFLAGS) -L/opt/local/lib/ # needed for macports jansson
else
	LDFLAGS:=$(LDFLAGS) -lrt
endif

all: snpBridge

$(LIBSDSL): $(LIBVG)

$(LIBPROTOBUF): $(LIBVG)

$(LIBVG):
	cd vg && $(MAKE) libvg.a

$(LIBXG): $(LIBVG)
	cd vg && $(MAKE) xg/libxg.a

# Needs XG to be built for the protobuf headers
main.o: main.cpp snpbridge.h graphvariant.h $(LIBXG)
	$(CXX) main.cpp -c $(CXXFLAGS)

graphvariant.o: graphvariant.h graphvariant.cpp
	$(CXX) graphvariant.cpp -c $(CXXFLAGS)

snpbridge.o: snpbridge.h snpbridge.cpp graphvariant.h
	$(CXX) snpbridge.cpp -c $(CXXFLAGS)

snpBridge: main.o snpbridge.o graphvariant.o $(VGLIBS)
	$(CXX) main.o snpbridge.o graphvariant.o $(VGLIBS) -o snpBridge $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f snpBridge
	rm -f *.o
#	cd vg && $(MAKE) clean
