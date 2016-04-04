# this Makefile is largely derived from https://github.com/adamnovak/corg/blob/master/Makefile
.PHONY: all clean

CXX=g++
INCLUDES=-Ivg/src -Ivg/include
CXXFLAGS=-O3 -std=c++11 -fopenmp -g $(INCLUDES)
VGLIBDIR=vg/lib
LDSEARCH=-Lvg/src -L$(VGLIBDIR)
LDFLAGS=-lm -lpthread -lz -lbz2 -lsnappy -ldivsufsort -ldivsufsort64 -ljansson $(LDSEARCH)
LIBVG=$(VGLIBDIR)/libvg.a
LIBXG=$(VGLIBDIR)/libxg.a
LIBPROTOBUF=$(VGLIBDIR)/libprotobuf.a
LIBSDSL=$(VGLIBDIR)/libsdsl.a
LIBGSSW=$(VGLIBDIR)/libgssw.a
LIBSNAPPY=$(VGLIBDIR)/libsnappy.a
LIBROCKSDB=$(VGLIBDIR)/librocksdb.a
LIBHTS=$(VGLIBDIR)/libhts.a
LIBGCSA2=$(VGLIBDIR)/libgcsa2.a
LIBVCFLIB=$(VGLIBDIR)/libvcflib.a
LIBRAPTOR=$(VGLIBDIR)/libraptor2.a
LIBGFAK=$(VGLIBDIR)/libgfakluge.a
VGLIBS=$(LIBVG) $(LIBXG) $(LIBVCFLIB) $(LIBGSSW) $(LIBSNAPPY) $(LIBROCKSDB) $(LIBHTS) $(LIBGCSA2) $(LIBSDSL) $(LIBPROTOBUF) $(LIBRAPTOR) $(LIBGFAK)

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
	cd vg && . ./source_me.sh && $(MAKE) lib/libvg.a

$(LIBXG): $(LIBVG)
	cd vg && . ./source_me.sh && $(MAKE) lib/libxg.a

# Needs XG to be built for the protobuf headers
main.o: $(LIBXG)

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
