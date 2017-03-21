CXX := h5c++

PLATFORM := $(shell uname -s)
COMPILER := $(shell ($(CXX) -v 2>&1) | tr A-Z a-z )

ifdef DEBUG
	OPT = -O0 -DDEBUG=1 --debug -g -ggdb
else
	ifneq (,$(findstring gcc,$(COMPILER)))
		OPT = -O4
		TGTFLAGS = -fwhole-program
	else
		OPT = -O3
	endif
endif

CPPFLAGS = -Wall -pthreads -Wextra -std=c++11 -pedantic -I. $(OPT) -mfma 

ifeq ($(PLATFORM),Darwin)
	AVX2 := $(shell sysctl -a | grep -c AVX2)
else
	AVX2 := $(shell grep "^flags" | head -n 1 | grep -c avx2)
endif

ifeq ($(AVX2),1)
	CPPFLAGS += -mavx2
else
	CPPFLAGS += -mavx
endif


test: tree.o test_su.cpp biom.o unifrac.o
	$(CXX) $(CPPFLAGS) -Wno-unused-parameter test_su.cpp -o test_su tree.o biom.o unifrac.o 

main: tree.o biom.o unifrac.o cmd.o
	$(CXX) $(CPPFLAGS) su.cpp -o ssu tree.o biom.o unifrac.o cmd.o
	cp ssu ../

ocltest: opencl.o
	clang++ -framework OpenCL opencl.cpp -c
	clang++ -framework OpenCL opencltest.cpp opencl.o -o opencltest

%.o: %.cpp %.hpp
	$(CXX) $(CPPFLAGS) -static -c $< -o $@

clean:
	rm *.o
