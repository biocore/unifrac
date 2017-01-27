CXX := h5c++

PLATFORM := $(shell uname)
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

CPPFLAGS = -Wall -Wextra -std=c++11 -pedantic -I. $(OPT)

test: tree.o test_su.cpp biom.o unifrac.o
	$(CXX) $(CPPFLAGS) -Wno-unused-parameter test_su.cpp -o test_su tree.o biom.o unifrac.o 
	./test_su

main: tree.o biom.o unifrac.o cmd.o
	$(CXX) $(CPPFLAGS) su.cpp -o su tree.o biom.o unifrac.o cmd.o
	cp su ../

ocltest: opencl.o
	clang++ -framework OpenCL opencl.cpp -c
	clang++ -framework OpenCL opencltest.cpp opencl.o -o opencltest

%.o: %.cpp %.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	rm *.o
