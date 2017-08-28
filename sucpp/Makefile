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

CPPFLAGS = -Wall -Wextra -std=c++11 -pedantic -I. $(OPT) -mfma -march=native -fPIC

ifeq ($(PLATFORM),Darwin)
	AVX2 := $(shell sysctl -a | grep -c AVX2)
else
	AVX2 := $(shell grep "^flags" /proc/cpuinfo | head -n 1 | grep -c avx2)
endif

ifeq ($(AVX2),1)
	CPPFLAGS += -mavx2
else
	CPPFLAGS += -mavx
endif


test: tree.o test_su.cpp biom.o unifrac.o unifrac_task.o api.o
	$(CXX) $(CPPFLAGS) -Wno-unused-parameter test_su.cpp -o test_su tree.o biom.o unifrac.o unifrac_task.o api.o

main: tree.o biom.o unifrac.o cmd.o unifrac_task.o api.o
	$(CXX) $(CPPFLAGS) su.cpp -o ssu tree.o biom.o unifrac.o cmd.o unifrac_task.o api.o

api: tree.o biom.o unifrac.o cmd.o unifrac_task.o
	$(CXX) $(CPPFLAGS) api.cpp -c -o api.o -fPIC
	$(CXX) -shared -fPIC -Wl -o libssu.so tree.o biom.o unifrac.o cmd.o unifrac_task.o api.o -lc

%.o: %.cpp %.hpp
	$(CXX) $(CPPFLAGS) -static -c $< -o $@

clean:
	-rm *.o ssu
