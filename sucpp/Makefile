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

ifeq ($(PREFIX),)
	PREFIX := $(CONDA_PREFIX)
endif

ifeq ($(PLATFORM),Darwin)
	AVX2 := $(shell sysctl -a | grep -c AVX2)
	LDDFLAGS = -dynamiclib -install_name @rpath/libssu.so
else
	AVX2 := $(shell grep "^flags" /proc/cpuinfo | head -n 1 | grep -c avx2)
	LDDFLAGS = -shared
endif

ifeq ($(PERFORMING_CONDA_BUILD),True)
	CPPFLAGS += -mtune=generic -fPIC
	CPPFLAGS=$(echo $CPPFLAGS | sed 's/-fno-plt //')
	LDDFLAGS += -fPIC
else
	CPPFLAGS += -mfma -march=native
endif

CPPFLAGS += -Wall -Wextra -std=c++11 -pedantic -I. $(OPT) -fPIC

test: tree.o test_su.cpp biom.o unifrac.o unifrac_task.o api.o
	$(CXX) $(CPPFLAGS) -Wno-unused-parameter test_su.cpp -o test_su tree.o biom.o unifrac.o unifrac_task.o api.o -pthread
	$(CXX) $(CPPFLAGS) -Wno-unused-parameter test_api.cpp -o test_api tree.o biom.o unifrac.o unifrac_task.o api.o -pthread

main: tree.o biom.o unifrac.o cmd.o unifrac_task.o api.o
	$(CXX) $(CPPFLAGS) su.cpp -o ssu tree.o biom.o unifrac.o cmd.o unifrac_task.o api.o -lhdf5_cpp -pthread
	cp ssu ${PREFIX}/bin/

rapi_test: main
	mkdir -p ~/.R
	if [ -a ~/.R/Makevars ] ; \
	then \
		echo "WARNING: OVERWRITING ~/.R/Makevars" ; \
		echo "The original Makevars file has been copied to ~/.R/Makevars" ;\
		cp ~/.R/Makevars Makevars-original ; \
	fi;
	echo CXX1X=h5c++ > ~/.R/Makevars
	echo CXX=h5c++ >> ~/.R/Makevars 
	echo CC=h5c++ >> ~/.R/Makevars
	pushd R_interface && Rscript rapi_test.R; popd
	
api: tree.o biom.o unifrac.o cmd.o unifrac_task.o
	$(CXX) $(CPPFLAGS) -v api.cpp -c -o api.o -fPIC
	$(CXX) $(LDDFLAGS) -v -o libssu.so tree.o biom.o unifrac.o cmd.o unifrac_task.o api.o -lc -lhdf5_cpp -L$(PREFIX)/lib
	cp libssu.so ${PREFIX}/lib/

capi_test: api
	gcc -std=c99 capi_test.c -lssu -L${PREFIX}/lib -o capi_test
	LD_LIBRARY_PATH="${PREFIX}/lib":"./capi_test"

%.o: %.cpp %.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f *.o ssu

