CXX := h5c++

PLATFORM := $(shell uname -s)
COMPILER := $(shell ($(CXX) -v 2>&1) | tr A-Z a-z )

ifdef DEBUG
	ifneq (,$(findstring pgi,$(COMPILER)))
		OPT = -g
	else
		OPT = -O0 -DDEBUG=1 --debug -g -ggdb
	endif
else
	ifneq (,$(findstring pgi,$(COMPILER)))
		OPT = -fast
	else
	  ifneq (,$(findstring gcc,$(COMPILER)))
		OPT = -O4
		TGTFLAGS = -fwhole-program
	  else
		OPT = -O3
	  endif
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

EXEFLAGS =

ifneq (,$(findstring pgi,$(COMPILER)))
	MPFLAG =  -mp
else
	MPFLAG = -fopenmp
endif

LDDFLAGS += $(MPFLAG)
CPPFLAGS += $(MPFLAG)

UNIFRAC_FILES = unifrac_internal.o unifrac_cmp_cpu.o

ifndef NOGPU
	ifneq (,$(findstring pgi,$(COMPILER)))
		CPPFLAGS += -DUNIFRAC_ENABLE_ACC=1
                UNIFRAC_FILES += unifrac_cmp_acc.o
		ACCCPPFLAGS += -acc
	        ifeq ($(PERFORMING_CONDA_BUILD),True)
	            ACCCPPFLAGS += -ta=tesla:ccall
		else
	            ACCCPPFLAGS += -ta=tesla
                endif
		# optional info
		ACCCPPFLAGS += -Minfo=accel
	        LDDFLAGS += -shlib -acc -Bstatic_pgi
	        EXEFLAGS += -acc -Bstatic_pgi
	endif
endif

ifneq (,$(findstring pgi,$(COMPILER)))
	ifeq ($(PERFORMING_CONDA_BUILD),True)
		CPPFLAGS += -tp=px
	endif
else
	ifeq ($(PERFORMING_CONDA_BUILD),True)
		CPPFLAGS += -mtune=generic
	else
		CPPFLAGS += -mfma -march=native
	endif
endif

ifeq (,$(findstring pgi,$(COMPILER)))
	# basically, not gcc
	CPPFLAGS += -Wextra -Wno-unused-parameter
endif

ifeq ($(PLATFORM),Darwin)
        BLASLIB=-llapacke -lcblas
else
        BLASLIB=-lcblas
endif


CPPFLAGS += -Wall  -std=c++11 -pedantic -I. $(OPT) -fPIC -L$(CONDA_PREFIX)/lib

test: tree.o test_su.cpp biom.o unifrac.o skbio_alt.o api.o $(UNIFRAC_FILES)
	$(CXX) $(CPPFLAGS) $(EXEFLAGS) test_su.cpp -o test_su tree.o biom.o $(UNIFRAC_FILES) unifrac.o skbio_alt.o api.o -llz4 $(BLASLIB) -lpthread
	$(CXX) $(CPPFLAGS) $(EXEFLAGS) test_ska.cpp -o test_ska skbio_alt.o tree.o biom.o $(UNIFRAC_FILES) unifrac.o api.o -llz4 $(BLASLIB) -lpthread
	$(CXX) $(CPPFLAGS) $(EXEFLAGS) test_api.cpp -o test_api tree.o biom.o $(UNIFRAC_FILES) unifrac.o skbio_alt.o api.o -llz4 $(BLASLIB) -lpthread

main: tree.o biom.o unifrac.o cmd.o skbio_alt.o api.o $(UNIFRAC_FILES)
	$(CXX) $(CPPFLAGS) $(EXEFLAGS) su.cpp -o ssu tree.o biom.o $(UNIFRAC_FILES) unifrac.o cmd.o skbio_alt.o api.o -lhdf5_cpp -llz4 $(BLASLIB) -lpthread
	$(CXX) $(CPPFLAGS) $(EXEFLAGS) faithpd.cpp -o faithpd tree.o biom.o $(UNIFRAC_FILES) unifrac.o cmd.o skbio_alt.o api.o -lhdf5_cpp -llz4 $(BLASLIB) -lpthread
	cp ssu ${PREFIX}/bin/
	cp faithpd ${PREFIX}/bin/

rapi_test: main
	mkdir -p ~/.R
	if [ -e ~/.R/Makevars ] ; \
	then \
		echo "WARNING: OVERWRITING ~/.R/Makevars" ; \
		echo "The original Makevars file has been copied to ~/.R/Makevars" ;\
		cp ~/.R/Makevars Makevars-original ; \
	fi;
	echo CXX1X=h5c++ > ~/.R/Makevars
	echo CXX=h5c++ >> ~/.R/Makevars 
	echo CC=h5c++ >> ~/.R/Makevars
	echo LDFLAGS= -llz4 $(BLASLIB) >> ~/.R/Makevars
	rm -f *.o
	Rscript R_interface/rapi_test.R
	rm -f *.o
	
api: tree.o biom.o unifrac.o cmd.o skbio_alt.o api.o $(UNIFRAC_FILES)
	$(CXX) $(LDDFLAGS) -o libssu.so tree.o biom.o $(UNIFRAC_FILES) unifrac.o cmd.o skbio_alt.o api.o -lc -lhdf5_cpp -llz4 $(BLASLIB) -L$(PREFIX)/lib
	cp libssu.so ${PREFIX}/lib/

capi_test: api
	gcc -std=c99 capi_test.c -lssu -L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib -o capi_test
	export LD_LIBRARY_PATH="${PREFIX}/lib":"./capi_test"

api.o: api.cpp api.hpp unifrac.hpp skbio_alt.hpp biom.hpp tree.hpp
	$(CXX) $(CPPFLAGS) api.cpp -c -o api.o -fPIC

unifrac_cmp_cpu.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.cpp unifrac_task.hpp biom_interface.hpp tree.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

unifrac_cmp_acc.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.cpp unifrac_task.hpp biom_interface.hpp tree.hpp
	$(CXX) $(CPPFLAGS) $(ACCCPPFLAGS) -c $< -o $@


%.o: %.cpp %.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f *.o ssu

