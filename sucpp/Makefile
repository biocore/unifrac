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

rapi_test:
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
	
capi_test: api
	gcc -std=c99 capi_test.c -lssu -L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib -o capi_test
	export LD_LIBRARY_PATH="${PREFIX}/lib":"./capi_test"


%.o: %.cpp %.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f *.o ssu

