H5CXX := h5c++

PLATFORM := $(shell uname -s)
COMPILER := $(shell ($(H5CXX) -v 2>&1) | tr A-Z a-z )

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

ifndef NOGPU
	ifneq (,$(findstring pgi,$(COMPILER)))
		CPPFLAGS += -acc
	        ifeq ($(PERFORMING_CONDA_BUILD),True)
	            CPPFLAGS += -ta=tesla:ccall
		else
	            CPPFLAGS += -ta=tesla
                endif
		# optional info
		CPPFLAGS += -Minfo=accel
	        LDDFLAGS += -shlib -acc -Bstatic_pgi
	        EXEFLAGS += -Bstatic_pgi
	endif	
	ifdef SMALLGPU
		CPPFLAGS += -DSMALLGPU
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

test: test_su.cpp test_ska.cpp test_api.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) test_su.cpp -o test_su -lssu -lpthread
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) test_ska.cpp -o test_ska -lssu -lpthread
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) test_api.cpp -o test_api -lssu -lpthread

main: su.cpp faithpd.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) su.cpp -o ssu -lssu -lpthread
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) faithpd.cpp -o faithpd -lssu -lpthread
	cp ssu ${PREFIX}/bin/
	cp faithpd ${PREFIX}/bin/

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
	echo LDFLAGS=-llz4 $(BLASLIB) >> ~/.R/Makevars
	Rscript R_interface/rapi_test.R
	
libssu.so: tree.o biom.o unifrac_cmp.o unifrac.o cmd.o skbio_alt.o api.o
	$(H5CXX) $(LDDFLAGS) -o libssu.so tree.o biom.o unifrac_cmp.o unifrac.o cmd.o skbio_alt.o api.o -lc -lhdf5_cpp -llz4 $(BLASLIB) -L$(PREFIX)/lib
	cp libssu.so ${PREFIX}/lib/

api: libssu.so

capi_test: api
	gcc -std=c99 capi_test.c -lssu -L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib -o capi_test
	export LD_LIBRARY_PATH="${PREFIX}/lib":"./capi_test"

api.o: api.cpp api.hpp unifrac.hpp skbio_alt.hpp biom.hpp tree.hpp
	$(H5CXX) $(CPPFLAGS) api.cpp -c -o api.o -fPIC

unifrac.o: unifrac.cpp unifrac.hpp unifrac_internal.hpp unifrac_cmp.hpp biom_interface.hpp tree.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

unifrac_cmp.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.cpp unifrac_task.hpp biom_interface.hpp tree.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@


%.o: %.cpp %.hpp
	$(H5CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f *.o ssu

