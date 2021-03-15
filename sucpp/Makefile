H5CXX := h5c++

PLATFORM := $(shell uname -s)
COMPILER := $(shell ($(H5CXX) -v 2>&1) | tr A-Z a-z )

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

EXEFLAGS =

MPFLAG = -fopenmp

LDDFLAGS += $(MPFLAG)
CPPFLAGS += $(MPFLAG)

ifeq ($(PERFORMING_CONDA_BUILD),True)
	CPPFLAGS += -mtune=generic
else
       	CPPFLAGS += -mfma -march=native
endif

CPPFLAGS += -Wextra -Wno-unused-parameter

ifeq ($(PLATFORM),Darwin)
        BLASLIB=-llapacke -lcblas
else
        BLASLIB=-lcblas
endif


LDDFLAGS += -L$(CONDA_PREFIX)/lib
CPPFLAGS += -Wall  -std=c++11 -pedantic -I. $(OPT) -fPIC -L$(CONDA_PREFIX)/lib

ifeq ($(PLATFORM),Darwin)
  LDDFLAGS += -Wl,-rpath,$(PREFIX)/lib
else
  LDDFLAGS += -Wl,-rpath-link,$(PREFIX)/lib
endif
BASE_LDDFLAGS = $(LDDFLAGS)
 
R_LDFLAGS = -llz4 $(BLASLIB)

UFCMP_LIBS=libssu_cpu.so
UFCMP_LINK=-lssu_cpu
ifdef ACC_CXX
	# Tell the generic code we will be building the ACC code, too
	CPPFLAGS += -DUNIFRAC_ENABLE_ACC=1

	UFCMP_LIBS+= libssu_acc.so
	UFCMP_LINK+= -lssu_acc

        R_LDFLAGS += -lssu_acc

	ACC_CPPFLAGS += -mp -acc
	ACC_CPPFLAGS += -Wall  -std=c++11 -pedantic -I. -fPIC -L$(CONDA_PREFIX)/lib

	ifdef DEBUG
                ACC_OPT = -g
	else
                ACC_OPT = -fast
	endif
        ACC_CPPFLAGS += $(ACC_OPT)

        ifeq ($(PERFORMING_CONDA_BUILD),True)
          ACC_CPPFLAGS += -ta=tesla:ccall
        else
          ACC_CPPFLAGS += -ta=tesla
        endif
        # optional info
        ACC_CPPFLAGS += -Minfo=accel

        # use the GNU OMP library to avoid conflicts
        ACC_LDDFLAGS = -shared -mp -acc -Wl,-rpath-link,$(PREFIX)/lib -L$(CONDA_PREFIX)/lib -lgomp -Bstatic_pgi

        ifeq ($(PERFORMING_CONDA_BUILD),True)
                ACC_CPPFLAGS += -tp=px
        endif
endif

ifeq ($(PLATFORM),Darwin)
  TEST_DEPS = -lssu
else
  TEST_DEPS = -lssu -lssu_internal
endif

test_su: test_su.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) test_su.cpp -o test_su $(TEST_DEPS) -lpthread

test_ska: test_ska.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) test_ska.cpp -o test_ska $(TEST_DEPS) -lpthread

test_api: test_api.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) test_api.cpp -o test_api $(TEST_DEPS) -lpthread

test: test_su test_ska test_api
	# test = (su,ska,api)

ssu: su.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) su.cpp -o ssu -lssu -lpthread
	cp ssu ${PREFIX}/bin/

faithpd: faithpd.cpp libssu.so
	$(H5CXX) $(CPPFLAGS) $(EXEFLAGS) faithpd.cpp -o faithpd -lssu -lpthread
	cp faithpd ${PREFIX}/bin/

main: ssu faithpd
	# main == (ssu,faithpd)

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
	echo LDFLAGS=$(R_LDFLAGS) >> ~/.R/Makevars
	Rscript R_interface/rapi_test.R
	
ifeq ($(PLATFORM),Darwin)

# We never use ACC under MacOS, so keep it simple

libssu.so: tree.o biom.o unifrac.o unifrac_internal.o unifrac_cmp_cpu.o cmd.o skbio_alt.o api.o 
	$(H5CXX) $(LDDFLAGS) -o libssu.so tree.o biom.o unifrac.o unifrac_internal.o unifrac_cmp_cpu.o cmd.o skbio_alt.o api.o -lc -lhdf5_cpp -llz4 $(BLASLIB)
	cp libssu.so ${PREFIX}/lib/

else

libssu.so: biom.o unifrac.o cmd.o api.o $(UFCMP_LIBS) libssu_internal.so
	$(H5CXX) $(LDDFLAGS) -o libssu.so biom.o unifrac.o cmd.o api.o $(UFCMP_LINK) -lssu_internal -lc -llz4 
	cp libssu.so ${PREFIX}/lib/

libssu_internal.so: tree.o skbio_alt.o unifrac_internal.o
	$(H5CXX) $(LDDFLAGS) -o libssu_internal.so tree.o skbio_alt.o unifrac_internal.o -lc -lhdf5_cpp $(BLASLIB)
	cp libssu_internal.so ${PREFIX}/lib/

libssu_cpu.so: unifrac_cmp_cpu.o libssu_internal.so
	$(CXX) $(BASE_LDDFLAGS) -o libssu_cpu.so unifrac_cmp_cpu.o -lssu_internal -lc
	cp libssu_cpu.so ${PREFIX}/lib/

libssu_acc.so: unifrac_cmp_acc.o libssu_internal.so
	$(ACC_CXX) $(ACC_LDDFLAGS) -o libssu_acc.so unifrac_cmp_acc.o -lssu_internal -lc
	cp libssu_acc.so ${PREFIX}/lib/

endif

api: libssu.so
	# api == libssu.so

capi_test: api
	gcc -std=c99 capi_test.c -lssu -L${PREFIX}/lib -Wl,-rpath,${PREFIX}/lib -o capi_test
	export LD_LIBRARY_PATH="${PREFIX}/lib":"./capi_test"

api.o: api.cpp api.hpp unifrac.hpp skbio_alt.hpp biom.hpp tree.hpp
	$(H5CXX) $(CPPFLAGS) api.cpp -c -o api.o -fPIC

unifrac.o: unifrac.cpp unifrac.hpp unifrac_internal.hpp unifrac_cmp.hpp biom_interface.hpp tree.hpp
	$(CXX) $(CPPFLAGS) -c $< -o $@

unifrac_cmp_cpu.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.cpp unifrac_task.hpp biom_interface.hpp tree.hpp
	$(CXX) $(CPPFLAGS) -Wno-unknown-pragmas -c $< -o $@

unifrac_cmp_acc.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.cpp unifrac_task.hpp biom_interface.hpp tree.hpp
	$(ACC_CXX) $(ACC_CPPFLAGS) -c $< -o $@


%.o: %.cpp %.hpp
	$(H5CXX) $(CPPFLAGS) -c $< -o $@

clean:
	-rm -f *.o *.so ssu faithpd test_su test_api test_ska

