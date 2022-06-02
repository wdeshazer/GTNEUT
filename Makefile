# Multi-platform Makefile for the GTNEUT neutral transport code. 
# This Makefile has been written assuming use of the GNU Make. It may
# or may not work with other implementations of the Make utility.
# Written by John Mandrekas, GIT, 08/15/96, to replace the old Makefile.
# Currently supports CRAY, SUN and HP systems.

# June 10, 2003, Sparse matrix version (UMFPACK)

# Set the desired system:

# HP   : General Atomics HP
# SUN  : SUN SparcStation 2 running SUNOS
# SOL  : SUN ULTRA 10 running Solaris 7
# CRAY : NERSC CRAY
# WSL2 : Windows Subsystem Layer (Linux VM)

SYS = WSL2

ifeq ($(SYS), WSL2)
	F = .f
	O = .o
	E = 
	FF = gfortran
	LD = gfortran
	FFLAGS = -O0 -c -g -ggdb
	LDFLAGS = -O0 -LUMFPACK2 -Ljson-fortran/lib
	LIBS =  -lumfpack -llapack -lblas -ljsonfortran
endif

SOURCES= main.f		\
     calctransm.f	\
	 transmcoeff.f	\
	 rectinp.f      \
	 checkinp.f 	\
	 calcrparms.f	\
	 calcRect.f 	\
	 tij.f		\
	 qgauss.f	\
	 calcmfp.f	\
	 svjanev.f	\
	 degasread.f    \
	 svdegas.f	\
	 calcrefln.f	\
	 reflect.f	\
	 escape.f	\
	 setup.f	\
	 solvers.f	\
	 output.f	\
	 calcxswms.f	\
	 wmsdata.f	\
	 simpson.f	\
	 bickley.f	\
	 ndata.f	\
	 zstop.f        \
	 fem.f          \
	 pbalance.f

OBJ = $(SOURCES:$F=$O)

gtneut$E : $(OBJ) UMFPACK2/libumfpack.a json-fortran/lib/libjsonfortran.a
	echo Makefile for GTNEUT     20081010 tbt
	echo SYS = $(SYS)
	echo 
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

%$O : %$F
	$(FF) $(FFLAGS) $<

installumfpack:
	$(MAKE) -C UMFPACK2 install

# Include file dependencies:

UMFPACK2/libumfpack.a:
	$(MAKE) -C UMFPACK2

json-fortran/lib/libjsonfortran.a:
	pip install Fobis.py
	pip install ford
	cd json-fortran; ./build.sh

main.o : \
	consts.inc \
	neutGlob.inc \
	comiou.inc 

rectinp.o :	\
	neutGlob.inc \
	comiou.inc

calcRect.o : 	\
	locGeom.inc \
    neutGlob.inc

calcmfp.o : \
	consts.inc	\
	locGeom.inc	\
	neutGlob.inc

calcrparms.o : \
	consts.inc	\
	locGeom.inc	\
	neutGlob.inc

calctransm.o : \
	neutGlob.inc

transmcoeff.o : \
	consts.inc	\
	locGeom.inc	\
	neutGlob.inc

calcxswms.o : \
	wmsdata.inc

checkinp.o : \
	consts.inc \
	neutGlob.inc \
	comiou.inc

degasread.o : degdat.inc \
		comiou.inc

svdegas.o   : degdat.inc

calcrefln.o : neutGlob.inc

escape.o : \
	consts.inc \
	locGeom.inc \
	neutGlob.inc \

output.o : \
	consts.inc \
	neutGlob.inc \
	comiou.inc 

setup.o : \
	consts.inc \
	neutGlob.inc \


solvers.o : \
	consts.inc \
	neutGlob.inc

tij.o : \
	locGeom.inc \

wmsdata.o : \
	wmsdata.inc \

ndata.o  :	\
	comiou.inc \
	consts.inc

pbalance.o: \
	neutGlob.inc

fem.o: \
	 neutGlob.inc \
	 esc.inc

clean:
	rm -f gtneut$E $(OBJ)

realclean: \
	clean
	$(MAKE) -C UMFPACK2 purge
	cd json-fortran; ./build.sh --clean

