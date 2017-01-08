SHELL=/bin/bash

# Compiler, I only test this using gcc/5.3
CC=gcc-5

# Compiler optimization level. (O3 is recommended)
#OPT_LEVEL=-O3 #-qopt-report

# OpenMP flag, -fopenmp for gcc, -openmp for icc. Leave it blank if use clang
OMPFLAG=-fopenmp

# Header and lib path for openBLAS, change them to ur install path
HPATH=/Users/Linling/Documents/OpenBLAS/include/
LIBPATH=/Users/Linling/Documents/OpenBLAS/lib


################
# OTHER INPUTS #
################
# Lenard-Jones parameter
LJ_EPSILON=0.31
# Spring Constant
SPR_CONST=200.0
# Delta T
DELTA_T=0.001
# Total time steps
TS_TOTAL=10000
# Ewald sum parameter
K_MAX=4
# Overlap concentration
C_OVERLAP=0.05102
# Bin size
BINSIZE=2.0
# Max iterations you want to run
ITER_MAX=10
################

# Common program arguments.
COMMON_PROG_ARGS= \
								 -std=c99 \
								 $(OMPFLAG) \
								 $(OPT_LEVEL) \
								 #-lrt \

# Input arguments.
INPUT_ARGS = \
		-DEPSILON=${LJ_EPSILON} \
		-DKAPPA=${SPR_CONST} \
		-DDT=${DELTA_T} \
		-DTMAX=${TS_TOTAL} \
		-DKMAX=${K_MAX} \
		-DC_STAR=${C_OVERLAP} \
		-DBIN_SIZE=${BINSIZE} \
		-DMAXITER=${ITER_MAX} \
								
# Program arguments.
OMP_CC = $(CC) $(COMMON_PROG_ARGS) $(INPUT_ARGS)

SRC=*.c
TARGETS=smdlt_static

all: $(TARGETS)

smdlt_static: $(SRC)
	$(OMP_CC) $(SRC) -lopenblas -I$(HPATH) -L$(LIBPATH) -o $@
	if [ ! -d ./output ];then \
	mkdir ./output;	\
	fi

clean:
	rm -f *.o *.optrpt $(TARGETS)
