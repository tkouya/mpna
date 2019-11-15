# ------------------------------------------------
# ----- Setting for Intel Compiler environment
# ----- Copyright (c) 2019 Tomonori Kouya
# ------------------------------------------------
CPP = icpc -std=c++11 -fp-model precise
CC = icc -fp-model precise

INC = -I/opt/intel/include -I/usr/local/include
LIB = -L/opt/intel/lib/intel64 -L/usr/local/lib

#IMKL_INC=-I/opt/intel/include -I/opt/intel/mkl/include 
IMKL_INC=-DUSE_IMKL -I/opt/intel/include -I/opt/intel/mkl/include 
#IMKL_LIB = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -L/opt/intel/lib/intel64 -lifcore
IMKL_LIB = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -L/opt/intel/lib/intel64 -lifcore -liomp5

CBLAS_INC = $(IMKL_INC)
#CBLAS_LIB = $(LIB) -lcblas -lrefblas -lm
CBLAS_LIB = $(IMKL_LIB)

LAPACKE_INC = $(CBLAS_INC)
LAPACKE_LIB = $(CBLAS_LIB)

OPENMP = -openmp
OPENMP_LIB =

all: quadratic_eq get_sec num_size matmul cg iterative_ref gmp rsa
