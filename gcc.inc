# ------------------------------------------------
# ----- Setting for GCC environment
# ----- Copyright (c) 2019 Tomonori Kouya
# ------------------------------------------------
#CPP = g++
CPP = g++ -std=c++11
#CC = gcc -O1
CC = gcc -O3

INC = -I/usr/local/include
#LIB = -L/usr/local/lib -lm
LIB = -L/usr/local/lib -lm -lstdc++

CBLAS_INC = $(INC)
#CBLAS_LIB = $(LIB) -lcblas -lrefblas -lm
CBLAS_LIB = $(LIB) -lblas -lm

#LAPACKE_INC = -I/usr/local/include/lapacke $(CBLAS_INC)
LAPACKE_INC = $(CBLAS_INC)
LAPACKE_LIB = -L/usr/local/lib -llapacke -llapack $(CBLAS_LIB)

OPENMP = -fopenmp
OPENMP_LIB = -lgomp

all: quadratic_eq get_sec matmul cg iterative_ref mpfr qd lu gcc_float128 gmp

# GCC Quadruple precision library (__float128, __complex128)
gcc_float128: float128.c complex128.c
	$(CC) test_quadmath.c -lquadmath -o test_quadmath
	$(CC) float128.c -o float128 -lquadmath
	$(CC) complex128.c -o complex128 -lquadmath

