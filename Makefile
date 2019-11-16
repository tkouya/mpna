# ------------------------------------------------
# ----- Makefile for "Multiple Precision Numerical Computation"
# ----- Copyright (c) 2019 Tomonori Kouya
# ------------------------------------------------
# GCC
include gcc.inc

# Intel C/C++
#include icc.inc

GMP_INC = $(INC)
GMP_LIB =  -lgmp $(LIB)
GMPXX_LIB = -lgmpxx $(GMP_LIB)
INNER_GMP_INC = -I/home/tkouya/pool/gmp/gmp-6.1.2 -I/home/tkouya/pool/gmp/gmp-6.1.2/mpn
#INNER_GMP_LIB = /home/tkouya/local/gmp/noassembly/lib/libgmp.la
INNER_GMP_LIB = /home/tkouya/local/gmp/x64/lib/libgmp.a
MPFR_INC = $(INC)
MPFR_LIB = $(LIB) -lmpfr $(GMP_LIB)
ARB_INC = $(INC)
ARB_LIB = $(LIB) -larb -lflint $(MPFR_LIB)
MPFI_INC = $(INC)
MPFI_LIB = $(LIB) -lmpfi $(MPFR_LIB)
QD_INC = $(INC)/qd
QD_LIB = $(LIB) -lqd
MPC_LIB = $(LIB) -lmpc

rsa: rsa.c
	$(CC) $(GMP_INC) rsa.c -o rsa $(GMP_LIB)

quadratic_eq: quadratic_eq.c
	$(CC) $(INC) quadratic_eq.c -o quadratic_eq $(LIB)

logistic: logistic_f.c logistic.c logistic_dd.cpp logistic_qd.cpp logistic_mpreal.cpp logistic_mpreal_dd.cpp logistic_mpreal_dd_qd.cpp logistic_mpfr.c logistic_err.c logistic_f_err.c
	$(CC) $(INC) logistic_f.c -o logistic_f $(LIB)
	$(CC) $(INC) logistic_f_err.c -o logistic_f_err $(LIB)
	$(CC) $(INC) logistic.c -o logistic $(LIB)
	$(CC) $(INC) logistic_err.c -o logistic_err $(LIB)
	$(CPP) $(MPFR_INC) logistic_mpfr.c -o logistic_mpfr $(MPFR_LIB)
	$(CPP) $(MPFR_INC) logistic_mpreal.cpp -o logistic_mpreal $(MPFR_LIB)
	$(CPP) $(QD_INC) $(INC) logistic_dd.cpp -o logistic_dd $(QD_LIB)
	$(CPP) $(QD_INC) $(INC) logistic_qd.cpp -o logistic_qd $(QD_LIB)
	$(CPP) $(QD_INC) $(INC) $(MPFR_INC) logistic_mpreal_dd.cpp -o logistic_mpreal_dd $(MPFR_LIB) $(QD_LIB)
	$(CPP) $(QD_INC) $(INC) $(MPFR_INC) logistic_mpreal_dd_qd.cpp -o logistic_mpreal_dd_qd $(MPFR_LIB) $(QD_LIB)

get_sec: get_sec.cpp get_secv.h
	$(CPP) -c get_sec.cpp -o get_sec.o

gmp: mpn_sample_full.c  mpz_test.c mpz_input.c mpz_mersenne.c mpf_template.c mpz_factorial.c mpz_input_nextprime.c mpz_prime_factorization.c mpz_binomial.c mpz_input_gcd_lcm.c mpf_relerr.c mpq_input.c mpq_input_convert.c
	$(CC) $(GMP_INC) mpn_sample_full.c -o mpn_sample_full $(GMP_LIB)
	$(CC) $(GMP_INC) -DUSE_STR_INPUT mpz_input.c -o mpz_input_str $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_test.c -o mpz_test $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_input.c -o mpz_input $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_mersenne.c -o mpz_mersenne $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_factorial.c -o mpz_factorial $(GMP_LIB)
	$(CC) $(GMP_INC) mpf_template.c -o mpf_template $(GMP_LIB)
	$(CC) $(GMP_INC) mpf_relerr.c -o mpf_relerr $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_input_nextprime.c -o mpz_input_nextprime $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_prime_factorization.c -o mpz_prime_factorization $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_mersenne.c -o mpz_mersenne $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_binomial.c -o mpz_binomial $(GMP_LIB)
	$(CC) $(GMP_INC) mpz_input_gcd_lcm.c -o mpz_input_gcd_lcm $(GMP_LIB)
	$(CC) $(GMP_INC) mpq_input.c -o mpq_input $(GMP_LIB)
	$(CC) $(GMP_INC) mpq_input_convert.c -o mpq_input_convert $(GMP_LIB)

gmp_class: mpf_template.cpp mpz_test.cpp mpz_mersenne.cpp mpz_input_nextprime.cpp mpz_prime_factorization.cpp mpz_input_gcd_lcm.cpp mpq_test.cpp
	$(CPP) $(GMP_INC) mpf_template.cpp -o mpf_template_cxx $(GMPXX_LIB)
	$(CPP) $(GMP_INC) mpz_mersenne.cpp -o mpz_mersenne_cxx $(GMPXX_LIB)
	$(CPP) $(GMP_INC) mpz_input_nextprime.cpp -o mpz_input_nextprime_cxx $(GMPXX_LIB)
	$(CPP) $(GMP_INC) mpz_prime_factorization.cpp -o mpz_prime_factorization_cxx $(GMPXX_LIB)
	$(CPP) $(GMP_INC) mpz_input_gcd_lcm.cpp -o mpz_input_gcd_lcm_cxx $(GMPXX_LIB)
	$(CPP) $(GMP_INC) mpq_test.cpp -o mpq_test_cxx $(GMPXX_LIB)
	$(CPP) $(GMP_INC) mpz_test.cpp -o mpz_test_cxx $(GMPXX_LIB)

mpfr: mpfr_template.c mpfr_template.cpp mpfr_newton_inverse.c mpfr_newton_sqrt.c mpfr_pi_simple.c mpfr_exp.c  mpfr_printf.c mpfr_dec2bin.c mpfr_relerr.c mpfr_relerr.cpp
	$(CC) $(MPFR_INC) mpfr_template.c -o mpfr_template $(MPFR_LIB)
	$(CPP) $(MPFR_INC) mpfr_template.cpp -o mpfr_template_cxx $(MPFR_LIB)
	$(CC) $(MPFR_INC) mpfr_newton_inverse.c -o mpfr_newton_inverse $(MPFR_LIB)
	$(CC) $(MPFR_INC) mpfr_newton_sqrt.c -o mpfr_newton_sqrt $(MPFR_LIB)
	$(CC) $(MPFR_INC) mpfr_pi_simple.c -o mpfr_pi_simple $(MPFR_LIB)
	$(CC) $(MPFR_INC) mpfr_exp.c -o mpfr_exp $(MPFR_LIB)
	$(CC) $(MPFR_INC) mpfr_printf.c -o mpfr_printf $(MPFR_LIB)
	$(CC) $(MPFR_INC) mpfr_relerr.c -o mpfr_relerr $(MPFR_LIB)
	$(CPP) $(MPFR_INC) mpfr_relerr.cpp -o mpfr_relerr_cxx $(MPFR_LIB)

matmul: matmul_simple.cpp matmul_simple_omp.cpp matmul_block.cpp matmul_winograd.cpp matmul_block.h get_sec
	$(CPP) matmul_simple.cpp get_sec.o -o matmul_simple
	$(CPP) -DBLOCK matmul_block.cpp get_sec.o -o matmul_block
	$(CPP) -DSTRASSEN matmul_block.cpp get_sec.o -o matmul_strassen
	$(CPP) matmul_winograd.cpp get_sec.o -o matmul_winograd
	$(CPP) $(OPENMP) matmul_simple_omp.cpp get_sec.o -o matmul_simple_omp $(OPENMP_LIB)
	$(CPP) $(OPENMP) -DBLOCK matmul_block.cpp get_sec.o -o matmul_block_omp $(OPENMP_LIB)

arb: test_arb.c quadratic_eq_arb.c logistic_arb.c
	$(CC) $(ARB_INC) test_arb.c -o test_arb $(ARB_LIB)
	$(CC) $(ARB_INC) quadratic_eq_arb.c -o quadratic_eq_arb $(ARB_LIB)
	$(CC) $(ARB_INC) logistic_arb.c -o logistic_arb $(ARB_LIB)

mpfi: test_mpfi.c quadratic_eq_mpfi.c
	$(CC) $(MPFI_INC) test_mpfi.c -o test_mpfi $(MPFI_LIB)
	$(CC) $(MPFI_INC) quadratic_eq_mpfi.c -o quadratic_eq_mpfi $(MPFI_LIB)
	$(CC) $(MPFI_INC) logistic_mpfi.c -o logistic_mpfi $(MPFI_LIB)

linear: linear_c.h linear_c_omp.h template_linear.h sample_linear.c sample_linear.cpp
	$(CC) $(MPFR_INC) sample_linear.c -o sample_linear $(MPFR_LIB)
	$(CPP) $(MPFR_INC) sample_linear.cpp -o sample_linear_cxx $(MPFR_LIB)

complex: template_complex_efunc.cpp complex_mpc.c complex_dd.cpp complex_qd.cpp complex_mpreal.cpp complex_d.cpp complex_d.c test_mpc.c
	$(CPP) $(MPFR_INC) template_complex_efunc.cpp -o template_complex_efunc $(MPFR_LIB) $(QD_LIB)
	$(CC) $(MPFR_INC) test_mpc.c -o test_mpc $(MPC_LIB) $(MPFR_LIB)
	$(CPP) $(QD_INC) complex_dd.cpp -o complex_dd $(QD_LIB)
	$(CPP) $(QD_INC) complex_qd.cpp -o complex_qd $(QD_LIB)
	$(CPP) $(MPFR_INC) complex_mpreal.cpp -o complex_mpreal $(MPFR_LIB)
	$(CPP) complex_d.cpp -o complex_d_cpp $(LIB)
	$(CC) complex_d.c -o complex_d_c $(LIB)

qd: dd_test.cpp dd_print.cpp qd_test.cpp
	$(CPP) $(QD_INC) dd_test.cpp -o dd_test $(QD_LIB)
	$(CPP) $(QD_INC) dd_print.cpp -o dd_print $(QD_LIB)
	$(CPP) $(QD_INC) qd_test.cpp -o qd_test $(QD_LIB)

cg: cg_d_omp.c cg_mpfr_omp.c get_sec cg_d.cpp cg_dd.cpp cg_qd.cpp cg_mpreal.cpp
	$(CC) cg_d_omp.c get_sec.o -o cg_double $(LIB)
	$(CC) $(OPENMP) cg_d_omp.c get_sec.o -o cg_d_omp $(LIB) $(OPENMP_LIB)
	$(CC) $(MPFR_INC) cg_mpfr_omp.c -o cg_mpfr get_sec.cpp $(MPFR_LIB)
	$(CC) $(MPFR_INC) $(OPENMP) cg_mpfr_omp.c get_sec.o -o cg_mpfr_omp $(MPFR_LIB) $(OPENMP_LIB)
	$(CPP) $(INC) cg_d.cpp get_sec.o -o cg_d_cxx $(LIB)
	$(CPP) $(QD_INC) cg_dd.cpp get_sec.o -o cg_dd_cxx $(QD_LIB)
	$(CPP) $(QD_INC) cg_qd.cpp get_sec.o -o cg_qd_cxx $(QD_LIB)
	$(CPP) $(MPFR_INC) cg_mpreal.cpp get_sec.o -o cg_mpreal_cxx $(MPFR_LIB)

lu: lu_d.c template_lu_dd.cpp template_lu.cpp template_lu_qd.cpp template_lu_mpreal.cpp
	$(CC) $(INC) lu_d.c get_sec.o -o lu_d $(LIB)
	$(CC) $(INC) template_lu.cpp get_sec.o -o lu_d_cxx $(LIB)
	$(CPP) $(QD_INC) template_lu_dd.cpp get_sec.o -o lu_dd $(QD_LIB)
	$(CPP) $(QD_INC) template_lu_qd.cpp get_sec.o -o lu_qd $(QD_LIB)
	$(CPP) $(MPFR_INC) template_lu_mpreal.cpp get_sec.o -o lu_mpreal $(MPFR_LIB)

power: power_mpreal_hilbert.cpp cond_mpreal_hilbert.cpp
	$(CPP) $(MPFR_INC) power_mpreal_hilbert.cpp get_sec.o -o power_mpreal_hilbert $(MPFR_LIB)
	$(CPP) $(MPFR_INC) inverse_power_mpreal_hilbert.cpp get_sec.o -o inverse_power_mpreal_hilbert $(MPFR_LIB)
	$(CPP) $(MPFR_INC) cond_mpreal_hilbert.cpp get_sec.o -o cond_mpreal_hilbert $(MPFR_LIB)

iterative_ref: test_iterative_ref_mpreal.cpp test_iterative_ref_qd.cpp test_iterative_ref_dd.cpp get_sec
	$(CPP) $(LAPACKE_INC) -DUSE_MPFR $(MPFR_INC) test_iterative_ref_mpreal.cpp get_sec.o -o test_iterative_ref_mpreal $(LAPACKE_LIB) $(MPFR_LIB) $(QD_LIB)
	$(CPP) $(LAPACKE_INC) test_iterative_ref_dd.cpp get_sec.o -o test_iterative_ref_dd $(LAPACKE_LIB) $(QD_LIB)
	$(CPP) $(LAPACKE_INC) test_iterative_ref_qd.cpp get_sec.o -o test_iterative_ref_qd $(LAPACKE_LIB) $(QD_LIB)

zip:
	-rm *.tar
	-rm *.tar.gz
	tar zcvf mpna.tar.gz *.cpp *.txt *.c *.h Makefile *.inc

clean:
	-rm *.o
	-rm rsa
	-rm quadratic_eq
	-rm complex128
	-rm float128
	-rm mpfr_template
	-rm mpfr_template1
	-rm mpfr_template_cxx
	-rm mpfr_newton_inverse
	-rm mpfr_newton_sqrt
	-rm mpfr_pi_simple
	-rm mpfr_exp
	-rm mpfr_printf
	-rm mpfr_relerr
	-rm mpfr_relerr_cxx
	-rm test_quadmath
	-rm test_arb
	-rm test_mpfi
	-rm quadratic_eq
	-rm logistic
	-rm logistic_arb
	-rm logisitc_mpfi
	-rm logisitc_mpfr
	-rm logisitc_2f
	-rm logistic_dd
	-rm logistic_qd
	-rm logistic_mpreal
	-rm logistic_mpreal_dd
	-rm logistic_mpreal_dd_qd
	-rm matmul_simple
	-rm matmul_simple_omp
	-rm matmul_block
	-rm matmul_block_omp
	-rm matmul_strassen
	-rm matmul_winograd
	-rm cg_mpfr_omp
	-rm mpn_test
	-rm mpn_sample_full
	-rm mpz_test
	-rm mpz_mersenne
	-rm mpz_mersenne_cxx
	-rm mpz_factorial
	-rm mpz_input
	-rm mpz_input_str
	-rm mpq_input
	-rm mpfr_pi_simple
	-rm template_complex
	-rm template_complex_efunc
	-rm mpz_mersenne
	-rm test_mpc
	-rm mpz_factorial
	-rm complex_dd
	-rm mpf_template
	-rm complex_qd
	-rm mpf_template2
	-rm complex_mpreal
	-rm mpf_template3
	-rm complex_d_cpp
	-rm mpf_relerr
	-rm complex_d_c
	-rm mpz_input_nextprime
	-rm dd_test
	-rm mpz_prime_factorization
	-rm dd_print
	-rm mpz_mersenne
	-rm qd_test
	-rm mpz_binomial
	-rm cg_double
	-rm mpz_input_gcd_lcm
	-rm cg_d_omp
	-rm mpq_input
	-rm cg_mpfr
	-rm mpq_input_convert
	-rm cg_mpfr_omp
	-rm cg_d_cxx
	-rm cg_dd_cxx
	-rm cg_qd_cxx
	-rm cg_mpreal_cxx
	-rm lu_d
	-rm lu_d_cxx
	-rm lu_dd
	-rm lu_qd
	-rm lu_mpreal
	-rm power_mpreal_hilbert
	-rm inverse_power_mpreal_hilbert
	-rm cond_mpreal_hilbert
	-rm test_iterative_ref_mpreal
	-rm test_iterative_ref_dd
	-rm test_iterative_ref_qd
	-rm sample_linear
	-rm sample_linear_cxx
