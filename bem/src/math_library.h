#ifndef math_library_h
#define math_library_h

/* header file for routines using the math library */
/* one letter prefixes in routine names have the following meanings

      --  no prefix indicates single precision real
    c --  complex data being used 
    d --  double sized data being used
    d_c -- double complex data is being used

*/

#ifndef complex_numbers_h
#include "complex_numbers.h"
#endif

#ifndef math_h
#include <math.h>
#endif

/* Create aliases for NSWC FORTRAN routines so C programs can call them */

//--------------------------------------


#define FFT F77_FUNC(fft,FFT)
#define CEIGV F77_FUNC(ceigv,CEIGV)
#define CMSLV1 F77_FUNC(cmslv1,CMSLV1)
#define CMTMS F77_FUNC(cmtms,CMTMS)
#define DCMSLV F77_FUNC(dcmslv,DCMSLV)
#define MSLV F77_FUNC(mslv,MSLV)
#define DMSLV F77_FUNC(dmslv,DMSLV)
#define MTMS F77_FUNC(mtms,MTMS)
//#define HBRD F77_FUNC(hbrd,HBRD)
//#define HYBRD F77_FUNC(hybrd,HYBRD)
//#define LMDIFF F77_FUNC(lmdiff,LMDIFF)
//#define LMDIF F77_FUNC(lmdif,LMDIF)
#define SGEFA F77_FUNC(sgefa,SGEFA)
#define SGECO F77_FUNC(sgeco,SGECO)
#define SGESL F77_FUNC(sgesl,SGESL)
#define DGEFA F77_FUNC(dgefa,DGEFA)
#define DGESL F77_FUNC(dgesl,DGESL)

/* And create aliases so C routines in math_library.c can be called from
   FORTRAN */

#define c_fft F77_FUNC(c_fft,C_FFT)
#define c_fft_inv F77_FUNC(c_fft_inv,C_FFT_INV)
#define c_set_calc_eigenvalues F77_FUNC(c_set_calc_eigenvalues,C_SET_CALC_EIGENVALUES)
#define c_init_calc_eigenvalues F77_FUNC(c_init_calc_eigenvalues,C_INIT_CALC_EIGENVALUES)
#define c_calc_eigenvalues F77_FUNC(c_calc_eigenvalues,C_CALC_EIGENVALUES)
#define c_mult_matricies F77_FUNC(c_mult_matricies,C_MULT_MATRICIES)
#define mult_matricies F77_FUNC(mult_matricies,MULT_MATRICIES)
#define c_set_invert_matrix F77_FUNC(c_set_invert_matrix,C_SET_INVERT_MATRIX)
#define c_init_invert_matrix F77_FUNC(c_init_invert_matrix,C_INIT_INVERT_MATRIX)
#define c_invert_matrix F77_FUNC(c_invert_matrix,C_INVERT_MATRIX)
#define set_invert_matrix F77_FUNC(set_invert_matrix,SET_INVERT_MATRIX)
#define init_invert_matrix F77_FUNC(init_invert_matrix,INIT_INVERT_MATRIX)
#define invert_matrix F77_FUNC(invert_matrix,INVERT_MATRIX)
#define invert_matrix_cond F77_FUNC(invert_matrix_cond,INVERT_MATRIX_COND)
#define d_set_invert_matrix F77_FUNC(d_set_invert_matrix,D_SET_INVERT_MATRIX)
#define d_init_invert_matrix F77_FUNC(d_init_invert_matrix,D_INIT_INVERT_MATRIX)
#define d_invert_matrix F77_FUNC(d_invert_matrix,D_INVERT_MATRIX)
#define c_set_solve_linear F77_FUNC(c_set_solve_linear,C_SET_SOLVE_LINEAR)
#define c_init_solve_linear F77_FUNC(c_init_solve_linear,C_INIT_SOLVE_LINEAR)
#define c_solve_linear F77_FUNC(c_solve_linear,C_SOLVE_LINEAR)
#define d_c_set_solve_linear F77_FUNC(d_c_set_solve_linear,D_C_SET_SOLVE_LINEAR)
#define d_c_init_solve_linear F77_FUNC(d_c_init_solve_linear,D_C_INIT_SOLVE_LINEAR)
#define d_c_solve_linear F77_FUNC(d_c_solve_linear,D_C_SOLVE_LINEAR)
//#define set_solve_nonlinear F77_FUNC(set_solve_nonlinear,SET_SOLVE_NONLINEAR)
//#define init_solve_nonlinear F77_FUNC(init_solve_nonlinear,INIT_SOLVE_NONLINEAR)
//#define solve_nonlinear F77_FUNC(solve_nonlinear,SOLVE_NONLINEAR)
#define lu_factor F77_FUNC(lu_factor,LU_FACTOR)
#define lu_factor_cond F77_FUNC(lu_factor_cond,LU_FACTOR_COND)
#define lu_solve_linear F77_FUNC(lu_solve_linear,LU_SOLVE_LINEAR)


// C function definitions


extern "C" void c_fft(int *n,COMPLEX *a,COMPLEX *c,int *status);

extern "C" void c_fft_inv(int *n,COMPLEX *a,COMPLEX *c,int *status);


extern "C" void c_set_calc_eigenvalues(int *n);

extern "C" void c_init_calc_eigenvalues(int * status);

extern "C" void c_calc_eigenvalues(COMPLEX *c,int *n,COMPLEX *eval,COMPLEX *evect,
			      int *ldc, int *ldevect, int *status); 


extern "C" void c_mult_matricies(int *ra,int *carb, int *cb,
			   COMPLEX *a, COMPLEX *b, COMPLEX *c,
			   int *lda, int *ldb, int *ldc,int *status);

extern "C" void mult_matricies(int *ra,int *carb, int *cb,
			   float *a, float *b, float *c,
			   int *lda, int *ldb, int *ldc,int *status);


extern "C" void c_set_invert_matrix(int *n);

extern "C" void c_init_invert_matrix(int *status);

extern "C" void c_invert_matrix(int *n,COMPLEX *a,COMPLEX *b,
			   int *lda, int *ldb, int *status);

extern "C" void set_invert_matrix(int *n);

extern "C" void init_invert_matrix(int *status);

extern "C" void invert_matrix(int *n,float *a,float *b,
			   int *lda, int *ldb, int *status);

extern "C" void invert_matrix_cond(int *n,float *a,float *b,
				   int *lda, int *ldb, float *rcond,
				   int *status);

extern "C" void d_set_invert_matrix(int *n);

extern "C" void d_init_invert_matrix(int *status);

extern "C" void d_invert_matrix(int *n,double *a,double *b,
			   int *lda, int *ldb, int *status);


extern "C" void c_set_solve_linear(int *n);

extern "C" void c_init_solve_linear(int *status);

extern "C" void c_solve_linear(int *n, COMPLEX *a, COMPLEX *b,
		      COMPLEX *c, int *lda, int *status);


extern "C" void d_c_set_solve_linear(int *n);

extern "C" void d_c_init_solve_linear(int * status);

extern "C" void d_c_solve_linear(int *n, DOUBLE_COMPLEX *a, DOUBLE_COMPLEX *b,
		      DOUBLE_COMPLEX *c, int *lda, int *status);


//extern "C" void set_solve_nonlinear(int *n);

//extern "C" void init_solve_nonlinear(int *status);

//extern "C" void solve_nonlinear(void (*f)(),void (*f2)(), int *n,float *eps,
//               float *tol, float *x, float *fvec, int *status);


extern "C" void lu_factor(int *n, float *a, float *lu, int *lda, 
		 int *ipvt, int *status);

extern "C" void lu_factor_cond(int *n, float *a, float *lu, int *lda, 
		 int *ipvt, float *rcond, int *status);

extern "C" void lu_solve_linear(int *n, float *a, float *x, float *b, int *lda,
		 int *ipvt, int *status);

extern "C" void dlu_factor(int *n, double *a, double *lu, int *lda, 
		 int *ipvt, int *status);

extern "C" void dlu_solve_linear(int *n, double *a, double *x, double *b,
		 int *lda, int *ipvt, int *status);

/* Declarations of NSWC routines */

extern "C" void FFT(COMPLEX *c, int *n, int *sign, int *ierr);

extern "C" void CEIGV(int *bal,
		 float *c_calc_eigenvalues_ar,float *c_calc_eigenvalues_ai,
		 int *n1, int *n2,
		 float *c_calc_eigenvalues_wr,float *c_calc_eigenvalues_wi,
		 float *c_calc_eigenvalues_zr,float *c_calc_eigenvalues_zi,
		 int *ierr,
		 float *c_calc_eigenvalues_temp);

extern "C"  void CMTMS(int *ra, int *carb, int *cb, COMPLEX *a, int *lda,
			   COMPLEX *b, int *ldb, COMPLEX *c, int *ldc);

extern "C"  void MTMS(int *ra,int *carb,int *cb,float *a,int *lda,
			  float *b,int *ldb,float *c,int *ldc);

extern "C"  void CMSLV1(int *calc_inv,int *n,int *zero_dim1,
			    COMPLEX *c_solve_linear_a, int *lda,
			    COMPLEX *c_solve_linear_b, 
			    int *n2,
			    int *ierr,
			    int *c_invert_matrix_ipvt,
			    COMPLEX *c_invert_matrix_wrk);

extern "C"  void MSLV(int *calc_inv,int *n,int *zero_dim1,
			  float *b,int *ldb,int *dum,int *zero_dim2,
			  float *t1, float *rcond,int *ierr,
			  int *invert_matrix_ipvt,float *invert_matrix_wrk);

extern "C"  void DMSLV(int *calc_inv,int *n,int *zero_dim1,
			  double *b,int *ldb,int *dum,int *zero_dim2,
			  double *t1, double *rcond,int *ierr,
			  int *invert_matrix_ipvt,double *invert_matrix_wrk);

extern "C"  void DCMSLV(int *calc_inv,
                   int *n1,
			    int *one_dim,
			    double *d_c_solve_linear_ar,
			    double *d_c_solve_linear_ai,
			    int *lda,
			    double *d_c_solve_linear_br,
			    double *d_c_solve_linear_bi,
			    int *n2,
			    int *ierr,
			    int *d_c_solve_linear_ipvt,
			    double *d_c_solve_linear_wrk);

extern "C" void SGEFA(float *lu, int *lda, int *n, int *ipvt, int *info);

extern "C" void SGECO(float *lu, int *lda, int *n, int *ipvt, float *rcond, 
                 int *info);

extern "C" void SGESL(float *a, int *lda, int *n, int *ipvt, float *x, int *job);

extern "C" void DGEFA(double *lu, int *lda, int *n, int *ipvt, int *info);

extern "C" void DGESL(double *a, int *lda, int *n, int *ipvt, double *x, int *job);


#endif
