#ifndef __CMINPACK_H
#define __CMINPACK_H

#if	defined(__cplusplus)
extern "C" {
#endif

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)	(((a) < (b)) ? (a) : (b))

	double enorm(int n, double x[]);
	double rownorm(int m, int r, int c, double **x);
	double colnorm(int m, int r, int c, double **x);

	void fdjac2(void f(int,int, double*,double*,int*),int m,int n,double x[],double fvec[],double **fjac,
		int *iflag,double epsfcn,double wa[]);

	void lmpar(int n,double **fjac,int ipvt[],double diag[],
		double qtf[],double delta,double *par,double wa1[],
		double wa2[],double wa3[],double wa4[]);

	void qrfac(int m,int n,double **a,int pivot,int ipvt[],
		double rdiag[],double acnorm[],double wa[]);

	void lmdif(void f(int,int, double*,double*,int*),int m,int n,double x[],int msk[],double fvec[],
		double ftol,double xtol,double gtol,int maxfev,double epsfcn,
		double diag[],int mode,double factor,int *info,int *nfev,
		double **fjac,int ipvt[],double qtf[],double wa1[],double wa2[],
		double wa3[],double wa4[]);

	void qrsolv(int n,double **r,int ipvt[],double diag[],
		double qtb[],double x[],double sdiag[],double wa[]);

	int lmdif0(void f(int,int, double*,double*,int*),int m,int n,double x[],int msk[],double fvec[],
		double tol,int *info,int *nfev);


#if defined(__cplusplus)
}
#endif

#endif  /* __CMINPACK_H */

/* CMINPACK.H */
