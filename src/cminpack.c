
#include <stdlib.h>
#include <math.h>
#include "cminpack.h"

double dpmpar[] = { 2.220446049250313e-16,
									2.225073858507201e-308,
									1.797693134862315e308 };

double enorm(int n, double x[])
{
	int i;
	double sum;

	sum = x[0] * x[0];
	for (i=1;i<n;i++)
		sum += x[i]*x[i];
	return sqrt(sum);
}

/* compute row norm for column c from row r to row m */
double rownorm(int m, int r, int c, double **x)
{
	int i;
	double sum;
	sum = x[r][c] * x[r][c];
	for (i = r+1;i < m; i++)
		sum += x[i][c] * x[i][c];

	return sqrt(sum);
}

double colnorm(int m, int r, int c, double **x)
{
	int i;
	double sum;
	sum = x[r][c] * x[r][c];
	for (i = c+1;i < m; i++)
		sum += x[r][i] * x[r][i];

	return sqrt(sum);
}


void fdjac2(void f(int,int,double *,double *,int *),int m,int n,double x[],double fvec[],double **fjac,
						int *iflag,double epsfcn,double wa[])
{
	int i,j;
	double eps,epsmch,h,temp;

	epsmch = (epsfcn > dpmpar[0]) ? epsfcn : dpmpar[0];
	eps = sqrt(epsmch);

	for (j = 0;j < n; j++) {
		temp = x[j];
		if (temp == 0.0) h = eps;
		else h = eps * fabs(temp);
		x[j] = temp + h;
		f(m,n,x,wa,iflag);
		if (*iflag < 0) break;
		x[j] = temp;
		for (i = 0;i < m; i++)
			fjac[j][i] = (wa[i] - fvec[i]) / h;
	}
}

void lmpar(int n,double **r,int ipvt[],double diag[],double qtb[],
					 double delta,double *par,double x[],double sdiag[],double wa1[],
					 double wa2[])
{
	int i,iter,j,jp1,k,l,nsing;
	double dxnorm,dwarf,fp,gnorm,parc,parl,paru;
	double sum,temp;

	dwarf = dpmpar[1];
	nsing = n;
	for (j = 0; j < n; j++) {
		wa1[j] = qtb[j];
		if ((r[j][j] == 0.0) && (nsing == n))
			nsing = j;
		if (nsing < n) wa1[j] = 0.0;
	}
	if (nsing >= 1) {
		for (k = 0;k < nsing; k++) {
			j = nsing - k - 1;
			wa1[j] /= r[j][j];
			temp = wa1[j];
			if (j < 1) continue;
			for (i = 0; i < j; i++)
				wa1[i] -= r[j][i] * temp;
		}
	}
	for (j = 0; j < n; j++) {
		l = ipvt[j];
		x[l] = wa1[j];
	}
	iter = 0;
	for (j = 0; j < n; j++)
		wa2[j] = diag[j] * x[j];
	dxnorm = enorm(n,wa2);
	fp = dxnorm - delta;
	if (fp <= 0.1*delta) {
		if (iter == 0)
			*par = 0.0;
		return;
	}
	parl = 0.0;
	if (nsing >= n) {
		for (j = 0; j < n; j++) {
			l = ipvt[j];
			wa1[j] = diag[l] * wa2[l] / dxnorm;
		}
		for (j = 0; j < n; j++) {
			sum = 0.0;
			if (j >= 1) {
				for (i = 0; i < j; i++)
					sum += r[j][i] * wa1[i];
			}
			wa1[j] = (wa1[j] - sum) / r[j][j];
		}
		temp = enorm(n,wa1);
		parl = ((fp / delta) / temp) / temp;
	}
	for (j = 0;j < n; j++) {
		sum = 0.0;
		for (i = 0; i <= j; i++)
			sum += r[j][i] * qtb[i];
		l = ipvt[j];
		wa1[j] = sum / diag[l];
	}
	gnorm = enorm(n,wa1);
	paru = gnorm / delta;
	if (paru == 0.0)
		paru = dwarf / min(delta,0.1);
	*par = max(*par,parl);
	*par = min(*par,paru);
	if (*par == 0.0)
		*par = gnorm / dxnorm;
	while (1) {
		iter++;
		if (*par == 0.0)
			*par = max(dwarf,0.001 * paru);
		temp = sqrt(*par);
		for (j = 0; j < n; j++)
			wa1[j] = temp * diag[j];
		qrsolv(n,r,ipvt,wa1,qtb,x,sdiag,wa2);
		for (j = 0; j < n; j++)
			wa2[j] = diag[j] * x[j];
		dxnorm = enorm(n,wa2);
		temp = fp;
		fp = dxnorm - delta;

		if ((fabs(fp) <= 0.1*delta) || (parl == 0.0) && (fp <= temp)
			&& (temp > 0.0) || iter == 10) {
				if (iter == 0)
					*par = 0.0;
				return;
			}
			for (j = 0; j < n; j++) {
				l = ipvt[j];
				wa1[j] = diag[l] * wa2[l] / dxnorm;
			}
			for (j = 0; j < n; j++) {
				wa1[j] /= sdiag[j];
				temp = wa1[j];
				jp1 = j + 1;
				if (jp1 < n)
					for (i = jp1; i < n; i++)
						wa1[i] -= r[j][i] * temp;
			}
			temp = enorm(n,wa1);
			parc = ((fp/delta) / temp) / temp;
			if (fp > 0.0)
				parl = max(parl,*par);
			if (fp < 0.0)
				paru = min(paru,*par);
			*par = max(parl,*par+parc);
	}
}

void qrfac(int m,int n,double **a,int pivot,int ipvt[],
					 double rdiag[],double acnorm[],double wa[])
{
	int i,j,jp1,k,kmax,minmn;
	double ajnorm,epsmch,sum,temp;

	/* get machine precision */
	epsmch = dpmpar[0];
	/* compute the initial column norms and initialize several arrays */
	for (j = 0;j < n; j++) {
		acnorm[j] = colnorm(m,j,0,a);
		rdiag[j] = acnorm[j];
		wa[j] = rdiag[j];
		if (pivot) ipvt[j] = j;
	}
	/* reduce a to r with householder transformations */
	minmn = (m < n) ? m : n;
	for (j = 0;j < minmn; j++) {
		if (pivot) {
			/* bring column with largest norm into the pivot position */
			kmax = j;
			for  (k = j;k < n; k++)
				if (rdiag[k] > rdiag[kmax]) kmax = k;
			if (kmax != j) {
				for (i = 0;i < m; i++) {
					temp = a[j][i];
					a[j][i] = a[kmax][i];
					a[kmax][i] = temp;
				}
				rdiag[kmax] = rdiag[j];
				wa[kmax] = wa[j];
				k = ipvt[j];
				ipvt[j] = ipvt[kmax];
				ipvt[kmax] = k;
			}
		}
		/* compute the householder transformation */
		ajnorm = colnorm(m,j,j,a);
		if (ajnorm != 0.0) {
			if (a[j][j] < 0.0) ajnorm = -ajnorm;
			for (i = j;i < m; i++)
				a[j][i] /= ajnorm;
			a[j][j] += 1.0;
			jp1 = j + 1;
			if (n > jp1) {
				for (k = jp1;k < n; k++) {
					sum = 0.0;
					for (i = j;i < m; i++)
						sum += a[j][i]*a[k][i];
					temp = sum / a[j][j];
					for (i = j; i < m; i++)
						a[k][i] -=temp*a[j][i];
					if (!pivot || !rdiag[k]) continue;
					temp = a[k][j] / rdiag[k];
					rdiag[k] *= sqrt(max(0.0,1.0-temp*temp));
					if (0.5 * (rdiag[k] * rdiag[k] / (wa[k] * wa[k])) > epsmch) continue;
					rdiag[k] = colnorm(m,k,jp1,a);
					wa[k] = rdiag[k];
				}
			}
		}
		rdiag[j] = -ajnorm;
	}
}

void lmdif(void f(int,int, double*,double*,int*),int m,
					 int n,double x[],int msk[],double fvec[],double ftol,
					 double xtol,double gtol,int maxfev,double epsfcn,double diag[],
					 int mode,double factor,int *info,int *nfev,double **fjac,
					 int ipvt[],double qtf[],double wa1[],double wa2[],
					 double wa3[],double wa4[])
{
	int i,iflag,iter,j,l;
	double actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm;
	double par,pnorm,prered,ratio,sum,temp,temp1,temp2,xnorm;

	/* initialize */
	epsmch = dpmpar[0];

	*info = 0;
	iflag = 0;
	*nfev = 0;

	/* check for input parameter errors */
	if ((n <= 0) || (m < n) || (maxfev <= 0)
		|| (factor <= 0)) return;
	if (mode == 2) {
		for (j = 0; j < n; j++)
			if (diag[j] <= 0) return;
	}

	/* evaluate the function at the starting point and calculate its norm */   
	iflag = 1;
	f(m,n,x,fvec,&iflag);
	*nfev = 1;
	if (iflag < 0) {
		*info = iflag;
		return;
	}
	fnorm = enorm(m,fvec);

	/* initialize levenberg-marquardt counters */
	par = 0;
	iter = 1;

	/* outer loop */
	while(1) {
		/* calculate jacobian matrix */
		iflag = 2;
		fdjac2(f,m,n,x,fvec,fjac,&iflag,epsfcn,wa4);
		*nfev += n;
		if (iflag < 0) {
			*info = iflag;
			return;
		} 
		f(m,n,x,fvec,&iflag);
		/* compute the qr factorization of the jacobian */
		qrfac(m,n,fjac,1,ipvt,wa1,wa2,wa3);                
		if (iter == 1) {
			if (mode != 2) {
				for (j = 0;j < n; j++) {
					diag[j] = wa2[j];
					if (wa2[j] == 0.0) diag[j] = 1.0;
				}
			}
			for (j = 0;j < n; j++)
				wa3[j] = diag[j] * x[j];
			xnorm = enorm(n,wa3);
			delta = factor * xnorm;
			if (delta == 0) delta = factor;
		}
		for (i = 0; i < m; i++)
			wa4[i] = fvec[i];
		for (j = 0;j < n; j++) {
			if (fjac[j][j] != 0.0) {
				sum = 0.0;
				for (i = j;i < m; i++)
					sum += fjac[j][i] * wa4[i];
				temp = -sum / fjac[j][j];
				for (i = j; i < m; i++)
					wa4[i] += fjac[j][i] * temp;
			}
			fjac[j][j] = wa1[j];
			qtf[j] = wa4[j];
		}
		/* compute the norm of the scaled gradient */ 
		gnorm = 0.0;
		if (fnorm != 0.0) {
			for (j = 0; j < n; j++) {
				l = ipvt[j];
				if (wa2[l] == 0.0) continue;
				sum = 0.0;
				for (i = l; i <= j; i++)
					sum += fjac[j][i] * qtf[i] / fnorm;
				gnorm = max(gnorm,fabs(sum/wa2[l]));
			}
		}
		/* test for convergence of the gradient norm */
		if (gnorm <= gtol) *info = 4;
		if (*info != 0) {
			*info = iflag;
			return;
		}
		/* rescale if necessary */
		if (mode != 2) {
			for (j = 0; j < n; j++)
				diag[j] = max(diag[j],wa2[j]);
		}
		/* beginning of inner loop */
		do {
			/* determine the levenberg-marquardt parameter */
			lmpar(n,fjac,ipvt,diag,qtf,delta,&par,wa1,wa2,wa3,wa4);
			for (j = 0;j < n; j++) {
				wa1[j] = -wa1[j];
				wa2[j] = x[j] + wa1[j];
				wa3[j] = diag[j] * wa1[j];
			}
			pnorm = enorm(n,wa3);
			if (iter == 1) delta = min(delta,pnorm);
			iflag = 1;
			f(m,n,wa2,wa4,&iflag);
			(*nfev)++;
			if (iflag < 0) {
				*info = iflag;
				return;
			}
			fnorm1 = enorm(m,wa4);
			actred = -1.0;
			if (0.1 * fnorm1 < fnorm)
				actred = 1.0 - (fnorm1*fnorm1/(fnorm*fnorm));
			for (j = 0;j < n; j++) {
				wa3[j] = 0.0;
				l = ipvt[j];
				temp = wa1[l];
				for (i = 0; i <= j; i++)
					wa3[i] += fjac[j][i] * temp;
			}
			temp1 = enorm(n,wa3) / fnorm;
			temp2 = sqrt(par) * pnorm / fnorm;
			prered = temp1*temp1 + temp2*temp2 / 0.5;
			dirder = -(temp1*temp1 + temp2*temp2);
			ratio = 0.0;
			if (prered != 0.0) ratio = actred/prered;
			if (ratio <= 0.25) {
				if (actred > 0.0) temp = 0.5;
				if (actred < 0.0) temp = 0.5*dirder/(dirder+0.5*actred);
				delta = temp * min(delta,pnorm/0.1);
				par /= temp;
			}
			else {
				if ((par == 0.0) || (ratio >= 0.75)) {
					delta = pnorm / 0.5;
					par *= 0.5;
				}
			}
			if (ratio >= 0.0001) {
				for (j = 0; j < n; j++) {
					if (msk[j]) {           /* handle masked variables */
						x[j] = wa2[j];
					}
					wa2[j] = diag[j] * x[j];
				}
				for (i = 0; i < m; i++)
					fvec[i] = wa4[i];
				xnorm = enorm(n,wa2);
				fnorm = fnorm1;
				iter++;
			}
			if ((fabs(actred) <= ftol) && (prered <= ftol) &&
				(0.5*ratio <= 1.0)) *info = 1;
			if (delta <= xtol*xnorm) *info = 2;
			if ((fabs(actred) <= ftol) && (prered <= ftol) &&
				(0.5*ratio <= 1.0) && (*info == 2)) *info = 3;
			if (*nfev >= maxfev) *info = 5;
			if ((fabs(actred) <= epsmch) && (prered <= epsmch) &&
				(0.5*ratio <= 1.0)) *info = 6;
			if (delta <= epsmch*xnorm) *info = 7;
			if (gnorm <= epsmch) *info = 8;
			if (*info != 0) {
				*info = iflag;
				return;
			}
		} while (ratio <= 0.0001);
			
	}
}


void qrsolv(int n,double **r,int ipvt[],double diag[],
						double qtb[],double x[],double sdiag[],double wa[])
{
	int i,j,jp1,k,kp1,l,nsing;
	double qtbpj,sum,temp,dsin,dcos,dtan,dcotan;

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++)
			r[j][i] = r[i][j];
		x[j] = r[j][j];
		wa[j] = qtb[j];
	}
	for (j = 0; j < n; j++) {
		l = ipvt[j];
		if (diag[l] != 0.0) {
			for (k = j; k < n; k++)
				sdiag[k] = 0.0;
			sdiag[j] = diag[l];
			qtbpj = 0.0;
			for (k = j; k < n; k++) {
				if (sdiag[k] == 0.0) continue;
				if (fabs(r[k][k]) < fabs(sdiag[k])) {
					dcotan = r[k][k] / sdiag[k];
					dsin = 1.0 / sqrt(1.0 + dcotan * dcotan);
					dcos = dsin * dcotan;
				}
				else {
					dtan = sdiag[k] / r[k][k];
					dcos = 1.0 / sqrt(1.0 + dtan * dtan);
					dsin = dcos * dtan;
				}
				r[k][k] = dcos * r[k][k] + dsin * sdiag[k];
				temp = dcos * wa[k] + dsin * qtbpj;
				qtbpj = -dsin * wa[k] + dcos * qtbpj;
				wa[k] = temp;
				kp1 = k + 1;
				if (n <= kp1) continue;
				for (i = kp1; i < n; i++) {
					temp = dcos * r[k][i] + dsin * sdiag[i];
					sdiag[i] = -dsin * r[k][i] + dcos * sdiag[i];
					r[k][i] = temp;
				}
			}
		}
		sdiag[j] = r[j][j];
		r[j][j] = x[j];
	}
	nsing = n;
	for (j = 0; j < n; j++) {
		if ((sdiag[j] == 0.0) && (nsing == n))
			nsing = j;
		if (nsing < n)
			wa[j] = 0.0;
	}
	if (nsing >= 1) {
		for (k = 0; k < nsing; k++) {
			j = nsing - k - 1;
			sum = 0.0;
			jp1 = j + 1;
			if (nsing > jp1) {
				for (i = jp1; i < nsing; i++)
					sum += r[j][i] * wa[i];
			}
			wa[j] = (wa[j] - sum) / sdiag[j];
		}
	}
	for (j = 0; j < n; j++) {
		l = ipvt[j];
		x[l] = wa[j];
	}
}

int lmdif0(void fcn(int,int, double*,double*,int*),int m, int n,double x[],int msk[],
					 double fvec[],double tol,int *info,int *nfev)
{
	int j,maxfev,mode;
	int *ipvt;
	double ftol,xtol,gtol,epsfcn,factor;
	double *diag,**fjac,*qtf,*wa1,*wa2,*wa3,*wa4;

	/* Check input parameters */
	if (n <= 0 || m < n || tol < 0.0) {
		*info = 0;
		return(1);
	}
	/* Allocate memory for working arrays. */
	ipvt = (int *)malloc(n*sizeof(int));
	diag = (double *)malloc(n*sizeof(double));
	qtf = (double *)malloc(n*sizeof(double));
	wa1 = (double *)malloc(n*sizeof(double));
	wa2 = (double *)malloc(n*sizeof(double));
	wa3 = (double *)malloc(n*sizeof(double));
	wa4 = (double *)malloc(m*sizeof(double));


	/* Create 2d matrix for Jacobian */
	fjac = (double **)malloc(n*sizeof(double *));
	for (j=0;j<n;j++)
		fjac[j] = (double *)malloc(m*sizeof(double));

	/* Set convergence tolerances */
	ftol = 1.0E-5; 
	xtol = 1.01E-7;
	gtol = 0.0;

	maxfev = n*1000;
	epsfcn = 1.01E-16;
	mode = 1;
	factor = 100;
	*nfev = 0;

	lmdif(fcn,m,n,x,msk,fvec,ftol,xtol,gtol,maxfev,epsfcn,diag,mode,
		factor,info,nfev,fjac,ipvt,qtf,wa1,wa2,wa3,wa4);

	if (*info == 8) *info = 4;
	for (j=0;j<n;j++)
		//free(fjac[j]);
	//free(fjac);
	//free(wa4);
	//free(wa3);
	//free(wa2);
	//free(wa1);
	//free(qtf);
	//free(diag);
	//free(ipvt);
	return(0);
}
