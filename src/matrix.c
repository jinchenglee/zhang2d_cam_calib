/**
* @file matrix.c
* Simple Math		
*					
* @author Marcelo Gattass	
*
* @author Manuel E. L. Fernandez	 
*
* @date Jul06,2006
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

/* Utitlity Macros and Functions */
#define SIGN(a,b) ((b) >= 0. ? fabs(a) : -fabs(a))
#define SQR(a) ((a)*(a))
#define ABS(a) (((a)<0)?-(a):(a))
#define MIN(a,b) (((a)<(b))? (a): (b))
#define MAX(a,b) (((a)>(b))? (a): (b))

#define TOL 1e-10

static void swap(double* a, double* b) 
{
	double tmp=*a;
	*a=*b;
	*b=tmp;
}


/* lib functions */
int mtxGaussAxb (double* a, int n, double* b) 
{
	int    i, j, k;
	int    imax;             /* pivot line */
	double amax, rmax; 

	for (j=0; j<n-1; j++) {  /*  Loop in the columns of [a] */
		rmax = 0.;
		imax = j;
		for (i=j; i<n; i++) {   /* Loop to find the best ration a[i-1)*n-1+j]/a[i-1)*n-1+k] */
			amax = 0.0f;
			for (k=j; k<n; k++)    /* Loop to find largest element of line i */
				if (ABS(a[i*n+k]) > amax)
					amax = ABS(a[i*n+k]);
			if (amax < TOL)        /* Check if all elements are null */
				return 0;             /* no solution */
			else if ((ABS(a[i*n+j]) > rmax*amax) && (ABS(a[i*n+j]) >= TOL)) { /* teste current line */
				rmax = ABS(a[i*n+j]) / amax;
				imax = i;             /* find pivot line */
			}
		}
		if (ABS(a[imax*n+j])<TOL) {         /* Check if pivot is null */
			for (k=imax+1; k<n; k++)          /* Search for a line with a no-null pivot */
				if (ABS(a[k*n+j]) < TOL)
					imax = k;                       /* exchange line j with k */
			if (ABS(a[imax*n+j]) < TOL)
				return 0;                        /* no unique soluition */
		}
		if (imax != j) {                   /* Exchange j by line imax */
			for (k=j; k<n; k++)
				swap(&a[imax*n+k], &a[j*n+k]);
			swap(&b[imax], &b[j]);
		}
		for (i=j+1; i<n; i++) {            /* Clear elements under the diagonal */
			double aux = a[i*n+j] / a[j*n+j];
			for (k=j+1; k<n; k++)             /* Transforms the rest of the elements of the line */
				a[i*n+k] -= aux * a[j*n+k];
			b[i] -= aux * b[j];
		}
	}
	if (ABS(a[(n-1)*n+n-1]) <= TOL)        /* Check the unicity of the solution */
		return 0;                          /* no solution */
	else {
		b[n-1] /= a[(n-1)*n+n-1];          /* back substitution */
		for (i=n-2; i>=0; i--) {           
			for (j=i+1; j<n; j++)
				b[i] -= a[i*n+j] * b[j];
			b[i] /= a[i*n+i];
		}
	}
	return 1;     /* solution ok */                      
}

int mtxSVD(double*a, int m, int n, double* u, double* d, double* v, double* tmp)
{
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;

	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			u[i*n+j]=a[i*n+j];

	g=scale=anorm=0.;
	for (i=0;i<n;i++) {
		l=i+2;
		tmp[i]=scale*g;
		g=s=scale=0.;
		if (i < m) {
			for (k=i;k<m;k++) scale +=  fabs(u[k*n+i]);
			if (scale != 0.) {
				for (k=i;k<m;k++) {
					u[k*n+i] /= scale;
					s += u[k*n+i]*u[k*n+i];
				}
				f=u[i*n+i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i*n+i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.,k=i;k<m;k++) s += u[k*n+i]*u[k*n+j];
					f=s/h;
					for (k=i;k<m;k++) u[k*n+j] += f*u[k*n+i];
				}
				for (k=i;k<m;k++) u[k*n+i] *= scale;
			}
		}
		d[i]=scale *g;
		g=s=scale=0.;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += fabs(u[i*n+k]);
			if (scale!=0.) {
				for (k=l-1;k<n;k++) {
					u[i*n+k] /= scale;
					s += u[i*n+k]*u[i*n+k];
				}
				f=u[i*n+l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				u[i*n+l-1]=f-g;
				for (k=l-1;k<n;k++) tmp[k]=u[i*n+k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.,k=l-1;k<n;k++) s += u[j*n+k]*u[i*n+k];
					for (k=l-1;k<n;k++) u[j*n+k] += s*tmp[k];
				}
				for (k=l-1;k<n;k++) u[i*n+k] *= scale;
			}
		}
		anorm=(anorm>(fabs(d[i])+fabs(tmp[i]))?anorm:(fabs(d[i])+fabs(tmp[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g!=0.) {
				for (j=l;j<n;j++)
					v[j*n+i]=(u[i*n+j]/u[i*n+l])/g;
				for (j=l;j<n;j++) {
					for (s=0.,k=l;k<n;k++) s += u[i*n+k]*v[k*n+j];
					for (k=l;k<n;k++) v[k*n+j] += s*v[k*n+i];
				}
			}
			for (j=l;j<n;j++) v[i*n+j]=v[j*n+i]=0.;
		}
		v[i*n+i]=1.;
		g=tmp[i];
		l=i;
	}
	for (i=(m<n?m:n)-1;i>=0;i--) {
		l=i+1;
		g=d[i];
		for (j=l;j<n;j++) u[i*n+j]=0.;
		if (g != 0.) {
			g=1./g;
			for (j=l;j<n;j++) {
				for (s=0.,k=l;k<m;k++) s += u[k*n+i]*u[k*n+j];
				f=(s/u[i*n+i])*g;
				for (k=i;k<m;k++) u[k*n+j] += f*u[k*n+i];
			}
			for (j=i;j<m;j++) u[j*n+i] *= g;
		} else for (j=i;j<m;j++) u[j*n+i]=0.;
		++u[i*n+i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((fabs(tmp[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((fabs(d[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.;
				s=1.;
				for (i=l;i<k+1;i++) {
					f=s*tmp[i];
					tmp[i]=c*tmp[i];
					if ((fabs(f)+anorm) == anorm) break;
					g=d[i];
					h=sqrt(f*f+g*g);
					d[i]=h;
					h=1./h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=u[j*n+nm];
						z=u[j*n+i];
						u[j*n+nm]=y*c+z*s;
						u[j*n+i]=z*c-y*s;
					}
				}
			}
			z=d[k];
			if (l == k) {
				if (z < 0.) {
					d[k] = -z;
					for (j=0;j<n;j++) v[j*n+k] = -v[j*n+k];
				}
				break;
			}
			if (its == 49) return 0;
			x=d[l];
			nm=k-1;
			y=d[nm];
			g=tmp[nm];
			h=tmp[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0f*h*y);
			g=sqrt(f*f+1.);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=tmp[i];
				y=d[i];
				h=s*g;
				g=c*g;
				z=sqrt(f*f+h*h);
				tmp[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj*n+j];
					z=v[jj*n+i];
					v[jj*n+j]=x*c+z*s;
					v[jj*n+i]=z*c-x*s;
				}
				z=sqrt(f*f+h*h);
				d[j]=z;
				if (z) {
					z=1.0f/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g; 
				for (jj=0;jj<m;jj++) {
					y=u[jj*n+j];
					z=u[jj*n+i];
					u[jj*n+j]=y*c+z*s;
					u[jj*n+i]=z*c-y*s;
				}
			}
			tmp[l]=0.;
			tmp[k]=f;
			d[k]=x;
		}
	}
	return 1;
}

/* 
* Add the tensor product of the vector {v} (i.e., {v}{v}T) 
* to the matrix [A].  
* [A]+={v}{v}T
*
*/
void mtxAddMatVecTensor(double* a, double* v, int n)
{
	int i,j;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) { 
			a[i*n+j]+=v[i]*v[j];
		}
	}
}

/* {x}=[A]{b}    Dimensions: [A]=mxn, {b}=n and {x}=m. --modificado*/
void mtxAb(double* a, double* b, int m, int n, double* x)
{
	int i,j;
	for (i=0;i<m;i++) {
		x[i]=0.;
		for (j=0;j<n;j++)
			x[i]+=a[i*n+j]*b[j];
	}
}

/* {x}=[A]T{b}    Dimensions: [A]=mxn, {b}=m and {x}=n. */
void mtxAtb(double* a, double* b, int m, int n, double* x)
{
	int i,j;
	for (i=0;i<n;i++) {
		x[i]=0.;
		for (j=0;j<m;j++)
			x[i]+=a[j*n+i]*b[j];
	}
}


/* [X]=[A-1)*n-1+B]    Dimensions: [A]=mxp, [B]=pxn and [X]=mxn. */
void mtxAB(double* a, double* b, int m, int p, int n, double* x)
{
	int i,j,k;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			x[i*n+j]=0.;
			for (k=0;k<p;k++)
				x[i*n+j]+=a[i*p+k]*b[k*n+j];
		}
	}
}

/* [X]=[A-1)*n-1+B]T    Dimensions: [A]=mxp, [B]=nxp and [X]=mxn. */
void mtxABt(double* a, double* b, int m, int p, int n, double* x)
{
	int i,j,k;
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++) {
			x[i*n+j]=0.;
			for (k=0;k<p;k++)
				x[i*n+j]+=a[i*p+k]*b[j*p+k];
		}
	}
}

/*  [X]=[A]T[B]    Dimensions: [A]=mxp, [B]=mxn and [X]=pxn. */
void mtxAtB(double* a, double* b, int m, int p, int n, double* x)
{
	int i,j,k;
	for (i=0;i<p;i++) {
		for (j=0;j<n;j++) {
			x[i*n+j]=0.;
			for (k=0;k<m;k++)
				x[i*n+j]+=a[k*p+i]*b[k*n+j];
		}
	}
}

/*  [X]=[A]+[B]    Dimensions: [A]=mxn, [B]=mxn and [X]=mxn. */
void mtxAddMat(double* a, double* b, int m, int n, double* x)
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			x[i*n+j]=a[i*n+j]+b[i*n+j];
}

/*  [X]=[A]-[B]    Dimensions: [A]=mxn, [B]=mxn and [X]=mxn. */
void mtxSubMat(double* a, double* b, int m, int n, double* x)
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			x[i*n+j]=a[i*n+j]-b[i*n+j];
}


/*  [X]=s[A]    Dimensions: [A]=mxn and [X]=mxn. */
void mtxScaMat(double* a, double s,int m, int n, double* x)
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			x[i*n+j]=s*a[i*n+j];
}


/*  [X]=[A]T    Dimensions: [A]=mxn and [X]=nxm. */
void mtxAt(double* a, int m, int n, double* x)
{
	int i,j;
	for (i=0;i<m;i++)
		for (j=0;j<n;j++)
			x[j*m+i]=a[i*n+j];
}

/*  [X]=[A]    Dimensions: [A]=mxn and [X]=mxn. */
void mtxMatCopy(double* a, int m, int n, double* x)
{
	int i,j;
	for (i=0;i<m;i++) 
		for (j=0;j<n;j++) 
			x[i*n+j]=a[i*n+j];
}

/*  {x}={v}+{u}  Dimensions: {v}=m, {u}=m and {x}=m. */
void mtxAddVec(double* u, double* v, int m, double* x)
{
	int i;
	for (i=0;i<m;i++) x[i]=u[i]+v[i];
}

/*  {x}={v}-{u}  Dimensions: {v}=m, {u}=m and {x}=m. */
void mtxSubVec(double* u, double* v, int m, double* x)
{
	int i;
	for (i=0;i<m;i++) x[i]=u[i]-v[i];
}

/*  {x}=s{u}  Dimensions: {u}=m and {x}=m. */
void mtxScaVec(double* u, int m,  double s, double* x)
{
	int i;
	for (i=0;i<m;i++) x[i]=s*u[i];
}

/*  s={v}.{u}  Dimensions: {v}=m and  {u}=m. */
double mtxDotProd(double* u, double* v, int m)
{
	double tmp=0;
	int i;
	for (i=0;i<m;i++) 
		tmp+=u[i]*v[i];
	return tmp;
}

/*  {x}={u}x{v}  Dimensions: {v}=m and  {u}=m. */
void mtxProdVec3(double* u, double* v, double*x)
{
	x[0] = u[1]*v[2] - u[2]*v[1] ;
	x[1] = u[2]*v[0] - u[0]*v[2] ;
	x[2] = u[0]*v[1] - u[1]*v[0] ;
}

double mtxNormalizeVector(int n, double* v)
{
	double norm=mtxDotProd(v,v,n);
	int i;
	if (norm>TOL) {
		norm = sqrt(norm);
		for (i=0;i<n;i++) v[i]/=norm;
	}
	return norm;
}

void mtxCol(double* A,int col,  int m, int n, double* x)
{
	int i;
	for (i=0;i<m;i++){
		x[i]=A[i*n+col];
	}
}

/* {x}={v} Dimensions: {v}={x}=m. */
void mtxVecCopy(double* v, int m, double* x)
{
	int i;
	for(i=0;i<m;i++) x[i]=v[i];
}



/* Print matrix */
void mtxShowMat(char*  title, double* a, int m, int n)
{
	int i,j;
	printf("\n%s",title);
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			if ((j%6)==0) printf("\n");
			printf(" %lf ", a[i*n+j]);
		}
	}
	printf("\n");
}

/* Print vector */
void mtxShowVec(char*  title, double* v, int n)
{
	int j;
	printf("\n%s",title);
	for (j=0;j<n;j++){
		if ((j%6)==0) printf("\n");
		printf("\t%lf ", v[j]);
	}
	printf("\n");
}



int mtxDecompLU ( double* a, int n, int *rp, double* d, double* scale )
{
	int i, imax, j, k;
	double max, tmp, sum;

	rp--;
	scale--;

	*d=1.;                   /* No row interchanges (yet). */
	for (i=0;i<n;i++) {     /* Loop over rows to get the implicit scaling */
		max=0.;          
		for (j=0;j<n;j++)
			max = (max>fabs(a[i*n+j]))?max:fabs(a[i*n+j]);
		if (max == 0.) {
			printf( "\nLUdmcp: Singular Matrix");
			return 0;
		}
		/* Largest nonzero largest element. */
		scale[i+1]=1./max;    /* Save the scaling. */
	}
	for (j=1;j<=n;j++) {  /* This is the loop over columns of Crout's method */
		for (i=1;i<j;i++) {  /* Sum form of a triangular matrix except for i=j */
			sum=a[(i-1)*n-1+j];     
			for (k=1;k<i;k++) sum -= a[(i-1)*n-1+k]*a[(k-1)*n-1+j];
			a[(i-1)*n-1+j]=sum;
		}
		max=0.;  /* Initialize for the search for largest pivot element. */
		imax = 0;        /* Set default value for imax */
		for (i=j;i<=n;i++) {  /* This is i=j of previous sum and i=j+1...N */
			sum=a[(i-1)*n-1+j];      /* of the rest of the sum */
			for (k=1;k<j;k++)
				sum -= a[(i-1)*n-1+k]*a[(k-1)*n-1+j];
			a[(i-1)*n-1+j]=sum;
			if ( (tmp=scale[i]*fabs(sum)) >= max) {
				/* Is the figure of merit for the pivot better than the best so far? */
				max=tmp;
				imax=i;
			}
		}
		if (j != imax) {          /* Do we need to interchange rows? */
			for (k=1;k<=n;k++) {  /* Yes, do so... */
				tmp=a[(imax-1)*n-1+k];
				a[(imax-1)*n-1+k]=a[(j-1)*n-1+k];
				a[(j-1)*n-1+k]=tmp;
			}
			*d = -(*d);           /* ...and change the parity of d */
			scale[imax]=scale[j];       /* Also interchange the scale factor */
		}
		rp[j]=imax;
		if (a[(j-1)*n-1+j] == 0.) a[(j-1)*n-1+j]=TOL;
		/* If the pivot element is zero the matrix is singular (at least */
		/* to the precision of the algorithm). For some applications on */
		/* singular matrices, it is desiderable to substitute TOL for zero */
		if (j != n) {           /* Now, finally divide by pivot element */
			tmp=1./(a[(j-1)*n-1+j]);
			for ( i=j+1 ; i<=n ; i++ ) a[(i-1)*n-1+j] *= tmp;
		}
	}   /* Go back for the next column in the reduction. */

	return 1;
}

void mtxBackSubLU ( double* a, int n, int *rp, double* b )
{
	int i, ii=0, ip, j;
	double sum;
	b--;
	rp--;

	for (i=1;i<=n;i++) {  
		ip=rp[i];     
		sum=b[ip];  
		b[ip]=b[i]; 
		if (ii)  
			for ( j=ii ; j<=i-1 ; j++ ) sum -= a[(i-1)*n-1+j]*b[j];
		else if (sum)      
			ii = i; 
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {   /* Backsubstitution */
		sum=b[i+1];
		for (j=i+1;j<n;j++) sum -= a[i*n+j]*b[j+1];
		b[i+1]=sum/a[i*n+i];  
		if (fabs(b[i])<1.e-8) b[i] = 0.;   
	}
}

double mtxDetLU( double* a, double d, int n)
{
	double det=d;
	int i;
	for (i=0;i<n;i++)
		det*=a[i*n+i];
	return det;
}

/*
* Finds a non trivial solution for the system Ax=0
* A mxn, m>n and rank(A)< n (better if rank(A)= n-1).
*
* The return value indicates the error. He is the ratio from the smallest  
* to the largest singular value of A. If rank(A) < n
* this value should be zero (~=TOL for the implementation).
*/ 
double mtxSVDAx0(double* a, int m, int n, double* x, double* u, double* d, double* v, double* tmp)
{
	double wmax,wmin,wmin2;
	int i,j,jmin;

	/* perform decomposition */
	mtxSVD(a,m,n,u,d,v,tmp);

	/* test if A has a non-trivial solution */
	wmax=d[0]; wmin = d[0];
	for (j=1;j<n;j++){
		if (d[j] < wmin) { wmin=d[j]; jmin =j; }
		if (d[j] > wmax) wmax=d[j];
	}

	/* test for the second smallest singular value */
	wmin2=wmax;
	for (j=0;j<n;j++){
		if (j==jmin) continue;
		if (d[j] < wmin2) wmin2=d[j];
	}

	/* copy the column of V correspondent to the smallest singular value of A */
	for (i=0;i<n;i++)
		x[i] = v[i*n+jmin];

	return (wmin/wmax);
}



#define TINY 1.0e-20

void mtxSVDbacksub(double* u, double* d, double* v, int m, int n, double* b, double* x, double* tmp)
{
	int j,i;
	double s;

	for (j=0;j<n;j++) {
		s=0.0;
		if (d[j]<TINY) {  
			for (i=0;i<m;i++) s += u[i*n+j]*b[i]; /* computes [U]T{b} */
			s /= d[j];   /* multiply by [d]-1 */
		}
		tmp[j]=s;
	}
	for (i=0;i<n;i++) {
		s=0.0;
		for (j=0;j<n;j++) s += v[i*n+j]*tmp[j];  /* computes [V]{tmp} */
		x[i]=s;
	}
}

/*
* Finds a solution for the system Ax=b
* A mxn, m>n.
* 
*/ 
void mtxSVDAxb(double* a, int m, int n, double* x, double* b, double* u, double* d, double* v, double* tmp)
{
	double wmax,wmin;
	int j;

	/* perform decomposition */
	mtxSVD(a,m,n,u,d,v,tmp);

	/* find the larger single value */
	wmax=d[0];
	for (j=1;j<n;j++)
		if (d[j] > wmax) wmax = d[j];

	/* cutoff small values */
	wmin=1.e-6*wmax;
	for (j=0;j<n;j++)
		if (d[j]<wmin) d[j]=0;

	/* backsubstitution */
	mtxSVDbacksub(u,d,v,m,n,b,x,tmp);
}


/*
* Computes the square matrix nxn [X] such that [A][X]=[I],
* then [X] is the inverse of [A]nxn.  
*(uses LU decomposition).
*/ 
void mtxAxInxn(double *a, int n, double *x)
{
	int i,j;
	int *rp = 0;
	double d, *lu = 0, *scale = 0,*col = 0;

	lu = (double*) malloc(sizeof(double)*(n*n));
	rp = (int*) malloc(sizeof(int)*n);
	scale = (double*) malloc(sizeof(double)*n);
	col = (double*) malloc(sizeof(double)*n);

	/* Copy the original matrix inside lu */
	for (i=0;i<n*n;i++)	lu[i]=a[i];

	/* LU decomposition */	
	mtxDecompLU ( lu, n, rp, &d, scale );

	for(j=0;j<n;j++) { //Find inverse by columns.
		for(i=0;i<n;i++) 
			col[i]=0.0;
		col[j]=1.0;
		mtxBackSubLU ( lu, n, rp, col );
		for(i=0;i<n;i++) 
			x[i*n + j]=col[i];
	}

	//free(lu);	//free(rp);	//free(scale); //free(col);	
	return;
}

/*
* Computes the nxm matrix [X] such that [X][A]=[I],
* then [X] is the inverse of [A]mxn.
* (uses LU decomposition).  
*/
void mtxAxImxn(double *A, int m, int n, double *X)
{
	double *At, *AtA, *invAtA, *invAtA_At;

	At = (double*) malloc(sizeof(double)*n*m);
	AtA = (double*) malloc(sizeof(double)*n*n);
	invAtA = (double*) malloc(sizeof(double)*n*n);
	invAtA_At = (double*) malloc(sizeof(double)*n*m);

	mtxAt(A, m, n, At);
	mtxAB(At,A,n,m,n,AtA);
	mtxAxInxn(AtA,n,invAtA);
	mtxAB(invAtA,At,n,n,m,X);

	return;
}
/*
* Get a solution x for the system Ax=b
*										[A]x =				b
*								[At][A]x =	[At]b
*	 inv([At][A])	[At][A]x =	inv([At][A])[At]b	
*										[I]x =  inv([At][A])[At]b	
*/ 
void mtxAxb(double *A, int m, int n, double *b, double *x)
{
	double *At, *AtA, *invAtA, *invAtA_At;

	At = (double*) malloc(sizeof(double)*n*m);
	AtA = (double*) malloc(sizeof(double)*m*m);
	invAtA = (double*) malloc(sizeof(double)*n*n);
	invAtA_At = (double*) malloc(sizeof(double)*n*m);

	mtxAt(A, m, n, At);
	mtxAB(At,A,n,m,n,AtA);
	mtxAxInxn(AtA,n, invAtA);
	mtxAB(invAtA,At,n,n,m,invAtA_At);
	mtxAb(invAtA_At,b,n,m,x);

	return;
}
