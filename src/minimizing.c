#include "matrix.h"
#include <memory.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "cminpack.h"

#define SQR(a) ((a)*(a))
#define ABS(a) (((a)<0)?-(a):(a))
#define MIN(a,b) (((a)<(b))? (a): (b))
#define MAX(a,b) (((a)>(b))? (a): (b))

#define PRINT 0

double* GB_modelPoints ;
double* GB_photoPoints ;

double erroQuadratico(double* photo,double* proj)
{
	return SQR(photo[0] - proj[0]/proj[2]) + SQR(photo[1] - proj[1]/proj[2]);
}

double calcErro(double* matrix, double* model, double* photo,int numPoints)
{
	int i;	
	double res[3];
	double sumError=0;
	for(i=0;i<numPoints;i++)
	{			
		mtxAb(matrix,&(model[3*i]),3,3,res);	
		sumError +=SQR(erroQuadratico(&(photo[2*i]),res));
	}	
	return sumError;
}




void erro_ajuste_cminpack (int m, int n,	double* params, double* err,int *iflag)
/* *m_ptr -->  m : pointer to number of points to fit */
/* *n_ptr -->  n : pointer to number of parameters */
/* *params -->  params :  vector of parameters */
/* *err -->  vector of error from data */
{
	int       k;	
	double res[3];
	for(k=0; k< m; k++) {		
		mtxAb(params,&(GB_modelPoints[3*k]),3,3,res);	
		err[k] =erroQuadratico(&(GB_photoPoints[2*k]),res);

	}
}

void _lmHomography_cminpack(double* matrix,int numPoints,double *ret)
{

	int     i;
	int     NPARAMS  = 9;
	// Parameters needed by MINPACK's lmdif() 
	int     m = numPoints;
	int     n = NPARAMS;

	double  *x;    //	x[NPARAMS];
	int     *msk;  // msk[NPARAMS];
	double  *fvec; // fvec[m];
	double  tol = 1.0e-10;
	int     info;
	int     nfev;

	// allocate some workspace 
	if ( (x = (double*)malloc(NPARAMS * sizeof(double))) == NULL ) 
	{
		printf("malloc: Cannot allocate workspace x\n");
		exit(-1);
	}

	if ( (msk = (int*)malloc(NPARAMS * sizeof(int))) == NULL ) 
	{
		printf("malloc: Cannot allocate workspace msk\n");
		exit(-1);
	}

	if ( (fvec = (double*)malloc(m * sizeof(double))) == NULL ) 
	{
		printf("malloc: Cannot allocate workspace fvec\n");
		exit(-1);
	}


	for(i = 0; i < NPARAMS; i++) x[i] = matrix[i];

	if(lmdif0(erro_ajuste_cminpack,
		m,n,x,msk,fvec,tol,&info,&nfev)){
			printf("\nParameters error lmdif process \n");
		}

		for(i = 0; i < NPARAMS; i++) ret[i] = x[i];

		// release allocated workspace 
		free (fvec);
		free (msk);
		free (x);

		if(PRINT)	// print the number of function calls during iteration 
			printf("\ninfo: %d nfev: %d\n",info,nfev);


}

void lmHomography(double *photoPoints,double *modelPoints,double *initial,int numPoints,double *lm)
{
	GB_photoPoints = photoPoints;
	GB_modelPoints = modelPoints;

	if(PRINT)	mtxShowMat("homografia inicial",initial,3,3);

	_lmHomography_cminpack(initial,numPoints,lm);

	if(PRINT)	mtxShowMat("homografia lm",lm,3,3);

}

void erro_ajuste_ARt_cminpack (int m, int n,	double* params, double* err,int *iflag)
{
	/* *m_ptr -->  m : pointer to number of points to fit */
	/* *n_ptr -->  n : pointer to number of parameters */
	/* *params -->  params :  vector of parameters */
	/* *err -->  vector of error from data */
	int    j,k;	
	double res[3];
	double model[4];

	double A[9],Rt[12];
	double ARt[12];

	j=0;
	A[0] = params[j++]; A[1] = params[j++]; A[2] = params[j++];
	A[3] = 0.0;    A[4] = params[j++]; A[5] = params[j++];
	A[6] = 0.0;    A[7] = 0.0;    A[8] = 1.0;
	
	Rt[0] = params[j++]; Rt[1] = params[j++]; Rt[2] = params[j++];	
	Rt[4] = params[j++]; Rt[5] = params[j++];	Rt[6] = params[j++];       
	Rt[8] = params[j++]; Rt[9] = params[j++];	Rt[10] = params[j++];       

	Rt[3] = params[j++];    
	Rt[7] = params[j++];
	Rt[11] = params[j++];

	mtxAB(A,Rt,3,3,4,ARt);	

	for(k=0; k< m; k++) {

		model[0] = GB_modelPoints[3*k];
		model[1] = GB_modelPoints[(3*k)+1];
		model[2] = 0.0;
		model[3] = 1.0;

		mtxAb(ARt,model,3,4,res);	
		
		err[k] =erroQuadratico(&(GB_photoPoints[2*k]),res);
	}
}

void _lmARt_cminpack(double *Ai,double *Alm,
										 double *Ri,double *Rlm,
										 double *ti,double *tlm,
										 int numPoints)
{

	int     i,j;
	int     NPARAMS  = 17;
	// Parameters needed by MINPACK's lmdif() 
	int     m = numPoints;
	int     n = NPARAMS;

	double  *x;    //	x[NPARAMS];
	int     *msk;  // msk[NPARAMS];
	double  *fvec; // fvec[m];
	double  tol = 1.0e-10;
	int     info;
	int     nfev;

	// allocate some workspace 
	if ( (x = (double*)malloc(NPARAMS * sizeof(double))) == NULL ) 
	{
		printf("malloc: Cannot allocate workspace x\n");
		exit(-1);
	}

	if ( (msk = (int*)malloc(NPARAMS * sizeof(int))) == NULL ) 
	{
		printf("malloc: Cannot allocate workspace msk\n");
		exit(-1);
	}

	if ( (fvec = (double*)malloc(m * sizeof(double))) == NULL ) 
	{
		printf("malloc: Cannot allocate workspace fvec\n");
		exit(-1);
	}


	j=0;
	for(i = 0; i < 5; i++) 		x[j++] = Ai[i];
	for(i = 0; i < 9; i++) 		x[j++] = Ri[i];
	for(i = 0; i < 3; i++) 		x[j++] = ti[i];
		
	if(lmdif0(erro_ajuste_ARt_cminpack,
		m,n,x,msk,fvec,tol,&info,&nfev)){
			printf("\nParameters error lmdif process \n");
		}

		j=0;
		for(i = 0; i < 5; i++) 		Alm[i] = x[j++];
		for(i = 0; i < 9; i++) 		Rlm[i] = x[j++];
		for(i = 0; i < 3; i++) 		tlm[i] = x[j++];

		// release allocated workspace 
		free (fvec);
		free (msk);
		free (x);

#ifdef _DEBUG
		// print the number of function calls during iteration 
		printf("\ninfo: %d nfev: %d\n",info,nfev);
#endif

}


void lmARt(double *photoPoints,double *modelPoints,
										double *Ai,double *Alm,double *Ri,double *Rlm,double *ti,double *tlm,
										int numPoints)
{
	GB_photoPoints = photoPoints;
	GB_modelPoints = modelPoints;
	_lmARt_cminpack(Ai,Alm,Ri,Rlm,ti,tlm,numPoints);
}

