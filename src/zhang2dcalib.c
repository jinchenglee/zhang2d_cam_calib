#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minimizing.h"
#include "matrix.h"
#include "m3.h"
#include "cminpack.h"

#define PRINT 0
#define TOL 1.E-8
#define NORMALIZE 1
 

/* Solve equation system Ax = 0 */
void solveAx0(int m, int n, double* A, double* x) 
{
	/* uses the Singular Value Decomposition of A, i.ie, A=UDVt */
	double* U   =(double*) malloc(sizeof(double)*m*n); /* [U]mxn */
	double* D   =(double*) malloc(sizeof(double)*n);   /* [D]nxn, but only the diagonal elements are stored */
	double* V   =(double*) malloc(sizeof(double)*n*n); /* [V]nxn */
	double* tmp =(double*) malloc(sizeof(double)*m);   /*  temporary area 1xm */

	mtxSVDAx0(A,m,n,x,U,D,V,tmp);

	if (1) {
	mtxScaVec(D,n,1/D[0],D);
	mtxShowVec("\nSingular Values:",D,n);
	}
	free(tmp);
	free(V);
	free(D);
	free(U);
}

/* Normalize image points that belong to calibration pattern. */
void normalizeImagePoints(double w, double h, int n,int views,double* imagePoints, double* N)
{
	double sx=2.0/w;
	double sy=2.0/h;
	double x0=w/2.0;
	double y0=h/2.0;

	int i,k;
	for (k=0;k<views;k++) {
		for (i=0;i<n;i++) {
			imagePoints[k*2*n+2*i  ] = sx*(imagePoints[k*2*n+2*i  ]-x0);
			imagePoints[k*2*n+2*i+1] = sy*(imagePoints[k*2*n+2*i+1]-y0);
		}
	}
	N[0]=sx; N[1]=0;  N[2]=-sx*x0;
	N[3]=0;  N[4]=sy; N[5]=-sy*y0;
	N[6]=0;  N[7]=0;  N[8]=1;
}

/* Get homography between model and image points. */
int  homography(int nPoints, double* modelPoints, double* imagePoints, double* H)
{
	int k;	
	double lm_cminpack[9];

	double* L=(double*) malloc(2*nPoints*9*sizeof(double));  /* L is a 2nx9 matrix where Lij is in L[9*i+j] */

	/* Assembles coefficients matrix L */
	for(k=0; k<nPoints; k++) {
		double X=modelPoints[3*k+0];  /* X coord of model point k */
		double Y=modelPoints[3*k+1];  /* Y coord of model point k */
		double W=modelPoints[3*k+2];  /* W coord of model point k */
		double u=imagePoints[2*k+0];  /* u coord of image point k */
		double v=imagePoints[2*k+1];  /* v coord of image point k */
		int    i=2*k;                 /* line number in matrix L  */
		L[9*i+0] =    X; L[9*i+1] =    Y; L[9*i+2] =    W; 
		L[9*i+3] =    0; L[9*i+4] =    0; L[9*i+5] =    0; 
		L[9*i+6] = -u*X; L[9*i+7] = -u*Y; L[9*i+8] = -u*W; 
		i=2*k+1;         
		L[9*i+0] =    0; L[9*i+1] =    0; L[9*i+2] =    0; 
		L[9*i+3] =    X; L[9*i+4] =    Y; L[9*i+5] =    W; 
		L[9*i+6] = -v*X; L[9*i+7] = -v*Y; L[9*i+8] = -v*W;
	}

	if (PRINT) mtxShowMat("[L]",L,2*nPoints,9);
	solveAx0(2*nPoints,9,L,H); /* solves the system [L]{h}={0} */
	if (PRINT) mtxShowMat("[H]",H,3,3);


	lmHomography(imagePoints,modelPoints,H,nPoints,lm_cminpack);
	if (PRINT) mtxShowMat("\n[H - lm_cminpack]",lm_cminpack,3,3);

	
	mtxMatCopy(lm_cminpack,3,3,H);
	

	free(L);
	return 0;
}

/**/
double hTransform(double* h, double X, double Y, double W, double* u, double* v)
{
	double z=h[6]*X+h[7]*Y+h[8]*W;
	*u = h[0]*X+h[1]*Y+h[2]*W;
	*v = h[3]*X+h[4]*Y+h[5]*W;
	if (fabs(z)>1.E-6) {
		*u = (*u)/z;
		*v = (*v)/z;
	}
	return z;
}

/**/
void calcB(int nH, double* H, double* B)
{	

	int m = 2*nH;
	int n = 6;
	double* V =(double*) calloc(sizeof(double),m*n);
	double* x =(double*) malloc(sizeof(double)*n); 
	int i;

	for(i=0;i<nH;i++)
	{
		double* h = &H[9*i];
		int line1 = 2*n*i;
		int line2 = line1+6;

		V[line1+0]=h[0]*h[1];
		V[line1+1]=h[0]*h[4] + h[3]*h[1];
		V[line1+2]=h[3]*h[4];
		V[line1+3]=h[0]*h[7]+h[6]*h[1];
		V[line1+4]=h[3]*h[7]+h[6]*h[4];
		V[line1+5]=h[6]*h[7];


		V[line2+0]=h[0]*h[0] - h[1]*h[1];
		V[line2+1]=2*(h[0]*h[3] - h[1]*h[4]);
		V[line2+2]=h[3]*h[3] - h[4]*h[4];
		V[line2+3]=2*(h[0]*h[6] - h[1]*h[7]);
		V[line2+4]=2*(h[3]*h[6] - h[4]*h[7]);
		V[line2+5]=h[6]*h[6] - h[7]*h[7];
		if (NORMALIZE) {
			mtxNormalizeVector(6,&V[line1]);
			mtxNormalizeVector(6,&V[line2]);
		}
	}

	if(PRINT) mtxShowMat("[V]",V,m,n);
	solveAx0(m,n,V,x); /* solves the system [V]{x}={0} */
	i=0;
	B[i++]=x[0]; B[i++]=x[1]; B[i++]=x[3];
	B[i++]=x[1]; B[i++]=x[2]; B[i++]=x[4];
	B[i++]=x[3]; B[i++]=x[4]; B[i++]=x[5];

	mtxShowMat("[B]",B,3,3);
	free(x);
	free(V);
}

/**/
void calcB_UoVo_Known(int nH, double* H, double* B)
{	

	int m = 2*nH;
	int n = 4; // only 4 parameters needed to be calculated
	double* V =(double*) calloc(sizeof(double),m*n);
	double* x =(double*) malloc(sizeof(double)*n); 
	int i;

	double Uo = 0.0;
	double Vo = 0.0;
	double a,b,c,d,e,f;

	for(i=0;i<nH;i++)
	{
		double* h = &H[9*i];
		int line1 = 2*n*i;
		int line2 = line1+4;

		a =  h[0]*h[1];
		b = (h[0]*h[4] + h[3]*h[1]);
		c =  h[3]*h[4];
		d = (h[0]*h[7]+h[6]*h[1]);
		e = (h[3]*h[7]+h[6]*h[4]);
		f =  h[6]*h[7];

		V[line1+0]= a - d*Uo + f*Uo*Uo;
		V[line1+1]= b - d*Vo - e*Uo + 2.0*f*Uo*Vo;
		V[line1+2]= c - e*Vo + f*Vo*Vo;
		V[line1+3]= f;
		
		a = h[0]*h[0] - h[1]*h[1];
		b = 2*(h[0]*h[3] - h[1]*h[4]);
		c = h[3]*h[3] - h[4]*h[4];
		d = 2*(h[0]*h[6] - h[1]*h[7]);
		e = 2*(h[3]*h[6] - h[4]*h[7]);
		f = h[6]*h[6] - h[7]*h[7];

		V[line2+0]= a - d*Uo + f*Uo*Uo;
		V[line2+1]= b - d*Vo - e*Uo + 2.0*f*Uo*Vo;
		V[line2+2]= c - e*Vo + f*Vo*Vo;
		V[line2+3]= f;
		
		if (NORMALIZE) {
			mtxNormalizeVector(4,&V[line1]);
			mtxNormalizeVector(4,&V[line2]);
		}
	}

	if(PRINT) mtxShowMat("[V]",V,m,n);
	solveAx0(m,n,V,x); /* solves the system [V]{x}={0} */
	
	x[0]/=x[3];
	x[1]/=x[3];
	x[2]/=x[3];
	x[3]/=x[3];

	i=0;
	B[i++]=x[0];               B[i++]=x[1];               B[i++]= -Vo*x[1] - Uo*x[0];
	B[i++]=x[1];               B[i++]=x[2];               B[i++]= -Vo*x[2] - Uo*x[1];
	B[i++]=-Vo*x[1] - Uo*x[0]; B[i++]=-Vo*x[2] - Uo*x[1]; B[i++]= Vo*Vo*x[2] + 2.0*Vo*Uo*x[1] + Uo*Uo*x[0] + x[3];

	mtxShowMat("[B]",B,3,3);
	free(x);
	free(V);
}

/* Get intrinsic parameters to decomposition matrix B*/
int  calcA(double* B, double* A)
{  
	double alpha,betha,gamma,u0,v0,lambda;
	double den=B[0]*B[4]-B[1]*B[1];

	if  (fabs(den)< TOL ) { printf("den=B[0]*B[4]-B[1]*B[1]=%f\n",den); return 0; }
	v0 = (B[1]*B[2]-B[0]*B[5])/den;
	if (fabs(B[0])<TOL){ printf("B[0]=%f\n",B[0]); return 0; }
	lambda = B[8]-(B[2]*B[2]+v0*(B[1]*B[2]-B[0]*B[5]))/B[0];
	if (lambda/B[0]<0) { printf("lambda/B[0]=%f\n",lambda/B[0]); return 0; }
	alpha=sqrt(lambda/B[0]);
	if ((lambda*B[0]/den)<0) { printf("lambda*B[0]/den=%f\n",lambda*B[0]/den); return 0; }
	betha = sqrt(lambda*B[0]/den);
	gamma = - B[1]*alpha*alpha*betha/lambda;
	u0=gamma*v0/betha-B[2]*alpha*alpha/lambda;

	A[0]=alpha; A[1]=gamma; A[2]=u0;
	A[3]=0;     A[4]=betha; A[5]=v0;
	A[6]=0;     A[7]=0;     A[8]=1;

	return 1;
}

/* Get intrinsic parameters to decomposition matrix B*/
void correctA(double* AL, double* N, double* A)
{
	double Ninv[3*3];
	double det=m3Inv(N,Ninv);
	m3CopyAB(Ninv,N);
	m3MultAB(Ninv,AL,A);
}

/* Get intrinsic parameters processing homographies */
int Zhang2DInternal(int nHomographies, double* H, double N[9], double A[9])
{
	double B[3*3], AL[3*3];
	calcB(nHomographies,H,B);
	if(calcA(B,AL))
	{
		if (PRINT)mtxShowMat("[A']",AL,3,3);
		correctA(AL,N,A);
		return 1;
	}
	else
		printf("\nCan´t calculated [A']\n");

	return 0;
}

/* Get intrinsic parameters processing homographies and guessing initial values for (Uo,Vo)*/
int Zhang2DInternal_UoVo_Know(int nHomographies, double* H, double N[9], double A[9])
{
	double B[3*3], AL[3*3];
	calcB_UoVo_Known(nHomographies,H,B);
	if(calcA(B,AL))
	{
		if (PRINT)mtxShowMat("[A']",AL,3,3);
		correctA(AL,N,A);
		mtxShowMat("\nIntrinsic parameters [A]",A,3,3);
		return 1;
	}
	else
		printf("\nCan´t calculated [A']\n");

	return 0;
}

/* Get extrinsic parameters for each view using in calibration process */
int Zhang2DExternal(double H[9], double A[9], double K[12], double N[9])
{
	double Q[9],R[9],H_[9];
	double lambda1,lambda2,lambda3;
	double Ainv[3*3],hi[3], Ninv[3*3];
	double r1[3],r2[3],r3[3],t[3];
	int i;

	double lbda1,lbda2,lbda3; 

	m3Inv(N,Ninv);
	m3Inv(A,Ainv);
		
	//				m´ =	[H´]*M
	//				m´ =	[N]*m and [H´] = [N]*[H]*M
	// [Ninv]*m´ =	[Ninv]*[H´]*M
	//				m =		[H]*M
	mtxAB(Ninv,H,3,3,3,H_);

	if (PRINT) mtxShowMat("[H_]:",H_,3,3);

	mtxCol(H_,0,3,3,hi);//h1
	m3MultAb(Ainv,hi,r1);
	if (PRINT)mtxShowMat("[r1]:",r1,1,3);
	
	//lambda1 = 1/mtxNormalizeVector(3,r1);
	lbda1 = mtxNormalizeVector(3,r1);
	lambda1 = 1/lbda1;

	if (PRINT)printf("lambda %lf",lambda1);
	mtxCol(H_,1,3,3,hi);//h2
	m3MultAb(Ainv,hi,r2);
	if (PRINT)mtxShowMat("[r2]:",r2,1,3);
	
	//lambda2 = 1/mtxNormalizeVector(3,r2);
	lbda2 = mtxNormalizeVector(3,r2);
	lambda2 = 1/lbda2;
	
	//mtxScaVec(r2,3,1/lambda1,r2);
	m3Cross(r1,r2,r3);
	if (PRINT)mtxShowMat("[r3]:",r3,1,3);

	mtxCol(H_,2,3,3,hi);//h3
	m3MultAb(Ainv,hi,t);

	lbda3 = (lbda1 + lbda2)/2;
	lbda3 = 1./lbda3;

	lambda3 = (lambda1+lambda2)/2;
	mtxScaVec(t,3,lambda3,t);

	i=0;
	Q[i++]=r1[0]; Q[i++]=r2[0]; Q[i++]=r3[0];
	Q[i++]=r1[1]; Q[i++]=r2[1]; Q[i++]=r3[1];
	Q[i++]=r1[2]; Q[i++]=r2[2]; Q[i++]=r3[2];

	if (PRINT)mtxShowMat("\nOrientation initial [Q]",Q,3,3);

	// This function don´t work well
	//BestRotationMatrix(Q,R);
	for(i=0;i<9;i++){
		R[i]= Q[i];
	};
	
	if (PRINT)mtxShowMat("\nOrientation normalize [R]:",R,3,3);

	i=0;
	K[i++]=R[0]; K[i++]=R[1]; K[i++]=R[2]; K[i++]=t[0];
	K[i++]=R[3]; K[i++]=R[4]; K[i++]=R[5]; K[i++]=t[1];
	K[i++]=R[6]; K[i++]=R[7]; K[i++]=R[8]; K[i++]=t[2];

	//TODO falta minimizacao

	if (PRINT) mtxShowMat("\nExtrinsic parameters externos [K]:",K,3,4);
	return 1;
}

/* Get distortion paramters */
void Zhang2D_get_distortion(int nPts, double * imgPts, double * imgPtsNormalized, double* imgPtsProjected, double A[9], double* k){

	int i,line1,line2;
	double Uo,Vo,u_uo,v_vo,x2_y2;
	double *_imgPts_uv,*_imgPts_xy_Norm,*_imgPts_uv_proj;
	double* D =(double*) malloc(sizeof(double)*((nPts)*2 *2));
	double* d =(double*) malloc(sizeof(double)*((nPts)*2)); 
	
	
	Uo =  A[2];  Vo =  A[5]; 

	for (i=0;i<nPts;i++) {
		line1 = 2*2*i;
		line2 = line1+2;
		_imgPts_uv = &imgPts[i*2];
		_imgPts_xy_Norm = &imgPtsNormalized[i*2];
		_imgPts_uv_proj = &imgPtsProjected[i*2];

		x2_y2 = (_imgPts_xy_Norm[0]*_imgPts_xy_Norm[0] + _imgPts_xy_Norm[1]*_imgPts_xy_Norm[1]);
		u_uo = _imgPts_uv_proj[0] - Uo;
		v_vo = _imgPts_uv_proj[1] - Vo;
		D[line1+0]= (u_uo)*(x2_y2); 
		D[line1+1]= (u_uo)*(x2_y2)*(x2_y2); 
		
		D[line2+0]= (v_vo)*(x2_y2); 
		D[line2+1]= (v_vo)*(x2_y2)*(x2_y2);

		d[2*i] = _imgPts_uv[0] - _imgPts_uv_proj[0];
		d[2*i+1] = _imgPts_uv[1] - _imgPts_uv_proj[1];
	};

	mtxAxb(D,((nPts)*2), 2, d, k);
	free(D);
	free(d);

}

