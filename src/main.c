#include <stdio.h>
//#include <conio.h>
#include <ncurses.h>
#include <math.h>
 
#include "zhang2dcalib.h"
#include "matrix.h"

#define PRINT 0
#define MAX_POINTS 512
#define MAX_IMAGES 5

void saveIntrinsicData(char* name, double* A)
{ 
	FILE* fpcsv=fopen(name,"wt");

	fprintf(fpcsv,"Intrinsic parameters:\n");
	fprintf(fpcsv,"	%f \t	%f	\t %f \n",A[0],A[1],A[2]);
	fprintf(fpcsv,"	%f \t	%f	\t %f \n",A[3],A[4],A[5]);
	fprintf(fpcsv,"	%f \t	%f	\t %f \n",A[6],A[7],A[8]);      
	fclose(fpcsv);

}

void saveData(char* name, int n, double* modelPoints, double* imagePoints)
{
   FILE* fpcsv=fopen(name,"wt");
   int i;

   fprintf(fpcsv,"n,xm,ym,zm,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5\n");
   for (i=0; i<n; i++ ) {
         fprintf(fpcsv,"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",i,
         modelPoints[3*i],modelPoints[3*i+1],modelPoints[3*i+2],
         imagePoints[0*2*n+2*i],imagePoints[0*2*n+2*i+1], 
         imagePoints[1*2*n+2*i],imagePoints[1*2*n+2*i+1], 
         imagePoints[2*2*n+2*i],imagePoints[2*2*n+2*i+1], 
         imagePoints[3*2*n+2*i],imagePoints[3*2*n+2*i+1], 
         imagePoints[4*2*n+2*i],imagePoints[4*2*n+2*i+1] );
      }
      fclose(fpcsv);
   }

int loadModel(double* modelPoints, double* imagePoints)
{
   int i,n=0;
   FILE* fpm  = fopen("model.txt","rt");
   FILE* fpi1 = fopen("data1.txt","rt");
   FILE* fpi2 = fopen("data2.txt","rt");
   FILE* fpi3 = fopen("data3.txt","rt");
   FILE* fpi4 = fopen("data4.txt","rt");
   FILE* fpi5 = fopen("data5.txt","rt");

   if (fpi1==NULL ||fpi2==NULL ||fpi3==NULL ||fpi4==NULL ||fpi5==NULL || fpm==NULL) { printf("Arq error\n"); return 1; }

   for (n=0; !feof(fpm); n++ ) {
      fscanf(fpm,"%lf %lf ",&modelPoints[3*n],&modelPoints[3*n+1]);
      modelPoints[3*n+2]=1;
   }

   fclose(fpm); 

   for (i=0; i<n; i++ ) {
      fscanf(fpi1,"%lf %lf ",&imagePoints[0*2*n+2*i],&imagePoints[0*2*n+2*i+1]);
      fscanf(fpi2,"%lf %lf ",&imagePoints[1*2*n+2*i],&imagePoints[1*2*n+2*i+1]);
      fscanf(fpi3,"%lf %lf ",&imagePoints[2*2*n+2*i],&imagePoints[2*2*n+2*i+1]);
      fscanf(fpi4,"%lf %lf ",&imagePoints[3*2*n+2*i],&imagePoints[3*2*n+2*i+1]);
      fscanf(fpi5,"%lf %lf ",&imagePoints[4*2*n+2*i],&imagePoints[4*2*n+2*i+1]);

   }
   fclose(fpi1); fclose(fpi2); fclose(fpi3); fclose(fpi4); fclose(fpi5);
   saveData("original.csv",n,modelPoints,imagePoints);
   return n;
}

double erroSqrt(double* photo,double* proj)
{
	double X,Y;
	double xProj,yProj;
	xProj = proj[0]/proj[2];
	yProj = proj[1]/proj[2];

	X = photo[0] - xProj;
	Y = photo[1] - yProj;

	return sqrt(X*X + Y*Y);
}


void   erro_ARt_calculated (int m,	double* A,double* K,  
														double* err, double* errAVG,  
														double* GB_modelPoints, 
														double* GB_photoPoints, 
														double* GB_cameraPoints,
														double* GB_photoPoints_proj)
														/* m : number of points to fit */
														/* err :  vector of error from data */
{
	int    k;	
	double res[3],model[4], camPt[3];
	double ARt[12];

	mtxAB(A,K,3,3,4,ARt);	
	if (PRINT)mtxShowMat("\n[ARt] ----\n",ARt,3,4);

	for(k=0; k< m; k++) {

		model[0] = GB_modelPoints[3*k];
		model[1] = GB_modelPoints[(3*k)+1];
		model[2] = 0.0;
		model[3] = 1.0;

		mtxAb(K,model,3,4,camPt);	
		mtxAb(ARt,model,3,4,res);	

		GB_cameraPoints[3*k]		 = camPt[0];
		GB_cameraPoints[(3*k)+1] = camPt[1];
		GB_cameraPoints[(3*k)+2] = camPt[2];

		GB_photoPoints_proj[2*k] = res[0]/res[2];
		GB_photoPoints_proj[(2*k)+1] = res[1]/res[2];

		err[k] =erroSqrt(&(GB_photoPoints[2*k]),res);
		*errAVG+=err[k];
	}
	*errAVG/=m;
}



int main( )
{
   double modelPoints[3*MAX_POINTS];
	 double modelPointsCopy[3*MAX_POINTS];
   double imagePoints[2*MAX_POINTS*MAX_IMAGES];
	 double imagePointsCopy[2*MAX_POINTS*MAX_IMAGES];
	 double imagePointsProjected[2*MAX_POINTS*MAX_IMAGES];
	 double cameraPoints[3*MAX_POINTS*MAX_IMAGES];

   double H[5*3*3],A[3*3], N[3*3],K[5*3*4];
   int i,j,k,n;
   int pts[4]={3,31,255,227};
   double X[4],Y[4];
   double u[4],v[4];
   double Ncopy[3*3];
	 double K_distortion[2] = {0.0,0.0};

	 double error_AVG = 0.0;
	 double ErrorTotal= 0.0;
	 double err[MAX_POINTS];

	 /* Auxiliary pointers */
		double* h, *imPts, *_Hpt, *_Kpt, *camCoord, *imgPtsProj, * imgPtsCpy;

		/* This values have been copied of original site of Zhengyou Zhang and was used here 
			 to show the accuracy of our implementation of calibration 2D using a planar pattern 
		*/ 
		double Acalib[9]={832.5, 0.204494, 303.959, 0.0, 832.53, 206.585,0.0,0.0,1.0};
		double Kcalib[5*12]={
			0.992759 ,-0.026319, 0.117201, -3.84019,
				0.0139247, 0.994339, 0.105341, 3.65164,
				-0.11931, -0.102947, 0.987505, 12.791,

				0.997397,   -0.00482564, 0.0719419, -3.71693, 
				0.0175608,   0.983971,   -0.17746,   3.76928, 
				-0.0699324,  0.178262,   0.981495,  13.1974,

				0.915213, -0.0356648, 0.401389, -2.94409, 
				-0.00807547, 0.994252, 0.106756, 3.77653, 
				-0.402889, -0.100946, 0.909665, 14.2456,

				0.986617, -0.0175461, -0.16211, -3.40697, 
				0.0337573, 0.994634, 0.0977953, 3.6362, 
				0.159524, -0.101959, 0.981915, 12.4551,

				0.967585, -0.196899, -0.158144, -4.07238, 
				0.191542, 0.980281, -0.0485827, 3.21033, 
				0.164592, 0.0167167, 0.98622, 14.3441
		};



   n=loadModel(modelPoints,imagePoints);   
	 
   for (j=0;j<4;j++) { k=pts[j]; X[j]=modelPoints[3*k]; Y[j]=modelPoints[3*k+1]; }
	 printf("\n ModelPts (%.2f,%.2f) (%.2f,%.2f) (%.2f,%.2f) (%.2f,%.2f) \n",X[0],Y[0],X[1],Y[1],X[2],Y[2],X[3],Y[3]);
	 
	 /* Get a copy of image points */
	 for ( j = 0; j < n*5; j++ ) { 
		 imagePointsCopy[2*j]   = imagePoints[2*j]; 
		 imagePointsCopy[2*j+1] = imagePoints[2*j+1];
	}

	 for ( j = 0; j < n; j++ ) { 
		 modelPointsCopy[3*j]   = modelPoints[3*j]; 
		 modelPointsCopy[3*j+1] = modelPoints[3*j+1]; 
		 modelPointsCopy[3*j+2] = modelPoints[3*j+2];
	 }

   normalizeImagePoints(640,480,n,5,imagePoints,N);
	 for (i=0;i<9;i++) Ncopy[i] = N[i]; 
   saveData("normalized.csv",n,modelPoints,imagePoints);
	 
   for (i=0;i<5;i++) {
      h     = &H[9*i];        
			imPts = &imagePoints[i*2*n];
      for (j=0;j<4;j++) 
			{ k=pts[j]; u[j]=imPts[2*k]; v[j]=imPts[2*k+1]; }
			printf("\n ImagePts (%.2f,%.2f) (%.2f,%.2f) (%.2f,%.2f) (%.2f,%.2f) \n",u[0],v[0],u[1],v[1],u[2],v[2],u[3],v[3]);      
      homography(n,modelPoints,imPts,h);
      mtxScaMat(h,  1./h[8]  ,3,3,h  );      
      printf("\nHomography = %d ",i);
      mtxShowMat("\n[H]",h,3,3);
      
   }
   printf("\n\n");

	 printf("\n****************** TEST WITH OUR IMPLEMENTATION *****************\n");
	 if(Zhang2DInternal(5, H, N, A)){
		 mtxShowMat("\nIntrinsic parameters [A]:",A,3,3);
		 for (i=0; i<5; i++){
			 _Hpt     = &H[9*i];
			 _Kpt     = &K[12*i];
			 camCoord = &cameraPoints[3*n*i]; 
			 imgPtsProj = &imagePointsProjected[2*n*i] ;
			 Zhang2DExternal(_Hpt, A,_Kpt, Ncopy);
			 
			 printf("\n");
			 printf("\nViews %d",i);
			 mtxShowMat("\nExtrinsic parameters [R |t]:",_Kpt,3,4);
			 error_AVG = 0.0;

			 imgPtsCpy = &imagePointsCopy[i*2*n];
			 
			 /* We calculated the error of re-reprojection of model points inside each pattern view 
					used to calibrate our camera, for that we used matrixs that have been calculated like
				  matrix of Intrinsic parameters [A] and Extrinsic parameters [K], remember that we calculated 
					one matrix of extrinsic parameters [Ki] for each view "i" of the planar pattern. */
			 erro_ARt_calculated (n,	A, _Kpt, err, 
				 &error_AVG, modelPointsCopy, 
				 imgPtsCpy, 
				 camCoord, 
				 imgPtsProj);
			 printf("\nAverage error of reprojection points: %.5lf",error_AVG);
			 ErrorTotal+=error_AVG;
			};
	 }

	 /* Save values of Intrinsic parameters calculated with our implementation*/
	 saveIntrinsicData("Zhang2D_calibration.txt", A);

	 printf("\n*************************************************************");
	 printf("\nTotal error with our implementation : %.5lf",ErrorTotal);
	 printf("\n*************************************************************");

	 Zhang2D_get_distortion(n*5, imagePointsCopy, imagePoints, 
		 imagePointsProjected, A, K_distortion);
	 mtxShowMat("\nDistortion parameters: ",K_distortion,1,2);

	 /* Test error of re-projection model points to image points using values in calib.txt that we found in 
			Zhengyou Zhang microsoft site */
	 printf("\n\n\n*************** TEST WITH PARAMETERS FOUND IN calib.txt ***************\n");
	 ErrorTotal = 0.0;
	 mtxShowMat("\nCalib.txt Intrinsic parameters [A]:",Acalib,3,3);
	 for (i=0; i<5; i++){
			 _Kpt     = &Kcalib[12*i];
			 camCoord = &cameraPoints[3*n*i]; 
			 imgPtsProj = &imagePointsProjected[2*n*i] ;
			 imgPtsCpy = &imagePointsCopy[i*2*n];

			 printf("\n");
			 printf("\nViews %d",i);
			 mtxShowMat("\nCalib.txt Extrinsic parameters [R |t]:",_Kpt,3,4);
			 
			 error_AVG = 0.0;
			 
			 erro_ARt_calculated (n,	Acalib, _Kpt, err, 
				 &error_AVG, modelPointsCopy, 
				 imgPtsCpy, 
				 camCoord, 
				 imgPtsProj);
			 printf("\nCalib.txt Average error of reprojection points: %.5lf",error_AVG);
			 ErrorTotal+=error_AVG;

	};
	 printf("\n*************************************************************");
	 printf("\nTotal error with parameters in calib.txt: %.5lf",ErrorTotal);
	 printf("\n*************************************************************");
   printf("\n\n");
	 //getch();

   return 0;
}





