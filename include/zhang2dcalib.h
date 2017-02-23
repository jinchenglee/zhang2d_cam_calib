#ifndef _ZHANG2DCALIB_H
#define _ZHANG2DCALIB_H

/** 	
 * @file zhang2dcalib.h
 * Zhang's camera calibration 2D	
 * @author Lucas Teixeira
 * @author Marcelo Gattass
 * @author Manuel E. L. Fernandez	
 * @date july 06,2006
 */

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	 * @brief Computes the homography [H]3x3 given nPoints pairs of Model \n
	 * and Image Points, i.e, computes [H]                         
	 *
	 * Model point k is obtained in:                               \n
	 *   double X=modelPoints[3*k+0];   X coord of model point k   \n
	 *   double Y=modelPoints[3*k+1];   Y coord of model point k   \n
	 *   double W=modelPoints[3*k+2];   W coord of model point k   \n
	 * Image point k is obtained in:                               \n
	 *   double u=imagePoints[2*k+0];   u coord of image point k   \n
	 *   double v=imagePoints[2*k+1];   v coord of image point k   \n
	 * @param[in] nPoints number of points pairs
	 * @param[in] modelPoints points model array 
	 * @param[in] imagePoints points image array
	 * @param[out] H homography
	 * @return 0 = no-error or 1 = error                           
	 */
int homography(int nPoints, double* modelPoints, double* imagePoints, double* H);

	/** 
	 * @brief Computes the projective transformation of the point (X,Y,W) according     
	 * with the homography [H]                                                  \n
	 *
	 * @param[in] H homography
	 * @param[in] X x model points
	 * @param[in] Y y model points
	 * @param[in] W w model points
	 * @param[out] u u image points
	 * @param[out] v v image points
	 * @return the w componente
	 */
double hTransform(double* H, double X, double Y, double W, double* u, double* v);

	/**
	 * @brief Normalizes	the	points in the image.  I.e.,	maps from [0,w]x[0,h]->[-1,1]x[-1,1].
	 * 
	 * Also computes	the	inverse	of this	transformation to be obtain	the	internal 
	 * Camera matrix	in the original	image system.
	 * \n
	 * imagePoints[2*i+0] - u coordinate	of image point i.
	 * imagePoints[2*i+1] - v coordinate	of image point i.
	 * @param[in] w width image resolutions.
	 * @param[in] h height image	resolutions.
	 * @param[in] n number of image points.
	 * @param[in] views number of views .	 
	 * @param[in] imagePoints points
	 * @param[out] N	 matrix	that transforms	[-1,1]x[-1,1]->[0,w]x[0,h] stored by rows.
	 */
void normalizeImagePoints(double w, double h, int n,int views,double* imagePoints, double* N);

	/**
	 * @brief Computes the internal camera Matrix [A] (see eq. 1 from Zhang 1998 paper) 
	 * from a set of homographies.
	 *
	 * @param[in] nHomographies - Number of Homographies.
	 * @param[in] H  homographies.
	 * @param[in] N - matrix that transforms [-1,1]x[-1,1]->[0,w]x[0,h] stored by rows.
	 * 
	 * @param[out] A - 3x3 matrix of internal camera parameters stored by rows.
	 * @return 0 (zero) means fail to compute [A] and {kt]. 1 (one) means ok. 
	 */
int Zhang2DInternal(int nHomographies, double* H, double N[9], double A[9]);

	/**
	 * @brief Computes the internal camera Matrix [A] (see eq. 1 from Zhang 1998 paper) 
	 * from a set of homographies.
	 *
	 * @param[in] nHomographies - Number of Homographies.
	 * @param[in] H  homographies.
	 * @param[in] N - matrix that transforms [-1,1]x[-1,1]->[0,w]x[0,h] stored by rows.
	 * 
	 * @param[in,out] A - 3x3 matrix of internal camera parameters stored by rows.
	 * @return 0 (zero) means fail to compute [A] and {kt]. 1 (one) means ok. 
	 */
int Zhang2DInternal_UoVo_Know(int nHomographies, double* H, double N[9], double A[9]);

	/**
	 * @brief Computes the external camera matrix [Kt] given [A] and a  homography.
	 *
	 * @param[in] H 3x3 homography matrix stores row by row.
	 * @param[in] A  3x3 matrix of internal camera parameters stored by rows.
	 * @param[in] N - matrix that transforms [-1,1]x[-1,1]->[0,w]x[0,h] stored by rows.
	 * @param[out] K 3x4 matrix of external camera parameters stored by rows.
	 * @return 0 (zero) means fail to compute [A] and {kt]. 1 (one) means ok. 
	 */
int Zhang2DExternal(double H[9],double A[9], double K[12], double N[9]);

	/**
	 * Computes the radial distortion camera parameters [ k1, k2 ] (see eq. 1 from Zhang 1998 paper) 
	 * from a set of image and image reprojected points.
	 *
	 * @param[in] nPts Number of points
	 * @param[in] imgPts  Image coordinates of points 
	 * @param[in] imgPtsNormalized Image coordinates normalized
	 * @param[in] imgPtsProjected Image coordinates of points after reprojected model points using Intrinsic and Extrinsic matrix found
	 * @param[in] A Intrinsic parameters matrix
	 * @param[out] k 2x1 matrix of radial distortion parameters.
	 */
void Zhang2D_get_distortion(int nPts, double * imgPts, double * imgPtsNormalized, 
														double* imgPtsProjected, double A[9], double* k);

#ifdef __cplusplus
}
#endif

#endif 
