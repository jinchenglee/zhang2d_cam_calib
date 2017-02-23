#ifndef _MINIMIZING_H
#define _MINIMIZING_H

/** 	
* @file minimizing.h
* nonlinear minimization 
*/

#ifdef __cplusplus
extern "C" {
#endif


void lmHomography( double *photoPoints,double *modelPoints,double *initial,int numPoints,double *lm);

void lmARt( double *photoPoints,double *modelPoints,double *Ai,double *Alm,double *Ri,double *Rlm,
		    double *ti,double *tlm,int numPoints);

#ifdef __cplusplus
}
#endif

#endif 