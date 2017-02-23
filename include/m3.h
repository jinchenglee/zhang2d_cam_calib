#ifndef _M3_H
#define _M3_H


/**
 * @file m3.h
 * TAD 3x3 Matrix		
 * @author Marcelo Gattass	 
 * @date Jul10,2006
 */

#ifdef __cplusplus
extern "C" {
#endif

	/**
	* computes the determinant of the matrix [A]
	* @param A [A] matrix
	* @return determinant  of the matrix [A]	     
	*/
	double m3Det(double* A);                           

	/**
	* computes the trace of the matrix [A]
	* @param[in] A matrix [A] 
	* @return trace  of the matrix [A]	     
	*/
	double m3Trace(double* A);                            


	/**
	* computes the inverse,[Ainv], of the matrix [A]
	* @param[in] A matrix [A] 
	* @param[out] Ainv matrix [Ainv]
	* @return determinant  of the matrix [A]
	*/
	double m3Inv( double* A, double* Ainv );


	/**
	* computes the cross product {c} = {a}x{b}
	* @param[in] a vector {a}
	* @param[in] b vector {b}
	* @param[out] c vector {c}     
	*/
	void m3Cross(double *a,double *b,double *c);


	/**
	* computes the matrix product [AB]=[A][B]
	* @param[in] A matrix [A]
	* @param[in] B matrix [B]
	* @param[out] AB matrix [AB]
	*/
	void m3MultAB(double* A, double* B, double* AB); 


	/**
	* computes the transformation {x}=[A]{b}
	* @param[in] A matrix [A]
	* @param[in] b vector {b}
	* @param[out] x vector {x}
	*/
	void m3MultAb(double* A, double* b, double* x);   

	/**
	* computes the solution of [A]{x}={b}            
	* @param[in] A matrix [A]
	* @param[in] b vector {b}
	* @param[out] x vector {x}
	*/
	double m3SolvAxb(double* A, double* b, double* x);

	/**
	* copy [A] to [B]
	* @param[in] A matrix [A]
	* @param[out] B matrix [B]
	*/
	void m3CopyAB(double* A, double* B);

	/**
	* print in the console the matrix [A]
	* @param[in] text text
	* @param[in] A matrix [A]
	*/
	void m3PrintMat(char* text,double* A);

	/**
	* print in the console the vector {v} 
	* @param[in] text text
	* @param[in] v vector {v}
	*/
	void m3PrintVet(char* text,double* v);  

#ifdef __cplusplus
}
#endif

#endif 