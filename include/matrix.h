#ifndef _MATRIX_H_
#define _MATRIX_H_

/**
 * @file matrix.h
 * Simple Math		
 *					
 * @author Marcelo Gattass
 *
 * @author Manuel E. L. Fernandez	 
 *
 * @date Jul06,2006
 *
 * Vectors are represented in ANSI C as double*
 * Matrices are stored as vectors in a line by line fashion.
 * i.e., the ij element of a matrix Amxn (m rowns and n coluns) is 
 * stored in the i*n+j  position of the vector that stores the 
 * matrix.
 *
 * The ideia is to enable simple matrix creations such as:
 * double a[2*3] = {  1., 2., 7.,
 *                    -1., 5., 4.};
 * This ANSI C statment declares and stores the 2x3 matrix:\n
 *     |  1. 2. 7. |\n
 *     | -1  5. 4. |\n
 */

#if	defined(__cplusplus)
extern "C" {
#endif


	/**
	* Computes the solution for the linear equations [a]{x}={b} using Gauss method.
	* @param[in] a Coefficient matrix (stored as a vector in a line by line mode).
	* @param[in] n number of equations.
	* @param[in,out] b input = Constant vector, output = Solution {x}.
	* @return 0 => no solution or infinity; 1=> one solution.
	*/
	int mtxGaussAxb (double* a, int n, double* b);

	/**
	* Computes the LU decomposition of a matrix [a]nxn, [A]=[L][U].
	* @param[in] a input = a nxn square matrix stored as a vector in a line by line mode, output = represented as [L] and [U} in the same matrix.
	* @param[in] n dimension of [a].
	* @param[in] scale scale.
	* @param[out] rp row permutation.
	* @param[out] d +1 ,number of row interchanges is even (-1 = odd).
	*/
	int mtxDecompLU ( double* a, int n, int *rp, double* d, double* scale );

	/**
	* Solves the set of n linear equations [A]{x} = {b} where the matrix for this function is de LU decomposition determined by the function mtxDecompLU.
	* @param[in] a [a] is not as the matrix A but rather as its LU decomposition.
	* @param[in] n dimension of [a].
	* @param[in] rp {rp} = row permutation (determined by the function mtxDecompLU).
	* @param[in,out] b {b} is the right hand side vector {b}.
	*/
	void mtxBackSubLU( double* a, int n, int* rp, double* b );

	/**
	* Computes the determinant of a matrix [A] where the matrix for this function is de LU decomposition determined by the function mtxDecompLU.
	* @param[in] a is not as the matrix A but rather as its LU decomposition.
	* @param[in] d +1 ,number of row permutation in the decomposition waseven (-1 = odd).
	* @param[in] n dimension of [a].
	* @return determinant of [a].
	*/
	double mtxDetLU( double* a, double d, int n);

	/**
	* Computes the Singlar Value Decomposition of a matrix [a], i.e.: [a]=[u][d][v]T
	* @param[in] a a mxn matrix (stored as a vector in a line by line mode).
	* @param[in] m number of rows of matrix [a].
	* @param[in] n number of coluns of matrix [a].
	* @param[out] u a mxn matrix (columns are orthogonal vectors)(stored as a vector in a line by line mode).
	* @param[out] d 1xn matrix that contains the elements of the diagonal matrix [d](nxn).
	* @param[out] v a nxn orthogonal matrix(stored as a vector in a line by line mode).
	* @param[out] tmp temporary 1xn vector (space used by the function).
	*/
	int mtxSVD(double* a, int m, int n, double* u, double* d, double* v, double* tmp);

	/**
	* Finds a non trivial solution for the system Ax=0. A mxn, m>n and rank(A)< n (better if = n-1).
	* @param[in] a a mxn matrix (stored as a vector in a line by line mode).
	* @param[in] m number of rows of matrix [a].
	* @param[in] n number of coluns of matrix [a].
	* @param[out] x 1xm solution vector.
	* @param[out] u a mxn matrix (columns are orthogonal vectors)(stored as a vector in a line by line mode).
	* @param[out] d 1xn matrix that contains the elements of the diagonal matrix [d](nxn).
	* @param[out] v a nxn orthogonal matrix(stored as a vector in a line by line mode).
	* @param[out] tmp a temporary 1xm vector (space used by the function).
	* @return  The return value indicates the error. He is the ratio from the smalest to the largest singular value of A. If rank(A) < n this value should be zero (~=TOL for the implementation).
	*/
	double mtxSVDAx0(double* a, int m, int n, double* x, double* u, double* d, double* v, double* tmp);

	/**
	* Finds a solution for the system Ax=b, A mxn, m>n. 
	* @param[in] a a mxn matrix (stored as a vector in a line by line mode).
	* @param[in] m number of rows of matrix [a].
	* @param[in] n number of coluns of matrix [a].
	* @param[out] x 1xm solution vector.
	* @param[in] b 1xm right hand side vector.
	* @param[out] u a mxn matrix (columns are orthogonal vectors)(stored as a vector in a line by line mode).
	* @param[out] d 1xn matrix that contains the elements of the diagonal matrix [d](nxn).
	* @param[out] v a nxn orthogonal matrix(stored as a vector in a line by line mode).
	* @param[out] tmp a temporary 1xm vector (space used by the function)
	*/
	void mtxSVDAxb(double* a, int m, int n, double* x, double* b, double* u, double* d, double* v, double* tmp);

	/**
	* Add the tensor product of the vector {v} (i.e., {v}{v}T) to the matrix [A] -> [A]+={v}{v}T
	* @param[in,out] a matrix [A]
	* @param[in] v vector {v}
	* @param[in] n number of elements of v
	*/
	void mtxAddMatVecTensor(double* a, double* v, int n);

	/**
	* computes the matrix product {x}=[A]{b}
	* @param[in] a matrix [A]
	* @param[in] b vector {b}
	* @param[in] m number of lines of [A] and number of elements of {b}
	* @param[in] n number of columns of [A]
	* @param[out] x vector {x}=[A]{b}, {x}=n
	*/
	void mtxAb(double* a, double* b, int m, int n, double* x);

	/**
	* computes the matrix product {x}=[A]T{b}
	* @param[in] a matrix [A]
	* @param[in] b vector {b}
	* @param[in] m number of lines of [A] and number of elements of {b}
	* @param[in] n number of columns of [A]
	* @param[out] x vector {x}=[A]T{b}, {x}=n
	*/
	void mtxAtb(double* a, double* b, int m, int n, double* x);

	/**
	* computes the matrix product [X]=[A][B]
	* @param[in] a matrix [A], [A]=mxp
	* @param[in] b matrix [B], [B]=pxn
	* @param[in] m number of lines of [A]
	* @param[in] p number of columns of [A] and number of lines of [B]
	* @param[in] n number of columns of [B]
	* @param[out] x matrix [X]=[A][B], [X]=mxn
	*/
	void mtxAB(double* a, double* b, int m, int p, int n, double* x);

	/**
	* computes the matrix product [X]=[A][B]T
	* @param[in] a matrix [A], [A]=mxp
	* @param[in] b matrix [B], [B]=nxp
	* @param[in] m number of lines of [A]
	* @param[in] p number of columns of [A] and [B]
	* @param[in] n number of lines of [B]
	* @param[out] x matrix [X]=[A][B]T, [X]=mxn
	*/
	void mtxABt(double* a, double* b, int m, int p, int n, double* x);

	/**
	* computes the matrix product [X]=[A]T[B]
	* @param[in] a matrix [A], [A]=mxp
	* @param[in] b matrix [B], [B]=mxn
	* @param[in] m number of lines of [A] and [B]
	* @param[in] p number of columns of [A]
	* @param[in] n number of columns of [B]
	* @param[out] x matrix [X]=[A]T[B], [X]=pxn
	*/
	void mtxAtB(double* a, double* b, int m, int p, int n, double* x);

	/**
	* computes the sum [A]+[B]
	* @param[in] a matrix [A]
	* @param[in] b matrix [B]
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] x matrix [X]=[A]+[B]
	*/
	void mtxAddMat(double* a, double* b, int m, int n, double* x);

	/**
	* computes the subtration [A]-[B]
	* @param[in] a matrix [A]
	* @param[in] b matrix [B]
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] x matrix [X]=[A]-[B]
	*/
	void mtxSubMat(double* a, double* b, int m, int n, double* x);

	/**
	* multiplies a matrix by a scalar
	* @param[in] a matrix
	* @param[in] s scalar 
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] x matrix [X]=s[A]
	*/
	void mtxScaMat(double* a, double s,int m, int n, double* x);

	/**
	* computes the transpose [X]=[A]T
	* @param[in] a matrix 
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] x transpose of [A]
	*/
	void mtxAt(double* a, int m, int n, double* x);

	/**
	* copy [A] to [X]
	* @param[in] a matrix 
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] x copy of matrix [A]
	*/
	void mtxMatCopy(double* a, int m, int n, double* x);

	/**
	* computes the sum {x}={v}+{u}
	* @param[in] u vector 
	* @param[in] v vector
	* @param[in] m dimension of u and v
	* @param[out] x new vector {x}={v}+{u}
	*/
	void mtxAddVec(double* u, double* v, int m, double* x);

	/**
	* computes the Subtraction {x}={v}-{u}
	* @param[in] u vector 
	* @param[in] v vector
	* @param[in] m dimension of u and v
	* @param[out] x new vector {x}={v}-{u}
	*/
	void mtxSubVec(double* u, double* v, int m, double* x);

	/**
	* multiplies a vector by a scalar -> {x}=s{u}  Dimensions: {u}=m and {x}=m
	* @param[in] u vector
	* @param[in] m dimension of u
	* @param[in] s scalar
	* @param[out] x new vector {x}=s{u}
	*/
	void mtxScaVec(double* u, int m,  double s, double* x);

	/**
	* dot product of vectors -> s={v}.{u}  Dimensions: {v}=m and  {u}=m
	* @param[in] u vector
	* @param[in] v vector
	* @param[in] m dimension of u and v
	* @return dot product of vectors
	*/
	double mtxDotProd(double* u, double* v, int m);

	/**
	* vectorial product of vectors with dimension 3 -> x={u}x{v}  Dimensions: {v}=3 and  {u}=3
	* @param[in] u vector
	* @param[in] v vector
	* @param[out] x vector
	* @return vectorial product of vectors
	*/
	void mtxProdVec3(double* u, double* v, double*x);

	/**
	* normalize {v} and returns the norm
	* @param[in] n number of elements
	* @param[in] v vector to be normalized
	* @return norm of v
	*/
	double mtxNormalizeVector(int n, double* v);

	/**
	* copies a vector -> {x}={v} Dimensions: {v}={x}=m
	* @param[in] v vector
	* @param[in] m number of elements
	* @param[out] x copy of v
	*/
	void mtxVecCopy(double* v, int m, double* x);

	/**
	* gets a column -> {x}=col(A)
	* @param[in] A matrix
	* @param[in] col number of the desired column  
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] x desired column
	*/
	void mtxCol(double* A,int col,  int m, int n, double* x);

	/**
	* print in the console the matrix
	* @param[in] title text
	* @param[in] a matrix to be impressed
	* @param[in] m number of lines
	* @param[in] n number of columns
	*/
	void mtxShowMat(char*  title, double* a, int m, int n);

	/**
	* print in the console the vector
	* @param[in] title text
	* @param[in] v vector to be impressed
	* @param[in] n number of elements
	*/
	void mtxShowVec(char*  title, double* v, int n);

	/**
	* Computes the nxn matrix [X] such that [A][X]=[I]
	* , then [X] is the inverse of [A].
	*(uses LU decomposition).
	* @param[in] a matrix [A]
	* @param[in] n number of columns
	* @param[out] x inverse of [A]
	*/ 
	void mtxAxInxn(double *a, int n, double *x);

	/**
	* Computes the mxn matrix [X] such that [A][X]=[I]
	* , then [X] is the inverse of [A].
	* @param[in] A matrix [A]
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[out] X inverse of [A]
	*/
	void mtxAxImxn(double *A, int m, int n, double *X);

	/**
	* @brief find solution  x to the system  Ax=b	
	*
	*						[A]x = b
	*				    [At][A]x = [At]b
	*	 inv([At][A])	[At][A]x =	inv([At][A])[At]b	
	*						[I]x =  inv([At][A])[At]b	
	* @param[in] A matrix nxn [A]
	* @param[in] m number of lines
	* @param[in] n number of columns
	* @param[in] b vector
	* @param[out] x solution
	*/ 
	void mtxAxb(double *A, int m, int n, double *b, double *x);

#if defined(__cplusplus)
}
#endif

#endif

