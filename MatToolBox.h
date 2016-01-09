#ifndef __MATTOOLBOX_H__
#define __MATTOOLBOX_H__

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\vectortemplate.h"
#include "F:\ASU Courses\CEE526 Finite Elements\computational tools\Library\Library\matrixtemplate.h"

const int ELEMENTSPERLINE = 5;  // # of vector/matrix elements per line
const int FW = 9;              // field width
const int NUM_ELEMENTS_PER_LINE = 4;
// ---------------------------------------------------------------
// ------------------------ vector functions ---------------------
// ---------------------------------------------------------------

template <class T>
bool Add (const CVector<T>& A, const CVector<T>& B, 
                          CVector<T>& C)
// ==================================================================
// Function: adds two vectors and stores the result in the
//           third vector C = A + B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
    // check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        return false;
    // add
    for (int i=1; i <= n; i++)
        C(i) = A(i) + B(i);
    return true;
}

template <class T>
bool Subtract (const CVector<T>& A, 
                               const CVector<T>& B, CVector<T>& C)
// ==================================================================
// Function: subtracts one vector from another and stores the result
//           in the third vector C = A - B
//    Input: vectors A and B 
//   Output: vector C
// ==================================================================
{
	// check for incompatible vectors
    int n = A.GetSize();
    if (n != B.GetSize() || n != C.GetSize())
        return false;
    // subtract
    for (int i=1;i<=n;i++)  
		C(i)=A(i)-B(i);
    return true;
}

template <class T>
bool DotProduct (const CVector<T>& A,
                                 const CVector<T>& B, T& product)
// ==================================================================
// Function: computes the dot product of two vectors such that
//           product = A dot B
//    Input: vectors A and B 
//   Output: product 
// ==================================================================
{
    int n = A.GetSize();
	// check for incompatible vectors
    if (n != B.GetSize())
        return false;
	//dot product
	product=0.0;
   	for (int i=1;i<=n;i++)
		product += A(i)*B(i);
		 return true;
}

template <class T>
bool Normalize (CVector<T>& A)
// ==================================================================
// Function: normalizes a vector
//    Input: vector A 
//   Output: normalized vector A 
// ==================================================================
{
    T sum(0);
    int n = A.GetSize();
    for (int i=1; i <= n; i++)
        sum += A(i)*A(i);
    if (sum == T(0))
        return false;
    else
    {
        sum = sqrt(sum);
        for (int i=1; i <= n; i++)
            A(i) /= sum;
    }
    return true;
}

template <class T>
void Scale (CVector<T>& A, T c)
// ==================================================================
// Function: scales a vector by a constant c such that A = c A
//    Input: vector A and constant c 
//   Output: scaled vector A
// ==================================================================
{
	int n = A.GetSize();
	//Scale
    for (int i=1;i<=n;i++)
        A(i)=c*A(i);
}

template <class T>
T MaxValue (const CVector<T>& A)
// ==================================================================
// Function: finds the largest value among all the elements in A
//    Input: vector A 
//   Output: return value is the largest element in A
// ==================================================================
{
    int n=A.GetSize();
	T fa=A(1);
    for (int i=1;i<=n;i++)    
		if(fa<A(i))
			fa=A(i);
	return fa;
}

template <class T>
T MinValue (const CVector<T>& A)
// ==================================================================
// Function: finds the smallest value among all the elements in A
//    Input: vector A 
//   Output: return value is the smallest element in A
// ==================================================================
{
	int n = A.GetSize();
	T fa = A(1);
    for (int i=1;i<=n;i++)    
		if(fa>A(i))
			fa=A(i);
	return fa;
}

template <class T>
T TwoNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the two norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{

    int n = A.GetSize();
	T fa=0.0f;
    for (int i=1;i<=n;i++)
       fa+=pow(abs(A(i)),2);
	fa=sqrt(fa);
	return fa;
}

template <class T>
T MaxNorm (const CVector<T>& A)
// ==================================================================
// Function: computes the max norm of vector A
//    Input: vector A 
//   Output: return value is the max-norm
// ==================================================================
{
	int n = A.GetSize();
	T fa=0.0f;
    for (int i=1;i<n;i++)
		if(abs(fa)<abs(A(i+1)))
	fa=abs(fa);
	return fa;
}

template <class T>
bool CrossProduct (const CVector<T>& A,
                                   const CVector<T>& B,
                                   CVector<T>& C)
// ==================================================================
// Function: computes the cross-product of two vectors and stores the
//           result in the third vector such that C = A x B
//           (3-dimensional space)
//    Input: vectors A, B and C
//   Output: vector C
// ==================================================================
{

    int n = A.GetSize();
	// check for incompatible vectors
    if (n != B.GetSize() || n != C.GetSize())
        return false;
	//Cross product of Vectors
    C(1) = A(2)*B(3) - B(2)*A(3);
	C(2) = B(1)*A(3) - A(1)*B(3);
	C(3) = A(1)*B(2) - B(1)*A(2);
    return true;
}

// ---------------------------------------------------------------
// ------------------------ matrix functions ---------------------
// ---------------------------------------------------------------

template <class T>
bool Add (const CMatrix<T>& A, const CMatrix<T>& B, 
                          CMatrix<T>& C)
// ==================================================================
// Function: adds two matrices and stores the result in the
//           third matrix C = A + B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
	// check for incompatible of Matrix for Addition
	if((A.GetColumns()!=B.GetColumns())&&(A.GetRows()!=B.GetRows())&&(A.GetColumns()!=C.GetColumns())&&(A.GetRows()!=C.GetRows()))
		return false;
	for(int i=1;i<=A.GetRows();i++)
		for(int j=1;j<=A.GetColumns();j++)
			C(i,j)=A(i,j)+B(i,j);
	return true;
}

template <class T>
bool Subtract (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: subtracts one matrix from another and stores the result
//           in the third matrix C = A - B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
	// check for incompatible of Matrix for Subtraction
	if((A.GetColumns()!=B.GetColumns())&&(A.GetRows()!=B.GetRows())&&(A.GetColumns()!=C.GetColumns())&&(A.GetRows()!=C.GetRows()))
		return false;
	for(int i=1;i<=A.GetRows();i++)
		for(int j=1;j<=A.GetColumns();j++)
			C(i,j)=A(i,j)-B(i,j);
	return true;

}

template <class T>
bool Multiply (const CMatrix<T>& A,
                               const CMatrix<T>& B, CMatrix<T>& C)
// ==================================================================
// Function: multiplies two matrices and stores the result
//           in the third matrix C = A * B
//    Input: matrices A and B 
//   Output: matrix C
// ==================================================================
{
	// check for incompatible of Matrix multiplication
	if(A.GetColumns()!=B.GetRows())
		return false;
	for(int i=1;i<=A.GetRows();i++)
	{
		for(int j=1;j<=B.GetColumns();j++)
		{
			C(i,j)=0.00;
			for(int k=1;k<=B.GetRows();k++)
				C(i,j)+=A(i,k)*B(k,j);
		}
	}
		return true;
}

template <class T>
bool Determinant (const CMatrix<T>& A, T& det)
// ==================================================================
// Function: computes the determinant of matrix A
//    Input: matrix A and variable to hold the determinant
//   Output: determinant
// ==================================================================
{
	return false;
}

template <class T>
void Scale (CMatrix<T>& A, T c)
// ==================================================================
// Function: scales all the elements of a matrix by a constant c
//           such that A = c A
//    Input: matrix A and constant c
//   Output: scaled matrix A
// ==================================================================
{
    for (int i=1;i<=A.GetRows();i++)
		for (int j=1;j<=A.GetColumns();j++)
			A(i,j)=c*A(i,j);
}

template <class T>
T MaxNorm (const CMatrix<T>& A)
// ==================================================================
// Function: computes the max norm of matrix A
//    Input: matrix A 
//   Output: return value is the max norm
// ==================================================================
{
	int n=A.GetRows();
	T fa=A(1,1);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			if(abs(fa)<abs(A(i,j)))
				fa=A(i,j);
	return abs(fa);
}
template <class T>
T TwoNorm (const CMatrix<T>& A)
// ==================================================================
// Function: computes the two norm of vector A
//    Input: vector A 
//   Output: return value is the two-norm
// ==================================================================
{
	T fa=0.0f;
    for (int i=1;i<=A.GetRows();i++)
	{
		for(int j=1;j<=A.GetColumns();j++)
			fa+=pow(A(i,j),2);
	}
	return sqrt(fa);
}

template <class T>
bool Transpose (const CMatrix<T>& A,CMatrix<T>& B)
// ==================================================================
// Function: computes the transpose of a matrix and stores the result
//           in another matrix B = A(T)
//    Input: matrices A and B
//   Output: matrix B
// ==================================================================
{
	// check for incompatible of Matrix Transpose
  if((A.GetRows()!=B.GetColumns())&&(A.GetColumns()!=B.GetRows()))
	  return false;
	for(int i=1;i<=A.GetColumns();i++)
		for(int j=1;j<=A.GetRows();j++)
			B(i,j)=A(j,i);
	return true;
}

template <class T>
bool MatMultVec (const CMatrix<T>& A,
                                 const CVector<T>& x,
                                 CVector<T>& b)
// ==================================================================
// Function: multiplies a matrix and a vector and stores the result
//           in a vector b = A * x
//    Input: vectors A and x 
//   Output: vector b
// ==================================================================
{
	//check for compatibililty of matrix and vector
   if(A.GetColumns()!= x.GetSize())
	   return false;
   for(int i=1;i<=A.GetRows();i++)
	{
		b(i) = 0;
		for(int j=1;j<=A.GetColumns();j++)
			b(i) += A(i,j)*x(j);
    }
    return true;
}

template <class T>
bool LUFactorization (CMatrix<T>& A, T TOL)
// ==================================================================
// Function: carries out LU factorization of matrix A
//           A is replaced with L and U
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{
    return false;
}

template <class T>
bool LUSolve (const CMatrix<T>& A,
                              CVector<T>& x,
                              const CVector<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and U terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
    return false;
}

template <class T>
bool AxEqb (CMatrix<T>& A,CMatrix<T>& x,CMatrix<T>& b,T TOL)
// ==================================================================
// Function: solves A x = b using Gaussian Elimination Technique
//           (this version only for one rhs vector)
//    Input: Matrices A and b
//   Output: Matrix x
//           A false return value indicates a singular A matrix.
// ==================================================================
{
    // solves A x = b
    int i, j, k, ii;
    double c;

    // number of equations to solve
    int n = A.GetRows();
    if (n != A.GetColumns() || n != x.GetRows() || n != b.GetRows() ||
        x.GetColumns() != b.GetColumns())
        return 1;

    // x initially contains b
    x = b;

    // forward elimination
    for (k=1; k <= n-1; k++)
    {
        for (i=k+1; i <= n; i++)
        {
            // singular matrix?
            if (fabs(A(k,k)) <= TOL)
            {
                return k;
            }
            c = A(i,k)/A(k,k);
            for (j=k+1; j <= n; j++)
            {
                A(i,j) -= c * A(k,j);
            }
            x(i,1) -= c * x(k,1);
        }
    }

    // back substitution
    x(n,1) /= A(n,n);

    for (ii=1; ii <= n-1; ii++)
    {
        i = n - ii;
        double sum = 0.0;
        for (j=i+1; j <= n; j++)
        {
            sum += A(i,j) * x(j,1);
        }
        x(i,1) = (x(i,1) - sum)/A(i,i);
    }
    return true;
}

template <class T>
bool LDLTFactorization (CMatrix<T>& A,T TOL)
// ==================================================================
// Function: carries out LDL(T) factorization of matrix A
//           A is replaced with L and D. A is a symmetric matrix.
//    Input: matrix A and tolerance value to detect singular A
//   Output: matrix A 
// ==================================================================
{
	int n= A.GetRows();
	CVector<T> fva(n),fvb(n);
	double da=0.00;
	// check for whether Matrix is Symmetric A=AT
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			if(A(i,j)!=A(j,i))
				return false;
	// check for whether Matrix is positive definite x*m*xt
	for(int i=1;i<=n;i++)
		fva(i)=i;
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			fvb(i)=fva(i)*A(i,j);
	for(int i=1;i<=n;i++)
		da+=fvb(i)*fva(i);
	if(da<0)
		return false;
	//LULT Factorization
	for(int i=1;i<=n;i++)
		{
			if (A(i,i) <= TOL)
				return false;
			for(int j=1;j<=i-1;j++)
				A(i,i)-=A(j,j)*pow(A(i,j),2);	//Step 1
			for (int j=i+1;j<=n;j++)		//Step 2
			{
				A(j,i)=(A(j,i)/A(i,i));			
				for(int k=1;k<=i-1;k++)
					A(j,i)-=((A(j,k)*A(k,k)*A(i,k))/A(i,i));
			}
		}	
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			A(i,j)=A(j,i);
    return true;
}

template <class T>
bool LDLTSolve (const CMatrix<T>& A,CMatrix<T>& x,const CMatrix<T>& b)
// ==================================================================
// Function: carries out forward and backward substitution so as to
//           solve A x = b. A contains L and D terms.
//    Input: matrix A, vectors x and b
//   Output: vector x 
// ==================================================================
{
	//check for compatibililty of matrix and vector
	if((A.GetColumns()!= x.GetRows())||(A.GetRows()!=b.GetRows())||(x.GetRows()!=b.GetRows()))
	  return false;
	int n = A.GetRows();
	//Forward Substitution
	x(1,1)=b(1,1);
	for (int i=2;i<=n;i++)
	{
		x(i,1)=b(i,1);
		for(int j=1;j<=i-1;j++)
			x(i,1)-=A(i,j)*x(j,1);
	}
	//Backward Substitution
	x(n,1)=x(n,1)/A(n,n);
	for (int i=(n-1);i>=1;i--)
	{
		x(i,1)=(x(i,1)/A(i,i));
		for (int j=(i+1);j<=n;j++)
			x(i,1)-=(A(i,j)*x(j,1));
	}
	return true;
}
template <class T>
void PrintMatrixRowWise (CMatrix<T>& A, const std::string& heading,
                         std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: outputs a matrix into stream Out
// Input:    matrix, a heading, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, j, k, nRows, nColumns;
    
    Out << '\n' << heading << '\n';
    Out << setiosflags (std::ios::left);

    nRows = A.GetRows(); nColumns = A.GetColumns();
    for (i=1; i <= nRows; i++)
    {
        Out << "Row No: " << i << '\n';
        for (j=1; j <= nColumns; j=j+NUM_ELEMENTS_PER_LINE)
        {
            for (k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nColumns); k++)
            {
                Out << "[" << std::setw (4) << k << "]";
                Out << std::setw (15) << A(i,k) << " ";
            }
            Out << '\n';
        }
    }
}

template <class T>
void PrintMatrixColumnWise (CMatrix<T>& A, const std::string& heading,
                            std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: outputs a matrix into stream Out
// Input:    matrix, a heading, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    int i, j, k, nRows, nColumns;
    
    Out << '\n' << heading << '\n';
    Out << setiosflags (std::ios::left);

    nRows = A.GetRows(); nColumns = A.GetColumns();
    for (i=1; i <= nColumns; i++)
    {
        Out << "Column No: " << i << '\n';
        for (j=1; j <= nRows; j=j+NUM_ELEMENTS_PER_LINE)
        {
            for (k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nRows); k++)
            {
                Out << "[" << std::setw (4) << k << "]";
                Out << std::setw (15) << A(k,i) << " ";
            }
            Out << '\n';
        }
    }
}

template <class T>
void PrintMatrixColumn (CMatrix<T>& A, const std::string& heading, int i,
                        std::ostream& Out)
// ---------------------------------------------------------------------------
// Function: outputs a column of the matrix into stream Out
// Input:    matrix, a heading, column index, output stream object
// Output:   none
// ---------------------------------------------------------------------------
{
    int j, k, nRows;
    
    Out << '\n' << heading << '\n';
    Out << setiosflags (std::ios::left);

    nRows = A.GetRows();
    for (j=1; j <= nRows; j=j+NUM_ELEMENTS_PER_LINE)
    {
        for (k=j; k <= std::min(j+NUM_ELEMENTS_PER_LINE-1, nRows); k++)
        {
            Out << "[" << std::setw (4) << k << "]";
            Out << std::setw (15) << A(k,i) << " ";
        }
        Out << '\n';
    }
}

#endif
