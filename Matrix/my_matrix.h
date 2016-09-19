#ifndef MY_MATRIX_H
#define MY_MATRIX_H

#include <iostream>
#include "Utils.h"

using std::cout;
using std::endl;

template< class T >
inline void SWAP(T &a, T &b)
{
	T temp = a;
	a = b;
	b = temp;
}

inline float RandomFloat(float min, float max)
{
	float r = (float)rand() / (float)RAND_MAX;
	return min + r * (max - min);
}

inline float RandomInt(int min, int max)
{
	float r = (float)rand() / (float)RAND_MAX;
	return (int)((float)min + r * float(max - min));
}

class vector2d
{
public:

	float v[2];

	float operator[](unsigned int idx) { return (idx < 2 ? v[idx] : 0.0f); }
	void operator=(vector2d b){ v[0] = b.v[0]; v[1] = b.v[1]; }
};


class LINALG
{

private:
template<class T>
class matrix
{
public:
	matrix<T>(){
		this->m_sizeX = 0; m_sizeY = 0;
		data = 0;
	}


	matrix<T>(matrix *p){
		this->m_sizeX = p->m_sizeX; m_sizeY = p->m_sizeY;
		this->create();
		for (int i = 0; i < m_sizeX; i++)
		{
			for (int j = 0; j < m_sizeY; j++)
				data[i][j] = (*p)[i][j];
		}
		is_transposed = false;
	}

	matrix<T>(matrix &p){
		this->m_sizeX = p.m_sizeX; m_sizeY = p.m_sizeY;
		this->create();
		for (int i = 0; i < m_sizeX; i++)
		{
			for (int j = 0; j < m_sizeY; j++)
				data[i][j] = p[i][j];
		}
		is_transposed = false;
	}

	matrix<T>(unsigned int n, unsigned int m)
	{
		m_sizeX = n;
		m_sizeY = m;
		create();
	}
	~matrix<T>()
	{
		destroy();
		if (out) delete out;
	}

	void destroy()
	{
		if (data != 0)
		{
			for (int i = 0; i < m_sizeX; i++)
			{
				delete[] data[i];
			}
			delete[] data;
		}
		data = 0;
		is_transposed = false;
	}
	void create()
	{
		if (m_sizeY > 0 && m_sizeX > 0)
		{
			data = new T*[m_sizeX];
			for (int i = 0; i < m_sizeX; i++)
			{
				data[i] = new T[m_sizeY];
				for (int j = 0; j < m_sizeY; j++)
					data[i][j] = 0.0f;
			}
		}
		is_transposed = false;
	}

	void operator=(matrix &b)
	{
		this->destroy();
		this->m_sizeX = b.NumRows();
		this->m_sizeY = b.NumColumns();
		this->create();
		for (int i = 0; i < m_sizeX; i++)
		{
			for (int j = 0; j < m_sizeY; j++)
				data[i][j] = b(i, j);
		}


	}

	void operator=(matrix *b)
	{
		this->destroy();
		this->m_sizeX = b->NumRows();
		this->m_sizeY = b->NumColumns();
		this->create();
		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumColumns(); j++)
				data[i][j] = (*b)(i, j);
		}
	}





	T& operator()(unsigned int i, unsigned int j)
	{
		return get(i, j);
	}

	T& get(unsigned int i, unsigned int j)
	{
		T null_return = 0.0f;
		if (i < NumRows() && j < NumColumns())
			return is_transposed ? data[j][i] : data[i][j];
		else return null_return;
	}
	// Hadamard element wise product
	matrix& operator | (matrix &b)
	{
		if (this->NumColumns() == b.NumColumns() && this->NumRows() == b.NumRows())
		{
			if (out) delete out;
			out = new matrix(this->NumRows(), this->NumColumns());

			for (int i = 0; i < this->NumRows(); i++)
			{

				for (int j = 0; j < this->NumColumns(); j++)
				{
					//for (int k = 0; k < this->NumColumns(); k++)
					{
						(*out)(i, j) = get(i, j) * b(i, j);
					}
				}
			}

			return *out;
		}

		return matrix(0, 0);
	}

	matrix& operator*(matrix &b)
	{
		if (b.NumColumns() == 1 && b.NumRows() == 1)
		{
			return (*this) * b(0, 0);
		}
		if (this->NumColumns() == b.NumRows())
		{
			if (out) delete out;
			out = new matrix(this->NumRows(), b.NumColumns());

			for (int i = 0; i < this->NumRows(); i++)
			{

				for (int j = 0; j < b.NumColumns(); j++)
				{
					for (int k = 0; k < this->NumColumns(); k++)
					{
						(*out)(i, j) += get(i, k) * b(k, j);
					}
				}
			}

			return *out;
		}
		return matrix(0, 0);
	}

	matrix& operator*(T s)
	{
		if (out) delete out;
		out = new matrix(this->NumRows(), this->NumColumns());

		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumColumns(); j++)
			{
				(*out)(i, j) = get(i, j) * s;
			}
		}
		return *out;
	}

	matrix& operator/(T s)
	{
		if (out) delete out;
		out = new matrix(this->NumRows(), this->NumColumns());

		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumColumns(); j++)
			{
				(*out)(i, j) = get(i, j) / s;
			}
		}
		return *out;
	}
	matrix& operator+(T s)
	{
		if (out) delete out;
		out = new matrix(this->NumRows(), this->NumColumns());

		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumColumns(); j++)
			{
				(*out)(i, j) = get(i, j) + s;
			}
		}
		return *out;
	}
	matrix& operator+(matrix &b)
	{
		if (this->NumColumns() != b.NumColumns() || this->NumRows() != b.NumRows())
			return matrix(0, 0);
		else
		{
			if (out) delete out;
			out = new matrix(this->NumRows(), this->NumColumns());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumColumns(); j++)
				{
					(*out)(i, j) = this->get(i, j) + b.get(i, j);
				}
			}

			return *out;
		}
		return matrix(0, 0);
	}

	matrix& operator-(matrix &b)
	{
		if (this->NumColumns() != b.NumColumns() || this->NumRows() != b.NumRows())
			return matrix(0, 0);
		else
		{
			if (out) delete out;
			out = new matrix(this->NumRows(), this->NumColumns());

			for (int i = 0; i < this->NumRows(); i++)
			{

				for (int j = 0; j < this->NumColumns(); j++)
				{
					(*out)(i, j) = this->get(i, j) - b.get(i, j);
				}
			}

			return *out;
		}
		return matrix(0, 0);
	}

	void print()
	{
		if (this->NumRows() == 0 || this->NumColumns() == 0 || this->data == 0)
		{
			cout << "Empty Matrix" << endl;
			return;
		}
		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumColumns(); j++)
			{
				cout << get(i, j) << "  ";
			}
			cout << endl;
		}
	}



	T trace()
	{
		T sum = 0.0f;
		if (m_sizeX != m_sizeY) return 0.0f;

		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumColumns(); j++)
			{
				if (i == j) sum += get(i, j);
			}

		}
		return sum;
	}



	bool ContainsNAN()
	{
		for (int i = 0; i < m_sizeX; i++)
		{
			for (int j = 0; j < m_sizeY; j++)
			{
				if (data[i][j] != data[i][j])
					return true;
			}
		}
		return false;
	}

	void transpose()
	{
		is_transposed = !is_transposed;
	}
	//void T()
	//{
	//	is_transposed = !is_transposed;
	//}


	inline unsigned int NumRows()
	{
		return  (is_transposed ? m_sizeY : m_sizeX);
	}
	inline unsigned int NumColumns()
	{
		return (is_transposed ? m_sizeX : m_sizeY);
	}

	inline bool IsSquare()
	{
		return m_sizeX == m_sizeY;
	}

	bool IsSymmetric()
	{
		if (!this->IsSquare())
		{
			return false;
		}
		
		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j <= i; j++)
			{
				if (get(i, j) != get(j, i)) return false;
			}
		}

		return true;
	}

	bool IsPositiveDefinite()
	{
		if (!this->IsSquare())
		{
			cout << "Error: Positive Definite Matrixces must be square" << endl;
			return false;
		}

		matrix X(this->NumRows(), 1);

		for (int i = 0; i < this->NumRows(); i++)
			X(i, 0) = RandomFloat(1, 5);

		matrix XT = X;
		XT.transpose();

		matrix res = XT * (*this) * X;

		for (int i = 0; i < res.NumColumns(); i++)
		{
			for (int j = 0; j < res.NumRows(); j++)
			{
				if (res(j, i) <= 0) return false;
			}
		}
		return true;
	}

	//============================================================================
	// For solution of the system Ax = b, 
	// relies on the condition that this matrix A is square and 
	// that b has the same number of rows as A
	//
	// inputs b, outputs x, 
	//
	// For more information please consult Krezig: Advanced Engineering Mathematics, sec 19.1
	//
	//============================================================================
	matrix& Gauss(matrix &b)
	{
		if (!this->IsSquare())
		{
			cout << "Error in Gauss calculation: System Matrix must be square." << endl;
			return b;
		}

		if (b.NumRows() != this->NumRows())
		{
			cout << "Error in Gauss calculation: b must have the same number of rows as A" << endl;
			return b;
		}

		unsigned int n = b.NumRows();
		
		for (int k = 0; k < n - 1; k++)
		{
			bool bSolutionExists = false;
			unsigned int j = k + 1;
			for (j = k + 1; j < n; j++)
			{
				if (get(j, k) != 0)
				{
					bSolutionExists = true;
					break;
				}
			}

			if (!bSolutionExists)
			{
				cout << "No Unique Solution Exists" << endl;
				return b; // no solution
			}

			for (int i = 0; i < n; i++)
			{
				SWAP(get(j, i), get(k, i));

				SWAP(b(j, 0), b(k, 0));
			}
			for (j = k + 1; j < n; j++)
			{
				T mjk = get(j, k) / get(k, k);
				for (int p = k ; p < n; p++)
				{
					get(j, p) = get(j, p) - mjk * get(k, p);
				}
				b(j, 0) = b(j, 0) - mjk * b(k, 0);
			}
		}
		if (get(n - 1, n - 1) == 0)
		{
			cout << "No Unique Solution Exists" << endl;
			return b; // no solution
		}

		if (out) delete out;
		out = new matrix(n, 1);

		(*out)(n - 1, 0) = b(n - 1, 0) / get(n - 1, n - 1);

		for (int i = n - 2; i > -1; i--)
		{
			T the_sum = 0;
			for (int j = i + 1; j < n; j++) 
				the_sum += get(i, j)*(*out)(j, 0);
			
			(*out)(i, 0) = (1 / get(i, i)) * (b(i, 0) - the_sum);
		}

		return *out;
	}

	// generalization of the reduction to triangular form above
	bool ReduceToUpperTriangularForm()
	{
		if (!this->IsSquare())
		{
			cout << "Error (ReduceToLowerTriangularForm) : System Matrix must be square." << endl;
			return false;
		}

		unsigned int n = this->NumRows();

		// this part is almost exactly the same as the section from the Gauss method 
		for (int k = 0; k < n - 1; k++)
		{
			bool bSolutionExists = false;
			unsigned int j = n - 1;
			for (j = k + 1; j < n; j++)
			{
				if (get(j, k) != 0)
				{
					bSolutionExists = true;
					break;
				}
			}

			if (!bSolutionExists) return false;

			for (int i = 0; i < n; i++)
			{
				SWAP(get(j, i), get(k, i));	
			}

			for (j = k + 1; j < n; j++)
			{
				T mjk = get(j, k) / get(k, k);
				for (int p = k; p < n; p++)
				{
					get(j, p) = get(j, p) - mjk * get(k, p);
				}
			}
		}	
	}

	// interesting way of evaluating the determinant
	T Determinant()
	{
		if( this->ReduceToUpperTriangularForm() )
			return this->trace();
		
		return 0;
	}

	//============================================================================
	// Matrix Inversion via Gauss-Jordan elimination
	// The matrix should be square, returns the inverse
	//
	//============================================================================
	matrix& Gauss_Jordan()
	{
		if (!this->IsSquare())
		{
			cout << "Error in Gauss calculation: System Matrix must be square." << endl;
			return (*this);
		}


		unsigned int n = this->NumRows();

		if (out) delete out;
		out = new matrix(n,n);

		out->Identity();


		// this part is almost exactly the same as the section from the Gauss method 
		for (int k = 0; k < n-1 ; k++)
		{
			bool bSolutionExists = false;
			unsigned int j = n-1;
			for (j = k + 1; j < n; j++)
			{
				if (get(j, k) != 0)
				{
					bSolutionExists = true;
					break;
				}
			}

			if (!bSolutionExists)
			{
				cout << "Error (Gauss_Jordan): No Unique Solution Exists." << endl;
				return (*this);
			}
	
			for (int i = 0; i < n; i++)
			{
				SWAP(get(j, i), get(k, i));

				SWAP((*out)(j, i), (*out)(k, i));
			}
			for (j = k + 1; j < n; j++)
			{
				T mjk = get(j, k) / get(k, k);
				
				for (int p = k; p < n; p++)
				{
					
					get(j, p) = get(j, p) - mjk * get(k, p);
					
					(*out)(j, p) = (*out)(j, p) - mjk * (*out)(k, p);
				}
				//b(j, 0) = b(j, 0) - mjk * b(k, 0);
			}


		}

		// make sure all the elements on the diagonal are 1 by dividing all 
		// rows by the elements in the diagonals of that row
		// need this for the next loop
		for (int k = 0; k < n; k++)
		{
			T val = get(k, k);
			for (int c = 0; c < n; c++)
			{
				get(k, c) = get(k, c) / val;
				(*out)(k, c) = (*out)(k, c) / val;
			}
		}

		// starting at the n-1 th row, we have
		// a matrix in echelon form with ones on the main diagonal
		// the 
		for (int r = n - 2; r >-1; r--) // n-2 here because we index from zero instead of one
		{


			// the next element after the diagonal is Row+1
			for (int c = r+1; c < n; c++)
			{
				// get the value in the row / column position, we will make the element in (r,c) into zero
				// by multiplying by the *next* row, elemnts before the diagonal are already zero, and 
				// the next elements in this row will be handled next ...
				// ...
				T val = get(r, c); 
				for (int p =0; p < n; p++)
				{
					// we use c to use a mulitiple of c's row, subtracting a multiple of the 1's on the diagonal
					get(r, p) = get(r, p) - get(c, p)*val; 

					// all operations are mirrored on the other *output* matrix, or the inverse
					(*out)(r, p) = (*out)(r, p) - (*out)(c,p)*val; 
				}


			}

		}
		

		return *out;

	}

	//============================================================================
	// Solve a matrix system of the form Ly = b by back substitution
	// L, is Lower triangular, y is the solution vector, b is the output
	// inputs L, y, b
	//============================================================================
	matrix& Solve_Lower_TriangularSystem(matrix& L, matrix & y, matrix& b)
	{
		if (!L.IsSquare())
		{
			cout << "Error (Solve_Lower_TriangularSystem): Lower Triangular Matrix should be square" << endl;
			return y;
		}

		if (!y.NumRows() == L.NumRows()) {
			cout << "Error (Solve_Lower_TriangularSystem): vector y should have the same number of rows as L" << endl;
			return y;
		}
		if (!b.NumRows() == y.NumRows())
		{
			cout << "Error (Solve_Lower_TriangularSystem): vector y should have the same number of rows as b" << endl;
			return y;
		}

		y(0, 0) = b(0, 0) / L(0, 0);

		for (int i = 1; i < L.NumRows(); i++)
		{
			T sum0 = 0.0f;
			for (int s = 0; s < i; s++)
			{
				sum0 += L(i, s) * y(s,0);
			}
			y(i, 0) = (1 / L(i, i)) * (b(i, 0) - sum0);
		}

		return y;
	}

	//============================================================================
	// Solve a system of the form Ux = y by back substitution
	// U is upper triangular, x is the solution vector, y is the output
	// inputs U,x,y
	//============================================================================
	matrix& Solve_Upper_TriangularSystem(matrix& U, matrix & x, matrix& y)
	{
		if (!U.IsSquare()) {
			cout << "Error (Solve_Upper_TriangularSystem): matrix should be square" << endl;
			return y;
		}
		if (!y.NumRows() == U.NumRows())
		{
			cout << "Error (Solve_Upper_TriangularSystem): vector y should have same number of rows as matrix" << endl;
			return y;
		}
		if (!x.NumRows() == y.NumRows()) {
			cout << "Error (Solve_Upper_TriangularSystem): x should be same shape as y" << endl;
			return y;
		}
		
		int n = x.NumRows();

		x(n - 1, 0) = y(n - 1, 0) / U(n - 1, n - 1);

		for (int i = n - 2; i > -1; i--)
		{
			T sum0 = 0;
			for (int j = i + 1; j < n; j++)
				sum0 += U(i, j)*x(j, 0);

			x(i, 0) = (1 / U(i, i)) * (y(i, 0) - sum0);
		}

		return x;
	}

	//============================================================================
	// Solves a system Ax = b via the Cholesky method
	// this matrix must be symmetric and positive definite
	// inputs b, output vector
	// returns x, solution vector
	//============================================================================
	matrix& Cholesky(matrix& b)
	{
		if (!this->IsSymmetric()) {
			cout << "Error (Cholesky): Matrix should be symetric" << endl;
			return b;
		}
		if (!this->IsPositiveDefinite())
		{
			cout << "Error (Cholesky): Matrix should be positive definite" << endl;
			return b;
		}

		int n = this->NumRows();
		matrix M(n, n);

		M(0, 0) = sqrt(get(0, 0));


		for (int j = 1; j < n; j++)
		{
			M(j, 0) = get(j, 0) / M(0, 0);
		}

		for (int k = 1; k < this->NumColumns(); k++)
		{
			for (int j = k; j < n; j++)
			{
				if (j == k)
				{
					T sum1 = 0;
					for (int s = 0; s < j ; s++) sum1 += M(j, s) * M(j, s);

					M(j, j) = sqrt(get(j, j) - sum1);
				}
				else
				{
					T sum2 = 0;
					for (int s = 0; s < k ; s++) sum2 += M(j, s) * M(k, s);

					M(j, k) = (1 / M(k, k)) * (get(j, k) - sum2);
				}
			}
		}

		matrix y(n, 1);
		Solve_Lower_TriangularSystem(M, y, b);
		M.transpose();

		if (out) delete out;
		out = new matrix(n, 1);

		Solve_Upper_TriangularSystem(M, *out, y);

		return *out;
	}

	inline bool EqualSize(matrix& b){ return (this->NumRows() == b.NumRows()) && (this->NumColumns() == b.NumColumns());  }

	void Identity()
	{
		if (!this->IsSquare())
		{
			cout << "Error (Identity): Identity Matrix must be square" << endl;
			return;
		}
		for (int i = 0; i < this->NumColumns(); i++)
		{
			for (int j = 0; j < this->NumRows(); j++)
			{
				if (i == j)
				{
					get(j, i) = 1.0f;
				}
				else
				{
					get(j, i) = 0.0f;
				}
				
			}
		}
	}

	//============================================================================
	// The LU decomposition using doolittles method, this function
	// does not solve a system, it simple decompses this matrix into two
	//matrices such that A = LU
	//============================================================================
	int LU_Decomposition_Doolittle(matrix& L, matrix& U)
	{
		if (!L.EqualSize(U))
		{
			cout << "Error (LU_Decomposition_Doolittle): L and U must be equal size" << endl;
			return 0;
		}
		if (!this->EqualSize(L))
		{
			cout << "Error (LU_Decomposition_Doolittle):  L and U must be the same size as the source matrix" << endl;
			return 0;
		}
		if (!this->IsSquare())
		{
			cout << "Error (LU_Decomposition_Doolittle):  source matrix must be square" << endl;
			return 0;
		}



		int n = this->NumColumns();

		L.Identity();

		for (int k = 0; k < n; k++)
		{
			U(0, k) = get(0, k);

			L(k, 0) = get(k, 0) / U(0, 0); // loop 1 with rows indexed with k, ignore cols

			L(k, k) = 1.0f; // loop 1 as L(k,k)
		}

		for (int j = 1; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				if (k >= j)
				{
					T sum0 = 0.0f;
					for (int s = 0; s < j; s++) sum0 += L(j, s) * U(s, k);

					U(j, k) = get(j, k) - sum0;
  				}
				else
				{
					T sum1 = 0.0f;
					for (int s = 0; s < k; s++) sum1 += L(j, s) * U(s, k);

					L(j, k) = (1 / U(k, k))*(get(j, k) - sum1);
				}
				
			}
		}

		return 1;
	}

	//============================================================================
	// Solves the system via doolittles method, for the equation Ax = b
	// it first factorizes A into LU such that A=LU, then solves the 
	// lower triangular system Ly = b and then the upper triangular system Ux=y
	// returns x
	//
	//============================================================================
	matrix& Solve_System_Doolittle(matrix& b)
	{

		if (!this->IsSquare())
		{
			cout << "Error (Solve_System_Doolittle):  source matrix must be square" << endl;
			return b;
		}
		if (!(b.NumRows() == this->NumRows()))
		{
			cout << "Error (Solve_System_Doolittle):  matrix b must have same number of rows as source matrix" << endl;
			return b;
		}

		
		int n = this->NumColumns();

		matrix L(n, n);
		matrix U(n, n);
		

		this->LU_Decomposition_Doolittle(L, U);

		matrix y(n, 1);
		Solve_Lower_TriangularSystem(L, y, b);

		if (out) delete out;
		out = new matrix(n, 1);

		Solve_Upper_TriangularSystem(U, *out, y);

		return *out;
	}


	//============================================================================
	// Very similar to the Doolittle method, solves A = LU for L and U
	// inputs require L and U to be of the same dimensionality as A
	//============================================================================
	int LU_Decomposition_Crout(matrix& L, matrix& U)
	{
		if (!L.EqualSize(U))
		{
			cout << "Error (LU_Decomposition_Doolittle): L and U must be equal size" << endl;
			return 0;
		}
		if (!this->EqualSize(L))
		{
			cout << "Error (LU_Decomposition_Doolittle):  L and U must be the same size as the source matrix" << endl;
			return 0;
		}
		if (!this->IsSquare())
		{
			cout << "Error (LU_Decomposition_Doolittle):  source matrix must be square" << endl;
			return 0;
		}



		int n = this->NumColumns();

		U.Identity();

		for (int j = 0; j < n; j++)
		{
			L(j, 0) = get(j, 0);

			U(0, j) = get(0, j) / L(0, 0); 

			U(j, j) = 1.0f; 
			
		}

		for (int k = 1; k < n; k++)
		{
			for (int j = 0; j < n; j++)
			{
				if (j >= k)
				{
					T sum0 = 0.0f;
					for (int s = 0; s < k; s++) sum0 += L(j, s) * U(s, k);

					L(j, k) = get(j, k) - sum0;
  				}
				else
				{
					T sum1 = 0.0f;
					for (int s = 0; s < j; s++) sum1 += L(j, s) * U(s, k);

					U(j, k) = (1 / L(j,j))*(get(j, k) - sum1);
				}		
			}
		}

		return 1;
	}

	//============================================================================
	// Solves the system via Crout method, for the equation Ax = b
	// it first factorizes A into LU such that A=LU, then solves the
	// lower triangular system Ly = b and then the upper triangular system Ux=y
	// returns x
	//
	//============================================================================
	matrix& Solve_System_Crout(matrix& b)
	{

		if (!this->IsSquare())
		{
			cout << "Error (Solve_System_Doolittle):  source matrix must be square" << endl;
			return b;
		}
		if (!(b.NumRows() == this->NumRows()))
		{
			cout << "Error (Solve_System_Doolittle):  matrix b must have same number of rows as source matrix" << endl;
			return b;
		}


		int n = this->NumColumns();

		matrix L(n, n);
		matrix U(n, n);


		this->LU_Decomposition_Crout(L, U);

		matrix y(n, 1);
		Solve_Lower_TriangularSystem(L, y, b);

		if (out) delete out;
		out = new matrix(n, 1);

		Solve_Upper_TriangularSystem(U, *out, y);

		return *out;
	}

	//============================================================================
	// Solves the system via Doolittles or Crouts method, for the equation Ax = b
	// If the matrices L and U satisfy A = LU (where A is the system matix)
	// it solves the lower triangular system Ly = b and then the upper triangular system Ux=y
	// returns x
	//
	//============================================================================
	matrix& Solve_System_LU(matrix& L, matrix& U, matrix& b)
	{

		if (!this->IsSquare())
		{
			cout << "Error (Solve_System):  source matrix must be square" << endl;
			return b;
		}
		if (!L.EqualSize(U))
		{
			cout << "Error (Solve_System): L and U must be equal size" << endl;
			return b;
		}
		if (!this->EqualSize(L))
		{
			cout << "Error (Solve_System):  L and U must be the same size as the source matrix" << endl;
			return b;
		}

		int n = this->NumColumns();

		//this->LU_Decomposition_Crout(L, U);

		matrix y(n, 1);
		Solve_Lower_TriangularSystem(L, y, b);

		if (out) delete out;
		out = new matrix(n, 1);

		Solve_Upper_TriangularSystem(U, *out, y);

		return *out;
	}

	//============================================================================
	// For solution of the system Ax = b, relies on the coondition
	// that this matrix A is square and that b has the same number of rows as A
	//
	// inputs b, initial estimation x0, tolerance (eps > 0) small number, MAX_ITERATIONS for number of loops,
	// outputs x or reports failure and outputs ,
    //
	// For more information please consult Krezig: Advanced Engineering Mathematics, 
	//
	//============================================================================
	matrix & Gauss_Seidel(matrix& b, matrix& x0, float tolerance, unsigned int MAX_ITERATIONS)
	{
		if (x0.NumRows() != this->NumRows())
		{
			cout << "Gauss Seidel error: x0 must have the same number of rows as this matrix" << endl;
			return x0;
		}
		
		if (b.NumRows() != this->NumRows())
		{
			cout << "Gauss Seidel error: b must have the same number of rows as this matrix" << endl;
			return x0;
		}

		int n = this->NumRows();

		if (out) delete out;
		out = new matrix(n, 1);

		//for (int i = 0; i < n; i++) (*out)(i, 0) = x0(i, 0);
		
		for ( int m = 0; m < MAX_ITERATIONS; m++)
		{
			for ( int j = 0; j < n; j++)
			{
				T sum1 = 0.0f;
				for (int k = 0; k < j; k++)
				{
					sum1 += get(j, k) * (*out)(k, 0);
				}
				T sum2 = 0.0f;
				for (int k = j+1; k < n; k++)
				{
					sum2 += get(j, k) * x0(k, 0);
				}
				
				(*out)(j, 0) = (-1 / get(j, j))*(sum1 + sum2 - b(j, 0));
			}
			float max_magnitude = 0.0f;
			for (int j = 0; j < n; j++)
			{
				float magnitude = abs((*out)(j, 0) - x0(j, 0));
				if (magnitude > max_magnitude) max_magnitude = magnitude;

				x0(j, 0) = (*out)(j, 0);
			}
			if (max_magnitude < tolerance) return *out;

		//	for (int j = 0; j < n; j++) 
		}
		cout << "No output satisfying the tolerance condition obtained after N iteration steps." << endl;
		return x0;
	}

	
	bool is_transposed = false;
	unsigned int m_sizeX = 0;
	unsigned int m_sizeY = 0;
	
	T** data;


//private:
	T* operator[](unsigned int a)
	{
		if (a < m_sizeX)
			return data[a];
		else return 0;
	}
	

	matrix *out = 0;


};


public:

	typedef matrix<float> matrixf;
	typedef matrix<double> matrixd;
};


#endif