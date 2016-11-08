#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

#include "my_matrix.h"
#include "my_complex_matrix.h"

//============================================================================
// 
//============================================================================
namespace LINALG
{


	typedef matrix < float >  matrixf;
	typedef  matrix < double > matrixd;
	//typedef  matrix < long double > matrixld;

	typedef matrix_complex < float >  matrix_cf;
	typedef  matrix_complex < double > matrix_cd;
	//typedef  matrix_complex < long double > matrix_cld;


	// operators that must be defined outside the class
	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix_complex<T>& matrix_complex<T>::operator | (matrix<T> &b)
	{
		if (this->NumRows() == b.NumRows() && this->NumCols() == b.NumCols())
		{

			out = Find_out(this->NumRows(), this->NumCols());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumCols(); j++)
				{
					(*out)(i, j) = get(i, j) * b(i, j);
				}
			}

			return *out;
		}

		return matrix_complex(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix_complex<T>& matrix_complex<T>::operator*(matrix<T> &b)
	{
		if (this->NumCols() == b.NumRows())
		{
			out = Find_out(this->NumRows(), b.NumCols());

			for (int r = 0; r < this->NumRows(); r++)
			{

				for (int c = 0; c < b.NumCols(); c++)
				{
					(*out)(r, c) = complex < T >(0.0, 0.0);
					for (int k = 0; k < this->NumCols(); k++)
					{
						(*out)(r, c) += get(r, k) * b(k, c);
					}
				}
			}

			return *out;
		}
		return matrix_complex(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix_complex<T>& matrix_complex<T>::operator+(matrix<T> &b)
	{
		if (this->NumRows() == b.NumRows() && this->NumCols() == b.NumCols())
		{
			out = Find_out(this->NumRows(), this->NumCols());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumCols(); j++)
				{
					(*out)(i, j) = this->get(i, j) + b.get(i, j);
				}
			}

			return *out;
		}
		return matrix_complex(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix_complex<T>& matrix_complex<T>::operator-(matrix<T> &b)
	{
		if (this->NumRows() == b.NumRows() && this->NumCols() == b.NumCols())
		{
			out = Find_out(this->NumRows(), this->NumCols());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumCols(); j++)
				{
					(*out)(i, j) = this->get(i, j) - b.get(i, j);
				}
			}

			return *out;
		}
		return matrix_complex(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix_complex<T>& matrix_complex<T>::operator = (matrix<T> &b)
	{
		if (this->NumRows()!= b.NumRows() || this->NumCols() != b.NumCols())
		{
			this->destroy();
			this->create(b.NumRows(), b.NumCols());
		}

		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumCols(); j++)
				data[i * SY + j] = complex<T> ( (b)(i, j), 0.0 );
		}
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix<T>& matrix<T>::operator = (matrix_complex<T> &b)
	{
		if (this->NumRows() != b.NumRows() || this->NumCols() != b.NumCols())
		{
			this->destroy();
			this->create(b.NumRows(), b.NumCols());
		}

		for (int i = 0; i < this->NumRows(); i++)
		{
			for (int j = 0; j < this->NumCols(); j++)
				data[i * SY + j] = =(b)(i, j).real();
		}
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix<T>& matrix<T>::operator | (matrix_complex<T> &b)
	{
		if (this->NumRows() == b.NumRows() && this->NumCols() == b.NumCols())
		{

			out = Find_out(this->NumRows(), this->NumCols());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumCols(); j++)
				{
					(*out)(i, j) = get(i, j) * b(i, j).real();
				}
			}

			return *out;
		}

		return matrix(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix<T>& matrix<T>::operator*(matrix_complex<T> &b)
	{
		if (this->NumCols() == b.NumRows())
		{
			out = Find_out(this->NumRows(), b.NumCols());

			for (int r = 0; r < this->NumRows(); r++)
			{

				for (int c = 0; c < b.NumCols(); c++)
				{
					(*out)(r, c) = 0.0;
					for (int k = 0; k < this->NumCols(); k++)
					{
						(*out)(r, c) += get(r, k) * b(k, c).real();
					}
				}
			}

			return *out;
		}
		return matrix(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix<T>& matrix<T>::operator+(matrix_complex<T> &b)
	{
		if (this->NumRows() == b.NumRows() && this->NumCols() == b.NumCols())
		{
			out = Find_out(this->NumRows(), this->NumCols());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumCols(); j++)
				{
					(*out)(i, j) = this->get(i, j) + b.get(i, j).real();
				}
			}

			return *out;
		}
		return matrix(0, 0);
	}

	//============================================================================
	//
	//============================================================================
	template <class T>
	matrix<T>& matrix<T>::operator-(matrix_complex<T> &b)
	{
		if (this->NumRows() == b.NumRows() && this->NumCols() == b.NumCols())
		{
			out = Find_out(this->NumRows(), this->NumCols());

			for (int i = 0; i < this->NumRows(); i++)
			{
				for (int j = 0; j < this->NumCols(); j++)
				{
					(*out)(i, j) = this->get(i, j) - b.get(i, j).real();
				}
			}

			return *out;
		}
		return matrix(0, 0);
	}


	//============================================================================
	//
	//============================================================================
	template< class T >
	inline matrix<T> operator*(T s, matrix<T> a)
	{
		return a * s;
	}


	//============================================================================
	//
	//============================================================================
	template< class T >
	inline matrix_complex<T> operator*(T s, matrix_complex<T> a)
	{
		return a * s;
	}

	//============================================================================
	//
	//============================================================================
	template< class T >
	inline matrixd operator*(T s, matrixd a)
	{
		return a * s;
	}
	//============================================================================
	//
	//============================================================================
	template< class T >
	inline matrixf operator*(T s, matrixf a)
	{
		return a * s;
	}
	//============================================================================
	//
	//============================================================================
	template< class T >
	inline matrix_cf operator*(T s, matrix_cf a)
	{
		return a * s;
	}
	//============================================================================
	//
	//============================================================================
	template< class T >
	inline matrix_cd operator*(T s, matrix_cd a)
	{
		return a * s;
	}

}


#endif