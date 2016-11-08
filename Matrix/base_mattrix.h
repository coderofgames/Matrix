#ifndef BASE_MATRIX
#define BASE_MATRIX

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "Utils.h"
#include <complex>
#include <limits>

using std::complex;

using std::cout;
using std::endl;

using namespace std;


//============================================================================
//
//============================================================================
namespace LINALG
{

	template <class T> class matrix;
	template <class T> class matrix_complex;

	template<class T>
	bool almost_equal(T a, T b)
	{
		return (std::abs(a - b) < std::numeric_limits<T>::epsilon());
	}



	// stack overflow
	template<class T>
	bool almost_equal(T a, T b, int ulps)
	{
		return (std::abs(a - b) < std::numeric_limits<T>.epsilon() * std::abs(a + b) * ulps)
			|| std::abs(x - y) < std::numeric_limits<T>::min();
	}

	//============================================================================
	//
	//============================================================================
	template< class T >
	static inline void SWAP(T &a, T &b)
	{
		T temp = a;
		a = b;
		b = temp;
	}






	//============================================================================
	//
	//============================================================================
	template<class T>
	static inline T sgn(T x)
	{
		if (x > 0.0 + numeric_limits<T>::epsilon())
			return 1.0;
		if (x < 0.0 - numeric_limits<T>::epsilon())
			return -1.0;

		return x;
	}


	//============================================================================
	// https://en.wikipedia.org/wiki/Sign_function
	//============================================================================
	template< class T >
	static inline complex <T> sgn(complex <T> z)
	{
		if (almost_equal<T>(z.real(), 0.0) && almost_equal<T>(z.imag(), 0.0))
			return z;

		T abs_z_inv = 1 / std::abs(z);

		return z * abs_z_inv;
	}


	//============================================================================
	//
	//============================================================================
	template<class T>
	static T CSGN(complex<T> x)
	{
		if (x.real() > 0.0 + numeric_limits<T>::epsilon() )
		{
			return 1.0;

		}
		else if (x.real() < 0.0 - numeric_limits<T>::epsilon())
		{
			return -1.0;
		}
		else if ( almost_equal <T> ( x.real(), 0.0) )
		{
			if (x.imag() > 0.0 + numeric_limits<T>::epsilon())
			{
				return 1.0;
			}
			else if (x.imag() < 0.0 - numeric_limits<T>::epsilon())
			{
				return -1.0;
			}
		}
		return 0.0;
	}

	//============================================================================
	//
	//============================================================================
	static inline double round_to_n_digits(double x, int n)
	{
		double scale = pow(10.0, ceil(log10(fabs(x))) + n);

		return round(x * scale) / scale;
	}





	//============================================================================
	//
	//============================================================================
	static inline float RandomFloat(float min, float max)
	{
		float r = (float)rand() / (float)RAND_MAX;
		return min + r * (max - min);
	}

	//============================================================================
	//
	//============================================================================
	static inline float RandomInt(int min, int max)
	{
		float r = (float)rand() / (float)RAND_MAX;
		return (int)((float)min + r * float(max - min));
	}

	//============================================================================
	//
	//============================================================================
	template<class Scalar>
	static inline Scalar RandomScalar(Scalar min, Scalar max)
	{
		Scalar r = (Scalar)rand() / (Scalar)RAND_MAX;
		return min + r * (max - min);
	}

	//============================================================================
	//
	//============================================================================
	template< class Scalar >
	static inline complex<Scalar> RandomComplex(Scalar min, Scalar max)
	{
		return complex<Scalar>(RandomScalar(min, max), RandomScalar(min, max));
	}






};

#endif