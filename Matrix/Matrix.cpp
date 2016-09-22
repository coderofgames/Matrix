// matrixf.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "my_matrix.h"
#include <iostream>

using std::cout;
using std::endl;

typedef LINALG::matrixf matrixf;
typedef LINALG::matrixd matrixd;




void TestHessenburg_2()
{

	matrixf H_(5,5);

	// using the wisdom from
	// http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0258255-3/S0025-5718-1969-0258255-3.pdf
	// I have decided to only use Housedholder methods
	// this matrix shows the evidence of the stability issue ...

	H_(0, 0) = 6.0;
	H_(0, 1) = 1.0;
	H_(0, 2) = -2.0;
	H_(0, 3) = 19.0;
	H_(0, 4) = 4.0;


	H_(1, 0) = 1.0;
	H_(1, 1) = 4.0;
	H_(1, 2) = 2.0;
	H_(1, 3) = 1.0;
	H_(1, 4) = 3.0;


	H_(2, 0) = 7.0;
	H_(2, 1) = 12.0;
	H_(2, 2) = 5.0;
	H_(2, 3) = 11.0;
	H_(2, 4) = -1.0;


	H_(3, 0) = -3.0;
	H_(3, 1) = 21.0;
	H_(3, 2) = 4.0;
	H_(3, 3) = 1.0;
	H_(3, 4) = 3.0;


	H_(4, 0) = 4.0;
	H_(4, 1) = 6.0;
	H_(4, 2) = -1.0;
	H_(4, 3) = 3.0;
	H_(4, 4) = 7.0;


	matrixf H__3 = H_;

	H__3.Householder_Tridiagonalize();
	H__3.print(2);
	cout << endl;
	cout << "Eigenvalues are on the diagonal or in complex conjugate pair" << endl;
	matrixf eigen_values(H__3.NumColumns(), 2);
	H__3.QR_algorithm(eigen_values);
	H__3.print(2);
	eigen_values.print(3);
}


void TestHessenburg_QR_real()
{

	matrixf H_(4,4);

	// using the wisdom from
	// http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0258255-3/S0025-5718-1969-0258255-3.pdf
	// I have decided to only use Housedholder methods
	// this matrix shows the evidence of the stability issue ...

	H_(0, 0) = 6.0;
	H_(0, 1) = -sqrt(18.0);
	H_(0, 2) = 0.0;
	H_(0, 3) = 0.0;



	H_(1, 0) = -sqrt(18.0);
	H_(1, 1) = 7.0;
	H_(1, 2) = sqrt(2.0);
	H_(1, 3) = 0.0;



	H_(2, 0) = 0.0;
	H_(2, 1) = sqrt(2.0);
	H_(2, 2) = 6.0;
	H_(2, 3) = 0.0;



	H_(3, 0) = 0.0;
	H_(3, 1) = 0.0;
	H_(3, 2) = 0.0;
	H_(3, 3) = 3.0;






	matrixf H__3 = H_;

//	H__3.Householder_Tridiagonalize();
//	H__3.print(2);
	cout << endl;
	cout << "Eigenvalues are on the diagonal or in complex conjugate pair" << endl;
	matrixf eigen_values(H__3.NumColumns(), 2);
	H__3.QR_algorithm(eigen_values);
	H__3.print(2);
	eigen_values.print(3);

}
int main(int argc, char* argv[])
{
	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination" << endl;
	matrixf m(3, 3);
	m(0, 0) = 0;
	m(0, 1) = 8;
	m(0, 2) = 2;
	
	m(1, 0) = 3;
	m(1, 1) = 5;
	m(1, 2) = 2;

	m(2, 0) = 6;
	m(2, 1) = 2;
	m(2, 2) = 8;

	cout << endl;
	cout << "Printing A" << endl;
	m.print(4);
	
	matrixf b(3, 1);
	b(0, 0) = -7;
	b(1, 0) = 8;
	b(2, 0) = 26;

	cout << endl;
	cout << "Printing b" << endl;
	m.print(4);
	
	cout << endl;
	cout << "Solvng the system A x = b for x" <<endl;
	cout << endl;

	matrixf sol = m.Gauss_Elimination(b);
	
	cout << "printing solution" << endl;
	
	sol.print(4);
	
	cout << endl;


	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination 2" << endl;

	matrixf m2(3, 3);
	m2(0, 0) = 4;
	m2(0, 1) = 2;
	m2(0, 2) = 14;

	m2(1, 0) = 2;
	m2(1, 1) = 17;
	m2(1, 2) = -5;

	m2(2, 0) = 14;
	m2(2, 1) = -5;
	m2(2, 2) = 83;

	cout << endl;
	cout << "printing A2" << endl;

	m2.print(4);

	matrixf b2(3, 1);
	b2(0, 0) = 14;
	b2(1, 0) = -101;
	b2(2, 0) = 155;

	cout << endl;
	cout << "Printing b" << endl;

	b2.print(4);


	cout << endl;
	cout << "Solvng the system A x = b" << endl;
	
	matrixf sol2 = m2.Gauss_Elimination(b2);
	
	cout << endl;
	cout << "printing solution (x) " << endl;
	sol2.print(4);


	cout << endl;
	cout << "=================================================" << endl;
	cout << "matrixf Multiplication A * A2" << endl;
	matrixf D = m * m2;
	D.print(4);

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;

	cout << "tesing Gauss Seidel" << endl;

	matrixf M(4, 4);
	M(0, 0) = 1;
	M(0, 1) = -0.25;
	M(0, 2) = -0.25;
	M(0, 3) = 0;

	M(1, 0) = -0.25;
	M(1, 1) = 1;
	M(1, 2) = 0;
	M(1, 3) = -0.25;

	M(2, 0) = -0.25;
	M(2, 1) = 0;
	M(2, 2) = 1;
	M(2, 3) = -0.25;

	M(3, 0) = 0;
	M(3, 1) = -0.25;
	M(3, 2) = -0.25;
	M(3, 3) = 1;

	cout << endl;
	cout << "Printing M" << endl;
	M.print(4);

	cout << endl;

	if (M.IsSymmetric()) cout << "M is Symmetric" << endl;
	if (M.IsPositiveDefinite()) cout << "M is positive definite" << endl;

	matrixf B(4, 1);

	B(0, 0) = 50;
	B(1, 0) = 50;
	B(2, 0) = 25;
	B(3, 0) = 25;

	cout << endl;
	cout << "Printing B" << endl;
	B.print(4);

	cout << endl;
	cout << "Printing Initial Approximation x0" << endl;

	matrixf initial_approx(4, 1);
	initial_approx(0, 0) = 100;
	initial_approx(1, 0) = 100;
	initial_approx(2, 0) = 100;
	initial_approx(3, 0) = 100;
	
	cout << endl;
	initial_approx.print(4);
	cout << endl;
	
	cout << "Solvng the system A x = b with A = M, b = B"<<endl<<" tolerance = 0.1, MAX_ITERATIONS = 200" << endl;
	M.transpose();

	matrixf SOL = M.Gauss_Seidel(B, initial_approx,  0.1, 200);

	cout << endl;
	cout << "solution" << endl;
	SOL.print(4);


	cout << endl;
	cout << "=================================================" << endl;
	cout << "Lower triangular system test." << endl;

	matrixf L(3, 3);
	L(0, 0) = 2;
	L(0, 1) = 0;
	L(0, 2) = 0;

	L(1, 0) = 1;
	L(1, 1) = 4;
	L(1, 2) = 0;

	L(2, 0) = 7;
	L(2, 1) = -3;
	L(2, 2) = 5;

	cout << endl;
	cout << "Printing L" << endl;
	L.print(4);

	cout << endl;


	matrixf y(3, 1);

	matrixf b3(3, 1);

	b3(0, 0) = 14;
	b3(1, 0) = -101;
	b3(2, 0) = 155;

	cout << endl;
	cout << "Printing B" << endl;
	b3.print(4);

	cout << endl;


	M.Solve_Lower_TriangularSystem(L, y, b3);
	cout << "printing solution" << endl;
	y.print(4);


	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Upper Triangular Solver" << endl;
	
	L.transpose();
	
	cout << endl;
	cout << "Printing U" << endl;
	L.print(4);
	
	cout << endl;
	
	matrixf x(3, 1);

	M.Solve_Upper_TriangularSystem(L, x, y);


	cout << "printing solution" << endl;
	x.print(4);

	cout << endl;
	cout << "=================================================" << endl;
	cout << "testing cholesky" << endl;

	matrixf A(3, 3);
	A(0, 0) = 4;
	A(0, 1) = 2;
	A(0, 2) = 14;

	A(1, 0) = 2;
	A(1, 1) = 17;
	A(1, 2) = -5;

	A(2, 0) = 14;
	A(2, 1) = -5;
	A(2, 2) =83;
	
	cout << endl;
	cout << "Printing A" << endl;
	A.print(4);
	cout << endl;

	cout << "Printing b" << endl;
	b3.print(4);
	cout << endl;

	cout << "To Solve Ax = b" << endl;
	matrixf c = A.Cholesky(b3);
	c.print(4);
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing matrixd Inversion via Gauss-Jordan Elimination " << endl;

	matrixd G(3, 3);
	G(0, 0) = -1;
	G(0, 1) = 1;
	G(0, 2) = 2;

	G(1, 0) = 3;
	G(1, 1) = -1;
	G(1, 2) = 1;

	G(2, 0) = -1;
	G(2, 1) = 3;
	G(2, 2) = 4;

	cout << endl;
	cout << "printing G" << endl;
	G.print(4);
	

	matrixd old_G = G;
	matrixd G2 = G.Gauss_Jordan();
	cout << endl;
	cout << "printing G2" << endl;
	G2.print(4);
	
	cout << endl;
	cout << "printing G" << endl;
	G.print(4);
	cout << endl;


	cout << endl;
	cout << "printing G2 * old_G" << endl;
	matrixd what = old_G * G2;
	what.print(4);
	cout << endl;
	cout << "Close enough ignoring numbers below 1e-010" << endl;

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Doolittle LU decomposition" << endl;
	cout << endl;

	matrixf L_2(3, 3);
	L_2.Identity();

	matrixf U_2(3, 3);

	matrixf M_(3, 3);
	M_(0, 0) =3;
	M_(0, 1) =5;
	M_(0, 2) =2;

	M_(1, 0) =0;
	M_(1, 1) =8;
	M_(1, 2) = 2;
	M_(2, 0) = 6;
	M_(2, 1) = 2;
	M_(2, 2) =8;

	cout << "Printing M_" << endl;
	M_.print(4);
	cout << endl;

	M_.LU_Decomposition_Doolittle(L_2, U_2);

	cout << endl;
	cout << "Printing L" << endl;
	L_2.print(4);


	cout << endl;
	cout << "Printing U" << endl;
	U_2.print(4);
	cout << endl;
	matrixf LU = L_2 * U_2;

	cout << "Printing L * U" << endl;
	LU.print(4);

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Doolittle system solver" << endl;
	cout << endl;

	matrixf B_lu(3, 1);
	B_lu(0, 0) = 8;
	B_lu(1, 0) = -7;
	B_lu(2, 0) = 26;

	matrixf X_lu = M_.Solve_System_LU(L_2, U_2, B_lu);
	cout << "1. Using previously computed LU matrices"<< endl;
	cout << "Printing solution" << endl;
	X_lu.print(4);

	B_lu(0, 0) = 8;
	B_lu(1, 0) = -7;
	B_lu(2, 0) = 26;
	matrixf X_lu2 = M_.Solve_System_Doolittle(B_lu);
	cout << endl;
	cout << "2. Using internally computed LU matrices" << endl;
	cout << "Printing solution" << endl;
	X_lu2.print(4);

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing LU Decomposition Crout's Method" << endl;
	cout << endl;
	matrixf L_3(3, 3);
	matrixf U_3(3, 3);
	M_.LU_Decomposition_Crout(L_3, U_3);

	cout << "Printing L" << endl;
	L_3.print(4);

	cout << endl;
	cout << "Printing U" << endl;
	U_3.print(4);

	cout << endl;
	cout << "Printing L*U" << endl;
	matrixf crout = L_3 * U_3;
	crout.print(4);

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Crout system solver" << endl;
	cout << endl;
	cout << "1. Using previously computed LU matrices" << endl;

	matrixf X_crout_lu = M_.Solve_System_LU(L_3, U_3, B_lu);
	cout << "Printing Solution" << endl;
	X_crout_lu.print(4);
	
	cout << endl;
	cout << "2. Using internally computed LU matrices" << endl;
	matrixf X_crout_lu2 = M_.Solve_System_Crout(B_lu);
	cout << "Printing Solution" << endl;
	X_crout_lu2.print(4);


	matrixf A_(4,4);
	A_(0, 0) = 6;
	A_(0, 1) = 4;
	A_(0, 2) = 1;
	A_(0, 3) = 1;

	A_(1, 0) = 4;
	A_(1, 1) = 6;
	A_(1, 2) = 1;
	A_(1, 3) = 1;

	A_(2, 0) = 1;
	A_(2, 1) = 1;
	A_(2, 2) = 5;
	A_(2, 3) = 2;

	A_(3, 0) = 1;
	A_(3, 1) = 1;
	A_(3, 2) = 2;
	A_(3, 3) = 5;

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Housholder Tridiagonalization" << endl;
	cout << endl;
	cout << "printing A_" << endl;
	A_.print(4);

	A_.Householder_Tridiagonalize();


	cout << endl;
	cout << "printing A_ after procedure" << endl;
	A_.print(4);


	cout << endl;
	cout << "printing A_ after clipping to 2e-06 " << endl;
	A_.ClipToZero(2e-06);
	A_.print(4);


		
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Housholder algorithm on Non-Symmetric Matrix (Hessenburg)" << endl;
	TestHessenburg_2();

	cout << endl;
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing TestHessenburg_QR_real algorithm on Hpuseholder matrix" << endl;
	TestHessenburg_QR_real();

	// results are in accordance with 
	// http://mathfaculty.fullerton.edu/mathews/n2003/hessenberg/HessenbergMod/Links/HessenbergMod_lnk_9.html
	// so the other procedures for Hessenburg will be kept until I can prove they are useless.
	return 0;
}

