// matrix<float>.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "my_matrix.h"
#include <iostream>

using std::cout;
using std::endl;



int main(int argc, char* argv[])
{
	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination" << endl;
	matrix<float> m(3, 3);
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
	m.print();
	
	matrix<float> b(3, 1);
	b(0, 0) = -7;
	b(1, 0) = 8;
	b(2, 0) = 26;

	cout << endl;
	cout << "Printing b" << endl;
	m.print();
	
	cout << endl;
	cout << "Solvng the system A x = b for x" <<endl;
	cout << endl;

	matrix<float> sol = m.Gauss(b);
	
	cout << "printing solution" << endl;
	
	sol.print();
	
	cout << endl;


	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination 2" << endl;

	matrix<float> m2(3, 3);
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

	m2.print();

	matrix<float> b2(3, 1);
	b2(0, 0) = 14;
	b2(1, 0) = -101;
	b2(2, 0) = 155;

	cout << endl;
	cout << "Printing b" << endl;

	b2.print();


	cout << endl;
	cout << "Solvng the system A x = b" << endl;
	
	matrix<float> sol2 = m2.Gauss(b2);
	
	cout << endl;
	cout << "printing solution (x) " << endl;
	sol2.print();


	cout << endl;
	cout << "=================================================" << endl;
	cout << "matrix<float> Multiplication A * A2" << endl;
	matrix<float> D = m * m2;
	D.print();

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;

	cout << "tesing Gauss Seidel" << endl;

	matrix<float> M(4, 4);
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
	M.print();

	cout << endl;

	if (M.IsSymmetric()) cout << "M is Symmetric" << endl;
	if (M.IsPositiveDefinite()) cout << "M is positive definite" << endl;

	matrix<float> B(4, 1);

	B(0, 0) = 50;
	B(1, 0) = 50;
	B(2, 0) = 25;
	B(3, 0) = 25;

	cout << endl;
	cout << "Printing B" << endl;
	B.print();

	cout << endl;
	cout << "Printing Initial Approximation x0" << endl;

	matrix<float> initial_approx(4, 1);
	initial_approx(0, 0) = 100;
	initial_approx(1, 0) = 100;
	initial_approx(2, 0) = 100;
	initial_approx(3, 0) = 100;
	
	cout << endl;
	initial_approx.print();
	cout << endl;
	
	cout << "Solvng the system A x = b with A = M, b = B"<<endl<<" tolerance = 0.1, MAX_ITERATIONS = 200" << endl;
	M.transpose();

	matrix<float> SOL = M.Gauss_Seidel(B, initial_approx,  0.1, 200);

	cout << endl;
	cout << "solution" << endl;
	SOL.print();


	cout << endl;
	cout << "=================================================" << endl;
	cout << "Lower triangular system test." << endl;

	matrix<float> L(3, 3);
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
	L.print();

	cout << endl;


	matrix<float> y(3, 1);

	matrix<float> b3(3, 1);

	b3(0, 0) = 14;
	b3(1, 0) = -101;
	b3(2, 0) = 155;

	cout << endl;
	cout << "Printing B" << endl;
	b3.print();

	cout << endl;


	M.Solve_Lower_TriangularSystem(L, y, b3);
	cout << "printing solution" << endl;
	y.print();


	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Upper Triangular Solver" << endl;
	
	L.transpose();
	
	cout << endl;
	cout << "Printing U" << endl;
	L.print();
	
	cout << endl;
	
	matrix<float> x(3, 1);

	M.Solve_Upper_TriangularSystem(L, x, y);


	cout << "printing solution" << endl;
	x.print();

	cout << endl;
	cout << "=================================================" << endl;
	cout << "testing cholesky" << endl;

	matrix<float> A(3, 3);
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
	A.print();
	cout << endl;

	cout << "Printing b" << endl;
	b3.print();
	cout << endl;

	cout << "To Solve Ax = b" << endl;
	matrix<float> c = A.Cholesky(b3);
	c.print();
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing matrix<double> Inversion via Gauss-Jordan Elimination " << endl;

	matrix<double> G(3, 3);
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
	G.print();
	

	matrix<double> old_G = G;
	matrix<double> G2 = G.Gauss_Jordan();
	cout << endl;
	cout << "printing G2" << endl;
	G2.print();
	
	cout << endl;
	cout << "printing G" << endl;
	G.print();
	cout << endl;


	cout << endl;
	cout << "printing G2 * old_G" << endl;
	matrix<double> what = old_G * G2;
	what.print();
	cout << endl;
	cout << "Close enough ignoring numbers below 1e-010" << endl;

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Doolittle LU decomposition" << endl;
	cout << endl;

	matrix<float> L_2(3, 3);
	L_2.Identity();

	matrix<float> U_2(3, 3);

	matrix<float> M_(3, 3);
	M_(0, 0) =3;
	M_(0, 1) =5;
	M_(0, 2) =2;

	M_(1, 0) =0;
	M_(1, 1) =8;
	M_(1, 2) = 2;
	M_(2, 0) = 6;
	M_(2, 1) = 2;
	M_(2, 2) =8;

	cout << "Prining M_" << endl;
	M_.print();
	cout << endl;

	M_.LU_Decomposition_Doolittle(L_2, U_2);

	cout << endl;
	cout << "Prining L" << endl;
	L_2.print();


	cout << endl;
	cout << "Prining U" << endl;
	U_2.print();
	cout << endl;
	matrix<float> LU = L_2 * U_2;

	cout << "Printing L * U" << endl;
	LU.print();

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Doolittle system solver" << endl;
	cout << endl;

	matrix<float> B_lu(3, 1);
	B_lu(0, 0) = 8;
	B_lu(1, 0) = -7;
	B_lu(2, 0) = 26;

	matrix<float> X_lu = M_.Solve_System_LU(L_2, U_2, B_lu);
	cout << "1. Using previously computed LU matrices"<< endl;
	cout << "Prining solution" << endl;
	X_lu.print();

	B_lu(0, 0) = 8;
	B_lu(1, 0) = -7;
	B_lu(2, 0) = 26;
	matrix<float> X_lu2 = M_.Solve_System_Doolittle(B_lu);
	cout << endl;
	cout << "2. Using internally computed LU matrices" << endl;
	cout << "Prining solution" << endl;
	X_lu2.print();

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing LU Decomposition Crout's Method" << endl;
	cout << endl;
	matrix<float> L_3(3, 3);
	matrix<float> U_3(3, 3);
	M_.LU_Decomposition_Crout(L_3, U_3);

	cout << "Printing L" << endl;
	L_3.print();

	cout << endl;
	cout << "Printing U" << endl;
	U_3.print();

	cout << endl;
	cout << "Printing L*U" << endl;
	matrix<float> crout = L_3 * U_3;
	crout.print();

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Crout system solver" << endl;
	cout << endl;
	cout << "1. Using previously computed LU matrices" << endl;

	matrix<float> X_crout_lu = M_.Solve_System_LU(L_3, U_3, B_lu);
	cout << "Printing Solution" << endl;
	X_crout_lu.print();
	
	cout << endl;
	cout << "2. Using internally computed LU matrices" << endl;
	matrix<float> X_crout_lu2 = M_.Solve_System_Crout(B_lu);
	cout << "Printing Solution" << endl;
	X_crout_lu2.print();
	return 0;
}

