// Matrix.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "my_matrix.h"
#include <iostream>

using std::cout;
using std::endl;



int main(int argc, char* argv[])
{
	matrix m(3, 3);
	m(0, 0) = 0;
	m(0, 1) = 8;
	m(0, 2) = 2;
	
	m(1, 0) = 3;
	m(1, 1) = 5;
	m(1, 2) = 2;

	m(2, 0) = 6;
	m(2, 1) = 2;
	m(2, 2) = 8;

	matrix b(3, 1);
	b(0, 0) = -7;
	b(1, 0) = 8;
	b(2, 0) = 26;

	matrix sol = m.Gauss(b);

	sol.print();

	cout << endl;
	m.print();


	matrix m2(3, 3);
	m2(0, 0) = 4;
	m2(0, 1) = 2;
	m2(0, 2) = 14;

	m2(1, 0) = 2;
	m2(1, 1) = 17;
	m2(1, 2) = -5;

	m2(2, 0) = 14;
	m2(2, 1) = -5;
	m2(2, 2) = 83;

	matrix b2(3, 1);
	b2(0, 0) = 14;
	b2(1, 0) = -101;
	b2(2, 0) = 155;

	matrix sol2 = m2.Gauss(b2);
	cout << endl;
	sol2.print();

	cout << endl;
	m2.print();


	cout << "tesing Gauss Seidel" << endl;

	matrix M(4, 4);
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

	if (M.IsSymmetric()) cout << "M is Symmetric" << endl;
	if (M.IsPositiveDefinite()) cout << "M is positive definite" << endl;

	matrix B(4, 1);

	B(0, 0) = 50;
	B(1, 0) = 50;
	B(2, 0) = 25;
	B(3, 0) = 25;

	matrix initial_approx(4, 1);
	initial_approx(0, 0) = 93;
	initial_approx(1, 0) = 90;
	initial_approx(2, 0) = 65;
	initial_approx(3, 0) = 6;
	M.transpose();
	matrix SOL = M.Gauss_Seidel(B, initial_approx,  0.1, 200);

	cout << "solution" << endl;
	SOL.print();


	cout << endl;
	cout << "Lower triangular system test." << endl;
	matrix L(3, 3);
	L(0, 0) = 2;
	L(0, 1) = 0;
	L(0, 2) = 0;

	L(1, 0) = 1;
	L(1, 1) = 4;
	L(1, 2) = 0;

	L(2, 0) = 7;
	L(2, 1) = -3;
	L(2, 2) = 5;

	matrix y(3, 1);

	matrix b3(3, 1);

	b3(0, 0) = 14;
	b3(1, 0) = -101;
	b3(2, 0) = 155;

	M.Solve_Lower_TriangularSystem(L, y, b3);
	y.print();


	cout << endl;
	cout << "Testing Upper Triangular Solver" << endl;
	L.transpose();
	matrix x(3, 1);

	M.Solve_Upper_TriangularSystem(L, x, y);

	x.print();

	cout << endl;
	cout << "testing cholesky" << endl;

	matrix A(3, 3);
	A(0, 0) = 4;
	A(0, 1) = 2;
	A(0, 2) = 14;

	A(1, 0) = 2;
	A(1, 1) = 17;
	A(1, 2) = -5;

	A(2, 0) = 14;
	A(2, 1) = -5;
	A(2, 2) =83;
	
	cout << "Printing A" << endl;
	A.print();
	cout << "Printing b" << endl;
	b3.print();
	cout << "To Solve Ax = b" << endl;
	matrix c = A.Cholesky(b3);
	c.print();
	return 0;
}

