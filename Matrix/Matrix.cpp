// Matrix.cpp : Defines the entry point for the console application.
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

	cout << endl;
	cout << "Printing A" << endl;
	m.print();
	
	matrix b(3, 1);
	b(0, 0) = -7;
	b(1, 0) = 8;
	b(2, 0) = 26;

	cout << endl;
	cout << "Printing b" << endl;
	m.print();
	
	cout << endl;
	cout << "Solvng the system A x = b for x" <<endl;
	cout << endl;

	matrix sol = m.Gauss(b);
	
	cout << "printing solution" << endl;
	
	sol.print();
	
	cout << endl;


	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination 2" << endl;

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

	cout << endl;
	cout << "printing A2" << endl;

	m2.print();

	matrix b2(3, 1);
	b2(0, 0) = 14;
	b2(1, 0) = -101;
	b2(2, 0) = 155;

	cout << endl;
	cout << "Printing b" << endl;

	b2.print();


	cout << endl;
	cout << "Solvng the system A x = b" << endl;
	
	matrix sol2 = m2.Gauss(b2);
	
	cout << endl;
	cout << "printing solution (x) " << endl;
	sol2.print();


	cout << endl;
	cout << "=================================================" << endl;
	cout << "Matrix Multiplication A * A2" << endl;
	matrix D = m * m2;
	D.print();

	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;

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

	cout << endl;
	cout << "Printing M" << endl;
	M.print();

	cout << endl;

	if (M.IsSymmetric()) cout << "M is Symmetric" << endl;
	if (M.IsPositiveDefinite()) cout << "M is positive definite" << endl;

	matrix B(4, 1);

	B(0, 0) = 50;
	B(1, 0) = 50;
	B(2, 0) = 25;
	B(3, 0) = 25;

	cout << endl;
	cout << "Printing B" << endl;
	B.print();

	cout << endl;
	cout << "Printing Initial Approximation x0" << endl;

	matrix initial_approx(4, 1);
	initial_approx(0, 0) = 100;
	initial_approx(1, 0) = 100;
	initial_approx(2, 0) = 100;
	initial_approx(3, 0) = 100;
	
	cout << endl;
	initial_approx.print();
	cout << endl;
	
	cout << "Solvng the system A x = b with A = M, b = B"<<endl<<" tolerance = 0.1, MAX_ITERATIONS = 200" << endl;
	M.transpose();

	matrix SOL = M.Gauss_Seidel(B, initial_approx,  0.1, 200);

	cout << endl;
	cout << "solution" << endl;
	SOL.print();


	cout << endl;
	cout << "=================================================" << endl;
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

	cout << endl;
	cout << "Printing L" << endl;
	L.print();

	cout << endl;


	matrix y(3, 1);

	matrix b3(3, 1);

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
	
	matrix x(3, 1);

	M.Solve_Upper_TriangularSystem(L, x, y);


	cout << "printing solution" << endl;
	x.print();

	cout << endl;
	cout << "=================================================" << endl;
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
	
	cout << endl;
	cout << "Printing A" << endl;
	A.print();
	cout << endl;

	cout << "Printing b" << endl;
	b3.print();
	cout << endl;

	cout << "To Solve Ax = b" << endl;
	matrix c = A.Cholesky(b3);
	c.print();
	cout << endl;
	cout << "=================================================" << endl;

	return 0;
}

