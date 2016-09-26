// matrixf.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "my_matrix.h"
#include <iostream>

using std::cout;
using std::endl;

typedef LINALG::matrixf matrixf;
typedef LINALG::matrixd matrixd;


matrixf  dat_M_LU = {
	{ 3, 5, 2 },
	{ 0, 8, 2 },
	{ 6, 2, 8 }
};
matrixf  dat_M_tri = {
	{ 1, -0.25, -0.25, 0 },
	{ -0.25, 1, 0, -0.25 },
	{ -0.25, 0, 1, -0.25 },
	{ 0, -0.25, -0.25, 1 }
};
matrixf  dat_L = {
	{ 2, 0, 0 },
	{ 1, 4, 0 },
	{ 7, -3, 5 }
};

matrixf  dat_G = {
	{ -1, 1, 2 },
	{ 3, -1, 1 },
	{ -1, 3, 4 }
};

matrixf  dat_A_1 = {
	{ 6, 4, 1, 1 },
	{ 4, 6, 1, 1 },
	{ 1, 1, 5, 2 },
	{ 1, 1, 2, 5 }
};
matrixf  dat_A = {
	{ 4, 2, 14 },
	{ 2, 17, -5 },
	{ 14, -5, 83 }
};

matrixf  dat_Hess_QR2 = {
	{ 5.0f, 1.0f, 0.0f, 4.0f, 0.0f },
	{ 1.0f, 4.0f, 2.0f, 1.0f, 3.0f },
	{ 2.0f, 2.0f, 5.0f, 4.0f, 0.0f },
	{ 0.0, 1.0, 4.0, 1.0, 3.0 },
	{ 4.0, 3.0, 0.0, 3.0, 4.0 }
};
float  dat_Hess_QR1_2[4][4] = {
	{ 6.0f, -sqrt(18.0), 0.0f, 0.0f },
	{ -sqrt(18), 7.0, sqrt(2.0), 0.0f },
	{ 0.0f, sqrt(2.0), 6.0f, 0.0f },
	{ 0.0f, 0.0, 0.0f, 3.0f }
};
matrixf dat_Hess_QR1 = dat_Hess_QR1_2;
matrixf  dat_Hess_3 = {
	{ 5.0, 4.0, 3.0, 2.0, 1.0 },
	{ 1.0, 4.0, 0.0, 3.0, 3.0 },
	{ 2.0, 0.0, 3.0, 0.0, 0.0 },
	{ 3.0, 2.0, 1.0, 2.0, 5.0 },
	{ 4.0, 2.0, 1.0, 2.0, 1.0 }
};
matrixf  dat_Hess_2 = {
	{ 6.0, 1.0, -2.0, 19.0, 4.0 },
	{ 1.0, 4.0, 2.0, 1.0, 3.0 },
	{ 7.0, 12.0, 5.0, 11.0, -1.0 },
	{ -3.0, 21.0, 4.0, 1.0, 3.0 },
	{ 4.0, 6.0, -1.0, 3.0, 7.0 }
};
matrixf  dat_Hess_0 = {
	{ 4.0, 1.0, -1.0, 2.0 },
	{ 1.0, 4.0, 1.0, -1.0 },
	{ -1.0, 1.0, 4.0f, 1.0f },
	{ 2.0f, -1.0, 1.0f, 4.0f }
};
matrixf mat_QR_2 = {
	{ 4, 2, 2, 1 },
	{ 2, -3, 1, 1 },
	{ 2, 1, 3, 1 },
	{ 1, 1, 1, 2 }
};

matrixf mat_QR_1 = {
	{ 1, 3, 0, 0 },
	{ 3, 2, 1, 0 },
	{ 0, 1, 3, 4 },
	{ 0, 0, 4, 1 }
};


matrixf list_of_matrices[] = { mat_QR_1, mat_QR_2, dat_Hess_0, dat_Hess_2, dat_Hess_3, dat_Hess_QR1, dat_Hess_QR2, dat_A, dat_A_1, dat_G, dat_M_tri };












void TestingProcedure_Householder_Hessenburg_QR(matrixf &H_)
{
	cout << endl;
	cout << "===================================================" << endl;
	// using the wisdom from
	// http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0258255-3/S0025-5718-1969-0258255-3.pdf
	// I have decided to only use Housedholder methods
	// this matrix shows the evidence of the stability issue ...


	cout << endl;
	cout << "Printing the initial matrix" << endl;

	H_.print(3);

	cout << endl;

	matrixf H__3 = H_;

	cout << endl;
	cout << "Applying Householder Algorithm (required before QR algorithm" << endl;
	cout << endl;
	
	std::cin.clear();
	char sel1;
	std::cin.get(sel1);

	
//	if (sel1 == 'y' || sel1 == 'Y')
//	{
		H__3.Householder_Tridiagonalize();

		cout << endl;
		cout << "Printing the matrix afer Householder reduction" << endl;
		cout << endl;

		H__3.print(4);
//	}

	cout << endl;
	
	cout << endl;
	cout << "Apply QR Algorithm ? (Y/N) " << endl;
	cout << endl;

	std::cin.clear();
	char sel;
	std::cin.get( sel);

	if (sel == 'y' || sel == 'Y')
	{
		matrixf eigen_values(H__3.NumColumns(), 2);

		H__3.QR_algorithm(eigen_values);

		cout << endl;
		cout << "Printing the matrix afer QR algorithm" << endl;
		cout << endl;

		H__3.print(4);

		cout << endl;
		cout << "printing the eigen values" << endl;
		cout << endl;

		eigen_values.print(3);
	}

	cout << endl;
	cout << endl;
}



void TestGaussElimination_1()
{
	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination" << endl;
	matrixf m = dat_M_LU;
	
	m.SwapRow(0, 1); // these were swapped

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
	cout << "Solvng the system A x = b for x" << endl;
	cout << endl;

	matrixf sol = m.Gauss_Elimination(b);

	cout << "printing solution" << endl;

	sol.print(4);

	cout << endl;
}

void TestGaussElimination_2()
{
	cout << "=================================================" << endl;
	cout << "Testing Gauss Elimination 2" << endl;

	matrixf m2 = dat_A;

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
}

void TestGaussSeidel()
{
	cout << endl;
	cout << "=================================================" << endl;
	cout << endl;

	cout << "tesing Gauss Seidel" << endl;

	matrixf M = dat_M_tri;

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

	cout << "Solvng the system A x = b with A = M, b = B" << endl << " tolerance = 0.1, MAX_ITERATIONS = 200" << endl;
	M.transpose();

	matrixf SOL = M.Gauss_Seidel(B, initial_approx, 0.1, 200);

	cout << endl;
	cout << "solution" << endl;
	SOL.print(4);
}




void Test_Doolittle_LU_System()
{
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Doolittle LU decomposition" << endl;
	cout << endl;

	matrixf L_2(3, 3);
	L_2.Identity();

	matrixf U_2(3, 3);





	matrixf M_ = dat_M_LU;



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
	cout << "1. Using previously computed LU matrices" << endl;
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

}
void TestCrout_LU_Method()
{
	cout << "=================================================" << endl;
	cout << "Testing LU Decomposition Crout's Method" << endl;
	cout << endl;
	matrixf L_3(3, 3);
	matrixf U_3(3, 3);

	matrixf M_ = dat_M_LU;

	matrixf A3 = M_;
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

	matrixf B_lu(3, 1);
	B_lu(0, 0) = 8;
	B_lu(1, 0) = -7;
	B_lu(2, 0) = 26;

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
}
void Testing_Determinant()
{
	matrixf M_ = dat_M_LU;

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Det" << endl;

	matrixf M_2 = M_;
	cout << M_.Det_3x3(0, 0) << "    " << M_2.Determinant();
	cout << endl;

}



void Test_TriangularSystemSolvers()
{
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Lower triangular system test." << endl;


	matrixf L = dat_L;


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



	matrixf M = dat_M_tri;


	M(3, 0) = 0;
	M(3, 1) = -0.25;
	M(3, 2) = -0.25;
	M(3, 3) = 1;
	M.Solve_Lower_TriangularSystem(L, y, b3);
	cout << "printing solution" << endl;
	y.print(4);
}

void TestUpperTriangularSolver(){

	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Upper Triangular Solver" << endl;

	matrixf L = dat_L;
	matrixf M = dat_M_tri;

	L.transpose();

	cout << endl;
	cout << "Printing U" << endl;
	L.print(4);

	cout << endl;

	matrixf y(3, 1);

	matrixf b3(3, 1);

	b3(0, 0) = 14;
	b3(1, 0) = -101;
	b3(2, 0) = 155;

	matrixf x(3, 1);

	M.Solve_Upper_TriangularSystem(L, x, y);


	cout << "printing solution" << endl;
	x.print(4);
}
void TestCholesky(){


	cout << endl;
	cout << "=================================================" << endl;
	cout << "testing cholesky" << endl;

	matrixf y(3, 1);

	matrixf b3(3, 1);

	b3(0, 0) = 14;
	b3(1, 0) = -101;
	b3(2, 0) = 155;

	matrixf A = dat_A;
	matrixf A2 = A;
	
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


}


void Testing_GaussJordan_Elimination()
{
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing matrixd Inversion via Gauss-Jordan Elimination " << endl;


	matrixf G = dat_G;


	cout << endl;
	cout << "printing G" << endl;
	G.print(4);


	matrixf old_G = G;
	matrixf G2 = G.Gauss_Jordan();
	cout << endl;
	cout << "printing G2" << endl;
	G2.print(4);

	cout << endl;
	cout << "printing G" << endl;
	G.print(4);
	cout << endl;


	cout << endl;
	cout << "printing G2 * old_G" << endl;
	matrixf what = old_G * G2;
	what.print(4);
	cout << endl;
	cout << "Close enough ignoring numbers below 1e-010" << endl;
}


void Testing_Householder_Tridiagonalize()
{
	matrixf A_ = dat_A_1;
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

}

void Test_Transposition()
{

	matrixf m = dat_M_tri;
	m.print(2);
	m.transpose();
	cout << endl << endl;
	m.print(2);
	cout << endl << endl;

	float dat__[4][3] = {
		{ 2, 1, 3 },
		{ 4, 5, 6 },
		{ 6, 5, 4 },
		{ 3, 2, 1 }
	};

	matrixf b = {
		{ 2, 1, 3 },
		{ 4, 5, 6 },
		{ 6, 5, 4 },
		{ 3, 2, 1 }
	};
	//for (int i = 0; i < 4; i++)
	//	for (int j = 0; j < 3; j++)
	//		b(i, j) = dat__[i][j];

	b.print(2);
	//b.transpose();
	cout << endl << endl;


	//	matrixf c = m * b;

	//	c.print(2);

	cout << endl << endl;
	b.transpose();
	b.print(2);
	cout << endl << endl;

	matrixf d = m*b;

	d.print(2);
	cout << endl << endl;
}

int main(int argc, char* argv[])
{



	//TestGaussElimination_1();
	//TestGaussElimination_2();

//	TestCholesky();
//	TestCrout_LU_Method();

//	TestGaussSeidel();




		
	cout << endl;
	cout << "=================================================" << endl;
	cout << "Testing Housholder algorithm on Non-Symmetric Matrix (Hessenburg)" << endl;
	
	Testing_Householder_Tridiagonalize();

	

	int matrix_selection;

	cout << "printing list of matrices " << endl;
	for (int i = 0; i < 11; i++)
	{
		cout << (char)('A' + i) << " = " << endl << endl;
		
		list_of_matrices[i].print(2);
		cout << endl << endl;
	}

	cout << "Select a matrix from the above (A - K)" << endl;

	char sel;
	std::cin.clear();

	std::cin.get(sel);

	if (sel >= 'A' && sel <= 'A' + 11)
	{
		cout << "You Selected: " << sel << endl;
		matrix_selection = sel - 'A';

		//list_of_matrices[matrix_selection].print(2);
		cout << "Testing Housholder algorithm and QR algorithm "<< endl;

		TestingProcedure_Householder_Hessenburg_QR(list_of_matrices[matrix_selection]);
	}









//	cout << endl;
//	cout << "=================================================" << endl;
//	cout << "matrixf Multiplication A * A2" << endl;


//	Test_LU_System();

	//Test_TriangularSystemSolvers();




	//Test_Transposition();
		return 0;
}

