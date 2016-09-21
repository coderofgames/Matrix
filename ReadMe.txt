========================================================================
    CONSOLE APPLICATION : Matrix Project Overview
========================================================================

Running a couple of tests on Systematic solution of linear equations.

Completed Procedures
1) Gaussian Elimination
2) Gauss-Jordan Elimination
3) Determinant of a triangular matrix
4) Solve lower triangular system (back substitution)
5) Solve upper triangular system (back substitution)
6) LU Decomposition Doolittle
7) LU Decomposition Crout
8) Solve LU system (gemeric)
9) LU Decomposition Cholesky
10) Gauss-Siedel iteration
11) Householder tridiagonalization


had problems with the process of constructing a Hessenburg matrix for the general case
of an n*n square matrix because I thought the Householder algorithm only worked for 
symmetric matrices. I implemented a method that performs a similarity transform and returns
H and S from http://www.mymathlib.com/matrices/eigen/, however I found a paper
http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0258255-3/S0025-5718-1969-0258255-3.pdf
that proves numerical instability of the row reduction methods. The methods I used for the 
Hessenburg showed indeterminate numbers using the example matrix. I implemented both Hessenburg 
matrix routines from the website however there may be mistakes in the code because the memory
is accessed differently ...

Anyway I discovered that the Householder method can be used to transform a general N*N matrix
into the Hessenburg form (with maybe a couple of snag cases) ... The matrix in the example is from
this page http://mathfaculty.fullerton.edu/mathews/n2003/hessenberg/HessenbergMod/Links/HessenbergMod_lnk_9.html
and for more information
http://mathfaculty.fullerton.edu/mathews/n2003/HessenbergMod.html


