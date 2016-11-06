========================================================================<br>
    CONSOLE APPLICATION : Matrix Project Overview<br>
========================================================================<br>
<br>
Running a couple of tests on Systematic solution of linear equations.<br>
<br>
Completed Procedures<br>
1) Gaussian Elimination<br>
2) Gauss-Jordan Elimination<br>
3) Determinant of a triangular matrix<br>
4) Solve lower triangular system (back substitution)<br>
5) Solve upper triangular system (back substitution)<br>
6) LU Decomposition Doolittle<br>
7) LU Decomposition Crout<br>
8) Solve LU system (gemeric)<br>
9) LU Decomposition Cholesky<br>
10) Gauss-Siedel iteration<br>
11) Householder tridiagonalization<br>
<br>
<br>
had problems with the process of constructing a Hessenburg matrix for the general case<br>
of an n*n square matrix because I thought the Householder algorithm only worked for <br>
symmetric matrices. <br>
<br>
I implemented a method that performs a similarity transform and returns<br>
H and S from http://www.mymathlib.com/matrices/eigen/, <br>
(these are removed, I used the Krezig method, see much older versions)<br>
<br>
I found a paper<br>
http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0258255-3/S0025-5718-1969-0258255-3.pdf
that proves numerical instability of the row reduction methods.(this showed up again in the complex matrix).
<br><br>
Anyway I discovered that the Householder method can be used to transform a general N*N matrix<br>
into the Hessenburg form ... The matrix in the example is from this page <br>
http://mathfaculty.fullerton.edu/mathews/n2003/hessenberg/HessenbergMod/Links/HessenbergMod_lnk_9.html<br>
and for more information<br>
http://mathfaculty.fullerton.edu/mathews/n2003/HessenbergMod.html<br>
<br>
<br>
12) implemented QR algorithm for finding real eigen-values. <br>
13) the QR procedure is working well for real valued solutions and especially symmetric matrices<br>
14) the QR procedure solved for the complex number pairs<br>
15) There are now constructors for various static initializer types<br> 
(see examples at top for information, the max is a 6x6 matrix)<br>
16) memory optimized for real matrices<br>
17) Implemented a Complex matrix type - still experimental<br> 
18) Set up matrix inversion using LU decomposition thanks to <br>
http://mathworld.wolfram.com/MatrixInverse.html<br>
and <br>
http://www.gamedev.net/page/resources/_/technical/math-and-physics/matrix-inversion-using-lu-decomposition-r3637<br>
for completely explaining it<br>
19) Added some functions for finding a determinant up to 5x5 however the input function isn't implemented.<br>
20) Added method for newton iteration to improve bad inverses.<br>
21) first complex matrix version of the householder algorithm is implemented<br>
22) implemented QR factorisation and tested here - the output produces the same R matrix<br>
23) Reduction to Hessenburg is *exactly* the same as the python version, the solution was to use CSGN instead of SGN<br> 
24) Eigenvalues for the Hermitian matrix are *exactly* the same as the python version<br> 
25) Eigenvalues for general complex matrix are working (recognizing the fact that the Hermitian is its own conjugate transpose)
    
<img src="Untitled.jpg"></img><br>

