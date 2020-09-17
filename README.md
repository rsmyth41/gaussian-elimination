# gaussian-elimination

Programme to calculate the inverse of an NxN square matrix using the Gaussian elimination method. Determinants can be calculated. The matrix equation Ax=b can be solved, giving the solution vector x. A print of the LU decomposition i.e. L and U matrices, L^-1 matrix and the modified vector b is given. Sanity checks included to ensure accuracy and consistency.

As a SUBROUTINE - Input of (size)N, (matrix)A and (column vector)b required.

TO DO: A check of the tolerances may be needed. Check print out of matrices/vectors.

To compile:
    gfortran gauss_elim.f95 -mcmodel=medium -O3 -o gauss_elim.x

To run:
    ./gauss_elim.x