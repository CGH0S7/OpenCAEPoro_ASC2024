#include <cstdio>
#include <iostream>
#include <string>
#include <mpi.h>

#include "cblas.h"
#include "lapacke.h"

// #include "PETScBSolverPS.h"


using namespace std;

// int main(int argc, char* argv[])
// {
//     // int myrank, numproc, blockdim, allBegin, allEnd, iA, jA, A, b, x;
//     // FIM_solver_p(myrank, numproc, blockdim, allBegin, allEnd, iA, jA, A, b, x);
    
// }

// extern "C" {
//     // LU decomoposition of a general matrix
//     void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

//     // generate inverse of a matrix given its LU decomposition
//     void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
// }

void inverse(double* A, int N)
{
    int *IPIV = new int[N];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    delete[] IPIV;
    delete[] WORK;
}

int main(){

    double A [4] = {1,2,3,4};

    inverse(A, 2);

    printf("%f %f\n", A[0], A[1]);
    printf("%f %f\n", A[2], A[3]);

    double B [3*3] = {
        1,0,0,
        0,1,0,
        0,0,0.5
    };

    inverse(B, 3);

    printf("%f %f %f\n", B[0], B[1], B[2]);
    printf("%f %f %f\n", B[3], B[4], B[5]);
    printf("%f %f %f\n", B[6], B[7], B[8]);

    return 0;
}