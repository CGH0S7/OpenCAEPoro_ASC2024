#ifndef PETScBSolverPS_h
#define PETScBSolverPS_h

// #ifdef __cplusplus
// extern "C" {
// #endif

#include "PETScSolver.h"

// #ifdef __cplusplus
// }
// #endif

// #ifdef __cplusplus
// extern "C" {
// #endif

class PETScBSolverPS{
public:
    dBSRmat_ Absr;
    dvector_ b;
    dvector_ x;

    void allocate(int row, int col, int nb, int nnz, int storage_manner);
    void destroy();
    int petscsolve();

};

extern "C" int FIM_solver(const int commRoot, int myid, int num_procs, int nrow, int nb, int *rpt, int *cpt, double *val, double *rhs, double *sol);

extern "C" int fim_solver_(const int *commRoot, int *myid, int *num_procs, int *nrow, int *nb, int *rpt, int *cpt, double *val, double *rhs, double *sol);

// #ifdef __cplusplus
// }
// #endif

#endif
