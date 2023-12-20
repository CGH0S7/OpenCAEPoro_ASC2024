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

class PETScBSolverPS
{
public:
    dBSRmat_ Absr;
    dvector_ b;
    dvector_ x;

    void allocate(int row, int col, int nb, int nnz, int storage_manner);
    void destroy();
    int petscsolve();
};

#ifdef __cplusplus
extern "C"
{
#endif

    int FIM_solver(const int commRoot, int myid, int num_procs, int nrow, int nb, int *rpt, int *cpt, double *val, double *rhs, double *sol);

    // int FIM_solver_p(int precondID, int is_thermal, int myid, int num_procs, int nb, int *allLower, int *allUpper, int *rpt, int *cpt, double *val, double *rhs, double *sol);
    int FIM_solver_p(int precondID, int is_thermal, int myid, int num_procs, int nb, int *allLower, int *allUpper, int *rpt, int *cpt, double *val, double *rhs, double *sol, double *mct);

    int FIM_solver_p_cpr(int is_thermal, int myid, int num_procs, int nb, int *allLower, int *allUpper, int *rpt, int *cpt, double *val, double *rhs, double *sol);

    int FIM_solver_p_msp(int is_thermal, int myid, int num_procs, int nb, int *allLower, int *allUpper, int *rpt, int *cpt, double *val, double *rhs, double *sol);

    int FIM_solver_p_bamg(int is_thermal, int myid, int num_procs, int nb, int *allLower, int *allUpper, int *rpt, int *cpt, double *val, double *rhs, double *sol);

    int fim_solver_(const int *commRoot, int *myid, int *num_procs, int *nrow, int *nb, int *rpt, int *cpt, double *val, double *rhs, double *sol);

#ifdef __cplusplus
}
#endif

#endif
