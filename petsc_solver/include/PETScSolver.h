#ifndef PETScSolver_h
#define PETScSolver_h

// #ifdef __cplusplus
// extern "C" {
#include "petscksp.h"
// #endif

// #ifdef __cplusplus
// }
// #endif
#define PrintLvl 0   // Print Information

#define PRE_AMG 1
#define PRE_JACOBI 2
#define PRE_BJACOBI 3
#define PRE_ILU 4
#define PRE_EU 5
#define PRE_PILUT 6

#define PS_MIN(a, b) (((a) < (b)) ? (a) : (b))

typedef struct dCSRmat_
{
    int row;
    int col;
    int nnz;
    int *IA;
    int *JA;
    double *val;
} dCSRmat_;

typedef struct dBSRmat_
{
    int ROW;
    int COL;
    int NNZ;
    int nb;
    int storage_manner; // 0: row-major order, 1: column-major order
    double *val;
    int *IA;
    int *JA;
} dBSRmat_;

typedef struct dvector_
{
    int row;
    double *val;
} dvector_;

typedef struct dBlockDiag
{
    int ROW;
    int nb;
    double *val;
} dBlockDiag;

typedef struct shellContext
{
    Mat BMat;
    Mat App;
    Mat Ass;
    int nb;
    int *lower;
    int *upper;
} shellContext;

void calBLowerUpper(int myid, int num_procs, int nrow, int nb, int *rpt, int &local_size, int *allLower, int *allUpper, int *allDisp, int *allSize);

int get_Prhs(Vec globalVec, Vec &localVec, int nBlockRows, int iStart, int nb);
int get_Srhs(Vec globalVec, Vec &localVec, int nBlockRows, int iStart, int nb);

int get_PP(dBSRmat_ *A, int Istart, int Iend, int matrixDim, Mat &localApp);
int get_SS(dBSRmat_ *A, int Istart, int Iend, int matrixDim, Mat &localAss);
int combine_PS(Vec Psol, Vec Ssol, Vec sol, int nBlockRows, int iStart, int nb);
int combine_P(Vec Psol, Vec sol, int nBlockRows, int iStart, int nb);

// PETSc shell preconditioner
PetscErrorCode precondApplyMSP(PC pc, Vec xin, Vec xout);
PetscErrorCode precondApplyCPR(PC pc, Vec xin, Vec xout);
int preSolver(Mat &A, Vec &b, Vec &x, bool zeroGuess, int solverType);

// diagonal scaling (ABF)
void smat_inv_4x4(double *A);
void smat_identity_4x4(double *A);
void smat_mul_4x4(double *A, double *B, double *C);
void smat_vec_mul_4(double *A, double *b, double *c);
void decoup_abf_4x4(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow);

void smat_inv_2x2(double *A);
void smat_identity_2x2(double *A);
void smat_mul_2x2(double *A, double *B, double *C);
void smat_vec_mul_2(double *A, double *b, double *c);
void decoup_abf_2x2(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow);

void decouple_abf(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow);
void decoup_abf_nb(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow);

// Analytical decoupling method
void decouple_anl(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int is_thermal);
// Semi-analytical decoupling method
void decouple_sem(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int is_thermal);

// decoupling method
void decoup(double *val, double *rhs, int *rpt, int *cpt, int nb, int nrow, int decoup_type, int is_thermal);
#endif
