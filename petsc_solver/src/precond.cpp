#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "PETScSolver.h"

using namespace std;

PetscErrorCode precondApplyMSP(PC pc, Vec xin, Vec xout)
{

    FILE *fp;
    // Get A from PC
    shellContext *context;
    PCShellGetContext(pc, (void **)&context);
    Mat Atmp = (*context).BMat;
    Mat App = (*context).App;
    Mat Ass = (*context).Ass;

    int *allLower = (*context).lower;
    int *allUpper = (*context).upper;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int iStart = allLower[myid];
    int iEnd = allUpper[myid];
    int local_size = iEnd - iStart + 1;

    int nb = (*context).nb;

    double start, finish;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    // Solve App*Xp = Bp
    Vec Xp;
    Vec Bp;
    get_Prhs(xin, Bp, local_size, iStart, nb);
    get_Prhs(xout, Xp, local_size, iStart, nb);

    preSolver(App, Bp, Xp, true, PRE_AMG);

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fp = fopen("linear_solver.log", "a");
        fprintf(fp, "App presolver time: %lf Sec\n", finish - start);
    }

    start = MPI_Wtime();
    // Solve Ass*Xs = Bs
    Vec Xs;
    Vec Bs;
    get_Srhs(xin, Bs, local_size, iStart, nb);
    get_Srhs(xout, Xs, local_size, iStart, nb);

    // preSolver( Ass, Bs, Xs, true, PRE_BJACOBI);
    preSolver(Ass, Bs, Xs, true, PRE_BJACOBI);

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0)
        fprintf(fp, "Ass PRE_BJACOBI time: %lf Sec\n", finish - start);

    start = MPI_Wtime();
    // Combine Psol and Ssol to sol
    combine_PS(Xp, Xs, xout, local_size, iStart, nb);

    // Smooth xout
    // preSolver( Atmp, xin, xout, false, PRE_BJACOBI);
    preSolver(Atmp, xin, xout, false, PRE_BJACOBI);

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fprintf(fp, "Atmp PRE_BJACOBI time: %lf Sec\n", finish - start);
        fclose(fp);
    }

    PetscErrorCode ierr = VecDestroy(&Bp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Xp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Bs);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Xs);
    CHKERRQ(ierr);

    

    return 0;
}

PetscErrorCode precondApplyCPR(PC pc, Vec xin, Vec xout)
{

    FILE *fp;
    // Get A from PC
    shellContext *context;
    PCShellGetContext(pc, (void **)&context);
    Mat Atmp = (*context).BMat;
    Mat App = (*context).App;
    // Mat Ass = (*context).Ass;

    int *allLower = (*context).lower;
    int *allUpper = (*context).upper;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int iStart = allLower[myid];
    int iEnd = allUpper[myid];
    int local_size = iEnd - iStart + 1;

    int nb = (*context).nb;

    double start, finish;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    // Solve App*Xp = Bp
    Vec Xp;
    Vec Bp;
    get_Prhs(xin, Bp, local_size, iStart, nb);
    get_Prhs(xout, Xp, local_size, iStart, nb);

    preSolver(App, Bp, Xp, true, PRE_AMG);

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fp = fopen("linear_solver.log", "a");
        fprintf(fp, "App PRE_AMG time: %lf Sec\n", finish - start);
    }

    // start = MPI_Wtime();
    // Solve Ass*Xs = Bs
    // Vec Xs;
    // Vec Bs;
    // get_Srhs(xin, Bs, local_size, iStart, nb);
    // get_Srhs(xout, Xs, local_size, iStart, nb);

    // // preSolver( Ass, Bs, Xs, true, PRE_BJACOBI);
    // preSolver(Ass, Bs, Xs, true, PRE_BJACOBI);

    // finish = MPI_Wtime();
    // MPI_Barrier(MPI_COMM_WORLD);
    // if (myid == 0)
    //     fprintf(fp, "Ass PRE_BJACOBI time: %lf Sec\n", finish - start);

    start = MPI_Wtime();
    // prolongation Psol to sol
    combine_P(Xp, xout, local_size, iStart, nb);

    // Smooth xout
    // preSolver( Atmp, xin, xout, false, PRE_BJACOBI);
    preSolver(Atmp, xin, xout, false, PRE_BJACOBI); // PRE_PILUT, PRE_ILU, PRE_EU

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fprintf(fp, "Atmp PRE_BJACOBI time: %lf Sec\n", finish - start);
        fclose(fp);
    }

    PetscErrorCode ierr = VecDestroy(&Bp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Xp);
    CHKERRQ(ierr);
    // ierr = VecDestroy(&Bs);
    // CHKERRQ(ierr);
    // ierr = VecDestroy(&Xs);
    // CHKERRQ(ierr);

    return 0;
}

// ********************************************
// parameter:
//     myid: 当前进程号
//     num_procs: 进程总数
//     nrow、nb、rpt: BSR格式矩阵
//     allLower: 每一个分段的起始位置
//     allUpper: 每一个分段的终止位置
//     allDisp: 每一个分段的偏移量
//     allSize: 每一个分段的分段长度
// ********************************************

void calBLowerUpper(int myid, int num_procs, int nrow, int nb, int *rpt, int &local_size, int *allLower, int *allUpper, int *allDisp, int *allSize)
{
    int i;

    // 计算每个进程上分到的块矩阵的行数
    local_size = nrow / num_procs;

    // 因无法整除而余下的行数
    int tmpExtra = nrow - local_size * num_procs;

    // 定义数组，用来表示每一个分段的起始位置、终止位置、偏移量、分段长度
    for (i = 0; i < num_procs; i++)
    {
        allLower[i] = local_size * i;
        allLower[i] += PS_MIN(i, tmpExtra);

        allUpper[i] = local_size * (i + 1);
        allUpper[i] += PS_MIN(i + 1, tmpExtra) - 1;

        allDisp[i] = rpt[allLower[i]];
        allSize[i] = rpt[allUpper[i] + 1] - rpt[allLower[i]];
    }

    local_size = allUpper[myid] - allLower[myid] + 1;
}

int combine_P(Vec Psol, Vec sol, int nBlockRows, int iStart, int nb)
{
    // VecRestoreSubVector
    int *pIndex = (int *)malloc(nBlockRows * sizeof(int));
    for (int i = 0; i < nBlockRows; i++)
    {
        pIndex[i] = i + iStart;
    }
    double *Pvec = (double *)malloc(nBlockRows * sizeof(double));
    VecGetValues(Psol, nBlockRows, pIndex, Pvec);

    for (int i = 0; i < nBlockRows; i++)
    {
        pIndex[i] = (i + iStart) * nb;
    }
    VecSetValues(sol, nBlockRows, pIndex, Pvec, INSERT_VALUES);

    PetscErrorCode ierr = VecAssemblyBegin(sol);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(sol);
    CHKERRQ(ierr);

    free(pIndex);
    free(Pvec);
    return 0;
}

int combine_PS(Vec Psol, Vec Ssol, Vec sol, int nBlockRows, int iStart, int nb)
{
    // VecRestoreSubVector
    int *pIndex = (int *)malloc(nBlockRows * sizeof(int));
    for (int i = 0; i < nBlockRows; i++)
    {
        pIndex[i] = i + iStart;
    }
    double *Pvec = (double *)malloc(nBlockRows * sizeof(double));
    VecGetValues(Psol, nBlockRows, pIndex, Pvec);

    for (int i = 0; i < nBlockRows; i++)
    {
        pIndex[i] = (i + iStart) * nb;
    }
    VecSetValues(sol, nBlockRows, pIndex, Pvec, INSERT_VALUES);

    PetscErrorCode ierr = VecAssemblyBegin(sol);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(sol);
    CHKERRQ(ierr);

    //-------------------------------------------------------------------
    int numS = nBlockRows * (nb - 1);
    int *sIndex = (int *)malloc(numS * sizeof(int));
    for (int i = 0; i < numS; i++)
    {
        sIndex[i] = i + iStart * (nb - 1);
    }
    double *Svec = (double *)malloc(numS * sizeof(double));
    VecGetValues(Ssol, numS, sIndex, Svec);

    int index = 0;

    for (int i = 0; i < nBlockRows; i++)
    {
        for (int j = 1; j < nb; j++)
        {
            sIndex[index] = (iStart + i) * nb + j;
            index++;
        }
    }

    VecSetValues(sol, numS, sIndex, Svec, INSERT_VALUES);

    ierr = VecAssemblyBegin(sol);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(sol);
    CHKERRQ(ierr);

    free(pIndex);
    free(Pvec);
    free(sIndex);
    free(Svec);
    return 0;
}

int get_Prhs(Vec globalVec, Vec &localVec, int nBlockRows, int iStart, int nb)
{
    IS is;
    PetscErrorCode ierr = ISCreateStride(PETSC_COMM_WORLD, nBlockRows, iStart * nb, nb, &is);
    CHKERRQ(ierr);
    ierr = VecGetSubVector(globalVec, is, &localVec);
    CHKERRQ(ierr);

    ISDestroy(&is);
    return 0;
}

int get_Srhs(Vec globalVec, Vec &localVec, int nBlockRows, int iStart, int nb)
{
    int *idx = (int *)malloc(sizeof(int) * nBlockRows * (nb - 1));
    int index = 0;
    for (int i = 0; i < nBlockRows; i++)
    {
        for (int j = 1; j < nb; j++)
        {
            idx[index] = iStart + i * nb + j;
            index++;
        }
    }
    IS is;
    PetscErrorCode ierr = ISCreateGeneral(PETSC_COMM_WORLD, nBlockRows * (nb - 1), idx, PETSC_COPY_VALUES, &is);
    CHKERRQ(ierr);
    ierr = VecGetSubVector(globalVec, is, &localVec);
    CHKERRQ(ierr);

    ISDestroy(&is);
    free(idx);
    return 0;
}

int get_PP(dBSRmat_ *A, int Istart, int Iend, int matrixDim, Mat &localApp)
{
    FILE *fp;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double start, finish;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    dCSRmat_ App;

    const int row = A->ROW;
    const int col = A->COL;
    const int nnz = A->NNZ;
    const int nb = A->nb;
    const int nb2 = nb * nb;
    double *val = A->val;
    int *IA = A->IA;
    int *JA = A->JA;

    App.row = row;
    App.col = col;
    App.nnz = nnz;

    App.IA = (int *)malloc(sizeof(int) * (row + 1));
    App.JA = (int *)malloc(sizeof(int) * nnz);
    App.val = (double *)malloc(sizeof(double) * nnz);

    double *Pval = App.val;

    memcpy(App.IA, IA, (row + 1) * sizeof(int));
    memcpy(App.JA, JA, nnz * sizeof(int));

    int i, j, Ii;
    for (i = 0; i < nnz; ++i)
    {
        Pval[i] = val[i * nb2];
    }

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fp = fopen("linear_solver.log", "a");
        fprintf(fp, "get_PP_fasp time: %lf Sec\n", finish - start);
    }
    //-------------------------------------------

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    Mat matApp;

    int blockSize = 1;
    int rowWidth = App.row * blockSize;
    int nBlockRows = App.row;

    int *rpt = App.IA;
    int *cpt = App.JA;
    val = App.val;

    int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
    int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

    for (i = 0; i < nBlockRows; i++)
    {
        nDCount[i] = 0;
        for (j = rpt[i]; j < rpt[i + 1]; j++)
        {
            if (cpt[j] >= Istart && cpt[j] <= Iend)
            {
                nDCount[i]++;
            }
        }
    }

    for (i = 0; i < nBlockRows; i++)
    {
        nNDCount[i] = rpt[i + 1] - rpt[i] - nDCount[i];
    }

    int dim = matrixDim / nb;
    PetscErrorCode ierr = MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, dim, dim, 0, nDCount, 0, nNDCount, &matApp);
    CHKERRQ(ierr);

    ierr = MatSetFromOptions(matApp);
    CHKERRQ(ierr);

    ierr = MatSetUp(matApp);
    CHKERRQ(ierr);

    double *valpt = val;
    int *globalx = (int *)malloc(blockSize * sizeof(int));
    int *globaly = (int *)malloc(blockSize * sizeof(int));
    int b2 = blockSize * blockSize;

    for (Ii = 0; Ii < nBlockRows; Ii++)
    {
        for (i = 0; i < blockSize; i++)
        {
            globalx[i] = (Ii + Istart) * blockSize + i;
        }

        for (i = rpt[Ii]; i < rpt[Ii + 1]; i++)
        {

            for (j = 0; j < blockSize; j++)
            {
                globaly[j] = cpt[i] * blockSize + j;
            }
            ierr = MatSetValues(matApp, blockSize, globalx, blockSize, globaly, valpt, INSERT_VALUES);
            CHKERRQ(ierr);
            valpt += b2;
        }
    }

    ierr = MatAssemblyBegin(matApp, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matApp, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    localApp = matApp;

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fprintf(fp, "get_PP_petsc time: %lf Sec\n", finish - start);
        fclose(fp);
    }

    free(App.IA);
    free(App.JA);
    free(App.val);
    return 0;
}

int get_SS(dBSRmat_ *A, int Istart, int Iend, int matrixDim, Mat &localAss)
{
    FILE *fp;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double start, finish;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    dBSRmat_ Ass;

    const int ROW = A->ROW;
    const int COL = A->COL;
    const int NNZ = A->NNZ;
    const int nb = A->nb;
    const int nb2 = nb * nb;
    double *val = A->val;
    int *IA = A->IA;
    int *JA = A->JA;

    int nbs = nb - 1;
    int nbs2 = nbs * nbs;

    Ass.ROW = ROW;
    Ass.COL = COL;
    Ass.NNZ = NNZ;
    Ass.nb = nbs;

    Ass.IA = (int *)malloc(sizeof(int) * (ROW + 1));
    Ass.JA = (int *)malloc(sizeof(int) * NNZ);
    Ass.val = (double *)malloc(sizeof(double) * nbs2 * NNZ);
    double *Sval = Ass.val;

    memcpy(Ass.IA, IA, (ROW + 1) * sizeof(int));
    memcpy(Ass.JA, JA, NNZ * sizeof(int));

    int i, j, k, Ii;
    for (i = 0; i < NNZ; ++i)
    {
        for (j = 0; j < nbs; ++j)
        {
            for (k = 0; k < nbs; ++k)
            {
                Sval[i * nbs2 + j * nbs + k] = val[i * nb2 + (j + 1) * nb + k + 1];
            }
        }
    }
    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fp = fopen("linear_solver.log", "a");
        fprintf(fp, "get_SS_fasp time: %lf Sec\n", finish - start);
    }
    
    //-------------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    Mat matAss;

    int blockSize = Ass.nb;
    int rowWidth = Ass.ROW * blockSize;
    int nBlockRows = Ass.ROW;

    int *rpt = Ass.IA;
    int *cpt = Ass.JA;
    val = Ass.val;

    int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
    int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

    for (i = 0; i < nBlockRows; i++)
    {
        nDCount[i] = 0;
        for (j = rpt[i]; j < rpt[i + 1]; j++)
        {
            if (cpt[j] >= Istart && cpt[j] <= Iend)
            {
                nDCount[i]++;
            }
        }
    }

    for (i = 0; i < nBlockRows; i++)
    {
        nNDCount[i] = rpt[i + 1] - rpt[i] - nDCount[i];
    }

    int dim = matrixDim / nb * (nb - 1);
    PetscErrorCode ierr = MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, dim, dim, 0, nDCount, 0, nNDCount, &matAss);
    CHKERRQ(ierr);

    ierr = MatSetFromOptions(matAss);
    CHKERRQ(ierr);

    ierr = MatSetUp(matAss);
    CHKERRQ(ierr);

    double *valpt = val;
    int *globalx = (int *)malloc(blockSize * sizeof(int));
    int *globaly = (int *)malloc(blockSize * sizeof(int));
    int b2 = blockSize * blockSize;

    for (Ii = 0; Ii < nBlockRows; Ii++)
    {
        for (i = 0; i < blockSize; i++)
        {
            globalx[i] = (Ii + Istart) * blockSize + i;
        }

        for (i = rpt[Ii]; i < rpt[Ii + 1]; i++)
        {

            for (j = 0; j < blockSize; j++)
            {
                globaly[j] = cpt[i] * blockSize + j;
            }
            ierr = MatSetValues(matAss, blockSize, globalx, blockSize, globaly, valpt, INSERT_VALUES);
            CHKERRQ(ierr);
            valpt += b2;
        }
    }

    ierr = MatAssemblyBegin(matAss, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matAss, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    localAss = matAss;

    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0 && PrintLvl > 0){
        fprintf(fp, "get_SS_petsc time: %lf Sec\n", finish - start);
        fclose(fp);
    }

    free(Ass.IA);
    free(Ass.JA);
    free(Ass.val);

    return 0;
}

/**
 * A =
 * [PP  PN1  PN2  ... PT ]
 * [N1P N1N1 N1N2 ... N1T]
 * [N2P N2N1 N2N2 ... N2T]
 *         ...
 * [TP  TN1  TN2  ... TT]
 *
 *
 * */
// int get_SS_thermal(dBSRmat_ *A, int Istart, int Iend, int matrixDim, Mat &localAss)
// {
//     FILE *fp;
//     fp = fopen("linear_solver.log", "a");
//     int myid;
//     MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//     double start, finish;
//     MPI_Barrier(MPI_COMM_WORLD);
//     start = MPI_Wtime();

//     dBSRmat_ Ass;

//     const int ROW = A->ROW;
//     const int COL = A->COL;
//     const int NNZ = A->NNZ;
//     const int nb = A->nb;
//     const int nb2 = nb * nb;
//     double *val = A->val;
//     int *IA = A->IA;
//     int *JA = A->JA;

//     int nbs = nb - 1;
//     int nbs2 = nbs * nbs;

//     Ass.ROW = ROW;
//     Ass.COL = COL;
//     Ass.NNZ = NNZ;
//     Ass.nb = nbs;

//     Ass.IA = (int *)malloc(sizeof(int) * (ROW + 1));
//     Ass.JA = (int *)malloc(sizeof(int) * NNZ);
//     Ass.val = (double *)malloc(sizeof(double) * nbs2 * NNZ);
//     double *Sval = Ass.val;

//     memcpy(Ass.IA, IA, (ROW + 1) * sizeof(int));
//     memcpy(Ass.JA, JA, NNZ * sizeof(int));

//     int i, j, k, Ii;
//     for (i = 0; i < NNZ; ++i)
//     {
//         for (j = 0; j < nbs; ++j)
//         {
//             for (k = 0; k < nbs; ++k)
//             {
//                 Sval[i * nbs2 + j * nbs + k] = val[i * nb2 + (j + 1) * nb + k + 1];
//             }
//         }
//     }
//     finish = MPI_Wtime();
//     MPI_Barrier(MPI_COMM_WORLD);
//     if (myid == 0)
//         fprintf(fp, "get_SS_fasp time: %lf Sec\n", finish - start);

//     //-------------------------------------------
//     MPI_Barrier(MPI_COMM_WORLD);
//     start = MPI_Wtime();

//     Mat matAss;

//     int blockSize = Ass.nb;
//     int rowWidth = Ass.ROW * blockSize;
//     int nBlockRows = Ass.ROW;

//     int *rpt = Ass.IA;
//     int *cpt = Ass.JA;
//     val = Ass.val;

//     int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
//     int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

//     for (i = 0; i < nBlockRows; i++)
//     {
//         nDCount[i] = 0;
//         for (j = rpt[i]; j < rpt[i + 1]; j++)
//         {
//             if (cpt[j] >= Istart && cpt[j] <= Iend)
//             {
//                 nDCount[i]++;
//             }
//         }
//     }

//     for (i = 0; i < nBlockRows; i++)
//     {
//         nNDCount[i] = rpt[i + 1] - rpt[i] - nDCount[i];
//     }

//     int dim = matrixDim / nb * (nb - 1);
//     PetscErrorCode ierr = MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, dim, dim, 0, nDCount, 0, nNDCount, &matAss);
//     CHKERRQ(ierr);

//     ierr = MatSetFromOptions(matAss);
//     CHKERRQ(ierr);

//     ierr = MatSetUp(matAss);
//     CHKERRQ(ierr);

//     double *valpt = val;
//     int *globalx = (int *)malloc(blockSize * sizeof(int));
//     int *globaly = (int *)malloc(blockSize * sizeof(int));
//     int b2 = blockSize * blockSize;

//     for (Ii = 0; Ii < nBlockRows; Ii++)
//     {
//         for (i = 0; i < blockSize; i++)
//         {
//             globalx[i] = (Ii + Istart) * blockSize + i;
//         }

//         for (i = rpt[Ii]; i < rpt[Ii + 1]; i++)
//         {

//             for (j = 0; j < blockSize; j++)
//             {
//                 globaly[j] = cpt[i] * blockSize + j;
//             }
//             ierr = MatSetValues(matAss, blockSize, globalx, blockSize, globaly, valpt, INSERT_VALUES);
//             CHKERRQ(ierr);
//             valpt += b2;
//         }
//     }

//     ierr = MatAssemblyBegin(matAss, MAT_FINAL_ASSEMBLY);
//     CHKERRQ(ierr);
//     ierr = MatAssemblyEnd(matAss, MAT_FINAL_ASSEMBLY);
//     CHKERRQ(ierr);

//     localAss = matAss;

//     finish = MPI_Wtime();
//     MPI_Barrier(MPI_COMM_WORLD);
//     if (myid == 0)
//         fprintf(fp, "get_SS_petsc time: %lf Sec\n", finish - start);
//     fclose(fp);

//     free(Ass.IA);
//     free(Ass.JA);
//     free(Ass.val);

//     return 0;
// }
