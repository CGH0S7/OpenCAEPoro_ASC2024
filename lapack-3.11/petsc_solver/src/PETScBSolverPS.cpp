#include <stdlib.h>
#include <time.h>

#include "PETScBSolverPS.h"
// #ifdef __cplusplus
// extern "C" {
// #endif

// #include "petscksp.h"
// #include "PETScSolver.h"

void PETScBSolverPS::allocate(int row, int col, int nb, int nnz, int storage_manner)
{

    if (row > 0)
    {
        Absr.IA = (int *)calloc(row + 1, sizeof(int));
        b.val = (double *)calloc(row * nb, sizeof(double));
    }
    else
    {
        Absr.IA = NULL;
        b.val = NULL;
    }

    if (nnz > 0)
    {
        Absr.JA = (int *)calloc(nnz, sizeof(int));
        x.val = (double *)calloc(col * nb, sizeof(double));
    }
    else
    {
        Absr.JA = NULL;
        x.val = NULL;
    }

    if (nb > 0 && nnz > 0)
    {
        Absr.val = (double *)calloc(nnz * nb * nb, sizeof(double));
    }
    else
    {
        Absr.val = NULL;
    }

    Absr.ROW = row;
    Absr.COL = col;
    Absr.NNZ = nnz;
    Absr.nb = nb;
    Absr.storage_manner = storage_manner;
    b.row = row * nb;
    x.row = col * nb;
}

void PETScBSolverPS::destroy()
{
    free(Absr.IA);
    free(Absr.JA);
    free(Absr.val);
    free(b.val);
    free(x.val);
    Absr.ROW = 0;
    Absr.COL = 0;
    Absr.NNZ = 0;
    Absr.nb = 0;
    Absr.storage_manner = 0;
    b.row = 0;
    x.row = 0;
}

int PETScBSolverPS::petscsolve()
{
    int myid, num_procs;
    int parallel_done = 0;
    const int commRoot = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Bcast(&parallel_done, 1, MPI_INT, commRoot, MPI_COMM_WORLD);
    return FIM_solver(commRoot, myid, num_procs, Absr.ROW, Absr.nb, Absr.IA, Absr.JA, Absr.val, b.val, x.val);
}

// #ifdef __cplusplus
// extern "C" {
// #endif

int FIM_solver(const int commRoot, int myid, int num_procs, int nrow, int nb, int *rpt, int *cpt, double *val, double *rhs, double *sol)
{
    // Partitioning A and b

    int local_size;

    MPI_Bcast(&nrow, 1, MPI_INT, commRoot, MPI_COMM_WORLD);
    MPI_Bcast(&nb, 1, MPI_INT, commRoot, MPI_COMM_WORLD);

    int matrixDim = nrow * nb;
    int nb2 = nb * nb;
    int rhs_size = nrow * nb;

    if (myid != commRoot)
    {
        rhs = (double *)malloc(sizeof(double) * rhs_size);
        rpt = (int *)malloc(sizeof(double) * (nrow + 1));
    }

    MPI_Bcast(rhs, rhs_size, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);
    MPI_Bcast(rpt, nrow + 1, MPI_INT, commRoot, MPI_COMM_WORLD);

    int *allLower = (int *)malloc(sizeof(int) * num_procs);
    int *allUpper = (int *)malloc(sizeof(int) * num_procs);
    int *allDisp = (int *)malloc(sizeof(int) * num_procs);
    int *allSize = (int *)malloc(sizeof(int) * num_procs);

    calBLowerUpper(myid, num_procs, nrow, nb, rpt, local_size, allLower, allUpper, allDisp, allSize);

    int iStart = rpt[allLower[myid]];
    int *local_rpt = (int *)malloc(sizeof(int) * (local_size + 1));
    double *local_rhs = (double *)malloc(sizeof(double) * local_size * nb);
    double *local_sol = (double *)malloc(sizeof(double) * local_size * nb);

    int i, j;
    for (i = 0; i < local_size + 1; ++i)
    {
        local_rpt[i] = rpt[allLower[myid] + i] - iStart;
    }
    for (i = 0; i < local_size; ++i)
    {
        for (j = 0; j < nb; ++j)
        {
            local_rhs[i * nb + j] = rhs[(allLower[myid] + i) * nb + j];
            local_sol[i * nb + j] = 0.0;
        }
    }

    int local_nnz = local_rpt[local_size] - local_rpt[0];
    int *local_cpt = (int *)malloc(sizeof(int) * local_nnz);
    double *local_val = (double *)malloc(sizeof(double) * local_nnz * nb2);

    MPI_Scatterv(cpt, allSize, allDisp, MPI_INT, local_cpt, local_nnz, MPI_INT, commRoot, MPI_COMM_WORLD);
    for (i = 0; i < num_procs; ++i)
    {
        allSize[i] *= nb2;
        allDisp[i] *= nb2;
    }
    MPI_Scatterv(val, allSize, allDisp, MPI_DOUBLE, local_val, local_nnz * nb2, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);

    decoup_abf_4x4(local_val, local_rhs, local_rpt, local_cpt, nb, local_size);
    //----------------------------------------------------------------------
    // Initializing MAT and VEC

    Vec x, b; /* approx solution, RHS, exact solution */
    Mat A;    /* linear system matrix */
    KSP ksp;  /* linear solver context */
    PC pc;
    // PetscReal      norm;     /* norm of solution error */
    PetscInt Ii, its;
    PetscErrorCode ierr;

    int Istart = allLower[myid];
    int Iend = allUpper[myid];
    int blockSize = nb;

    int rowWidth = (Iend - Istart + 1) * blockSize;
    int nBlockRows = Iend - Istart + 1;
    int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
    int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

    for (i = 0; i < nBlockRows; i++)
    {
        nDCount[i] = 0;
        for (j = local_rpt[i]; j < local_rpt[i + 1]; j++)
        {
            if (local_cpt[j] >= Istart && local_cpt[j] <= Iend)
            {
                nDCount[i]++;
            }
        }
    }

    for (i = 0; i < nBlockRows; i++)
    {
        nNDCount[i] = local_rpt[i + 1] - local_rpt[i] - nDCount[i];
    }

    MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, matrixDim, matrixDim, 0, nDCount, 0, nNDCount, &A);

    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    ierr = MatSetUp(A);
    CHKERRQ(ierr);

    double *valpt = local_val;
    int *globalx = (int *)malloc(blockSize * sizeof(int));
    int *globaly = (int *)malloc(blockSize * sizeof(int));
    for (Ii = 0; Ii < nBlockRows; Ii++)
    {
        for (i = 0; i < blockSize; i++)
        {
            globalx[i] = (Ii + Istart) * blockSize + i;
        }

        for (i = local_rpt[Ii]; i < local_rpt[Ii + 1]; i++)
        {

            for (j = 0; j < blockSize; j++)
            {
                globaly[j] = local_cpt[i] * blockSize + j;
            }
            ierr = MatSetValues(A, blockSize, globalx, blockSize, globaly, valpt, INSERT_VALUES);
            CHKERRQ(ierr);
            valpt += blockSize * blockSize;
        }
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, rowWidth, matrixDim);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);

    int *bIndex = (int *)malloc(rowWidth * sizeof(int));
    for (i = 0; i < rowWidth; i++)
    {
        bIndex[i] = Istart * blockSize + i;
    }

    VecSetValues(b, rowWidth, bIndex, local_rhs, INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    dBSRmat_ BMat;
    BMat.ROW = local_size;
    BMat.COL = nrow;
    BMat.NNZ = local_nnz;
    BMat.nb = nb;
    BMat.IA = local_rpt;
    BMat.JA = local_cpt;
    BMat.val = local_val;

    Mat App;
    Mat Ass;
    get_PP(&BMat, Istart, Iend, matrixDim, App);
    get_SS(&BMat, Istart, Iend, matrixDim, Ass);

    // Set preconditoner
    shellContext context;
    context.BMat = A;
    context.App = App;
    context.Ass = Ass;
    context.nb = blockSize;
    context.lower = allLower;
    context.upper = allUpper;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);

    ierr = KSPSetType(ksp, KSPGMRES);
    CHKERRQ(ierr);

    ierr = KSPSetTolerances(ksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRQ(ierr);

    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);

    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCSHELL);
    PCShellSetApply(pc, precondApplyMSP);

    PCShellSetContext(pc, (void *)&context);

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);

    /*
     Check the error
     */
    if (myid == 0)
    {
        ierr = KSPGetIterationNumber(ksp, &its);
        CHKERRQ(ierr);
        PetscReal rnorm;
        ierr = KSPGetResidualNorm(ksp, &rnorm);
        CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", its);
    }

    int *xIndex = (int *)malloc(rowWidth * sizeof(int));
    for (i = 0; i < rowWidth; i++)
    {
        xIndex[i] = i + Istart * blockSize; // !!!!
    }
    VecGetValues(x, rowWidth, xIndex, local_sol);

    /////////////////////////////////////////////////////////

    for (i = 0; i < num_procs; ++i)
    {
        allSize[i] = (allUpper[i] - allLower[i] + 1) * nb;
        allDisp[i] = allLower[i] * nb;
    }
    MPI_Gatherv(local_sol, local_size * nb, MPI_DOUBLE, sol, allSize, allDisp, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);

    free(nDCount);
    free(nNDCount);
    free(xIndex);

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);

    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&App);
    CHKERRQ(ierr);
    ierr = MatDestroy(&Ass);
    CHKERRQ(ierr);

    free(globalx);
    free(globaly);

    if (myid != commRoot)
    {
        free(rhs);
        free(rpt);
    }
    free(allLower);
    free(allUpper);
    free(allDisp);
    free(allSize);

    free(local_cpt);
    free(local_rpt);
    free(local_val);
    free(local_sol);
    free(local_rhs);

    return 0;
}

#if 0
int FIM_solver_p(const int commRoot, int myid, int num_procs, int nrow, int ncol, int nnz, int nb, int *rpt, int *cpt, double *val, double *rhs, double *sol)
{

    int matrixDim = nrow*nb;
    int nb2=nb*nb;
    int rhs_size=nrow*nb;

    decoup_4x4(val, rhs, rpt, cpt, nb, nrow, nnz);
    //----------------------------------------------------------------------
    // Initializing MAT and VEC
    
    Vec            x,b;  /* approx solution, RHS, exact solution */
    Mat            A;        /* linear system matrix */
    KSP            ksp;     /* linear solver context */
    PC             pc;
    //PetscReal      norm;     /* norm of solution error */
    PetscInt       Ii,its;
    PetscErrorCode ierr;
    int i;
    int Istart = 0;
    int Iend = nrow-1;
    int blockSize = nb;

    int nBlockRows = Iend-Istart+1;
    int rowWidth = (Iend-Istart+1)*blockSize;
    int colWidth = ncol*blockSize;
    int* nDCount = (int*)malloc(sizeof(int)*nBlockRows);
    int* nNDCount = (int*)malloc(sizeof(int)*nBlockRows);
    
    for (int i=0; i<nBlockRows; i++) {
        nDCount[i] = 0;
        for (int j=rpt[i]; j<rpt[i+1]; j++) {
            if(cpt[j]>=Istart && cpt[j]<=Iend)
            {
                nDCount[i]++;
            }
        }
    }
    
    for (int i=0; i<nBlockRows; i++) {
        nNDCount[i] = rpt[i+1]-rpt[i]-nDCount[i];
    }
    
    MatCreateBAIJ(PETSC_COMM_WORLD,blockSize,rowWidth,colWidth,matrixDim,matrixDim,0,nDCount,0, nNDCount, &A);
    
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
    
    double* valpt = val;
    int* globalx = (int*)malloc(blockSize*sizeof(int));
    int* globaly = (int*)malloc(blockSize*sizeof(int));
    for (Ii=0; Ii<nBlockRows; Ii++) {
        for(i=0; i<blockSize; i++)
        {
            globalx[i] = (Ii+Istart)*blockSize + i;
        }
        
        for (int i=rpt[Ii]; i<rpt[Ii+1]; i++) {
            
            for (int j = 0; j<blockSize; j++) {
                globaly[j] = cpt[i]*blockSize + j;
            }
            ierr = MatSetValues(A,blockSize,globalx,blockSize,globaly,valpt,INSERT_VALUES);CHKERRQ(ierr);
            valpt += blockSize * blockSize;
        }
    }
    
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,colWidth,matrixDim);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
    
    int* bIndex = (int*)malloc(colWidth*sizeof(int));
    for (i=0; i<rowWidth; i++) {
        bIndex[i] = Istart*blockSize+i;
    }
    
    VecSetValues(b,colWidth,bIndex,rhs,INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    
    dBSRmat_ BMat;
    BMat.ROW = nrow;
    BMat.COL = nrow;
    BMat.NNZ = nnz;
    BMat.nb = nb;
    BMat.IA = rpt;
    BMat.JA = cpt;
    BMat.val = val;
    
    Mat App;
    Mat Ass;
    get_PP(&BMat, Istart, Iend, matrixDim, App);
    get_SS(&BMat, Istart, Iend, matrixDim, Ass);
    
    // Set preconditoner
    shellContext context;
    context.BMat = A;
    context.App = App;
    context.Ass = Ass;
    context.nb = blockSize;
    context.lower=allLower;
    context.upper=allUpper;
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    
    ierr = KSPSetType(ksp,KSPGMRES);CHKERRQ(ierr);
    
    ierr = KSPSetTolerances(ksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCSHELL);
    PCShellSetApply(pc,precondApply);
    
    PCShellSetContext(pc,(void*)&context);
    
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    
    /*
     Check the error
     */
    // if(myid==0)
    // {
    //     ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    //     PetscReal rnorm;
    //     ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
    //     ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
    //     ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", its);
    // }

    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    PetscReal rnorm;
    ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", its);

    
    int* xIndex = (int*)malloc(colWidth * sizeof(int));
    for (i=0; i<colWidth; i++) {
        xIndex[i]=i+Istart*blockSize;  // !!!!
    }
    VecGetValues(x,colWidth,xIndex,sol);
    
    /////////////////////////////////////////////////////////
    
    // for (i=0; i<num_procs; ++i)
    // {
    //     allSize[i]=(allUpper[i]-allLower[i]+1)*nb;
    //     allDisp[i]=allLower[i]*nb;
    // }
    // MPI_Gatherv (local_sol, local_size*nb, MPI_DOUBLE, sol, allSize, allDisp, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);
    
    free(nDCount);
    free(nNDCount);
    free(xIndex);
    
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&App);CHKERRQ(ierr);
    ierr = MatDestroy(&Ass);CHKERRQ(ierr);
    
    free(globalx);
    free(globaly);
    
    // if(myid!=commRoot)
    // {
    //     free(rhs);
    //     free(rpt);
    // }
    // free(allLower);
    // free(allUpper);
    // free(allDisp);
    // free(allSize);
    
    // free(local_cpt);
    // free(local_rpt);
    // free(local_val);
    // free(local_sol);
    // free(local_rhs);
    
    return 0;
}
#endif

#if 0
int fim_solver_(const int *commRoot, int *myid, int *num_procs, int *nrow, int *nb, int *rpt, int *cpt, double *val, double *rhs, double *sol)
{
    // Partitioning A and b
    //std::cout<<"fim_solver_ get in "<<std::endl;
    printf("fim_solver_ get in\n");

    int local_size;
    clock_t start, end;
    double Runtime;

        
    int matrixDim = *nrow**nb;
    int nb2=*nb**nb;
    int rhs_size=*nrow**nb;
    int rpt_size=*nrow + 1;

    printf("matrixDim = %d",matrixDim);
    printf("matrixDim = %d",nb2);
    
    if (*myid!=*commRoot){
        rhs=(double *)malloc(sizeof(double)*rhs_size);
        rpt=(int *)malloc(sizeof(int)*(*nrow+1));
    }
    
    start = clock();
    MPI_Bcast(&*nrow, 1, MPI_INT, *commRoot, MPI_COMM_WORLD);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MPI_Bcast1 done in %.16lf s\n", Runtime);

    start = clock();
    MPI_Bcast(&*nb, 1, MPI_INT, *commRoot, MPI_COMM_WORLD);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MPI_Bcast2 done in %.16lf s\n", Runtime);

    
    start = clock();
    MPI_Bcast(rhs, rhs_size, MPI_DOUBLE, *commRoot, MPI_COMM_WORLD);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MPI_Bcast3 done in %.16lf s\n", Runtime);
    
    start = clock();
    printf("MPI_Bcast4 beginnnnn");
    MPI_Bcast(rpt, rpt_size, MPI_INT, *commRoot, MPI_COMM_WORLD);
    printf("MPI_Bcast4 done");
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MPI_Bcast4 done in %.16lf s\n", Runtime);
    
    int *allLower = (int *)malloc(sizeof(int)**num_procs);
    int *allUpper = (int *)malloc(sizeof(int)**num_procs);
    int *allDisp = (int *)malloc(sizeof(int)**num_procs);
    int *allSize = (int *)malloc(sizeof(int)**num_procs);
    
    start = clock();
    calBLowerUpper(*myid, *num_procs, *nrow, *nb, rpt, local_size, allLower, allUpper, allDisp, allSize);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("calBLowerUpper1 done in %.16lf s\n", Runtime);
    
    int iStart = rpt[allLower[*myid]];
    int *local_rpt = (int *)malloc(sizeof(int)*(local_size+1));
    double *local_rhs = (double *)malloc(sizeof(double)*local_size**nb);
    double *local_sol = (double *)malloc(sizeof(double)*local_size**nb);
    
    int i, j;
    for (i=0; i<local_size+1; ++i)
    {
        local_rpt[i] = rpt[allLower[*myid]+i]-iStart;
    }
    for (i=0; i<local_size; ++i)
    {
        for (j=0; j<*nb; ++j)
        {
            local_rhs[i**nb+j]=rhs[(allLower[*myid]+i)**nb+j];
            local_sol[i**nb+j]=0.0;
        }
    }
    
    int local_nnz = local_rpt[local_size]-local_rpt[0];
    int *local_cpt = (int *)malloc(sizeof(int)*local_nnz);
    double *local_val = (double *)malloc(sizeof(double)*local_nnz*nb2);
    
    start = clock();
    MPI_Scatterv(cpt, allSize, allDisp, MPI_INT, local_cpt, local_nnz, MPI_INT, *commRoot, MPI_COMM_WORLD);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MPI_Scatterv1 done in %.16lf s\n", Runtime);
    
    for (i=0; i<*num_procs; ++i)
    {
        allSize[i]*=nb2;
        allDisp[i]*=nb2;
    }

    start = clock();
    MPI_Scatterv(val, allSize, allDisp, MPI_DOUBLE, local_val, local_nnz*nb2, MPI_DOUBLE, *commRoot, MPI_COMM_WORLD);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MPI_Scatterv2 done in %.16lf s\n", Runtime);
    decoup_2x2(local_val, local_rhs, local_rpt, local_cpt, *nb, local_size, local_nnz);
    printf("decoup done\n");
    //----------------------------------------------------------------------
    // Initializing MAT and VEC

    start = clock();    
    // if(*myid==0)
    //     PetscInitializeNoArguments();
    
    Vec            x,b;  /* approx solution, RHS, exact solution */
    Mat            A;        /* linear system matrix */
    KSP            ksp;     /* linear solver context */
    PC             pc;
    //PetscReal      norm;     /* norm of solution error */
    PetscInt       Ii,its;
    PetscErrorCode ierr;
    
    int Istart = allLower[*myid];
    int Iend = allUpper[*myid];
    int blockSize = *nb;
    
    int rowWidth = (Iend-Istart+1)*blockSize;
    int nblockRows = Iend-Istart+1;
    int* nDCount = (int*)malloc(sizeof(int)*nblockRows);
    int* nNDCount = (int*)malloc(sizeof(int)*nblockRows);

    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PetscInitialize done in %.16lf s\n", Runtime);
    

    start = clock();
    for (i=0; i<nblockRows; i++) {
        nDCount[i] = 0;
        for (j=local_rpt[i]; j<local_rpt[i+1]; j++) {
            if(local_cpt[j]>=Istart && local_cpt[j]<=Iend)
            {
                nDCount[i]++;
            }
        }
    }
    
    for (i=0; i<nblockRows; i++) {
        nNDCount[i] = local_rpt[i+1]-local_rpt[i]-nDCount[i];
    }
    
    MatCreateBAIJ(PETSC_COMM_WORLD,blockSize,rowWidth,rowWidth,matrixDim,matrixDim,0,nDCount,0, nNDCount, &A);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MatCreateBAIJ done in %.16lf s\n", Runtime);
    

    start = clock();
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
    
    double* valpt = local_val;
    int* globalx = (int*)malloc(blockSize*sizeof(int));
    int* globaly = (int*)malloc(blockSize*sizeof(int));
    for (Ii=0; Ii<nblockRows; Ii++) {
        for(i=0; i<blockSize; i++)
        {
            globalx[i] = (Ii+Istart)*blockSize + i;
        }
        
        for (i=local_rpt[Ii]; i<local_rpt[Ii+1]; i++) {
            
            for (j = 0; j<blockSize; j++) {
                globaly[j] = local_cpt[i]*blockSize + j;
            }
            ierr = MatSetValues(A,blockSize,globalx,blockSize,globaly,valpt,INSERT_VALUES);CHKERRQ(ierr);
            valpt += blockSize * blockSize;
        }
    }
    
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,rowWidth,matrixDim);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
    
    int* bIndex = (int*)malloc(rowWidth*sizeof(int));
    for (i=0; i<rowWidth; i++) {
        bIndex[i] = Istart*blockSize+i;
    }
    
    VecSetValues(b,rowWidth,bIndex,local_rhs,INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("Mat/VEC Assembly done in %.16lf s\n", Runtime);


    start = clock();
    dBSRmat_ BMat;
    BMat.ROW = local_size;
    BMat.COL = *nrow;
    BMat.NNZ = local_nnz;
    BMat.nb = *nb;
    BMat.IA = local_rpt;
    BMat.JA = local_cpt;
    BMat.val = local_val;
    
    Mat App;
    Mat Ass;
    get_PP(&BMat, Istart, Iend, matrixDim, App);
    get_SS(&BMat, Istart, Iend, matrixDim, Ass);
    
    // Set preconditoner
    shellContext context;
    context.BMat = A;
    context.App = App;
    context.Ass = Ass;
    context.nb = blockSize;
    context.lower=allLower;
    context.upper=allUpper;

    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("Set preconditoner done in %.16lf s\n", Runtime);
    
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    
    ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
    
    ierr = KSPSetTolerances(ksp, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    
    start = clock();
    KSPGetPC(ksp,&pc);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("KSPGetPC done in %.16lf s\n", Runtime);

    start = clock();
    PCSetType(pc,PCSHELL);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PCSetType done in %.16lf s\n", Runtime);

    start = clock();
    PCShellSetApply(pc,precondApply);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PCShellSetApply done in %.16lf s\n", Runtime);

    start = clock();
    PCShellSetContext(pc,(void*)&context);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PCShellSetContext done in %.16lf s\n", Runtime);

    ierr = KSPMonitorSet(ksp,KSPMonitorTrueResidualNorm,NULL,NULL);CHKERRQ(ierr);

    start = clock();
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("KSPSolve done in %.16lf s\n", Runtime);
    
    /*
     Check the error
     */
    if(*myid==0)
    {
        ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
        PetscReal rnorm;
        ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", its);
    }

    printf("Check the error done\n");
    
    int* xIndex = (int*)malloc(rowWidth * sizeof(int));
    for (i=0; i<rowWidth; i++) {
        xIndex[i]=i+Istart*blockSize;  // !!!!
    }
    VecGetValues(x,rowWidth,xIndex,local_sol);

    printf("VecGetValues done\n");
    
    /////////////////////////////////////////////////////////
    
    for (i=0; i<*num_procs; ++i)
    {
        allSize[i]=(allUpper[i]-allLower[i]+1)**nb;
        allDisp[i]=allLower[i]**nb;
    }
    MPI_Gatherv (local_sol, local_size**nb, MPI_DOUBLE, sol, allSize, allDisp, MPI_DOUBLE, *commRoot, MPI_COMM_WORLD);
    
    printf("MPI_Gatherv done\n");

    free(nDCount);
    free(nNDCount);
    free(xIndex);

    printf("free1 done\n");
    
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&App);CHKERRQ(ierr);
    ierr = MatDestroy(&Ass);CHKERRQ(ierr);
    
    free(globalx);
    free(globaly);

    printf("free2 done\n");
    
    if(*myid!=*commRoot)
    {
        free(rhs);
        free(rpt);
    }
    free(allLower);
    free(allUpper);
    free(allDisp);
    free(allSize);

    printf("free3 done\n");
    
    free(local_cpt);
    free(local_rpt);
    free(local_val);
    free(local_sol);
    free(local_rhs);

    printf("free4 done\n");
    
    // if(*myid==0)
    //     PetscFinalize();
    return 0;

}

// #ifdef __cplusplus
// }
#endif

#if 0
int fim_solver_(const int *commRoot, int *myid, int *num_procs, int *nrow, int *nb, int *rpt, int *cpt, double *val, double *rhs, double *sol)
{
    // Partitioning A and b
    //std::cout<<"fim_solver_ get in "<<std::endl;
    printf("fim_solver_ get in\n");

    int local_size;
    clock_t start, end;
    double Runtime;

        
    int matrixDim = *nrow**nb;
    int nb2=*nb**nb;
    int rhs_size=*nrow**nb;
    int rpt_size=*nrow + 1;

    printf("matrixDim = %d",matrixDim);
    printf("matrixDim = %d",nb2);
    
    
    // rhs=(double *)malloc(sizeof(double)*rhs_size);
    // rpt=(int *)malloc(sizeof(int)*(*nrow+1));

    
    // start = clock();
    // MPI_Bcast(&*nrow, 1, MPI_INT, *commRoot, MPI_COMM_WORLD);
    // end = clock();
    // Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    // printf("MPI_Bcast1 done in %.16lf s\n", Runtime);

    // start = clock();
    // MPI_Bcast(&*nb, 1, MPI_INT, *commRoot, MPI_COMM_WORLD);
    // end = clock();
    // Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    // printf("MPI_Bcast2 done in %.16lf s\n", Runtime);

    
    // start = clock();
    // MPI_Bcast(rhs, rhs_size, MPI_DOUBLE, *commRoot, MPI_COMM_WORLD);
    // end = clock();
    // Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    // printf("MPI_Bcast3 done in %.16lf s\n", Runtime);
    
    // start = clock();
    // printf("MPI_Bcast4 beginnnnn");
    // MPI_Bcast(rpt, rpt_size, MPI_INT, *commRoot, MPI_COMM_WORLD);
    // printf("MPI_Bcast4 done");
    // end = clock();
    // Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    // printf("MPI_Bcast4 done in %.16lf s\n", Runtime);
    
    int *allLower = (int *)malloc(sizeof(int)**num_procs);
    int *allUpper = (int *)malloc(sizeof(int)**num_procs);
    int *allDisp = (int *)malloc(sizeof(int)**num_procs);
    int *allSize = (int *)malloc(sizeof(int)**num_procs);
    
    start = clock();
    calBLowerUpper(*myid, *num_procs, *nrow, *nb, rpt, local_size, allLower, allUpper, allDisp, allSize);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("calBLowerUpper1 done in %.16lf s\n", Runtime);

    printf("local_size= %d \n", local_size);
    
    int iStart = rpt[allLower[*myid]];
    int *local_rpt = (int *)malloc(sizeof(int)*(local_size+1));
    double *local_rhs = (double *)malloc(sizeof(double)*local_size**nb);
    double *local_sol = (double *)malloc(sizeof(double)*local_size**nb);
    
    int i, j;
    for (i=0; i<local_size+1; ++i)
    {
        local_rpt[i] = rpt[i];
    }
    for (i=0; i<local_size; ++i)
    {
        for (j=0; j<*nb; ++j)
        {
            local_rhs[i**nb+j]=rhs[(i)**nb+j];
            local_sol[i**nb+j]=0.0;
        }
    }
    
    int local_nnz = local_rpt[local_size]-local_rpt[0];
    int *local_cpt = (int *)malloc(sizeof(int)*local_nnz);
    double *local_val = (double *)malloc(sizeof(double)*local_nnz*nb2);

    for (i=0; i<local_nnz; ++i)
    {
        local_cpt[i] = cpt[i];
    }
    for (i=0; i<local_nnz; ++i)
    {
        for (j=0; j<nb2; ++j)
        {
            local_val[i*nb2+j]=val[(i)*nb2+j];
        }
    }
    
    // start = clock();
    // MPI_Scatterv(cpt, allSize, allDisp, MPI_INT, local_cpt, local_nnz, MPI_INT, *commRoot, MPI_COMM_WORLD);
    // end = clock();
    // Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    // printf("MPI_Scatterv1 done in %.16lf s\n", Runtime);
    
    for (i=0; i<*num_procs; ++i)
    {
        allSize[i]*=nb2;
        allDisp[i]*=nb2;
    }

    // start = clock();
    // MPI_Scatterv(val, allSize, allDisp, MPI_DOUBLE, local_val, local_nnz*nb2, MPI_DOUBLE, *commRoot, MPI_COMM_WORLD);
    // end = clock();
    // Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    // printf("MPI_Scatterv2 done in %.16lf s\n", Runtime);

    for (i=0; i<10; ++i)
    {
        printf("rpt[i]=%d \n", rpt[i]);
    }
    for (i=0; i<10; ++i)
    {
        printf("local_rpt[i]=%d \n", local_rpt[i]);
    }

    for (i=0; i<10; ++i)
    {
        printf("val[i]=%.16lf \n", val[i]);
    }
    for (i=0; i<10; ++i)
    {
        printf("local_val[i]=%.16lf \n", local_val[i]);
    }

    decoup_2x2(local_val, local_rhs, local_rpt, local_cpt, *nb, local_size, local_nnz);
    printf("decoup done\n");
    //----------------------------------------------------------------------
    // Initializing MAT and VEC

    start = clock();    
    // if(*myid==0)
    //     PetscInitializeNoArguments();
    
    Vec            x,b;  /* approx solution, RHS, exact solution */
    Mat            A;        /* linear system matrix */
    KSP            ksp;     /* linear solver context */
    PC             pc;
    //PetscReal      norm;     /* norm of solution error */
    PetscInt       Ii,its;
    PetscErrorCode ierr;
    
    int Istart = allLower[*myid];
    int Iend = allUpper[*myid];
    int blockSize = *nb;
    
    int rowWidth = (Iend-Istart+1)*blockSize;
    int nblockRows = Iend-Istart+1;
    int* nDCount = (int*)malloc(sizeof(int)*nblockRows);
    int* nNDCount = (int*)malloc(sizeof(int)*nblockRows);

    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PetscInitialize done in %.16lf s\n", Runtime);
    

    start = clock();
    for (i=0; i<nblockRows; i++) {
        nDCount[i] = 0;
        for (j=local_rpt[i]; j<local_rpt[i+1]; j++) {
            if(local_cpt[j]>=Istart && local_cpt[j]<=Iend)
            {
                nDCount[i]++;
            }
        }
    }
    
    for (i=0; i<nblockRows; i++) {
        nNDCount[i] = local_rpt[i+1]-local_rpt[i]-nDCount[i];
    }
    printf("MatCreateBAIJ begin");
    
    MatCreateBAIJ(PETSC_COMM_WORLD,blockSize,rowWidth,rowWidth,matrixDim,matrixDim,0,nDCount,0, nNDCount, &A);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("MatCreateBAIJ done in %.16lf s\n", Runtime);
    

    start = clock();
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
    
    double* valpt = local_val;
    int* globalx = (int*)malloc(blockSize*sizeof(int));
    int* globaly = (int*)malloc(blockSize*sizeof(int));
    for (Ii=0; Ii<nblockRows; Ii++) {
        for(i=0; i<blockSize; i++)
        {
            globalx[i] = (Ii+Istart)*blockSize + i;
        }
        
        for (i=local_rpt[Ii]; i<local_rpt[Ii+1]; i++) {
            
            for (j = 0; j<blockSize; j++) {
                globaly[j] = local_cpt[i]*blockSize + j;
            }
            ierr = MatSetValues(A,blockSize,globalx,blockSize,globaly,valpt,INSERT_VALUES);CHKERRQ(ierr);
            valpt += blockSize * blockSize;
        }
    }
    
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
    ierr = VecSetSizes(b,rowWidth,matrixDim);CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);CHKERRQ(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
    
    int* bIndex = (int*)malloc(rowWidth*sizeof(int));
    for (i=0; i<rowWidth; i++) {
        bIndex[i] = Istart*blockSize+i;
    }
    
    VecSetValues(b,rowWidth,bIndex,local_rhs,INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("Mat/VEC Assembly done in %.16lf s\n", Runtime);


    start = clock();
    dBSRmat_ BMat;
    BMat.ROW = local_size;
    BMat.COL = *nrow;
    BMat.NNZ = local_nnz;
    BMat.nb = *nb;
    BMat.IA = local_rpt;
    BMat.JA = local_cpt;
    BMat.val = local_val;
    
    Mat App;
    Mat Ass;
    get_PP(&BMat, Istart, Iend, matrixDim, App);
    get_SS(&BMat, Istart, Iend, matrixDim, Ass);
    
    // Set preconditoner
    shellContext context;
    context.BMat = A;
    context.App = App;
    context.Ass = Ass;
    context.nb = blockSize;
    context.lower=allLower;
    context.upper=allUpper;

    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("Set preconditoner done in %.16lf s\n", Runtime);
    
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    
    ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);
    
    ierr = KSPSetTolerances(ksp, 1.e-3, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    
    start = clock();
    KSPGetPC(ksp,&pc);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("KSPGetPC done in %.16lf s\n", Runtime);

    start = clock();
    PCSetType(pc,PCSHELL);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PCSetType done in %.16lf s\n", Runtime);

    start = clock();
    PCShellSetApply(pc,precondApply);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PCShellSetApply done in %.16lf s\n", Runtime);

    start = clock();
    PCShellSetContext(pc,(void*)&context);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("PCShellSetContext done in %.16lf s\n", Runtime);

    ierr = KSPMonitorSet(ksp,KSPMonitorTrueResidualNorm,NULL,NULL);CHKERRQ(ierr);

    start = clock();
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    end = clock();
    Runtime = (double)(end - start)/CLOCKS_PER_SEC;
    printf("KSPSolve done in %.16lf s\n", Runtime);
    
    /*
     Check the error
     */
    if(*myid==0)
    {
        ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
        PetscReal rnorm;
        ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", its);
    }

    printf("Check the error done\n");
    
    int* xIndex = (int*)malloc(rowWidth * sizeof(int));
    for (i=0; i<rowWidth; i++) {
        xIndex[i]=i+Istart*blockSize;  // !!!!
    }
    VecGetValues(x,rowWidth,xIndex,local_sol);

    printf("VecGetValues done\n");
    
    /////////////////////////////////////////////////////////
    
    for (i=0; i<*num_procs; ++i)
    {
        allSize[i]=(allUpper[i]-allLower[i]+1)**nb;
        allDisp[i]=allLower[i]**nb;
    }
    MPI_Gatherv (local_sol, local_size**nb, MPI_DOUBLE, sol, allSize, allDisp, MPI_DOUBLE, *commRoot, MPI_COMM_WORLD);
    
    printf("MPI_Gatherv done\n");

    free(nDCount);
    free(nNDCount);
    free(xIndex);

    printf("free1 done\n");
    
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = MatDestroy(&App);CHKERRQ(ierr);
    ierr = MatDestroy(&Ass);CHKERRQ(ierr);
    
    free(globalx);
    free(globaly);

    printf("free2 done\n");
    
    if(*myid!=*commRoot)
    {
        free(rhs);
        free(rpt);
    }
    free(allLower);
    free(allUpper);
    free(allDisp);
    free(allSize);

    printf("free3 done\n");
    
    free(local_cpt);
    free(local_rpt);
    free(local_val);
    free(local_sol);
    free(local_rhs);

    printf("free4 done\n");
    
    // if(*myid==0)
    //     PetscFinalize();
    return 0;

}

// #ifdef __cplusplus
// }
#endif
