int FIM_solver(const int commRoot, int myid, int num_procs, int nrow, int nb, int *rpt, int *cpt, double *val, double *rhs, double *sol)
{
    // Partitioning A and b
    
    int local_size;
    
    MPI_Bcast(&nrow, 1, MPI_INT, commRoot, MPI_COMM_WORLD);
    MPI_Bcast(&nb, 1, MPI_INT, commRoot, MPI_COMM_WORLD);
    
    int matrixDim = nrow*nb;
    int nb2=nb*nb;
    int rhs_size=nrow*nb;
    
    if (myid!=commRoot){
        rhs=(double *)malloc(sizeof(double)*rhs_size);
        rpt=(int *)malloc(sizeof(double)*(nrow+1));
    }
    
    MPI_Bcast(rhs, rhs_size, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);
    MPI_Bcast(rpt, nrow+1, MPI_INT, commRoot, MPI_COMM_WORLD);
    
    int *allLower = (int *)malloc(sizeof(int)*num_procs);
    int *allUpper = (int *)malloc(sizeof(int)*num_procs);
    int *allDisp = (int *)malloc(sizeof(int)*num_procs);
    int *allSize = (int *)malloc(sizeof(int)*num_procs);
    
    calBLowerUpper(myid, num_procs, nrow, nb, rpt, local_size, allLower, allUpper, allDisp, allSize);
    
    int iStart = rpt[allLower[myid]];
    int *local_rpt = (int *)malloc(sizeof(int)*(local_size+1));
    double *local_rhs = (double *)malloc(sizeof(double)*local_size*nb);
    double *local_sol = (double *)malloc(sizeof(double)*local_size*nb);
    
    int i, j;
    for (i=0; i<local_size+1; ++i)
    {
        local_rpt[i] = rpt[allLower[myid]+i]-iStart; // global index to local index
    }
    for (i=0; i<local_size; ++i)
    {
        for (j=0; j<nb; ++j)
        {
            local_rhs[i*nb+j]=rhs[(allLower[myid]+i)*nb+j];
            local_sol[i*nb+j]=0.0;
        }
    }
    
    int local_nnz = local_rpt[local_size]-local_rpt[0];
    int *local_cpt = (int *)malloc(sizeof(int)*local_nnz);
    double *local_val = (double *)malloc(sizeof(double)*local_nnz*nb2);
    
    MPI_Scatterv(cpt, allSize, allDisp, MPI_INT, local_cpt, local_nnz, MPI_INT, commRoot, MPI_COMM_WORLD);
    for (i=0; i<num_procs; ++i)
    {
        allSize[i]*=nb2;
        allDisp[i]*=nb2;
    }
    MPI_Scatterv(val, allSize, allDisp, MPI_DOUBLE, local_val, local_nnz*nb2, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);
    decoup(local_val, local_rhs, local_rpt, local_cpt, nb, local_size, local_nnz);
    //----------------------------------------------------------------------
    // Initializing MAT and VEC
    
    Vec            x,b;  /* approx solution, RHS, exact solution */
    Mat            A;        /* linear system matrix */
    KSP            ksp;     /* linear solver context */
    PC             pc;
    //PetscReal      norm;     /* norm of solution error */
    PetscInt       Ii,its;
    PetscErrorCode ierr;
    
    int Istart = allLower[myid];
    int Iend = allUpper[myid];
    int blockSize = nb;
    
    int rowWidth = (Iend-Istart+1)*blockSize;
    int nBlockRows = Iend-Istart+1;
    int* nDCount = (int*)malloc(sizeof(int)*nBlockRows);
    int* nNDCount = (int*)malloc(sizeof(int)*nBlockRows);
    
    for (i=0; i<nBlockRows; i++) {
        nDCount[i] = 0;
        for (j=local_rpt[i]; j<local_rpt[i+1]; j++) {
            if(local_cpt[j]>=Istart && local_cpt[j]<=Iend)
            {
                nDCount[i]++;
            }
        }
    }
    
    for (i=0; i<nBlockRows; i++) {
        nNDCount[i] = local_rpt[i+1]-local_rpt[i]-nDCount[i];
    }
    
    MatCreateBAIJ(PETSC_COMM_WORLD,blockSize,rowWidth,rowWidth,matrixDim,matrixDim,0,nDCount,0, nNDCount, &A);
    
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
    
    double* valpt = local_val;
    int* globalx = (int*)malloc(blockSize*sizeof(int));
    int* globaly = (int*)malloc(blockSize*sizeof(int));
    for (Ii=0; Ii<nBlockRows; Ii++) {
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
    context.lower=allLower;
    context.upper=allUpper;
    
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
    
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
    if(myid==0)
    {
        ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
        PetscReal rnorm;
        ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
        ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", its);
    }
    
    int* xIndex = (int*)malloc(rowWidth * sizeof(int));
    for (i=0; i<rowWidth; i++) {
        xIndex[i]=i+Istart*blockSize;  // !!!!
    }
    VecGetValues(x,rowWidth,xIndex,local_sol);
    
    /////////////////////////////////////////////////////////
    
    for (i=0; i<num_procs; ++i)
    {
        allSize[i]=(allUpper[i]-allLower[i]+1)*nb;
        allDisp[i]=allLower[i]*nb;
    }
    MPI_Gatherv (local_sol, local_size*nb, MPI_DOUBLE, sol, allSize, allDisp, MPI_DOUBLE, commRoot, MPI_COMM_WORLD);
    
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
    
    if(myid!=commRoot)
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

