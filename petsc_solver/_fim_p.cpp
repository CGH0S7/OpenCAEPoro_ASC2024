int FIM_solver_p(int nrow, int nnz, int nb, int *rpt, int *cpt, double *val, double *rhs, double *sol)
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
    
    MatCreateBAIJ(PETSC_COMM_WORLD,blockSize,rowWidth,rowWidth,matrixDim,matrixDim,0,nDCount,0, nNDCount, &A);
    
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