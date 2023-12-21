#include <stdlib.h>
#include <time.h>

#include "PETScBSolverPS.h"

#define TOL 1e-3     // tolerance for linear solver
#define MAXIT 200    // Maximum iteration number of linear solver
#define DecoupType 1 // 0: None, 1: ABF, 2: ANL, 3: SEM, 4: QI

int Call_Solver_times = 0;
int Is_Convergence = 1;
int NonConNum = 0; // Count the number of non convergence times
double Matrix_Conversion_Time = 0.0; // Global variables: matrix conversion time

int FIM_solver_p(int precondID, int is_thermal, int myid, int num_procs, int nb,
                 int *allLower, int *allUpper, int *rpt, int *cpt, double *val,
                 double *rhs, double *sol, double *mct) {
  int status = 0;
  precondID = -1; // 0: MSP, 1: BAMG, -1: CPR
  is_thermal = 0; // 0: non-thermal, 1: thermal

  switch (precondID) {
  case 0:
    status = FIM_solver_p_msp(is_thermal, myid, num_procs, nb, allLower,
                              allUpper, rpt, cpt, val, rhs, sol);
    break;

  case 1:
    status = FIM_solver_p_bamg(is_thermal, myid, num_procs, nb, allLower,
                               allUpper, rpt, cpt, val, rhs, sol);
    break;

  default: // case -1
    status = FIM_solver_p_cpr(is_thermal, myid, num_procs, nb, allLower,
                              allUpper, rpt, cpt, val, rhs, sol);
    break;
  }

  *mct = Matrix_Conversion_Time;

  return status;
}

// CPR preconditioner
int FIM_solver_p_cpr(int is_thermal, int myid, int num_procs, int nb,
                     int *allLower, int *allUpper, int *rpt, int *cpt,
                     double *val, double *rhs, double *sol) {
  // printf("%s %s %d begin !\n", __FILE__, __FUNCTION__, __LINE__);
  // printf("PETSc_FIM_Solver begin !\n");
  FILE *fp;
  if (myid == 0 && PrintLvl > 0) {
    fp = fopen("linear_solver.log", "a");
    fprintf(fp, "Call_Solver_times: %d\n", Call_Solver_times);
  }

  //! Construct command line parameters for petsc
  int argc = 21, len = 128, i;
  char **args = (char **)malloc(argc * sizeof(char *));
  for (i = 0; i < argc; i++)
    args[i] = (char *)calloc(len, sizeof(char));

  strcpy(args[0], "program");
  strcpy(args[1],
         "-pc_hypre_boomeramg_coarsen_type"); // {"CLJP", "Ruge-Stueben", "",
                                              // "modifiedRuge-Stueben", "", "",
                                              // "Falgout", "", "PMIS", "",
                                              // "HMIS"};
  strcpy(args[2], "HMIS");
  strcpy(args[3], "-pc_hypre_boomeramg_interp_type");
  strcpy(args[4], "ext+i"); // default: classical, ext+i,  block block-wtd:
                            // no-ok/no-convergence
  // -pc_hypre_boomeramg_interp_type. Available options:
  // classical   direct multipass multipass-wts ext+i ext+i-cc standard
  // standard-wts block block-wtd ff ff1 ext ad-wts ext-mm ext+i-mm ext+e-mm
  strcpy(args[5], "-pc_hypre_boomeramg_relax_type_all"); // Relax type for the
                                                         // up and down cycles
  strcpy(args[6], "symmetric-SOR/Jacobi");               // Jacobi,
  // {"Jacobi", "sequential-Gauss-Seidel", "seqboundary-Gauss-Seidel",
  // "SOR/Jacobi", "backward-SOR/Jacobi", "" /* [5] hybrid chaotic Gauss-Seidel
  // (works only with OpenMP) */, "symmetric-SOR/Jacobi", "" /* 7 */,
  // "l1scaled-SOR/Jacobi", "Gaussian-elimination", "" /* 10 */, "" /* 11 */, ""
  // /* 12 */, "l1-Gauss-Seidel" /* nonsymmetric */, "backward-l1-Gauss-Seidel"
  // /* nonsymmetric */, "CG" /* non-stationary */, "Chebyshev", "FCF-Jacobi",
  // "l1scaled-Jacobi"};
  strcpy(args[7],
         "-pc_hypre_boomeramg_relax_type_coarse"); // Relax type on coarse grid
  strcpy(args[8], "Gaussian-elimination");         // Gaussian-elimination, CG
  strcpy(args[9],
         "-pc_hypre_boomeramg_P_max"); // Max elements per row for interpolation
                                       // operator (0=unlimited)
  strcpy(args[10], "4");               // 4
  strcpy(args[11],
         "-pc_hypre_boomeramg_truncfactor"); // Truncation factor for
                                             // interpolation (0=no truncation)
  strcpy(args[12], "0.01");                  // 0, 0.01, 0.05, 0.1
  strcpy(args[13],
         "-pc_hypre_boomeramg_min_coarse_size"); // Minimum size of coarsest
                                                 // grid, default: 1
  strcpy(args[14], "1");
  strcpy(args[15],
         "-pc_hypre_boomeramg_max_coarse_size"); // Maximum size of coarsest
                                                 // grid, default: 9
  strcpy(args[16], "100");
  strcpy(args[17],
         "-pc_hypre_boomeramg_grid_sweeps_coarse"); // Number of sweeps for the
                                                    // coarse level, default: 1
  strcpy(args[18], "1");
  strcpy(args[19], "-pc_hypre_boomeramg_print_statistics"); // print level
  strcpy(args[20], "0");
  // sprintf(args[5], "-ksp_monitor");
  //! end

  // diag-check
  double start, finish;
  // PetscInitializeNoArguments();
  PetscCall(PetscInitialize(&argc, &args, NULL, NULL)); // 有命令行参数
  int nrow_global = allUpper[num_procs - 1] + 1;
  int matrixDim = nrow_global * nb;
  int Istart = allLower[myid];
  int Iend = allUpper[myid];
  int nBlockRows = Iend - Istart + 1;
  int rowWidth = nBlockRows * nb;
  int nb2 = nb * nb;
  int rhs_size = nBlockRows * nb;
  int nnz = rpt[nBlockRows] - rpt[0];

  int blockSize = nb;
  int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
  int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

  MPI_Barrier(MPI_COMM_WORLD);

  ////////////////////////////////////////////////////////////////////
  // Modifiable Area  |
  // Modifiable Area  v
  /////////////////////////////////////////////////////////////////////
  start = MPI_Wtime();
  decoup(val, rhs, rpt, cpt, nb, nBlockRows, DecoupType, is_thermal);

  //----------------------------------------------------------------------
  // Initializing MAT and VEC

  Vec x, b; /* approx solution, RHS, exact solution */
  Mat A;    /* linear system matrix */
  KSP ksp;  /* linear solver context */
  PC pc;
  // PetscReal      norm;     /* norm of solution error */
  PetscInt Ii, iters;
  PetscErrorCode ierr;

#pragma acc parallel loop gang vector
  for (int i = 0; i < nBlockRows; i++) {
    nDCount[i] = 0;
    for (int j = rpt[i]; j < rpt[i + 1]; j++) {
      if (cpt[j] >= Istart && cpt[j] <= Iend) {
        nDCount[i]++;
      }
    }
  }

#pragma acc parallel loop
  for (int i = 0; i < nBlockRows; i++) {
    nNDCount[i] = rpt[i + 1] - rpt[i] - nDCount[i];
  }

  MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, matrixDim,
                matrixDim, 0, nDCount, 0, nNDCount, &A);

  ierr = MatSetFromOptions(A);
  CHKERRQ(ierr);
  ierr = MatSetUp(A);
  CHKERRQ(ierr);

  double *valpt = val;
  int *globalx = (int *)malloc(blockSize * sizeof(int));
  int *globaly = (int *)malloc(blockSize * sizeof(int));

  // #pragma acc parallel
  for (Ii = 0; Ii < nBlockRows; Ii++) {
#pragma acc parallel loop gang vector
    for (i = 0; i < blockSize; i++) {
      globalx[i] = (Ii + Istart) * blockSize + i;
    }
    for (int i = rpt[Ii]; i < rpt[Ii + 1]; i++) {

#pragma acc parallel loop vector
      for (int j = 0; j < blockSize; j++) {
        globaly[j] = cpt[i] * blockSize + j;
      }
      ierr = MatSetValues(A, blockSize, globalx, blockSize, globaly, valpt,
                          INSERT_VALUES);
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
#pragma acc parallel loop
  for (i = 0; i < rowWidth; i++) {
    bIndex[i] = Istart * blockSize + i;
  }

  VecSetValues(b, rowWidth, bIndex, rhs, INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  dBSRmat_ BMat;
  BMat.ROW = nBlockRows;
  BMat.COL = nrow_global;
  BMat.NNZ = nnz;
  BMat.nb = nb;
  BMat.IA = rpt;
  BMat.JA = cpt;
  BMat.val = val;

  Mat App;
  // Mat Ass;

  finish = MPI_Wtime();
  Matrix_Conversion_Time +=
      finish - start; //! matrix conversion time (Global variables)

  ////////////////////////////////////////////////////////////////////
  // Modifiable Area  ^
  // Modifiable Area  |
  /////////////////////////////////////////////////////////////////////

  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "petsc A_tmp assemble time: %lf Sec\n", finish - start);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  get_PP(&BMat, Istart, Iend, matrixDim, App);
  // get_SS(&BMat, Istart, Iend, matrixDim, Ass);

  finish = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "get_PP time: %lf Sec\n", finish - start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  // Set preconditoner
  shellContext context;
  context.BMat = A;
  context.App = App;
  // context.Ass = Ass;
  context.nb = blockSize;
  context.lower = allLower;
  context.upper = allUpper;

  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"shell Set preconditoner time: %lf
  // Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);
  CHKERRQ(ierr);

  ierr = KSPSetType(ksp, KSPFGMRES);
  CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp, TOL, PETSC_DEFAULT, PETSC_DEFAULT, MAXIT);
  CHKERRQ(ierr);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  KSPGetPC(ksp, &pc);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"KSPGetPC time: %lf Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  PCSetType(pc, PCSHELL);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"PCSetType time: %lf Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  PCShellSetApply(pc, precondApplyCPR);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"PCShellSetApply time: %lf Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  PCShellSetContext(pc, (void *)&context);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"PCShellSetContext time: %lf Sec\n",finish-start);

  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"KSPSet time: %lf Sec\n",finish-start);

  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  ierr = KSPSolve(ksp, b, x);
  CHKERRQ(ierr);

  finish = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "KSPSolve time: %lf Sec\n", finish - start);

  /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  if (Call_Solver_times == 0 && PrintLvl > 0)
    KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // view ksp parameters

  /*
   Check the error
   */
  // if(myid==0)
  // {
  //     ierr = KSPGetIterationNumber(ksp,&iters);CHKERRQ(ierr);
  //     PetscReal rnorm;
  //     ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
  //     ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
  //     ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n",
  //     iters);
  // }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check the solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPGetIterationNumber(ksp, &iters);
  CHKERRQ(ierr);
  PetscReal rnorm, bnorm;
  PetscCall(VecNorm(b, NORM_2, &bnorm));
  ierr = KSPGetResidualNorm(ksp, &rnorm);
  CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD, "PETSc_FIM_Solver done! number
  // iterations = %d, residual norm = %f\n", iters, rnorm);
  if (myid == 0 && PrintLvl > 0) {
    printf("iter = %d, residual norm = %e, RRN = %e\n", iters, rnorm,
           rnorm / bnorm);
    // write files
    fprintf(fp,
            "Call_Solver_times = %d, iter = %d, residual norm = %e, RRN = %e\n",
            Call_Solver_times, iters, rnorm, rnorm / bnorm);
    fprintf(fp, "----------------------------------------------\n");
    fclose(fp);
  }

  // ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", iters);
  int iter_done = iters;

  int *xIndex = (int *)malloc(rowWidth * sizeof(int));
  for (i = 0; i < rowWidth; i++) {
    xIndex[i] = i + Istart * blockSize; // !!!!
  }
  VecGetValues(x, rowWidth, xIndex, sol);

  /////////////////////////////////////////////////////////

  // for (i=0; i<num_procs; ++i)
  // {
  //     allSize[i]=(allUpper[i]-allLower[i]+1)*nb;
  //     allDisp[i]=allLower[i]*nb;
  // }
  // MPI_Gatherv (local_sol, local_size*nb, MPI_DOUBLE, sol, allSize, allDisp,
  // MPI_DOUBLE, commRoot, MPI_COMM_WORLD);

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
  // ierr = MatDestroy(&Ass);
  // CHKERRQ(ierr);

  free(globalx);
  free(globaly);

  // free command line parameters
  for (i = 0; i < argc; i++)
    free(args[i]);
  free(args);

  PetscFinalize();
  Call_Solver_times++;
  return iter_done;
}

// MSP preconditioner
int FIM_solver_p_msp(int is_thermal, int myid, int num_procs, int nb,
                     int *allLower, int *allUpper, int *rpt, int *cpt,
                     double *val, double *rhs, double *sol) {
  // printf("%s %s %d begin !\n", __FILE__, __FUNCTION__, __LINE__);
  // printf("PETSc_FIM_Solver begin !\n");
  FILE *fp;
  if (myid == 0 && PrintLvl > 0) {
    fp = fopen("linear_solver.log", "a");
    fprintf(fp, "Call_Solver_times: %d\n", Call_Solver_times);
  }

  // diag-check
  double start, finish;
  PetscInitializeNoArguments();
  int nrow_global = allUpper[num_procs - 1] + 1;
  int matrixDim = nrow_global * nb;
  int Istart = allLower[myid];
  int Iend = allUpper[myid];
  int nBlockRows = Iend - Istart + 1;
  int rowWidth = nBlockRows * nb;
  int nb2 = nb * nb;
  int rhs_size = nBlockRows * nb;
  int nnz = rpt[nBlockRows] - rpt[0];

  int blockSize = nb;
  int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
  int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  decoup(val, rhs, rpt, cpt, nb, nBlockRows, DecoupType, is_thermal);

  finish = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "decoup time: %lf Sec\n", finish - start);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  //----------------------------------------------------------------------
  // Initializing MAT and VEC

  Vec x, b; /* approx solution, RHS, exact solution */
  Mat A;    /* linear system matrix */
  KSP ksp;  /* linear solver context */
  PC pc;
  // PetscReal      norm;     /* norm of solution error */
  PetscInt Ii, iters;
  PetscErrorCode ierr;
  int i;

  for (int i = 0; i < nBlockRows; i++) {
    nDCount[i] = 0;
    for (int j = rpt[i]; j < rpt[i + 1]; j++) {
      if (cpt[j] >= Istart && cpt[j] <= Iend) {
        nDCount[i]++;
      }
    }
  }

  for (int i = 0; i < nBlockRows; i++) {
    nNDCount[i] = rpt[i + 1] - rpt[i] - nDCount[i];
  }

  MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, matrixDim,
                matrixDim, 0, nDCount, 0, nNDCount, &A);

  ierr = MatSetFromOptions(A);
  CHKERRQ(ierr);
  ierr = MatSetUp(A);
  CHKERRQ(ierr);

  double *valpt = val;
  int *globalx = (int *)malloc(blockSize * sizeof(int));
  int *globaly = (int *)malloc(blockSize * sizeof(int));
  for (Ii = 0; Ii < nBlockRows; Ii++) {
    for (i = 0; i < blockSize; i++) {
      globalx[i] = (Ii + Istart) * blockSize + i;
    }

    for (int i = rpt[Ii]; i < rpt[Ii + 1]; i++) {

      for (int j = 0; j < blockSize; j++) {
        globaly[j] = cpt[i] * blockSize + j;
      }
      ierr = MatSetValues(A, blockSize, globalx, blockSize, globaly, valpt,
                          INSERT_VALUES);
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
  for (i = 0; i < rowWidth; i++) {
    bIndex[i] = Istart * blockSize + i;
  }

  VecSetValues(b, rowWidth, bIndex, rhs, INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  dBSRmat_ BMat;
  BMat.ROW = nBlockRows;
  BMat.COL = nrow_global;
  BMat.NNZ = nnz;
  BMat.nb = nb;
  BMat.IA = rpt;
  BMat.JA = cpt;
  BMat.val = val;

  Mat App;
  Mat Ass;

  finish = MPI_Wtime();
  Matrix_Conversion_Time +=
      finish - start; //! matrix conversion time (Global variables)
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "petsc A_tmp assemble time: %lf Sec\n", finish - start);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  get_PP(&BMat, Istart, Iend, matrixDim, App);
  get_SS(&BMat, Istart, Iend, matrixDim, Ass);

  finish = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "get_PP_SS time: %lf Sec\n", finish - start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  // Set preconditoner
  shellContext context;
  context.BMat = A;
  context.App = App;
  context.Ass = Ass;
  context.nb = blockSize;
  context.lower = allLower;
  context.upper = allUpper;

  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"shell Set preconditoner time: %lf
  // Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);
  CHKERRQ(ierr);

  ierr = KSPSetType(ksp, KSPFGMRES);
  CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp, TOL, PETSC_DEFAULT, PETSC_DEFAULT, MAXIT);
  CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  KSPGetPC(ksp, &pc);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"KSPGetPC time: %lf Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  PCSetType(pc, PCSHELL);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"PCSetType time: %lf Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  PCShellSetApply(pc, precondApplyMSP);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"PCShellSetApply time: %lf Sec\n",finish-start);

  // MPI_Barrier(MPI_COMM_WORLD);
  // start = MPI_Wtime();
  PCShellSetContext(pc, (void *)&context);
  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"PCShellSetContext time: %lf Sec\n",finish-start);

  // finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  // if(myid==0)fprintf(fp,"KSPSet time: %lf Sec\n",finish-start);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  ierr = KSPSolve(ksp, b, x);
  CHKERRQ(ierr);

  finish = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "KSPSolve time: %lf Sec\n", finish - start);

  /*
   Check the error
   */
  // if(myid==0)
  // {
  //     ierr = KSPGetIterationNumber(ksp,&iters);CHKERRQ(ierr);
  //     PetscReal rnorm;
  //     ierr = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(ierr);
  //     ierr = PetscPrintf(PETSC_COMM_WORLD, "residual norm = %f\n", rnorm);
  //     ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n",
  //     iters);
  // }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check the solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPGetIterationNumber(ksp, &iters);
  CHKERRQ(ierr);
  PetscReal rnorm, bnorm;
  PetscCall(VecNorm(b, NORM_2, &bnorm));
  ierr = KSPGetResidualNorm(ksp, &rnorm);
  CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD, "PETSc_FIM_Solver done! number
  // iterations = %d, residual norm = %f\n", iters, rnorm);
  if (myid == 0 && PrintLvl > 0) {
    printf("iter = %d, residual norm = %e, RRN = %e\n", iters, rnorm,
           rnorm / bnorm);
    // write files
    fprintf(fp,
            "PETSc_FIM_Solver done! iter = %d, residual norm = %e, RRN = %e\n",
            iters, rnorm, rnorm / bnorm);
    fprintf(fp, "----------------------------------------------\n");
    fclose(fp);
  }

  // ierr = PetscPrintf(PETSC_COMM_WORLD, "number iterations = %d\n", iters);
  int iter_done = iters;

  int *xIndex = (int *)malloc(rowWidth * sizeof(int));
  for (i = 0; i < rowWidth; i++) {
    xIndex[i] = i + Istart * blockSize; // !!!!
  }
  VecGetValues(x, rowWidth, xIndex, sol);

  /////////////////////////////////////////////////////////

  // for (i=0; i<num_procs; ++i)
  // {
  //     allSize[i]=(allUpper[i]-allLower[i]+1)*nb;
  //     allDisp[i]=allLower[i]*nb;
  // }
  // MPI_Gatherv (local_sol, local_size*nb, MPI_DOUBLE, sol, allSize, allDisp,
  // MPI_DOUBLE, commRoot, MPI_COMM_WORLD);

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
  PetscFinalize();

  return iter_done;
}

// BAMG preconditioner for Hypre software
int FIM_solver_p_bamg(int is_thermal, int myid, int num_procs, int nb,
                      int *allLower, int *allUpper, int *rpt, int *cpt,
                      double *val, double *rhs, double *sol) {
  if (myid == 0 && PrintLvl > 0)
    printf("Call_Solver_times: %d\n", Call_Solver_times);

  //! Construct command line parameters for petsc
  int argc = 29, len = 128, i;
  char **args = (char **)malloc(argc * sizeof(char *));
  for (i = 0; i < argc; i++)
    args[i] = (char *)calloc(len, sizeof(char));

  strcpy(args[0], "program");
  strcpy(args[1], "-pc_hypre_boomeramg_numfunctions");
  sprintf(args[2], "%d", nb);
  strcpy(args[3], "-pc_hypre_boomeramg_nodal_coarsen"); // 0-6
  if (Is_Convergence)
  // if (Is_Convergence)
  {
    // if (Call_Solver_times < 40)
    //     strcpy(args[4], "10"); // nodal coarsening strategy,
    //     nodal_coarsen_diag = 2
    // else
    // strcpy(args[4], "10"); // nodal coarsening strategy, nodal_coarsen_diag =
    // 2
    strcpy(args[4], "0");
    // strcpy(args[4], "11");
    //  * The default is 0 (unknown-based coarsening,
    //  *                   only coarsens within same function).
    //  * For the remaining options a nodal matrix is generated by
    //  * applying a norm to the nodal blocks and applying the coarsening
    //  * algorithm to this matrix.
    //  *    - 1 : Frobenius norm
    //  *    - 2 : sum of absolute values of elements in each block
    //  *    - 3 : largest element in each block (not absolute value)
    //  *    - 4 : row-sum norm, i.e., inf. norm
    //  *    - 5 : col-sum norm, i.e., 1 norm, Li Zhao 09/04/2023
    //  *    - 6 : sum of all values in each block
    //  *    user defined weights in each block for reservior, Li Zhao
    //  07/09/2023
    //  *    - 7 : |PP|
    //  *    - 8 : |TT|, (n,n)
    //  *    - 9 : |PP| + |TT|
    //  *    - 10 : PP + PT + TP + TT = F22
    //  *    - 11 : |TT|, (2,2)
    //  *    - 17 : PP, sign(PP)*|PP|
  } else
    strcpy(args[4], "0");

  strcpy(
      args[5],
      "-pc_hypre_boomeramg_nodal_coarsen_diag"); // {"CLJP", "Ruge-Stueben", "",
                                                 // "modifiedRuge-Stueben", "",
                                                 // "", "Falgout", "", "PMIS",
                                                 // "", "HMIS"};
  strcpy(args[6], "0"); // The default is 0; 1: the diagonal entry is set to the
                        // negative sum of all off diagonal entries; 2: the
                        // signs of all diagonal entries are inverted.

  // sprintf(args[5], "-ksp_monitor");
  strcpy(args[7],
         "-pc_hypre_boomeramg_coarsen_type"); // {"CLJP", "Ruge-Stueben", "",
                                              // "modifiedRuge-Stueben", "", "",
                                              // "Falgout", "", "PMIS", "",
                                              // "HMIS"};
  strcpy(args[8], "HMIS");                    // HMIS
  strcpy(args[9], "-pc_hypre_boomeramg_interp_type"); // ext+i
  strcpy(args[10], "ext+i"); // default: classical, ext+i,  block block-wtd:
                             // no-ok/no-convergence
  // -pc_hypre_boomeramg_interp_type. Available options:
  // classical   direct multipass multipass-wts ext+i ext+i-cc standard
  // standard-wts block block-wtd ff ff1 ext ad-wts ext-mm ext+i-mm ext+e-mm
  // strcpy(args[5], "-pc_hypre_boomeramg_print_debug");
  // strcpy(args[9], "-pc_hypre_boomeramg_vec_interp_qmax");
  // strcpy(args[10], "4");
  strcpy(args[11], "-pc_hypre_boomeramg_relax_type_all"); // Relax type for the
                                                          // up and down cycles
  strcpy(args[12], "symmetric-SOR/Jacobi"); // Jacobi, symmetric-SOR/Jacobi
  // {"Jacobi", "sequential-Gauss-Seidel", "seqboundary-Gauss-Seidel",
  // "SOR/Jacobi", "backward-SOR/Jacobi", "" /* [5] hybrid chaotic Gauss-Seidel
  // (works only with OpenMP) */, "symmetric-SOR/Jacobi", "" /* 7 */,
  // "l1scaled-SOR/Jacobi", "Gaussian-elimination", "" /* 10 */, "" /* 11 */, ""
  // /* 12 */, "l1-Gauss-Seidel" /* nonsymmetric */, "backward-l1-Gauss-Seidel"
  // /* nonsymmetric */, "CG" /* non-stationary */, "Chebyshev", "FCF-Jacobi",
  // "l1scaled-Jacobi"};
  strcpy(args[13],
         "-pc_hypre_boomeramg_relax_type_coarse"); // Relax type on coarse grid
  strcpy(args[14], "Gaussian-elimination");        // Gaussian-elimination, CG
  strcpy(args[15],
         "-pc_hypre_boomeramg_P_max"); // Max elements per row for interpolation
                                       // operator (0=unlimited)
  strcpy(args[16], "4");               // 4
  strcpy(args[17],
         "-pc_hypre_boomeramg_truncfactor"); // Truncation factor for
                                             // interpolation (0=no truncation)
  strcpy(args[18], "0.00");                  // 0, 0.01, 0.05, 0.1
  strcpy(args[19],
         "-pc_hypre_boomeramg_min_coarse_size"); // Minimum size of coarsest
                                                 // grid, default: 1
  strcpy(args[20], "1");
  strcpy(args[21],
         "-pc_hypre_boomeramg_max_coarse_size"); // Maximum size of coarsest
                                                 // grid, default: 9
  strcpy(args[22], "100");
  strcpy(args[23],
         "-pc_hypre_boomeramg_grid_sweeps_coarse"); // Number of sweeps for the
                                                    // coarse level, default: 1
  strcpy(args[24], "1");

  /* complex smoothers */
  // strcpy(args[25], "-pc_hypre_boomeramg_smooth_type"); // Enable more complex
  // smoothers: {"Schwarz-smoothers", "Pilut", "ParaSails", "Euclid"};
  // strcpy(args[26], "Euclid");
  // strcpy(args[27], "-pc_hypre_boomeramg_smooth_num_levels"); // Number of
  // levels on which more complex smoothers are used, num - 1 strcpy(args[28],
  // "1"); strcpy(args[29], "-pc_hypre_boomeramg_eu_level"); // Number of levels
  // for ILU(k) for Euclid, k = ... strcpy(args[30], "0"); // ILU(0)
  // strcpy(args[31], "-pc_hypre_boomeramg_eu_bj");    // Use Block Jacobi for
  // ILU in Euclid smoother? strcpy(args[32], "1");
  strcpy(args[25], "-pc_hypre_boomeramg_strong_threshold"); // theta
  strcpy(args[26], "0.25");
  strcpy(args[27], "-pc_hypre_boomeramg_print_statistics"); // print level
  strcpy(args[28], "0");
  // strcpy(args[23], "-pc_hypre_boomeramg_nodal_relaxation"); // Nodal
  // relaxation via Schwarz
  //! end

  // printf("%s %s %d nb = %d, begin !\n", __FILE__, __FUNCTION__, __LINE__,
  // nb);
  FILE *fp;
  if (myid == 0 && PrintLvl > 0) {
    fp = fopen("linear_solver_BAMG.log", "a");
    fprintf(fp, "Call_Solver_times: %d\n", Call_Solver_times);
  }

  // diag-check
  double start, finish;
  // PetscInitializeNoArguments(); // 这没有命令行参数
  PetscCall(PetscInitialize(&argc, &args, NULL, NULL)); // 有命令行参数
  int nrow_global = allUpper[num_procs - 1] + 1;
  int matrixDim = nrow_global * nb;
  int Istart = allLower[myid];
  int Iend = allUpper[myid];
  int nBlockRows = Iend - Istart + 1;
  int rowWidth = nBlockRows * nb;
  int nb2 = nb * nb;
  int rhs_size = nBlockRows * nb;
  int nnz = rpt[nBlockRows] - rpt[0];

  int blockSize = nb;
  int *nDCount = (int *)malloc(sizeof(int) * nBlockRows);
  int *nNDCount = (int *)malloc(sizeof(int) * nBlockRows);

  // MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  decoup(val, rhs, rpt, cpt, nb, nBlockRows, DecoupType, is_thermal);

  finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "decoup time: %lf Sec\n", finish - start);

  // MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  //----------------------------------------------------------------------
  // Initializing MAT and VEC
  Vec x, b; /* approx solution, RHS, exact solution */
  Mat A;    /* linear system matrix */
  KSP ksp;  /* linear solver context */
  PC pc;
  // PetscReal      norm;     /* norm of solution error */
  PetscInt Ii, iters;
  PetscErrorCode ierr;

  for (int i = 0; i < nBlockRows; i++) {
    nDCount[i] = 0;
    for (int j = rpt[i]; j < rpt[i + 1]; j++) {
      if (cpt[j] >= Istart && cpt[j] <= Iend) {
        nDCount[i]++;
      }
    }
  }

  for (int i = 0; i < nBlockRows; i++) {
    nNDCount[i] = rpt[i + 1] - rpt[i] - nDCount[i];
  }

  MatCreateBAIJ(PETSC_COMM_WORLD, blockSize, rowWidth, rowWidth, matrixDim,
                matrixDim, 0, nDCount, 0, nNDCount, &A);

  ierr = MatSetFromOptions(A);
  CHKERRQ(ierr);
  ierr = MatSetUp(A);
  CHKERRQ(ierr);

  double *valpt = val;
  int *globalx = (int *)malloc(blockSize * sizeof(int));
  int *globaly = (int *)malloc(blockSize * sizeof(int));
  for (Ii = 0; Ii < nBlockRows; Ii++) {
    for (i = 0; i < blockSize; i++) {
      globalx[i] = (Ii + Istart) * blockSize + i;
    }

    for (int i = rpt[Ii]; i < rpt[Ii + 1]; i++) {

      for (int j = 0; j < blockSize; j++) {
        globaly[j] = cpt[i] * blockSize + j;
      }
      ierr = MatSetValues(A, blockSize, globalx, blockSize, globaly, valpt,
                          INSERT_VALUES);
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
  for (i = 0; i < rowWidth; i++) {
    bIndex[i] = Istart * blockSize + i;
  }

  VecSetValues(b, rowWidth, bIndex, rhs, INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  // dBSRmat_ BMat;
  // BMat.ROW = nBlockRows;
  // BMat.COL = nrow_global;
  // BMat.NNZ = nnz;
  // BMat.nb = nb;
  // BMat.IA = rpt;
  // BMat.JA = cpt;
  // BMat.val = val;

  finish = MPI_Wtime();
  Matrix_Conversion_Time +=
      finish - start; //! matrix conversion time (Global variables)
  // MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "petsc A_tmp assemble time: %lf Sec\n", finish - start);

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A);
  CHKERRQ(ierr);

  ierr = KSPSetType(ksp, KSPFGMRES);
  CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp, TOL, PETSC_DEFAULT, PETSC_DEFAULT, MAXIT);
  CHKERRQ(ierr);

  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCHYPRE);
  PCHYPRESetType(pc, "boomeramg");

  ierr = KSPSetFromOptions(ksp);
  CHKERRQ(ierr);

  // MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  ierr = KSPSolve(ksp, b, x);
  CHKERRQ(ierr);

  finish = MPI_Wtime();
  // MPI_Barrier(MPI_COMM_WORLD);
  if (myid == 0 && PrintLvl > 0)
    fprintf(fp, "KSPSolve time: %lf Sec\n", finish - start);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check the solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPGetIterationNumber(ksp, &iters);
  CHKERRQ(ierr);
  PetscReal rnorm, bnorm;
  PetscCall(VecNorm(b, NORM_2, &bnorm));
  ierr = KSPGetResidualNorm(ksp, &rnorm);
  CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD, "PETSc_FIM_Solver done! number
  // iterations = %d, residual norm = %f\n", iters, rnorm);

  /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  if (Call_Solver_times == 0 && PrintLvl > 0)
    KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // view ksp parameters

  int iter_done = iters;
  if (iter_done < MAXIT) {
    Is_Convergence = 1;
  } else {
    Is_Convergence = 0;
    NonConNum++;
  }

  if (myid == 0 && PrintLvl > 0) {
    printf("iter = %d, residual norm = %e, RRN = %e, nodal_coarsen = %s, "
           "NonConNum = %d\n",
           iters, rnorm, rnorm / bnorm, args[4], NonConNum);
    // write files
    fprintf(fp,
            "PETSc_BAMG: iter = %d, residual norm = %e, RRN = %e, "
            "nodal_coarsen = %s, NonConNum = %d\n",
            iters, rnorm, rnorm / bnorm, args[4], NonConNum);
    fprintf(fp, "----------------------------------------------\n");
    fclose(fp);
  }

  int *xIndex = (int *)malloc(rowWidth * sizeof(int));
  for (i = 0; i < rowWidth; i++) {
    xIndex[i] = i + Istart * blockSize; // !!!!
  }
  VecGetValues(x, rowWidth, xIndex, sol);

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

  free(globalx);
  free(globaly);

  PetscFinalize();

  // free command line parameters
  for (i = 0; i < argc; i++)
    free(args[i]);
  free(args);

  Call_Solver_times++;
  // exit(0);
  return iter_done;
}