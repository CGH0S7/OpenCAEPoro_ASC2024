//
//  preSolver.cpp
//
//
//  Created by Wenchao Guan on 6/26/14.
//
//

#include "PETScSolver.h"

int preSolver(Mat &A, Vec &b, Vec &x, bool zeroGuess, int solverType)
{
    KSP ksp; /* linear solver context */
    PC pc;

    PetscInt its;
    PetscErrorCode ierr;

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
     Create linear solver context
     */
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    if (!zeroGuess)
        KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
     */
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);

    // precondition only
    ierr = KSPSetType(ksp, KSPRICHARDSON);
    CHKERRQ(ierr);

    /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
     we can then directly call any KSP and PC routines to set
     various options.
     - The following two statements are optional; all of these
     parameters could alternatively be specified at runtime via
     KSPSetFromOptions().  All of these defaults can be
     overridden at runtime, as indicated below.
     */
    ierr = KSPSetTolerances(ksp, 1.e-2, PETSC_DEFAULT, PETSC_DEFAULT, 1);
    CHKERRQ(ierr);

    /*
     Set runtime options, e.g.,
     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
     These options will override those specified above as long as
     KSPSetFromOptions() is called _after_ any other customization
     routines.
     */
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);

    KSPGetPC(ksp, &pc);

    switch (solverType)
    {
    case PRE_JACOBI:
        PCSetType(pc, PCJACOBI);
        break;
    case PRE_BJACOBI:
        PCSetType(pc, PCBJACOBI);
        break;
    case PRE_ILU:
        PCSetType(pc, PCILU);
        break;
    case PRE_AMG:
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "boomeramg");
        break;
    case PRE_EU:
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "euclid");
        break;
    case PRE_PILUT:
        PCSetType(pc, PCHYPRE);
        PCHYPRESetType(pc, "pilut");
        break;
    default:
        PCSetType(pc, PCJACOBI);
        break;
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);

    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
     Check the error
     */
    PetscReal rnorm;
    ierr = KSPGetResidualNorm(ksp, &rnorm);
    CHKERRQ(ierr);

    /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
     */
    // int myid;
    // MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    // if (myid==0) {
    // ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);

    // if (zeroGuess) {
    //     ierr = PetscPrintf(PETSC_COMM_WORLD,"Use zero guess. ");
    //     CHKERRQ(ierr);
    // }
    // else{
    //     ierr = PetscPrintf(PETSC_COMM_WORLD,"Use user's guess. ");
    //     CHKERRQ(ierr);
    // }

    // switch (solverType) {
    //     case PRE_JACOBI:
    //         ierr = PetscPrintf(PETSC_COMM_WORLD,"JACOBI: Norm of residual %.12lf iterations %D\n",rnorm,its);
    //         CHKERRQ(ierr);
    //         break;
    //     case PRE_BJACOBI:
    //         ierr = PetscPrintf(PETSC_COMM_WORLD,"BJACOBI: Norm of residual %.12lf iterations %D\n",rnorm,its);
    //         CHKERRQ(ierr);
    //         break;
    //     case PRE_EU:
    //         ierr = PetscPrintf(PETSC_COMM_WORLD,"euclid: Norm of residual %.12lf iterations %D\n",rnorm,its);
    //         CHKERRQ(ierr);
    //         break;
    //     case PRE_AMG:
    //         ierr = PetscPrintf(PETSC_COMM_WORLD,"AMG: Norm of residual %.12lf iterations %D\n",rnorm,its);
    //         CHKERRQ(ierr);
    //         break;
    //     default:
    //         ierr = PetscPrintf(PETSC_COMM_WORLD,"BJACOBI: Norm of residual %.12lf iterations %D\n",rnorm,its);
    //         CHKERRQ(ierr);
    //         break;
    // }
    // }

    /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
     */
    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);

    return 0;
}
