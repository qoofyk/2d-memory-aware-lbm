#include "lb.h"
#include "boundaries.h"
#include <stdio.h>
#include <stdlib.h>

void setConstants(int argc, char *argv[]) {
  if (argc != 6) {
   printf("%s: Wrong parameters. The syntax is x, y, warmUpIter, numIter, tile\n", argv[0]); 
   exit(1);
  }

  lx        = atoi(argv[1]);     // channel lenghth
  ly        = atoi(argv[2]);     // channel width
  warmUpIter= atoi(argv[3]);     // warm up iterations
  numIter   = atoi(argv[4]);     // total number of iterations
  tile      = atoi(argv[5]);

  obst_r = 8;   // radius of the cylinder // test performance
  // obst_r = ly/10+1;   // radius of the cylinder // for visualization
  obst_x = lx/4;      // position of the cylinder; the cylinder is
  obst_y = ly/2;      // offset from the center to break symmetry

  uMax  = 0.02;       // maximum velocity of the Poiseuille inflow
  Re    = 100;        // Reynolds number
  nu    = uMax * 2.*obst_r / Re;  // kinematic fluid viscosity
  omega = 1. / (3*nu+1./2.);      // relaxation parameter

  tSave  = 2;          // frequency of periodic saves to disk

  printf("lx=%d, ly=%d, Memory=%f MB, omega=%f, warmUpIter=%d, numIter=%d, tile=%d\n", 
    lx, ly, (lx * ly * 80) / 1024.0 / 1024, omega, warmUpIter, numIter, tile);

  #ifdef _OPENMP
    NUM_THREADS = atoi(getenv("OMP_NUM_THREADS"));
    if (lx % NUM_THREADS != 0) {
      printf("lx % NUM_THREADS != 0\n");
      exit(1);
    }
    my_domain_H = lx / NUM_THREADS;
    printf("my_domain_H=%d, NUM_THREADS=%d\n", my_domain_H, NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
      Hello();
  #endif
}

  // Memory allocation and default initialisation of the simulation 
  //   and the left/right boundaries
void iniData() {
    poiseuilleBoundary =
        (VelocityBCData*) calloc(ly+2, sizeof(VelocityBCData));
    pressureBoundary =
        (PressureBCData*) calloc(ly+2, sizeof(PressureBCData));

    leftBoundary  = (Dynamics*) calloc(ly+2, sizeof(Dynamics));
    rightBoundary = (Dynamics*) calloc(ly+2, sizeof(Dynamics));

#ifdef ZGB
  //add by Yuankun
  myrho1 = (double *) calloc(ly, sizeof(double));
  myrho2 = (double *) calloc(ly, sizeof(double));
#endif

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for (int iY=2; iY<=ly-1; ++iY) {
        poiseuilleBoundary[iY].bulkDynamics   = &bulkDynamics;
        poiseuilleBoundary[iY].uy             = 0.;
        pressureBoundary[iY].bulkDynamics = &bulkDynamics;
        pressureBoundary[iY].rho          = 1.;
        pressureBoundary[iY].uPar         = 0.;

        leftBoundary[iY].dynamicsFun = &leftRegularized;
        leftBoundary[iY].selfData
            = (void*)&poiseuilleBoundary[iY];
        rightBoundary[iY].dynamicsFun = &rightPressureRegularized;
        rightBoundary[iY].selfData
            = (void*)&pressureBoundary[iY];
    }
}

  // De-allocation of the memory
void freeData() {
    free(rightBoundary);
    free(leftBoundary);
    free(poiseuilleBoundary);
    free(pressureBoundary);

    //add by Yuankun
    free(myrho1);
    free(myrho2);
}

  // compute parabolic Poiseuille profile
double computePoiseuille(int iY) {
    double y = (double)(iY-1);
    double L = (double)(ly-1);
    return 4.*uMax / (L*L) * (L*y-y*y);
}

  // Specify the geometry of the simulation
void iniGeometry() {
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for(int iX=1; iX<=lx; ++iX) {
        for(int iY=1; iY<=ly; ++iY) {
              // All bulk nodes are initialized at equilibrium with constant
              // density and a velocity determined by a y-dependend Poiseuille
              // profile.
            double uPoiseuille = computePoiseuille(iY);
            iniEquilibrium(&sim.lattice[iX][iY], 1., uPoiseuille, 0.);
              // on the obstacle, set bounce back dynamics
            if ( (iX-obst_x)*(iX-obst_x) +
                 (iY-obst_y)*(iY-obst_y) <= obst_r*obst_r )
            {
                setDynamics(&sim, iX, iY, &bounceBackDynamics);
            }
              // elsewhere, use lbgk dynamics
            else {
                setDynamics(&sim, iX, iY, &bulkDynamics);
            }
        }
    }

    // upper and lower boundary: u=0
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for (int iX=1; iX<=lx; ++iX) {
        setDynamics(&sim, iX, 1, &lowerBoundary);
        setDynamics(&sim, iX, ly, &upperBoundary);
    }

      // left boundary: Poiseuille profile, constant through time
      // right boundary: initially Poiseuille profile, then outlet
      //   condition grad_x u = 0
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for (int iY=2; iY<=ly-1; ++iY) {
        double uPoiseuille = computePoiseuille(iY);
        poiseuilleBoundary[iY].ux   = uPoiseuille;
        setDynamics(&sim, 1, iY, &leftBoundary[iY]);
        setDynamics(&sim, lx, iY, &rightBoundary[iY]);
    }
}

  // Compute a second order extrapolation on the right boundary to
  // ensure a zero-gradient boundary condition on the pressure.
  // This must be recomputed at every time step. The velocity is
  // constrained to be perpendicular to the outflow surface.
void updateZeroGradientBoundary() {
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for (int iY=2; iY<=ly-1; ++iY) {
        double rho1, ux1, uy1, rho2, ux2, uy2;
        computeMacros(sim.lattice[lx-1][iY].fPop, &rho1, &ux1, &uy1);
        computeMacros(sim.lattice[lx-2][iY].fPop, &rho2, &ux2, &uy2);
        pressureBoundary[iY].rho = 4./3.*rho1 - 1./3.*rho2;
        pressureBoundary[iY].uPar = 0.;
    }
}

//OpenMP test hello
void Hello(){

#ifdef _OPENMP
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
#else
  int my_rank = 0;
  int thread_count = 1;
#endif

  printf("Hello from thread %d of %d\n", my_rank, thread_count);
}