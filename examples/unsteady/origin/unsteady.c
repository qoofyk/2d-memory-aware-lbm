/* unsteady.c:
 * This example examines an unsteady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, where as the outlet implements an outflow condition:
 * grad_x u = 0. At Reynolds numbers around 100, an unstable periodic pattern is
 * created, the Karman vortex street.
 */

#include "lb.h"
#include "boundaries.h"
#include <stdio.h>
#include <stdlib.h>
#include "eval_tools.h"
#include "ittnotify.h"

// #define SAVE

  // These constants define the flow geometry and are commented in
  //   the function setConstants()
int lx, ly;
int obst_x, obst_y, obst_r;
double uMax, Re, nu, omega;
int numIter, tSave, warmUpIter;

/*-------------------------- added by Yuankun ----------------------*/
int NUM_THREADS;
int tile;
double *myrho1, *myrho2;
int iT=0, count;
int my_domain_H;
/*-------------------------- End------------ ----------------------*/

  // The dynamics that are to be specified on different regions of 
  //   the flow: LBGK in the bulk, bounce-back on the cylinder, 
  //   regularized boundary condition on the four domain boundaries.
  //   Alternatively, the Zou/He boundary condition can be used:
  //   replace upperRegularized by upperZouHe etc.
Dynamics       bulkDynamics         = { &bgk, (void*)&omega };
Dynamics       bounceBackDynamics   = { &bounceBack, 0 };

  // The velocity to be imposed on the upper/lower boundary: u=0
VelocityBCData zeroVelocityBoundary = { &bulkDynamics, 0., 0. };

Dynamics upperBoundary = { &upperRegularized,
                           (void*)&zeroVelocityBoundary };
Dynamics lowerBoundary = { &lowerRegularized,
                           (void*)&zeroVelocityBoundary };
Dynamics* leftBoundary;   // Those two objects are initialized in the
Dynamics* rightBoundary;  //   function iniData()

  // These arrays contain the velocities that are to be imposed on 
  //   the inlet (poiseuilleBoundary) and the outlet 
  //   (pressureBoundary) of the channel
VelocityBCData* poiseuilleBoundary;
PressureBCData* pressureBoundary;

  // The main object containing the simulation
Simulation sim;

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

  printf("\nlx=%d, ly=%d, omega=%f, warmUpIter=%d, numIter=%d, tile=%d\n", 
    lx, ly, omega, warmUpIter, numIter, tile);
}

  // Memory allocation and default initialisation of the simulation 
  //   and the left/right boundaries
void iniData() {
    int iX, iY;
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
    for (iY=2; iY<=ly-1; ++iY) {
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
    int iX, iY;

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for(iX=1; iX<=lx; ++iX) {
        for(iY=1; iY<=ly; ++iY) {
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
    for (iX=1; iX<=lx; ++iX) {
        setDynamics(&sim, iX, 1, &lowerBoundary);
        setDynamics(&sim, iX, ly, &upperBoundary);
    }

      // left boundary: Poiseuille profile, constant through time
      // right boundary: initially Poiseuille profile, then outlet
      //   condition grad_x u = 0
    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for (iY=2; iY<=ly-1; ++iY) {
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
    int iY;
    double rho1, ux1, uy1, rho2, ux2, uy2;

    #ifdef _OPENMP
    #pragma omp parallel for default(shared) schedule(static, my_domain_H)
    #endif
    for (iY=2; iY<=ly-1; ++iY) {
        computeMacros(sim.lattice[lx-1][iY].fPop, &rho1, &ux1, &uy1);
        computeMacros(sim.lattice[lx-2][iY].fPop, &rho2, &ux2, &uy2);
        pressureBoundary[iY].rho = 4./3.*rho1 - 1./3.*rho2;
        pressureBoundary[iY].uPar = 0.;
    }
}

//OpenMP test hello
void Hello(void){

#ifdef _OPENMP
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
#else
  int my_rank = 0;
  int thread_count = 1;
#endif

  printf("Hello from thread %d of %d\n", my_rank, thread_count);
}

int main(int argc, char *argv[]) {
  
  __itt_pause();

#ifdef _OPENMP
  NUM_THREADS = atoi(getenv("OMP_NUM_THREADS"));
  my_domain_H = atoi(getenv("OMP_MY_DOMAIN_HEIGHT"));
  printf("my_domain_H=%d, NUM_THREADS=%d\n", my_domain_H, NUM_THREADS);
  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel
    Hello();
#endif

    char case_name[32], filename[64];
    void (*stream_func)(Simulation *) = NULL;
    void (*collision_func)(Simulation *) = NULL;
    
    sprintf(case_name, "origin");
    collision_func=&collide;
    stream_func = &propagate;

    /*----------------Part I. Initialization----------------------*/
    // initialisation of a lx * ly simulation
    setConstants(argc, argv);
    int interval = numIter / 10;

    // allocate spaces for boundries
    iniData();

    // allocae space for nodes and lattice
    // set configurations
    constructSim(&sim, lx, ly);

    // bounce back dynamics in obstacles, else use lbgk(bulk dynamics)
    iniGeometry();
    /*----------------End Part I. -------------------------------*/


    /*----------------Part II. Run Benchmark----------------------*/
    // Run the benchmark once "to warm up the machine".
    for (int iT = 0; iT < warmUpIter; ++iT) {
      #ifdef SAVE
        if (iT % tSave == 0) {
          printf("iT=%d, save before computing\n", iT);
          sprintf(filename, "vel_%s_%d.dat", case_name, iT);
          saveVel(&sim, filename);
        }
      #endif

      #ifdef ZGB
        // on the right boundary, outlet condition grad_x u = 0
        updateZeroGradientBoundary();
      #endif

        //save after ZGB
      #ifdef SAVE
        if (iT % tSave == 0) {
          printf("iT=%d, save after ZGB\n", iT);
          sprintf(filename, "vel_%s_zgb_%d.dat", case_name, iT);
          saveVel(&sim, filename);
        }
      #endif

        // step 3,4,5:compute rho, u, get f^eq, then update f.
        collision_func(&sim);

      //save after collide
      #ifdef SAVE
        if (iT % tSave == 0) {
          sprintf(filename, "vel_%s_collide_%d.dat", case_name, iT);
          saveVel(&sim, filename);
        }
      #endif

        stream_func(&sim);

      //save after stream
      #ifdef SAVE
        if (iT % tSave == 0) {
          printf("iT=%d, save after stream, count=%d\n", iT);
          sprintf(filename, "vel_%s_stream_%d.dat", case_name, count);
          saveVel(&sim, filename);
        }
      #endif

    }

    __itt_resume();
    double t[2];
    t[0] = get_cur_time();

    // the main loop over time steps
    for (int iT = 0; iT < numIter; ++iT) {

        #ifdef ZGB
        // on the right boundary, outlet condition grad_x u = 0
        updateZeroGradientBoundary();
        #endif

        // step 3,4,5:compute rho, u and  update fp
        collision_func(&sim);

        // step 2: streamming step
        stream_func(&sim);

        // By default: periodic boundary conditions. In this case,
        //   this is important, because the upper and lower
        //   boundaries are horizontally periodic, so that no
        //   special corner nodes are needed.
        // makePeriodic(&sim);
    }

    t[1] = get_cur_time();
    __itt_pause();

    
    printf("After %d iterations: %f Mega site updates per second, Running time (s) = %f\n", 
      numIter, (lx * ly * numIter) / (t[1] - t[0]) / 1.e6, t[1] - t[0]);

    destructSim(&sim);
    freeData();
}