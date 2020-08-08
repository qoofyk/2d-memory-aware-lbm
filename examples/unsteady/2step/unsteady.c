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

int main(int argc, char *argv[]) {
  
    __itt_pause();

    char case_name[32], filename[64];
    void (*collision_func)(Simulation *) = NULL;
    void (*stream_func)(Simulation *) = NULL;
    
    sprintf(case_name, "step2-line");
    collision_func=&step2CollideStream;
    stream_func = NULL;

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
    for (int iT = 0; iT < warmUpIter; iT += 2) {
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

        // stream_func(&sim);

      //save after stream
      #ifdef SAVE
        if (iT % tSave == 0) {
          printf("iT=%d, save after stream\n", iT);
          sprintf(filename, "vel_%s_stream_%d.dat", case_name, iT);
          saveVel(&sim, filename);
        }
      #endif

    }

    __itt_resume();
    double t[2];
    t[0] = get_cur_time();

    // the main loop over time steps
    for (int iT = 0; iT < numIter; iT += 2) {

        #ifdef ZGB
        // on the right boundary, outlet condition grad_x u = 0
        updateZeroGradientBoundary();
        #endif

        // step 3,4,5:compute rho, u and  update fp
        collision_func(&sim);

        // step 2: streamming step, No need to swap lattice for 2 steps
        // stream_func(&sim);

        // By default: periodic boundary conditions. In this case,
        //   this is important, because the upper and lower
        //   boundaries are horizontally periodic, so that no
        //   special corner nodes are needed.
        // makePeriodic(&sim);
    }

    t[1] = get_cur_time();
    __itt_pause();

    
    printf("After %d iterations: %f Mega site updates per second, Running time (s) = %f\n", 
      numIter, (1.e-6 * lx * ly * numIter) / (t[1] - t[0]), t[1] - t[0]);

    destructSim(&sim);
    freeData();
}