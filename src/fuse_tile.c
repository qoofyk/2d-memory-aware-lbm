#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"

#if 0 // This function can satisfy any tile size, however because of "min", the speed will become a little slower
void collideStreamTile(Simulation* sim) {
  int lx = sim->lx, ly = sim->ly;
  // Outer loops.
  for (int outerX = 1; outerX <= lx; outerX += tile) {
    for (int outerY = 1; outerY <= ly ; outerY += tile) {
      // Inner loops.
      int innerX_max = MIN(outerX + tile - 1, lx);
      for (int innerX = outerX; innerX <= innerX_max; ++innerX) {
        
        int innerY_max = MIN(outerY + tile - 1, ly);
        for (int innerY = outerY; innerY <= innerY_max; ++innerY) {

          collideNode(&(sim->lattice[innerX][innerY]));

          for (int iPop = 0; iPop < 9; ++iPop) {
            int nextX = innerX + c[iPop][0];
            int nextY = innerY + c[iPop][1];

            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[innerX][innerY].fPop[iPop];
          }
        }
      }
    }
  }
}
#endif

#if 1
void collideStreamTile(Simulation* sim) {
  int lx = sim->lx, ly = sim->ly;
  for (int iX = 1; iX <= lx; iX += tile) {
    for (int iY = 1; iY <= ly; iY += tile) {
      for (int iix = 0; iix < tile; iix++){
        for (int iiy = 0; iiy < tile; iiy++){

          // step1: collision on this line y
          collideNode(&(sim->lattice[iX+iix][iY+iiy]));

          // step 2: stream from line x-1 to x
          for (int iPop = 0; iPop < 9; ++iPop) {
            int nextX = iX+iix + c[iPop][0];
            int nextY = iY+iiy + c[iPop][1];
            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[iX+iix][iY+iiy].fPop[iPop];
          }
        }
      }
    }
  }
}
#endif

// support tile_size doesn't need to divide exactly to lx or ly
void collideStreamTileOMP(Simulation* sim) {

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  int iX, iY, iPop, iix, iiy;
  int nextX, nextY;

  int tid = omp_get_thread_num();
  int my_lx[2];
  my_lx[0] = 1 + tid * my_domain_H;
  my_lx[1] = (tid + 1) * my_domain_H;

  for (int outerX = my_lx[0]; outerX <= my_lx[1]; outerX += tile) {
    for (int outerY = 1; outerY <= ly ; outerY += tile) {
      // Inner loops.
      int innerX_max = MIN(outerX + tile - 1, my_lx[1]);
      for (int innerX = outerX; innerX <= innerX_max; ++innerX) {

        int innerY_max = MIN(outerY + tile - 1, ly);
        for (int innerY = outerY; innerY <= innerY_max; ++innerY) {
          // fused collision and streaming
          collideNode(&(sim->lattice[innerX][innerY]));
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = innerX + c[iPop][0];
            nextY = innerY + c[iPop][1];

            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[innerX][innerY].fPop[iPop];
          }
        }
      }
    }
  }
}
#else
    printf("No OPENMP used");
#endif
}

#if 0
void collideStreamTileOMP(Simulation* sim) {

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  #pragma omp for schedule(static, my_domain_H / tile)
  for (int iX = 1; iX <= lx; iX += tile) {
    for (int iY = 1; iY <= ly; iY += tile) {
      for (int iix = 0; iix < tile; iix++){
        for (int iiy = 0; iiy < tile; iiy++){

          // step1: collision on this line y
          collideNode(&(sim->lattice[iX+iix][iY+iiy]));

          // step 2: stream from line x-1 to x
          for (int iPop = 0; iPop < 9; ++iPop) {
            int nextX = iX+iix + c[iPop][0];
            int nextY = iY+iiy + c[iPop][1];
            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[iX+iix][iY+iiy].fPop[iPop];
          }
        }
      }
    }
  }
}
#else
    printf("No OPENMP used");
#endif
}
#endif