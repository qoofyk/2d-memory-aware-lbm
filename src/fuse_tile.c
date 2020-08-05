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
      for (int innerX=outerX;
           innerX <= MIN(outerX + tile - 1, lx);
           ++innerX)
      {
        for (int innerY = outerY;
             innerY <= MIN(outerY + tile - 1, ly);
             ++innerY)
        {
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
          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
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
          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
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