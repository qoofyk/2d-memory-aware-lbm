#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "assert.h"

void collideStream(Simulation* sim) {
  // #pragma ivdep
  // #pragma vector always
  // #pragma vector nontemporal
  for (int iX = 1; iX <= sim->lx; ++iX) {
    for (int iY = 1; iY <= sim->ly; ++iY) {
      collideNode(&(sim->lattice[iX][iY]));
      // streamming imediately once we got the updated f
      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];
        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }
    }
  }
}

void collideStreamOMP(Simulation* sim) {
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  #pragma omp for schedule(static, my_domain_H)
  for (int iX = 1; iX <= sim->lx; ++iX) {
    for (int iY = 1; iY <= sim->ly; ++iY) {
      collideNode(&(sim->lattice[iX][iY]));
      // streamming imediately once we got the updated f
      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];
        sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX][iY].fPop[iPop];
      }
    }
  }
}
#else
    printf("No OPENMP used");
#endif

}

void swapLattice(Simulation* sim){
  Node** swapLattice = sim->lattice;
  sim->lattice = sim->tmpLattice;
  sim->tmpLattice = swapLattice;
}