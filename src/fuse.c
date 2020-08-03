#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "assert.h"

void collideStream(Simulation* sim) {
  int iX, iY, iPop;
  int nextX, nextY;

  // #pragma ivdep
  // #pragma vector always
  // #pragma vector nontemporal
  for (iX = 1; iX <= sim->lx; ++iX) {
    for (iY = 1; iY <= sim->ly; ++iY) {
      collideNode(&(sim->lattice[iX][iY]));
      // streamming imediately once we got the updated f
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }
    }
  }
}

void collideStreamOMP(Simulation* sim) {
  int iX, iY, iPop;
  int nextX, nextY;

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: total_values)
{
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, thread_block)
  for (iX = 1; iX <= sim->lx; ++iX) {
    for (iY = 1; iY <= sim->ly; ++iY) {
      collideNode(&(sim->lattice[iX][iY]));
      // streamming imediately once we got the updated f
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];
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