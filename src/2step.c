#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "boundaries.h"
#include "assert.h"

// immediately compute iX-1, iY-1
void step2CollideStream(Simulation* sim) {
  int iX, iY;
  int lx=sim->lx, ly=sim->ly;

  for (iX = 1; iX <= sim->lx; ++iX) {
    for (iY = 1; iY <= sim->ly; ++iY) {

      // step1: collision on this line y
      collideNode(&(sim->lattice[iX][iY]));

      // step 2: stream from line x-1 to x
      // #pragma ivdep
      // #pragma vector always
      // #pragma vector nontemporal
      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];
        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }

      if (iX > 1 && iY > 1){

#ifdef ZGB
        double ux1, uy1, ux2, uy2;
        //save rho
        if (iX == (lx-1) ){
          //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
          computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
        }
        if (iX == lx){
          computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
        }
#endif
        // step 3: second collision on line x-1, y-1
        // should be based on the result of first stream
        // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
        collideNode(&(sim->tmpLattice[iX-1][iY-1]));

        // another branch for iX=sim->lx-1 and iY=sim-lx-2

        // step 4: second stream from  line y-1
        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (int iPop = 0; iPop < 9; ++iPop) {
          int nextX = iX-1 + c[iPop][0];
          int nextY = iY-1 + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX-1][iY-1].fPop[iPop];
        }
      }
    }// end of iY loop
  }// end of iX loop

  // Line iX=1~lx-1, y=ly need to compute one more time
  // iY=sim->ly;

  for (iX = 1; iX < sim->lx; ++iX){
    collideNode(&(sim->tmpLattice[iX][ly]));

    // #pragma vector always
    // #pragma vector nontemporal
    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = ly + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][ly].fPop[iPop];
    }
  }

  // Line iY=1~ly, iX=lx need to compute one more time
  // iX=sim->lx;
  //simple optimize
  iY = 1;
  collideNode(&(sim->tmpLattice[lx][iY]));

  // #pragma vector always
  // #pragma vector nontemporal
  #pragma ivdep
  for (int iPop = 0; iPop < 9; ++iPop) {
    int nextX = lx + c[iPop][0];
    int nextY = iY + c[iPop][1];
    sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[lx][iY].fPop[iPop];
  }

  for (iY = 2; iY < sim->ly; ++iY){

#ifdef ZGB
    // Compute a second order extrapolation on the right boundary
    pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
    pressureBoundary[iY].uPar = 0.;
#endif
    collideNode(&(sim->tmpLattice[lx][iY]));

    // #pragma vector always
    // #pragma vector nontemporal
    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = lx + c[iPop][0];
      int nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[lx][iY].fPop[iPop];
    }
  }

  // compute lx, ly point
  collideNode(&(sim->tmpLattice[lx][ly]));

  // #pragma vector always
  // #pragma vector nontemporal
  #pragma ivdep
  for (int iPop=0; iPop<9; ++iPop) {
    int nextX = lx + c[iPop][0];
    int nextY = ly + c[iPop][1];
    sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[lx][ly].fPop[iPop];
  }

}// end of func

//go through whole line, then compute lower line
void step2CollideStream2(Simulation* sim) {
  unsigned int iX, iY, iPop;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;

  for (iX = 1; iX <= sim->lx; ++iX) {

#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iY, iPop, nextX, nextY) \
    schedule(static, my_domain_H)
#endif
    for (iY = 1; iY <= sim->ly; ++iY) {
        collide_stream_buf1_to_buf2(sim, iX, iY);
    }// end of iY loop

#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iY, iPop, nextX, nextY) \
    schedule(static, my_domain_H)
#endif
    for (iY = 1; iY <= sim->ly; ++iY){
      if (iX > 1 && iY > 1){

#ifdef ZGB
        //save rho
        if( iX == (lx-1) ){
          //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
          computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
        }
        if( iX == lx ){
          computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
        }
#endif
        // 2nd fused collision and streaming
        collide_stream_buf2_to_buf1(sim, iX-1, iY-1);
      }
    }

  }// end of iX loop

  //Line iX=1~lx-1, y=ly need to compute one more time
  iY = ly;
#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iX, iPop, nextX, nextY) \
    schedule(static, my_domain_H)
#endif
  for (iX = 1; iX < sim->lx; ++iX){
    // 2nd collision and streaming
    collide_stream_buf2_to_buf1(sim, iX, iY);
  }


  //Line iY=1~ly, iX=lx need to compute one more time
  // iX=sim->lx;
  //simple optimize
  // #pragma omp critical
  // {
      iX = lx; iY = 1;
      // 2nd collision and streaming
      collide_stream_buf2_to_buf1(sim, iX, iY);
  // }

#ifdef _OPENMP
#pragma omp parallel for default(shared) \
  private(iY, iPop, nextX, nextY) \
  schedule(static, my_domain_H)
#endif
  for (iY = 2; iY < sim->ly; ++iY){

#ifdef ZGB
    //Compute a second order extrapolation on the right boundary
    pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
    pressureBoundary[iY].uPar = 0.;
#endif
    // 2nd collision and streaming
    collide_stream_buf2_to_buf1(sim, iX, iY);
  }

  //compute lx, ly point
  // #pragma omp critical
  // {
      iX = lx; iY = ly;
      // 2nd collision and streaming
      collide_stream_buf2_to_buf1(sim, iX, iY);
  // }

}// end of func

// explicit, no use inline
void step2CollideStreamOMP(Simulation* sim) {
  int iX, iY, iPop;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;

  //compute each thread upper boundary line at iX=my_domain_H 1st c+s

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  // compute thread boundaries
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = my_domain_H; iX <= lx; iX += my_domain_H){
    for (iY = 1; iY <= ly; ++iY) {

      #ifdef DEBUG_PRINT
        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
          fflush(stdout);
        #endif
      #endif

      // 1st fused collision and streaming
      // collide_stream_buf1_to_buf2(sim, iX, iY);

      collideNode(&(sim->lattice[iX][iY]));

      // #pragma vector always
      // #pragma vector nontemporal
      #pragma ivdep
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }
    }
  }

  // compute bulk
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX) {
    for (iY = 1; iY <= ly; ++iY) {

      if (iX % my_domain_H != 0){

        #ifdef DEBUG_PRINT
          #ifdef _OPENMP
            int my_rank = omp_get_thread_num();
            printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
            fflush(stdout);
          #endif
        #endif

        // 1st fused collision and streaming
        // collide_stream_buf1_to_buf2(sim, iX, iY);

        collideNode(&(sim->lattice[iX][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX][iY].fPop[iPop];
        }
      }

      if (iX > 1 && iY > 1){

#ifdef ZGB
        //save rho
        if(iX == (lx-1) ){
            //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
            computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
        }
        if(iX == lx){
            computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
        }
#endif
        if ( (iX-1)%my_domain_H != 0){

          #ifdef DEBUG_PRINT
            #ifdef _OPENMP
              int my_rank = omp_get_thread_num();
              printf("T%d: 2nd c+s on iX=%d, iY=%d will compute iX-1=%d, iY-1=%d\n", my_rank, iX, iY, iX-1, iY-1);
              fflush(stdout);
            #endif
          #endif

          // 2nd fused collision and streaming
          // collide_stream_buf2_to_buf1(sim, iX-1, iY-1);

          collideNode(&(sim->tmpLattice[iX-1][iY-1]));

          // #pragma vector always
          // #pragma vector nontemporal
          #pragma ivdep
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX-1 + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX-1][iY-1].fPop[iPop];
          }
        }
      }
    }// end of iY loop
  }// end of iX loop

  //compute thread boundary line at iX=my_domain_H 2nd c+s
  //NOTICE: 1~ly-1 !!! use tmpLattice !!!
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = my_domain_H; iX < lx; iX += my_domain_H){
    for (iY = 1; iY <= (ly - 1); ++iY) {

      #ifdef DEBUG_PRINT
        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          printf("T%d: 2nd c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
          fflush(stdout);
        #endif
      #endif

      // 2nd fused collision and streaming
      // collide_stream_buf2_to_buf1(sim, iX, iY);

      collideNode(&(sim->tmpLattice[iX][iY]));

      // #pragma vector always
      // #pragma vector nontemporal
      #pragma ivdep
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
          sim->tmpLattice[iX][iY].fPop[iPop];
      }
    }
  }


  // Line iX=1~lx-1, y=ly need to compute one more time
  iY = ly;
  #pragma omp for private(iX, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX < lx; ++iX){
    // 2nd fused collision and streaming
    // collide_stream_buf2_to_buf1(sim, iX, iY);

    collideNode(&(sim->tmpLattice[iX][iY]));

    // #pragma vector always
    // #pragma vector nontemporal
    #pragma ivdep
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }

  // -----------Start Epilog for boundary---------------------
  // Line iY=1~ly, iX=lx need to compute one more time
  // iX=sim->lx;
  #pragma omp single
  {
      iX = lx; iY = 1;
      collide_stream_buf2_to_buf1(sim, iX, iY);
  }

  #pragma omp for private(iY, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iY = 2; iY < ly; ++iY){

#ifdef ZGB
    //Compute a second order extrapolation on the right boundary
    pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
    pressureBoundary[iY].uPar = 0.;
#endif
    // collide_stream_buf2_to_buf1(sim, iX, iY);

    collideNode(&(sim->tmpLattice[iX][iY]));

    // #pragma vector always
    // #pragma vector nontemporal
    #pragma ivdep
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }

  #pragma omp single
  {
      //compute lx, ly point
      iX = lx; iY = ly;
      collide_stream_buf2_to_buf1(sim, iX, iY);
  }
  // -----------End Epilog for boundary---------------------
}
#else
    printf("No OPENMP used");
#endif

}// end of func