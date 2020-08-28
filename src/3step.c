#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "assert.h"

#if 0
// immediately compute (iX-1, iY-1), then (iX-2, iY-2)
void step3CollideStream(Simulation* sim) {
  unsigned int iX, iY, iPop;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;

  for (iX = 1; iX <= sim->lx; ++iX) {
    for (iY = 1; iY <= sim->ly; ++iY) {

      // 1st fused collision and streaming
      // collide_stream_buf1_to_buf2(sim, iX, iY);
      collideNode(&(sim->lattice[iX][iY]));

      // #pragma ivdep
      // #pragma vector always
      // #pragma vector nontemporal
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }

      if (iX > 1 && iY > 1){

        // 2nd fused collision and streaming
        // collide_stream_buf2_to_buf1(sim, iX-1, iY-1);

        collideNode(&(sim->tmpLattice[iX-1][iY-1]));

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY-1 + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX-1][iY-1].fPop[iPop];
        }

        if (iX > 2 && iY > 2){

          // 3rd fused collision and streaming
          // collide_stream_buf1_to_buf2(sim, iX-2, iY-2);

          collideNode(&(sim->lattice[iX-2][iY-2]));

          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX-2 + c[iPop][0];
            nextY = iY-2 + c[iPop][1];

            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[iX-2][iY-2].fPop[iPop];
          }

        }
      }
    }// end of iY loop
  }// end of iX loop

  //-----------------------------------------------------------------------
  // Line iX=[1~lx-1], y=ly need to compute the second collision & stream
  iY = ly;
  for(iX = 1; iX < sim->lx; ++iX){
    // collide_stream_buf2_to_buf1(sim, iX, iY);
    collideNode(&(sim->tmpLattice[iX][iY]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }

  // Line iY=[1~ly], iX=lx need to compute the second collision & stream
  // iX=sim->lx;
  //simple optimize
  iX = lx;
  for (iY = 1; iY <= sim->ly; ++iY){
    // collide_stream_buf2_to_buf1(sim, iX, iY);
    collideNode(&(sim->tmpLattice[iX][iY]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }

  //-----------------------------------------------------------------------
  // Line iX=[1~lx-2], y=ly-1 need to compute the third collision & stream
  iY = ly-1;
  for (iX = 1; iX < (sim->lx-1); ++iX){
    // collide_stream_buf1_to_buf2(sim, iX, iY);
    collideNode(&(sim->lattice[iX][iY]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  // Line iY=[1~ly-1], iX=lx-1 need to compute the third collision & stream
  iX = lx-1;
  for (iY = 1; iY < sim->ly; ++iY){
    // collide_stream_buf1_to_buf2(sim, iX, iY);
    collideNode(&(sim->lattice[iX][iY]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  //-----------------------------------------------------------------------
  // Line iX=[1~lx-1], y=ly need to compute the third collision&stream
  iY = ly;
  for(iX = 1; iX < sim->lx; ++iX){
    // collide_stream_buf1_to_buf2(sim, iX, iY);
    collideNode(&(sim->lattice[iX][iY]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  // Line iY=1~ly, iX=lx need to compute the third collision&stream
  // iX=sim->lx;
  // simple optimize
  iX = lx;
  for (iY = 1; iY <= sim->ly; ++iY){
    // collide_stream_buf1_to_buf2(sim, iX, iY);
    collideNode(&(sim->lattice[iX][iY]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

}// end of func
#endif

#if 1
// immediately compute (iX-1, iY-1), then (iX-2, iY-2) & unrolling
void step3CollideStream(Simulation* sim) {
  int iX, iY;
  int lx = sim->lx, ly = sim->ly;
  double ux1, uy1, ux2, uy2;

  /* ------ line iX = 1 -----*/
  iX = 1;
  for (iY = 1; iY <= ly; ++iY) {
    // 1st fused collision and streaming
    collideNode(&(sim->lattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  /* ------ line iX = 2 -----*/
  iX = 2; 
  // iY = 1
  iY = 1;
  // 1st fused collision and streaming
  collideNode(&(sim->lattice[iX][iY]));

  #pragma ivdep
  for (int iPop = 0; iPop < 9; ++iPop) {
    int nextX = iX + c[iPop][0];
    int nextY = iY + c[iPop][1];

    sim->tmpLattice[nextX][nextY].fPop[iPop] =
      sim->lattice[iX][iY].fPop[iPop];
  }

  // iY = [2, ly-1]
  for (iY = 2; iY <= ly; ++iY) {
    // 1st fused collision and streaming
    collideNode(&(sim->lattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }

    // 2nd fused collision and streaming
    collideNode(&(sim->tmpLattice[iX-1][iY-1]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX-1 + c[iPop][0];
      int nextY = iY-1 + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX-1][iY-1].fPop[iPop];
    }
  }

  iY = ly;
  // 2nd fused collision and streaming on (iX-1, iY)
  collideNode(&(sim->tmpLattice[iX-1][iY]));

  #pragma ivdep
  for (int iPop = 0; iPop < 9; ++iPop) {
    int nextX = iX-1 + c[iPop][0];
    int nextY = iY + c[iPop][1];
    sim->lattice[nextX][nextY].fPop[iPop] =
      sim->tmpLattice[iX-1][iY].fPop[iPop];
  }

  /* ------ line iX = [3, lx] -----*/
  for (iX = 3; iX <= lx; ++iX) {
    // iY = 1
    iY = 1;
    // 1st fused collision and streaming
    collideNode(&(sim->lattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }

    iY = 2;
    // 1st fused collision and streaming
    collideNode(&(sim->lattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }

    // 2nd fused collision and streaming
    collideNode(&(sim->tmpLattice[iX-1][iY-1]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX-1 + c[iPop][0];
      int nextY = iY-1 + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX-1][iY-1].fPop[iPop];
    }

    // iY = [3, ly-1]
    for (iY = 3; iY <= ly; ++iY) {
      // 1st fused collision and streaming
      collideNode(&(sim->lattice[iX][iY]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }

      // 2nd fused collision and streaming
      collideNode(&(sim->tmpLattice[iX-1][iY-1]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX-1 + c[iPop][0];
        int nextY = iY-1 + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
          sim->tmpLattice[iX-1][iY-1].fPop[iPop];
      }

      // 3rd fused collision and streaming
      collideNode(&(sim->lattice[iX-2][iY-2]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX-2 + c[iPop][0];
        int nextY = iY-2 + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX-2][iY-2].fPop[iPop];
      }
    }// end of iY loop

    iY = ly;
    // 2nd fused collision and streaming
    collideNode(&(sim->tmpLattice[iX-1][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX-1 + c[iPop][0];
      int nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX-1][iY].fPop[iPop];
    }

    // 3rd fused collision and streaming
    collideNode(&(sim->lattice[iX-2][iY-1]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX-2 + c[iPop][0];
      int nextY = iY-1 + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX-2][iY-1].fPop[iPop];
    }

    // 3rd fused collision and streaming
    collideNode(&(sim->lattice[iX-2][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX-2 + c[iPop][0];
      int nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX-2][iY].fPop[iPop];
    }
  }// end of iX loop

  /* ------ line iX = [lx-1, lx] -----*/
  // 2nd fused collision and streaming on iX = lx
  iX = lx;
  for (iY = 1; iY <= ly; ++iY) {
    collideNode(&(sim->tmpLattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  } // end of iY loop
  

  // 3rd fused collision and streaming on iX = lx - 1 & lx
  for (iX = lx - 1; iX <= lx; ++iX) {
    for (iY = 1; iY <= ly; ++iY) {
      collideNode(&(sim->lattice[iX][iY]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }
    } // end of iY loop
  }
  

}// end of func
#endif

void step3CollideStreamOMP(Simulation* sim) {
  unsigned int iX, iY, iPop;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;
  int row_index;

// #define DEBUG_PRINT

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  // ------------------ step1: Prepare ------------------------------//
  // 1st fused collision and streaming on last row of 1 thread block (0, 1, 2)
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = 0; iX <= lx; iX += my_domain_H){
    for (iY = 1; iY <= ly; ++iY) {
      // 1st fused c&s on me (x,y) and (x-1, y), except row x=0
      if (iX == 0){
        // my up + 1
        // collide_stream_buf1_to_buf2(sim, iX+1, iY);

        collideNode(&(sim->lattice[iX+1][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+1][iY].fPop[iPop];
        }

        // my up + 2
        // collide_stream_buf1_to_buf2(sim, iX+2, iY);

        collideNode(&(sim->lattice[iX+2][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+2 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+2][iY].fPop[iPop];
        }

        if (iY > 1){
          // 2nd fused c&s on my left & up point (x+1, y-1)
          // collide_stream_buf2_to_buf1(sim, iX+1, iY-1);

          collideNode(&(sim->tmpLattice[iX+1][iY-1]));

          // #pragma vector always
          // #pragma vector nontemporal
          #pragma ivdep
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX+1 + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX+1][iY-1].fPop[iPop];
          }

        }
      }
      else if (iX == lx){
        // my down - 1
        // collide_stream_buf1_to_buf2(sim, iX-1, iY);

        collideNode(&(sim->lattice[iX-1][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX-1][iY].fPop[iPop];
        }

        // mine
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

        if (iY > 1){
          // 2nd fused c&s on my left (x,y-1) on row x=lx
          // collide_stream_buf2_to_buf1(sim, iX, iY-1);

          collideNode(&(sim->tmpLattice[iX][iY-1]));

          // #pragma vector always
          // #pragma vector nontemporal
          #pragma ivdep
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX][iY-1].fPop[iPop];
          }
        }
      }
      else{
        // my down - 1
        // collide_stream_buf1_to_buf2(sim, iX-1, iY);

        collideNode(&(sim->lattice[iX-1][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX-1][iY].fPop[iPop];
        }

        // mine
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

        // my up + 1
        // collide_stream_buf1_to_buf2(sim, iX+1, iY);

        collideNode(&(sim->lattice[iX+1][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+1][iY].fPop[iPop];
        }

        // my up + 2
        // collide_stream_buf1_to_buf2(sim, iX+2, iY);

        collideNode(&(sim->lattice[iX+2][iY]));

        // #pragma vector always
        // #pragma vector nontemporal
        #pragma ivdep
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+2 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+2][iY].fPop[iPop];
        }

        if (iY > 1){
          // 2nd fused c&s on my left & up point (x+1, y-1)
          // collide_stream_buf2_to_buf1(sim, iX+1, iY-1);

          collideNode(&(sim->tmpLattice[iX+1][iY-1]));

          // #pragma vector always
          // #pragma vector nontemporal
          #pragma ivdep
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX+1 + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX+1][iY-1].fPop[iPop];
          }

          // 2nd c&s on my left (x, y-1)
          // collide_stream_buf2_to_buf1(sim, iX, iY-1);

          collideNode(&(sim->tmpLattice[iX][iY-1]));

          // #pragma vector always
          // #pragma vector nontemporal
          #pragma ivdep
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX][iY-1].fPop[iPop];
          }
        }
      }
    }
  }
  // ------------------ End step1 ------------------------------//

  // ------------------ step2: compute bulk --------------------//
  #pragma omp for private(iX, iY, row_index, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX) {
    row_index =  iX % my_domain_H;
    if ((row_index == 1) || (row_index == 2)) continue;
    // from row 3, 4, ..., my_domain_H-1, 0
    for (iY = 1; iY <= ly; ++iY) {

      if (row_index != 0 && row_index != (my_domain_H-1) ){
        // 1st fused c&s
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
      if (iY > 1){
        // 2nd fused c&s
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

        if (iY > 2){
            // 3rd fused collision and streaming
            // collide_stream_buf1_to_buf2(sim, iX-2, iY-2);

            collideNode(&(sim->lattice[iX-2][iY-2]));

            // #pragma vector always
            // #pragma vector nontemporal
            #pragma ivdep
            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = iX-2 + c[iPop][0];
              nextY = iY-2 + c[iPop][1];

              sim->tmpLattice[nextX][nextY].fPop[iPop] =
                sim->lattice[iX-2][iY-2].fPop[iPop];
            }
        }
      }
    }// end of iY loop
  }// end of iX loop
  // ------------------ End step2 ------------------------------//

  // ------------------ step3: handle boundary --------------------//
  // compute 2nd c&s on iY=ly, row_index=2~tb-1 
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX){

    // 2nd fused collision and streaming
    iY = ly;
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

  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = my_domain_H; iX <= lx; iX += my_domain_H){
    for (iY = 1; iY <= ly-2; iY++){
      // 3rd fused collision and streaming
      // collide_stream_buf1_to_buf2(sim, iX-1, iY);
      // collide_stream_buf1_to_buf2(sim, iX, iY);

      collideNode(&(sim->lattice[iX-1][iY]));

      // #pragma vector always
      // #pragma vector nontemporal
      #pragma ivdep
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX-1 + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX-1][iY].fPop[iPop];
      }

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

  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX){
    // 3rd fused collision and streaming
    iY = ly;

    // collide_stream_buf1_to_buf2(sim, iX, iY-1);
    // collide_stream_buf1_to_buf2(sim, iX, iY);

    collideNode(&(sim->lattice[iX][iY-1]));

    // #pragma vector always
    // #pragma vector nontemporal
    #pragma ivdep
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY-1 + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY-1].fPop[iPop];
    }

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
  // ------------------ End step3 ------------------------------//
}
#else
    printf("No OPENMP used");
#endif

}// end of func

#if 0 // use inline function will become slower!
// immediately compute (iX-1, iY-1), then (iX-2, iY-2)
void step3CollideStream(Simulation* sim) {
  unsigned int iX, iY, iPop;
  int lx=sim->lx, ly=sim->ly;

  for (iX = 1; iX <= sim->lx; ++iX) {
    for (iY = 1; iY <= sim->ly; ++iY) {

      // 1st fused collision and streaming
      collide_stream_buf1_to_buf2(sim, iX, iY);

      if (iX > 1 && iY > 1){
        // 2nd fused collision and streaming
        collide_stream_buf2_to_buf1(sim, iX-1, iY-1);

        if (iX > 2 && iY > 2){
          // 3rd fused collision and streaming
          collide_stream_buf1_to_buf2(sim, iX-2, iY-2);
        }
      }
    }// end of iY loop
  }// end of iX loop

  //-----------------------------------------------------------------------
  // Line iX=[1~lx-1], y=ly need to compute the second collision & stream
  iY = ly;
  for(iX = 1; iX < sim->lx; ++iX){
    collide_stream_buf2_to_buf1(sim, iX, iY);
  }

  // Line iY=[1~ly], iX=lx need to compute the second collision & stream
  // iX=sim->lx;
  //simple optimize
  iX = lx;
  for (iY = 1; iY <= sim->ly; ++iY){
    collide_stream_buf2_to_buf1(sim, iX, iY);
  }

  //-----------------------------------------------------------------------
  // Line iX=[1~lx-2], y=ly-1 need to compute the third collision & stream
  iY = ly-1;
  for (iX = 1; iX < (sim->lx-1); ++iX){
    collide_stream_buf1_to_buf2(sim, iX, iY);
  }

  // Line iY=[1~ly-1], iX=lx-1 need to compute the third collision & stream
  iX = lx-1;
  for (iY = 1; iY < sim->ly; ++iY){
    collide_stream_buf1_to_buf2(sim, iX, iY);
  }

  //-----------------------------------------------------------------------
  // Line iX=[1~lx-1], y=ly need to compute the third collision&stream
  iY = ly;
  for(iX = 1; iX < sim->lx; ++iX){
    collide_stream_buf1_to_buf2(sim, iX, iY);
  }

  // Line iY=1~ly, iX=lx need to compute the third collision&stream
  // iX=sim->lx;
  // simple optimize
  iX = lx;
  for (iY = 1; iY <= sim->ly; ++iY){
    collide_stream_buf1_to_buf2(sim, iX, iY);
  }

}// end of func
#endif