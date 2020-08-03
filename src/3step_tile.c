#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "assert.h"

void step3CollideStreamTile(Simulation* sim) {
  unsigned int iX, iY, iPop;
  // unsigned int tile = 32;
  int iix, iiy;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;

  for (iX = 1; iX <= sim->lx; iX += tile) {
    for (iY = 1; iY <= sim->ly; iY += tile) {
      for (iix = 0; iix < tile; ++iix){
        for (iiy = 0; iiy < tile; ++iiy){
          // 1st fused collision and streaming
          collide_stream_buf1_to_buf2(sim, iX+iix, iY+iiy);

          if ( (iX+iix) > 1 && (iY+iiy) > 1){
#ifdef ZGB
            //save rho
            if ( (iX+iix) == (lx-1) ){
              //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
            }
            if ( (iX+iix) == lx ){
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
            }
#endif
            // 2nd collision on line x-1, y-1
            collide_stream_buf2_to_buf1(sim, iX+iix-1, iY+iiy-1);

            if ( (iX+iix) > 2 && (iY+iiy) > 2){
              // 3rd collision on line x-2, y-2
              collide_stream_buf1_to_buf2(sim, iX+iix-2, iY+iiy-2);
            }
          }
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
}

void step3CollideStreamTileOMP(Simulation* sim) {
  unsigned int iX, iY, iPop, iix, iiy;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;
  int row_index;
  //compute each thread upper boundary line at iX=thread_block 1st c+s

// #define DEBUG_PRINT

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: total_values)
{
  // ------------------ step1: Prepare ------------------------------//
  // 1st fused collision and streaming on last row of 1 thread block (0, 1, 2)
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = 0; iX <= lx; iX += thread_block){
    for (iY = 1; iY <= ly; ++iY) {
      // 1st fused c&s on me (x,y) and (x-1, y), except row x=0
      if (iX == 0){
        // my up + 1
        // collide_stream_buf1_to_buf2(sim, iX+1, iY);

        collideNode(&(sim->lattice[iX+1][iY]));

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+1][iY].fPop[iPop];
        }

        // my up + 2
        // collide_stream_buf1_to_buf2(sim, iX+2, iY);

        collideNode(&(sim->lattice[iX+2][iY]));

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
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

          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
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

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX-1][iY].fPop[iPop];
        }

        // mine
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

        if (iY > 1){
          // 2nd fused c&s on my left (x,y-1) on row x=lx
          // collide_stream_buf2_to_buf1(sim, iX, iY-1);

          collideNode(&(sim->tmpLattice[iX][iY-1]));

          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
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

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX-1][iY].fPop[iPop];
        }

        // mine
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

        // my up + 1
        // collide_stream_buf1_to_buf2(sim, iX+1, iY);

        collideNode(&(sim->lattice[iX+1][iY]));

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+1][iY].fPop[iPop];
        }

        // my up + 2
        // collide_stream_buf1_to_buf2(sim, iX+2, iY);

        collideNode(&(sim->lattice[iX+2][iY]));

        // #pragma ivdep
        // #pragma vector always
        // #pragma vector nontemporal
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

          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX+1 + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX+1][iY-1].fPop[iPop];
          }

          // 2nd c&s on my left (x, y-1)
          // collide_stream_buf2_to_buf1(sim, iX, iY-1);

          collideNode(&(sim->tmpLattice[iX][iY-1]));

          // #pragma ivdep
          // #pragma vector always
          // #pragma vector nontemporal
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
  #pragma omp for private(iix, iiy, iX, iY, row_index, iPop, nextX, nextY) schedule(static, thread_block/tile)
  for (iX = 1; iX <= lx; iX += tile) {
    for (iY = 1; iY <= ly; iY += tile) {
      for (iix = 0; iix < tile; ++iix){
        row_index =  (iX+iix) % thread_block;
        if ((row_index == 1) || (row_index == 2)) continue;

        // from row 3, 4, ..., thread_block-1, 0
        for (iiy = 0; iiy < tile; ++iiy) {

          if (row_index != 0 && row_index != (thread_block-1) ){
            // 1st fused c&s
            // collide_stream_buf1_to_buf2(sim, iX+iix, iY+iiy);

            collideNode(&(sim->lattice[iX+iix][iY+iiy]));

            // #pragma ivdep
            // #pragma vector always
            // #pragma vector nontemporal
            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = iX+iix + c[iPop][0];
              nextY = iY+iiy + c[iPop][1];

              sim->tmpLattice[nextX][nextY].fPop[iPop] =
                sim->lattice[iX+iix][iY+iiy].fPop[iPop];
            }
          }
          if ((iY+iiy) > 1){

            // 2nd fused c&s
            // collide_stream_buf2_to_buf1(sim, iX+iix-1, iY+iiy-1);

            collideNode(&(sim->tmpLattice[iX+iix-1][iY+iiy-1]));

            // #pragma ivdep
            // #pragma vector always
            // #pragma vector nontemporal
            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = iX+iix-1 + c[iPop][0];
              nextY = iY+iiy-1 + c[iPop][1];
              sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop[iPop];
            }

            if ((iY+iiy) > 2){
              // 3rd fused collision and streaming
              // collide_stream_buf1_to_buf2(sim, iX+iix-2, iY+iiy-2);

              collideNode(&(sim->lattice[iX+iix-2][iY+iiy-2]));

              // #pragma ivdep
              // #pragma vector always
              // #pragma vector nontemporal
              for (iPop = 0; iPop < 9; ++iPop) {
                nextX = iX+iix-2 + c[iPop][0];
                nextY = iY+iiy-2 + c[iPop][1];

                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                  sim->lattice[iX+iix-2][iY+iiy-2].fPop[iPop];
              }
            }
          }
        } // end of iiy
      } // end of iix
    }// end of iY loop
  }// end of iX loop
  // ------------------ End step2 ------------------------------//

  // ------------------ step3: handle boundary --------------------//
  // compute 2nd c&s on iY=ly, row_index=2~tb-1 
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, thread_block)
  for (iX = 1; iX <= lx; ++iX){

    // 2nd fused collision and streaming
    iY = ly;
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

  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = thread_block; iX <= lx; iX += thread_block){
    for (iY = 1; iY <= ly-2; iY++){
      // 3rd fused collision and streaming
      // collide_stream_buf1_to_buf2(sim, iX-1, iY);
      // collide_stream_buf1_to_buf2(sim, iX, iY);

      collideNode(&(sim->lattice[iX-1][iY]));

      // #pragma ivdep
      // #pragma vector always
      // #pragma vector nontemporal
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX-1 + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX-1][iY].fPop[iPop];
      }

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
  }

  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, thread_block)
  for (iX = 1; iX <= lx; ++iX){
    // 3rd fused collision and streaming
    iY = ly;

    // collide_stream_buf1_to_buf2(sim, iX, iY-1);
    // collide_stream_buf1_to_buf2(sim, iX, iY);

    collideNode(&(sim->lattice[iX][iY-1]));

    // #pragma ivdep
    // #pragma vector always
    // #pragma vector nontemporal
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY-1 + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY-1].fPop[iPop];
    }

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
  // ------------------ End step3 ------------------------------//
}
#else
    printf("No OPENMP used");
#endif

}// end of func