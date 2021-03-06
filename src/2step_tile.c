#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "boundaries.h"

#if 1
void step2CollideStreamTile(Simulation* sim) {
  int lx=sim->lx, ly=sim->ly;
  int iX, iY;

  for (iX = 1; iX <= sim->lx; iX += tile) {
    for (iY = 1; iY <= sim->ly; iY += tile) {
      for (int iix = 0; iix < tile; iix++){
        for (int iiy = 0; iiy < tile; iiy++){

          // step1: collision on this line y
          collideNode(&(sim->lattice[iX+iix][iY+iiy]));

          // step 2: stream from line x-1 to x
          // #pragma vector always
          // #pragma vector nontemporal
          #pragma ivdep
          for (int iPop = 0; iPop < 9; ++iPop) {
            int nextX = iX+iix + c[iPop][0];
            int nextY = iY+iiy + c[iPop][1];
            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[iX+iix][iY+iiy].fPop[iPop];
          }

          if (iX+iix>1 && iY+iiy>1){
#ifdef ZGB
            double ux1, uy1, ux2, uy2;
            //save rho
            if( (iX+iix)==(lx-1) ){
              //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
            }
            if( (iX+iix)==lx ){
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
            }
#endif
            // step 3: second collision on line x-1, y-1
            // should be based on the result of first stream
            // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
            collideNode(&(sim->tmpLattice[iX+iix-1][iY+iiy-1]));

            // another branch for iX=sim->lx-1 and iY=sim-lx-2

            // step 4: second stream from  line y-1
            // #pragma vector always
            // #pragma vector nontemporal
            #pragma ivdep
            for (int iPop = 0; iPop < 9; ++iPop) {
              int nextX = iX + iix - 1 + c[iPop][0];
              int nextY = iY + iiy - 1 + c[iPop][1];
              sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop[iPop];
            }
          }
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

  //Line iY=1~ly, iX=lx need to compute one more time
  // iX=sim->lx;
  //simple optimize
  iY=1;
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
    //Compute a second order extrapolation on the right boundary
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

  //compute lx, ly point
  collideNode(&(sim->tmpLattice[lx][ly]));

  // #pragma vector always
  // #pragma vector nontemporal
  #pragma ivdep
  for (int iPop = 0; iPop < 9; ++iPop) {
    int nextX = lx + c[iPop][0];
    int nextY = ly + c[iPop][1];
    sim->lattice[nextX][nextY].fPop[iPop] =
      sim->tmpLattice[lx][ly].fPop[iPop];
  }

}
#endif

#if 1 // use unrolling and remove %, support tile_size doesn't need to divide exactly to lx or ly
void step2CollideStreamTileOMP(Simulation* sim) {
  int lx = sim->lx, ly = sim->ly;

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  #ifdef ZGB
  double ux1, uy1, ux2, uy2;
  #endif
  
  int iX, iY;
  int tid = omp_get_thread_num();
  int my_lx[2];
  my_lx[0] = 1 + tid * my_domain_H;
  my_lx[1] = (tid + 1) * my_domain_H;

  // ------------------ step I : Preprocess intersection area ------------------------------//
  // 1st fused collision and streaming on last row of 1 thread's data domain
  iX = my_lx[1];
  for (iY = 1; iY <= ly; ++iY) {
    collideNode(&(sim->lattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }
  #pragma omp barrier
  // ------------------ End step I ------------------------------//

  // ------------------ step II: Main Area Computation --------------------//
  // Outer loops.
  for (int outerX = my_lx[0]; outerX <= my_lx[1]; outerX += tile) {
    for (int outerY = 1; outerY <= ly ; outerY += tile) {
      // Inner loops.
      int innerX_max = MIN(outerX + tile - 1, my_lx[1]);
      for (int innerX = outerX; innerX <= innerX_max; ++innerX) {

        int innerY_max = MIN(outerY + tile - 1, ly);

        if (innerX == my_lx[1]) {
          for (int innerY = outerY; innerY <= innerY_max; ++innerY) {
            if (innerY > 1){
              #ifdef ZGB
              // save rho
              if ( innerX == lx ){
                computeMacros(sim->tmpLattice[innerX-1][innerY-1].fPop, &myrho1[innerY-1], &ux1, &uy1);
              }
              #endif

              // 2nd fused c&s
              collideNode(&(sim->tmpLattice[innerX - 1][innerY - 1]));

              #pragma ivdep
              for (int iPop = 0; iPop < 9; ++iPop) {
                int nextX = innerX - 1 + c[iPop][0];
                int nextY = innerY - 1 + c[iPop][1];
                sim->lattice[nextX][nextY].fPop[iPop] =
                  sim->tmpLattice[innerX - 1][innerY - 1].fPop[iPop];
              }
            }
          }
        }
        else if (innerX == my_lx[0]) {
          for (int innerY = outerY; innerY <= innerY_max; ++innerY) {
            // 1st fused collision and streaming
            collideNode(&(sim->lattice[innerX][innerY]));

            #pragma ivdep
            for (int iPop = 0; iPop < 9; ++iPop) {
              int nextX = innerX + c[iPop][0];
              int nextY = innerY + c[iPop][1];

              sim->tmpLattice[nextX][nextY].fPop[iPop] =
                sim->lattice[innerX][innerY].fPop[iPop];
            }
          }
        }
        else {
          for (int innerY = outerY; innerY <= innerY_max; ++innerY) {
            // 1st fused collision and streaming
            collideNode(&(sim->lattice[innerX][innerY]));

            #pragma ivdep
            for (int iPop = 0; iPop < 9; ++iPop) {
              int nextX = innerX + c[iPop][0];
              int nextY = innerY + c[iPop][1];

              sim->tmpLattice[nextX][nextY].fPop[iPop] =
                sim->lattice[innerX][innerY].fPop[iPop];
            }

            if (innerY > 1) {
              #ifdef ZGB
              // save rho
              if ( innerX == (lx - 1) ){
                //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                computeMacros(sim->tmpLattice[innerX-1][innerY-1].fPop, &myrho2[innerY-1], &ux2, &uy2);
              }
              #endif

              // 2nd fused c&s
              collideNode(&(sim->tmpLattice[innerX - 1][innerY - 1]));

              #pragma ivdep
              for (int iPop = 0; iPop < 9; ++iPop) {
                int nextX = innerX - 1 + c[iPop][0];
                int nextY = innerY - 1 + c[iPop][1];
                sim->lattice[nextX][nextY].fPop[iPop] =
                  sim->tmpLattice[innerX - 1][innerY - 1].fPop[iPop];
              }
            }
          }
        }

      } // end of innerX loop
    }// end of iY loop
  }// end of iX loop
  #pragma omp barrier
  // ------------------ End step II ------------------------------//

  // ------------------ step III: handle intersection area --------------------//
  // III.1 Righmost Column iX=my_lx[0], iY=ly need to compute one more time
  iX = my_lx[0]; iY = ly;
  // 2nd collision and streaming
  collideNode(&(sim->tmpLattice[iX][iY]));

  #pragma ivdep
  for (int iPop = 0; iPop < 9; ++iPop) {
    int nextX = iX + c[iPop][0];
    int nextY = iY + c[iPop][1];
    sim->lattice[nextX][nextY].fPop[iPop] =
      sim->tmpLattice[iX][iY].fPop[iPop];
  }

  #pragma omp barrier
  
  // III.2 Other data are within each data domain: compute 2nd c+s
  if (my_lx[1] != lx) {
    // III.2.1.1 Row my_lx[1] at iX = my_domain_H 
    iX = my_lx[1];
    for (iY = 1; iY <= (ly-1); ++iY) {
      // 2nd collision and streaming
      collideNode(&(sim->tmpLattice[iX][iY]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
          sim->tmpLattice[iX][iY].fPop[iPop];
      }
    }

    // III.2.1.2 compute rightmost column
    iY = ly;
    for (iX = my_lx[0] + 1; iX <= my_lx[1]; ++iX) {
      // 2nd collision and streaming
      collideNode(&(sim->tmpLattice[iX][iY]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
          sim->tmpLattice[iX][iY].fPop[iPop];
      }
    }
  }
  else {
    // III.2.2.1 Row my_lx[1] at iX = lx 2nd c&s
    iX = lx; iY = 1;
    collideNode(&(sim->tmpLattice[iX][iY]));

    #pragma ivdep
    for (int iPop = 0; iPop < 9; ++iPop) {
      int nextX = iX + c[iPop][0];
      int nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }

    for (iY = 2; iY < ly; ++iY){
      #ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
      #endif
        // 2nd collision and streaming
        collideNode(&(sim->tmpLattice[iX][iY]));

        #pragma ivdep
        for (int iPop = 0; iPop < 9; ++iPop) {
          int nextX = iX + c[iPop][0];
          int nextY = iY + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX][iY].fPop[iPop];
        }
    }

    // III.2.2.1 rightmost column 2nd c&s
    iY = ly;
    for (iX = my_lx[0] + 1; iX <= my_lx[1]; ++iX) {
      collideNode(&(sim->tmpLattice[iX][iY]));

      #pragma ivdep
      for (int iPop = 0; iPop < 9; ++iPop) {
        int nextX = iX + c[iPop][0];
        int nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
          sim->tmpLattice[iX][iY].fPop[iPop];
      }
    }
  }
}
#else
    printf("No OPENMP used\n");
#endif
}
#endif

#if 0
void step2CollideStreamTileOMP(Simulation* sim) {
  unsigned int iX, iY, iPop;
  // unsigned int tile = 32;
  int iix, iiy;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;

  //compute each thread upper boundary line at iX=my_domain_H 1st c+s
#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = my_domain_H; iX <= lx; iX+=my_domain_H){
    for (iY = 1; iY <= ly; ++iY) {

      #ifdef DEBUG_PRINT
        int my_rank = omp_get_thread_num();
        printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
        fflush(stdout);
      #endif

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
    }
  }

  // int schedule_thread_chunk = sim->lx/tile/NUM_THREADS;
  // printf("sim->lx=%d, tile=%d, NUM_THREADS=%d, \n", sim->lx, tile, NUM_THREADS);
  // fflush(stdout);
  #pragma omp for private(iix, iiy, iX, iY, iPop, nextX, nextY) schedule(static, my_domain_H/tile)
  for (iX = 1; iX <= lx; iX+=tile) {
    for (iY = 1; iY <= ly; iY+=tile) {
      for (iix = 0; iix < tile; ++iix){
        for (iiy = 0; iiy < tile; ++iiy){

          if ( (iX+iix) % my_domain_H != 0){

              #ifdef DEBUG_PRINT
                #ifdef _OPENMP
                  int my_rank = omp_get_thread_num();
                  printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX+iix, iY+iiy);
                  fflush(stdout);
                #endif
              #endif

            // 1st fused collision and streaming
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

          if( iX+iix > 1 && iY+iiy > 1){
#ifdef ZGB
            //save rho
            if ( (iX+iix) == (lx-1) ){
              //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
            }
            if( (iX+iix)==lx ){
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
            }
#endif
            if ( (iX+iix-1) % my_domain_H != 0){

              #ifdef DEBUG_PRINT
                  int my_rank = omp_get_thread_num();
                  printf("T%d: 2nd c+s on iX=%d, iY=%d will compute iX-1=%d, iY-1=%d\n",
                      my_rank, iX+iix, iY+iiy, iX+iix-1, iY+iiy-1);
                  fflush(stdout);
              #endif

              // 2nd collision and streaming on line x-1, y-1
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
            }
          }
        }
      }
    }// end of iY loop
  }// end of iX loop

  //compute thread boundary line at iX=my_domain_H 2nd c+s
  //NOTICE: 1~ly-1 !!! use tmpLattice !!!
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = my_domain_H; iX < lx; iX += my_domain_H){
    for (iY = 1; iY <= (ly-1); ++iY) {

      #ifdef DEBUG_PRINT
        int my_rank = omp_get_thread_num();
        printf("T%d: 2nd c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
        fflush(stdout);
      #endif

      // 2nd collision and streaming
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
  }

  //Line iX=1~lx-1, y=ly need to compute one more time
  iY=ly;
  #pragma omp for private(iX, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX < lx; ++iX){
    // 2nd collision and streaming
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

  // Line iY=1~ly, iX=lx need to compute one more time
  // iX=sim->lx;
  // simple optimize
  #pragma omp single
  {
    iX = lx; iY = 1;
    // 2nd collision and streaming
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

  iX = lx;
  #pragma omp for private(iY, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iY = 2; iY < ly; ++iY){

#ifdef ZGB
    //Compute a second order extrapolation on the right boundary
    pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
    pressureBoundary[iY].uPar = 0.;
#endif
    // 2nd collision and streaming
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

  //compute lx, ly point
  #pragma omp single
  {
    iY = ly;
    // 2nd collision and streaming
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
}
#else
    printf("No OPENMP used");
#endif
}
#endif

#if 0 // use inline become slower
void step2CollideStreamTile(Simulation* sim) {
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
            if( (iX+iix) == (lx-1) ){
              //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
            }
            if( (iX+iix) == lx ){
              computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
            }
#endif
            // 2nd collision on line x-1, y-1
            collide_stream_buf2_to_buf1(sim, iX+iix-1, iY+iiy-1);
          }
        }
      }
    }// end of iY loop
  }// end of iX loop

  // -----------Start Epilog for boundary---------------------
  // Line iX=1~lx-1, y=ly need to compute one more time
  iY = ly;
  for (iX = 1; iX < sim->lx; ++iX){
    collide_stream_buf2_to_buf1(sim, iX, iY);
  }

  // Line iY=1~ly, iX=lx need to compute one more time
  // iX=sim->lx;
  //simple optimize
  iX = lx; iY = 1;
  collide_stream_buf2_to_buf1(sim, iX, iY);

  for (iY = 2; iY < sim->ly; ++iY){

#ifdef ZGB
    // Compute a second order extrapolation on the right boundary
    pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
    pressureBoundary[iY].uPar = 0.;
#endif
    collide_stream_buf2_to_buf1(sim, iX, iY);
  }

  // compute lx, ly point
  iX = lx; iY = ly;
  collide_stream_buf2_to_buf1(sim, iX, iY);
  // -----------End Epilog for boundary---------------------
}
#endif