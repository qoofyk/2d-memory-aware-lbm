#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"

#if 1
void step3CollideStreamTile(Simulation* sim) {
  int iX, iY, iPop;
  int iix, iiy;
  int lx = sim->lx, ly = sim->ly;
  int nextX, nextY;

  for (iX = 1; iX <= lx; iX += tile) {
    for (iY = 1; iY <= ly; iY += tile) {
      for (iix = 0; iix < tile; ++iix){
        for (iiy = 0; iiy < tile; ++iiy){
          // 1st fused collision and streaming
          collideNode(&(sim->lattice[iX+iix][iY+iiy]));

          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX+iix + c[iPop][0];
            nextY = iY+iiy + c[iPop][1];

            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[iX+iix][iY+iiy].fPop[iPop];
          }

          if ( (iX+iix) > 1 && (iY+iiy) > 1){

            // 2nd collision on line x-1, y-1
            collideNode(&(sim->tmpLattice[iX+iix-1][iY+iiy-1]));

            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = iX+iix-1 + c[iPop][0];
              nextY = iY+iiy-1 + c[iPop][1];
              sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop[iPop];
            }

            if ( (iX+iix) > 2 && (iY+iiy) > 2){
              // 3rd collision on line x-2, y-2
              collideNode(&(sim->lattice[iX+iix-2][iY+iiy-2]));

              for (iPop = 0; iPop < 9; ++iPop) {
                nextX = iX+iix-2 + c[iPop][0];
                nextY = iY+iiy-2 + c[iPop][1];

                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                  sim->lattice[iX+iix-2][iY+iiy-2].fPop[iPop];
              }
            }
          }
        }
      }
    }// end of iY loop
  }// end of iX loop

  //-----------------------------------------------------------------------
  // Column iX=[1~lx-1], y=ly need to compute the second collision & stream
  iY = ly;
  for(iX = 1; iX < sim->lx; ++iX){
    collideNode(&(sim->tmpLattice[iX][iY]));

    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }

  // Row iY=[1~ly], iX=lx need to compute the second collision & stream
  // iX=sim->lx;
  //simple optimize
  iX = lx;
  for (iY = 1; iY <= sim->ly; ++iY){
    collideNode(&(sim->tmpLattice[iX][iY]));

    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }

  //-----------------------------------------------------------------------
  // Column iX=[1~lx-2], y=ly-1 need to compute the third collision & stream
  iY = ly-1;
  for (iX = 1; iX < (sim->lx-1); ++iX){
    collideNode(&(sim->lattice[iX][iY]));

    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  // Row iY=[1~ly-1], iX=lx-1 need to compute the third collision & stream
  iX = lx-1;
  for (iY = 1; iY < sim->ly; ++iY){
    collideNode(&(sim->lattice[iX][iY]));

    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  //-----------------------------------------------------------------------
  // Column iX=[1~lx-1], y=ly need to compute the third collision&stream
  iY = ly;
  for(iX = 1; iX < sim->lx; ++iX){
    // collide_stream_buf1_to_buf2(sim, iX, iY);
    collideNode(&(sim->lattice[iX][iY]));

    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }

  // Row iY=1~ly, iX=lx need to compute the third collision&stream
  // iX=sim->lx;
  // simple optimize
  iX = lx;
  for (iY = 1; iY <= sim->ly; ++iY){
    collideNode(&(sim->lattice[iX][iY]));

    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY].fPop[iPop];
    }
  }
}
#endif

#if 1 // use unrolling and remove row_index %
void step3CollideStreamTileOMP(Simulation* sim) {
  int lx = sim->lx, ly = sim->ly;

#ifdef _OPENMP
#pragma omp parallel default(shared)
{
  int iX, iY, iPop, iix, iiy;
  int nextX, nextY;
  
  int tid = omp_get_thread_num();
  int my_lx[2];
  my_lx[0] = 1 + tid * my_domain_H;
  my_lx[1] = (tid + 1) * my_domain_H;

  // ------------------ step I : Preprocess intersection area ------------------------------//
  // 1st fused collision and streaming on last row of 1 thread block (0, 1, 2)
  if (my_lx[0] == 1) {
    iX = my_lx[0];
    for (iY = 1; iY <= ly; ++iY) {
      // my up + 1
      collideNode(&(sim->lattice[iX+1][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX+1 + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX+1][iY].fPop[iPop];
      }

      // mine
      collideNode(&(sim->lattice[iX][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }

      if (iY > 1){
        // 2nd fused c&s on my left & up point (x+1, y-1)
        collideNode(&(sim->tmpLattice[iX][iY-1]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX + c[iPop][0];
          nextY = iY-1 + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX][iY-1].fPop[iPop];
        }
      }
    }

    iX = lx;
    for (iY = 1; iY <= ly; ++iY) {
      // my down - 1
      collideNode(&(sim->lattice[iX-1][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX-1 + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX-1][iY].fPop[iPop];
      }

      // mine
      collideNode(&(sim->lattice[iX][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }

      if (iY > 1){
        // 2nd fused c&s on my left (x,y-1) on row x=lx-1
        collideNode(&(sim->tmpLattice[iX][iY-1]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX + c[iPop][0];
          nextY = iY-1 + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX][iY-1].fPop[iPop];
        }
      }
    }
    
  }
  else {
    iX = my_lx[0] - 2;
    for (iY = 1; iY <= ly; ++iY) {
      for (int i = 0; i < 4; ++i) {
        collideNode(&(sim->lattice[iX + i][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX + i + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+i][iY].fPop[iPop];
        }
      }

      if (iY > 1){
        // 2nd fused c&s on my left & up point (x+1, y-1)
        collideNode(&(sim->tmpLattice[iX+1][iY-1]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY-1 + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX+1][iY-1].fPop[iPop];
        }

        // 2nd c&s on my left (x, y-1)
        collideNode(&(sim->tmpLattice[iX+2][iY-1]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+2 + c[iPop][0];
          nextY = iY-1 + c[iPop][1];
          sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[iX+2][iY-1].fPop[iPop];
        }
      }
    }
  }
  #pragma omp barrier
  // ------------------ End step1 ------------------------------//

  // ------------------ step II: Main Area Computation --------------------//
  // Outer loops.
  for (int outerX = my_lx[0] + 2; outerX <= my_lx[1]; outerX += tile) {
    for (int outerY = 1; outerY <= ly ; outerY += tile) {
      // Inner loops.
      int innerX_max = MIN(outerX + tile - 1, my_lx[1]);
      for (int innerX = outerX; innerX <= innerX_max; ++innerX) {

        int innerY_max = MIN(outerY + tile - 1, ly);

        if (innerX == my_lx[1] || innerX == (my_lx[1] - 1) ) {
          for (int innerY = outerY; innerY <= innerY_max; ++innerY) {
            if (innerY > 1){           
              // 2nd fused c&s
              collideNode(&(sim->tmpLattice[innerX - 1][innerY - 1]));
              for (iPop = 0; iPop < 9; ++iPop) {
                nextX = innerX - 1 + c[iPop][0];
                nextY = innerY - 1 + c[iPop][1];
                sim->lattice[nextX][nextY].fPop[iPop] =
                  sim->tmpLattice[innerX - 1][innerY - 1].fPop[iPop];
              }

              if (innerY > 2){
                // 3rd fused c&s
                collideNode(&(sim->lattice[innerX - 2][innerY - 2]));
                for (iPop = 0; iPop < 9; ++iPop) {
                  nextX = innerX - 2 + c[iPop][0];
                  nextY = innerY - 2 + c[iPop][1];

                  sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[innerX - 2][innerY - 2].fPop[iPop];
                }
              }
            }
          }

          continue;
        }
        
        for (int innerY = outerY; innerY <= innerY_max; ++innerY) {
          // 1st fused c&s
          collideNode(&(sim->lattice[innerX][innerY]));
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = innerX + c[iPop][0];
            nextY = innerY + c[iPop][1];

            sim->tmpLattice[nextX][nextY].fPop[iPop] =
              sim->lattice[innerX][innerY].fPop[iPop];
          }
          
          if (innerY > 1){
            // 2nd fused c&s
            collideNode(&(sim->tmpLattice[innerX - 1][innerY - 1]));
            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = innerX - 1 + c[iPop][0];
              nextY = innerY - 1 + c[iPop][1];
              sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[innerX - 1][innerY - 1].fPop[iPop];
            }

            if (innerY > 2){
              // 3rd fused collision and streaming
              collideNode(&(sim->lattice[innerX - 2][innerY - 2]));
              for (iPop = 0; iPop < 9; ++iPop) {
                nextX = innerX - 2 + c[iPop][0];
                nextY = innerY - 2 + c[iPop][1];

                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                  sim->lattice[innerX - 2][innerY - 2].fPop[iPop];
              }
            }
          }
        } // end innerY
      } // end innerX
    } // end outerY
  } // end outerX
  #pragma omp barrier
  // ------------------ End step2 ------------------------------//

  // ------------------ step III: handle intersection area --------------------//
  // Each thread computes the 2nd c&s on the rightmost column iY = ly, row# = [1..my_domain_H]
  #pragma omp for schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX){
    // 2nd fused collision and streaming
    iY = ly;
    collideNode(&(sim->tmpLattice[iX][iY]));
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY + c[iPop][1];
      sim->lattice[nextX][nextY].fPop[iPop] =
        sim->tmpLattice[iX][iY].fPop[iPop];
    }
  }
  // This must has a barrierr here

  // The rest are within each data domain
  // Each thread computes the 3rd c&s on the upper row iX = my_domain_H,  iY = [1, ly - 2] 
  #pragma omp for nowait schedule(static)
  for (iX = my_domain_H; iX <= lx; iX += my_domain_H){
    for (iY = 1; iY <= ly-2; iY++){
      // 3rd fused collision and streaming
      collideNode(&(sim->lattice[iX-1][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX-1 + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX-1][iY].fPop[iPop];
      }

      // 3rd fused collision and streaming
      collideNode(&(sim->lattice[iX][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX][iY].fPop[iPop];
      }
    } 
  }

  #pragma omp for nowait schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX){
    // 3rd fused collision and streaming on rightmost - 1 column
    iY = ly;
    collideNode(&(sim->lattice[iX][iY-1]));
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY-1 + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY-1].fPop[iPop];
    }

    // 3rd fused collision and streaming on rightmost column
    collideNode(&(sim->lattice[iX][iY]));
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
#endif

#if 0 // use inline and become slower
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
#endif

#if 0
void step3CollideStreamTileOMP(Simulation* sim) {
  int iX, iY, iPop, iix, iiy;
  int lx = sim->lx, ly = sim->ly;
  int nextX, nextY;
  int row_index;
  //compute each thread upper boundary line at iX=my_domain_H 1st c+s

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
        collideNode(&(sim->lattice[iX+1][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+1][iY].fPop[iPop];
        }

        // my up + 2
        collideNode(&(sim->lattice[iX+2][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+2 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+2][iY].fPop[iPop];
        }

        if (iY > 1){
          // 2nd fused c&s on my left & up point (x+1, y-1)
          collideNode(&(sim->tmpLattice[iX+1][iY-1]));
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
        collideNode(&(sim->lattice[iX-1][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX-1][iY].fPop[iPop];
        }

        // mine
        collideNode(&(sim->lattice[iX][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX][iY].fPop[iPop];
        }

        if (iY > 1){
          // 2nd fused c&s on my left (x,y-1) on row x=lx

          collideNode(&(sim->tmpLattice[iX][iY-1]));
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
        collideNode(&(sim->lattice[iX-1][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX-1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX-1][iY].fPop[iPop];
        }

        // mine
        collideNode(&(sim->lattice[iX][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX][iY].fPop[iPop];
        }

        // my up + 1
        collideNode(&(sim->lattice[iX+1][iY]));
        for (iPop = 0; iPop < 9; ++iPop) {
          nextX = iX+1 + c[iPop][0];
          nextY = iY + c[iPop][1];

          sim->tmpLattice[nextX][nextY].fPop[iPop] =
            sim->lattice[iX+1][iY].fPop[iPop];
        }

        // my up + 2
        collideNode(&(sim->lattice[iX+2][iY]));
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
          for (iPop = 0; iPop < 9; ++iPop) {
            nextX = iX+1 + c[iPop][0];
            nextY = iY-1 + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
              sim->tmpLattice[iX+1][iY-1].fPop[iPop];
          }

          // 2nd c&s on my left (x, y-1)
          collideNode(&(sim->tmpLattice[iX][iY-1]));
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
  #pragma omp for private(iix, iiy, iX, iY, row_index, iPop, nextX, nextY) schedule(static, my_domain_H/tile)
  for (iX = 1; iX <= lx; iX += tile) {
    for (iY = 1; iY <= ly; iY += tile) {
      for (iix = 0; iix < tile; ++iix){
        row_index =  (iX+iix) % my_domain_H;
        if ((row_index == 1) || (row_index == 2)) continue;

        // from row 3, 4, ..., my_domain_H-1, 0
        for (iiy = 0; iiy < tile; ++iiy) {

          if (row_index != 0 && row_index != (my_domain_H-1) ){
            // 1st fused c&s
            collideNode(&(sim->lattice[iX+iix][iY+iiy]));
            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = iX+iix + c[iPop][0];
              nextY = iY+iiy + c[iPop][1];

              sim->tmpLattice[nextX][nextY].fPop[iPop] =
                sim->lattice[iX+iix][iY+iiy].fPop[iPop];
            }
          }
          if ((iY+iiy) > 1){

            // 2nd fused c&s
            collideNode(&(sim->tmpLattice[iX+iix-1][iY+iiy-1]));
            for (iPop = 0; iPop < 9; ++iPop) {
              nextX = iX+iix-1 + c[iPop][0];
              nextY = iY+iiy-1 + c[iPop][1];
              sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop[iPop];
            }

            if ((iY+iiy) > 2){
              // 3rd fused collision and streaming
              collideNode(&(sim->lattice[iX+iix-2][iY+iiy-2]));
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
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, my_domain_H)
  for (iX = 1; iX <= lx; ++iX){

    // 2nd fused collision and streaming
    iY = ly;
    collideNode(&(sim->tmpLattice[iX][iY]));
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
      collideNode(&(sim->lattice[iX-1][iY]));
      for (iPop = 0; iPop < 9; ++iPop) {
        nextX = iX-1 + c[iPop][0];
        nextY = iY + c[iPop][1];

        sim->tmpLattice[nextX][nextY].fPop[iPop] =
          sim->lattice[iX-1][iY].fPop[iPop];
      }

      collideNode(&(sim->lattice[iX][iY]));
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
    collideNode(&(sim->lattice[iX][iY-1]));
    for (iPop = 0; iPop < 9; ++iPop) {
      nextX = iX + c[iPop][0];
      nextY = iY-1 + c[iPop][1];

      sim->tmpLattice[nextX][nextY].fPop[iPop] =
        sim->lattice[iX][iY-1].fPop[iPop];
    }

    collideNode(&(sim->lattice[iX][iY]));
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
#endif