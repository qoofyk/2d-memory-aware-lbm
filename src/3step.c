#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "assert.h"

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

void step3CollideStreamOMP(Simulation* sim) {
  unsigned int iX, iY, iPop;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;
  int row_index;

// #define DEBUG_PRINT

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: total_values)
{
  #ifdef ADDPAPI
    long long local_values[NUM_EVENT];
    long long local_extra_values[NUM_EVENT];
    int event_codes[NUM_EVENT];
    int EventSet = PAPI_NULL;
    int retval;

    int EventCode;
    for (int i = 0; i < NUM_EVENT; ++i){
      /* Convert to integer */
      if (PAPI_event_name_to_code(native_name[i] , &EventCode) != PAPI_OK){
        PAPI_perror(errstring);
        ERROR_RETURN(retval);
      }

      event_codes[i] = EventCode;
    }

    /* Creating event set   */
    if ((retval = PAPI_create_eventset(&EventSet)) != PAPI_OK)
      ERROR_RETURN(retval);

    /* Add the array of events PAPI_TOT_INS and PAPI_TOT_CYC to the eventset*/
    if ((retval = PAPI_add_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK){
      PAPI_perror(errstring);
      ERROR_RETURN(retval);
    }

    /* Reset the counting events in the Event Set */
    if (PAPI_reset(EventSet) != PAPI_OK){
      fprintf(stderr, "PAPI_reset error!\n");
      exit(1);
    }

    /* Start counting */
    if ( (retval=PAPI_start(EventSet)) != PAPI_OK)
      ERROR_RETURN(retval);
  #endif

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
  #pragma omp for private(iX, iY, row_index, iPop, nextX, nextY) schedule(static, thread_block)
  for (iX = 1; iX <= lx; ++iX) {
    row_index =  iX % thread_block;
    if ((row_index == 1) || (row_index == 2)) continue;
    // from row 3, 4, ..., thread_block-1, 0
    for (iY = 1; iY <= ly; ++iY) {

      if (row_index != 0 && row_index != (thread_block-1) ){
        // 1st fused c&s
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
      if (iY > 1){
        // 2nd fused c&s
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

        if (iY > 2){
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

  #ifdef ADDPAPI
    /* read the counter values and store them in the values array */
    if ( (retval=PAPI_read(EventSet, local_values)) != PAPI_OK)
      ERROR_RETURN(retval);

    /* Stop counting, this reads from the counter as well as stop it. */
    if ( (retval=PAPI_stop(EventSet, local_extra_values)) != PAPI_OK)
      ERROR_RETURN(retval);

    if ( (retval=PAPI_remove_events(EventSet, event_codes, NUM_EVENT)) != PAPI_OK)
    ERROR_RETURN(retval);

    /* Free all memory and data structures, EventSet must be empty. */
    if ( (retval=PAPI_destroy_eventset(&EventSet)) != PAPI_OK)
      ERROR_RETURN(retval);

    #ifdef _OPENMP
      int my_rank = omp_get_thread_num();
      int thread_count = omp_get_num_threads();
    #else
      int my_rank = 0;
      int thread_count = 1;
    #endif

    for(int i = 0; i < NUM_EVENT; ++i){
      // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
      // fflush(stdout);
      total_values[i] += local_values[i];
    }
  #endif
}
#else
    printf("No OPENMP used");
#endif

}// end of func