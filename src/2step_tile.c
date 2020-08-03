#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "D2Q9.h"
#include "lb.h"
#include "assert.h"

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

void step2CollideStreamTileOMP(Simulation* sim) {
  unsigned int iX, iY, iPop;
  // unsigned int tile = 32;
  int iix, iiy;
  int lx=sim->lx, ly=sim->ly;
  double ux1, uy1, ux2, uy2;
  int nextX, nextY;

  //compute each thread upper boundary line at iX=thread_block 1st c+s
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

  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = thread_block; iX <= lx; iX+=thread_block){
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
  #pragma omp for private(iix, iiy, iX, iY, iPop, nextX, nextY) schedule(static, thread_block/tile)
  for (iX = 1; iX <= lx; iX+=tile) {
    for (iY = 1; iY <= ly; iY+=tile) {
      for (iix = 0; iix < tile; ++iix){
        for (iiy = 0; iiy < tile; ++iiy){

          if ( (iX+iix) % thread_block != 0){

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
            if ( (iX+iix-1) % thread_block != 0){

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

  //compute thread boundary line at iX=thread_block 2nd c+s
  //NOTICE: 1~ly-1 !!! use tmpLattice !!!
  #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
  for (iX = thread_block; iX < lx; iX += thread_block){
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
  #pragma omp for private(iX, iPop, nextX, nextY) schedule(static, thread_block)
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
  #pragma omp for private(iY, iPop, nextX, nextY) schedule(static, thread_block)
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
}