#!/usr/bin/env bash

# First created: 2020 Aug 27
# Last modified: 2020 Aug 27

# Author: Yuankun Fu
# email: qoofyk@gmail.com

rm -rf */*.dat

# This is correct results
cd origin
./unsteady 24 24 30 0 4
cd ..

# Now test and verify with different other code
for CODE in fuse fuse_tile 2step 3step 2step_tile 3step_tile; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady 24 24 30 0 4
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_24.dat vel_${CODE}_zgb_24.dat
  cd ..
done

# Verify OMP
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=1
for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady 24 24 30 0 4
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_24.dat vel_${CODE}_zgb_24.dat
  cd ..
done

export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=2
for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady 24 24 30 0 4
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_24.dat vel_${CODE}_zgb_24.dat
  cd ..
done

export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=4
for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady 24 24 30 0 4
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_24.dat vel_${CODE}_zgb_24.dat
  cd ..
done

rm -rf */*.dat