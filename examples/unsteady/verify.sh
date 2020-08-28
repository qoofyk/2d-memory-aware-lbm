#!/usr/bin/env bash

# First created: 2020 Aug 27
# Last modified: 2020 Aug 27

# Author: Yuankun Fu
# email: qoofyk@gmail.com

rm -rf */*.dat

Height=24
Width=24
Warm_up=60
Measure=0
Tile=4
Verify_timepoint=54

# This is correct results
cd origin
./unsteady $Height $Width $Warm_up $Measure $Tile
cd ..

# Now test and verify with different other code
for CODE in fuse fuse_tile 2step 3step 2step_tile 3step_tile; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady $Height $Width $Warm_up $Measure $Tile
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_${Verify_timepoint}.dat vel_${CODE}_zgb_${Verify_timepoint}.dat
  cd ..
  echo "---------------------------------------"
done

# Verify OMP
echo "Verify with 1 thread"
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=1
for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady $Height $Width $Warm_up $Measure $Tile
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_${Verify_timepoint}.dat vel_${CODE}_zgb_${Verify_timepoint}.dat
  cd ..
  echo "---------------------------------------"
done

echo "Verify with 2 thread"
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=2
for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady $Height $Width $Warm_up $Measure $Tile
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_${Verify_timepoint}.dat vel_${CODE}_zgb_${Verify_timepoint}.dat
  cd ..
  echo "---------------------------------------"
done

echo "Verify with 4 threads"
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=4
for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  echo "run ${CODE}"
  cd ${CODE}
  ./unsteady $Height $Width $Warm_up $Measure $Tile
  echo "verify ${CODE}"
  diff ../origin/vel_origin_zgb_${Verify_timepoint}.dat vel_${CODE}_zgb_${Verify_timepoint}.dat
  cd ..
  echo "---------------------------------------"
done

rm -rf */*.dat