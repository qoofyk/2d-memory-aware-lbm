#!/usr/bin/env bash

# First created: 2020 Aug 09
# Last modified: 2020 Aug 09

# Author: Yuankun Fu
# email: qoofyk@gmail.com

# 112 224 448 896 1792 3584 7168 14336
for d in 1792 3584 7168 14336; do
  DIM=$d CODE=origin_omp bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
  DIM=$d CODE=fuse_omp   bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
  DIM=$d CODE=2step_omp  bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
  DIM=$d CODE=3step_omp  bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
done

for t in 8 16 32 64 128 256; do
  for d in 1792 3584 7168 14336; do
    DIM=$d CODE=fuse_tile_omp  TILE=$t bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
    DIM=$d CODE=2step_tile_omp TILE=$t bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
    DIM=$d CODE=3step_tile_omp TILE=$t bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" --output="jobtime/omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
  done
done

# DIM=112 CODE=fuse_omp bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'
# DIM=7168 CODE=3step_omp_tile TILE=64 bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'