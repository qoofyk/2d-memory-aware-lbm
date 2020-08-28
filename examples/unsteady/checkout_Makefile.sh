#!/usr/bin/env bash

# First created: 
# Last modified: 2020 Aug 27

# Author: Yuankun Fu
# email: qoofyk@gmail.com

for CODE in origin fuse fuse_tile 2step 3step 2step_tile 3step_tile; do
  git checkout ${CODE}/Makefile
done

for CODE in origin_omp fuse_omp fuse_tile_omp 2step_omp 3step_omp 2step_tile_omp 3step_tile_omp; do
  git checkout ${CODE}/Makefile
done