#!/usr/bin/env bash

# First created: 
# Last modified: 2020 Aug 18

# Author: Yuankun Fu
# email: qoofyk@gmail.com
cd origin
yes | cp -f ../stampede2.seq.Makefile Makefile
cd ../fuse
yes | cp -f ../stampede2.seq.Makefile Makefile
cd ../fuse_tile
yes | cp -f ../stampede2.seq.Makefile Makefile
cd ../2step
yes | cp -f ../stampede2.seq.Makefile Makefile
cd ../2step_tile
yes | cp -f ../stampede2.seq.Makefile Makefile
cd ../3step
yes | cp -f ../stampede2.seq.Makefile Makefile
cd ../3step_tile
yes | cp -f ../stampede2.seq.Makefile Makefile

cd ../origin_omp
yes | cp -f ../stampede2.omp.Makefile Makefile
cd ../fuse_omp
yes | cp -f ../stampede2.omp.Makefile Makefile
cd ../fuse_tile_omp
yes | cp -f ../stampede2.omp.Makefile Makefile
cd ../2step_omp
yes | cp -f ../stampede2.omp.Makefile Makefile
cd ../2step_tile_omp
yes | cp -f ../stampede2.omp.Makefile Makefile
cd ../3step_omp
yes | cp -f ../stampede2.omp.Makefile Makefile
cd ../3step_tile_omp
yes | cp -f ../stampede2.omp.Makefile Makefile