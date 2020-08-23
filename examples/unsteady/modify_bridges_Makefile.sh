#!/usr/bin/env bash

# First created: 
# Last modified: 2020 Aug 18

# Author: Yuankun Fu
# email: qoofyk@gmail.com
cd origin
patch -p1 < ../bridges.seq.Makefile.patch
cd ../fuse
patch -p1 < ../bridges.seq.Makefile.patch
cd ../fuse_tile
patch -p1 < ../bridges.seq.Makefile.patch
cd ../2step
patch -p1 < ../bridges.seq.Makefile.patch
cd ../2step_tile
patch -p1 < ../bridges.seq.Makefile.patch
cd ../3step
patch -p1 < ../bridges.seq.Makefile.patch
cd ../3step_tile
patch -p1 < ../bridges.seq.Makefile.patch

cd ../origin_omp
patch -p1 < ../bridges.omp.Makefile.patch
cd ../fuse_omp
patch -p1 < ../bridges.omp.Makefile.patch
cd ../fuse_tile_omp
patch -p1 < ../bridges.omp.Makefile.patch
cd ../2step_omp
patch -p1 < ../bridges.omp.Makefile.patch
cd ../2step_tile_omp
patch -p1 < ../bridges.omp.Makefile.patch
cd ../3step_omp
patch -p1 < ../bridges.omp.Makefile.patch
cd ../3step_tile_omp
patch -p1 < ../bridges.omp.Makefile.patch