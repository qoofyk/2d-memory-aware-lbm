# sequential
cd origin
make clean;make
cd ../fuse
make
cd ../fuse_tile
make
cd ../2step
make
cd ../2step_tile
make
cd ../3step
make
cd ../3step_tile
make

# omp
cd ../origin_omp
make clean;make
cd ../fuse_omp
make
cd ../fuse_tile_omp
make
cd ../2step_omp
make
cd ../2step_tile_omp
make
cd ../3step_omp
make
cd ../3step_tile_omp
make