#test
echo "origin"
./origin/unsteady 400 400 0 120 20
echo "fuse"
./fuse/unsteady 400 400 0 120 20
echo "2step"
./2step/unsteady 400 400 0 120 20
echo "2step_tile"
./2step_tile/unsteady 400 400 0 120 20
echo "3step"
./3step/unsteady 400 400 0 120 20
echo "3step_tile"
./3step_tile/unsteady 400 400 0 120 20

export OMP_NUM_THREADS=2
echo "origin_omp"
./origin_omp/unsteady 400 400 0 120 20
echo "fuse_omp"
./fuse_omp/unsteady 400 400 0 120 20
echo "2step_omp"
./2step_omp/unsteady 400 400 0 120 20
echo "2step_tile_omp"
./2step_tile_omp/unsteady 400 400 0 120 20
echo "3step_omp"
./3step_omp/unsteady 400 400 0 120 20
echo "3step_tile_omp"
./3step_tile_omp/unsteady 400 400 0 120 20