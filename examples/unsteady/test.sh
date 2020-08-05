#test
echo "origin"
./origin/unsteady 400 400 60 60 20
echo
echo "fuse"
./fuse/unsteady 400 400 60 60 20
echo
echo "fuse_tile"
./fuse_tile/unsteady 400 400 60 60 20
echo
echo "2step"
./2step/unsteady 400 400 60 60 20
echo
echo "2step_tile"
./2step_tile/unsteady 400 400 60 60 20
echo
echo "3step"
./3step/unsteady 400 400 60 60 20
echo
echo "3step_tile"
./3step_tile/unsteady 400 400 60 60 20
echo

export OMP_NUM_THREADS=2
echo "origin_omp"
./origin_omp/unsteady 400 400 60 60 20
echo
echo "fuse_omp"
./fuse_omp/unsteady 400 400 60 60 20
echo
echo "fuse_tile_omp"
./fuse_tile_omp/unsteady 400 400 60 60 20
echo
echo "2step_omp"
./2step_omp/unsteady 400 400 60 60 20
echo
echo "2step_tile_omp"
./2step_tile_omp/unsteady 400 400 60 60 20
echo
echo "3step_omp"
./3step_omp/unsteady 400 400 60 60 20
echo
echo "3step_tile_omp"
./3step_tile_omp/unsteady 400 400 60 60 20
echo