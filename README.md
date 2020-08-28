# 2D memory aware LBM
Repo for paper **Yuankun Fu**, Feng Li, Fengguang Song, Luoding Zhu, **Designing a Parallel Memory-Aware Lattice Boltzmann Algorithm on Manycore Systems**, *Proceedings of 30th International Symposium on Computer Architecture and High Performance Computing* (**SBAC-PAD'18**), Lyon, France, September 2018. [[PDF](https://ieeexplore.ieee.org/abstract/document/8645909)]

nearest neighbors: D2Q9

# Implementation

We have `origin-(omp)`, `fuse-(omp)`, `two-steps-(tile)-(omp)`, `three steps-(tile)-(omp)` versions

# Build instruction

download scons
```bash
tar xf scons-3.1.1.tar.gz
cd scons-3.1.1
mkdir build
python3 setup.py install --prefix=./build/

source ~/intel/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 
```

# Verification of Correctness

Uncomment the `#define SAVE` in `src/lb.h` and each algorithm will generate output for every 1 step, 2 steps or 3steps during machine warm up stage.

The input parameters are `Height`, `Width`, `Warm_Up_Iter`, `Measure_Iter`, `Tile_size`.

You can verify with larger `Height` or `Width` or longer `warmUpIter`.

Use the `vel_${Algorithm}_zgb_#` data to verify, since this is the starting scalar velocity of each grid point before calling different algorithm function.
Always Use the data in `origin` as the true base data, then compare other algorithm's output with `origin`.

E.g., you can manually test with `origin`, `2step` and `3step`.
```bash
cd 2d-memory-aware-lbm/examples/unsteady/origin or 2step or 3step or 3step_tile
make

2d-memory-aware-lbm/examples/unsteady/2step$./unsteady 24 24 30 0 4
2d-memory-aware-lbm/examples/unsteady/3step$./unsteady 24 24 30 0 4
2d-memory-aware-lbm/examples/unsteady/origin$./unsteady 24 24 30 0 4

$diff ../origin/vel_origin_zgb_12.dat vel_step3-line_zgb_12.dat
$diff ../origin/vel_origin_zgb_12.dat vel_step2-line_zgb_12.dat
diff ../origin/vel_origin_zgb_24.dat vel_step3-tile_zgb_24.dat

# Batch to verify for sequential and OMP code using 1, 2, 4 threads
cd 2d-memory-aware-lbm/examples/unsteady/
sh verify.sh > v.log
```

# Run Benchmark

The input parameters are `Height`, `Width`, `warmUpIter`, `numIter`, `tile_size`
For OMP code, need `export OMP_NUMTHREADS=8`
```bash
cd 2d-memory-aware-lbm/examples/unsteady/origin
make
./unsteady 100 100 100 100 0
```

Compile & Test
```bash
cd 2d-memory-aware-lbm/examples/unsteady
sh compile.sh
sh test.sh
```

# Experiments
## Bridges
```bash
cd 2d-memory-aware-lbm/examples/unsteady/experiments/bridges/omp_square
DIM=14336 CODE=3step_tile_omp TILE=64 bash -c 'sbatch --nodes=1 --job-name="omp-square-dim-${DIM}-${CODE}${TILE:+_}${TILE}" mflups.sh'

Or
sh submit_all.sh
```

# Post operation on job output
```bash
less seq_tile=16.txt | grep 'Mega' | cut -d ":" -f 2 | cut -d " " -f 2
```

# Collect Data and generate graph
make sure install python3 and matplotlib
```bash
module load python3
pip install matplotlib
module load texlive

cd /home/qoofyk/2d-memory-aware-lbm/examples/unsteady/experiments
sh render_all_bridges.sh
```