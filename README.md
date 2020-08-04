# 2D memory aware LBM
Repo for paper **Yuankun Fu**, Feng Li, Fengguang Song, Luoding Zhu, **Designing a Parallel Memory-Aware Lattice Boltzmann Algorithm on Manycore Systems**, *Proceedings of 30th International Symposium on Computer Architecture and High Performance Computing* (**SBAC-PAD'18**), Lyon, France, September 2018. [[PDF](https://ieeexplore.ieee.org/abstract/document/8645909)]

nearest neighbors: D2Q9

# Implementation

We have `origin-(omp)`, `fuse-(omp)`, `two-steps-(tile)-(omp)`, `three steps-(tile)-(omp)` versions

# Build instruction

download scons
```
tar xf scons-3.1.1.tar.gz
cd scons-3.1.1
mkdir build
python3 setup.py install --prefix=./build/

source ~/intel/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 
```

# Run Benchmark

The input parameters are `Height`, `Width`, `warmUpIter`, `numIter`, `tile_size`

For OMP code, need `export OMP_NUMTHREADS=8`
```
cd 2d-memory-aware-lbm/examples/unsteady/origin
make
./unsteady 100 100 100 100 0
```

Compile & Test
```
cd 2d-memory-aware-lbm/examples/unsteady
sh compile.sh
sh test.sh
```