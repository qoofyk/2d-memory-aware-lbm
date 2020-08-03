# 2D memory aware LBM

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
