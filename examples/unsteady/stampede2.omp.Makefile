##########################################################################
## Makefile.
##
## The present Makefile is a pure configuration file, in which
## you can select compilation options. Compilation dependencies
## are managed automatically through the Python library SConstruct.
##
## If you don't have Python, or if compilation doesn't work for other
## reasons, consult the Palabos user's guide for instructions on manual
## compilation.
##########################################################################

# USE: multiple arguments are separated by spaces.
#   For example: projectFiles = file1.cpp file2.cpp
#                optimFlags   = -O -finline-functions

# Leading directory of the Palabos source code
palabosRoot  = ../../..
# Name of source files in current directory to compile and link with Palabos
projectFiles = unsteady.c

# Set optimization flags on/off
optimize     = true
# Set debug mode and debug flags on/off
debug        = false
# Set profiling flags on/off
profile      = false
# Set MPI-parallel mode on/off (parallelism in cluster-like environment)
MPIparallel  = false
# Set SMP-parallel mode on/off (shared-memory parallelism)
SMPparallel  = true
# Decide whether to include calls to the POSIX API. On non-POSIX systems,
#   including Windows, this flag must be false, unless a POSIX environment is
#   emulated (such as with Cygwin).
usePOSIX     = true

# Path to external source files (other than Palabos)
srcPaths =
# Path to external libraries (other than Palabos)
libraryPaths =
# Path to external libraries (other than Palabos)
includePaths = -I/opt/intel/vtune_amplifier_2018.4.0.574913/include
# Dynamic and static libraries (other than Palabos)
libraries    = /opt/intel/vtune_amplifier_2018.4.0.574913/lib64/libittnotify.a

# Compiler to use without MPI parallelism
serialCXX    = icc
# Compiler to use with MPI parallelism
parallelCXX  = mpiicc
# General compiler flags (e.g. -Wall to turn on all warnings on g++) 
compileFlags = -Wall -Wextra -Wnon-virtual-dtor -Wno-deprecated-declarations -DPLB_MAC_OS_X -std=gnu99 -D_OPENMP -qopenmp -qopt-report5 -qopt-report-phase=vec -qopenmp-simd ${TACC_VEC_FLAGS} # -fno-strict-aliasing -no-fma

# General linker flags (don't put library includes into this flag)
linkFlags    = -lm -qopenmp
# Compiler flags to use when optimization mode is on
optimFlags   = -O3
#optimFlags   = -xHOST -O3 -ip -no-prec-div -static
# Compiler flags to use when debug mode is on
debugFlags   = -g -gdwarf-2 #-parallel-source-info=2 -debug inline-debug-info
# Compiler flags to use when profile mode is on
profileFlags = -pg


##########################################################################
# All code below this line is just about forwarding the options
# to SConstruct. It is recommended not to modify anything there.
##########################################################################

SCons     = $(palabosRoot)/scons/scons-3.1.1/build/bin/scons -j 16 -f $(palabosRoot)/SConstruct

SConsArgs = palabosRoot=$(palabosRoot) \
            projectFiles="$(projectFiles)" \
            optimize=$(optimize) \
            debug=$(debug) \
            profile=$(profile) \
            MPIparallel=$(MPIparallel) \
            SMPparallel=$(SMPparallel) \
            usePOSIX=$(usePOSIX) \
            serialCXX=$(serialCXX) \
            parallelCXX=$(parallelCXX) \
            compileFlags="$(compileFlags)" \
            linkFlags="$(linkFlags)" \
            optimFlags="$(optimFlags)" \
            debugFlags="$(debugFlags)" \
            profileFlags="$(profileFlags)" \
            srcPaths="$(srcPaths)" \
            libraryPaths="$(libraryPaths)" \
            includePaths="$(includePaths)" \
            libraries="$(libraries)"

compile:
	python $(SCons) $(SConsArgs)

clean:
	python $(SCons) -c $(SConsArgs)
	/bin/rm -vf `find $(palabosRoot) -name '*~'`
