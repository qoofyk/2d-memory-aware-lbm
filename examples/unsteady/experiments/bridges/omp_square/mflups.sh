#!/bin/bash
#SBATCH --partition=RM
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL

export OMP_PROC_BIND=spread
total_cores=$(( $(echo $SLURM_JOB_CPUS_PER_NODE | cut -d'(' -f 1) ))
cores=$(( $total_cores))

mybin=../../../${CODE}/unsteady

omp_square () {
  # module list
  echo $total_cores #bridges give 28
  echo $cores

  # Start Simulation
  Height=$DIM
  Width=$DIM
    # for ((k=0; k<${#threads[@]}; k++)); do
    for ((k=${#threads[@]}-1; k>=0; k--)); do
      export OMP_NUM_THREADS=${threads[k]}

      for repeat in 0 1 2 3 4; do
        echo "spread $mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${TILE:-1}"
        $mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${TILE:-1}
        echo "---------------------------------------------------------------------"
      echo
      done
    done
}

threads=(1 2 4 8 14 16 28)

if [ $DIM == 112 ]; then
  warmup_steps=(600 600 600 600 600 600 600)
  steps=(600 600 600 600 600 600 600)

elif [ $DIM == 224 ]; then
  warmup_steps=(600 600 600 600 600 600 600)
  steps=(600 600 600 600 600 600 600)

elif [ $DIM == 448 ]; then
  warmup_steps=(600 600 600 600 600 600 600)
  steps=(600 600 600 600 600 600 600)

elif [ $DIM == 896 ]; then
  warmup_steps=(900 900 900 900 900 900 900)
  steps=(600 600 600 600 600 600 600)

elif [ $DIM == 1792 ]; then
  warmup_steps=(768 768 768 768 768 768 768)
  steps=(600 600 600 600 600 600 600)

elif [ $DIM == 3584 ]; then
  warmup_steps=(48 48 96 96 192 192 192)
  steps=(96 96 192 192 384 384 384)

elif [ $DIM == 7168 ]; then
  warmup_steps=(24 24 48 48 96 96 96)
  steps=(48 48 96 96 192 192 192)

elif [ $DIM == 14336 ]; then
  warmup_steps=(6 6 12 24 24 24 48)
  steps=(6 6 12 48 48 48 96)

else
  echo "Wrong parameters"
fi

time omp_square > omp_dim_${DIM}_${CODE}${TILE:+_}${TILE}.log