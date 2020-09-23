#!/bin/bash
#SBATCH --partition=RM
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=04:00:00

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

    if (( $OMP_NUM_THREADS == 0 )); then
      continue
    fi

    for repeat in 0 1 2 3 4; do
      echo "spread $mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${TILE:-1}"
      $mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${TILE:-1}
      echo "---------------------------------------------------------------------"
    echo
    done
  done
}

# threads=(1 2 4 8 14 16 28)
threads=(0 0 0 0 14 16 28)

if [ $DIM == 112 ]; then
  warmup_steps=(900 900 900 900 900 900 900)
  steps=(1800 1800 1800 1800 1800 1800 1800)

elif [ $DIM == 224 ]; then
  warmup_steps=(900 900 900 900 900 900 900)
  steps=(1800 1800 1800 1800 1800 1800 1800)

elif [ $DIM == 448 ]; then
  warmup_steps=(900 900 900 900 900 900 900)
  steps=(1800 1800 1800 1800 1800 1800 1800)

elif [ $DIM == 896 ]; then
  warmup_steps=(900 900 900 900 900 900 900)
  steps=(900 900 900 900 900 900 900)

elif [ $DIM == 1792 ]; then
  warmup_steps=(768 768 768 768 768 768 768)
  steps=(600 600 600 600 600 600 600)

elif [ $DIM == 3584 ]; then
  warmup_steps=(48 48 96 96 192 192 192)
  steps=(96 96 192 192 384 384 384)

  if [[ ${CODE} == *"_tile"* ]]; then
    echo "contains _tile"
    warmup_steps=(48 48 96 96 768 768 768)
    steps=(96 96 192 192 600 600 600)
  fi

elif [ $DIM == 7168 ]; then
  warmup_steps=(24 24 48 48 96 96 96)
  steps=(48 48 96 96 192 192 192)

  if [[ ${CODE} == *"_tile"* ]]; then
    echo "contains _tile"
    warmup_steps=(24 24 48 48 768 768 768)
    steps=(48 48 96 96 600 600 600)
  fi

elif [ $DIM == 14336 ]; then
  warmup_steps=(6 6 12 24 48 48 48)
  steps=(6 6 12 48 96 96 96)

  if [[ ${CODE} == *"_tile"* ]]; then
    echo "contains _tile"
    warmup_steps=(6 6 12 24 768 768 768)
    steps=(6 6 12 48 600 600 600)
  fi

elif [ $DIM == 20720 ]; then
  warmup_steps=(6 6 12 24 48 48 48)
  steps=(6 6 12 48 96 96 96)

  if [[ ${CODE} == *"_tile"* ]]; then
    echo "contains _tile"
    warmup_steps=(6 6 12 24 192 192 192)
    steps=(6 6 12 48 384 384 384)
  fi

elif [ $DIM == 27104 ]; then
  warmup_steps=(6 6 12 24 48 48 48)
  steps=(6 6 12 48 96 96 96)

  if [[ ${CODE} == *"_tile"* ]]; then
    echo "contains _tile"
    warmup_steps=(6 6 12 24 96 96 96)
    steps=(6 6 12 48 192 192 192)
  fi

else
  echo "Wrong parameters"
fi

time omp_square > omp_dim_${DIM}_${CODE}${TILE:+_}${TILE}.log