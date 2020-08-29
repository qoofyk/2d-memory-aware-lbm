#!/bin/bash
#SBATCH --partition=RM
#SBATCH --qos=regular
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL

total_cores=$(( $(echo $SLURM_JOB_CPUS_PER_NODE | cut -d'(' -f 1) ))
cores=$(( $total_cores))

mybin=../../../${CODE}/unsteady

seq_square () {
  # module list
  echo $total_cores #bridges give 28
  echo $cores

  # Start Simulation
  for ((k=0; k<${#dim[@]}; k++)); do
  # for ((k=${#dim[@]}-1; k>=0; k--)); do
    Height=$(( 1 << dim[k] ))
    Width=$(( 1 << dim[k] ))
    
    if [[ $Height < ${TILE:-1} ]]
      continue
    fi

    for repeat in 0 1 2 3 4; do
      echo "$mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${TILE:-1}"
      $mybin $Height $Width ${warmup_steps[k]} ${steps[k]} ${TILE:-1}
      echo "---------------------------------------------------------------------"
    echo
    done
  done
}

#128 256 512 1024 2048 4096 8192 16384
dim=(7 8 9 10 11 12 13 14)
warmup_steps=(150 270 540 1050 210 30 12 6)
steps=(300 300 300 300 300 90 48 12)
seq_square > seq_${CODE}${TILE:+_}${TILE}.log