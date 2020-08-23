#!/usr/bin/env bash

# First created: 2020 Aug 09
# Last modified: 2020 Aug 09

# Author: Yuankun Fu
# email: qoofyk@gmail.com

CODE=origin bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}" --output="jobtime/seq-square-${CODE}.out" mflups.sh'
CODE=fuse   bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}" --output="jobtime/seq-square-${CODE}.out" mflups.sh'
CODE=2step  bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}" --output="jobtime/seq-square-${CODE}.out" mflups.sh'
CODE=3step  bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}" --output="jobtime/seq-square-${CODE}.out" mflups.sh'

for t in 8 16 32 64 128 256; do
  CODE=fuse_tile  TILE=$t bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}${TILE:+_}${TILE}" --output="jobtime/seq-square-${CODE}${TILE:+_}${TILE}.out" mflups.sh'
  CODE=2step_tile TILE=$t bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}${TILE:+_}${TILE}" --output="jobtime/seq-square-${CODE}${TILE:+_}${TILE}.out" mflups.sh'
  CODE=3step_tile TILE=$t bash -c 'sbatch --nodes=1 --job-name="seq-square-${CODE}${TILE:+_}${TILE}" --output="jobtime/seq-square-${CODE}${TILE:+_}${TILE}.out" mflups.sh'
done