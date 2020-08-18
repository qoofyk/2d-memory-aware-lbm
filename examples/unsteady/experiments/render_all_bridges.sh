#!/usr/bin/env bash

# First created: 2020 Aug 17
# Last modified: 2020 Aug 17

# Author: Yuankun Fu
# email: qoofyk@gmail.com

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo $root_dir

if [[ $1 = crop ]]; then
    function crop {
        ./pdfcrop "$1" >/dev/null && mv "$(basename "$1" .pdf)"-crop.pdf "$1"
    }
else
    function crop { true; }
fi

python3 ${root_dir}/seq_mflups.py -m "bridges" > bridges_seq_square.csv
python3 ${root_dir}/render_mflups.py ./bridges_seq_square.csv --xdata "dims" --xlabel "Length of a square grid"  --ylim "(0,35)" --no-ylog
crop bridges_seq_square.pdf

for DIM in 112 224 448 896 1792 3584 7168 14336; do
  python3 ${root_dir}/omp_mflups.py -m "bridges" -d ${DIM} > omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./omp_square_dim_${DIM}.csv --ylim "(0,600)" --no-ylog
  crop omp_square_dim_${DIM}.pdf
done