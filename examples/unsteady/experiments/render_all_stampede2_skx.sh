#!/usr/bin/env bash

# First created: 2020 Aug 17
# Last modified: 2020 Aug 17

# Author: Yuankun Fu
# email: qoofyk@gmail.com

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo $root_dir

python3 ${root_dir}/seq_mflups.py -m "stampede2_skx" > stampede2_skx_seq_square.csv
python3 ${root_dir}/render_mflups.py ./stampede2_skx_seq_square.csv --xdata "dims" --xlabel "Side length of a square lattice"  --ylim "(0,65)" --no-ylog
${root_dir}/pdfcrop stampede2_skx_seq_square.pdf stampede2_skx_seq_square.pdf

for DIM in 192; do
  python3 ${root_dir}/omp_mflups.py -m "stampede2_skx" -d ${DIM} > stampede2_skx_omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./stampede2_skx_omp_square_dim_${DIM}.csv --ylim "(0,1200)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop stampede2_skx_omp_square_dim_${DIM}.pdf stampede2_skx_omp_square_dim_${DIM}.pdf
done

for DIM in 768; do
  python3 ${root_dir}/omp_mflups.py -m "stampede2_skx" -d ${DIM} > stampede2_skx_omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./stampede2_skx_omp_square_dim_${DIM}.csv --ylim "(0,1700)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop stampede2_skx_omp_square_dim_${DIM}.pdf stampede2_skx_omp_square_dim_${DIM}.pdf
done

for DIM in 384 1536 3072 6144 12288 24576 28800 33600; do
  python3 ${root_dir}/omp_mflups.py -m "stampede2_skx" -d ${DIM} > stampede2_skx_omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./stampede2_skx_omp_square_dim_${DIM}.csv --ylim "(0,1600)" --no-xlog --no-ylog
  ${root_dir}/pdfcrop stampede2_skx_omp_square_dim_${DIM}.pdf stampede2_skx_omp_square_dim_${DIM}.pdf
done