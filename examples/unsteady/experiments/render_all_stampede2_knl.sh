#!/usr/bin/env bash

# First created: 2020 Aug 17
# Last modified: 2020 Aug 17

# Author: Yuankun Fu
# email: qoofyk@gmail.com

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo $root_dir

python3 ${root_dir}/seq_mflups.py -m "stampede2_knl" > stampede2_knl_seq_square.csv
python3 ${root_dir}/render_mflups.py ./stampede2_knl_seq_square.csv --xdata "dims" --xlabel "Side Length of a square lattice"  --ylim "(0,15)" --no-ylog
${root_dir}/pdfcrop stampede2_knl_seq_square.pdf stampede2_knl_seq_square.pdf

# --no-xlog
# for DIM in 272; do
#   python3 ${root_dir}/omp_mflups.py -m "stampede2_knl" -d ${DIM} > stampede2_knl_omp_square_dim_${DIM}.csv
#   python3 ${root_dir}/render_mflups.py ./stampede2_knl_omp_square_dim_${DIM}.csv --width 15 --ylim "(0,400)" --no-xlog --no-ylog --xlabel-size 8
#   ${root_dir}/pdfcrop stampede2_knl_omp_square_dim_${DIM}.pdf stampede2_knl_omp_square_dim_${DIM}.pdf
# done

# for DIM in 544; do
#   python3 ${root_dir}/omp_mflups.py -m "stampede2_knl" -d ${DIM} > stampede2_knl_omp_square_dim_${DIM}.csv
#   python3 ${root_dir}/render_mflups.py ./stampede2_knl_omp_square_dim_${DIM}.csv --width 15 --ylim "(0,500)" --no-xlog --no-ylog --xlabel-size 8
#   ${root_dir}/pdfcrop stampede2_knl_omp_square_dim_${DIM}.pdf stampede2_knl_omp_square_dim_${DIM}.pdf
# done

# for DIM in 1088 2176 4352 8704 17408; do
#   python3 ${root_dir}/omp_mflups.py -m "stampede2_knl" -d ${DIM} > stampede2_knl_omp_square_dim_${DIM}.csv
#   python3 ${root_dir}/render_mflups.py ./stampede2_knl_omp_square_dim_${DIM}.csv --width 15 --ylim "(0,600)" --no-xlog --no-ylog --xlabel-size 8
#   ${root_dir}/pdfcrop stampede2_knl_omp_square_dim_${DIM}.pdf stampede2_knl_omp_square_dim_${DIM}.pdf
# done

# --xlog
for DIM in 272; do
  python3 ${root_dir}/omp_mflups.py -m "stampede2_knl" -d ${DIM} > stampede2_knl_omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./stampede2_knl_omp_square_dim_${DIM}.csv --ylim "(0,400)" --no-xlog --no-ylog --xlabel-size 8.75
  ${root_dir}/pdfcrop stampede2_knl_omp_square_dim_${DIM}.pdf stampede2_knl_omp_square_dim_${DIM}.pdf
done

for DIM in 544; do
  python3 ${root_dir}/omp_mflups.py -m "stampede2_knl" -d ${DIM} > stampede2_knl_omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./stampede2_knl_omp_square_dim_${DIM}.csv --ylim "(0,500)" --no-xlog --no-ylog --xlabel-size 8.75
  ${root_dir}/pdfcrop stampede2_knl_omp_square_dim_${DIM}.pdf stampede2_knl_omp_square_dim_${DIM}.pdf
done

for DIM in 1088 2176 4352 8704 17408 20400 21760; do
  python3 ${root_dir}/omp_mflups.py -m "stampede2_knl" -d ${DIM} > stampede2_knl_omp_square_dim_${DIM}.csv
  python3 ${root_dir}/render_mflups.py ./stampede2_knl_omp_square_dim_${DIM}.csv --ylim "(0,600)" --no-xlog --no-ylog --xlabel-size 8.75
  ${root_dir}/pdfcrop stampede2_knl_omp_square_dim_${DIM}.pdf stampede2_knl_omp_square_dim_${DIM}.pdf
done