# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 19:30:57 2020

@author: qoofyk
"""

# ./seq_square.csv --xdata "dims" --xlabel "Length of a square grid"  --ylim "(0,35)" --no-ylog
# ./omp_square_dim_14336.csv --ylim "(0,570)" --no-ylog
# ./omp_square/omp_dim_14336_fuse_tile_omp_16.csv --xbase 10 --no-xticks --ylabel "Efficiency" --ylim "(0,32)" --no-ylog --y-percent --highlight-column "metg"

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import argparse
import ast
import csv
import math
import numpy as np
import os

def csv2rec(filename):
    # print(filename)
    return np.recfromtxt(filename, dtype=None, delimiter=',', names=True)

def str_csv2rec(filename):
    # print(filename)
    return np.genfromtxt(filename, dtype=str, delimiter=',')

parser = argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('--title')
parser.add_argument('--width', type=float, default=9)
parser.add_argument('--height', type=float, default=5)
parser.add_argument('--legend', default='./legend.csv')
parser.add_argument('--legend-ncol', type=int, default=1)
parser.add_argument('--legend-fontsize', type=int, default=12)
parser.add_argument('--legend-position', default='best')
parser.add_argument('--legend-base', type=int, default=0)
parser.add_argument('--legend-suffix', action='append', default=[])
parser.add_argument('--limit-suffix')
parser.add_argument('--limit-intersection-filename')
parser.add_argument('--limit-intersection-system')
parser.add_argument('--filter-legend-even-powers', action='store_true')
parser.add_argument('--xlabel', default='Number of cores')
parser.add_argument('--ylabel', default='MFLUPS')
parser.add_argument('--xdata', default='cores')
parser.add_argument('--xscale', type=float, default=0)
parser.add_argument('--yscale', type=float, default=0)
parser.add_argument('--xlim', type=ast.literal_eval, default='None')
parser.add_argument('--ylim', type=ast.literal_eval, default='None')
parser.add_argument('--x-invert', action='store_true')
parser.add_argument('--y-invert', action='store_true')
parser.add_argument('--x-percent', action='store_true')
parser.add_argument('--y-percent', action='store_true')
parser.add_argument('--xbase', type=int, default=2)
parser.add_argument('--no-xlog', action='store_false', dest='xlog')
parser.add_argument('--no-ylog', action='store_false', dest='ylog')
parser.add_argument('--no-xticks', action='store_false', dest='xticks')
parser.add_argument('--connect-missing', action='store_true')
parser.add_argument('--highlight-column')
parser.add_argument('--ideal-column')
parser.add_argument('--xticks-rotate', dest='xticks_rotate', type=int, default=0)
parser.add_argument('--xlabel-size', dest='xlabel_size', type=float, default=12)
args = parser.parse_args()

markers = [
    '*',
    '^',
    'h',
    'D',
    'o',
    '<',
    's',
    '1',
    'p',
    'x',
    'D',
    '8',
    '4',
    'v',
    '+',
    'H',
    '2',
    '>'
] * 10

colors = [
    # Built in colors
    'gold', # no use
    (126 / 255.0, 0 / 255.0, 33 / 255.0), #'tab:red', # origin
    'orangered', # (255 / 255.0, 66 / 255.0, 14 / 255.0), # fuse_tile = google orange
    'C1', # fuse_tile = google yellow (255 / 255.0, 211 / 255.0, 32 / 255.0),
    'tab:green', # 2step
    'mediumorchid', #(49 / 255.0, 64 / 255.0, 4 / 255.0), # 2step_tile = google dark green
    'tab:blue', # 3step = google blue (0 / 255.0, 69 / 255.0, 134 / 255.0) / 'C0' / 
    'navy', #3step_tile
    
    # 'tab:brown',
    # 'mediumorchid',
    # 'orangered',
] * 10

# matplotlib.rcParams["font.family"] = "STIXGeneral"
# matplotlib.rcParams["mathtext.fontset"] = "stix"
matplotlib.rcParams["mathtext.fontset"] = "stixsans"

matplotlib.rc('xtick', labelsize=args.xlabel_size)
matplotlib.rc('ytick', labelsize=12)

fig = plt.figure(figsize=(args.width, args.height))
ax = fig.add_subplot(111)
plt.subplots_adjust(bottom=0.16, left=0.10)

ax.spines['top'].set_linewidth(1.0)
ax.spines['bottom'].set_linewidth(1.0)
ax.spines['left'].set_linewidth(1.0)
ax.spines['right'].set_linewidth(1.0)
ax.tick_params(axis='x', width=0.5)
ax.tick_params(axis='y', width=0.5)

if args.x_percent:
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: '{:.0%}'.format(x)))
if args.y_percent:
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda y, _: '{:.0%}'.format(y)))

if args.x_invert:
    ax.invert_xaxis()
if args.y_invert:
    ax.invert_yaxis()

if args.xlog and args.ylog:
    plt.loglog(basex=args.xbase)
elif args.xlog:
    plt.semilogx(basex=args.xbase)
elif args.ylog:
    plt.semilogy()

data = csv2rec(args.filename)
# print(data)
nodes = getattr(data, args.xdata)
if args.xscale:
    nodes = nodes * args.xscale
columns = list(data.dtype.names)
columns.remove(args.xdata)

if args.legend:
    # print(args.legend)
    legend_raw = str_csv2rec(args.legend)
    # print(legend_raw)
    legend_label = dict(zip(legend_raw[:, 0], legend_raw[:, 1]))
    legend_visible = dict(zip(legend_raw[:, 0], legend_raw[:, 2]))
    legend_idx = {}
    next_idx = 0
    for name in legend_raw[:, 0]:
        if legend_visible[name]:
            legend_idx[name] = next_idx
            next_idx += 1
        else:
            legend_idx[name] = None
else:
    next_idx = 0

for column in columns:
    visible = True
    limit = False
    if args.legend:
        colname = column
        colsuffix = ''
        for suffix in args.legend_suffix:
            if colname.endswith(suffix):
                colname = colname[:-len(suffix)]
                colsuffix = suffix.replace('_', ' ')
                if args.limit_suffix == suffix:
                    limit = True
                break
        if colname in legend_label:
            label = legend_label[colname] + colsuffix
            visible = legend_visible[colname]
            idx = legend_idx[colname]
        else:
            label = column.replace('_', ' ')
            idx = next_idx
            next_idx += 1
    else:
        label = column.replace('_', ' ')
        idx = next_idx
        next_idx += 1

    if args.legend_base > 0 and args.highlight_column != column:
        exponent = int(math.log(int(label), args.legend_base))
        if args.filter_legend_even_powers and exponent % 2 != 0:
            continue
        label = '$%s^{%.0f}$' % (args.legend_base, exponent)

    if not visible:
        continue

    if args.highlight_column == column:
        color = 'red'
        marker = None
        linetype = '--'
        linewidth = 3
        label = None
    elif args.ideal_column == column:
        color = 'black'
        marker = None
        linetype = ':'
        linewidth = 1
        label = None
    elif limit:
        color = colors[idx]
        marker = None
        linetype = '--'
        linewidth = 1
        label = '%s 50%%' % label
    else:
        color = colors[idx]
        marker = markers[idx]
        linetype = '-'
        linewidth = 1

    column_data = getattr(data, column)
    if args.yscale:
        column_data = column_data * args.yscale

    mask = np.isfinite(column_data)
    plt.plot(nodes[mask] if args.connect_missing else nodes, column_data[mask] if args.connect_missing else column_data, linestyle=linetype, color=color, marker=marker, markerfacecolor='none', linewidth=linewidth, label=label)

if args.limit_intersection_system:
    assert args.limit_intersection_filename is not None
    with open(args.limit_intersection_filename, newline='') as f:
        intersections = list(csv.DictReader(f, dialect='excel'))

    row = next(x for x in intersections if x['system'] == args.limit_intersection_system)

    assert not args.yscale, 'unimplemented'

    plt.plot(float(row['limit_actual_nodes']), float(row['limit_actual_time']), color='black', marker='s', markerfacecolor='none', markersize=12, markeredgewidth=3)
    plt.plot(float(row['limit_ideal_nodes']), float(row['limit_ideal_time']), color='black', marker='o', markerfacecolor='none', markersize=12, markeredgewidth=3)

if args.xticks:
    plt.xticks(nodes, nodes, rotation=args.xticks_rotate)
if args.xlim:
    plt.xlim(*args.xlim)
if args.ylim:
    plt.ylim(*args.ylim)

plt.xlabel(args.xlabel, fontsize=12)
plt.ylabel(args.ylabel, fontsize=12)
if args.title:
    plt.title(args.title, fontsize=14)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.80, box.height])

if args.legend:
    kwargs = {'loc': args.legend_position}
    if args.legend_position == 'center left':
        # Put a legend to the right of the current axis
        kwargs['bbox_to_anchor'] = (1, 0.5)
    plt.legend(
        ncol=args.legend_ncol, fontsize=args.legend_fontsize,
        # Square corners, disable transparency, set color to black
        fancybox=False, framealpha=1, edgecolor='black',
        **kwargs)

plt.grid(True, color='black', linestyle='--', linewidth=0.5, dashes=(1, 5))
output_filename = '%s.pdf' % os.path.splitext(args.filename)[0]
print('Generating %s' % output_filename)
plt.savefig(output_filename)
# plt.show()
plt.close()