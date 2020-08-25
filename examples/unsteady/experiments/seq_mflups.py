# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 02:06:58 2020

@author: Yuankun Fu
"""

import argparse
import collections
import csv
import numpy
import os
import re
import sys
import glob
import traceback

_columns = collections.OrderedDict([
    ('mflups', (re.compile(r'^.*iterations: ([0-9.]+) Mega.*$', re.MULTILINE), float)),
    ('dims', (re.compile(r'^lx=([0-9]+).*$', re.MULTILINE), int)),
    ('iterations', (re.compile(r'^.*numIter=([0-9]+).*$', re.MULTILINE), int)),
    ('warmup', (re.compile(r'^.*warmUpIter=([0-9]+).*$', re.MULTILINE), int)),
])

def group_by(keys, values):
    last_key = None
    last_group = None
    for key, value in zip(keys, values):
        if key != last_key:
            if last_group is not None:
                yield (last_key, last_group)
            last_group = []
        last_key = key
        last_group.append(value)

    if last_group is not None:
        yield (last_key, last_group)

class Parser:
    def __init__(self, machine, tile):
        self.tile = tile
        self.machine = machine

        self.header = ['origin', 'fuse', 'fuse tile', '2step', '2step tile', '3step', '3step tile']
        self.table = collections.defaultdict(lambda: collections.defaultdict(lambda: float('-inf')))
        self.max_mflups_tile = collections.defaultdict(lambda: collections.defaultdict(lambda: 1))
    
    def error_value(self):
        raise Exception('error_value() must be customized by the subclass')
    
    def parse_filename(self, filename):
        #splittext give filename[0] + extension [1]
        fields = os.path.splitext(os.path.basename(filename))[0].split('_')
        
        try:
            self.tile = int(fields[-1])
        except ValueError:
            self.tile = None
        
        return {
            'exp_kind': fields[0],
            'name': ' '.join(fields[1:]),
            'tile': self.tile,
        }
    
    def analyze(self, filename):
        # print("analysis", filename)

        # compute = collections.OrderedDict([
        #     ('mflups_per_core', lambda t: t['mflups'] / t['cores']),
        # ])

        # Parse input columns:
        with open(filename) as f:
            text = f.read()
        table = collections.OrderedDict(
            (k, numpy.asarray([t(m.group(1)) for m in re.finditer(p, text)]))
            for k, (p, t) in _columns.items())
        
        # print(table)
        # table['mflups'] = numpy.delete(table['mflups'], table['mflups'].argmin()

        assert table['dims'].size > 0, "logs are empty"
        # Group by cores count and compute statistics:
        table['dims'], table['mflups'], table['std'], table['reps'], table['iterations'], table['warmup'], table['tile'] = list(map(
            numpy.asarray,
            zip(*[(k, numpy.mean(mflups), numpy.std(mflups), len(mflups), iteration[0], warmup[0], self.tile)
                  for k, vs in group_by(table['dims'], zip(table['mflups'], table['iterations'], table['warmup']))
                  for mflups, iteration, warmup in [zip(*vs)]])))

        for column in ('iterations', 'warmup'):
            table[column] = numpy.resize(table[column], table['dims'].shape)

        # Compute derived columns:
        # for k, f in compute.items():
        #     table[k] = f(table)
        
        for outlier in table['std']:
            if (outlier > 1):
                print(filename, "has a potential outlier", outlier, file=open("omp_square_outlier.txt", "a"))

        # Post-process table for output:
        result = collections.OrderedDict()
        for k, c in table.items():
            if any(isinstance(x, float) for x in c):
                result[k] = ['%f' % x for x in c]
            else:
                result[k] = table[k]

        # print(result)

        out_filename = os.path.splitext(filename)[0] + '.csv'
        with open(out_filename, 'w') as f:
            out = csv.writer(f)
            out.writerow(result.keys())
            out.writerows(zip(*list(result.values())))
        
        return (result)
    
    def process(self, row, data, mflups=None):
        for values in zip(*list(data.values())):
            items = dict(zip(data.keys(), values))
            name = row['name']
            if 'tile' in name:
                name = name[0 : name.index('tile') + 4]

            # if name not in self.header:
            #     self.header.append(name)

            # self.table[items['cores']][name] = max(
            #     items['mflups'],
            #     self.table[items['cores']][name],
            #     key=float)
        
            if (float(items['mflups']) > float(self.table[items['dims']][name])):
                self.table[items['dims']][name] = items['mflups']
                self.max_mflups_tile[items['dims']][name] = items['tile']
                
    
    def complete(self):
        # FIXME: This isn't actually the criteria we'd like to sort on,
        # we'd prefer to sort so that the list of names roughly parallels
        # the order of the bars in the graph.
        # self.header.sort()
        self.header.insert(0, 'dims')
        # print(self.header)

        out = csv.DictWriter(sys.stdout, self.header, dialect="excel")
        out.writeheader()
        
        # print(self.table.keys())
        for d in sorted(self.table.keys()):
            row = self.table[d]
            # print(row)
            row = {k: None if v == float('-inf') else v for k, v in row.items()}
            row['dims'] = d #add a key-value pair
            out.writerow(row)
        
        with open(self.machine + 'seq_square_best_tile.csv', 'w') as f:
            out = csv.DictWriter(f, self.header, dialect="excel")
            out.writeheader()
            for d in sorted(self.max_mflups_tile.keys()):
                row = self.max_mflups_tile[d]
                # print(row)
                row = {k : v for k, v in row.items()}
                row['dims'] = d #add a key-value pair
                out.writerow(row)

    def parse(self, tile):
        has_exception = False
        log_filenames = glob.glob(self.machine + '/seq_square/*.log', recursive=False)
        # print(log_filenames)
        
        for filename in log_filenames:
            row = self.parse_filename(filename)

            try:
                data = self.analyze(filename)
            except Exception as e:
                print("get exception")
                
                print('%s: %s: %s' % (filename, type(e).__name__, e), file=sys.stderr)
                data = (self.error_value(),)
                has_exception = True
            self.process(row, data)

        self.complete()
        
        if has_exception:
            print('Errors were encountered while parsing. Run with -v to see full error messages.', file=sys.stderr)
            print(file=sys.stderr)
            
def driver(machine, tile):
    parser = Parser(machine, tile)
    parser.parse(tile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--machine', required=True)
    # parser.add_argument('-e', '--exp', dest='exp_kind', required=True, help='specify experiment type')
    parser.add_argument('-t', '--tile', type=int, default=0)

    args = parser.parse_args()
    driver(**vars(args))
