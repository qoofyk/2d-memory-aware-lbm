# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 19:03:33 2020

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
    ('cores', (re.compile(r'^.*NUM_THREADS=([0-9]+)$', re.MULTILINE), int)),
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
    def __init__(self, machine, dim, tile):
        self.dim = dim
        self.tile = tile
        self.machine = machine

        self.header = ['origin', 'fuse', 'fuse tile', '2step', '2step tile', '3step', '3step tile']
        self.table = collections.defaultdict(lambda: collections.defaultdict(lambda: float('-inf')))
        self.max_mflups_tile = collections.defaultdict(lambda: collections.defaultdict(lambda: 1))
    
    def filter(self, row):
        return row['dim'] == self.dim
    
    def error_value(self):
        raise Exception('error_value() must be customized by the subclass')
    
    def parse_filename(self, filename):
        #splittext give filename[0] + extension [1]
        fields = os.path.splitext(os.path.basename(filename))[0].split('_')
        
        try:
            dim_idx = fields.index('dim')
        except ValueError:
            dim_idx = None
        
        try:
            self.tile = int(fields[-1])
        except ValueError:
            self.tile = None
        
        return {
            'exp_kind': fields[0],
            'dim': dim_idx and int(fields[dim_idx+1]),
            'name': ' '.join(fields[dim_idx+2:]),
            'tile': self.tile,
        }
    
    def analyze(self, filename):
        # print("analysis", filename)
        
        compute = collections.OrderedDict([
            ('mflups_per_core', lambda t: t['mflups'] / t['cores']),
        ])
    
        # Parse input columns:
        with open(filename) as f:
            text = f.read()
        table = collections.OrderedDict(
            #group(1) is the first value of group
            (k, numpy.asarray([t(m.group(1)) for m in re.finditer(p, text)])) # a key-value pair (k, numpy_array). p='mflups', text=a re
            for k, (p, t) in _columns.items()) #outermost loop to iterate _columns
        
        # print(table)
    
        assert table['cores'].size > 0, "logs are empty"
        
        # table['mflups'] = numpy.delete(table['mflups'], table['mflups'].argmin())
        
        # Group by cores count and compute statistics:
        # multi variable assignment
        table['cores'], table['mflups'], table['std'], table['reps'], table['dims'], table['iterations'], table['warmup'], table['tile']= list(map( #list()) will combine all data into list
            numpy.asarray,
            zip(*[(k, numpy.mean(mflups), numpy.std(mflups), len(mflups), dims[0], iteration[0], warmup[0], self.tile) #this is a tuple
                  for k, vs in group_by(table['cores'], zip(table['mflups'], table['dims'], table['iterations'], table['warmup']))
                  for mflups, dims, iteration, warmup in [zip(*vs)]])))
    
        for column in ('dims', 'iterations', 'warmup'):
            # print(table['cores'].shape)
            table[column] = numpy.resize(table[column], table['cores'].shape)
    
        # Compute derived columns:
        for k, f in compute.items():
            table[k] = f(table)
        
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
        # print(*list(data.values()))
        for values in zip(*list(data.values())):
            # print(values)
            items = dict(zip(data.keys(), values))
            # print(items)
            name = row['name']
            name = name[0 : name.index('omp') - 1]

            # if name not in self.header:
            #     self.header.append(name)

            # self.table[items['cores']][name] = max(
            #     items['mflups'],
            #     self.table[items['cores']][name],
            #     key=float)
        
            if (float(items['mflups']) > float(self.table[items['cores']][name])):
                self.table[items['cores']][name] = items['mflups']
                self.max_mflups_tile[items['cores']][name] = items['tile']
                
    
    def complete(self):
        # FIXME: This isn't actually the criteria we'd like to sort on,
        # we'd prefer to sort so that the list of names roughly parallels
        # the order of the bars in the graph.

        self.header.insert(0, 'cores')
        # print(self.header)

        out = csv.DictWriter(sys.stdout, self.header, dialect="excel")
        out.writeheader()
        
        # print(self.table.keys())
        for c in sorted(self.table.keys()):
            row = self.table[c]
            # print(row)
            row = {k: None if v == float('-inf') else v for k, v in row.items()}
            row['cores'] = c #add a key-value pair
            out.writerow(row)
        
        with open("omp_square_dim_" + str(self.dim) + "_best_tile.csv", 'w') as f:
            out = csv.DictWriter(f, self.header, dialect="excel")
            out.writeheader()
            for c in sorted(self.max_mflups_tile.keys()):
                row = self.max_mflups_tile[c]
                # print(row)
                row = {k : v for k, v in row.items()}
                row['cores'] = c #add a key-value pair
                out.writerow(row)

    def parse(self, dim, tile):
        has_exception = False
        log_filenames = glob.glob(self.machine + '/omp_square/*.log', recursive=False)
        # print(log_filenames)
        
        for filename in log_filenames:
            row = self.parse_filename(filename)
            # print(len(row), row['name'])
            if not self.filter(row): #filter a specific dim of logs
                continue
            
            try:
                data = self.analyze(filename)
                # print(len(data))
                # print(*data['mflups'])
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
            
def driver(machine, dim, tile):
    parser = Parser(machine, dim, tile)
    parser.parse(dim, tile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--machine', required=True)
    parser.add_argument('-d', '--dim',  type=int, required=True)
    parser.add_argument('-t', '--tile', type=int, default=0)

    args = parser.parse_args()
    driver(**vars(args)) #all parameters combined into dict (hash table)
