#!/usr/bin/env python

import os
import sys
import argparse

__description__ = ""
__author__ = "Malvika Sharan <malvika.sharan@uni-wuerzburg.de>"
__email__ = "malvika.sharan@uni-wuerzburg.de"
__version__ = ""

def main():
    '''arguments'''
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("in_file")
    parser.add_argument("out_file")
    args = parser.parse_args()

    create_mir_gff = CreateMirGff(args.in_file,
                             args.out_file)
    create_mir_gff.parse_in_data()
    create_mir_gff.create_out_files()
    
class CreateMirGff(object):
    def __init__(self, in_file,
                 out_file):
        self._in_file = in_file
        self._out_file = out_file
        
        self._mir_seed = {}
        self._stemloop_list = []
        
    def parse_in_data(self):
        ''''''
        with open(self._in_file, 'r') as in_fh:
            for entry in in_fh.read().split('>')[1:]:
                #print(entry)
                mir_id = entry.split('\n')[0].split(' ')[0]
                primary = entry.split('\n')[0].split(' ')[1]
                seq = entry.strip().split('\n')[1]
                if 'hsa' in mir_id:
                    seed = ''.join(list(seq)[1:7])
                    self._mir_seed.setdefault(seed, []).append("%s\t%s\t%s\t%s" % (seed, primary, mir_id, seq))
        return self._mir_seed

    def create_out_files(self):
        with open(self._out_file, 'w') as out_fh:
            for seeds in self._mir_seed.keys():
                out_fh.write('\n'.join(self._mir_seed[seeds])+'\n//\n')

if __name__ == '__main__':
    main()
