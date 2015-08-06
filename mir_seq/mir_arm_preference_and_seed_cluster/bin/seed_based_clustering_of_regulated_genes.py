#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict

__description__ = ""
__author__ = "Malvika Sharan <malvika.sharan@uni-wuerzburg.de>"
__email__ = "malvika.sharan@uni-wuerzburg.de"
__version__ = ""

def main():
    '''arguments'''
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("reference_file")
    parser.add_argument("in_file")
    parser.add_argument("outpath")
    args = parser.parse_args()

    classify_mir_families = ClassifyMirFamilies(args.reference_file,
                            args.in_file, args.outpath)
    classify_mir_families.parse_ref_data()
    classify_mir_families.parse_in_data()
    classify_mir_families.create_out_file()
    
class ClassifyMirFamilies(object):
    def __init__(self, reference_file,
                 in_file,
                 outpath):
        self._reference_file = reference_file
        self._in_file = in_file
        self._outpath = outpath
        
        self._mir_family_data = {}
        self._in_data = {}
        self._family_stat = {}
        self._mir_ref_list = {}
        self._compiled_data = {}
        
    def parse_ref_data(self):
        with open(self._reference_file, 'r') as ref_fh:
            for entry_list in ref_fh.read().split('//'):
                for mirs in entry_list.split('\n'):
                    if not mirs == '':
                        #print(mirs)
                        seed = mirs.split('\t')[0]
                        primary = mirs.split('\t')[1]
                        mir_id = mirs.split('\t')[2]
                        self._mir_family_data[mir_id] = seed
        return self._mir_family_data

    def parse_in_data(self):
        with open(self._in_file, 'r') as in_fh:
            for entry in in_fh:
                if not 'Libraries' in entry.split('\t')[0]:
                    mir = entry.split('\t')[0].split(' - ')[0]
                    log_fold = entry.strip().split('\t')[-3]
                    if float(log_fold) >= 0.99 or float(log_fold) == 'Inf':
                        regulation_type = 'up_regulated'
                        self._in_data[mir] = "%s\t%s\t%s" % (
                                mir, log_fold, regulation_type)
                    elif float(log_fold) <= -0.99 or float(log_fold) == '-Inf':
                        regulation_type = 'down_regulated'
                        self._in_data[mir] = "%s\t%s\t%s" % (
                                mir, log_fold, regulation_type)
        return self._in_data

    def create_out_file(self):
        filename = self._in_file.split('/')[-1]
        for mir_ref in self._mir_family_data.keys():
            for mir_set in self._in_data.keys():
                if mir_set == mir_ref:
                    self._mir_ref_list.setdefault(
                        self._mir_family_data[mir_set], []).append(
                        self._in_data[mir_set])
                    self._compiled_data.setdefault(
                        self._mir_family_data[mir_set], []).append(self._in_data[mir_set])
        with open(self._outpath+'/seed_cluster_'+filename, 'w') as out_fh:
            for mir_ref_collected in self._mir_ref_list.keys():
                total_mirs = len(self._mir_ref_list[mir_ref_collected])
                for seed_based_entry in self._compiled_data[mir_ref_collected]:
                    out_fh.write("%s\t%s\t%s\n" % (mir_ref_collected,
                    total_mirs, seed_based_entry))

if __name__ == '__main__':
    main()