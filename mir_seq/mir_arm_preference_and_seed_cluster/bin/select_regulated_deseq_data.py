#!/usr/bin/env python 

'''Filter all the rna related data.'''

#'''FUNCTION & USAGE'''

import os
import argparse
from collections import defaultdict

__description__ = ""
__author__ = "Malvika Sharan <malvika.sharan@uni-wuerzburg.de>"
__email__ = "malvika.sharan@uni-wuerzburg.de"
__version__ = ""

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("infile")
    parser.add_argument("threshold")
    parser.add_argument("outpath")
    args = parser.parse_args()
    
    select_enriched_deseq_data = SelectEnrichedDeseqData(
        args.infile, args.threshold, args.outpath)
    select_enriched_deseq_data.select_regulated_gene()
    select_enriched_deseq_data.create_outfile()
    
class SelectEnrichedDeseqData(object):
    '''Filter the RPS blast result'''
    def __init__(self, infile,
                 threshold,
                 outpath):
        self._infile = infile
        self._threshold = threshold
        self._outpath = outpath
        
        self._header = ''
        self._selected_entries = set()
        
    def select_regulated_gene(self):
        ''''''
        with open(self._infile, 'r') as in_fh:
            for entry in in_fh:
                if 'Libraries' in entry:
                    self._header = entry
                else:
                    read1 = entry.split('\t')[1]
                    if float(read1) >= 10:
                        log2fc = entry.split('\t')[-3]
                        if float(log2fc) >= float(self._threshold) or float(log2fc) <= float(self._threshold)*-1:
                            self._selected_entries.add(
                                entry.strip())
        return self._header, self._selected_entries
    
    def create_outfile(self):
        ''''''
        filename = self._infile.split('/')[-1].split(
            'deseq_')[1].split('_extended')[0]
        if 'old' in filename:
            filename = filename.replace('old_', '')
        with open(self._outpath+'/'+filename+'.csv', 'w') as out_fh:
            out_fh.write(self._header)
            out_fh.write(('\n').join(list(self._selected_entries)))
        
        
if __name__ == '__main__':
    main()
