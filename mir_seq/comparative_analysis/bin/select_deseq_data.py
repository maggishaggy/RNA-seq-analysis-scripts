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
    parser.add_argument("inpath")
    parser.add_argument("threshold")
    parser.add_argument("outpath")
    args = parser.parse_args()
    
    select_enriched_deseq_data = SelectEnrichedDeseqData(
        args.inpath, args.threshold, args.outpath)
    select_enriched_deseq_data.select_regulated_gene()
    
class SelectEnrichedDeseqData(object):
    def __init__(self, inpath,
                 threshold,
                 outpath):
        self._inpath = inpath
        self._threshold = threshold
        self._outpath = outpath
        
    def select_regulated_gene(self):
        ''''''
        for files in os.listdir(self._inpath):
            print(self._inpath, files)
            selected_entries = set()
            header = ''
            outfile = self._outpath+'/'+files.split('.')[0]+'.csv'
            with open(self._inpath+'/'+files, 'r') as in_fh:
                for entry in in_fh:
                    if 'log2FoldChange' in entry:
                        header = entry
                    else:
                        read1 = entry.split('\t')[1]
                        read2 = entry.split('\t')[-10]
                        if float(read1) >= 10 or float(read2) >= 10:
                            log2fc = entry.split('\t')[-3]
                            if not log2fc == 'NA':
                                if log2fc == '#NAME?':
                                    log2fc = 8
                                if float(log2fc) >= float(
                                    self._threshold) or float(
                                    log2fc) <= float(self._threshold)*-1:
                                    selected_entries.add(entry.strip())
            self._create_outfile(outfile, header, selected_entries)
    
    def _create_outfile(self, outfile, header, selected_entries):
        ''''''
        with open(outfile, 'w') as out_fh:
            out_fh.write(header)
            out_fh.write(('\n').join(list(selected_entries)))
        
        
if __name__ == '__main__':
    main()