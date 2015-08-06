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
    parser.add_argument("mirna_path")
    parser.add_argument("full_transcriptome_path")
    parser.add_argument("outpath")
    args = parser.parse_args()
    
    compare_enriched_deseq_data = compareEnrichedDeseqData(
        args.mirna_path, args.full_transcriptome_path, args.outpath)
    compare_enriched_deseq_data.compare_regulated_gene()
    
class compareEnrichedDeseqData(object):
    '''Filter the RPS blast result'''
    def __init__(self, mirna_path,
                 full_transcriptome_path,
                 outpath):
        self._mirna_path = mirna_path
        self._full_transcriptome_path = full_transcriptome_path
        self._outpath = outpath
        
    def compare_regulated_gene(self):
        ''''''
        for files in os.listdir(self._mirna_path):
            outfile = self._outpath+'/common_'+files.split(
                '.')[0]+'.csv'
            header = "Libraries\tmiRNA:BM-A\tmiRNA:BM-B\tmiRNA:Log2FC\tTotal:BM-A\tTotal:BM-B\tTotal:Log2FC" 
            gene_dict1 = self._compile_mir_data(self._mirna_path+'/'+files)
            gene_dict2 = self._compile_mir_data(self._full_transcriptome_path+'/'+files)
            self._create_outfile(outfile, header,
                                 gene_dict1, gene_dict2)
    
    def _compile_mir_data(self, infile):
        ''''''
        data_dict = {}
        with open(infile, 'r') as in_fh:
            for entry in in_fh:
                mirna = entry.split('\t')[0]
                bm_a = entry.split()[-6]
                bm_b = entry.split()[-5]
                log2fc = entry.split()[-3]
                data_dict[mirna] = "%s\t%s\t%s" % (bm_a, bm_b, log2fc
                )
        return data_dict
        
    def _create_outfile(self, outfile, header,
                        gene_dict1, gene_dict2):
        ''''''
        with open(outfile, 'w') as out_fh:
            out_fh.write(header)
            for genes in set(gene_dict1.keys()
                             ).intersection(set(gene_dict2)):
                out_fh.write('%s\t%s\t%s\n' % (
                    genes, gene_dict1[genes], gene_dict2[genes]))
        
        
if __name__ == '__main__':
    main()