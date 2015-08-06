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
    parser.add_argument("outpath")
    args = parser.parse_args()
    
    regulated_mir_arm_preference = RegulatedMirArmPreference(
        args.infile, args.outpath)
    regulated_mir_arm_preference.regulated_gene_info()
    regulated_mir_arm_preference.create_out_file()
    
class RegulatedMirArmPreference(object):
    '''Filter the RPS blast result'''
    def __init__(self, infile,
                 outpath):
        self._infile = infile
        self._outpath = outpath
        
        self._arm_dict = {}
        
    def regulated_gene_info(self):
        ''''''
        with open(self._infile, 'r') as in_fh:
            for entry in in_fh:
                if not 'Libraries' in entry:
                    mir = entry.split('\t')[0].split(' - ')[0]
                    if '-3p' in mir:
                        self._arm_dict.setdefault('3p', []).append(mir)
                    elif '-5p' in mir:
                        self._arm_dict.setdefault('5p', []).append(mir)
                    else:
                        self._arm_dict.setdefault('UK', []).append(mir)
        return self._arm_dict
    
    def create_out_file(self):
        ''''''
        filename = self._infile.split('/')[-1]
        with open(self._outpath+'/arm_preference_'+filename, 'w') as out_fh:
            out_fh.write("Arm preference quantification:\n")
            out_fh.write("3p Arm preference: %i\n" % (int(len(self._arm_dict['3p']))))
            out_fh.write("5p Arm preference: %i\n" % (int(len(self._arm_dict['5p']))))
            out_fh.write("Undefined Arm: %i\n" % (int(len(self._arm_dict['UK']))))
            out_fh.write("\nmiRNAs with 3p Arm preference: %s\n" % ('\n'.join(self._arm_dict['3p'])))
            out_fh.write("\nmiRNAs with 5p Arm preference: %s\n" % ('\n'.join(self._arm_dict['5p'])))
            out_fh.write("\nmiRNAs with undefined Arms: %s\n" % ('\n'.join(self._arm_dict['UK'])))

if __name__ == '__main__':
    main()