#!/usr/bin/python
# Programmer: Wei Guifeng, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 11 Jun 2020 22:24:20

import sys, os, argparse, string, re

'''
subset of GTF, this Program will take less than 3 minutes for all the analysis

'''

def main():
    ''' main scripts '''
    gene_pattern = re.compile('gene_id \"([A-Za-z0-9\-\.\(\)]+)\"')
    trans_pattern= re.compile('transcript_id \"([A-Za-z0-9\-\.\(\)]+)\"')
    Transcript2Gene = {}
    geneid=[]
    for line in open(sys.argv[1], 'r'):
        geneid.append(line.strip())
    for line in open(sys.argv[2], 'r'):
        if line.startswith('#'): continue
        line = line.strip()
        #print line
        gene= re.findall(gene_pattern, line)[0]
        transcript = re.findall(trans_pattern, line)[0]
        if gene in geneid:
            print(line)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        print('Usage: This_script.py genelist *.gtf', file=sys.stderr)

