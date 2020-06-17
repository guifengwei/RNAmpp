#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 12 Jun 2020 14:30:15

import sys, os, argparse, string
from collections import Counter 

def main():
    ''' main scripts '''
    m6A2anno={}
    for line in open(sys.argv[1], "r"):
        if not line.startswith("#"):
            line = line.strip().split("\t")
            if line[3] not in m6A2anno:
                m6A2anno[line[3]]=[ line[6]+"\t"+line[7] ]
            else:
                m6A2anno[line[3]].append(  line[6]+"\t"+line[7]  )
    for k,v in m6A2anno.items():
        if len(v) == 1:
            print("{0}\t{1}".format(k,v[0]))
        else:
            if "protein_coding\texon" in v:
                print ("{0}\t{1}".format(k, "protein_coding\texon") )
            elif "ncRNA\texon" in v:
                print ("{0}\t{1}".format(k, "ncRNA\texon") )
            elif "protein_coding\tintron" in v:
                print ("{0}\t{1}".format(k, "protein_coding\tintron") )
            elif "ncRNA\tintron" in v:
                print ("{0}\t{1}".format(k, "ncRNA\tintron") )
            else:
                print ("{0}\t{1}".format(k, "NA\tNA") )



if __name__ == "__main__":
    main()

