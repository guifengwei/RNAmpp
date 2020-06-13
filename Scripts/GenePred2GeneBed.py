#!/usr/bin/python
# Wei Guifeng, <guifengwei@gmail.com>
# -- coding: utf-8 --

import sys,string
#from itertools import izip

''' convert genepred format to standard genebed format '''

class genebed():
    '''genebed '''
    ''' 0-5: chrom, chromstart, chromend, name, score, strand,\
        6-11: cdsStart, cdsEnd, itemRGB, exonCount, exonSizes, exonStarts '''
    def __init__(self,x):
        try:
            self.bin=int(x[0])
            if(self.bin<10000):
                pass
                #x=x[1:]
        except:
            pass
        self.chr = x[0].rstrip()
        self.start = int(x[1])
        self.end = int(x[2])
        self.id  = x[3]
        self.score = x[4]
        self.strand = x[5]
        self.cdsStart = int(x[6])
        self.cdsEnd = int(x[7])
        self.itemRGB= x[8]
        self.exonCount = int(x[9])
        self.exonSizes = x[10].strip(",").split(",")
        self.exonSizes = list(map(int, self.exonSizes))
        # coordinates relative to the trnascropt start(start)  self.exonStarts
        self.exonStarts = x[11].strip(",").split(",")
        self.exonStarts = list(map(int, self.exonStarts))
        # coordinates in the chromosome
        # self.exon_Starts self.exon_Stops
        self.exon_Starts = [each + self.start for each in self.exonStarts]
        self.exon_Stops=[]
        for k,v in zip(self.exonStarts, self.exonSizes):
            self.exon_Stops.append(self.start + k + v)
    def __str__(self):
        return "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s,\t%s," % (self.chr, self.start, self.end, self.id, self.score, self.strand, self.cdsStart, self.cdsEnd, self.itemRGB, self.exonCount, ','.join(list(map(str,self.exonSizes))), ','.join(list(map(str, self.exonStarts))) )

class genepred():
    ''' genepred'''
    '''0-6: name(id), chrom, strand, txStart, txEnd, cdsStart, cdsEnd, \
       7-11: exonCount, exonStarts, exonEnds, score, name2, other(desc) '''
    def __init__(self,x):
        try:
            self.bin=int(x[0])
            if(self.bin < 10000):
                x=x[1:]
        except:
            pass
        self.id=x[0]
        self.chr=x[1]
        self.strand=x[2]
        if(self.strand == 1):
            self.strand="+"
        elif(self.strand == -1):
            self.strand="-"
        elif(self.strand == 0):
            self.strand="+"
        self.start=int(x[3])
        self.stop=int(x[4])
        self.cds_start=int(x[5])
        self.cds_stop=int(x[6])
        self.exon_count=int(x[7])
        self.exon_starts=x[8].split(",")
        self.exon_stops=x[9].split(",")
        for i in range(self.exon_count):
            try:
                self.exon_starts[i]=int(self.exon_starts[i])
                self.exon_stops[i]=int(self.exon_stops[i])
            except:
                pass
        self.exonSizes   = [self.exon_stops[i] - self.exon_starts[i] for i in range(self.exon_count) ]
        self.exonStartsR = [self.exon_starts[i] -self.exon_starts[0] for i in range(self.exon_count)]
        self.score=0
        try:
            self.protein_id=x[11]
            self.name2=x[11]
        except:
            self.name2="None"
        try:
            self.align_id=x[12]
        except:
            pass
    def __str__(self):
        return "%s\t%s\t%c\t%d\t%d\t%d\t%d\t%d\t%s,\t%s,\t%s" % (self.id,self.chr,self.strand,self.start,self.stop,self.cds_start,self.cds_stop,self.exon_count,",".join([str(self.exon_starts[j]) for j in range(self.exon_count)]),",".join([str(self.exon_stops[j]) for j in range(self.exon_count)]),self.name2)

def main():
    with open(sys.argv[1],'r') as fh:
        for line in fh:
            if not line.startswith("#"):
                line=line.strip().split("\t")
                g=genepred(line)
                print ('%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s,\t%s,' %(g.chr, g.start, g.stop, g.id, g.score, g.strand, g.cds_start, g.cds_stop, '0', g.exon_count, ','.join(list(map(str,g.exonSizes))), ','.join(list(map(str,g.exonStartsR)))  ))

if __name__=="__main__":
    if len(sys.argv) == 1:
        print ("python Genebed2Genepred.py genetab(genePred).file", file=sys.stderr)
    else:
        main()

