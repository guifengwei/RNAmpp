#!/usr/bin/python
# Programmer: Guifeng Wei, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 11 Jun 2020 16:58:41

import sys, os, argparse, string
from Select_RepresentativeTranscript_For_Gene import ORF_UTR_size_of_v2 as ORF_UTR_size_of
from UCSC_Class import Bed,genebed
import numpy as np
#import warnings
#warnings.filterwarnings("ignore")

def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --gene *.genebed --intersect *bed_*genebed.intersect')
    p.add_argument('-g','--gene',dest='genefile',metavar='Gene',type=str,  required=True, help="The genebed file. Format: GeneBed")
    p.add_argument('--intersect',dest='intersect', metavar='intersectfile', required=True, type=str, help="The intersect file between query.bed(6-col) and genebed(12-col). Total is 19-col")
    p.add_argument('-o','--output',dest='output',metavar='output', type=str, default="output.txt", help="The output files keeps original file but also report the relative result to annotation, defalut=output.txt")
    if len(sys.argv) == 1 :
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def gene_utr_cds_length(genebed_file):
    UTR5 = []
    ORF = []
    UTR3 = []
    ncRNA = []
    for line in open(genebed_file, "r"):
        if not line.startswith("#"):
            g = genebed( line.strip().split("\t") )
            utr5, cds, utr3 = ORF_UTR_size_of(g)
            if cds < 10:
                ncRNA.append( utr5+cds+utr3 )
            elif cds > 10:
                UTR5.append(utr5)
                ORF.append(cds)
                UTR3.append(utr3)
            else:
                pass
    mean_utr5 = float( np.where( len(UTR5) == 0, 0, np.mean(UTR5) ) )
    mean_utr3 = float( np.where( len(UTR3) == 0, 0, np.mean(UTR3) ) )
    mean_cds  = float( np.where( len(ORF)  == 0, 0, np.mean(ORF) ) )
    mean_ncRNA= float( np.where( len(ncRNA)== 0, 0, np.mean(ncRNA) ) )
    print("##########################")
    print("Average length for features: UTR5, CDS, UTR3, and ncRNA")
    print(mean_utr5, mean_cds, mean_utr3, mean_ncRNA)
    return mean_utr5, mean_cds, mean_utr3, mean_ncRNA

def bed_in_bedlist(bed, bedlist):
    ''' '''
    length = sum([exon.length() for exon in bedlist])
    for i, exon in enumerate(bedlist):
        if bed.overlap(bed, exon):
            if i == 0:
                f = 1.0*(bed.end - exon.start)/length
            else:
                pre_length = sum([bedlist[x].length() for x in range(i) ])
                f = 1.0*(pre_length + bed.end - exon.start)/length
            ### strand
            if bed.strand == exon.strand and bed.strand == "+":
                return f
            elif bed.strand == exon.strand and bed.strand == "-":
                return 1-f
            else:
                pass
            #### using for test
            if not f:
                print("#Something is wrong", file=sys.stderr)
                print(bed)
                print(exon)

def bed_in_Intron(bed, Introns):
    '''return the relative position into individual intron, without cumulation '''
    for intron in Introns:
        if bed.overlap(bed, intron):
            f = 1.0*(bed.end - intron.start)/intron.length()
            if bed.strand == intron.strand and bed.strand == "+":
                return f
            elif bed.strand == intron.strand and bed.strand == "-":
                return 1-f
            else:
                return -1
    return -1

def utr5_cds_utr3_rPos(intersect_file):
    ''' if the overlapped genebed is protein_coding gene and return the relative position to the feature start '''
    ''' for example: UTR5 is [0,1220], and the single-point is 125, so it returns 125/1220=0.10245902'''
    UTR5_rpos = []
    CDS_rpos  = []
    UTR3_rpos = []
    ncRNA_rpos= []
    Intron_rpos=[]
    for line in open(intersect_file, "r"):
        if not line.startswith("#"):
            line = line.strip().split("\t")
            b = Bed(line[0:6])
            g = genebed(line[6:18])
            ## print("\t".join(line)) ## this is test line
            utr5, cds, utr3 = ORF_UTR_size_of(g)
            ## print(utr5, cds, utr3) ## this is test line
            #### Introns
            if g.Introns():
                f = bed_in_Intron(b, g.Introns())
                ## print("{0}\t{1}".format("test ", f)) ## this si test line
                if f>=0:
                    Intron_rpos.append(f)
                    if cds<10:
                        fo.write("{0}\t{1}\t{2}\n".format("\t".join(line), "ncRNA", "intron"))
                    else:
                        fo.write("{0}\t{1}\t{2}\n".format("\t".join(line), "protein_coding", "intron"))
                    continue
            ########
            if cds < 10:
                ncRNA_rpos.append( bed_in_bedlist(b, g.Exons()) )
                fo.write("{0}\t{1}\t{2}\n".format("\t".join(line), "ncRNA", "exon"))
            else:
                try:
                    if g.CDS().start <= b.start < b.end <= g.CDS().end:
                        CDS_rpos.append(  bed_in_bedlist( bed=b, bedlist=g.CDS().Exons()  ) )
                    elif utr5 > 0 and g.UTR5().start <= b.start < b.end <= g.UTR5().end:
                        UTR5_rpos.append( bed_in_bedlist( bed=b, bedlist=g.UTR5().Exons() ) )
                    elif utr3 > 0 and g.UTR3().start <= b.start < b.end <= g.UTR3().end:
                        UTR3_rpos.append( bed_in_bedlist( bed=b, bedlist=g.UTR3().Exons() ) )
                    else:
                        pass
                    fo.write("{0}\t{1}\t{2}\n".format("\t".join(line), "protein_coding", "exon"))
                except AttributeError:
                    pass
    print("####################")
    print("# Intron peaks:")
    print(len(Intron_rpos))
    print("# ncRNA peaks: ")
    print(len(ncRNA_rpos))
    print("# 5UTR peaks:")
    print(len(UTR5_rpos))
    print("# CDS peaks:")
    print(len(CDS_rpos))
    print("# 3UTR peaks:")
    print(len(UTR3_rpos))
    return UTR5_rpos, CDS_rpos, UTR3_rpos, ncRNA_rpos, Intron_rpos

def bin_fraction(nBin, frac_list):
    ''' according to the step_frac, and return the frequency occurring within this step size '''
    freq = []
    step = 1.0/nBin
    for i in range(nBin):
        freq.append( len([ j for j in frac_list if i*step<= j< (i+1)*step ]) )
    return freq

def main():
    ''' main scripts '''
    args = parse_argument()
    genebed_file = args.genefile
    intersect_file = args.intersect
    output = args.output
    global fo
    fo = open(output, "w")
    Routput = open("plot_frac.R", "w")
    ## 1. get the mean size of feature length
    mean_utr5, mean_cds, mean_utr3, mean_ncRNA = gene_utr_cds_length(genebed_file)
    Routput.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format("#", "average_length: UTR5, CDS, UTR3, ncRNA ", mean_utr5, mean_cds, mean_utr3, mean_ncRNA))
    ## 2. get the relative position
    UTR5_rpos, CDS_rpos, UTR3_rpos, ncRNA_rpos, Intron_rpos = utr5_cds_utr3_rPos(intersect_file)
    
    common_step_width = 0
    for i in range(20):
        x = mean_utr5 * mean_cds * mean_utr3 * mean_ncRNA
        if 10<= x/10**i < 100:
            common_step_width = x/10**i
    Routput.write("{0} {1} {2}\n".format("#", "common_step_width acorss features: ", common_step_width))
    print("{0} {1} {2}".format("#", "common_step_width acorss features: ", common_step_width))

    nBin_utr5 = int(round(1.0*mean_utr5/common_step_width))
    nBin_cds  = int(round(1.0*mean_cds/common_step_width))
    nBin_utr3 = int(round(1.0*mean_utr3/common_step_width))
    nBin_ncRNA =int(round(1.0*mean_ncRNA/common_step_width))

    Routput.write( "{0}{1}{2}\n".format( "UTR5<-c(",  ",".join(list(map(str, bin_fraction(nBin=nBin_utr5,  frac_list=UTR5_rpos)))), ")"))
    Routput.write( "{0}{1}{2}\n".format( "CDS <-c(",  ",".join(list(map(str, bin_fraction(nBin=nBin_cds,   frac_list=CDS_rpos)))), ")"))
    Routput.write( "{0}{1}{2}\n".format( "UTR3<-c(",  ",".join(list(map(str, bin_fraction(nBin=nBin_utr3,  frac_list=UTR3_rpos)))), ")"))
    Routput.write( "{0}{1}{2}\n".format( "ncRNA<-c(", ",".join(list(map(str, bin_fraction(nBin=nBin_ncRNA, frac_list=ncRNA_rpos)))), ")"))
    Routput.write( "{0}{1}{2}\n".format( "Intron<-c(",",".join(list(map(str, bin_fraction(nBin=40, frac_list=Intron_rpos)))), ")"))
    for line in open("Scripts/plot_script.R", "r"):
        Routput.write( line.strip()+"\n" )
#    3. calculate the exon, intron, antisense, and intergenic fraction as piechart
#    4. calculate the protein_coding and ncRNA, fraction as barplot
    fo.close()
    Routput.close()

if __name__ == "__main__":
    main()


