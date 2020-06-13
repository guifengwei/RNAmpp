#!/usr/bin/python
# Programmer: Wei Guifeng, <Guifeng.Wei@bioch.ox.ac.uk>
#-- coding: utf-8 --
#Last-modified: 12 Jun 2020 17:50:37

import sys, os, argparse, string, random
#from itertools import izip
from UCSC_Class import genebed

'''
    Aim: Fetch the representive of the Genes, which contains multiple isoforms
    for noncoding genes, take the longest isoform
    for protein-coding genes, take the longest-ORF isoforms, if the ORFs are equal, take the one with longest size 

    two files are needed.
    file1: Genes.Genebed, which is a Genebed-format file, describes the gene structure
    file2: geneid_genename_transcriptid_genetype.txt, in the third column, the transcript ID were seperated by "|"
    
'''
def Iter_filehandle(filename=""):
    ''' split every lines into iterable lines'''
    for line in open(filename, 'r'):
        if not line.startswith("#"):
            line = line.strip().split("\t")
            yield line

def ORF_UTR_size_of(gene):
    ''' main scripts '''
    g = gene
    #print g
    if g == []:
        print("Transcript input error ! ", file=sys.stderr)
    ## assume that positive strand
    ## find CDS locating which exon
    cds_start_exon = 0
    for idx, start in enumerate(g.exon_starts):
        if g.cdsStart >= start:
            cds_start_exon += 1
        else:
            break
    UTR_f = sum(g.exonSizes[0:(cds_start_exon-1)]) + g.cdsStart - g.exon_starts[cds_start_exon-1]
    #print cds_start_exon
    cds_end_exon = 0
    ## count from the last exon end
    for idx, End in enumerate(g.exon_stops[::-1]):
        if g.cdsEnd <= End:
            cds_end_exon += 1
        else:
            break
    #print g.exon_stops
    UTR_b = sum(g.exonSizes[::-1][0:(cds_end_exon - 1)]) +  g.exon_stops[::-1][cds_end_exon - 1] - g.cdsEnd
    ORF = sum(g.exonSizes) - UTR_f - UTR_b
    UTR5 = 0
    UTR3 = 0
    if g.strand == '+':
        UTR5, UTR3 = UTR_f, UTR_b
    elif g.strand == '-':
        UTR5, UTR3 = UTR_b, UTR_f
    else:
        print('Strand Error ! ', file=sys.stderr)
    #print g.id, UTR5, ORF, UTR3, g.strand
    return UTR5, ORF, UTR3

def ORF_UTR_size_of_v2(gene):
    ''' main scripts '''
    g = gene
    if g == []:
        print("Transcript input error ! ", file=sys.stderr)
    try:
        UTR5 = sum(g.UTR5().exonSizes)
    except:
        UTR5 = 0
    try:
        UTR3 = sum(g.UTR3().exonSizes)
    except:
        UTR3 = 0
    try:
        ORF = sum(g.CDS().exonSizes)
    except:
        ORF = 0
    #print(g.id, UTR5, ORF, UTR3)
    return UTR5, ORF, UTR3

def length_of_transcript(transcript=[]):
    '''
    return the length of the transcript, pre-mRNA length
    '''
    L = sum(transcript.exonSizes)
    return L

def parse_argument():
    ''' argument parser'''
    p=argparse.ArgumentParser(description='example: %(prog)s --gene *.GeneBed -m [m, r]')
    p.add_argument('-g', '--gene', dest='genebed', metavar='Gene', type=str, required=True, help="The Gene annotation file. Format: Genebed")
    p.add_argument('-m', '--mode', dest='mode',  choices=["m", "r"], default = "m", type=str, help="Isoform selection mode. m=MaxORF_LongestNcRNA, r=randomTranscript. default: m")
    if len(sys.argv) == 1:
        sys.exit(p.print_help())
    args=p.parse_args()
    return args

def main():
    ''' '''
    size_of_ORF = {}
    size_of_RNA = {}
    transcript_bed = {}

    args = parse_argument()
    mode=args.mode
    genefile=args.genebed

    if mode=="m":
        fo=open(genefile.strip(".GeneBed")+".MaxORF_LongestNcRNA.GeneBed", "w")
        print("## Searching the maxORF and Longest Transcript ... ", file=sys.stderr)
    elif mode=="r":
        fo=open(genefile.strip(".GeneBed")+".randomTranscript.GeneBed", "w")
        print("## Searching the random Transcript ... ", file=sys.stderr)
    else:
        print("## Error in the isoform selection mode ", file=sys.stderr)

    print("## Calculating the size of ORF and gene length ... ", file=sys.stderr)
    ## loading the GeneBed file
    for x in Iter_filehandle(filename=genefile):
        gene = genebed(x)
        ## gene.id in genebed class is actually transcript_id
        transcript_bed[gene.id] = x
        #### for the gene without "NM_" or "NR_" in name
        size_of_ORF[gene.id] = ORF_UTR_size_of_v2(gene)[1]
        size_of_RNA[gene.id] = length_of_transcript(gene)
    ##
    ### geneid_genename_transcriptid_genetype.txt file 
    for x in Iter_filehandle(filename = "geneid_genename_transcriptid_genetype.txt"):
        #### Geneid, GeneName, Transcripts, GeneType
        genename = x[1]
        transcripts = x[2].split('|')
        longest_one = ''
        orf = []
        #### gene_type
        try:
            if x[3]:
                pass
        except IndexError:
            x.append("")
        ##
        for t in transcripts:
            try:
                orf.append(size_of_ORF[t])
            except:
                print("## Error for calculating the ORF size for, ", t, file=sys.stderr)
        ## max ORF, if equal orf, long transcript
        m = max(orf)
        index = []
        for idx, size in enumerate(orf):
            if size == m:
                index.append(idx)
        if len(index) == 1:
            longest_one = transcripts[orf.index(m)]
        else:
            RNALength = []
            for i, t in enumerate(transcripts):
                if i in index:
                    RNALength.append(size_of_RNA[t])
            longest_one = transcripts[ index[RNALength.index(max(RNALength)) ] ]
        random_index = random.randint(0,len(transcripts)-1)
        random_one = transcripts[random_index]
        if mode=="m":
            select_one = longest_one
        elif mode=="r":
            select_one = random_one
        else:
            print("## Error in the isoform selection mode ", file=sys.stderr)
        select_genebed = transcript_bed[select_one]
        if x[0] == x[1] and not x[3]:
            fo.write("{0}\t{1}\t{2}\n".format( "\t".join(select_genebed[0:3]), select_one+"|"+genename,"\t".join(select_genebed[4:]) ))
        elif x[0] == x[1] and x[3]:
            fo.write("{0}\t{1}\t{2}\n".format( "\t".join(select_genebed[0:3]), select_one+"|"+genename+"|"+x[3],"\t".join(select_genebed[4:]) ))
        elif x[0] != x[1] and not x[3]:
            fo.write("{0}\t{1}\t{2}\n".format( "\t".join(select_genebed[0:3]), select_one+"|"+genename+"|"+x[0],"\t".join(select_genebed[4:]) ))
        else:
            fo.write("{0}\t{1}\t{2}\n".format( "\t".join(select_genebed[0:3]), select_one+"|"+genename+"|"+x[0]+"|"+x[3],"\t".join(select_genebed[4:]) ))
    print("## Done ",  file=sys.stderr)
    fo.close()

if __name__ == "__main__":
    main()


