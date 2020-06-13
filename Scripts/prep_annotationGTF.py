#!/usr/bin/python
# Programmer: Wei Guifeng, <guifengwei@gmail.com>
#-- coding: utf-8 --
#Last-modified: 12 Jun 2020 11:55:41

import sys, os, string, re
from collections import defaultdict

'''
subset of GTF based on gene_type if exsits, this Program should take less than 1 min for all kinds of preparation
'''
def prep_anno(gtf):
    ''' main scripts '''
    ### patterns
    gene_type_pattern = re.compile('gene_(bio)?type \"([A-Za-z0-9\_\-\.\(\)]+)\"')
    trans_id_pattern = re.compile('transcript_id \"([A-Za-z0-9\_\-\.\(\)]+)\"')
    gene_id_pattern = re.compile('gene_id \"([A-Za-z0-9:\'\_\-\.\(\)\[\]]+)\"')
    gene_name_pattern = re.compile('gene_name \"([A-Za-z0-9:\'\_\-\.\(\)\[\]]+)\"')
    ###
    geneid2name={}
    geneid2genetype={}
    geneid2transcripts=defaultdict(list)
    ## Types filtration
    gene_type=""
    unwanted_type=["pseudogene", "IG_", "TR_", "miRNA", "rRNA", "TEC", "tRNA", "sRNA", "scRNA"]
    ## output the file, genetype is optional
    fo1 = open("geneid_genename_transcriptid_genetype.txt", "w")
    fo2 = open("filtered_"+gtf, "w")
    for line in open(gtf, 'r'):
        if line.startswith('#'): continue
        line = line.strip()
        # If it has gene_type, it will pass through this gene_type filter
        if "gene_type" in line or "gene_biotype" in line:
            gene_type = re.findall(gene_type_pattern, line)[0][1]
            if any([i in gene_type for i in unwanted_type]):
                continue
        fo2.write(line+"\n")
        items = line.split("\t")
        if items[2] == "exon":
            gene_id = re.findall(gene_id_pattern, line)[0]
            gene_name = re.findall(gene_name_pattern, line)[0]
            transcript_id = re.findall(trans_id_pattern, line)[0]
            if gene_id not in geneid2name:
                geneid2name[gene_id] = gene_name
            if gene_id not in geneid2genetype:
                geneid2genetype[gene_id] = gene_type
            geneid2transcripts[gene_id].append(transcript_id)
    #### output the genename file
    for gid, gname in geneid2name.items():
        tid = "|".join(set(geneid2transcripts[gid]))
        gtype=geneid2genetype[gid]
        fo1.write("{0}\t{1}\t{2}\t{3}\n".format(gid, gname, tid, gtype))

if __name__ == "__main__":
    if len(sys.argv) >= 1:
        prep_anno(gtf=sys.argv[1])
    else:
        print('\nUsage: prep_annotationGTF.py *.gtf\n', file=sys.stderr)


