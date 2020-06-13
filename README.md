# RNAmpp
RNA metaProfile plot (RNAmpp), the pipeline to generate the distribution on transcriptome (RNA). Written in Shell, Python3, and R language. Working with human, mouse, and drosophila genome annotations.

##### Dependency:
1. running in Linux with Python3 and R.
2. Bedtools >=v2.25.0
3. R package: ggplot2, moonBook, webr
4. Python package: numpy

##### How to run:

1. Isoform selection. 
   mRNA or lncRNAs have multiple isoforms, the prep step aims to select the representative isoform. This can be the MaxORF_LongestNcRNA, or randomTranscript by setting the mode m or r in the shell script. In addition, the user can also skip this step, with their own representativeTranscript for the next step.
   ```
   sh RNAmpp_prep.sh GENCODE_vM22.annotation.gtf
   ```
  This step usually takes less than 100 seconds, tested in GENCODE and UCSC annoated gtf files.
  
 2. Calculate the relative distribution and generate plots.
    This step requires two files as input, (1) query_bed must be a standard 6-column bed-format file (2) genebed file which is the output from Step1, either the randomTranscript or MaxORF_LongestNcRNA, also can be user-defined Transcript.
    ```
    sh RNAmpp_stat.sh query.bed MaxORF_LongestNcRNA.genebed
    ```
 
 ##### Result
 
 
 
 
 
 
 ##### Citation
 
 Please contact me at this monment.
  
