# RNAmpp
RNA metaProfile plot (RNAmpp), the pipeline to generate the distribution on transcriptome (RNA). Written in Shell, Python3, and R language. Working with human, mouse, and drosophila genome annotations.

##### Dependency:
1. running in Linux with Python3 and R.
2. Bedtools >=v2.25.0
3. R package: ggplot2, moonBook, webr
4. Python package: numpy

##### How to run:

1. Isoform selection. (mRNA or lncRNAs have multiple isoforms, the prep step aims to select the representative isoform. This can be the MaxORF_LongestNcRNA, or randomTranscript by setting the mode m or r in the shell script. In addition, the user can also skip this step, with their own representativeTranscript for the next step)

  ```
  sh RNAmpp_prep.sh GENCODE_vM22.annotation.gtf
  ```
  This step usually takes less than 100 seconds, tested in GENCODE and UCSC annoated gtf files.
