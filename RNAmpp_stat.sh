
echo "\n###  Calcuate the distribution of query sites across the unified RNA transcripts ### \n"
echo "### This shell script takes 2 input files: "
echo "### 1. query_bed file is standard 6-column bed-format file "
echo "### 2. gene_bed file is output file from Step1, either the randomTranscript or MaxORF_LongestNcRNA "
echo "### Optional for 2: custom gene annotation in 12-column genebed-format \n"

query_bed=$1
gene_bed=$2

sort -k1,1 -k2,2n $query_bed > ${query_bed}.sorted
sort -k1,1 -k2,2n $gene_bed > ${gene_bed}.sorted

mv ${query_bed}.sorted ${query_bed}
mv ${gene_bed}.sorted ${gene_bed}

echo "## Number of total query sites: "
num1=`more ${query_bed} | grep -v "^#"|wc -l`
echo $num1

echo "## perform the intersection and statistics ... \n"
intersectBed -a $query_bed -b $gene_bed -s -wo > ${query_bed}_${gene_bed}.intersect

echo "## Number of intergenic query sites: "
intersectBed -v -a $query_bed -b $gene_bed -s -wa | grep -v "^#" > ${query_bed}_${gene_bed}.intergenic
num2=`more ${query_bed}_${gene_bed}.intergenic | wc -l`
echo $num2

echo "## Calculating the intragneic overlap ... "
python Scripts/stat_bed_in_genebed.py -g $gene_bed --intersect ${query_bed}_${gene_bed}.intersect -o output.txt

echo "## preparing the xls for piechart ... "
more output.txt | cut -f1-6,20,21 > output.temp.txt
python Scripts/combine_anno.py output.temp.txt > output_1.txt
awk '{FS=OFS="\t"}{print $4, "intergenic", "intergenic"}' ${query_bed}_${gene_bed}.intergenic > output_2.txt
cat output_1.txt output_2.txt > output_forpie.txt
rm output_1.txt output_2.txt output.temp.txt

## plot
echo "\n### making enrichment plot for protein_coding gene, ncRNA, and intron"
Rscript plot_PieDonut.R
Rscript plot_frac.R

