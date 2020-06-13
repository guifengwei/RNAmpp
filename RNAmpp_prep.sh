
echo "\n### Gene Annotation preparation, filtration, and selection ###"
echo "### This step usually takes <= 100 sec"
date +"%T"

gtf=$1
name="$(echo $gtf | awk -F "/" '{print $NF}')"
echo "\n### The filtered gtf file is: "
echo " filtered_$gtf"

python Scripts/prep_annotationGTF.py $gtf

echo "\n### Converting the filtered GTF to GeneBed file\n"

Utlities/gtfToGenePred "filtered_"$gtf ${name}.genePred
python Scripts/GenePred2GeneBed.py  ${name}.genePred > ${name}.GeneBed
rm ${name}.genePred

echo "### Representative isoform selecting ...\n "
python Scripts/Select_RepresentativeTranscript_For_Gene.py -g ${name}.GeneBed -m m

date +"%T"
