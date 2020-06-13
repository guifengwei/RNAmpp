


data <-read.table("output_forpie.txt", header=F)
colnames(data)<- c("peak_name", "gene_annotation", "Location")


require(ggplot2)
require(moonBook)
require(webr)

pdf("plot_PieDonut.pdf")
PieDonut(data,aes(pies=Location, donuts=gene_annotation), ratioByGroup=FALSE, explode=3, showRatioThreshold=0)
dev.off()
