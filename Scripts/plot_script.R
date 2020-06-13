
library(ggplot2)
Protein_coding <- append(UTR5, CDS, after=length(UTR5))
Protein_coding <- append(Protein_coding, UTR3, after=length(Protein_coding))

PcG <- data.frame("Coord"=1:length(Protein_coding), "Fraction"=Protein_coding/sum(Protein_coding))
noncoding_RNA <- data.frame("Coord"=1:length(ncRNA), "Fraction"=ncRNA/sum(ncRNA))
INTRON <- data.frame("Coord"=1:length(Intron), "Fraction"=Intron/sum(Intron))

pdf("plot_frac.pdf", width=12.5, height=5.25)
qplot(Coord, Fraction, data=PcG,   geom=c("line")) + geom_vline(xintercept = c(length(UTR5), length(UTR5)+length(CDS)), colour="red",linetype = "longdash") 
qplot(Coord, Fraction, data=noncoding_RNA, geom=c("line"))
qplot(Coord, Fraction, data=INTRON,        geom=c("line"))
dev.off()
