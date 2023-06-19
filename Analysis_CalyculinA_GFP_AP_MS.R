
library(stringr)
library(beeswarm)



setwd('N:/VirionProject/Caly300nM/txt')

ULX <- read.csv("OutputPerseus.txt", sep="\t", stringsAsFactors = F)

#GeneNameExtractMQ3(ULX)$Gene.names

ULX$Gene.names[grepl(x=ULX$Fasta.headers, pattern = 'herpesvirus')] <- str_remove_all(str_extract(pattern = " [URI][A-Z0-9.,]{0,5} ",string=ULX$Fasta.headers)[grepl(pattern='herpesvirus', ULX$Fasta.headers)], pattern = " ")

Imp <- read.csv(file = 'OutputPerseusNoImp.txt', sep="\t")
Imp <- GeneNameExtractMQ(Imp)
Imp$Gene.names <- unlist(lapply(str_split(Imp$Gene.names, ";"), "[[", 1))


#pdf(pdfname <- "Stoichiometry_Interactors.pdf", width=30, height=10)
#layout(matrix(ncol=3, nrow=1, c(1,2,2)))

plot(Imp$Difference[Imp$X.Log.P.value.!=0], Imp$X.Log.P.value.[Imp$X.Log.P.value.!=0], cex=2,
     frame=F, col="darkgrey", main="pUL32 / pp150 AP- MS", xlab="log2 LFQ (300 nM Calyculin A vs DMSO)", ylab="-log 10 t-test p-value", xlim=c(-5.5,5.5))
#points(Imp$Student.s.T.test.Difference.GFP_WT[Imp$Gene.names %in% IT], Imp$X.Log.Student.s.T.test.p.value.GFP_WT[Imp$Gene.names %in% IT], col="orange", pch=20)
#points(Imp$Student.s.T.test.Difference.GFP_WT[Imp$Gene.names %in% c(OT,BMOT,Diffuse)], Imp$X.Log.Student.s.T.test.p.value.GFP_WT[Imp$Gene.names %in% c(OT,BMOT,Diffuse)], col="mediumseagreen", pch=20)
abline(h=2,lty=3,lwd=3)
abline(v=2,lty=3,lwd=3)
abline(v=-2,lty=3,lwd=3)
HitsUL32 <- Imp$Difference > 2 & Imp$X.Log.P.value.>2 | Imp$Gene.names=="UL32"
wordcloud::textplot(Imp$Difference[HitsUL32], Imp$X.Log.P.value.[HitsUL32], words = Imp$Gene.names[HitsUL32], new=F)
grid(NULL,NULL)


## write.csv(file="SourceDataExtDataFig9c.csv", Imp)


#iBAQs <- read.csv("N:/VirionProject/IP_GFPPP150_14_3_3_PP1/txt/proteinGroups.txt", sep = "\t")


Phos <- read.csv("N:/VirionProject/Caly300nM/txt/Phospho (STY)Sites.txt", sep = "\t")
evi <- read.csv("N:/VirionProject/Caly300nM/txt/peptides.txt", sep = "\t")

NPhos <- subset(Phos,!Phos$Potential.contaminant == "+" )

eviUL32 <- evi[grepl("UL32", evi$Proteins),]
PhosUL32 <- Phos[grepl("UL32",Phos$Proteins),]

eviUL32_min <- rowSums(cbind.data.frame(eviUL32$Intensity.Minus_1, eviUL32$Intensity.Minus_2, eviUL32$Intensity.Minus_3), na.rm=T)
eviUL32_min<- eviUL32_min[eviUL32_min != 0]
eviUL32_plus <- rowSums(cbind.data.frame(eviUL32$Intensity.Plus_1, eviUL32$Intensity.Plus_2, eviUL32$Intensity.Plus_3), na.rm=T)
eviUL32_plus<- eviUL32_plus[eviUL32_plus != 0]

PhosUL32_min <- rowSums(cbind.data.frame(PhosUL32$Intensity.Minus_1, PhosUL32$Intensity.Minus_2, PhosUL32$Intensity.Minus_3), na.rm=T)
PhosUL32_min <- PhosUL32_min[PhosUL32_min != 0]
PhosUL32_plus <- rowSums(cbind.data.frame(PhosUL32$Intensity.Plus_1, PhosUL32$Intensity.Plus_2, PhosUL32$Intensity.Plus_3), na.rm=T)
PhosUL32_plus <- PhosUL32_plus[PhosUL32_plus != 0]


pdf(pdfname <- "Fig5d.pdf", width=7, height=9)
b <- boxplot(list(log2(PhosUL32_min), log2(PhosUL32_plus), log2(eviUL32_min), log2(eviUL32_plus)), main="pp150/UL32 peptides",
        notch=T, ylab="log2 sum peptide intensity", ylim=c(18,36),
        names = c("phospho \n DMSO ","phospho \n Calyculin A ", "unmodified \n DMSO ", "unmodified \n Calyculin A "), horizontal =F, frame=F)
beeswarm(list(log2(PhosUL32_min), log2(PhosUL32_plus), log2(eviUL32_min), log2(eviUL32_plus)) , corral = "wrap", add=T, col=c("#899DA490","#DC863B90"), pch=20, cex=2)
text(y=b$stats[5,]+0.75, x=c(1:4),labels = paste("n=\n",b$n, sep=""))
graphics.off()
system(paste("open", pdfname))


HitsUL32 <- Imp$Difference > 2 & Imp$X.Log.P.value.>2 | Imp$Gene.names %in% GOI

plot(Imp$Difference[Imp$X.Log.P.value.!=0], Imp$X.Log.P.value.[Imp$X.Log.P.value.!=0], cex=2,
     frame=F, col="darkgrey", main="pUL32 / pp150 AP- MS", xlab="log2 difference (300 nM Calyculin A vs DMSO)", ylab="-log 10 P (t-test)", xlim=c(-3,5.5))
points(Imp$Difference[Imp$X.Log.P.value.!=0 & HitsUL32], Imp$X.Log.P.value.[Imp$X.Log.P.value.!=0 &HitsUL32], cex=2,
     col="red")
#points(Imp$Student.s.T.test.Difference.GFP_WT[Imp$Gene.names %in% c(OT,BMOT,Diffuse)], Imp$X.Log.Student.s.T.test.p.value.GFP_WT[Imp$Gene.names %in% c(OT,BMOT,Diffuse)], col="mediumseagreen", pch=20)
abline(h=2,lty=3,lwd=3)
abline(v=2,lty=3,lwd=3)
abline(v=-2,lty=3,lwd=3)

wordcloud::textplot(Imp$Difference[HitsUL32], Imp$X.Log.P.value.[HitsUL32], words = Imp$Gene.names[HitsUL32], new=F)
grid(NULL,NULL)

graphics.off()
system(paste("open", pdfname))




wilcox.test(log2(eviUL32_min), log2(eviUL32_plus))
wilcox.test(log2(PhosUL32_min), log2(PhosUL32_plus))


