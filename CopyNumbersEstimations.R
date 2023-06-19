
setwd("N:/VirionProject/NewDataset/")
sapply(paste('./functions/',list.files(pattern="[.]R$", path="./functions"),sep='/'), source);


library(plotly)
library(stringr)
library(RColorBrewer)

setwd("N:/VirionProject/iBAQs_Virion_Prep/txt")


### first biological replicate
MQ <- read.csv("proteinGroups.txt", sep="\t", stringsAsFactors = F, header = T)
MQ <- subset(MQ, !(MQ$Potential.contaminant == "+" | MQ$Reverse == "+" | MQ$Only.identified.by.site == "+"))


### second biological replicate
MQ2 <- read.csv("N:/VirionProject/iBAQs_Virion_Prep/txt_triplicates/proteinGroups.txt", sep="\t", stringsAsFactors = F, header = T)
MQ2 <- subset(MQ2, !(MQ2$Potential.contaminant == "+" | MQ2$Reverse == "+" | MQ2$Only.identified.by.site == "+"))

MQ <- GeneNameExtractMQ(PGs=MQ)

MQ$Gene.names <- unlist(lapply(str_split(MQ$Gene.names, ";"), "[[",1))

MQ2 <- GeneNameExtractMQ(PGs=MQ2)

MQ2$Gene.names <- unlist(lapply(str_split(MQ2$Gene.names, ";"), "[[",1))

MQ$Gene.names[grepl('ABV71575.1',MQ$Fasta.headers)] <- "UL45"
MQ2$Gene.names[grepl('ABV71575.1',MQ2$Fasta.headers)] <- "UL45"


### get HighConfidenceViral 


MQ$iBAQ[MQ$iBAQ == 0] <- NA


iBAQ1 <- MQ[,grepl("iBAQ",names(MQ))][,4:5]
iBAQ2 <- MQ2[,grepl("iBAQ",names(MQ2))][,2:4]



Cap <- c("UL32", "UL46", "UL48A", "UL77", "UL85", "UL86", "UL93", "UL104")

CapVal <- c(960, 320, 960, 24, 650,960, 12,12)


col1 <- brewer.pal(n=5, name = 'Set2')

CopyMaker <- function(iBAQsVir, MQ=MQ, add=F, col1=col1){
  In <- MQ$Gene.names %in% Cap
if (add == F){
  plot(y=log10(iBAQsVir[In]), x=CapVal, log = "x", xlab="protein copies per virion", xlim=c(1,1000), ylim=c(6.5,10), ylab="log10 iBAQ", frame=F, col=col1, pch=20, cex=1.5)
  grid(NULL, NULL)
}
  else {
    
    points(y=log10(iBAQsVir[In]), x=CapVal,col=col1, pch=20, cex=1.5)
    
  }
    
  
  lx <- lm(log10(iBAQsVir[In])~log10(CapVal))
  abline(lx, lwd=2, lty=3,col=col1)
  
 # wordcloud::textplot(y=log10(iBAQsVir)[In], x=CapVal, words=MQ$Gene.names[In], new=F)
#  text(x=c(10), y=c(9.25), labels =paste("log10 iBAQ = ",round(lx$coefficients[2],3),"* log10 copy number + ", round(lx$coefficients[1],2), sep=""), cex=1.5)
  Copies <- 10^((log10(iBAQsVir)-lx$coefficients[1])/lx$coefficients[2])
  
  return(list(Copies, lx, iBAQsVir[In], log10(CapVal)))
  
}


MQ2$Gene.names[MQ2$Gene.names %in% Cap]

C1 <- CopyMaker(iBAQsVir = iBAQ1$iBAQ.Vir_1, MQ=MQ, col1=col1[1], add=F)
C2 <- CopyMaker(iBAQsVir = iBAQ1$iBAQ.Vir_3, MQ=MQ, col1=col1[2], add=T)
C3 <- CopyMaker(iBAQsVir = iBAQ2$iBAQ.V1, MQ=MQ2, col1=col1[3], add=T)
C4 <- CopyMaker(iBAQsVir = iBAQ2$iBAQ.V2, MQ=MQ2, col1=col1[4], add=T)
C5 <- CopyMaker(iBAQsVir = iBAQ2$iBAQ.V3, MQ=MQ2, col1=col1[5], add=T)

legend('topleft', legend = c("rep1", "rep2", "rep3", "rep4", "rep5"), col=col1, lwd = 2,pch = 20, lty=3, bty='n', cex=2.5)

