


library(wordcloud)
library(stringr)


PGs <- read.csv('N:/VirionProject/XLb4Gradient/5mM_DSSO/SILAC_Virions/txt_14_3_3//proteinGroups.txt', sep='\t')
PGs <- GeneNameExtractMQ3(PGs)
PGs <- subset(PGs ,!(PGs$Only.identified.by.site == '+' | PGs$Potential.contaminant == "+" | PGs$Reverse == '+'))

PGs$Gene.names <- unlist(lapply(str_split(PGs$Gene.names,";") , "[[",1))

#

a <- median(log2(PGs$Ratio.H.L.normalized.WTH_1433L), na.rm = T)
b <- median(-log2(PGs$Ratio.H.L.normalized.WTL_1433H), na.rm=T)


pdf(pdfname <- "FigureS4BC.pdf", width=14, height=7, useDingbats = F)

layout(matrix(ncol=2, nrow=1,c(1,2)))

plot(log2(PGs$Ratio.H.L.normalized.WTH_1433L)-a, -log2(PGs$Ratio.H.L.normalized.WTL_1433H)-b, xlim=c(-2,2), ylim=c(-2,2),frame=F ,
     xlab='log2 SILAC fold-change virions H / L \n WT / 4_5_6mut',ylab='log2 SILAC fold-change virions L / H \n WT / 4_5_6mut', col="lightgrey", cex=2) 
points((log2(PGs$Ratio.H.L.normalized.WTH_1433L)-a)[host], (-log2(PGs$Ratio.H.L.normalized.WTL_1433H)-b)[host], xlim=c(-2,2), ylim=c(-2,2),
       xlab='log2 SILAC fold-change virions H / L \n WT / 4_5_6mut',ylab='log2 SILAC fold-change virions \n L / H ; WT / 4_5_6mut', col="#DC863B90", cex=2, pch=20) 
points((log2(PGs$Ratio.H.L.normalized.WTH_1433L)-a)[viral], (-log2(PGs$Ratio.H.L.normalized.WTL_1433H)-b)[viral], xlim=c(-2,2), ylim=c(-2,2),
       xlab='log2 SILAC fold-change virions H / L \n WT / 4_5_6mut',ylab='log2 SILAC fold-change virions \n L / H ; WT / 4_5_6mut', col="#899DA490", cex=2, pch=20) 
grid(NULL,NULL)

hits <- (log2(PGs$Ratio.H.L.normalized.WTH_1433L) > 0.67 & log2(PGs$Ratio.H.L.normalized.WTL_1433H) < -0.67) 
ISNAS <- is.na(log2(PGs$Ratio.H.L.normalized.WTH_1433L)) | is.na(log2(PGs$Ratio.H.L.normalized.WTL_1433H))
wordcloud::textplot(x = log2(PGs$Ratio.H.L.normalized.WTH_1433L)[hits&!ISNAS], y=-log2(PGs$Ratio.H.L.normalized.WTL_1433H)[hits&!ISNAS], words=PGs$Gene.names[hits&!ISNAS], new=F)

legend("bottomleft", legend = c("host protein", "viral protein"), pch = 20, col=c("#DC863B90", "#899DA490"), bty="n", cex=1.5)


write.csv(cbind(PGs[,1:7], log2(PGs$Ratio.H.L.normalized.WTH_1433L)-a, -log2(PGs$Ratio.H.L.normalized.WTL_1433H)-b), file="SourceDataExtDataFig5e.csv")

names(PGs)

PGs2 <- read.csv("N:/VirionProject/XLb4Gradient/5mM_DSSO/SILAC_Virions/txt_old_PP1_SILAC/proteinGroups.txt", sep="\t", header=T)
PGs2 <- GeneNameExtractMQ6(PGs2)
PGs2 <- subset(PGs2, !(PGs2$Potential.contaminant == "+" | PGs2$Reverse == "+" | PGs2$Only.identified.by.site =="+"))
PGs2$Gene.names <- unlist(lapply(str_split(PGs2$Gene.names,";") , "[[",1))


a <- median(log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions), na.rm=T)
b <- median(-log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions), na.rm=T)
  
  
plot(log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions)-a, -log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions)-b, ylim=c(-2,2), xlim=c(-2,2), frame=F, col="lightgrey", cex=2,
     xlab='log2 SILAC fold-change virions H / L \n WT / SILK+RVXFmut',ylab='log2 SILAC fold-change virions L / H \n WT / SILK+RVXFmut')
points((log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions)-a)[host], (-log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions)-b)[host], ylim=c(-2,2), xlim=c(-2,2), cex=2,pch=20,
     xlab='log2 SILAC fold-change virions H / L \n WT / SILK+RVXFmut',ylab='log2 SILAC fold-change virions \n L / H ; WT / SILK+RVXFmut',col="#DC863B90")
points((log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions)-a)[viral], (-log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions)-b)[viral], ylim=c(-2,2), xlim=c(-2,2), cex=2,pch=20,
       xlab='log2 SILAC fold-change virions H / L \n WT / SILK+RVXFmut',ylab='log2 SILAC fold-change virions \n L / H ; WT / SILK+RVXFmut',col="#899DA490")

hits <- (log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions)-a > 0.5 & -log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions)-b > 0.67) 
ISNAS <- is.na(log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions)) | is.na(-log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions))
wordcloud::textplot(x = (log2(PGs2$Ratio.H.L.normalized.WTH_524L_Virions)-a)[hits&!ISNAS], y=(-log2(PGs2$Ratio.H.L.normalized.WTL_524H_Virions)-b)[hits&!ISNAS], words=PGs2$Gene.names[hits&!ISNAS], new=F)
grid(NULL,NULL)

legend("bottomleft", legend = c("host protein", "viral protein"), pch = 20, col=c("#DC863B90", "#899DA490"), bty="n", cex=1.5)




graphics.off()
system(paste("open", pdfname))






