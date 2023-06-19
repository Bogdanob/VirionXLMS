library(stringr)
library(beeswarm)
#### this is checking PhosphoPeptide Data PP1 mutant



PGs <- read.csv('//agms/M/Fusion/2022/BB/SUBMISSION_PRIDE_Virions/Virion_Phosphoprotein_Abundance_RVxF+SILKmut_vs_WT/combined/txt/proteinGroups.txt', sep='\t')
PGs <- GeneNameExtractMQ3(PGs)
PGs <- subset(PGs ,!(PGs$Only.identified.by.site == '+' | PGs$Potential.contaminant == "+" | PGs$Reverse == '+'))

PGs$Gene.names <- unlist(lapply(str_split(PGs$Gene.names,";") , "[[",1))

Phos <- read.csv('//agms/M/Fusion/2022/BB/SUBMISSION_PRIDE_Virions/Virion_Phosphoprotein_Abundance_RVxF+SILKmut_vs_WT/combined/txt/Phospho (STY)Sites.txt', sep='\t')
Phos <- GeneNameExtractMQ3(Phos)

Phos <- subset(Phos ,!(Phos$Potential.contaminant == "+" | Phos$Reverse == '+'))
Phos <- subset(Phos, Phos$Localization.prob.WTH_524L_phos > 0.75 | Phos$Localization.prob.WTL_524H_phos > 0.75)


### remove mutated phosphosites
Phos <- subset(Phos, !Phos$id %in% c(3542,3539,170,185,184,166,153,152,183,182,181,180))

plot(log2(Phos$Ratio.H.L.normalized.WTH_524L_phos), -log2(Phos$Ratio.H.L.normalized.WTL_524H_phos))
#identify(log2(Phos$Ratio.H.L.normalized.WTH_524L_phos_vir), -log2(Phos$Ratio.H.L.normalized.WTL_524H_phos_vir),labels =  paste(Phos$Gene.names, Phos$Positions.within.proteins))


MERGIS <- merge(Phos, PGs, by="Gene.names", all=T)

H.L.WTvs524 <- log2(MERGIS$Ratio.H.L.WTH_524L_phos.x) - log2(MERGIS$Ratio.H.L.WTH_524L_WVL.y)
H.L.524vsWT <- log2(MERGIS$Ratio.H.L.WTL_524H_phos.x) - log2(MERGIS$Ratio.H.L.WTL_524H_WVL.y)

AVIS <- rowMeans(cbind(H.L.WTvs524, -H.L.524vsWT), na.rm=F)

host <- !grepl(pattern="lcl", x = MERGIS$Fasta.headers.x)
viral <- grepl(pattern="lcl", x = MERGIS$Fasta.headers.x)

pdf(pdfname <- "Compare_UL32_PhosphoSites.pdf", width=6, height=9, useDingbats = F)

boxplot(list(AVIS[host],AVIS[viral & !MERGIS$Gene.names == "UL32"], AVIS[MERGIS$Gene.names == "UL32"]), frame=F,
        notch=T, ylab="log2 Phospho fold-change WT/SILK+RVxFmut; Virions", col="white", ylim=c(-1.6,1.6), outline = F,
        names = c("phosphosites in \n host proteins", "phosphosites in \n viral proteins", "phosphosites in \n pp150/UL32"))
beeswarm(list(AVIS[host],AVIS[viral & !MERGIS$Gene.names == "UL32"], 
              AVIS[MERGIS$Gene.names == "UL32"]) , corral = "wrap", add=T, col=c("#DC863B90", "#899DA490", "grey"), pch=20)
graphics.off()
system(paste("open", pdfname))


wilcox.test(AVIS[host], AVIS[MERGIS$Gene.names == "UL32"])


MERGIS$viral_exclUL32 <- viral & !MERGIS$Gene.names == "UL32"
MERGIS$host <- host
MERGIS$UL32 <- MERGIS$Gene.names == "UL32"
MERGIS$log2_average_phosphosite_fold_change <- AVIS

MERGIS2 <- MERGIS[!is.na(AVIS),]
