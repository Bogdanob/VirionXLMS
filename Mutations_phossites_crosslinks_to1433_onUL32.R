library(stringr)
library(wesanderson)


setwd("N:/VirionProject/NewDataset/DataInput")


## cross-links
DataS <- read.csv(file="N:/VirionProject/XLb4Gradient/5mM_DSSO/Manuscript/SupplementaryTables/SupplementaryTable2.csv", stringsAsFactors = F)

## predicted 14-3-3 binding
DataP <- read.csv(file="N:/VirionProject/NewDataset/14-3-3/1433pred_s9634275465.csv", stringsAsFactors = F)

## phosphosites
DataPhos <- read.csv(file="N:/VirionProject/PP1/txt_PrideUpload/Phospho (STY)Sites.txt", stringsAsFactors = F, sep="\t")


DataPhos <- subset(x=DataPhos, subset = DataPhos$Localization.prob.WTH_524L_phos > 0.75 & DataPhos$Localization.prob.WTL_524H_phos > 0.75)

DataPhos <- subset(x=DataPhos, subset = grepl(x = DataPhos$Fasta.headers,pattern =  "UL32"))

LISTG <- c("UL32","YWHAZ","YWHAH","YWHAG","YWHAE","YWHAQ","YWHAB")
Example <- DataS[DataS$gene_a %in% LISTG & DataS$gene_b %in% LISTG,]

Example_cln <- Example[!Example$gene_a == Example$gene_b,]

Example_cln <- Example_cln[!(grepl(Example_cln$gene_a, pattern="YWH") & grepl(Example_cln$gene_b, pattern="YWH")),]



Residue <- c(1:1049)
R1 <- cbind.data.frame(Residue,rep(NA, length(Residue)))

Pred <- (cbind.data.frame(DataP$Site, DataP$Consensus))

names(Pred)[1] <- "Residue"

A <- merge(Pred,R1, all.y = T)[,-3]

LINKS <- as.data.frame(table(Example_cln$LinkPos1))
names(LINKS)[1] <- "Residue"
Ph <- unique(unlist(str_split(DataPhos$Positions.within.proteins, pattern = ";")))
Ph1 <- cbind.data.frame(Ph, rep(1,length(Ph)))

SET <- merge(A,LINKS, all.x = T)

names(Ph1)[1] <- "Residue"

SET2 <- merge(SET,Ph1,by = "Residue",all.x=T)[,-5]
#SET2$`DataP$Consensus` <- scale(SET2$`DataP$Consensus`)
#SET2$Freq <- scale(SET2$Freq)
SET2$`rep(1, length(Ph))`[!is.na(SET2$`rep(1, length(Ph))`)] <- 2

pdf(pdfname <- "14-3-3Mutations.pdf", width=20, height=6)
layout(matrix(c(1,2,3), ncol=1, nrow=3))

sels <- SET2[,2]>0.6
sels[sels==F] <- NA

image(as.matrix(SET2[sels,2]), col = "tomato", axes=F, main="14-3-3 binding prediction", )
axis(1,at = c(0,(c(1:10)*100)/1049,1), labels = c(c(0:10)*100,1049))
polygon(x = c(982/1049,982/1049,993/1049,993/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(931/1049,931/1049,942/1049,942/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(1006/1049,1006/1049,1017/1049,1017/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(915/1049,915/1049,926/1049,926/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(862/1049,862/1049,874/1049,874/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(513/1049,513/1049,524/1049,524/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(473/1049,473/1049,484/1049,484/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(494/1049,494/1049,505/1049,505/1049),y=c(-1,1,1,-1), col = my_col)


image(as.matrix(SET2[,3]), col = c("tomato"), axes=F, main="cross-links")
axis(1,at = c(0,(c(1:10)*100)/1049,1), labels = c(c(0:10)*100,1049))
polygon(x = c(982/1049,982/1049,993/1049,993/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(931/1049,931/1049,942/1049,942/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(1006/1049,1006/1049,1017/1049,1017/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(915/1049,915/1049,926/1049,926/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(862/1049,862/1049,874/1049,874/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(513/1049,513/1049,524/1049,524/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(473/1049,473/1049,484/1049,484/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(494/1049,494/1049,505/1049,505/1049),y=c(-1,1,1,-1), col = my_col)

image(as.matrix(SET2[,4]), col = c("mediumseagreen"), axes=F, main="Phosphosites")
axis(1,at = c(0,(c(1:10)*100)/1049,1), labels = c(c(0:10)*100,1049))
polygon(x = c(982/1049,982/1049,993/1049,993/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(931/1049,931/1049,942/1049,942/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(1006/1049,1006/1049,1017/1049,1017/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(915/1049,915/1049,926/1049,926/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(862/1049,862/1049,874/1049,874/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(513/1049,513/1049,524/1049,524/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(473/1049,473/1049,484/1049,484/1049),y=c(-1,1,1,-1), col = my_col)
polygon(x = c(494/1049,494/1049,505/1049,505/1049),y=c(-1,1,1,-1), col = my_col)

graphics.off()
system(paste("open", pdfname))


