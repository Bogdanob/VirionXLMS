#### this contains many different things

## cytoscape export
library(RColorBrewer)
library("viridis")           # Load
library(wesanderson)
library(gplots)


VirTable <- read.csv("N:/VirionProject/XLb4Gradient/5mM_DSSO/Manuscript/SupplementaryTables/SupplementaryTable1.csv")
VirTable <- GeneNameExtractXlinkX(VirTable = VirTable)

VirTable$gene_a[grepl('ABV71575.1',VirTable$protein_a)] <- "UL45"
VirTable$gene_b[grepl('ABV71575.1',VirTable$protein_b)] <- "UL45"

VirTable$protein_a_cln[VirTable$gene_a == "UL45"] <- ">lcl|EF999921.1_prot_ABV71575.1_78 [protein=UL45] [protein_id=ABV71575.1] [location=complement(90411..93131)] [gbkey=CDS]"
VirTable$protein_b_cln[VirTable$gene_b == "UL45"] <- ">lcl|EF999921.1_prot_ABV71575.1_78 [protein=UL45] [protein_id=ABV71575.1] [location=complement(90411..93131)] [gbkey=CDS]"

# virus-virus links
VirPosPos <- grepl(x=VirTable$protein_a, pattern="lcl") & grepl(x=VirTable$protein_b, pattern="lcl")

# viral proteins of interest
des1<-c('UL46','UL85', 'UL86', 'UL80', 'UL104', 'UL93','UL48A', 'UL48','UL32', 'UL47','UL88', 'UL82','UL83', 'UL25', 'UL35', 'UL71', 'UL24', 'UL97', 'UL69','UL122', 'UL99','UL45','UL94','UL100' ,'RL10','UL33','UL132','US27', 'UL119', 'UL41A', 'UL148D', 'UL55', 'UL128', 'UL131A') 


KnownCapsid <- c("UL86", "UL85", "UL46", "UL48A", "UL93", "UL77", "UL80", "UL104")
KnownTegument <- c("UL38",'UL83', 'UL82', 'UL32', "UL94", 'UL99', "UL47", "UL48", "UL88", "UL25", "UL45", "UL103", "UL71", "UL26", "UL24", "UL35", "UL97", "UL69", "TRS1", 'IRS1', 'UL43', 'UL96', "US23")
KnownGlyco <- c('UL148D','UL41A','UL33', "UL55", "UL100", "UL73", "RL10", "UL75", "UL132", "UL115", "UL119", "UL118", "UL74", "UL116", "UL128", "UL130", "UL131A", "US27")


#
VirionVirTable <- VirTable[(VirTable$gene_a %in% des1 & VirPosPos) | (VirTable$gene_b %in% des1 & VirPosPos),]


IntNamesViral <- apply(cbind(str_extract(string = VirionVirTable$gene_a, pattern="[^ ]{2,6}"), 
                             str_extract(string = VirionVirTable$gene_b, pattern="[^ ]{2,6}")), 1, function(x) paste(sort(x), collapse=" "))

IntNamesViral<- IntNamesViral[!str_extract(string = VirionVirTable$gene_a, pattern="[^ ]{2,6}")==str_extract(string = VirionVirTable$gene_b, pattern="[^ ]{2,6}")]


HM_Viral <-  names(table(unlist(str_split(IntNamesViral," "))))



HighConfidenceViral <- (VirTable$gene_a %in% HM_Viral & VirPosPos) | (VirTable$gene_b %in% HM_Viral & VirPosPos)

HighConfidenceViral <- unique((c(VirTable$gene_b[HighConfidenceViral], VirTable$gene_a[HighConfidenceViral])))[unique((c(VirTable$gene_b[HighConfidenceViral], VirTable$gene_a[HighConfidenceViral]))) %in% HM_Viral]


HostVirMatrix <- matrix(ncol=length(HighConfidenceViral), nrow=length(HighConfidenceViral), dimnames = list(HighConfidenceViral, HighConfidenceViral))


colnames(HostVirMatrix)<- colnames(HostVirMatrix)[order(match(colnames(HostVirMatrix),des1))]
row.names(HostVirMatrix)<- row.names(HostVirMatrix)[order(match(row.names(HostVirMatrix),des1))]




### Clustering PPI specificities ###

pdf(pdfname <- "Figure2ab.pdf", width=10, height=10)


toInput <- table(apply(cbind(str_extract(string = VirTable$gene_a, pattern="[^ ]{2,6}"), 
                             str_extract(string = VirTable$gene_b, pattern="[^ ]{2,6}")), 1, function(x) paste(sort(x), collapse=" ")))

for (i in 1:length(colnames(HostVirMatrix))){
  for (j in 1:length(row.names(HostVirMatrix))){
    tek <-  paste(sort(c(colnames(HostVirMatrix)[i],row.names(HostVirMatrix)[j])), collapse=" ") 
    if (any(names(toInput) == tek) == F){
      tek2 <- 0}
    if (any(names(toInput) == tek) == T){
      tek2 <- toInput[names(toInput) == tek]}
    HostVirMatrix[j,i] <- tek2
  }  
}

dig <- diag(HostVirMatrix)

#diag(HostVirMatrix)<-0

### remove entires with fewer than 10 links
HostVirMatrix <- HostVirMatrix[!rowSums(HostVirMatrix) < 9,]
HostVirMatrix <- HostVirMatrix[,!colSums(HostVirMatrix) < 9]
HostVirMatrix <- HostVirMatrix[!rowSums(HostVirMatrix) < 9,]
HostVirMatrix <- HostVirMatrix[,!colSums(HostVirMatrix) < 9]



JaccMat <- HostVirMatrix
for (i in 1:dim(HostVirMatrix)[1]){
  for (j in 1:dim(HostVirMatrix)[2]){
  JaccMat[i,j]<- (HostVirMatrix[i,j])/(colSums(HostVirMatrix)[j])
  }
}

mywesCols <- wes_palette(name="Royal1", n=4, type="discrete")[c(2,1,3,4)]

JaccMat <- JaccMat[-32,-32]

col_breaks = c(seq(0,0.5,length=40))
my_x_cols <- colnames(JaccMat) 
my_x_cols[colnames(JaccMat) %in% KnownTegument] <- mywesCols[2]
my_x_cols[colnames(JaccMat) %in% KnownCapsid] <- mywesCols[1]
my_x_cols[colnames(JaccMat) %in% c(KnownGlyco, "UL41A", "UL148D", 'US27')] <- mywesCols[3]
my_x_cols[!(colnames(JaccMat) %in% c(KnownGlyco, "UL41A", "UL148D", 'US27')) & !(colnames(JaccMat) %in% KnownCapsid) & !(colnames(JaccMat) %in% KnownTegument)] <- "grey"

#install.packages("viridis")  # Install


my_y_cols <- my_x_cols
my_palette <- colorRampPalette(rev(viridis(39)))


#A <- wes_palette("Zissou1", 3, type = c("continuous"))

#wes_palette(name="Zissou1", n=4, type="discrete")

h1 <- heatmap.2(JaccMat, scale = "none", trace = "none", 
                breaks = col_breaks, col=my_palette, Rowv = T,Colv = F, key = T, sepcolor = "white", 
                sepwidth=c(0.005,0.005), ColSideColors = my_x_cols, RowSideColors = my_y_cols,reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean),
              distfun=function(x) as.dist(1-cor(t(x))),
                hclustfun=function(x) hclust(x, method="complete"), rowsep=c(7,13,20,31), colsep=c(7,23), key.title = "", key.ylab = "", 
              key.xlab = "PPI\nspecificity", density.info = 'none', xlab="interactor 1", ylab="interactor 2")

write.csv(JaccMat, file="sourcedataFig2a.csv")




#### Scaffold Indexes withi cluster ####


diag(JaccMat) <- 0

cl1_suba <- c("UL32", "UL48A", "UL48", "UL47", "UL82", "UL35", "UL88")
cl1_subb <- c("UL86", "UL80", "UL104", "UL46", "UL93", "UL85")
cl2_suba <- c("UL94", "UL83", "UL25", "UL24", "UL122", "UL97", "UL69")
cl2_subb <- c("UL45", "UL132", "UL33", "UL99", "UL100", "RL10", "UL55", "US27", "UL41A", "UL71", "UL119")



pdf(pdfname <- "Figure2B.pdf", width=6, height=12)
layout(matrix(c(1:4),ncol=1, nrow=4))



MyLayerScaffold <- function(cl1_suba, lim=2){
  
  y <- rowSums(JaccMat[,colnames(JaccMat) %in% cl1_suba])[h1$rowInd]
  
  x<- barplot(y, frame=F, ylab="Scaffold index = sum PPI specificity", axes=F, xlab="", col=my_x_cols[h1$rowInd], pch=20, las=2, xaxt="n", ylim=c(0,3))
  labs <- paste(names(rowSums(JaccMat)[h1$rowInd]))
  text(cex=1, x=x-.5, y=-.25, labs, xpd=TRUE, srt=90)
  
  #axis(1, at=c(1:dim(JaccMat)[1]), labels=names((rowSums(JaccMat))[h1$rowInd]), las=2)
  axis(2, at=c(0:3), labels=c(0:3), las=1)
  #abline(h=2, lty=2, lwd=2)
  grid(NA,NULL)
  
  return(y)
}


v1<-MyLayerScaffold(cl1_suba = cl1_suba, lim=2)
v2<-MyLayerScaffold(cl1_suba = cl1_subb, lim=1)
v3<-MyLayerScaffold(cl1_suba = cl2_suba, lim=1)
v4<-MyLayerScaffold(cl1_suba = cl2_subb, lim=2)






### this is doing the Lysine based calculations ###



VirTable <- read.csv("N:/VirionProject/XLb4Gradient/5mM_DSSO/Manuscript/SupplementaryTables/SupplementaryTable1.csv")
VirTable <- GeneNameExtractXlinkX(VirTable = VirTable)
#write.csv(file = "N:/VirionProject/XLb4Gradient/5mM_DSSO/Manuscript/SupplementaryTables/SupplementaryTable1.csv", x=VirTable)

VirTable$gene_a[grepl('ABV71575.1',VirTable$protein_a)] <- "UL45"
VirTable$gene_b[grepl('ABV71575.1',VirTable$protein_b)] <- "UL45"

VirTable$protein_a_cln[VirTable$gene_a == "UL45"] <- ">lcl|EF999921.1_prot_ABV71575.1_78 [protein=UL45] [protein_id=ABV71575.1] [location=complement(90411..93131)] [gbkey=CDS]"
VirTable$protein_b_cln[VirTable$gene_b == "UL45"] <- ">lcl|EF999921.1_prot_ABV71575.1_78 [protein=UL45] [protein_id=ABV71575.1] [location=complement(90411..93131)] [gbkey=CDS]"


VirPosPos <- grepl(x=VirTable$protein_a, pattern="lcl") & grepl(x=VirTable$protein_b, pattern="lcl")



VirionVirTable <- VirTable[(VirTable$gene_a %in% HM_Viral & VirPosPos) | (VirTable$gene_b %in% HM_Viral & VirPosPos),]



A<- sort(paste(VirionVirTable$gene_a,VirionVirTable$LinkPos1,VirionVirTable$gene_b, VirionVirTable$LinkPos2))


A2 <- cbind.data.frame(VirionVirTable$gene_a,VirionVirTable$LinkPos1,VirionVirTable$gene_b, VirionVirTable$LinkPos2)

A2 <- subset(A2, A2$`VirionVirTable$gene_a` != A2$`VirionVirTable$gene_b`)

A2 <- subset(A2, A2$`VirionVirTable$gene_a` == "UL32" |  A2$`VirionVirTable$gene_b` == "UL32") ### only include links to protein of interest

A2 <- subset(A2, !(A2$`VirionVirTable$gene_a` %in% c("UL32") &  A2$`VirionVirTable$gene_b` %in% c("UL32"))) ### remove intralinks

ClusterVM <- c("UL100", "RL10", "UL33", "UL132", "US27", "UL119", "UL41A", "UL148D", "UL131A", "UL55", "UL128", "UL73") ## "UL99", "UL45", "UL71"
ClusterOT <- c("UL24", "UL97", "UL69", "UL83", "UL94", "UL25", "UL122")
ClusterIT <- c("UL47", "UL48","UL82", "UL88", "UL35") ### remove the protein of interest here UL32
ClusterNC <- c("UL46", "UL86", "UL93", "UL85", "UL80", "UL104", "UL77", "UL48A")


A2$`VirionVirTable$gene_a`[A2$`VirionVirTable$gene_a` %in% ClusterVM] <- "VM"
A2$`VirionVirTable$gene_b`[A2$`VirionVirTable$gene_b` %in% ClusterVM] <- "VM"

A2$`VirionVirTable$gene_a`[A2$`VirionVirTable$gene_a` %in% ClusterOT] <- "OT"
A2$`VirionVirTable$gene_b`[A2$`VirionVirTable$gene_b` %in% ClusterOT] <- "OT"

A2$`VirionVirTable$gene_a`[A2$`VirionVirTable$gene_a` %in% ClusterIT] <- "IT"
A2$`VirionVirTable$gene_b`[A2$`VirionVirTable$gene_b` %in% ClusterIT] <- "IT"

A2$`VirionVirTable$gene_a`[A2$`VirionVirTable$gene_a` %in% ClusterNC] <- "NC"
A2$`VirionVirTable$gene_b`[A2$`VirionVirTable$gene_b` %in% ClusterNC] <- "NC"

A2.1 <- A2[A2$`VirionVirTable$gene_a` != "UL32",c(3,4,1,2)] 
A2.2 <- A2[A2$`VirionVirTable$gene_a` == "UL32",c(1,2,3,4)]

colnames(A2.1) <- c("gene_a", "LinkPos1", 'gene_b', 'LinkPos2')
colnames(A2.2) <- c("gene_a", "LinkPos1", 'gene_b', 'LinkPos2')

A2 <- rbind(A2.1,A2.2)

A3 <- A2[,c(2,3)]

Aggri <- aggregate(A3, by=list(A3$LinkPos1), table)$gene_b
Aggri_Pos <- aggregate(A3, by=list(A3$LinkPos1), table)$Group.1

ResiVM <- Aggri
ResiNC <- Aggri
ResiOT <- Aggri
ResiIT <- Aggri

for (i in 1:length(Aggri)){
  ResiVM[[i]]<- as.numeric(Aggri[[i]][names(Aggri[[i]])=="VM"]/sum(Aggri[[i]]))
  ResiNC[[i]]<- as.numeric(Aggri[[i]][names(Aggri[[i]])=="NC"]/sum(Aggri[[i]]))
  ResiOT[[i]]<- as.numeric(Aggri[[i]][names(Aggri[[i]])=="OT"]/sum(Aggri[[i]]))
  ResiIT[[i]]<- as.numeric(Aggri[[i]][names(Aggri[[i]])=="IT"]/sum(Aggri[[i]]))
  }

ResiVM[lengths(ResiVM) == 0] <- 0
ResiIT[lengths(ResiIT) == 0] <- 0
ResiOT[lengths(ResiOT) == 0] <- 0
ResiNC[lengths(ResiNC) == 0] <- 0


#pdf(pdfname <- "pp150_connectivity.pdf", width=8, height=8)
#layout(matrix(ncol=1, nrow=2, c(1:2)))

mycols <- wes_palette(name="Royal1", n=4, type="discrete")[c(2,1,3,4)]

###  "#C93312" "#899DA4" "#FAEFD1" "#DC863B"
par(mar=c(5, 4, 4, 6) + 0.1)


layout(matrix(ncol=1, nrow=2, c(1:2)))

plot( x= Aggri_Pos, (unlist(ResiIT)),axes = F, pch=20,xlab='sequence position UL32', 
        ylab="fraction of cross-links (Inter-links)", main="UL32/pp150 to inner tegument proteins",
        frame = F, col="#899DA4", cex=2)
axis(2, ylim=c(0,1),las=1)  ## las=1 makes horizontal labels
axis(1, xlim=c(1,1049))

plot( x= Aggri_Pos, unlist(ResiVM),axes = F, pch=20,xlab='sequence position UL32', 
      ylab="fraction of cross-links (Inter-links)", main="UL32/pp150 to viral membrane proteins", 
      ylim=c(0,1),frame = F, col="#F5DE91", cex=2)
axis(2, ylim=c(0,1),las=1)  ## las=1 makes horizontal labels
axis(1, xlim=c(1,1049))




write.csv(rbind(Aggri_Pos, unlist(ResiIT), unlist(ResiVM)), file="sourceDataFig2cd.csv") 


wilcox.test(unlist(ResiVM)[Aggri_Pos>774],unlist(ResiVM)[Aggri_Pos<774]) 
wilcox.test(unlist(ResiIT)[Aggri_Pos>340],unlist(ResiIT)[Aggri_Pos<340]) 

