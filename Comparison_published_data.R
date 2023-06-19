library(openxlsx)
library(VennDiagram)

S2 <- read.xlsx("N:/VirionProject/XLb4Gradient/5mM_DSSO/CouteS2/Table S2_mod.xlsx")

Coute  <- na.omit(S2$Gene.name[S2$Couté.et.al.== 1])
Rieder <- na.omit(S2$Gene.name[S2$Rieder.et.al.== 1])
Reyda  <- na.omit(S2$Gene.name[S2$Reyda.et.al.== 1])
Varnum <- na.omit(S2$Gene.name[S2$Varnum.et.al.== 1])


rownames(HostMatrix)[!rownames(HostMatrix) %in% Coute]
rownames(HostMatrix)[!rownames(HostMatrix) %in% Rieder]
rownames(HostMatrix)[!rownames(HostMatrix) %in% Reyda]
rownames(HostMatrix)[!rownames(HostMatrix) %in% Varnum]



drawMyVenn <- function(Coute, title="Coute et al."){
  A <- row.names(HostMatrix)
  B <- Coute
  
  MyVenn1 <-VennDiagram::draw.pairwise.venn(length(A), length(B), cross.area = length(intersect(A,B)), col = brewer.pal(3, name = "Set1")[2:3], 
                                            category = c("XL-MS PPI map", paste("associated with particles \n",title, sep="")), alpha=0.25, 
                                                         fill=brewer.pal(3, name = "Set1")[2:3], cex = 3, cat.cex = 3, cat.dist = c(-0.02,-0.02))
                                            
}

pdf(pdfname <- "Comparison_published.pdf",width=12,height=12)

drawMyVenn(Coute = Coute)
plot.new()
drawMyVenn(Coute = Rieder, title = "Rieder et al.")
plot.new()
drawMyVenn(Coute = Reyda, title = "Reyda et al.")
plot.new()
drawMyVenn(Coute = Varnum, title = "Varnum et al.")
plot.new()
drawMyVenn(Coute = S2$Gene.name, title = "Coute et al.\nRieder et al.\nReyda et al.\nVarnum et al.")



graphics.off()

system(paste("open", pdfname))