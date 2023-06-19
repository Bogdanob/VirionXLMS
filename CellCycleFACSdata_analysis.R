
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")

library(ggpubr)

##### This is making the analysis of FACS dotplots #####


MOCK <- cbind.data.frame(mock_G1_IE_minus <- c(61,60,62),
mock_G1_IE_plus <- c(0,0,0),
mock_S_IE_minus <- c(39,40,38),
mock_S_IE_plus <- c(0,0,0))

WT <- cbind.data.frame(WT_G1_IE_minus <- c(11,8,9),
WT_G1_IE_plus <- c(47,51,46),
WT_S_IE_minus <- c(40,38,42),
WT_S_IE_plus <- c(2,3,3))


SILK <- cbind.data.frame(SILK_G1_IE_minus <- c(13,12,15),
SILK_G1_IE_plus <- c(45,43,43),
SILK_S_IE_minus <- c(39,42,40),
SILK_S_IE_plus <- c(3,3,2))


RVXF <- cbind.data.frame(RVXF_G1_IE_minus <- c(13,11,13),
RVXF_G1_IE_plus <- c(44,46,44),
RVXF_S_IE_minus <- c(40,41,40),
RVXF_S_IE_plus <- c(3,2,3))



BOTH <- cbind.data.frame(BOTH_G1_IE_minus <- c(18,16,21),
BOTH_G1_IE_plus <- c(40,40,34),
BOTH_S_IE_minus <- c(40,41,43),
BOTH_S_IE_plus <- c(2,3,2))



RXL <- cbind.data.frame(RXL_G1_IE_minus <- c(5,2,7),
RXL_G1_IE_plus <- c(52,55,51),
RXL_S_IE_minus <- c(9,7,10),
RXL_S_IE_plus <- c(34,36,32))


BOTHRXL <- cbind.data.frame(BOTHRXL_G1_IE_minus <- c(9,7,14),
BOTHRXL_G1_IE_plus <- c(47,49,42),
BOTHRXL_S_IE_minus <- c(26,25,28),
BOTHRXL_S_IE_plus <- c(18,19,16))


DATS <- list((WT[,1:2]/rowSums(WT[,1:2]))[,2], (WT[,3:4]/rowSums(WT[,3:4]))[,2], 
          (RVXF[,1:2]/rowSums(RVXF[,1:2]))[,2], (RVXF[,3:4]/rowSums(RVXF[,3:4]))[,2],
          (SILK[,1:2]/rowSums(SILK[,1:2]))[,2], (SILK[,3:4]/rowSums(SILK[,3:4]))[,2],
          (BOTH[,1:2]/rowSums(BOTH[,1:2]))[,2], (BOTH[,3:4]/rowSums(BOTH[,3:4]))[,2], 
          (RXL[,1:2]/rowSums(RXL[,1:2]))[,2], (RXL[,3:4]/rowSums(RXL[,3:4]))[,2],
          (BOTHRXL[,1:2]/rowSums(BOTHRXL[,1:2]))[,2], (BOTHRXL[,3:4]/rowSums(BOTHRXL[,3:4]))[,2])

DATS1 <- cbind.data.frame((WT[,1:2]/rowSums(WT[,1:2]))[,2], (WT[,3:4]/rowSums(WT[,3:4]))[,2], 
             (RVXF[,1:2]/rowSums(RVXF[,1:2]))[,2], (RVXF[,3:4]/rowSums(RVXF[,3:4]))[,2],
             (SILK[,1:2]/rowSums(SILK[,1:2]))[,2], (SILK[,3:4]/rowSums(SILK[,3:4]))[,2],
             (BOTH[,1:2]/rowSums(BOTH[,1:2]))[,2], (BOTH[,3:4]/rowSums(BOTH[,3:4]))[,2], 
             (RXL[,1:2]/rowSums(RXL[,1:2]))[,2], (RXL[,3:4]/rowSums(RXL[,3:4]))[,2],
             (BOTHRXL[,1:2]/rowSums(BOTHRXL[,1:2]))[,2], (BOTHRXL[,3:4]/rowSums(BOTHRXL[,3:4]))[,2])

MEANS <- unlist(lapply(DATS, mean))



DATS2 <- unlist(DATS)
DATS3 <- cbind.data.frame(DATS2*100, c(rep("WT",6), rep("SILKmut",6), rep("RVXFmut",6), rep("SILKmut+RVXFmut",6), rep("RXLmut",6), rep("RXLmut+SILKmut+RVXFmut",6)),
                 rep(c(rep("G1",3), rep("S/G2",3)),6))

names(DATS3) <- c("percentage_IE_positive", "genotype","cell_cycle_stage")


pdf(pdfname <-"WT_PP1mutsRXLmut.pdf", width=10, height=7)

#layout(matrix(ncol=3, nrow=1, c(1:3)))


ggbarplot(DATS3, x = "genotype", y = "percentage_IE_positive", 
          add = c("mean_se", "jitter"), position = position_dodge(0.8), color = "cell_cycle_stage", palette = "npg", 
          ggtheme=theme_pubclean(), cex=1.5, lab.pos = "in")

DATS4 <- DATS3
DATS4$percentage_IE_positive <- DATS4$percentage_IE_positive/MEANS[9]

ggbarplot(DATS4, x = "genotype", y = "percentage_IE_positive", 
          add = c("mean_se", "jitter"), position = position_dodge(0.8), color = "cell_cycle_stage", palette = "npg", 
          ggtheme=theme_pubclean(), cex=1.5, lab.pos = "in", remove = c("WT","SILKmut", "RVXFmut", "SILKmut+RVXFmut"), ylab = "IE positive cells (RXLmut G1 = 100)")


DATS5 <- DATS3
DATS5$percentage_IE_positive <- DATS5$percentage_IE_positive/MEANS[1]

ggbarplot(DATS5, x = "genotype", y = "percentage_IE_positive", 
          add = c("mean_se", "jitter"), position = position_dodge(0.8), color = "cell_cycle_stage", palette = "npg", 
          ggtheme=theme_pubclean(), cex=1.5, lab.pos = "in", remove = c("RXLmut+\nSILKmut+RVXFmut", "RXLmut"), ylab = "IE positive cells (WT G1 = 100)")

graphics.off()

system(paste("open", pdfname))
