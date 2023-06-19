
library(ggpubr)
library(ggplot2)

### extended data figure ### wb quant ###

## each normalized to mcp

N1433_WT <- c(0.550257724, 1.429410194, 1.627201343, 1.262834194)
Npp150_WT <- c(1.972696404, 1.310977553, 2.361373607, 2.484556395)

N1433_3 <- c(0.208013039, 0.288012279,0.892381934,0.976169156)
Npp150_3 <- c(2.022926225, 1.416214776, 2.207979524, 1.970932351)

N1433_34 <- c(0.051744733, 0.020611878, 0.141852268, 0.006457134)
Npp150_34 <- c(1.80675717,1.465637636, 2.127054607, 1.846731456)

N1433_45 <- c(0.00443956, -0.001769838, 0.032209527, 0.009999852)
Npp150_45 <- c(2.295645893, 1.586939099, 2.529915283, 1.946321436)



df2<- as.data.frame(matrix(ncol=1, nrow=16,c(N1433_WT, N1433_3, N1433_34, N1433_45)))

df2[c(1,5,9,13),1] <- df2[c(1,5,9,13),1]/df2[1,1]
df2[c(2,6,10,14),1] <- df2[c(2,6,10,14),1]/df2[2,1]
df2[c(3,7,11,15),1] <- df2[c(3,7,11,15),1]/df2[3,1]
df2[c(4,8,12,16),1] <- df2[c(4,8,12,16),1]/df2[4,1]

df2[,2] <- rep("14-3-3 (pan)", 16)
df2[,3] <- c(rep("WT",4), rep("3mut",4), rep("3/4mut",4), rep("4/5mut",4))

df3<- as.data.frame(matrix(ncol=1, nrow=16,c(Npp150_WT, Npp150_3, Npp150_34, Npp150_45)))

df3[c(1,5,9,13),1] <- df3[c(1,5,9,13),1]/df3[1,1]
df3[c(2,6,10,14),1] <- df3[c(2,6,10,14),1]/df3[2,1]
df3[c(3,7,11,15),1] <- df3[c(3,7,11,15),1]/df3[3,1]
df3[c(4,8,12,16),1] <- df3[c(4,8,12,16),1]/df3[4,1]

df3[,2] <- rep("UL32 (pp150)", 16)
df3[,3] <- c(rep("WT",4), rep("3mut",4), rep("3/4mut",4), rep("4/5mut",4))

df2 <- rbind(df2,df3)


names(df2) <- c("normalized_levels", "protein", "genotype")

df2$normalized_levels[df2$normalized_levels < 0] <- 0


pdf(pdfname <- "WB_quant_1433muts.pdf", width=5, height=5)

ggbarplot(df2, x = "genotype", y = "normalized_levels", 
          add = c("mean_sd", "jitter"), position = position_dodge(0.8), color = "protein", palette = "npg", 
          ggtheme=theme_pubclean(), cex=1.5, lab.pos = "in")

graphics.off()
system(paste("open", pdfname))


### this is making quan for RXL vs RXL/SILK/RVXF ###

N1433_RXL <- c(0.316072809,0.13919807, 0.334770186, 0.233706472)
N1433_RXLPP1 <- c(0.788871741,0.271991919,0.864933449, 0.388749219)
NPP1_RXL <- c(0.532279332, 0.372208131, 0.24532924)
NPP1_RXLPP1 <- c(0.238947044, 0.168459874, 0.132183913)
NPP150_RXL <- c(1.881043124, 2.050458241, 2.135458335)
NPP150_RXLPP1 <- c(1.523312158, 1.972809001, 1.645898377)

df4<- as.data.frame(matrix(ncol=1, nrow=20,c(N1433_RXL, N1433_RXLPP1, NPP1_RXL, NPP1_RXLPP1, NPP150_RXL, NPP150_RXLPP1)))

df4<- as.data.frame(df4[c(1:20),]/df4[c(1,2,3,4,1,2,3,4,9,10,11,9,10,11,15,16,17,15,16,17),])

df4$"genotype" <- c(rep("RXLmut",4), rep("RXL/SILK/RVXFmut",4), rep("RXLmut",3), rep("RXL/SILK/RVXFmut",3), rep("RXLmut",3), rep("RXL/SILK/RVXFmut",3))

df4$"protein" <- c(rep("14-3-3(pan)",8), rep("PP-1",6), rep("UL32 (pp150)",6))
names(df4)[1] <- "normalized_levels"

pdf(pdfname <- "WB_quant_PP1_RXLmuts.pdf", width=5, height=5)
ggbarplot(df4, x = "genotype", y = "normalized_levels", 
          add = c("mean_sd", "jitter"), position = position_dodge(0.8), color = "protein", palette = "npg", 
          ggtheme=theme_pubclean(), cex=1.5, lab.pos = "in")

graphics.off()
system(paste("open", pdfname))

t.test(c(0.4489,0.4525,0.5488), c(0.81,0.96,0.77))

t.test(c(2.5,1.95,2.58,1.66), c(0.81,0.96,0.77))
