

#### High MOI growth PP1 mutant + single mutants

WT_d0 <- c(5.88E+04,8.14E+03,4.49E+04) 
WT_d2 <- c(3.28E+02,8.29E+01,3.41E+02)
WT_d4 <- c(4.80E+05,4.00E+05,1.04E+06)
WT_d6 <- c(2.80E+07,3.24E+07,5.35E+07)
WT_d8 <- c(2.44E+07,2.32E+07,3.28E+07)

RV522_d0 <- c(1.69E+04,4.05E+03,3.57E+04)
RV522_d2 <- c(2.55E+02,3.99E+01,1.15E+02)
RV522_d4 <- c(1.60E+05,1.60E+05,4.37E+05)
RV522_d6 <- c(1.52E+07,1.48E+07,2.59E+07)
RV522_d8 <- c(1.72E+07,1.76E+07,2.47E+07)

RV523_d0 <- c(4.00E+04,4.18E+03,7.94E+04)
RV523_d2 <- c(2.07E+02,8.36E+01,3.45E+02)
RV523_d4 <- c(4.00E+05,2.40E+05,6.67E+05)
RV523_d6 <- c(1.32E+07,1.48E+07,3.05E+07)
RV523_d8 <- c(1.64E+07,1.48E+07,2.30E+07)

RV524_d0 <- c(3.82E+04,	4.00E+03,	4.95E+04)
RV524_d2 <- c(8.59E+01,8.00E+01,3.45E+02)
RV524_d4 <- c(1.12E+05,9.60E+04,3.54E+05)
RV524_d6 <- c(6.4E+06,6.8E+06,1.15E+07)
RV524_d8 <- c(6.40E+06,6.4E+06,1.73E+07)


WT <- log10(rbind.data.frame(WT_d0, WT_d2, WT_d4, WT_d6, WT_d8))
RV522 <- log10(rbind.data.frame(RV522_d0, RV522_d2, RV522_d4, RV522_d6, RV522_d8))
RV523 <- log10(rbind.data.frame(RV523_d0, RV523_d2, RV523_d4, RV523_d6, RV523_d8))
RV524 <- log10(rbind.data.frame(RV524_d0, RV524_d2, RV524_d4, RV524_d6, RV524_d8))

colnames(WT) <- c("rep1", "rep2", "rep3")
row.names(WT) <- c("d0", "d2", "d4", "d6", "d8") 

colnames(RV522) <- c("rep1", "rep2", "rep3")
row.names(RV522) <-c("d0", "d2", "d4", "d6", "d8")  

colnames(RV523) <- c("rep1", "rep2", "rep3")
row.names(RV523) <- c("d0", "d2", "d4", "d6", "d8")  

colnames(RV524) <- c("rep1", "rep2", "rep3")
row.names(RV524) <- c("d0", "d2", "d4", "d6", "d8")  

WT[,4]<- apply(WT,1, FUN=sd, na.rm=T)
RV522[,4] <- apply(RV522,1, FUN=sd, na.rm=T)
RV523[,4]<- apply(RV523,1, FUN=sd, na.rm=T)
RV524[,4]<- apply(RV524,1, FUN=sd, na.rm=T)

WT[,5] <- rowMeans(WT[,1:3], na.rm = T)
RV522[,5] <- rowMeans(RV522[,1:3], na.rm=T)
RV523[,5] <- rowMeans(RV523[,1:3], na.rm=T)
RV524[,5] <- rowMeans(RV524[,1:3], na.rm=T)


#pdf(pdfname <- "HighMOI_growth.pdf", width=8, height=8)


plot(y=WT$V5, x=c(0,2,4,6,8), axes=F, type="l", ylab="log10 [IU / mL]", xlab="days post infection", col="orange", lwd=1, ylim=c(1,8))
axis(2, at=NULL)
axis(1, at=c(0,2,4,6,8))
grid(NULL, NULL)
lines(y=RV522$V5, x=c(0,2,4,6,8), type="l", col="mediumseagreen", lwd=1)
lines(y=RV523$V5, x=c(0,2,4,6,8), type="l", col="steelblue2", lwd=1)
lines(y=RV524$V5, x=c(0,2,4,6,8), type="l", col="red", lwd=1)
abline(h=log10(3.33E+01), col="darkgrey", lty=3, lwd=1)
text(x=5, y=2.05, "detection limit", cex=1.5)


arrows(x0=c(0,2,4,6,8), y0=RV522$V5,  x1=c(0,2,4,6,8), y1=RV522$V5-RV522$V4, angle = 90,length = 0.05, col="mediumseagreen", lwd=2)
arrows(x0=c(0,2,4,6,8), y0=RV522$V5,  x1=c(0,2,4,6,8), y1=RV522$V5+RV522$V4, angle = 90,length = 0.05, col="mediumseagreen", lwd=2)

arrows(x0=c(0,2,4,6,8), y0=RV523$V5,  x1=c(0,2,4,6,8), y1=RV523$V5-RV523$V4, angle = 90,length = 0.05, col="steelblue2", lwd=2)
arrows(x0=c(0,2,4,6,8), y0=RV523$V5,  x1=c(0,2,4,6,8), y1=RV523$V5+RV523$V4, angle = 90,length = 0.05, col="steelblue2", lwd=2)

arrows(x0=c(0,2,4,6,8), y0=WT$V5,  x1=c(0,2,4,6,8), y1=WT$V5-WT$V4, angle = 90,length = 0.05, col="orange", lwd=0.5)
arrows(x0=c(0,2,4,6,8), y0=WT$V5,  x1=c(0,2,4,6,8), y1=WT$V5+WT$V4, angle = 90,length = 0.05, col="orange", lwd=0.5)

arrows(x0=c(0,2,4,6,8), y0=RV524$V5,  x1=c(0,2,4,6,8), y1=RV524$V5-RV524$V4, angle = 90,length = 0.05, col="red", lwd=0.5)
arrows(x0=c(0,2,4,6,8), y0=RV524$V5,  x1=c(0,2,4,6,8), y1=RV524$V5+RV524$V4, angle = 90,length = 0.05, col="red", lwd=0.5)

legend("topleft", c("WT","SILKmut", "RVxFmut" ,"SILK+RVxFmut"), col=c("orange", "mediumseagreen", "steelblue2","red"), bty="n", cex=1.5, lwd=3)

points(x=jitter(rep(0,3), amount = 0.15),y=WT[1,1:3], col="orange", pch=20, cex=1.75)
points(x=jitter(rep(2,3), amount = 0.15),y=WT[2,1:3], col="orange", pch=20, cex=1.75)
points(x=jitter(rep(4,3), amount = 0.15),y=WT[3,1:3], col="orange", pch=20, cex=1.75)
points(x=jitter(rep(6,3), amount = 0.15),y=WT[4,1:3], col="orange", pch=20, cex=1.75)
points(x=jitter(rep(8,3), amount = 0.15),y=WT[5,1:3], col="orange", pch=20, cex=1.75)


points(x=jitter(rep(0,3), amount = 0.15),y=RV522[1,1:3], col="mediumseagreen", pch=20, cex=1.75)
points(x=jitter(rep(2,3), amount = 0.15),y=RV522[2,1:3], col="mediumseagreen", pch=20, cex=1.75)
points(x=jitter(rep(4,3), amount = 0.15),y=RV522[3,1:3], col="mediumseagreen", pch=20, cex=1.75)
points(x=jitter(rep(6,3), amount = 0.15),y=RV522[4,1:3], col="mediumseagreen", pch=20, cex=1.75)
points(x=jitter(rep(8,3), amount = 0.15),y=RV522[5,1:3], col="mediumseagreen", pch=20, cex=1.75)


points(x=jitter(rep(0,3), amount = 0.15),y=RV523[1,1:3], col="steelblue2", pch=20, cex=1.75)
points(x=jitter(rep(2,3), amount = 0.15),y=RV523[2,1:3], col="steelblue2", pch=20, cex=1.75)
points(x=jitter(rep(4,3), amount = 0.15),y=RV523[3,1:3], col="steelblue2", pch=20, cex=1.75)
points(x=jitter(rep(6,3), amount = 0.15),y=RV523[4,1:3], col="steelblue2", pch=20, cex=1.75)
points(x=jitter(rep(8,3), amount = 0.15),y=RV523[5,1:3], col="steelblue2", pch=20, cex=1.75)

points(x=jitter(rep(0,3), amount = 0.15),y=RV524[1,1:3], col="red", pch=20, cex=1.75)
points(x=jitter(rep(2,3), amount = 0.15),y=RV524[2,1:3], col="red", pch=20, cex=1.75)
points(x=jitter(rep(4,3), amount = 0.15),y=RV524[3,1:3], col="red", pch=20, cex=1.75)
points(x=jitter(rep(6,3), amount = 0.15),y=RV524[4,1:3], col="red", pch=20, cex=1.75)
points(x=jitter(rep(8,3), amount = 0.15),y=RV524[5,1:3], col="red", pch=20, cex=1.75)





#graphics.off()

#system(paste("open", pdfname))

#t.test(x = as.numeric(WT[3,1:3]), y = as.numeric(RV524[3,1:3]), exact=F, paired = F)$p.value
t.test(x = as.numeric(WT[4,1:3]), y = as.numeric(RV524[4,1:3]), exact=F, pared = F)$p.value
#t.test(x = as.numeric(WT[5,1:3]), y = as.numeric(RV524[5,1:3]), exact=F, paired = F)$p.value

text("P = 0.005", x=6, y=8)
#text("P = 0.0013 \n (WT - SILK+RVxFmut)", x=12, y=7.55)
#text("P = 0.00017 (WT - SILK+RVxFmut)", x=16, y=7.85)








