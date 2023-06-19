



#### This is for growth curve low MOI PP1-binding deficient viral mutant

WT_d0 <- c(1217.73,	375.18,	129.90) 
WT_d4 <- c(180.95,1320.00	,392.66)
WT_d8 <- c(14880.00,	2640.00,	7920.00)
WT_d12 <- c(16144800.00,	6600000.00,	8976000.00)
WT_d16 <- c(53940000.00,	44220000.00,	41580000.00)
WT_d20 <- c(9300000.00	,8448000.00,	9108000.00)

RV524_d0 <- c(915.99,	134.05,	1320.00)
RV524_d4 <- c(361.87,	131.11,	124.68)
RV524_d8 <- c(540.91,	403.83,	519.23)
RV524_d12 <- c(595200.00,	348480.00,396000.00)
RV524_d16 <- c(1860000.00,	2904000.00,2508000.00)
RV524_d20 <- c(8370000.00,	7260000.00,	8184000.0)


WT <- log10(rbind.data.frame(WT_d0, WT_d4, WT_d8, WT_d12, WT_d16, WT_d20))
RV524 <- log10(rbind.data.frame(RV524_d0, RV524_d4, RV524_d8, RV524_d12, RV524_d16, RV524_d20))

colnames(WT) <- c("rep1", "rep2", "rep3")
row.names(WT) <- c("d0", "d4", "d8", "d12", "d16", "d20") 

colnames(RV524) <- c("rep1", "rep2", "rep3")
row.names(RV524) <- c("d0", "d4", "d8", "d12", "d16", "d20") 

WT[,4]<- apply(WT,1, FUN=sd, na.rm=T)
RV524[,4]<- apply(RV524,1, FUN=sd, na.rm=T)

WT[,5] <- rowMeans(WT[,1:3], na.rm = T)
RV524[,5] <- rowMeans(RV524[,1:3], na.rm=T)



pdf(pdfname <- "Figure5a.pdf", width = 6, height = 6)

plot(y=WT$V5, x=c(0,4,8,12,16,20), axes=F, type="l", ylab="log10 [IU / mL]", xlab="days post infection", col="orange", lwd=1.5, ylim=c(1,8))
axis(2, at=NULL)
axis(1, at=c(0,4,8,12,16,20))
grid(NULL, NULL)
#lines(y=RV522$V5, x=c(0,4,8,12,16,20), type="b", col="mediumseagreen", lwd=2.5)
#lines(y=RV523$V5, x=c(0,4,8,12,16,20), type="b", col="steelblue2", lwd=2.5)
lines(y=RV524$V5, x=c(0,4,8,12,16,20), type="l", col="red", lwd=1.5)
abline(h=2, col="darkgrey", lty=3, lwd=2)
text(x=5, y=2.05, "detection limit", cex=1.5)

points(x=rep(0,3),y=WT[1,1:3], col="orange", cex=1.25)
points(x=rep(4,3),y=WT[2,1:3], col="orange", cex=1.25)
points(x=rep(8,3),y=WT[3,1:3], col="orange", cex=1.25)
points(x=rep(12,3),y=WT[4,1:3], col="orange", cex=1.25)
points(x=rep(16,3),y=WT[5,1:3], col="orange", cex=1.25)
points(x=rep(20,3),y=WT[6,1:3], col="orange", cex=1.25)

points(x=rep(0,3),y=RV524[1,1:3], col="red", cex=1.25)
points(x=rep(4,3),y=RV524[2,1:3], col="red", cex=1.25)
points(x=rep(8,3),y=RV524[3,1:3], col="red", cex=1.25)
points(x=rep(12,3),y=RV524[4,1:3], col="red", cex=1.25)
points(x=rep(16,3),y=RV524[5,1:3], col="red", cex=1.25)
points(x=rep(20,3),y=RV524[6,1:3], col="red", cex=1.25)



#arrows(x0=c(0,4,8,12,16,20), y0=RV522$V5,  x1=c(0,4,8,12,16,20), y1=RV522$V5-RV522$V4, angle = 90,length = 0.05, col="mediumseagreen", lwd=2)
#arrows(x0=c(0,4,8,12,16,20), y0=RV522$V5,  x1=c(0,4,8,12,16,20), y1=RV522$V5+RV522$V4, angle = 90,length = 0.05, col="mediumseagreen", lwd=2)

#arrows(x0=c(0,4,8,12,16,20), y0=RV523$V5,  x1=c(0,4,8,12,16,20), y1=RV523$V5-RV523$V4, angle = 90,length = 0.05, col="steelblue2", lwd=2)
#arrows(x0=c(0,4,8,12,16,20), y0=RV523$V5,  x1=c(0,4,8,12,16,20), y1=RV523$V5+RV523$V4, angle = 90,length = 0.05, col="steelblue2", lwd=2)

arrows(x0=c(0,4,8,12,16,20), y0=WT$V5,  x1=c(0,4,8,12,16,20), y1=WT$V5-WT$V4, angle = 90,length = 0.05, col="orange", lwd=0.75)
arrows(x0=c(0,4,8,12,16,20), y0=WT$V5,  x1=c(0,4,8,12,16,20), y1=WT$V5+WT$V4, angle = 90,length = 0.05, col="orange", lwd=0.75)

arrows(x0=c(0,4,8,12,16,20), y0=RV524$V5,  x1=c(0,4,8,12,16,20), y1=RV524$V5-RV524$V4, angle = 90,length = 0.05, col="red", lwd=0.75)
arrows(x0=c(0,4,8,12,16,20), y0=RV524$V5,  x1=c(0,4,8,12,16,20), y1=RV524$V5+RV524$V4, angle = 90,length = 0.05, col="red", lwd=0.75)

legend("topleft", c("WT", "SILK+RVxFmut"), col=c("orange", "red"), bty="n", cex=1.5, lwd=3)



t.test(x = as.numeric(WT[3,1:3]), y = as.numeric(RV524[3,1:3]), exact=F, paired = F)$p.value
t.test(x = as.numeric(WT[4,1:3]), y = as.numeric(RV524[4,1:3]), exact=F, pared = F)$p.value
t.test(x = as.numeric(WT[5,1:3]), y = as.numeric(RV524[5,1:3]), exact=F, paired = F)$p.value

text("P = 0.031", x=8, y=4.55)
text("P = 0.0013", x=12, y=7.55)
text("P = 0.00017", x=16, y=7.85)

graphics.off()

system(paste("open", pdfname))

names(WT)[4:5] <- c("sd","mean")

names(RV524)[4:5] <- c("sd","mean")


