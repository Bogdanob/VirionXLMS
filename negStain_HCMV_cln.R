# Load required libraries
library(dplyr)
library(data.table)
library(stringr)
library(lattice)
library(HistogramTools)
library(plotly)

# Set working directory
setwd("N:/R scripts/negative_stain_HCMV/")

file_names <- dir()[-length(dir())][grepl(".csv", dir())]


#read and combine anonymized tables (FIJI output)
df <- do.call(rbind, lapply(file_names, fread))

df$New_name <- str_extract(df$Label, "a[0-9]{3}")


#read key file for de-anonymization
key <- read.csv("./key/keyfile.csv")

df2 <- left_join(df, key, by="New_name")

df3 <- na.omit(df2)

#combine all results in a df according to their experimental group
dat1 <- aggregate(df3$Major, by=list(df3$Group), c)



#plot histograms of measured diameters per experimental group
histogram((dat1$x)[[2]], breaks=100, main=dat1$Group.1[[2]],col="green", xlim = c(0,500), ylim=c(0,15))
histogram((dat1$x)[[1]], breaks=30, main=dat1$Group.1[[1]],col="blue", xlim=c(0,500), ylim=c(0,15))


histogram((dat1$x)[[4]], breaks=100, main=dat1$Group.1[[4]],col="green", xlim = c(0,500), ylim=c(0,15))
histogram((dat1$x)[[3]], breaks=50, main=dat1$Group.1[[3]],col="blue", xlim=c(0,500), ylim=c(0,15))


plot(density((dat1$x)[[2]]), col="green", ylim=c(0,0.03), xlim=c(0,350))
lines(density((dat1$x)[[1]]), add=T)


plot(density((dat1$x)[[4]]), col="green", ylim=c(0,0.03), xlim=c(0,350))
lines(density((dat1$x)[[3]]), add=T)


plot(density((dat1$x)[[1]]), col="green", dat1$Group.1[[4]])
lines(density((dat1$x)[[3]]), add=T)


plot(density((dat1$x)[[2]]), col="green", dat1$Group.1[[4]])
lines(density((dat1$x)[[4]]), add=T)


#table(dat1$x[[3]]< median((dat1$x)[[1]])+2*sd((dat1$x)[[1]]))



#we define the diameter limits of a particle with 160-220 nm since it reflects the mean minus/plus two times sd for the purified, non cross-linked sample

lowlim <- 160

uplim <- 220


#calculate how many particles per experiment meet cut-off
under <- c(sum(dat1$x[[1]]< lowlim, nar.rm=T),
           sum(dat1$x[[2]]< lowlim, nar.rm=T),
           sum(dat1$x[[3]]< lowlim, nar.rm=T),
           sum(dat1$x[[4]]< lowlim, nar.rm=T))

over <- c(sum(dat1$x[[1]]> uplim, nar.rm=T),
          sum(dat1$x[[2]]> uplim, nar.rm=T),
          sum(dat1$x[[3]]> uplim, nar.rm=T),
          sum(dat1$x[[4]]> uplim, nar.rm=T))

within <- c(sapply(dat1$x, length))-under-over

wit <- within/(within+under+over)*100

ov <- over/(within+under+over)*100

un <- under/(within+under+over)*100


#plot bar-diagramm for fractions below, in-between and above diameter limits per experiment using plotly package
plotdf <- data.frame("condition" = str_remove_all(str_remove_all(dat1$Group.1, "[0-9]*_LM_negSt_HCMV_"), "_1min_2perc"), "under" = un, "within" = wit, "over" = ov)

toplabels = c(paste("below ",lowlim," nm", sep=""), paste(lowlim," - ",uplim," nm", sep = ""), paste("over ",uplim," nm", sep=""))


fig <- plot_ly(plotdf, x=~under, y=~condition, type="bar", orientation="h", 
               marker=list(color = 'rgba(38, 24, 74, 0.8)',
                           line = list(color = 'rgb(248, 248, 249)', width = 1)))
fig <- fig %>% add_trace(x=~within, marker = list(color = 'rgba(71, 58, 131, 0.8)'))
fig <- fig %>% add_trace(x=~over, marker = list(color = 'rgba(122, 120, 168, 0.8)'))

fig <- fig %>% layout(xaxis = list(title = "",
                                   
                                   showgrid = FALSE,
                                   
                                   showline = FALSE,
                                   
                                   showticklabels = FALSE,
                                   
                                   zeroline = FALSE,
                                   
                                   domain = c(0.15, 1)),
                      
                      yaxis = list(title = "",
                                   
                                   showgrid = FALSE,
                                   
                                   showline = FALSE,
                                   
                                   showticklabels = FALSE,
                                   
                                   zeroline = FALSE),
                      
                      barmode = 'stack',
                      
                      paper_bgcolor = 'rgb(248, 248, 255)', plot_bgcolor = 'rgb(248, 248, 255)',
                      
                      margin = list(l = 120, r = 10, t = 140, b = 80),
                      
                      showlegend = FALSE)


fig <- fig %>% add_annotations(xref = 'paper', yref = 'condition', x = 0.14, y = plotdf$condition,
                               
                               xanchor = 'right',
                               
                               text = plotdf$condition,
                               
                               font = list(family = 'Arial', size = 12,
                                           
                                           color = 'rgb(67, 67, 67)'),
                               
                               showarrow = FALSE, align = 'right')


fig <- fig %>% add_annotations(xref = 'x', yref = 'paper',
                               
                               x = c(plotdf$under[4]/ 2, plotdf$under[4] + plotdf$within[4] / 2, plotdf$under[4] + plotdf$within[4] + plotdf$over[4] / 2),
                               
                               y = 1.05,
                               
                               text = toplabels,
                               
                               font = list(family = 'Arial', size = 12,
                                           
                                           color = 'rgb(67, 67, 67)'),
                               
                               showarrow = FALSE)


fig <- fig %>% add_annotations(xref = 'x', yref = 'y',
                               
                               x = un/2, y = plotdf$condition,
                               
                               text = paste(round(plotdf[,"under"], digits = 2), '%'),
                               
                               font = list(family = 'Arial', size = 12,
                                           
                                           color = 'rgb(248, 248, 255)'),
                               
                               showarrow = FALSE)


fig <- fig %>% add_annotations(xref = 'x', yref = 'y',
                               
                               x = un+wit/2, y = plotdf$condition,
                               
                               text = paste(round(plotdf[,"within"], digits = 2), '%'),
                               
                               font = list(family = 'Arial', size = 12,
                                           
                                           color = 'rgb(248, 248, 255)'),
                               
                               showarrow = FALSE)


fig <- fig %>% add_annotations(xref = 'x', yref = 'y',
                               
                               x = un+wit+ov/2, y = plotdf$condition,
                               
                               text = paste(round(plotdf[,"over"], digits = 2), '%'),
                               
                               font = list(family = 'Arial', size = 12,
                                           
                                           color = 'rgb(248, 248, 255)'),
                               
                               showarrow = FALSE)


fig


