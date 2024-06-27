setwd("~/Desktop/Result/Result_ta")

library(grDevices)
library(colorhcplot)
library(vioplot)

# parameter

datanum <- c(3,4,5,7,8,9,10)
clustnum <- 4

# color palette

plotcolor_sum <- c()
plotcolor_rec <- c()
for(plotclust in 1:clustnum){
  plotcolor_sum <- c(adjustcolor(palette()[1+plotclust], alpha=0.4), plotcolor_sum)
  plotcolor_rec <- c(adjustcolor(palette()[1+plotclust], alpha=1), plotcolor_rec)
#  assign(paste0("data",plotclust), Continue[Continue[,2]==plotclust,1])
}

nameseq <- c()
plotcolor_change <- c()
for(m in 1:clustnum){
  for(n in 1:clustnum){
    if(m!=n){
      mixedrgb <- colorRamp(c(palette()[1+m], palette()[1+n]))
      mixedcolor <- rgb(mixedrgb(0.3),maxColorValue=255)
      mixedcolor <- adjustcolor(mixedcolor, alpha=0.5)
      plotcolor_change <- c(plotcolor_change, mixedcolor)  
      nameseq <- c(nameseq, paste0(m," → ",n," "))
    }
  }
}

# plot

Change_mat <- c()
Ratio_mat <- c()
laymat1 <- c(1,1,1,1,1,2,3,3,3,3,3,4,0,0,0,0,0,0)
laymat2 <- c(5,5,5,5,5,6,7,7,7,7,7,8,0,0,0,0,0,0)
laymat3 <- c(9,9,9,9,9,10,11,11,11,11,11,12,13,13,13,13,13,14)
laymat <- rbind(laymat1, laymat2,laymat3)

png("vioplot_all.png", width = 2600, height = 1800, res=120)
layout(laymat)

for(dataid in datanum){
    
  setwd("../Result_ca")
  info <- read.csv(paste0("cluster_",dataid,".csv"))
  Psep <- info$Psep
  cluster <- info$cluster
  dlen <- length(cluster) -1
  
  Continue <- c()
  Psep_dif <- Psep[2:(dlen+1)] - Psep[1:dlen]
  Change <- rep(0, 100*(clustnum+1))
  for(m in 1:clustnum){
    for(n in 1:clustnum){
      if(m!=n){
        inputnum <- as.numeric(paste0(m,n))
        Change[inputnum] <- 0.0000000001   
      }
    }
  }
  for(j in 1:dlen){
    if(cluster[j]!=cluster[j+1] | j==dlen){
      Continue <- rbind(Continue, c(Psep_dif[j],cluster[j]))
      if(j!=dlen){
        inputnum <- as.numeric(paste0(cluster[j],cluster[j+1]))
        Change[inputnum] <- Change[inputnum] + 1   
      }
    } else{
      Psep_dif[j+1] <- Psep_dif[j] + Psep_dif[j+1]
    }
  }
    
  for(k in 1:clustnum){
    assign(paste0("data",k), Continue[Continue[,2]==k,1])
  }
  Ratio <- c(sum(data1),sum(data2),sum(data3),sum(data4))/sum(c(data1,data2,data3,data4))
  Ratio_vline <- c(Ratio[1],sum(Ratio[1:2]),sum(Ratio[1:3]))
  
  Change <- Change[Change!=0]
  for(k in 1:clustnum){
    start <- (k-1)*(clustnum-1) + 1
    end <- k*(clustnum-1)
    Change[start:end] <- Ratio[k]*(Change[start:end]/sum(Change[start:end]))
  }
  Change[Change<0.00001] <- 0
  
  par(mar=c(6, 4, 3, 2 ))
  par(xpd=T)
  vioplot(data4,data3,data2,data1 ,col=plotcolor_sum,horizontal=TRUE, ylim=c(0,2600),
          main="", xlab="", ylab="",xaxt="n", yaxt="n", bty="n",
          cex=1.5, lwd=1.5, cex.main=1.5, rectCol=plotcolor_rec, pchMed=20, colMed="white")
  axis(1, at=seq(0,2600,500), labels=seq(0,2600,500), cex.axis=2, line=0)
  if(dataid == 3 | dataid == 5 | dataid == 8) axis(2, at=seq(1,clustnum),formatC(c(paste0("ptn ",seq(clustnum,1,-1)))), cex.axis=2.4)
  mtext("duration", side = 1, line = 4.5, cex=1.6)
  legend(1950, par()$usr[4], legend=paste0("ID",dataid), cex=2.8  , box.lwd = 0, bg="transparent") 
  
  par(mar=c(6, 0, 3, 3.5))
  par(xpd=NA)
  barplot(as.matrix(Ratio), horiz=FALSE, main="", col=rev(plotcolor_sum), ylim=c(1,0), yaxt="n")
  axis(4, line=-0.12 ,at = c(0,0.2,0.4,0.6,0.8,1.0), cex.axis=2, las=1, lwd=0.5)
}
dev.off()






