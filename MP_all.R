setwd("~/Desktop/Result")

library(progress)
library(grDevices)

# set parameters

dataid <- 10
datanum <- 10

plot_f <- 10
pc1min <- -1.5
pc1max <- 1.5
pc2min <- -1.5
pc2max <- 1.5
plot_scale <- 0.1
plotcolor <- "black"
color_palette <- colorRampPalette(c(plotcolor, "white"))(n = plot_f)

plotnum <- 200
T <- TRUE

# data read 

setwd("./Result_ce")
CT <- data.frame(read.csv(paste0("CT_",dataid,".csv")))

setwd("../Result_fa")
interpolated <- data.frame(read.csv(paste0("interpolated_",dataid,".csv")))
for(i in 1:datanum){
  if(i==1) coef <-  data.frame(read.csv(paste0("coef_",i,".csv")))
  else coef <- rbind(coef,data.frame(read.csv(paste0("coef_",i,".csv"))))
}

setwd("../Result_pca")
pc <- read.csv(paste0("Pc_timeseries_all_",dataid,".csv"))

setwd("../Result_ta")
Psep <- as.ts(read.csv(paste0("Psep_",dataid,".csv")))

# settings

pc1 <- pc[,1]
pc2 <- pc[,2]
pca_result <- prcomp(coef)

plotlen <- length(pc1)

xmin <- min(CT$x_head) - 10
xmax <- max(CT$x_head) + 10
ymin <- min(CT$y_head) - 10
ymax <- max(CT$y_head) + 10

# contour settings

pc1_int <- (pc1max-pc1min)/6
pc2_int <- (pc2max-pc2min)/6
len <- length(pca_result$center)/4
seqlength <- 300
plotseq <- seq(1,seqlength)
xlongseq <- vector()
ylongseq <- vector()
for(i in -2:2){
  for(j in -2:2){
    xseq <- rep(i*pc1_int,seqlength)
    yseq <- rep(j*pc2_int,seqlength)
    pc_scores <- c(i*pc1_int, j*pc2_int, rep(0, ncol(pca_result$x) - 2)) 
    re_coef <- pca_result$center + pc_scores %*% t(pca_result$rotation)
    for(k in 1:len){
      xseq <- xseq + (re_coef[k+0*len]*cos((k-1)*plotseq) +
                        re_coef[k+1*len]*sin((k-1)*plotseq))*plot_scale
      yseq <- yseq + (re_coef[k+2*len]*cos((k-1)*plotseq) + 
                        re_coef[k+3*len]*sin((k-1)*plotseq))*plot_scale
    }
    xlongseq <- c(xlongseq, xseq)
    ylongseq <- c(ylongseq, yseq)
  }
}

if (!dir.exists(paste0("../Result_mp_",dataid))) dir.create(paste0("../Result_mp_",dataid))
setwd(paste0("../Result_mp_",dataid))

pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = plotlen  # 総数を設定
)

# plot movement graph

print("Ploting movement graph...")

for(i in 1:plotlen){
  
  pb$tick()
  
  plot_part <- interpolated[interpolated$interpolated_coords==(i-1),]
  
  png(filename=paste0(sprintf("%05d", i),".png"), width = 800, height = 800)
  layout(matrix(c(1,1,1,2,2,2,3,4,4,3,4,4,5,4,4,5,4,4), nrow = 6, byrow = TRUE)) 
  
  plotinit <- i-plotnum
  plotend <- i+plotnum

  par(lwd=1.5)
  par(mar=c(0, 6, 5, 6))
  plot(pc1, xlim=c(plotinit,plotend), ylim=c(pc1min, pc1max),
       xlab="", ylab="", xaxt="n", main="Timeseries", cex.main=2,
       col="black", type="l", pch=20, lwd=1)
  abline(v=Psep, col="red", lwd=2, lty=2)
  points(i,pc1[i],col=plotcolor, cex=2, pch=19)
  mtext("pc1", side = 2, line = 2.5)

  par(mar=c(4, 6, 1, 6))
  plot(pc2, xlim=c(plotinit,plotend), ylim=c(pc1min, pc1max),
       xlab="", ylab="", main="", cex.main=2,
       col="black", type="l", pch=20, lwd=1)
  abline(v=Psep, col="red", lwd=2, lty=2)
  points(i,pc2[i],col=plotcolor, cex=2, pch=19)
  mtext("time", side = 1, line = 2.5)
  mtext("pc2", side = 2, line = 2.5)
  
  par(mar=c(6, 6, 4, 3))
  plot(CT$x_head[(i+1):(i+plot_f)],CT$y_head[(i+1):(i+plot_f)],
       xlim=c(xmin,xmax), ylim=c(ymin, ymax),
       xlab="", ylab="", main="Position", cex.main=2,
       col=color_palette, pch=20)
 par(new=TRUE)
  plot(CT$x_head[i],CT$y_head[i],
       xlim=c(xmin,xmax), ylim=c(ymin, ymax),
       xlab="", ylab="", xaxt="n", yaxt="n",
       col=plotcolor, pch=19, cex=2, bty="n")
  mtext("x(pixel)", side = 1, line = 2.5)
  mtext("y(pixel)", side = 2, line = 2.5)
  
  par(mar=c(6, 3, 4, 6))
  plot(xlongseq, ylongseq, xlim=c(pc1min,pc1max), ylim=c(pc2min,pc2max), 
       main="PC plot", cex.main=2, xlab="", ylab="", col=8, pch=20, cex=0.2)
  mtext("pc1", side = 1, line = 2.5)
  mtext("pc2", side = 2, line = 2.5)
 par(new=TRUE)
  plot(pc1[(i+1):(i+plot_f)], pc2[(i+1):(i+plot_f)],
       xlim=c(pc1min,pc1max), ylim=c(pc2min,pc2max),xlab="", ylab="",
       xaxt="n", yaxt="n", col=color_palette, pch=20, cex=2)
 par(new=TRUE)
  plot(pc1[i], pc2[i],
       xlim=c(pc1min,pc1max), ylim=c(pc2min,pc2max), xlab="", ylab="",
       xaxt="n", yaxt="n", col=plotcolor, pch=20, cex=5, bty="n")
  
  par(mar=c(6, 6, 4, 3))
  plot(plot_part$interpolated_x,plot_part$interpolated_y, 
       xlim=c(-2.5,2.5), ylim=c(-2.5,2.5), 
       xlab="", ylab="", xaxt="n", yaxt="n", type="l", lwd=2,
       main="Shape", cex.main=2, col=plotcolor, pch=20)
  
  dev.off()
}