setwd("~/Desktop/Result")

# parameters

datanum <- 10

pc1min <- -1.5
pc1max <- 1.5
pc2min <- -1.5
pc2max <- 1.5
plot_scale <- 0.1 

setwd("./Result_fa")

datalen <- c()
for(i in 1:datanum){
  if(i==1){
    coef <-  data.frame(read.csv(paste0("coef_",i,".csv")))
  }
  else{
    coef <- rbind(coef,data.frame(read.csv(paste0("coef_",i,".csv"))))
  }
  datalen[i] <- nrow(data.frame(read.csv(paste0("coef_",i,".csv"))))
}

# settings

if (!dir.exists("../Result_pca")) dir.create("../Result_pca")
setwd("../Result_pca")

coef_len <- length(coef$xcos_0)
pca_result <- prcomp(coef)
pc1_timeseries <- as.ts(pca_result$x[1:coef_len,1])
pc2_timeseries <- as.ts(pca_result$x[1:coef_len,2])
pc3_timeseries <- as.ts(pca_result$x[1:coef_len,3])
pc4_timeseries <- as.ts(pca_result$x[1:coef_len,4])

pc_timeseries <- cbind(pc1_timeseries, pc2_timeseries, pc3_timeseries, pc4_timeseries)

for(i in 1:datanum){
  if(i==1) init <- 1
  else init <- sum(datalen[1:(i-1)])
  end <- sum(datalen[1:i])
  write.csv(pc_timeseries[init:end,], paste0("pc_timeseries_all_",i,".csv"), row.names = FALSE)
}

pca_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pca_rotation <- pca_result$rotation
pca_explained_rotation <- rbind(pca_explained, pca_rotation)
write.csv(pca_explained_rotation, "pca_result_all.csv", row.names = TRUE)

pc1_int <- (pc1max-pc1min)/6
pc2_int <- (pc2max-pc2min)/6
len <- length(pca_result$center)/4
seqlength <- 1000
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

#plot pc contribution

png(filename="PC_contribution_all.png", width = 2000, height = 600, res=150)
par(mar=c(6, 5, 4, 4))
plot(pca_result$rotation[-c(1,11,21,31),1],  xlab = "", ylab = "", xaxt="n", type="h",
     lwd=3, col=rgb(1, 0, 0, alpha = 0.8), ylim=c(-0.7,0.7), cex.axis=1) 
axis(1, at=seq(1.1,36.1,1), las=2, cex.axis=1,
     labels=c("w_a1","w_a2","w_a3","w_a4","w_a5","w_a6","w_a7","w_a8","w_a9",
              "w_b1","w_b2","w_b3","w_b4","w_b5","w_b6","w_b7","w_b8","w_b9",
              "w_c1","w_c2","w_c3","w_c4","w_c5","w_c6","w_c7","w_c8","w_c9",
              "w_d1","w_d2","w_d3","w_d4","w_d5","w_d6","w_d7","w_d8","w_d9"))
abline(h = 0, col = "gray", lty = 2)
abline(v = seq(1.1,36.1,2), col = "gray", lty = 2)
legend("topright", legend = c("pc1 ", "pc2 "), col = c("red", "blue"), lty=1, lwd=2, cex=1)
par(new=T)
par(mar=c(6, 5.3, 4, 3.7))
plot(pca_result$rotation[-c(1,11,21,31),2],  xlab = "", ylab = "", bty="n", xaxt="n", yaxt="n", 
     type="h", lwd=3, col=rgb(0, 0, 1, alpha = 0.8), ylim=c(-0.7,0.7))
mtext("Coefficients", side = 1, line = 4, cex=1.2)
mtext("Contribution", side = 2, line = 3, cex=1.2)
dev.off()

# plot variance graph

png(file="PC_explain_all.png", width = 1800, height = 1200, res=150)
par(mfrow=c(1,1))
plot(pca_explained, type="h", lwd=2.5, xlab = "", ylab = "",col=adjustcolor("black", alpha=0.6))
mtext("pc axis", side = 1, line = 3, cex=1.2)
mtext("explanatory power", side = 2, line = 3, cex=1.2)
dev.off()

# plot PC plane

png(file="PC_plane.png", width = 1600, height = 1600, res=200)
plot(xlongseq, ylongseq, xlim=c(pc1min+0.2,pc1max-0.2), ylim=c(pc2min+0.2,pc2max-0.2), 
     main="", cex.main=2, xlab="", ylab="", col=1, pch=20, cex=0.2, cex.axis=1.2, lwd=0.5)
mtext("PC1", side = 1, line = 2.9, cex=1.4)
mtext("PC2", side = 2, line = 2.9, cex=1.4)
par(lwd=1.5)
box()
dev.off()
