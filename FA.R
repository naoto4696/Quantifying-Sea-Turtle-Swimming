setwd("~/Desktop/Result")

library(progress)

# parameters

cutoff <- 2.5
dataid <- 1

#settings

setwd("./Result_ce")

data <- read.csv(paste0("output_rechain_",dataid,".csv"))
data$norm_rotated_y <- -1*data$norm_rotated_y
coordlen <- max(data$coords)+1
newdata <- matrix()

if (!dir.exists("../Result_fa")) dir.create("./Result_fa")
setwd("./Result_fa")

# add cumulative difference row

pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = coordlen  
)
print("adding cummulative row...")
  
for (i in 0:(coordlen-1)){
  pb$tick()
  data_part <- data[data$coords==i,]
  newdata_part <- rbind(data_part,data_part[1,])
  data_part_dummy <- newdata_part[-1,]

  difference <- sqrt((data_part$norm_rotated_x - data_part_dummy$norm_rotated_x)**2 +
                     (data_part$norm_rotated_y - data_part_dummy$norm_rotated_y)**2)
  difference <- c(0,difference)
  difference <- cumsum(difference)
  
  newdata_part <- cbind(newdata_part,difference)
  newdata_part$difference <- newdata_part$difference*(2*pi/max(newdata_part$difference))
  if (i==0)  newdata <- newdata_part 
  else  newdata <- rbind(newdata,newdata_part) 
}

# interpolation

interpolated <- matrix()
interpolatelength <- 100
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = coordlen  # 総数を設定
)
print("interpolating data...")

for (l in 0:(coordlen-1)){
  
  pb$tick()
  
  interpolated_part <- newdata[newdata$coords==l,]
  
  interpolated_coords <- rep(l,interpolatelength)
  interpolated_diff <- approx(interpolated_part$difference, interpolated_part$norm_rotated_x, 
                              xout = seq(0,2*pi,length.out = interpolatelength))$x
  interpolated_x <- approx(interpolated_part$difference, interpolated_part$norm_rotated_x, 
                           xout = seq(0,2*pi,length.out = interpolatelength))$y
  interpolated_y <- approx(interpolated_part$difference, interpolated_part$norm_rotated_y, 
                           xout = seq(0,2*pi,length.out = interpolatelength))$y
  
  interpolated_x[interpolatelength]  <- interpolated_part$norm_rotated_x[1]
  interpolated_y[interpolatelength]  <- interpolated_part$norm_rotated_y[1]
  
  if(any(abs(interpolated_x) > cutoff) | any(abs(interpolated_y) > cutoff)){
    interpolated_part <- cbind(interpolated_coords, interpolated_dummy)
  }
  else {
    interpolated_part <- cbind(interpolated_coords, interpolated_x, interpolated_y, interpolated_diff)
    interpolated_dummy <- cbind(interpolated_x, interpolated_y, interpolated_diff)
  }
  
  if (l==0)  interpolated  <- interpolated_part 
  else  interpolated <- rbind(interpolated,interpolated_part) 
}

interpolated <- data.frame(interpolated)
write.csv(interpolated, paste0("interpolated_",dataid,".csv"), row.names = FALSE)

# plot interpolated data

maxcoef <- 10
coef <- data.frame()
pb <- progress_bar$new(
  format = "[:bar] :percent :current/:total (:eta)",
  total = coordlen  # 総数を設定
)
print("calculating fourier coefficients...")

for (m in 0:(coordlen-1)){
  
  pb$tick()
  
  plot_part <- interpolated[interpolated$interpolated_coords==m,]
  maxdif <- max(plot_part$interpolated_diff)
  datalength <- length(plot_part$interpolated_diff)
  
  fft_x <- fft(plot_part$interpolated_x)
  fft_y <- fft(plot_part$interpolated_y)
  coef_x_cos <- 2*Re(fft_x[1:maxcoef])/datalength
  coef_x_sin <- -2*Im(fft_x[1:maxcoef])/datalength
  coef_y_cos <- 2*Re(fft_y[1:maxcoef])/datalength
  coef_y_sin <- -2*Im(fft_y[1:maxcoef])/datalength
  coef_for_x <- c(coef_x_cos,coef_x_sin)
  coef_for_y <- c(coef_y_cos,coef_y_sin)
  coef_mix <- c(coef_for_x,coef_for_y)
  pch_coef <- c(0:(maxcoef-1))
  
  if(i == 0) coef <- coef_mix
  else coef <- rbind(coef, coef_mix)
  
}

colnames(coef) <- c(paste0("xcos_",c(0:(maxcoef-1))),paste0("xsin_",c(0:(maxcoef-1))),
                    paste0("ycos_",c(0:(maxcoef-1))),paste0("ysin_",c(0:(maxcoef-1))))
coef[, grepl("0", names(coef))] <- 0
write.csv(coef, paste0("coef_",dataid,".csv"), row.names = FALSE)




 
