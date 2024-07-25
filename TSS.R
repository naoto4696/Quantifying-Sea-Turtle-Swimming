setwd("~/Desktop/Result")

library(vars)

# parameters

datanum <- 10
pcnum <- 2
arorder <- 10

minint <- 100
dint <- c(-20,-10,0,10,20)

# data input

setwd("./Result_pca")

for(i in 1:datanum){
  assign(paste0("pc_data",i), as.matrix(read.csv(paste0("pc_timeseries_all_",i,".csv"))[,1:pcnum]))
}

if (!dir.exists(paste0("../Result_ta"))) dir.create(paste0("../Result_ta"))
setwd("../Result_ta")

# calculate Separation point

for(i in 1:datanum){
  
  pc_part <- get(paste0("pc_data",i))
  datalength <- nrow(pc_part)
  repeatnum <- trunc(datalength/minint) - 1
  Psep <- 1
  Start <- 1
  Nsep <- 1
  
  for(j in 1:repeatnum){

    Csep <- j*minint
    
    pcsep1 <- pc_part[Start:(Csep-1),]
    pcsep2 <- pc_part[Csep:(Csep+minint-1),]
    pcint <- pc_part[Start:(Csep+minint-1),]
    
    lagsep1 <- VARselect(pcsep1, lag.max = arorder, type = "const")
    lagsep2 <- VARselect(pcsep2, lag.max = arorder, type = "const")
    lagint <- VARselect(pcint, lag.max = arorder, type = "const")
    
    Resultsep1 <- VAR(pcsep1, p = lagsep1$selection["AIC(n)"], type = "const")
    Resultsep2 <- VAR(pcsep2, p = lagsep2$selection["AIC(n)"], type = "const")
    Resultint <- VAR(pcint, p = lagint$selection["AIC(n)"], type = "const")
    
    AICsep <- AIC(Resultsep1) + AIC(Resultsep2)
    AICint <- AIC(Resultint)
    
    if(AICsep<AICint) {
      Psep <- c(Psep, Csep)
      Start <- Csep
    }
  }
  Psep <- c(Psep, datalength)
  
  summarydsep <- replicate(length(Psep)-1,list())
  AICmin <- 10000
  
  for(k in 2:(length(Psep)-1)){
    
    Start <- Psep[k-1]
    Csep <- Psep[k]
    End <- Psep[k+1]
    
    for(l in 1:length(dint)){
      
      Cdsep <- Csep + dint[l]
      
      pcdsep1 <- pc_part[Start:(Cdsep-1),]
      pcdsep2 <- pc_part[Cdsep:(End-1),]
    
      lagdsep1 <- VARselect(pcdsep1, lag.max = arorder, type = "const")
      lagdsep2 <- VARselect(pcdsep2, lag.max = arorder, type = "const")
      
      Resultdsep1 <- VAR(pcdsep1, p = lagdsep1$selection["AIC(n)"], type = "const")
      Resultdsep2 <- VAR(pcdsep2, p = lagdsep2$selection["AIC(n)"], type = "const")
      
      sumdsep1 <- summary(Resultdsep1)
      sumdsep2 <- summary(Resultdsep2)
      
      AICdsep <-  AIC(Resultdsep1) + AIC(Resultdsep2)
      
      if(l==1 || AICdsep<AICmin){
        AICmin <- AICdsep
        summarymin <- sumdsep1
        Psep[k] <- Cdsep
      }
    }
    summarydsep[[k]] <- summarymin
  }
  summarydsep[[k+1]] <- sumdsep2
  write.csv(Psep, paste0("Psep_data_",i,".csv"), row.names = FALSE)
  #write.csv(summarydsep, paste0("summary_data",i,".csv"), row.names = FALSE)
}


