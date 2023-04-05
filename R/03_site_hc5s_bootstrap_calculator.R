library(RCPSrl)
library(tidyverse)
library(MASS)
library(reshape2)
library(here)

#code shameless purloined from https://edild.github.io/ssd/

myboot <- function(fit, p){
  # resample from fitted distribution
  xr <- rlnorm(fit$n, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
  #   # for a nonparametric boostrap resample from data
  #   xr <- sample(df$val, length(df$val), replace = TRUE)
  # fit distribition to new data
  fitr <- fitdistr(xr, 'lognormal')
  # return HCp
  hc5r <- qlnorm(p, meanlog = fitr$estimate[1], sdlog = fitr$estimate[2])
  return(hc5r)
}

myboot2 <- function(fit, newxs){
  # resample
  xr <- rlnorm(fit$n, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
  # fit to resample
  fitr <- fitdistr(xr, 'lognormal')
  # predict for new data
  pyr <- plnorm(newxs, meanlog = fitr$estimate[1], sdlog = fitr$estimate[2])
  return(pyr)
}

#locate output files
files <- list.files(here("output","predictions"), pattern = "predicted_tox_")

files <-  files[grepl("GETREAL", files)]

samples <- list.files(path =  here("output","sample"), pattern="*.csv")
#readin in the chemicals
chems <- read_csv(here("input","chems","names.csv"))

for(i in 1:length(chems[[1]])){
  
  # #Check if there is an output folder for that chemical, if not, create it
  # ifelse(!dir.exists(here("output","predictions",chems[[2]][i])), dir.create(here("output","predictions",chems[[2]][i])), FALSE)
  # ifelse(!dir.exists(here("figures","predictions",chems[[2]][i])), dir.create(here("figures","predictions",chems[[2]][i])), FALSE)
  # 
  
  pred.tox <- read.csv(here("output","predictions",files[i]))
  
  #Only needed if the column names and . are still left in the predicted tox files
  
  #if(!any(grepl("\\.",pred.tox$latin))){
    # colnames(pred.tox) <- c("latin","tox")
    # # 
    # pred.tox$latin <- gsub("\\."," ", pred.tox$latin)
    # pred.tox$latin <- gsub(" sp "," sp\\.", pred.tox$latin)
    # # 
    # write.csv(pred.tox,here("output","predictions",paste("predicted_runs_",chems[[2]][i],"_all_sites.csv", sep = "")), row.names = FALSE)
  #}

  #tox.sites <- list.files(path =  "./hSSDs/CPS/tox_sites/", pattern="*.csv")
  
  summary <- c()
  
  for(j in 1:length(samples)){
  
    newspecies <- read.csv(here("output","sample",samples[j]), as.is=TRUE)
    if(nrow(newspecies) < 1){
      print(samples[j])
      next
    }
    sample.tox <- pred.tox[pred.tox$latin %in% tolower(newspecies$Latin),]
    sample.tox <- sample.tox[order(sample.tox$latin),]
    
    sample.tox$frac <- ppoints(sample.tox$tox, 0.5)
    
    # ggplot(data = site.tox) +
    #   geom_point(aes(x = tox, y = frac), size = 2) +
    #   # geom_text(aes(x = x, y = frac, label = X), hjust = 1.1, size = 4) +
    #   theme_bw() +
    #   scale_x_log10(limits = c(0.0075, max(site.tox$tox))) +
    #   labs(x = expression(paste('Concentration of Chlorpyrifos [ ', mu, 'g ', L^-1, ' ]')), 
    #        y = 'Fraction of species affected')
    sample.tox$tox <- 10^sample.tox$tox
    fit <- fitdistr(sample.tox$tox, 'lognormal')
    hc5 <- qlnorm(0.05, meanlog = fit$estimate[1], sdlog = fit$estimate[2])
    
    #set.seed(1234)
    hc5_boot <- replicate(1000, myboot(fit, p = 0.05))
    # as.numeric(log10(quantile(hc5_boot, probs = c(0.025, 0.975))))
    # 
    # hc5_boot <- replicate(1000, myboot(fit, p = 0.05))
    # hc50_boot <- replicate(1000, myboot(fit, p = 0.5))
    # hc95_boot <- replicate(1000, myboot(fit, p = 0.95))
    # 
    # hc5_sum <- data.frame(
    #   hc5 = as.numeric(log10(quantile(hc5_boot, probs = c(0.025, 0.5, 0.975)))),
    #   hc50 = as.numeric(log10(quantile(hc50_boot, probs = c(0.025, 0.5, 0.975)))),
    #   hc95 = as.numeric(log10(quantile(hc95_boot, probs = c(0.025, 0.5, 0.975)))))
    # hc5_sum    
    # # newdata to predict
    # newxs <- 10^(seq(log10(0.01), log10(max(site.tox$tox)), length.out = 1000))
    # boots <- replicate(1000, myboot2(fit, newxs))
    
    summary <- rbind(summary,c(gsub(".csv","",gsub("sample ","",samples[j])),hc5,  quantile(hc5_boot, probs = c(0.025,0.5, 0.975))))
    colnames(summary) <- c("sample","HC5","Lower","50th","Upper")
    # fancy plot
    
    # # extract boostrap values
    # bootdat <- data.frame(boots)
    # bootdat$newxs <- newxs
    # bootdat <- melt(bootdat, id = 'newxs')
    # # extract CI
    # cis <- apply(boots, 1, quantile, c(0.025, 0.975))
    # rownames(cis) <- c('lwr', 'upr')
    # # add fitted values
    # pdat <- data.frame(newxs, py = plnorm(newxs, meanlog = fit$estimate[1], sdlog = fit$estimate[2]))
    # # add CI
    # pdat <- cbind(pdat, t(cis))
    # # x koordinates for species names (better use lower Ci instead of fit...)
    # site.tox$fit <- 10^(log10(qlnorm(site.tox$frac, meanlog = fit$estimate[1], sdlog = fit$estimate[2])) - 0.4)
    # # plot
    # site.tox
    # image <- ggplot() +
    #   #(data = bootdat, aes(x = newxs, y = value, group = variable), col = 'steelblue', alpha = 0.05) +
    #   geom_point(data = site.tox, aes(x = tox, y = frac)) +
    #   geom_line(data = pdat, aes(x = newxs, y = py), col = 'red') +
    #   geom_line(data = pdat, aes(x = newxs, y = lwr), linetype = 'dashed') +
    #   geom_line(data = pdat, aes(x = newxs, y = upr), linetype = 'dashed') +
    #   geom_text(data = site.tox, aes(x = fit, y = frac, label = latin), hjust = 1, size = 2.5) +
    #   theme_bw() +
    #   scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000,10000,100000), limits = c(0.003, max(site.tox$tox))) +
    #   labs(x = paste('Concentration of ',chems[[1]][i]," (ug/L)",sep = ""),
    #        y = 'Fraction of species affected')
    # ggsave(file = here("figures","predicted_tox",chems[[2]][i],paste("hSSD_",gsub(" ","-",gsub(".csv","",sites[j])),"_simple.jpeg",sep = "")), plot= image)
    print(length(samples)-j)
  }
  summary <- data.frame(summary)
  colnames(summary) <- c("sample","HC5_aj","lower","HC5_50th","upper")
  write_csv(summary,here("output","predictions",paste("hSSD_GETREAL_",chems[[1]][i], ".csv", sep = "")))
}
