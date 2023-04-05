
#####
#convert site HC5 values to single file


library(here)
library(tidyverse)

chems <- read_csv(here("input","chems","names.csv"))

files <- list.files(here("output","predictions"), pattern = "hSSD_GETREAL")

sum.tox <- data.frame(1,2,3,4,5)
colnames(sum.tox) <- c(colnames(read.csv(here("output","predictions",files[1])))[c(1,2,3,5)],"chem")


for(i in 1:length(chems[[1]])){
  pred.tox <- read.csv(here("output","predictions",files[i]))
  pred.tox <- pred.tox[,c(1,2,3,5)]
  pred.tox$chem <- rep(chems[[2]][i],length(pred.tox[,1]))
  sum.tox <- rbind(sum.tox,pred.tox)
}

sum.tox <- sum.tox[2:length(sum.tox$sample),]  

sum.tox$sample <- gsub("_tax","",sum.tox$sample)
unique(sum.tox$chem)
write.csv(sum.tox,here("output","hSSD_GETREAL_all_chems.csv"), row.names = FALSE)

#####
sum.tox <- read.csv(here("output","hSSD_GETREAL_all_chems.csv"))

#add the sample properties in
sample_properties <- read.csv(here("output", "sample_properties.csv"))

#add the additional property columns
sum.tox$site_id <- NA
sum.tox$date_id <- NA
sum.tox$latitude <- NA
sum.tox$longitude <- NA
sum.tox$brt20 <- NA
sum.tox$brt20_aggregate <- NA
sum.tox$brt12 <- NA
sum.tox$RT_AF <- NA
sum.tox$richness <- NA
sum.tox$ta7 <- NA

for(i in 1:length(unique(sum.tox$sample))){
  print(paste(i,unique(sum.tox$sample)[i]))
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$site_id <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$site_id
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$date_id <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$date_id
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$latitude <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$latitude
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$longitude <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$longitude
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$brt20 <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$brt20
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$brt20_aggregate <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$brt20_aggregate
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$brt12 <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$brt12
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$RT_AF <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$RT_AF
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$richness <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$richness
  sum.tox[sum.tox$sample == unique(sum.tox$sample)[i],]$ta7 <- sample_properties[sample_properties$sample_id == unique(sum.tox$sample)[i],]$ta7
}
unique(sum.tox$chem)
write.csv(sum.tox,here("output","hSSD_GETREAL_all_chems_properties.csv"), row.names = FALSE)
