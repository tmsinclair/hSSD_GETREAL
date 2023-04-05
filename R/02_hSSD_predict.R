#install.packages("tidyverse")
library(tidyverse)
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("rhdf5")
library(rhdf5)
library(Matrix)
library(here)

#readin in the chemicals
chems <- read.csv(here("input","chems","names.csv"))
#read in the tox tax data
i <- 1
print(chems[[1]][i])
print(i)


#Prime R with the hSSD function
source(here("R","hSSD_function_2-1.R"))

#read in the assemblage for toxicity to be predicted
newspecies <- read.csv(here("output", "unique_taxa_summer.csv"))

for(i in 1:length(chems[[2]])){
  print(chems[[1]][i])

  #read in the tox data
  newdata <- read.csv(here("input", "tox", paste(chems[[1]][i], "_tox.csv", sep = "")))
  newdata$latin <- tolower(newdata$latin)
  #read in the taxonomy of the tox data
  newspecies.testedonly <- read.csv(here("input", "tox", paste(chems[[1]][i], "_tox_tax.csv", sep = "")))

  result <- MCMC.newchemical(
    newdata,
    here("input","database","test","toms-model122b.h5"),
    numeric(0),
    newspecies=newspecies,
    newspecies.testedonly=newspecies.testedonly,
    N=10000, burn=2500, thin=10, shout=100,
    detailed.output=TRUE
  )
  
  #have a look at the predicted outputs from the runs
  result$mu0
  
  write_csv(data.frame(result$mu0), here("output","predictions",paste("predicted_runs_",chems[[1]][i],"_tox_tax_GETREAL.csv", sep = "")))
  
  #mean the runs (as they are logged, this is a geomean)
  pred.tox <- data.frame(colnames(result$mu0),apply(result$mu0, 2, mean))
  
  #save the files
  colnames(pred.tox) <- c("latin","tox")
  write.csv(pred.tox, here("output","predictions",paste("predicted_tox_",chems[[1]][i],"_GETREAL.csv", sep = "")), row.names = FALSE)
}

