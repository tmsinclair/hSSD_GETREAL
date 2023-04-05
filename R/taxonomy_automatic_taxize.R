library(taxize)
library(here)

#file location of the species/ taxonomic rank names to be filled in
file <- here("data","sites","site_sprn_0-5_family.csv")

#the relevant column titles for the hSSD model
tax_ranks <- c("Kingdom",	"Phylum",	"Subphylum",	"Superclass",	"Class",	"Order",	"Family",	"Genus",	"Latin")

#read in the species/ taxonomic rank names to be filled in
tax_raw <- read.csv(file,stringsAsFactors = FALSE)

#correct the column names for the hSSD model
colnames(tax_raw) <- tax_ranks

#remove unknown for when the species is unknown
tax_raw$Genus <- gsub("\\s.*", "", gsub("Unknown \\(","",tax_raw$Latin))

#only get the unique values for genus
tax_incomplete <- as.character(unique(tax_raw$Family)) #as.character(unique(tax_raw[is.na(tax_raw$Family),]$Genus))

tax_raw <- tax_raw[1:length(tax_incomplete),]
tax_raw$Family <- tax_incomplete

tax_raw <- tax_raw[!is.na(tax_raw$Family),]

#for ever taxa
for(i in tax_incomplete){
#get a classification from the ncbi database
  tax_database <- classification(i ,db = "ncbi")
  #if there is no entry, skip
  if(is.na(tax_database)) next
  #create a dataframe to store the data
  tax_database <- data.frame(tax_database[1])
  tax_database <- data.frame(name = tax_database[[1]],
                             rank = tax_database[[2]],
                             id = tax_database[[3]],
                             stringsAsFactors = FALSE)
						
	 #for every taxonomic rank of interest, save it into the database
  for(j in 1:length(tax_raw)){
    if(!any(tax_database$rank == tolower(colnames(tax_raw)[j]))) next
    tax_raw[tax_raw$Family == i,j] <- tax_database[tax_database$rank == tolower(colnames(tax_raw)[j]),]$name
  }
}

#process the data so it works with the hSSD model
for(k in 1:length(tax_raw)){
 tax_raw[,k] <- tolower(tax_raw[,k])
}

#correct the column titles
colnames(tax_raw) <- c("Kingdom",	"Phylum_division",	"Subphylum",	"Superclass",	"Class",	"Order",	"Family",	"Genus",	"Latin")

#correct kingdo and superclass m to align with hSSD model
tax_raw$Kingdom <- "animalia"
tax_raw$Superclass <-NA

#save the file
write.csv(tax_raw, here("output","sites", "site_sprn_0-5.csv"), row.names = FALSE)

