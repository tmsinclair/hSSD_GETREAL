library(here)

#read in the raw data
EU_sites <- read.csv(here("input", "220124_data_tom_javier.csv"))
sweden_sites <- readRDS(here("input", "sweden_macroinvertebrates.rds"))

#before doing anything with taxonomy, get a list of samples and their properties

sample_props <- NULL

for(i in 1:length(unique(EU_sites$sample_id))){
  print(i)
  #get the first row from each sample
  sample <- EU_sites[EU_sites$sample_id == unique(EU_sites$sample_id)[i],][1,]

  #extract the relevant parameters
  sample_props <- rbind(sample_props, data.frame(sample$sample_id, sprintf("%05d",sample$site_id), sprintf("%04d",sample$date_id), sample$latitude, sample$longitude, sample$brt20, sample$brt20_aggregate, sample$brt12, sample$richness))
}

#sort column names out
colnames(sample_props) <- c("sample_id", "site_id", "date_id", gsub("sample\\.", "", colnames(sample_props))[4:9])

#sweden props

sweden_props <- NULL

for(i in 1:length(unique(sweden_sites$gr_sample_id))){
  print(i)
  #get the first row from each sample
  sample <- sweden_sites[sweden_sites$gr_sample_id == unique(sweden_sites$gr_sample_id)[i],][1,]
  
  #extract the relevant parameters
  sweden_props <- rbind(sweden_props, data.frame(sample$gr_sample_id, sprintf("%05d",as.numeric(sample$site_id)), sprintf("%04d", as.numeric(sample$date_id)), sample$y.coord, sample$x.coord, NA, NA, sample$brt12, nrow(sweden_sites[sweden_sites$gr_sample_id == unique(sweden_sites$gr_sample_id)[i],])))
}

colnames(sweden_props) <- colnames(sample_props)

sample_props<- rbind(sample_props, sweden_props)

#rearrange the data by site
sample_props <- sample_props[order(sample_props$date_id),]
sample_props <- sample_props[order(sample_props$site_id),]

sample_props

write.csv(sample_props, here("output", "sample_properties.csv"), row.names = FALSE)

#remove the subclass level as it is empty
EU_sites <- EU_sites[,-6]

#make taxonomy the same as the hSSD model
for(i in 2:8){
  #all lowercase
  EU_sites[,i] <- tolower(EU_sites[,i])
  sweden_sites[,i] <- tolower(sweden_sites[,i])
}
sweden_sites <- sweden_sites[,-13]
for(i in 9:15){
  #all lowercase
  sweden_sites[,i] <- tolower(sweden_sites[,i])
}

sweden_sites <- sweden_sites[!is.na(sweden_sites$kingdom),]

#sort out issues with gastropods
EU_sites[EU_sites$class == "gastropoda" & EU_sites$order == "",]$order <-  paste("unknown (",   EU_sites[EU_sites$class == "gastropoda" & EU_sites$order == "",]$family,  " family)", sep = "")
sweden_sites[sweden_sites$class == "gastropoda" & is.na(sweden_sites$order) == "",]$order <-  paste("unknown (",   sweden_sites[sweden_sites$class == "gastropoda" & sweden_sites$order == "",]$family,  " family)", sep = "")

sweden_sites[sweden_sites$order == "tricladida" & is.na(sweden_sites$class),]$class <-  "turbellaria"
sweden_sites[sweden_sites$order == "prorhynchida" & is.na(sweden_sites$class) &  !is.na(sweden_sites$order),]$class <-  "prorhynchida"

  #fill in missing lower taxonomy with the higher taxonomic classification
for(i in 7:2){
  EU_sites[EU_sites[,i] == "",2:i] <- paste("unknown (",   EU_sites[EU_sites[,i] == "",i+1], " " ,colnames(EU_sites)[i+1],  ")", sep = "")
}

EU_sites <- EU_sites[EU_sites$species != "unknown (arthropoda phylum)",]

for(i in 14:9){
  sweden_sites[is.na(sweden_sites[,i]),9:i] <- paste("unknown (",   sweden_sites[is.na(sweden_sites[,i]),i+1], " " ,colnames(sweden_sites)[i+1],  ")", sep = "")
}

sweden_sites <- sweden_sites[sweden_sites$species != "unknown (arthropoda phylum)",]

#make a unique list of taxa from the taxonomic data
unique_taxa <- rbind(EU_sites[,2:7], sweden_sites[,9:14])

#remove duplicates present at multiple sites
unique_taxa <- unique_taxa[!duplicated(unique_taxa),]

#reorder and rename to be the same as the hSSD model
hSSD_taxa <- data.frame(Kingdom = "animalia", Phylum_division = unique_taxa[,6],Subphylum = NA, Superclass = NA, Class = unique_taxa[,5], Order = unique_taxa[,4], Family = unique_taxa[,3], Genus = unique_taxa[,2], Latin = unique_taxa[,1])

#####
hSSD_taxa[hSSD_taxa$Phylum_division == "nematoda",]$Subphylum <- "unknown (nematoda phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "nematoda",]$Superclass <- "unknown (nematoda phylum)"
#nemata are nematodes
hSSD_taxa[hSSD_taxa$Phylum_division == "nematoda",]$Phylum_division <- "nemata"

hSSD_taxa[hSSD_taxa$Phylum_division == "platyhelminthes",]$Subphylum <- "unknown (platyhelminthes phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "platyhelminthes",]$Superclass <- "unknown (platyhelminthes phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "nematomorpha",]$Subphylum <- "unknown (nematomorpha phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "nematomorpha",]$Superclass <- "unknown (nematomorpha phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "porifera",]$Subphylum <- "unknown (porifera phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "porifera",]$Superclass <- "unknown (porifera phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "bryozoa",]$Subphylum <- "unknown (bryozoa phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "bryozoa",]$Superclass <- "unknown (bryozoa phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "acanthocephala",]$Subphylum <- "unknown (acanthocephala phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "acanthocephala",]$Superclass <- "unknown (acanthocephala phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "nemertea",]$Subphylum <- "unknown (nemertea phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "nemertea",]$Superclass <- "unknown (nemertea phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "cnidaria",]$Subphylum <- "unknown (cnidaria phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "cnidaria",]$Superclass <- "unknown (cnidaria phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "rotifera",]$Subphylum <- "unknown (rotifera phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "rotifera",]$Superclass <- "unknown (rotifera phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "kamptozoa",]$Subphylum <- "unknown (kamptozoa phylum)"
hSSD_taxa[hSSD_taxa$Phylum_division == "kamptozoa",]$Superclass <- "unknown (kamptozoa phylum)"
#collembola not in database at all
hSSD_taxa[hSSD_taxa$Class == "collembola",]$Subphylum <- "unknown (collembola class)"
hSSD_taxa[hSSD_taxa$Class == "collembola",]$Superclass <- "unknown (collembola class)"
#adenophorea nematodes missing
hSSD_taxa[hSSD_taxa$Class == "adenophorea",]$Subphylum <- "Missing Subphylum 1"
hSSD_taxa[hSSD_taxa$Class == "adenophorea",]$Superclass <- "Missing Superclass 1"
#taxonomy of cnidaria is different missing subclass rather than blank?
hSSD_taxa[hSSD_taxa$Class == "anthozoa",]$Subphylum <- "Missing Subphylum 2"
hSSD_taxa[hSSD_taxa$Class == "anthozoa",]$Superclass <- "Missing Superclass 3"
#taxonomy of cnidaria is missing subclass rather than blank?
hSSD_taxa[hSSD_taxa$Class == "hydrozoa",]$Subphylum <- "Missing Subphylum 9"
hSSD_taxa[hSSD_taxa$Class == "hydrozoa",]$Superclass <- "Missing Superclass 15"
#higher taxonomy of turbellaria is missing subclass rather than blank?
hSSD_taxa[hSSD_taxa$Class == "turbellaria",]$Subphylum <- "Missing Subphylum 16"
hSSD_taxa[hSSD_taxa$Class == "turbellaria",]$Superclass <- "Missing Superclass 26"
#plumatellida are bryozoa not ectoprocta
hSSD_taxa[hSSD_taxa$Order == "plumatellida",]$Phylum_division <- "ectoprocta"
hSSD_taxa[hSSD_taxa$Order == "plumatellida",]$Subphylum <- "Missing Subphylum 12"
hSSD_taxa[hSSD_taxa$Order == "plumatellida",]$Superclass <- "Missing Superclass 22"
hSSD_taxa[hSSD_taxa$Order == "plumatellida",]$Class <- "phylactolaemata"
#oligochaets are a subclass of clitellata
hSSD_taxa[hSSD_taxa$Order == "haplotaxida",]$Class <- "oligochaeta"
#spelling mistake
hSSD_taxa[hSSD_taxa$Order == "unionida",]$Order <- "unionoida"
#spelling mistake
hSSD_taxa[hSSD_taxa$Order == "venerida",]$Order <- "veneroida"
#pisidiidae is now it's own family of sphaeriidae and order of sphaeriida 
hSSD_taxa[hSSD_taxa$Genus == "sphaerium",]$Family <- "pisidiidae"
hSSD_taxa[hSSD_taxa$Family == "sphaeriidae",]$Family <- "pisidiidae"
hSSD_taxa[hSSD_taxa$Family == "pisidiidae",]$Order <- "veneroida"
#hydridae is part of the anthoathecata order
hSSD_taxa[hSSD_taxa$Family == "hydridae",]$Order <- "hydroida"
#heteroptera is a suborder, actual hemiptera
hSSD_taxa[hSSD_taxa$Family == "gerridae",]$Order <- "heteroptera"
hSSD_taxa[hSSD_taxa$Family == "corixidae",]$Order <- "heteroptera"
hSSD_taxa[hSSD_taxa$Family == "naucoridae",]$Order <- "heteroptera"
hSSD_taxa[hSSD_taxa$Family == "notonectidae",]$Order <- "heteroptera"
#aphelocheirus is in aphelocheiridae family and hemiptera
hSSD_taxa[hSSD_taxa$Family == "aphelocheiridae",]$Order <- "heteroptera"
hSSD_taxa[hSSD_taxa$Genus == "aphelocheirus",]$Family <- "naucoridae"
#snails
hSSD_taxa[hSSD_taxa$Family == "planorbidae",]$Order <- "basommatophora"
hSSD_taxa[hSSD_taxa$Family == "lymnaeidae",]$Order <- "basommatophora"
hSSD_taxa[hSSD_taxa$Family == "physidae",]$Order <- "basommatophora"
#ancylus is actually a planorbidae 
hSSD_taxa[hSSD_taxa$Genus == "ancylus",]$Family <- "ancylidae"
#snails
hSSD_taxa[hSSD_taxa$Family == "hydrobiidae",]$Order <- "neotaenioglossa"
hSSD_taxa[hSSD_taxa$Family == "thiaridae",]$Order <- "neotaenioglossa"
hSSD_taxa[hSSD_taxa$Family == "calyptraeidae",]$Order <- "neotaenioglossa"
#snails, potamopyrgus is a littorininimorpha order and tateidae  family
hSSD_taxa[hSSD_taxa$Family == "tateidae",]$Order <- "neotaenioglossa"
hSSD_taxa[hSSD_taxa$Genus == "potamopyrgus",]$Family <- "hydrobiidae"
#snails
hSSD_taxa[hSSD_taxa$Family == "melanopsidae",]$Order <- "mesogastropoda"
#zygoptera is a suborder of odonata
hSSD_taxa[hSSD_taxa$Family == "lestidae",]$Order <- "zygoptera"
#mussel is in order myida 
hSSD_taxa[hSSD_taxa$Family == "dreissenidae",]$Order <- "veneroida"
#specific reclassification of mussels
hSSD_taxa[hSSD_taxa$Genus == "corbicula",]$Family <- "corbiculidae"
#tubificidae now split into naididae
hSSD_taxa[hSSD_taxa$Genus == "tubifex",]$Family <- "tubificidae"
hSSD_taxa[hSSD_taxa$Genus == "limnodrilus",]$Family <- "tubificidae"
hSSD_taxa[hSSD_taxa$Genus == "spirosperma",]$Family <- "tubificidae"
hSSD_taxa[hSSD_taxa$Genus == "rhyacodrilus",]$Family <- "tubificidae"
hSSD_taxa[hSSD_taxa$Genus == "branchiura",]$Family <- "tubificidae"
#naididaeare tubificida , oligochaeta
hSSD_taxa[hSSD_taxa$Family == "naididae",]$Order <- "haplotaxida"
#tubificidae are tubificida , oligochaeta
hSSD_taxa[hSSD_taxa$Family == "tubificidae",]$Order <- "haplotaxida"
hSSD_taxa[hSSD_taxa$Order == "haplotaxida",]$Class <- "oligochaeta"
#tricladida are rhabditophora not turbellaria
hSSD_taxa[!is.na(hSSD_taxa$Order) & hSSD_taxa$Order == "tricladida",]$Class <- "turbellaria"
#valvatidae are not classified at order
hSSD_taxa[hSSD_taxa$Family == "valvatidae",]$Order <- "lower heterobranchia"
#aeolosomatidae is an ol,igochaet
hSSD_taxa[hSSD_taxa$Family == "aeolosomatidae",]$Order <- "crassiclitellata"
hSSD_taxa[hSSD_taxa$Family == "aeolosomatidae",]$Class <- "clitellata"
#tachinidae is a chironomid
hSSD_taxa[hSSD_taxa$Genus == "cricotopus",]$Family <- "chironomidae"
hSSD_taxa[hSSD_taxa$Genus == "microtendipes",]$Family <- "chironomidae"

hSSD_taxa[hSSD_taxa$Genus == "ilyodrilus",]$Family <- "tubificidae"
#####
write.csv(hSSD_taxa, here("output", "unique_taxa_summer.csv"), row.names = FALSE)

#reformat the sample data into hSSD format
EU_taxa <- data.frame(sample = EU_sites$sample_id, Kingdom = "animalia", Phylum_division = EU_sites[,7],Subphylum = NA, Superclass = NA, Class = EU_sites[,6], Order = EU_sites[,5], Family = EU_sites[,4], Genus = EU_sites[,3], Latin = EU_sites[,2])

#####
EU_taxa[EU_taxa$Phylum_division == "nematoda",]$Subphylum <- "unknown (nematoda phylum)"
EU_taxa[EU_taxa$Phylum_division == "nematoda",]$Superclass <- "unknown (nematoda phylum)"
#nemata are nematodes
EU_taxa[EU_taxa$Phylum_division == "nematoda",]$Phylum_division <- "nemata"

EU_taxa[EU_taxa$Phylum_division == "platyhelminthes",]$Subphylum <- "unknown (platyhelminthes phylum)"
EU_taxa[EU_taxa$Phylum_division == "platyhelminthes",]$Superclass <- "unknown (platyhelminthes phylum)"
EU_taxa[EU_taxa$Phylum_division == "nematomorpha",]$Subphylum <- "unknown (nematomorpha phylum)"
EU_taxa[EU_taxa$Phylum_division == "nematomorpha",]$Superclass <- "unknown (nematomorpha phylum)"
EU_taxa[EU_taxa$Phylum_division == "porifera",]$Subphylum <- "unknown (porifera phylum)"
EU_taxa[EU_taxa$Phylum_division == "porifera",]$Superclass <- "unknown (porifera phylum)"
EU_taxa[EU_taxa$Phylum_division == "bryozoa",]$Subphylum <- "unknown (bryozoa phylum)"
EU_taxa[EU_taxa$Phylum_division == "bryozoa",]$Superclass <- "unknown (bryozoa phylum)"
EU_taxa[EU_taxa$Phylum_division == "acanthocephala",]$Subphylum <- "unknown (acanthocephala phylum)"
EU_taxa[EU_taxa$Phylum_division == "acanthocephala",]$Superclass <- "unknown (acanthocephala phylum)"
EU_taxa[EU_taxa$Phylum_division == "nemertea",]$Subphylum <- "unknown (nemertea phylum)"
EU_taxa[EU_taxa$Phylum_division == "nemertea",]$Superclass <- "unknown (nemertea phylum)"
EU_taxa[EU_taxa$Phylum_division == "cnidaria",]$Subphylum <- "unknown (cnidaria phylum)"
EU_taxa[EU_taxa$Phylum_division == "cnidaria",]$Superclass <- "unknown (cnidaria phylum)"
EU_taxa[EU_taxa$Phylum_division == "rotifera",]$Subphylum <- "unknown (rotifera phylum)"
EU_taxa[EU_taxa$Phylum_division == "rotifera",]$Superclass <- "unknown (rotifera phylum)"
EU_taxa[EU_taxa$Phylum_division == "kamptozoa",]$Subphylum <- "unknown (kamptozoa phylum)"
EU_taxa[EU_taxa$Phylum_division == "kamptozoa",]$Superclass <- "unknown (kamptozoa phylum)"
#collembola not in database at all
EU_taxa[EU_taxa$Class == "collembola",]$Subphylum <- "unknown (collembola class)"
EU_taxa[EU_taxa$Class == "collembola",]$Superclass <- "unknown (collembola class)"
#adenophorea nematodes missing
EU_taxa[EU_taxa$Class == "adenophorea",]$Subphylum <- "Missing Subphylum 1"
EU_taxa[EU_taxa$Class == "adenophorea",]$Superclass <- "Missing Superclass 1"
#taxonomy of cnidaria is different missing subclass rather than blank?
EU_taxa[EU_taxa$Class == "anthozoa",]$Subphylum <- "Missing Subphylum 2"
EU_taxa[EU_taxa$Class == "anthozoa",]$Superclass <- "Missing Superclass 3"
#taxonomy of cnidaria is missing subclass rather than blank?
EU_taxa[EU_taxa$Class == "hydrozoa",]$Subphylum <- "Missing Subphylum 9"
EU_taxa[EU_taxa$Class == "hydrozoa",]$Superclass <- "Missing Superclass 15"
#higher taxonomy of turbellaria is missing subclass rather than blank?
EU_taxa[EU_taxa$Class == "turbellaria",]$Subphylum <- "Missing Subphylum 16"
EU_taxa[EU_taxa$Class == "turbellaria",]$Superclass <- "Missing Superclass 26"
#plumatellida are bryozoa not ectoprocta
EU_taxa[EU_taxa$Order == "plumatellida",]$Phylum_division <- "ectoprocta"
EU_taxa[EU_taxa$Order == "plumatellida",]$Subphylum <- "Missing Subphylum 12"
EU_taxa[EU_taxa$Order == "plumatellida",]$Superclass <- "Missing Superclass 22"
EU_taxa[EU_taxa$Order == "plumatellida",]$Class <- "phylactolaemata"
#oligochaets are a subclass of clitellata
EU_taxa[EU_taxa$Order == "haplotaxida",]$Class <- "oligochaeta"
#spelling mistake
EU_taxa[EU_taxa$Order == "unionida",]$Order <- "unionoida"
#spelling mistake
EU_taxa[EU_taxa$Order == "venerida",]$Order <- "veneroida"
#pisidiidae is now it's own family of sphaeriidae and order of sphaeriida 
EU_taxa[EU_taxa$Genus == "sphaerium",]$Family <- "pisidiidae"
EU_taxa[EU_taxa$Family == "sphaeriidae",]$Family <- "pisidiidae"
EU_taxa[EU_taxa$Family == "pisidiidae",]$Order <- "veneroida"
#hydridae is part of the anthoathecata order
EU_taxa[EU_taxa$Family == "hydridae",]$Order <- "hydroida"
#heteroptera is a suborder, actual hemiptera
EU_taxa[EU_taxa$Family == "gerridae",]$Order <- "heteroptera"
EU_taxa[EU_taxa$Family == "corixidae",]$Order <- "heteroptera"
EU_taxa[EU_taxa$Family == "naucoridae",]$Order <- "heteroptera"
EU_taxa[EU_taxa$Family == "notonectidae",]$Order <- "heteroptera"
#aphelocheirus is in aphelocheiridae family and hemiptera
EU_taxa[EU_taxa$Family == "aphelocheiridae",]$Order <- "heteroptera"
EU_taxa[EU_taxa$Genus == "aphelocheirus",]$Family <- "naucoridae"
#snails
EU_taxa[EU_taxa$Family == "planorbidae",]$Order <- "basommatophora"
EU_taxa[EU_taxa$Family == "lymnaeidae",]$Order <- "basommatophora"
EU_taxa[EU_taxa$Family == "physidae",]$Order <- "basommatophora"
#ancylus is actually a planorbidae 
EU_taxa[EU_taxa$Genus == "ancylus",]$Family <- "ancylidae"
#snails
EU_taxa[EU_taxa$Family == "hydrobiidae",]$Order <- "neotaenioglossa"
EU_taxa[EU_taxa$Family == "thiaridae",]$Order <- "neotaenioglossa"
EU_taxa[EU_taxa$Family == "calyptraeidae",]$Order <- "neotaenioglossa"
#snails, potamopyrgus is a littorininimorpha order and tateidae  family
EU_taxa[EU_taxa$Family == "tateidae",]$Order <- "neotaenioglossa"
EU_taxa[EU_taxa$Genus == "potamopyrgus",]$Family <- "hydrobiidae"
#snails
EU_taxa[EU_taxa$Family == "melanopsidae",]$Order <- "mesogastropoda"
#zygoptera is a suborder of odonata
EU_taxa[EU_taxa$Family == "lestidae",]$Order <- "zygoptera"
#mussel is in order myida 
EU_taxa[EU_taxa$Family == "dreissenidae",]$Order <- "veneroida"
#specific reclassification of mussels
EU_taxa[EU_taxa$Genus == "corbicula",]$Family <- "corbiculidae"
#tubificidae now split into naididae
EU_taxa[EU_taxa$Genus == "tubifex",]$Family <- "tubificidae"
EU_taxa[EU_taxa$Genus == "limnodrilus",]$Family <- "tubificidae"
EU_taxa[EU_taxa$Genus == "spirosperma",]$Family <- "tubificidae"
EU_taxa[EU_taxa$Genus == "rhyacodrilus",]$Family <- "tubificidae"
EU_taxa[EU_taxa$Genus == "branchiura",]$Family <- "tubificidae"
#naididaeare tubificida , oligochaeta
EU_taxa[EU_taxa$Family == "naididae",]$Order <- "haplotaxida"
#tubificidae are tubificida , oligochaeta
EU_taxa[EU_taxa$Family == "tubificidae",]$Order <- "haplotaxida"
EU_taxa[EU_taxa$Order == "haplotaxida",]$Class <- "oligochaeta"
#tricladida are rhabditophora not turbellaria
EU_taxa[!is.na(EU_taxa$Order) & EU_taxa$Order == "tricladida",]$Class <- "turbellaria"
#valvatidae are not classified at order
EU_taxa[EU_taxa$Family == "valvatidae",]$Order <- "lower heterobranchia"
#aeolosomatidae is an ol,igochaet
EU_taxa[EU_taxa$Family == "aeolosomatidae",]$Order <- "crassiclitellata"
EU_taxa[EU_taxa$Family == "aeolosomatidae",]$Class <- "clitellata"
#tachinidae is a chironomid
EU_taxa[EU_taxa$Genus == "cricotopus",]$Family <- "chironomidae"
EU_taxa[EU_taxa$Genus == "microtendipes",]$Family <- "chironomidae"
#
EU_taxa[EU_taxa$Genus == "ilyodrilus",]$Family <- "tubificidae"
#####
EU_taxa <- EU_taxa[!is.na(EU_taxa$sample),]

#make a list of all the sites
sample_ids <- unique(EU_taxa$sample)

for(i in 1:length(sample_ids)){
 
  print(i)
  
  sample_taxa <- EU_taxa[EU_taxa$sample == sample_ids[i],]
  
  #remove duplicatte species (or latin names)
  sample_taxa <- sample_taxa[!duplicated(sample_taxa$Latin),]

  #there is still a problem with nested taxonomic ranks
 #if there are multiple taxa with the same genus... 
if(any(duplicated(sample_taxa$Genus))){
  #record the duplicates
  duplicates <- sample_taxa[sample_taxa$Genus %in% sample_taxa[duplicated(sample_taxa$Genus),]$Genus,]$Latin
  #and remove any generic ones at higher taxonomic ranks ie. those that are unknown
  sample_taxa <- sample_taxa[!(sample_taxa$Latin %in% duplicates[grepl(" genus)",duplicates)]),]
}
  #repeat for multiple taxa with the same family... 
  if(any(duplicated(sample_taxa$Family))){
    #record the duplicates
    duplicates <- sample_taxa[sample_taxa$Family %in% sample_taxa[duplicated(sample_taxa$Family),]$Family,]$Latin
    #and remove any generic ones at higher taxonomic ranks ie. those that are unknown
    sample_taxa <- sample_taxa[!(sample_taxa$Latin %in% duplicates[grepl(" family)",duplicates)]),]
  }
  
  write.csv(sample_taxa[,2:length(sample_taxa)], here("output", "sample", paste(sample_ids[i], "_tax.csv", sep = "")), row.names = FALSE)
}
