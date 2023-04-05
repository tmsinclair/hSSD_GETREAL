library(here)
library(tidyverse)
library(stringi)

chems <- read_csv(here("input","chems","names.csv"))

files <- list.files(here("output","sample"), pattern = ".csv")

for(i in 1:nrow(chems)){
  
  nearest <- data.frame(chem = NA, name = NA, sample = NA, nearest_taxa = NA, nearest_taxa_rank = NA, nearest_n = NA, nearest_rank = NA)
  
  tox_tax <- read.csv(here("input", "tox", paste(chems$name[i], "_tox_tax.csv", sep = "")))
  
  for(j in 1:length(files)){
    print(paste(chems$name[i], j, "of", length(files)))
    sample_tax <- read.csv(here("output", "sample", files[j]))
    
    nearest_sub <- data.frame(chem = NA, name = NA, sample = NA, nearest_taxa = NA, nearest_taxa_rank = NA, nearest_n = NA, nearest_rank = NA)
    
    for(k in 1:nrow(sample_tax)){
      for(l in (length(sample_tax)-1):1){
        if(any(na.omit(sample_tax[k,l]) %in% tox_tax[,l])){
          #print(paste("match at", colnames(sample_tax)[l]))
          nearest_sub[nrow(nearest_sub)+1,] <- c(chems$abbreviation[i], chems$name[i], gsub("_tax.csv", "", files[j]), 
                                         tox_tax[tox_tax[,l] %in% na.omit(sample_tax[k,l]),]$Latin[1],
                                         tox_tax[tox_tax[,l] %in% na.omit(sample_tax[k,l]),][1,l], 
                                         nrow(tox_tax[tox_tax[,l] %in% na.omit(sample_tax[k,l]),]), 
                                         colnames(sample_tax)[l])
          break
        }
      }
    }
    nearest <- rbind(nearest, nearest_sub[-1,])
  }
  nearest <- nearest[-1,]
  write.csv(nearest, here("output", paste("nearest_tox_data_", chems$abbreviation[i], ".csv", sep = "")), row.names = FALSE)
}


nearest <- lapply(here("output",list.files(here("output"), pattern = "nearest_tox_data")), read.csv)
nearest <- data.frame(Reduce(rbind, nearest))

nearest[nearest$nearest_rank == "Phylum_division",]$nearest_rank <- "Phylum" 

nearest$nearest_rank <- factor(nearest$nearest_rank, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))

nearest$nearest_taxa_rank <- stri_trans_totitle(nearest$nearest_taxa_rank)

nearest_plot <- ggplot(data = nearest, aes(x = nearest_rank, fill = chem))+
  geom_bar(stat = "count", col = "#666666")+
  labs(x = "Nearest taxonomic rank with toxicity data", y = "Number of predicted sensitivity values")+
  scale_fill_manual(values = c("#a0d699","#38bbf1","#f0eaba", "#d69b34"))+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~name)
nearest_plot
ggsave(here("figures","nearest_taxonomic_rank.png"), nearest_plot, width = 12, height = 8,  bg = "transparent")

nearest$nearest <- paste(nearest$nearest_taxa_rank,nearest$nearest_rank)

nearest_count <- nearest |>
  group_by(name, chem, nearest) |>
  summarise(
    n = n()
  ) 

nearest_count$prop <- nearest_count$n / nrow(filter(nearest, chem == "ATZ"))

nearest_count$percent <- nearest_count$prop * 100

nearest_plot <- ggplot(data = filter(nearest_count, percent > 1 & chem == "ATZ"), aes(x = reorder(nearest, -percent), y = percent))+
  geom_bar(stat = "identity", col = "#666666", fill  = "#a0d699")+
  labs(x = "Nearest toxicity data", y = "Percent of predicted values (%)", title = chems$name[1])+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20, angle = 90),
    axis.title.x = element_text(angle = 180),
    plot.title = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))
nearest_plot
ggsave(here("figures","nearest_taxonomic_data_ATZ.png"), nearest_plot, width = 8, height = 12,  bg = "transparent")


nearest_plot <- ggplot(data = filter(nearest_count, percent > 1 & chem == "Cu"), aes(x = reorder(nearest, -percent), y = percent))+
  geom_bar(stat = "identity", col = "#666666", fill  = "#38bbf1")+
  labs(x = "Nearest toxicity data", y = "Percent of predicted values (%)", title = chems$name[2])+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20, angle = 90),
    axis.title.x = element_text(angle = 180),
    plot.title = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))
nearest_plot
ggsave(here("figures","nearest_taxonomic_data_Cu.png"), nearest_plot, width = 8, height = 12,  bg = "transparent")


nearest_plot <- ggplot(data = filter(nearest_count, percent > 1 & chem == "IMD"), aes(x = reorder(nearest, -percent), y = percent))+
  geom_bar(stat = "identity", col = "#666666", fill  = "#f0eaba")+
  labs(x = "Nearest toxicity data", y = "Percent of predicted values (%)", title = chems$name[3])+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20, angle = 90),
    axis.title.x = element_text(angle = 180),
    plot.title = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))
nearest_plot
ggsave(here("figures","nearest_taxonomic_data_IMD.png"), nearest_plot, width = 8, height = 12,  bg = "transparent")


nearest_plot <- ggplot(data = filter(nearest_count, percent > 1 & chem == "LCH"), aes(x = reorder(nearest, -percent), y = percent))+
  geom_bar(stat = "identity", col = "#666666", fill  = "#d69b34")+
  labs(x = "Nearest toxicity data", y = "Percent of predicted values (%)", title = chems$name[4])+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 20, angle = 90),
    axis.title.x = element_text(angle = 180),
    plot.title = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))
nearest_plot
ggsave(here("figures","nearest_taxonomic_data_LCH.png"), nearest_plot, width = 8, height = 12,  bg = "transparent")
