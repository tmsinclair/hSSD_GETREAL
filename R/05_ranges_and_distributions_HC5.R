library(here)
library(tidyverse)
library(agricolae)
library(stringi)

data_summary <- function(x) {
  m <- mean(x, na.rm=TRUE)
  sem <-sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))
  ymin<-m-1.96*sem
  ymax<-m+1.96*sem
  return(c(y=m,ymin=ymin,ymax=ymax))
}

tox <- read.csv(here("output","hSSD_GETREAL_all_chems_properties.csv"))

#readin in the chemicals
chems <- read_csv(here("input","chems","names.csv"))

tox$name <- NA


for(i in 1:nrow(chems)){
  tox[tox$chem == chems$abbreviation[i],]$name <- chems$name[i]
}

sum_ta7 <- tox %>%
  filter(name == "Atrazine" & !is.na(ta7)) |>
  group_by(ta7) %>%
  summarise(
    n = length(ta7)
  )
sum_ta7
  
TA_desc <- c("Very large",
             "Small, lowland, silicious",
             "Large, lowland, silicious",
             "Lowland calcareous",
             "Large mid-altitude",
             "Small mid-altitude",
             "High altitude")

sum_ta7$ta_desc <- TA_desc

sum_ta7 <- sum_ta7 %>%
  mutate(ta_desc = fct_relevel(ta_desc, 
                           "Very large",
                           "Small, lowland, silicious",
                           "Large, lowland, silicious",
                           "Lowland calcareous",
                           "Large mid-altitude",
                           "Small mid-altitude",
                           "High altitude")) 

summary <- ggplot(data = sum_ta7, aes(x = ta_desc, y = n))+
  geom_bar(fill = "#BBBBBb", col = "#999999", stat = "identity")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  geom_text(aes(x = ta_desc, y =  n + 20,
                                           label = n), size = 8, col = "#222222")+
  labs(y = "Number of assemblages")+
  theme(
    axis.ticks.x=element_blank(),
    text = element_text(size=20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    axis.text.x = element_text(size=24),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#EEEEEE"))
summary

ggsave(here("figures","HC5_TA7_summary_chem_bar.png"), summary, width = 14, height = 8,  bg = "transparent")

sum_all <- tox %>%
  group_by(chem, name) %>%
  summarise(
    mean = mean(HC5_aj),
    var = var(HC5_aj),
    max = max(HC5_aj),
    min = min(HC5_aj),
    range = max(HC5_aj) - min(HC5_aj),
    range_prop = max(HC5_aj) / min(HC5_aj),
    lower = quantile(HC5_aj, probs = 0.025),
    middle = quantile(HC5_aj, probs = 0.5),
    upper = quantile(HC5_aj, probs = 0.975),
    quartile = (quantile(HC5_aj, probs = 0.75) - quantile(HC5_aj, probs = 0.25)) / mean(HC5_aj),
    CI = upper - lower,
    CI_prop = upper / lower
  )
sum_all$CI_prop
sum <- tox %>%
  group_by(chem, name) %>%
  summarise(
    mean = mean(HC5_aj),
    var = var(HC5_aj),
    max = max(HC5_aj),
    min = min(HC5_aj),
    CI_prop = max(HC5_aj) - min(HC5_aj),
    CI_prop = max(HC5_aj) / min(HC5_aj),
    lower = quantile(HC5_aj, probs = 0.025),
    middle = quantile(HC5_aj, probs = 0.5),
    upper = quantile(HC5_aj, probs = 0.975),
    quartile = (quantile(HC5_aj, probs = 0.75) - quantile(HC5_aj, probs = 0.25)) / mean(HC5_aj),
    CI = upper - lower,
    CI_prop = upper / lower
  )


sum$CI_prop <- sum$CI_prop / sum$mean

tox$HC5_adj <- NA

for(i in 1:nrow(chems)){
  tox[tox$chem == chems$abbreviation[i],]$HC5_adj <- tox[tox$chem == chems$abbreviation[i],]$HC5_aj / sum[sum$chem == chems$abbreviation[i],]$middle
}

sum <- tox %>%
  group_by(chem, name, ta7) %>%
  summarise(
    mean = mean(HC5_aj),
    var = var(HC5_aj),
    max = max(HC5_aj),
    min = min(HC5_aj),
    CI_prop = max(HC5_aj) - min(HC5_aj),
    CI_prop = max(HC5_aj) / min(HC5_aj),
    sd = sd(HC5_aj),
    length = length(HC5_aj),
    sem = sd/sqrt(length),
    lower = quantile(HC5_aj, probs = 0.025),
    middle = quantile(HC5_aj, probs = 0.5),
    upper = quantile(HC5_aj, probs = 0.975),
    quartile = (quantile(HC5_aj, probs = 0.75) - quantile(HC5_aj, probs = 0.25)) / mean(HC5_aj),
    CI = upper - lower,
    CI_prop = upper / lower
  )



summary <- ggplot(data = tox, aes(y = log10(HC5_aj)))+
  geom_boxplot()+
  scale_fill_manual(values = c("#cc6752","#c3cce3","#85a3a6"), name = "Chemical type")+
  labs(x = "Chemical", y = "log[Predicted HC5] (\u00B5g/L)")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=16),
        legend.position = "bottom",
        panel.background = element_rect(fill = "#FFFFFF"),
        axis.line = element_line(colour = "#333333"),
        strip.background = element_rect(fill = "#FFFFFF"))+
  facet_wrap(~name,  strip.position = "bottom")
summary

ggsave(here("figures","HC5_site_summary.png"), summary, width = 12, height = 8,  bg = "transparent")

#####
#Typologies/ assemblage types
tox <- tox[!is.na(tox$ta7),]
#hc5


#points not box plots
#sum$ta7 <- factor(sum$ta7, levels = c("RT1", "RT2", "RT3", "RT4_5", "RT6", "RT7", "RT8_10_11_18", "RT9", "RT12", "RT14_15_16", "RT17", "RT19", "RT20"))


#same but for CI_props
#sum$ta7 <- factor(sum$ta7, levels = c("RT1", "RT2", "RT3", "RT4_5", "RT6", "RT7", "RT8_10_11_18", "RT9", "RT12", "RT14_15_16", "RT17", "RT19", "RT20"))
sum <- sum[!is.na(sum$ta7),]


#####
#ANOVAs
#create summary
#sum[,2:11] <- log10(sum[,2:11])
#Mean tox
#tox$ta7 <- factor(tox$ta7, levels = c("RT1", "RT2", "RT3", "RT4_5", "RT6", "RT7", "RT8_10_11_18", "RT9", "RT12", "RT14_15_16", "RT17", "RT19", "RT20"))
tox <- tox[!is.na(tox$ta7),]

anova_middle <- glm(tox$HC5_adj ~ tox$chem * tox$ta7)
summary(anova_middle)

tox$treatment <- paste(tox$ta7, tox$chem)

HSD_AB_middle <- HSD.test(y = tox$HC5_adj,
                          trt =  tox$treatment,
                          DFerror = anova_middle$df.residual,
                          MSerror = deviance(anova_middle)/anova_middle$df.residual,
                          alpha = 0.05)
HSD_AB_middle

HSD_AB_middle <- cbind(HSD_AB_middle$groups,HSD_AB_middle$means$Max[order(HSD_AB_middle$means$`tox$HC5_adj`, decreasing = TRUE)])
HSD_AB_middle$treatment <- row.names(HSD_AB_middle)
colnames(HSD_AB_middle)[c(1,3)] <- c("param","max")
HSD_AB_middle

HSD_AB_middle$chem <- gsub(".* ", "", HSD_AB_middle$treatment)
HSD_AB_middle$ta7 <- gsub(" .*", "", HSD_AB_middle$treatment)

middle <-  ggplot(data = tox, aes(x = ta7, y = log10(HC5_adj), fill = chem))+
  geom_boxplot()+
  geom_text(data = HSD_AB_middle, aes(x = ta7, y = 2, label =  groups), size = 4)+
  scale_fill_manual(values = c("#98d191","#38bbf1","#faf5cc","#d69b34"), name = "")+
  labs(x = "", y = "log10(Mean predicted HC5) (ug/l)")+
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~chem)

middle

ggsave(here("figures","HC5_sample_summary_ANOVA_chem.png"), middle, width = 20, height = 8,  bg = "transparent")


middle <-  ggplot(data = tox, aes(x = chem, y = log10(HC5_adj), fill = chem))+
  geom_boxplot()+
  
  geom_text(data = HSD_AB_middle, aes(x = chem, y = 1.5, label =  groups), size = 6)+
  scale_fill_manual(values = c("#98d191","#38bbf1","#faf5cc","#d69b34"), name = "")+
  labs(x = "", y = "log10(Mean predicted HC5) (ug/l)")+
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~ta7)

middle

ggsave(here("figures","HC5_sample_summary_ANOVA_ta7.png"), middle, width = 12, height = 8,  bg = "transparent")


shapiro.test(log(sum$CI_prop))

#CI_prop tox
anova_middle <- glm(sum$HC5_adj ~ tox$chem * tox$ta7)
summary(anova_middle)

tox$treatment <- paste(tox$ta7, tox$chem)

HSD_AB_middle <- HSD.test(y = tox$HC5_adj,
                          trt =  tox$treatment,
                          DFerror = anova_middle$df.residual,
                          MSerror = deviance(anova_middle)/anova_middle$df.residual,
                          alpha = 0.05)
HSD_AB_middle

HSD_AB_middle <- cbind(HSD_AB_middle$groups,HSD_AB_middle$means$Max[order(HSD_AB_middle$means$`tox$HC5_adj`, decreasing = TRUE)])
HSD_AB_middle$treatment <- row.names(HSD_AB_middle)
colnames(HSD_AB_middle)[c(1,3)] <- c("param","max")
HSD_AB_middle

HSD_AB_middle$chem <- gsub(".* ", "", HSD_AB_middle$treatment)
HSD_AB_middle$ta7 <- gsub(" .*", "", HSD_AB_middle$treatment)

middle <-  ggplot(data = tox, aes(x = ta7, y = log10(HC5_adj), fill = chem))+
  geom_boxplot()+
  geom_text(data = HSD_AB_middle, aes(x = ta7, y = 2, label =  groups), size = 4)+
  scale_fill_manual(values = c("#98d191","#38bbf1","#faf5cc","#d69b34"), name = "")+
  labs(x = "", y = "log10(Mean predicted HC5) (ug/l)")+
  theme(
    text = element_text(size = 12),
    axis.text = element_text(size = 12, angle = -20),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~chem)

middle

ggsave(here("figures","HC5_sample_CI_prop_ANOVA_chem.png"), middle, width = 20, height = 8,  bg = "transparent")


middle <-  ggplot(data = tox, aes(x = chem, y = log10(HC5_adj), fill = chem))+
  geom_boxplot()+
  
  geom_text(data = HSD_AB_middle, aes(x = chem, y = 1.5, label =  groups), size = 6)+
  scale_fill_manual(values = c("#98d191","#38bbf1","#faf5cc","#d69b34"), name = "")+
  labs(x = "", y = "log10(Mean predicted HC5) (ug/l)")+
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~ta7)

middle

ggsave(here("figures","HC5_sample_CI_prop_ANOVA_ta7.png"), middle, width = 12, height = 8,  bg = "transparent")

#the same but with bar charts
HSD_AB_middle_ta7 <- NULL
for(j in 1:length(unique(tox$ta7))){
  print(unique(tox$ta7)[j])
  anova_middle_ta7_sep <- glm(tox[tox$ta7 == unique(tox$ta7)[j],]$HC5_adj ~ tox[tox$ta7 == unique(tox$ta7)[j],]$chem)
summary(anova_middle_ta7_sep)

tox[tox$ta7 == unique(tox$ta7)[j],]$treatment <- paste(tox[tox$ta7 == unique(tox$ta7)[j],]$ta7, tox[tox$ta7 == unique(tox$ta7)[j],]$chem)
HSD_AB_middle_ta7_sep <- HSD.test(y = tox[tox$ta7 == unique(tox$ta7)[j],]$HC5_adj,
                          trt =  tox[tox$ta7 == unique(tox$ta7)[j],]$treatment,
                          DFerror = anova_middle_ta7_sep$df.residual,
                          MSerror = deviance(anova_middle_ta7_sep)/anova_middle_ta7_sep$df.residual,
                          alpha = 0.05)
HSD_AB_middle_ta7_sep

HSD_AB_middle_ta7_sep <- cbind(HSD_AB_middle_ta7_sep$groups,HSD_AB_middle_ta7_sep$means$Max[order(HSD_AB_middle_ta7_sep$means$`tox[tox$ta7 == unique(tox$ta7)[j], ]$HC5_adj`, decreasing = TRUE)])
HSD_AB_middle_ta7_sep$treatment <- row.names(HSD_AB_middle_ta7_sep)
colnames(HSD_AB_middle_ta7_sep)[c(1,3)] <- c("param","max")
HSD_AB_middle_ta7 <- rbind(HSD_AB_middle_ta7, HSD_AB_middle_ta7_sep)
}
HSD_AB_middle_ta7


HSD_AB_middle_ta7$chem <- gsub(".* ", "", HSD_AB_middle_ta7$treatment)
HSD_AB_middle_ta7$ta7 <- gsub(" .*", "", HSD_AB_middle_ta7$treatment)

HSD_AB_middle_ta7$name <- NA
for(i in 1:length(unique(sum$chem))){
  HSD_AB_middle_ta7[HSD_AB_middle_ta7$chem == unique(sum$chem)[i],]$name <- unique(sum$name)[i]
}


summary <- ggplot(data = sum, aes(x = name, y = mean))+
  geom_bar(aes(x = name, y = mean, fill = name), stat = "identity", col = "#666666")+
  geom_errorbar(data = sum, ymin = sum$mean-1.96*sum$sem , ymax = sum$mean+1.96*sum$sem, lwd = 1.15, width = 0.5, col = "#666666")+
  scale_fill_manual(values = c("#98d191","#38bbf1","#f0ebc4","#d69b34"), name = "Chemical type")+
  #scale_y_continuous(breaks=seq(from = , to = ceiling(max(sum$mean+1.96*sum$sem)), by = 0.5), limits = c((min(sum$mean-1.96*sum$sem)), (max(sum$mean+1.96*sum$sem))))+
  geom_text(data = HSD_AB_middle_ta7, aes(x = name, y = round(param *1.4 +0.05,digits = 1), 
                                      label =  groups), size = 8)+
  labs(x = "Chemical", y = "Mean adjusted HC5 (\u00B5g/L, +- 2 s.d.)")+
  theme(
    axis.ticks.x=element_blank(),
    text = element_text(size=20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~ta7,  strip.position = "top", scales = "free")
summary

ggsave(here("figures","HC5_sample_summary_ta7_bar.png"), summary, width = 12, height = 8,  bg = "transparent")


#repeat but with chemical instead
tox <- tox[!is.na(tox$ta7),]
tox$treatment <- NA

A2G <- unique(tox$ta7)
A2G <-  A2G[order(A2G)]

for(i in 1:length(A2G)){
  tox[tox$ta7 == A2G[i],]$ta7 <- TA_desc[i]
}

sum$ta7 <- rep(TA_desc, nrow(chems))

sum <- sum %>%
  mutate(ta7 = fct_relevel(ta7, 
                           "Very large",
                           "Small, lowland, silicious",
                           "Large, lowland, silicious",
                           "Lowland calcareous",
                           "Large mid-altitude",
                           "Small mid-altitude",
                           "High altitude")) 
sum |>
  group_by(chem) |>
  summarise(
    range = max(mean)/min(mean)
  )

HSD_AB_middle_chem <- NULL
for(j in 1:length(unique(tox$chem))){
  print(unique(tox$chem)[j])
  anova_middle_chem_sep <- glm(tox[tox$chem == unique(tox$chem)[j],]$HC5_aj ~ tox[tox$chem == unique(tox$chem)[j],]$ta7)
  summary(anova_middle_chem_sep)
  
  tox[tox$chem == unique(tox$chem)[j],]$treatment <- paste(tox[tox$chem == unique(tox$chem)[j],]$ta7, tox[tox$chem == unique(tox$chem)[j],]$chem)
  HSD_AB_middle_chem_sep <- HSD.test(y = tox[tox$chem == unique(tox$chem)[j],]$HC5_aj,
                                      trt =  tox[tox$chem == unique(tox$chem)[j],]$treatment,
                                      DFerror = anova_middle_chem_sep$df.residual,
                                      MSerror = deviance(anova_middle_chem_sep)/anova_middle_chem_sep$df.residual,
                                      alpha = 0.05)
  HSD_AB_middle_chem_sep
  
  
  HSD_AB_middle_chem_sep <- cbind(HSD_AB_middle_chem_sep$groups,HSD_AB_middle_chem_sep$means$Max[order(HSD_AB_middle_chem_sep$means$`tox[tox$chem == unique(tox$chem)[j], ]$HC5_aj`, decreasing = TRUE)])
  HSD_AB_middle_chem_sep$treatment <- row.names(HSD_AB_middle_chem_sep)
  colnames(HSD_AB_middle_chem_sep)[c(1,3)] <- c("param","max")
  HSD_AB_middle_chem <- rbind(HSD_AB_middle_chem, HSD_AB_middle_chem_sep)
}
HSD_AB_middle_chem


HSD_AB_middle_chem$chem <- gsub(".* ", "", HSD_AB_middle_chem$treatment)
HSD_AB_middle_chem$ta7 <- NA
HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "ATZ",]$ta7 <- gsub(" ATZ", "", HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "ATZ",]$treatment)
HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "Cu",]$ta7 <- gsub(" Cu", "",   HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "Cu",]$treatment)
HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "IMD",]$ta7 <- gsub(" IMD", "", HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "IMD",]$treatment)
HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "LCH",]$ta7 <- gsub(" LCH", "", HSD_AB_middle_chem[HSD_AB_middle_chem$chem == "LCH",]$treatment)
HSD_AB_middle_chem$name <- NA

sum <- sum[!is.na(sum$ta7),]

for(i in 1:length(unique(sum$chem))){
  HSD_AB_middle_chem[HSD_AB_middle_chem$chem == unique(sum$chem)[i],]$name <- unique(sum$name)[i]
}

sum <- sum %>%
  mutate(ta7 = fct_relevel(ta7, 
                           "Very large",
                           "Small, lowland, silicious",
                           "Large, lowland, silicious",
                           "Lowland calcareous",
                           "Large mid-altitude",
                           "Small mid-altitude",
                           "High altitude")) 

HSD_AB_middle_chem <- HSD_AB_middle_chem%>%
  mutate(ta7 = fct_relevel(ta7, 
                           "Very large",
                           "Small, lowland, silicious",
                           "Large, lowland, silicious",
                           "Lowland calcareous",
                           "Large mid-altitude",
                           "Small mid-altitude",
                           "High altitude")) 


summary <- ggplot(data = sum, aes(x = ta7, y = mean))+
  geom_bar(aes(x = ta7, y = mean, fill = name), stat = "identity", col = "#666666")+
  geom_errorbar(data = sum, ymin = sum$mean-1.96*sum$sem , ymax = sum$mean+1.96*sum$sem, lwd = 1.15, width = 0.5, col = "#666666")+
  scale_fill_manual(values = c("#98d191","#38bbf1","#f0ebc4","#d69b34"), name = "Chemical type")+
  #scale_y_continuous(breaks=seq(from = , to = ceiling(max(sum$mean+1.96*sum$sem)), by = 0.5), limits = c((min(sum$mean-1.96*sum$sem)), (max(sum$mean+1.96*sum$sem))))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  geom_text(data = HSD_AB_middle_chem, aes(x = ta7, y =  round(param * 1.4 + 0.001, digits = 3),
                                      label = groups, pch = chem), size = 4, col = "#222222")+
  labs(y = "Mean HC5 (\u00B5g/L, \U00B1 2 s.e.)")+
  theme(
    axis.ticks.x=element_blank(),
    text = element_text(size=20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    axis.text.x = element_text(size=14),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~name,  strip.position = "top", scales = "free")
summary
ggsave(here("figures","HC5_sample_summary_chem_bar_handout.png"), summary, width = 12, height = 4.5,  bg = "transparent")
ggsave(here("figures","HC5_sample_summary_chem_bar.png"), summary, width = 14, height = 8,  bg = "transparent")

#####
#distributions of the HC5 values
hc5_tox_data <- read.csv(here("output","SSD_tox_data.csv"))

hc5_tox_data$chem <- unique(tox$chem)
hc5_tox_data$name <- unique(tox$name)

sum

hist_hc5 <- ggplot(data = tox, aes(x =log10(HC5_aj), fill = name))+
  #geom_boxplot(x ,outlier.shape = NA, fill = "#CCCCCC")+
  geom_histogram(bins = 64, col = "#666666")+
  geom_histogram(bins = 64, alpha = 0.9)+
  geom_segment(data = sum_all, aes(x = ((log10(upper) + log10(lower))/2) -0.15, y = c(300, 220, 350, 280), xend = log10(lower), yend = c(300, 220, 350, 280)), lwd = 1,
              arrow = arrow(length = unit(0.5, "cm")), col = "#444444")+
  geom_segment(data = sum_all, aes(x = ((log10(upper) + log10(lower))/2) +0.15, y = c(300, 220, 350, 280), xend = log10(upper), yend = c(300, 220, 350, 280)), lwd = 1,
             arrow = arrow(length = unit(0.5, "cm")), col = "#444444")+
  geom_vline(data = sum_all, aes(xintercept = log10(lower), fill = name), col = "#888888", lwd = 1.2, alpha = 0.6)+
  geom_vline(data = sum_all, aes(xintercept = log10(upper), fill = name), col = "#888888", lwd = 1.2, alpha = 0.6)+
  geom_text(data = sum_all, aes(x = ((log10(upper) + log10(lower))/2), y = c(300, 220, 350, 280), label = paste("x",round(CI_prop, digits = 1), sep = "")), size = 6, col = "#444444")+
  labs(x = "[Predicted HC5] (\u00B5g/L)", y = "Number of Assemblages")+
  #scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values = c("#98d191","#38bbf1","#f0ebc4","#d69b34"), name = "Chemical type")+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~name, scales = "free", nrow = 3)
hist_hc5
ggsave(plot = hist_hc5, here("figures","hc5_histogram_distibution.png"), width = 14, height = 8)


hist_hc5 <- ggplot(data = tox, aes(x = log10(HC5_aj), fill = name))+
  #geom_boxplot(x ,outlier.shape = NA, fill = "#CCCCCC")+
  geom_histogram(bins = 64, col = "#666666")+
  geom_histogram(bins = 64, alpha = 0.9)+
  geom_vline(data = hc5_tox_data, aes(xintercept = log10(HC5_aj), fill = name), col = "#AA0000", lwd = 1.2, alpha = 0.6)+
  #geom_text(data = sum_all, aes(x = log10(max)-0.2, y = c(150,120,150,150), label = paste("x",round(CI_prop, digits = 1), sep = "")), size = 8)+
  labs(x = "log10[Predicted HC5] (\u00B5g/L)", y = "Number of Assemblages")+
  #scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values = c("#98d191","#38bbf1","#f0ebc4","#d69b34"), name = "Chemical type")+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~name, scales = "free", nrow = 3)
hist_hc5
ggsave(plot = hist_hc5, here("figures","hc5_histogram_ssd.png"), width = 14, height = 8)


hist_hc5 <- ggplot(data = tox, aes(x = log10(HC5_aj), fill = name))+
  #geom_boxplot(x ,outlier.shape = NA, fill = "#CCCCCC")+
  geom_histogram(bins = 64, col = "#666666")+
  geom_histogram(bins = 64, alpha = 0.9)+
  geom_vline(data = hc5_tox_data, aes(xintercept = log10(HC5_aj), fill = name), col = "#AAAAAA", lwd = 1.2, alpha = 0)+
  #geom_text(data = sum_all, aes(x = log10(max)-0.2, y = c(150,120,150,150), label = paste("x",round(CI_prop, digits = 1), sep = "")), size = 8)+
  labs(x = "log10[Predicted HC5] (\u00B5g/L)", y = "Number of Assemblages")+
  #scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values = c("#98d191","#38bbf1","#f0ebc4","#d69b34"), name = "Chemical type")+
  theme(
    text = element_text(size = 20),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~name, scales = "free", nrow = 3)
hist_hc5
ggsave(plot = hist_hc5, here("figures","hc5_histogram_ssd_lineless.png"), width = 14, height = 8)


prop_unprotected <- c()
for(i in 1:length(unique(tox$chem))){
  prop_unprotected <- c(prop_unprotected,sum(tox[tox$chem == unique(tox$chem)[i],]$HC5_aj < hc5_tox_data[hc5_tox_data$chem == unique(tox$chem)[i],]$HC5_aj)/length(tox[tox$chem == unique(tox$chem)[i],]$HC5_aj))
}

prop_unprotected*100

#atrazine is the only one with sites more sensitive than the baseline
#is this proportion the same for all RTs?
tox_ATZ <- tox[tox$chem == "ATZ",]


prop_unprotected <- c()
for(i in order(unique(tox_ATZ$ta7))){
  print(unique(tox_ATZ$ta7)[i])
  prop_unprotected <- c(prop_unprotected,sum(tox_ATZ[tox_ATZ$ta7 == unique(tox_ATZ$ta7)[i],]$HC5_aj < hc5_tox_data[hc5_tox_data$chem == "ATZ",]$HC5_aj)/length(tox_ATZ[tox_ATZ$ta7 == unique(tox_ATZ$ta7)[i],]$HC5_aj))
}

prop_unprotected_ATZ <- data.frame(ta7 =  unique(tox_ATZ$ta7)[order(unique(tox_ATZ$ta7))], prop = round(prop_unprotected, digits = 2), Xs = rep(log10(hc5_tox_data[hc5_tox_data$chem == "ATZ",]$HC5_aj)+0.05, 6), Ys = c(1.75,30, 20, 25, 50, 30), name = rep("Atrazine", 6))


hist_hc5 <- ggplot(data = tox_ATZ, aes(x = log10(HC5_aj), fill = name))+
  #geom_boxplot(x ,outlier.shape = NA, fill = "#CCCCCC")+
  geom_histogram(bins = 64, col = "#666666")+
  geom_histogram(bins = 64, alpha = 0.9)+
  geom_vline(data = hc5_tox_data[hc5_tox_data$chem == "ATZ",], aes(xintercept = log10(HC5_aj)), col = "#AA0000", lwd = 1.2, alpha = 0.6)+
  #geom_text(data = prop_unprotected_ATZ, aes(x = Xs, y = Ys, label = prop), size = 8)+
  labs(x = "log10[Predicted HC5] (\u00B5g/L)", y = "Number of Assemblages")+
  #scale_y_continuous(trans = 'log10')+
  scale_fill_manual(values = c("#98d191","#38bbf1","#f0ebc4","#d69b34"), name = "Chemical type")+
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~ta7, scales = "free")
hist_hc5
ggsave(plot = hist_hc5, here("figures","hc5_histogram_ssd_Atrazine.png"), width = 12, height = 8)
