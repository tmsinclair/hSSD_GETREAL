#libraries
library(here)
library(tidyverse)
library(ggmap)
library(maps)
library(viridis)
library(sp)
library(sf)
library(agricolae)

pnt.in.poly <- function(pnts, polypnts){
  pip <- point.in.polygon(pnts[,1], pnts[,2], polypnts[,1], polypnts[,2])
  pip <- ifelse(pip > 0, 1, 0)
  return(data.frame(x = pnts[,1], y = pnts[,2], pip))
}
#####

###
#Map the site locations data
###

#read in the sample data
samples <- read.csv(here("output","hSSD_GETREAL_all_chems_properties.csv"))


chem_full <- c("Atrazine", "Copper", "Imidacloprid", "Lambda-cyhalothrin")
names(chem_full) <- unique(samples$chem)
#get a map of the world
world <- map_data('world')

#focus on Europe
europe <- c("Albania", "Austria", "Belarus", "Belgium", "Bosnia", "Bulgaria", "Croatia", "Cyprus", "Czech Republic", 
            "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Georgia",
            "Hungary","Iceland", "Ireland", "Italy", "Kosovo", "Latvia", "Lithuania", "Luxembourg", 
            "Malta", "Moldova", "Netherlands", "Montenegro", "North Macedonia", "Norway", "Poland", "Portugal", "Romania", "Russia", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland", "Turkey", "UK", "Ukraine")

world_europe <- map_data('world', region = europe)

#convert the co-ordinates into angular distances not measured to map european projections
convert <- SpatialPoints(data.frame(x = samples$longitude,y = samples$latitude), proj4string=CRS("+init=epsg:3035"))
lat_lon_points <- spTransform(convert, CRS("+init=epsg:4326"))

samples$lon_4326 <- lat_lon_points$x
samples$lat_4326 <- lat_lon_points$y

over(data.frame(samples$lon_4326, samples$lat_4326), world_europe)

pnts_sf <- st_as_sf(samples, coords = c('lat_4326', 'lon_4326'), crs = st_crs(world_europe))

pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, world_europe))
  , area = if_else(is.na(intersection), '', map$SUA_NAME16[intersection])
) 

pnts

st_intersects(pnts_sf$geometry, world_europe)

#####
plot_europe <- ggplot(world_europe, aes(x=long, y=lat, group = group)) +
  geom_polygon(fill="#EFEFCC", colour="#776655")+
  xlim(-12, 35) + ylim(35 ,72) +
  geom_point(data = samples,inherit.aes=FALSE, aes(x = lon_4326, y = lat_4326, col = RT_AF), size = 2.5, alpha = 0.75)+
  scale_colour_manual(values = c("#cc222e", "#f7d676", "#39c166", "#30b4f1", "#5041ac", "#db7ad1", "#AAAAAA"))+
  #scale_fill_viridis(name = "Distance (km)", option = "inferno", begin = 0.1, end = 0.9)+
  # scale_colour_gradientn(colours = c("#332d40", "#394373", "#327fa6", "#49d1ba"), name = "Typical assemblage", 
  #                        limits = c(149, 181), 
  #                        breaks = seq(150, 180, 10),
  #                        labels = seq(150, 180, 10),
  #                        guide = guide_colorbar(
  #                          direction = "horizontal",
  #                          barheight = unit(5, units = "mm"),
  #                          barwidth = unit(75, units = "mm"),
  #                          ticks.colour = "886633",
  #                          ticks.linewidth = 0,
  #                          title.position = 'top',
  #                          title.hjust = 0.5,
  #                          label.hjust = 0.5
  #                        ))+
  labs(x = "Longitude", y = "Latitude", title = "Typical assemblages of Europe")+
  theme(
    text = element_text(colour = "#FFEEDD"),
    plot.background = element_rect(fill = "#112233"),
    plot.title = element_text(colour = "#FFEEDD", hjust = 0.5, size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#112233"),
    legend.background = element_rect(fill = "#112233"),
    #   legend.text =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title.align = 1,
    legend.position = "bottom"
  )

plot_europe

ggsave(plot = plot_europe, here("figures", "typical_assemblages.png"), width = 8, height = 10, dpi = 300)

#
plot_europe <- ggplot(world_europe, aes(x=long, y=lat, group = group)) +
  geom_polygon(fill="#EFEFCC", colour="#776655")+
  xlim(-12, 35) + ylim(35 ,72) +
  geom_point(data = filter(samples, chem == "ATZ"),inherit.aes=FALSE, aes(x = lon_4326, y = lat_4326, col = log10(HC5_aj)), size = 2.5, alpha = 0.95)+
  scale_colour_viridis(name = "log10(HC5) (ug/L)", option = "inferno", begin = 0.1, end = 0.9, direction = -1)+
  # scale_colour_gradientn(colours = c("#332d40", "#394373", "#327fa6", "#49d1ba"), name = "Typical assemblage", 
  #                        limits = c(149, 181), 
  #                        breaks = seq(150, 180, 10),
  #                        labels = seq(150, 180, 10),
  #                        guide = guide_colorbar(
  #                          direction = "horizontal",
  #                          barheight = unit(5, units = "mm"),
  #                          barwidth = unit(75, units = "mm"),
  #                          ticks.colour = "886633",
  #                          ticks.linewidth = 0,
#                          title.position = 'top',
#                          title.hjust = 0.5,
#                          label.hjust = 0.5
#                        ))+
labs(x = "Longitude", y = "Latitude", title = "Assemblage sensitivity to atrazine")+
  theme(
    text = element_text(colour = "#FFEEDD"),
    plot.background = element_rect(fill = "#112233"),
    plot.title = element_text(colour = "#FFEEDD", hjust = 0.5, size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#112233"),
    legend.background = element_rect(fill = "#112233"),
    #   legend.text =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title.align = 1,
    legend.position = "bottom"
  )

plot_europe

ggsave(plot = plot_europe, here("figures", "spatial_atrazine.png"), width = 8, height = 10, dpi = 300)



#
plot_europe <- ggplot(world_europe, aes(x=long, y=lat, group = group)) +
  geom_polygon(fill="#EFEFCC", colour="#776655")+
  xlim(-12, 35) + ylim(35 ,72) +
  geom_point(data = filter(samples, chem == "Cu"),inherit.aes=FALSE, aes(x = lon_4326, y = lat_4326, col = log10(HC5_aj)), size = 2.5, alpha = 0.75)+
  scale_colour_viridis(name = "log10(HC5) (ug/L)", option = "inferno", begin = 0.1, end = 0.9, direction = -1)+
  # scale_colour_gradientn(colours = c("#332d40", "#394373", "#327fa6", "#49d1ba"), name = "Typical assemblage", 
  #                        limits = c(149, 181), 
  #                        breaks = seq(150, 180, 10),
  #                        labels = seq(150, 180, 10),
  #                        guide = guide_colorbar(
  #                          direction = "horizontal",
  #                          barheight = unit(5, units = "mm"),
  #                          barwidth = unit(75, units = "mm"),
  #                          ticks.colour = "886633",
  #                          ticks.linewidth = 0,
  #                          title.position = 'top',
#                          title.hjust = 0.5,
#                          label.hjust = 0.5
#                        ))+
labs(x = "Longitude", y = "Latitude", title = "Assemblage sensitivity to copper")+
  theme(
    text = element_text(colour = "#FFEEDD"),
    plot.background = element_rect(fill = "#112233"),
    plot.title = element_text(colour = "#FFEEDD", hjust = 0.5, size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#112233"),
    legend.background = element_rect(fill = "#112233"),
    #   legend.text =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title.align = 1,
    legend.position = "bottom"
  )

plot_europe

ggsave(plot = plot_europe, here("figures", "spatial_copper.png"), width = 8, height = 10, dpi = 300)

#
plot_europe <- ggplot(world_europe, aes(x=long, y=lat, group = group)) +
  geom_polygon(fill="#EFEFCC", colour="#776655")+
  xlim(-12, 35) + ylim(35 ,72) +
  geom_point(data = filter(samples, chem == "IMD"),inherit.aes=FALSE, aes(x = lon_4326, y = lat_4326, col = log10(HC5_aj)), size = 2.5, alpha = 0.75)+
  scale_colour_viridis(name = "log10(HC5) (ug/L)", option = "inferno", begin = 0.1, end = 0.9, direction = -1)+
  # scale_colour_gradientn(colours = c("#332d40", "#394373", "#327fa6", "#49d1ba"), name = "Typical assemblage", 
  #                        limits = c(149, 181), 
  #                        breaks = seq(150, 180, 10),
  #                        labels = seq(150, 180, 10),
  #                        guide = guide_colorbar(
  #                          direction = "horizontal",
  #                          barheight = unit(5, units = "mm"),
  #                          barwidth = unit(75, units = "mm"),
  #                          ticks.colour = "886633",
  #                          ticks.linewidth = 0,
  #                          title.position = 'top',
#                          title.hjust = 0.5,
#                          label.hjust = 0.5
#                        ))+
labs(x = "Longitude", y = "Latitude", title = "Assemblage sensitivity to imidacloprid")+
  theme(
    text = element_text(colour = "#FFEEDD"),
    plot.background = element_rect(fill = "#112233"),
    plot.title = element_text(colour = "#FFEEDD", hjust = 0.5, size = 22, face = "bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#112233"),
       legend.background = element_rect(fill = "#112233"),
    #   legend.text =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title.align = 1,
    legend.position = "bottom"
  )

plot_europe

ggsave(plot = plot_europe, here("figures", "spatial_imidacloprid.png"), width = 8, height = 10, dpi = 300)

plot_europe <- ggplot(world_europe, aes(x=long, y=lat, group = group)) +
  geom_polygon(fill="#EFEFCC", colour="#776655")+
  xlim(-12, 35) + ylim(35 ,72) +
  geom_point(data = filter(samples, chem == "LCH", HC5_aj<0.05),inherit.aes=FALSE, aes(x = lon_4326, y = lat_4326, col = log10(HC5_aj)), size = 2.5, alpha = 0.75)+
  scale_colour_viridis(name = "log10(HC5) (ug/L)", option = "inferno", begin = 0.1, end = 0.9, direction = -1)+
  # scale_colour_gradientn(colours = c("#332d40", "#394373", "#327fa6", "#49d1ba"), name = "Typical assemblage", 
  #                        limits = c(149, 181), 
  #                        breaks = seq(150, 180, 10),
  #                        labels = seq(150, 180, 10),
  #                        guide = guide_colorbar(
  #                          direction = "horizontal",
  #                          barheight = unit(5, units = "mm"),
  #                          barwidth = unit(75, units = "mm"),
  #                          ticks.colour = "886633",
  #                          ticks.linewidth = 0,
  #                          title.position = 'top',
#                          title.hjust = 0.5,
#                          label.hjust = 0.5
#                        ))+
labs(x = "Longitude", y = "Latitude", title = "Assemblage sensitivity to imidacloprid")+
  theme(
    text = element_text(colour = "#FFEEDD"),
    plot.background = element_rect(fill = "#112233"),
    plot.title = element_text(colour = "#FFEEDD", hjust = 0.5, size = 22, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#112233"),
    legend.background = element_rect(fill = "#112233"),
    #   legend.text =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title =  element_text(colour = "#886633", size = 14, face = "bold"),
    #   legend.title.align = 1,
    legend.position = "bottom"
  )

plot_europe

ggsave(plot = plot_europe, here("figures", "spatial_lambda-cyhalothrin.png"), width = 8, height = 10, dpi = 300)
#####
#look at this zonally

#Northern
N_lon_max <- max(filter(world_europe, region == "Finland" | region == "Sweden")$long)
N_lon_min <- min(filter(world_europe, region == "Finland" | region == "Sweden")$long)
N_lat_max <- max(filter(world_europe, region == "Finland" | region == "Sweden")$lat)
N_lat_min <- min(filter(world_europe, region == "Finland" | region == "Sweden")$lat)



N_samples <- samples[pnt.in.poly(select(samples, c(lon_4326, lat_4326)),filter(world_europe, region == "Norway" | region == "Denmark" | region == "Finland" | region == "Sweden"))$pip == 1,]
N_samples <- rbind(N_samples, samples[grepl("sweden", samples$sample),])
N_samples <- N_samples[!duplicated(N_samples),]

N_samples$zone <- "Northern"

#central
C_samples <- samples[pnt.in.poly(select(samples, c(lon_4326, lat_4326)),filter(world_europe, region == "Austria" | region == "Germany" | region == "Netherlands" | region == "UK"| region == "Czech Republic"| region == "Poland"| region == "Romania"))$pip == 1,]
C_samples$zone <- "Central"

#southern
S_samples <- samples[pnt.in.poly(select(samples, c(lon_4326, lat_4326)),filter(world_europe, region == "Potugal" | region == "Spain" | region == "France" | region == "Greece"| region == "Croatia"))$pip == 1,]

S_samples$zone <- "Southern"


samples_zone <- rbind(N_samples, C_samples, S_samples)

samples_zone_summary <- samples_zone |>
  group_by(zone, chem) |>
  summarise(
    mean = mean(HC5_aj),
    var = var(HC5_aj),
    max = max(HC5_aj),
    min = min(HC5_aj),
    range = max(HC5_aj) - min(HC5_aj),
    range_prop = max(HC5_aj) / min(HC5_aj),
    sd = sd(HC5_aj),
    length = length(HC5_aj),
    sem = sd/sqrt(length),
    quantile_2.5 = quantile(HC5_aj, probs = 0.025),
    quantile_50 = quantile(HC5_aj, probs = 0.5),
    quantile_97.5 = quantile(HC5_aj, probs = 0.975),
    quartile = (quantile(HC5_aj, probs = 0.75) - quantile(HC5_aj, probs = 0.25)) / mean(HC5_aj),
    range_95 = quantile_97.5 - quantile_2.5,
    range_95_prop = quantile_97.5 / quantile_2.5
  )

samples_zone_summary_all <- samples_zone |>
  group_by(chem) |>
  summarise(
    mean = mean(HC5_aj))

samples_zone_summary$mean_adj <- samples_zone_summary$mean

samples_zone_summary[samples_zone_summary$chem == "ATZ",]$mean_adj <- samples_zone_summary[samples_zone_summary$chem == "ATZ",]$mean_adj / samples_zone_summary_all[samples_zone_summary_all$chem == "ATZ",]$mean
samples_zone_summary[samples_zone_summary$chem == "Cu",]$mean_adj <- samples_zone_summary[samples_zone_summary$chem == "Cu",]$mean_adj / samples_zone_summary_all[samples_zone_summary_all$chem == "Cu",]$mean
samples_zone_summary[samples_zone_summary$chem == "IMD",]$mean_adj <- samples_zone_summary[samples_zone_summary$chem == "IMD",]$mean_adj / samples_zone_summary_all[samples_zone_summary_all$chem == "IMD",]$mean
samples_zone_summary[samples_zone_summary$chem == "LCH",]$mean_adj <- samples_zone_summary[samples_zone_summary$chem == "LCH",]$mean_adj / samples_zone_summary_all[samples_zone_summary_all$chem == "LCH",]$mean


dummy <- data.frame(zone = rep("Central", 4), 
                   mean = c(max(filter(samples_zone_summary, chem == "ATZ")$mean) + 1.96*max(filter(samples_zone_summary, chem == "ATZ")$sem),
                            max(filter(samples_zone_summary, chem == "Cu")$mean)  + 1.96*max(filter(samples_zone_summary, chem == "Cu")$sem),
                            max(filter(samples_zone_summary, chem == "IMD")$mean) + 1.96*max(filter(samples_zone_summary, chem == "IMD")$sem),
                            max(filter(samples_zone_summary, chem == "LCH")$mean) + 1.96*max(filter(samples_zone_summary, chem == "LCH")$sem)),
                   chem = unique(samples_zone_summary$chem)
                   )

zone_plot <- ggplot(data = samples_zone_summary, aes(x = zone, y = mean, fill = chem))+
  #geom_point(data=dummy, colour = "#FFFFFF")+
  geom_blank(data = dummy)+
  geom_bar(stat = "identity")+
  geom_errorbar(ymin = samples_zone_summary$mean-1.96*samples_zone_summary$sem , ymax = samples_zone_summary$mean+1.96*samples_zone_summary$sem, lwd = 1, width = 0.5, col = "#333333")+
  scale_fill_manual(values = c("#a0d699","#38bbf1","#f0eaba","#d69b34"), name = "Chemical type")+
  labs(x = "Chemical", y = "log[Predicted HC5] (\u00B5g/L  \U00B1 2 s.e.)")+
  theme(
    axis.ticks.x=element_blank(),
    text = element_text(size=16),
    legend.position = "bottom",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#FFFFFF"),
    strip.text = element_blank())+
  facet_wrap(~chem,  strip.position = "top", scales = "free")
zone_plot
ggsave(plot = zone_plot, here("figures", "zonal_comparison.png"), width = 12, height = 8, dpi = 300)

#####
#ANOVAs

anova_middle <- glm(samples_zone$HC5_aj ~ samples_zone$chem * samples_zone$zone)
summary(anova_middle)

samples_zone$treatment <- paste(samples_zone$zone, samples_zone$chem)

HSD_AB_middle <- HSD.test(y = samples_zone$HC5_aj,
                          trt =  samples_zone$treatment,
                          DFerror = anova_middle$df.residual,
                          MSerror = deviance(anova_middle)/anova_middle$df.residual,
                          alpha = 0.05)
HSD_AB_middle

HSD_AB_middle <- cbind(HSD_AB_middle$groups,HSD_AB_middle$means$Max[order(HSD_AB_middle$means$`samples_zone$HC5_aj`, decreasing = TRUE)])
HSD_AB_middle$treatment <- row.names(HSD_AB_middle)
colnames(HSD_AB_middle)[c(1,3)] <- c("param","max")
HSD_AB_middle

HSD_AB_middle$chem <- gsub(".* ", "", HSD_AB_middle$treatment)
HSD_AB_middle$zone <- gsub(" .*", "", HSD_AB_middle$treatment)


zone_plot <- ggplot(data = samples_zone_summary, aes(x = zone, y = mean, fill = chem))+
  #geom_point(data=dummy, colour = "#FFFFFF")+
  geom_blank(data = dummy)+
  geom_bar(stat = "identity", col = "#999999")+
  geom_errorbar(ymin = samples_zone_summary$mean-1.96*samples_zone_summary$sem , ymax = samples_zone_summary$mean+1.96*samples_zone_summary$sem, lwd = 1, width = 0.5, col = "#333333")+
  geom_text(data = HSD_AB_middle, aes(x = zone, y = c(270, 270, 200, 200, 200, 270,0.016, 0.02, 0.016, 0.016, 0.02, 0.02), label =  c("b","c","a","a","a","a", "a", "a", "a", "a", "a", "a")), size = 6)+
  scale_fill_manual(values = c("#a0d699","#38bbf1","#f0eaba","#d69b34"), name = "Chemical type")+
  labs(y = "Predicted HC5 (\u00B5g/L \U00B1 2 s.e.)")+
  theme(
    axis.ticks.x=element_blank(),
    text = element_text(size=17),
    legend.position = "none",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#EEEEEE"))+
  facet_wrap(~chem,  strip.position = "top", scales = "free", labeller = labeller(chem = chem_full), nrow = 1)#(ATZ = "Atrazine", Cu = "Copper", IMD = "Imidacloprid", LCH = "Lambda-cyhalothrin"))label_bquote(.(chem)-.("Atrazine").("Copper", "Imidacloprid", "Lambda-cyhalothrin")))
zone_plot
ggsave(plot = zone_plot, here("figures", "zonal_comparison.png"), width = 12, height = 4, dpi = 300)
ggsave(plot = zone_plot, here("figures", "zonal_comparison.png"), width = 14, height = 8, dpi = 300)

#
samples_zone_ATZ <- filter(samples_zone, chem == "IMD")

anova_middle <- glm(samples_zone_ATZ$HC5_aj ~ samples_zone_ATZ$zone)
summary(anova_middle)

samples_zone_ATZ$treatment <- paste(samples_zone_ATZ$zone)

HSD_AB_middle <- HSD.test(y = samples_zone_ATZ$HC5_aj,
                          trt =  samples_zone_ATZ$treatment,
                          DFerror = anova_middle$df.residual,
                          MSerror = deviance(anova_middle)/anova_middle$df.residual,
                          alpha = 0.05)
HSD_AB_middle

HSD_AB_middle <- cbind(HSD_AB_middle$groups,HSD_AB_middle$means$Max[order(HSD_AB_middle$means$`samples_zone_ATZ$HC5_aj`, decreasing = TRUE)])
HSD_AB_middle$treatment <- row.names(HSD_AB_middle)
colnames(HSD_AB_middle)[c(1,3)] <- c("param","max")
HSD_AB_middle

HSD_AB_middle$chem <- gsub(".* ", "", HSD_AB_middle$treatment)
HSD_AB_middle$zone <- gsub(" .*", "", HSD_AB_middle$treatment)


zone_plot <- ggplot(data = samples_zone_summary, aes(x = zone, y = mean, fill = chem))+
  #geom_point(data=dummy, colour = "#FFFFFF")+
  geom_blank(data = dummy)+
  geom_bar(stat = "identity")+
  geom_errorbar(ymin = samples_zone_summary$mean-1.96*samples_zone_summary$sem , ymax = samples_zone_summary$mean+1.96*samples_zone_summary$sem, lwd = 1, width = 0.5, col = "#333333")+
  #geom_text(data = HSD_AB_middle, aes(x = zone, y = c(rep(500, 3), rep(600, 3), rep(0.03, 3), rep(0.02, 3)), label =  groups), size = 4)+
  scale_fill_manual(values = c("#a0d699","#38bbf1","#f0eaba","#d69b34"), name = "Chemical type")+
  labs(x = "Chemical", y = "Predicted HC5 (\u00B5g/L)")+
  theme(
    axis.ticks.x=element_blank(),
    text = element_text(size=16),
    legend.position = "bottom",
    panel.background = element_rect(fill = "#FFFFFF"),
    axis.line = element_line(colour = "#333333"),
    strip.background = element_rect(fill = "#FFFFFF"),
    strip.text = element_blank())+
  facet_wrap(~chem,  strip.position = "top", scales = "free")
zone_plot
ggsave(plot = zone_plot, here("figures", "zonal_comparison.png"), width = 12, height = 8, dpi = 300)