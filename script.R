library(devtools)
library(dplyr)
library(metafor)
library(ggplot2)
library(cowplot)
make_pct <- function(x) (exp(x) - 1) * 100
CO2inc <- 616-372
data <- read.csv("AGB_effects2019.csv",na.strings=c("",NA))
data <- data[complete.cases(data$id),]
data <- data %>% dplyr::rename(Age=Age2,Nyears=nyears)
levels(data$Biome) <- list("Bo" = "Boreal_Forest","Cr"="Cropland", "Gr"="Grassland", "Sh"="Shrubland",
                           "Te"="Temperate_Forest","Tr"="Tropical_Forest")
trop <- filter(data, Biome=="Tropical_Forest") %>% mutate(obs = 1:nrow(.))
write.csv(trop,"tropical_review.csv")

overall <- rma.mv(yi, vi, data=trop, random = ~ 1 | Site.Name / obs) # overall effect
(make_pct(coef(summary(overall))[1])*100)/CO2inc
(make_pct(coef(summary(overall))[2])*100)/CO2inc

rma.mv(yi, vi, mods=~N, data=trop, random = ~ 1 | Site.Name / obs)
mod <- rma.mv(yi, vi, mods=~N-1, data=trop, random = ~ 1 | Site.Name / obs)
summary(mod)
make_pct(coef(summary(mod)))

Nm.n <- trop %>%  group_by(N) %>% summarise(n = n())
Nm.df <- coef(summary(mod)) %>% mutate(type="Nutrient fertilization", 
                                       factor=as.factor(c("Yes", "No")),
                                       size=Nm.n$n)
expNlow.mean <- (make_pct(Nm.df$estimate[2])*100)/CO2inc
expNlow.se <- (make_pct(Nm.df$se[2])*100)/CO2inc
expNhigh.mean <-  (make_pct(Nm.df$estimate[1])*100)/CO2inc
expNhigh.se <-  (make_pct(Nm.df$se[1])*100)/CO2inc

###### CMIP6 #######
library("readxl")
cmip <- read_excel("CMIP6.xlsx") %>%
  mutate(percentage=percentage*100, percentage100ppm = (percentage*100)/CO2inc) %>%
  group_by(N) %>%
  summarise(mean=mean(percentage100ppm), se=sd(percentage100ppm)/sqrt(n()))

###### Terrer et et al. 2019 #######
terrer.mean <- (12.5*100)/CO2inc
terrer.se <- (3*100)/CO2inc

###### Fleischer et al. 2019 #######
fleischer <- read.csv("biomass_per_change_original.csv") %>% 
  mutate(CBIO100=CBIO/2) %>% # Original effect is per 200ppm
  group_by(MODgr) %>% summarise(mean=mean(CBIO100), 
                                se=sd(CBIO100)/sqrt(n()))

###### ALL ######
all<- data.frame(id=c("CMIP6 - C", "CMIP6 - CN", "Experiments - NOFERT", "Experiments - FERT", 
                      "Experiments - Upscaled", "Process models - C", "Process models - CN", "Process models - CNP"), 
                 mean=c(cmip$mean[1], cmip$mean[2], expNlow.mean,expNhigh.mean, terrer.mean,
                        fleischer$mean[1], fleischer$mean[2], fleischer$mean[3]),
                 se=c(cmip$se[1], cmip$se[2], expNlow.se, expNhigh.se, terrer.se,
                      fleischer$se[1],fleischer$se[2],fleischer$se[3]),
                 fill_group=c("mechanistic","mechanistic","observational","observational","observational", "mechanistic", "mechanistic", "mechanistic"),
                 facet_group=c(rep("CMIP6 models",2), rep("Experiments",3),rep("Process models",3)),
                 color_group=as.factor(c("C","CN","meta-analysis","meta-analysis","CNP","C","CN","CNP")))
write.csv(all,"tropical_meta.csv")
figure <- ggplot(all, aes(x=reorder(id, -mean), y=mean, fill= fill_group, ymin=mean-se, ymax=mean+se)) + 
  geom_errorbar(size=.2,width=.3) + 
  #scale_y_continuous(expand = c(0, 0), limits=c(0,31)) +
  geom_col(width=.55, size=1, aes(col=color_group)) +
  xlab(NULL) + ylab (expression(paste("Biomass ", beta, " (% 100 ", ppm^-1,")", sep=" "))) +
  scale_fill_manual(values = c("#5B8EC1","#EBB851"), name="") +
  scale_color_manual(values=c("#ffab40", "#ffd100","#6cc24a", "black"), name="") +
  guides(color = guide_legend(override.aes = list(fill = "transparent"))) +
  facet_wrap(.~facet_group, scales="free_x") +
  theme_cowplot(font_size=12) + panel_border() +
  theme(strip.background = element_blank(),
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1)
  )
figure
save_plot("tropicalCO2.png", figure, dpi=1200, nrow=1, ncol=1, base_height = 6, base_width = 6,type = "cairo-png",bg = "white")

##################### Whittaker's Biomes #############################################
#devtools::install_github("valentinitnelav/plotbiomes")
library(dplyr)
library(ggplot2)
library(raster)
library(cowplot)
library(devtools)
library(metafor)
library(ggrepel)
library(plotbiomes)
library(rcartocolor)
make_pct <- function(x) (exp(x) - 1) * 100

effects <- read.csv("~/OneDrive - LLNL/FACEreview/data/AGB_effects2019.csv",na.strings=c("",NA)) %>%
  dplyr::rename(Age=Age2, Duration=nyears, Fertilized=N) %>% mutate(Fertilized=recode_factor(Fertilized,Nhigh="Yes",Nlow="No"), new_id = paste(Site.Name,Fertilized)) %>%
  dplyr::select(new_id,Age,Duration,Fertilized,MAT,MAP2, Biome) %>% distinct() %>% mutate(X = 1:n())

my_outliers <- get_outliers(tp = dplyr::select(effects, MAT,MAP2)%>%mutate(MAP2=MAP2*0.1)) # Outliers
effects$status <- ifelse(effects$X %in% my_outliers$row_idx, "out", "in")

biomes <- ggplot() + geom_polygon(data = Whittaker_biomes,aes(x = temp_c, y = precp_cm*10, fill = biome),
                                  colour = "gray98",size   = 1) +
  scale_fill_manual(name   = "",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors) +
  geom_point(data=effects,
             #data = dplyr::filter(effects, status=="in"), 
             aes(x = jitter(MAT,100), y = jitter(MAP2,100), size=Age, alpha=Duration, col=Fertilized),
             shape  = 21,
             fill   = "black",
             stroke = 1) +
  #geom_label_repel(data=effects,aes(x = MAT, y = MAP2, label=SITE), 
  #                               colour = "red", size = 2, nudge_x=-100, nudge_y=1000) +
  guides(fill = guide_legend(title="Ecosystem",order = 1),
         color = guide_legend(order = 4),
         alpha = guide_legend(order = 3),
         size = guide_legend(order = 2)) +
  scale_y_continuous(name = 'Precipitation (mm)', expand = c(0, 0)) +
  # - set range on OX axes and adjust the distance (gap) from OY axes
  scale_x_continuous(name = expression("Temperature " ( degree*C)),expand = c(0, 0), limits = c(-16,31)) +
  #coord_fixed(ratio = 1/10) + # aspect ratio, expressed as y / x
  theme_bw() + theme_classic(base_size = 8) + 
  theme(legend.justification = c(0, 1), # pick the upper left corner of the legend box and
        legend.position = c(0, 1), # adjust the position of the corner as relative to axis
        legend.background = element_rect(fill = NA), # transparent legend background
        #legend.box = "horizontal", # horizontal arrangement of multiple legends
        legend.spacing.x = unit(0.1, units = "cm"), # horizontal spacing between legends
        legend.spacing.y = unit(0.1, units = "cm"),
        panel.grid = element_blank(), # eliminate grids
        legend.margin=margin(t=0, r=0, b=0, l=0.1, unit="lines"),
        legend.key.width = unit(.3, "cm"), legend.key.height = unit(0.3, "cm"),
  )
biomes
save_plot("biomes.png", biomes, dpi=1200, nrow=1, ncol=1, base_height = 4, base_width = 4,type = "cairo-png")


