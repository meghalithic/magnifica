## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) creates trait measurements
# 2) investigates causes of measurement differences

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)
require(ggplot2)
require(reshape2)
require(lmodel2)
require(tidyverse)

#### SET UP ENV ----

col.form = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
             "#00A9FF", "#C77CFF", "#FF61CC")


col.traits = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
               "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")

#### LOAD DATA ----
df <- read.csv("Data/traits_25Sept2023.csv",
               header = TRUE)

form.meta <- read.csv("~/Documents/GitHub/bryozoa/stegino_metadata/newMetadata/formations.csv",
                      header = TRUE)

oxy.18 <- read.csv("Data/∂18O.csv",
                   header = TRUE)

locality.df <- read.csv("Data/All.NZ.Samples_EDM_31.07.2023_sheet1.csv",
                        header = TRUE)

##### MANIPULATE DATA ----

traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d",
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

##### REMOVE OUTLIER COLONIES -----
## bimodality in traits, driven by colonies with unusually small zooids
traits.melt.trim <- traits.melt[traits.melt$measurementType == "ln.zh" |
                                  traits.melt$measurementType == "ln.mpw.b" |
                                  traits.melt$measurementType == "ln.cw.m" |
                                  traits.melt$measurementType == "ln.cw.d" |
                                  traits.melt$measurementType == "ln.ow.m" |
                                  traits.melt$measurementType == "ln.oh" |
                                  traits.melt$measurementType == "ln.o.side" |
                                  traits.melt$measurementType == "ln.c.side",]

ggplot(traits.melt.trim) +
  geom_density(aes(x = measurementValue,
                   group = measurementType,
                   col = measurementType)) + #lots are bimodal
  ggtitle("Distribution of traits, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement")

## Find cutoff for frequencies of size bins
# frequency by size bin (quarter ln bins)
df.bins <- df %>% 
  mutate(zh.bin = cut(ln.zh, breaks = seq(5.5, 7.5, .1))) %>%
  as.data.frame()

df.bin.f <- df.bins %>%
  dplyr::group_by(zh.bin) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()
View(df.bin.f)
#(6.2,6.3]
#6.25 like I eyeballed

#write.csv(df.bin.f,
#          "./Results/zh.bin.frequency_8Sept2023.csv",
#          row.names = FALSE)

sm.traits <- df[df$ln.zh < 6.25,]
sm.colonies <- unique(sm.traits$colony.id) #41; was 95 images out of 891

bins <- c("(5.5,5.6]", "(5.6,5.7]", "(5.7,5.8]",
          "(5.8,5.9]", "(5.9,6]", "(6,6.1]", 
          "(6.1,6.2]", "(6.2,6.3]")

df.bins$zh.bin <- as.character(df.bins$zh.bin)

df.bins$sm <- FALSE
df.bins$sm[df.bins$zh.bin %in% bins] <- TRUE

#look at proportions
prop.sm <- df.bins %>% 
  dplyr::group_by(colony.id) %>%
  dplyr::summarise(n.zooid = length(zooid.id),
                   n.sm.zooid = sum(sm),
                   prop.sm = n.sm.zooid/n.zooid)
View(prop.sm)


#goes from 100% to 20%; perhaps make 20% the cut off
rm.col <- prop.sm$colony.id[prop.sm$prop.sm == 1]
length(rm.col) #31 colonies

reg.colonies <- df[!(df$colony.id %in% rm.col),]

##### WRITE OUT DATASETS ----
write.csv(reg.colonies,
          "./Results/colonies.traits_8Sept2023.csv",
          row.names = FALSE)
df <- reg.colonies

##### ABOUT TRAITS -----
###### RANGE -------
range(df$zh)
#266.3675 1701.3353
range(df$zh[df$formation == "NKLS"]) #450.4467 1701.3353; smallest are 26% of biggest
range(df$zh[df$formation == "NKBS"]) #266.3675 1495.0723; smallest are 17% size of biggest
range(df$zh[df$formation == "Tewkesbury"]) #330.5689 1668.3389; smallest are 19% size of biggest
range(df$zh[df$formation == "Waipuru"]) #326.9368 1361.0060; smallest are 24% size of biggest
range(df$zh[df$formation == "Upper Kai-Iwi"]) #331.6873 1362.1250; smallest are 24% size of biggest
range(df$zh[df$formation == "Tainui"]) #592.3886 1203.1921; smallest are 50% size of biggest
range(df$zh[df$formation == "SHCSBSB"]) #418.482 1333.873; smallest are 30% size of biggest

mean(df$zh)
median(df$zh)

mean(df$zh[df$formation == "NKLS"])
median(df$zh[df$formation == "NKLS"])

mean(df$zh[df$formation == "NKBS"])
median(df$zh[df$formation == "NKBS"])

mean(df$zh[df$formation == "Tewkesbury"])
median(df$zh[df$formation == "Tewkesbury"])

mean(df$zh[df$formation == "Waipuru"])
median(df$zh[df$formation == "Waipuru"])

mean(df$zh[df$formation == "Upper Kai-Iwi"])
median(df$zh[df$formation == "Upper Kai-Iwi"])

mean(df$zh[df$formation == "Tainui"])
median(df$zh[df$formation == "Tainui"])

mean(df$zh[df$formation == "SHCSBSB"])
median(df$zh[df$formation == "SHCSBSB"])

###### DISTRIBUTION ------
p.ln.zh <- ggplot(df) +
  geom_density(aes(x = ln.zh)) +
  ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

#ggsave(p.zh, 
#  file = "./Results/zooid.height_8Sept2023.png", width = 14, height = 10, units = "cm")

traits.melt <- melt(data = df,
                    id.vars = c("boxID", "zooid.id","imageName", "colony.id", "formation", "magnification"),
                    variable.name = "measurementType",
                    value.name = "measurementValue")
length(unique(traits.melt$colony.id)) #731 unique colonies [previously 742]

traits.stats <- traits.melt %>%
  group_by(measurementType) %>%
  summarise(avg = mean(measurementValue))

traits.stats.form <- traits.melt %>%
  group_by(measurementType, formation) %>%
  summarise(avg = mean(measurementValue))


p.dist <- ggplot(traits.melt) +
  geom_density(aes(x = log(measurementValue),
                   group = measurementType,
                   col = measurementType)) + #lots are bimodal
  ggtitle("Distribution of traits, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

ggsave(p.dist, 
       file = "./Results/trait.distribution_8Sept2023.png", width = 14, height = 10, units = "cm")

traits.melt.trim <- traits.melt[traits.melt$measurementType == "ln.zh" |
                                traits.melt$measurementType == "ln.mpw.b" |
                                traits.melt$measurementType == "ln.cw.m" |
                                traits.melt$measurementType == "ln.cw.d" |
                                traits.melt$measurementType == "ln.ow.m" |
                                traits.melt$measurementType == "ln.oh" |
                                traits.melt$measurementType == "ln.o.side" |
                                traits.melt$measurementType == "ln.c.side",]
  
ggplot(traits.melt.trim) +
  geom_density(aes(x = measurementValue,
                   group = measurementType,
                   col = measurementType)) + #lots are bimodal
  ggtitle("Distribution of traits, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement")

p.zh <- ggplot(df) +
  geom_density(aes(x = ln.zh,
                   group = formation,
                   col = formation)) + #lots are bimodal
  ggtitle(paste0("Distribution of traits, N zooids = ", length(unique(df$zooid.id)),
                 ", N colony = ", length(unique(df$colony.id)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.nkbs <- ggplot(df[df$formation == "NKBS",]) +
  geom_density(aes(x = ln.zh),
                   colour = col.form[1]) + #lots are bimodal
  ggtitle(paste0("NKBS, N zooids = ", length(unique(df$zooid.id[df$formation == "NKBS"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "NKBS"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.nkls <- ggplot(df[df$formation == "NKLS",]) +
  geom_density(aes(x = ln.zh), colour = col.form[2]) + #lots are bimodal
  ggtitle(paste0("NKLS, N zooids = ", length(unique(df$zooid.id[df$formation == "NKLS"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "NKLS"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.tewk <- ggplot(df[df$formation == "Tewkesbury",]) +
  geom_density(aes(x = ln.zh), colour = col.form[3]) + #lots are bimodal
  ggtitle(paste0("Tewkesbury, N zooids = ", length(unique(df$zooid.id[df$formation == "Tewkesbury"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "Tewkesbury"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.wai <- ggplot(df[df$formation == "Waipuru",]) +
  geom_density(aes(x = ln.zh), colour = col.form[4]) + #lots are bimodal
  ggtitle(paste0("Waipuru, N zooids = ", length(unique(df$zooid.id[df$formation == "Waipuru"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "Waipuru"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.uki <- ggplot(df[df$formation == "Upper Kai-Iwi",]) +
  geom_density(aes(x = ln.zh), colour = col.form[5]) + #lots are bimodal
  ggtitle(paste0("Upper Kai-Iwi, N zooids = ", length(unique(df$zooid.id[df$formation == "Upper Kai-Iwi"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "Upper Kai-Iwi"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.tai <- ggplot(df[df$formation == "Tainui",]) +
  geom_density(aes(x = ln.zh), colour = col.form[6]) + #lots are bimodal
  ggtitle(paste0("Tainui, N zooids = ", length(unique(df$zooid.id[df$formation == "Tainui"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "Tainui"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

p.zh.shcsbsb <- ggplot(df[df$formation == "SHCSBSB",]) +
  geom_density(aes(x = ln.zh), colour = col.form[7]) + #lots are bimodal
  ggtitle(paste0("SHCSBSB, N zooids = ", length(unique(df$zooid.id[df$formation == "SHCSBSB"])),
                 ", N colony = ", length(unique(df$colony.id[df$formation == "SHCSBSB"])))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

Fig = list(p.zh.nkbs, p.zh.nkls, 
           p.zh.tewk, p.zh.wai,
           p.zh.uki, p.zh.tai,
           p.zh.shcsbsb)

ml <- marrangeGrob(Fig, nrow = 7, ncol = 1)
ml





##### BIMODALITY -----
##explore bimodality, using zooid height as an example then see if it generalizes
p.ln.zh <- ggplot(df) +
  geom_density(aes(x = ln.zh)) +
  ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

#ggsave(p.zh, 
#  file = "./Results/zooid.height_8Sept2023.png", width = 14, height = 10, units = "cm")

ggplot(df) +
  geom_histogram(aes(x = ln.zh)) +
  ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Frequency") +
  scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

##driven by FORMATION?
#all but Punneki Limestone are bimodal
p.zh.form <- ggplot(df) +
  geom_density(aes(x = ln.zh,
                   group = formation,
                   col = formation)) + 
  ggtitle(paste0("Zooid height by formation, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

#ggsave(p.zh.form, 
# file = "./Results/zooid.height.by.formation_8Sept2023.png", width = 14, height = 10, units = "cm")

#only see bimodal in Waipurpu, NKBS, Upper Kai-Iwi

###### SAME INDIVIDUALS & COLONIES? ------
## really convince self that these are the same individuals
# make data more manageable by reducing it to the three formations
# make data even more manageable by reducing it to the small hump that is seen
nk.wa.uki <- df[df$formation == "NKBS" |
                df$formation == "Waipuru" |
                df$formation == "Upper Kai-Iwi",]
sm.zoo <- nk.wa.uki[nk.wa.uki$ln.zh <= 6.25,]
nrow(sm.zoo) #492 (was 943) zooids
length(unique(sm.zoo$colony.id)) #34 (was 64) colonies
table(sm.zoo$colony.id) #a lot of one offs, but some clusters
## look at a couple of these:
sm.zoo[sm.zoo$colony.id == "077CV",]

## are these individuals from the same colony or across colonies?

p.ow.zh <- ggplot(data = df) +
  geom_smooth(aes(x = ln.zh, y = ln.ow.m), method = "lm") +
  geom_point(aes(x = ln.zh, y = ln.ow.m)) + #two clusters
  ggtitle(paste0("Scaling of operculum mid-width with zooid height, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "LN Operculum mid-width") +
  scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
summary(lm(df$ln.ow.m ~ df$ln.zh)) #slope = 0.772476

#ggsave(p.ow.zh, 
#file = "./Results/ow.zh.scaling.png", width = 14, height = 10, units = "cm")

#look at 4 images:
#2 with the largest operculum width
#2 with the largest zooid height

slice_max(df, n = 2, order_by = ln.ow.m)

slice_max(df, n = 2, order_by = ln.zh)

zh.owm.model <- lm(ln.ow.m ~ ln.zh,
                   data = df)

###### LOOK AT RESULTS FROM IMG J -----
table(sm.zoo$colony.id)
unique(sm.zoo$imageName[sm.zoo$colony.id == "077CV"]) #from images 1, 2, 3
## image 077_CV
imgj.res_077_1 <- read.csv("~/Desktop/007_CV_1_Results_13Sept2023.csv.csv", header = TRUE)
imgj.res_077_2 <- read.csv("~/Desktop/007_CV_2_Results_13Sept2023.csv", header = TRUE)
imgj.res_077_3 <- read.csv("~/Desktop/077_CV_3_Results_13Sept2023.csv", header = TRUE)
## scale
imgj.res_077_1[1,] #first row is the scale to .1 mm; scale is 28 pixels per .1 mm; or 280 per 1mm; or 280000 per 1 micrometer
# 28 px/.1 mm or .1 mm / 28 px
imgj.res_077_1$scale.zh <- (imgj.res_077_1$Length/28)*1000
imgj.res_077_2[1,] #29.017
imgj.res_077_2$scale.zh <- (imgj.res_077_2$Length/29.017)*1000
imgj.res_077_3[1,] #0.109
imgj.res_077_3$scale.zh <- (imgj.res_077_3$Length/0.109)*1000

imgj_077 <- rbind(imgj.res_077_1, imgj.res_077_2, imgj.res_077_3)

## log
imgj_077$ln.zh <- log(imgj_077$scale.zh)
range(imgj_077$ln.zh) # all small!!
#-1.715532  6.075998
range(df$ln.zh[df$imageName == "077_CV_1_15v_x30_BSE" |
                 df$imageName == "077_CV_2_15v_x30_BSE" |
                 df$imageName == "077_CV_3_15v_x30_BSE"])
#5.789308 6.086177

imgj.res_696_1 <- read.csv("~/Desktop/696_CC_1_Results_13Sept2023.csv", header = TRUE)
imgj.res_696_2 <- read.csv("~/Desktop/696_CC_2_Results_13Sept2023.csv", header = TRUE)
imgj.res_696_1[1,] #60
imgj.res_696_1$scale.zh <- (imgj.res_696_1$Length/60)*1000
imgj.res_696_2[1,] #60
imgj.res_696_2$scale.zh <- (imgj.res_696_2$Length/60)*1000
imgj.res_696 <- rbind(imgj.res_696_1, imgj.res_696_2)
imgj.res_696$ln.zh <- log(imgj.res_696$scale.zh)
range(imgj.res_696$ln.zh)
#6.907755 9.148988
range(df$ln.zh[df$imageName == "696_CC_1_10v_x30_BSE" |
                 df$imageName == "696_CC_2_10v_x30_BSE"])
#6.599759 7.061558

##273S
unique(sm.zoo$imageName[sm.zoo$colony.id == "273S"]) #from images 1, 2, 3
imgj.res_273_1 <- read.csv("~/Desktop/273_S_1_Results_13Sept2023.csv", header = TRUE)
imgj.res_273_2 <- read.csv("~/Desktop/273_S_2_Results_13Sept2023.csv", header = TRUE)
imgj.res_273_3 <- read.csv("~/Desktop/273_S_3_Results_13Sept2023.csv", header = TRUE)
imgj_273 <- rbind(imgj.res_273_1, imgj.res_273_2, imgj.res_273_3)
## scale
imgj_273$scale.zh <- imgj_273$Length/.606
## log
imgj_273$ln.zh <- log(imgj_273$scale.zh)
imgj_273$ln.zh #all small!

##### SMALL HUMP -----
# know that these are all related, use zh as test
## test where hump is...
# frequency by size bin (quarter ln bins)
df.bins <- df %>% 
  mutate(zh.bin = cut(ln.zh, breaks = seq(5.5, 7.5, .1))) %>%
  as.data.frame()

df.bin.f <- df.bins %>%
  dplyr::group_by(zh.bin) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()
View(df.bin.f)
#(6.2,6.3]
#6.25 like I eyeballed

#write.csv(df.bin.f,
#          "./Results/zh.bin.frequency_8Sept2023.csv",
#          row.names = FALSE)

sm.traits <- df[df$ln.zh < 6.25,]
sm.colonies <- unique(sm.traits$colony.id) #41; was 95 images out of 891

bins <- c("(5.5,5.6]", "(5.6,5.7]", "(5.7,5.8]",
          "(5.8,5.9]", "(5.9,6]", "(6,6.1]", 
          "(6.1,6.2]", "(6.2,6.3]")

df.bins$zh.bin <- as.character(df.bins$zh.bin)

df.bins$sm <- FALSE
df.bins$sm[df.bins$zh.bin %in% bins] <- TRUE

#look at proportions
prop.sm <- df.bins %>% 
  dplyr::group_by(colony.id) %>%
  dplyr::summarise(n.zooid = length(zooid.id),
            n.sm.zooid = sum(sm),
            prop.sm = n.sm.zooid/n.zooid)
View(prop.sm)

#write.csv(prop.sm,
#          "./Results/proportion.small.colonies_8Sept2023.csv",
#          row.names = FALSE)

#goes from 100% to 20%; perhaps make 20% the cut off
rm.col <- prop.sm$colony.id[prop.sm$prop.sm == 1]
length(rm.col) #31 colonies

small.colonies <- df[df$colony.id %in% rm.col,]
reg.colonies <- df[!(df$colony.id %in% rm.col),]

#for those with 20%, see where the sizes are
sm.zooids <- prop.sm$colony.id[prop.sm$prop.sm < 1 &
                                 prop.sm$prop.sm > 0]
length(sm.zooids) #19

range(df$zh[df$colony.id %in% sm.zooids])
#330.5689 1701.3353
sort(df$ln.zh[df$colony.id %in% sm.zooids])
#only 11 below 6.25, so probably fine

##### WRITE OUT DATASETS ----
write.csv(reg.colonies,
          "./Results/colonies.traits_8Sept2023.csv",
          row.names = FALSE)

######## SUMMARY STATS ------
reg.mean_by_formation_colony = reg.colonies %>% #use this going forward
  dplyr::group_by(formation, colony.id) %>%
  dplyr::summarize(n.zooid = length(unique(zooid.id)),
                   
                   min.zh = min(ln.zh, na.rm = T),
                   max.zh = max(ln.zh, na.rm = T),
                   avg.zh = mean(ln.zh, na.rm = T),
                   sd.zh = sd(ln.zh, na.rm = T),
                   
                   min.mpw.b = min(ln.mpw.b, na.rm = T),
                   max.mpw.b = max(ln.mpw.b, na.rm = T),
                   avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                   sd.mpw.b = sd(ln.mpw.b, na.rm = T),
                   
                   min.cw.m = min(ln.cw.m, na.rm = T),
                   max.cw.m = max(ln.cw.m, na.rm = T),
                   avg.cw.m = mean(ln.cw.m, na.rm = T),
                   sd.cw.m = sd(ln.cw.m, na.rm = T),
                   
                   min.cw.d = min(ln.cw.d, na.rm = T),
                   max.cw.d = max(ln.cw.d, na.rm = T),
                   avg.cw.d = mean(ln.cw.d, na.rm = T),
                   sd.cw.d = sd(ln.cw.d, na.rm = T),
                   
                   min.ow.m = min(ln.ow.m, na.rm = T),
                   max.ow.m = max(ln.ow.m, na.rm = T),
                   avg.ow.m = mean(ln.ow.m, na.rm = T),
                   sd.ow.m = sd(ln.ow.m, na.rm = T),
                   
                   min.oh = min(ln.oh, na.rm = T),
                   max.oh = max(ln.oh, na.rm = T),
                   avg.oh = mean(ln.oh, na.rm = T),
                   sd.oh = sd(ln.oh, na.rm = T),
                   
                   min.o.side = min(ln.o.side, na.rm = T),
                   max.o.side = max(ln.o.side, na.rm = T),
                   avg.o.side = mean(ln.o.side, na.rm = T),
                   sd.o.side = sd(ln.o.side, na.rm = T),
                   
                   min.c.side = min(ln.c.side, na.rm = T),
                   max.c.side = max(ln.c.side, na.rm = T),
                   avg.c.side = mean(ln.c.side, na.rm = T),
                   sd.c.side = sd(ln.c.side, na.rm = T),
                   
                   min.area = min(ln.area, na.rm = T),
                   max.area = max(ln.area, na.rm = T),
                   avg.area = mean(ln.area, na.rm = T),
                   sd.area = sd(ln.area, na.rm = T)) %>%
  as.data.frame()
nrow(reg.mean_by_formation_colony) #541

reg.mean_by_formation = reg.colonies %>%
  dplyr::group_by(formation) %>%
  dplyr::summarize(num.col = length(unique(colony.id)),
                   num.zooid = length(unique(zooid.id)),
                   avg.zooid = ceiling(num.zooid/num.col), #round up to nearest integer
                   
                   min.zh = min(ln.zh, na.rm = T),
                   max.zh = max(ln.zh, na.rm = T),
                   avg.zh = mean(ln.zh, na.rm = T),
                   sd.zh = sd(ln.zh, na.rm = T),
                   
                   min.mpw.b = min(ln.mpw.b, na.rm = T),
                   max.mpw.b = max(ln.mpw.b, na.rm = T),
                   avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                   sd.mpw.b = sd(ln.mpw.b, na.rm = T),
                   
                   min.cw.m = min(ln.cw.m, na.rm = T),
                   max.cw.m = max(ln.cw.m, na.rm = T),
                   avg.cw.m = mean(ln.cw.m, na.rm = T),
                   sd.cw.m = sd(ln.cw.m, na.rm = T),
                   
                   min.cw.d = min(ln.cw.d, na.rm = T),
                   max.cw.d = max(ln.cw.d, na.rm = T),
                   avg.cw.d = mean(ln.cw.d, na.rm = T),
                   sd.cw.d = sd(ln.cw.d, na.rm = T),
                   
                   min.ow.m = min(ln.ow.m, na.rm = T),
                   max.ow.m = max(ln.ow.m, na.rm = T),
                   avg.ow.m = mean(ln.ow.m, na.rm = T),
                   sd.ow.m = sd(ln.ow.m, na.rm = T),
                   
                   min.oh = min(ln.oh, na.rm = T),
                   max.oh = max(ln.oh, na.rm = T),
                   avg.oh = mean(ln.oh, na.rm = T),
                   sd.oh = sd(ln.oh, na.rm = T),
                   
                   min.o.side = min(ln.o.side, na.rm = T),
                   max.o.side = max(ln.o.side, na.rm = T),
                   avg.o.side = mean(ln.o.side, na.rm = T),
                   sd.o.side = sd(ln.o.side, na.rm = T),
                   
                   min.c.side = min(ln.c.side, na.rm = T),
                   max.c.side = max(ln.c.side, na.rm = T),
                   avg.c.side = mean(ln.c.side, na.rm = T),
                   sd.c.side = sd(ln.c.side, na.rm = T),
                   
                   min.area = min(ln.area, na.rm = T),
                   max.area = max(ln.area, na.rm = T),
                   avg.area = mean(ln.area, na.rm = T),
                   sd.area = sd(ln.area, na.rm = T)) %>%
  as.data.frame()

#### CORRELATIONS & ALLOMETRIES ----
## are these coming from the same individuals??
## ask KLV for other metadata for sites
## OH hump is on the other side

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[2]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[2]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[2]) +
  scale_color_manual(values = col.form)
summary(lm(df[, traits[2]] ~ df[, traits[1]]))
#slope = 0.79460; p-value < 2.2e-16; r2 = 0.4722

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[3]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[3]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[3]) +
  scale_color_manual(values = col.form)
summary(lm(df[, traits[3]] ~ df[, traits[1]]))
#slope = 0.72434; p-value < 2.2e-16; r2 = 0.3684

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[4]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[4]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[4]) +
  scale_color_manual(values = col.form)
summary(lm(df[, traits[4]] ~ df[, traits[1]]))
#slope = 0.758859; p-value < 2.2e-16; r2 = 0.503

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[5]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[5]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[5]) +
  scale_color_manual(values = col.form)
summary(lm(df[, traits[5]] ~ df[, traits[1]]))
#slope = 0.772476; p-value < 2.2e-16; r2 = 0.642

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[6]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[6]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[6]) +
  scale_color_manual(values = col.form)
summary(lm(df[, traits[6]] ~ df[, traits[1]]))
#slope = 0.345976; p-value < 2.2e-16; r2 = 0.325

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[7]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[7]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[7]) +
  scale_color_manual(values = col.form)
#ln.c.side is not an issue really
summary(lm(df[, traits[7]] ~ df[, traits[1]]))
#slope = 1.092931; p-value < 2.2e-16; r2 = 0.8933

ggplot(data = df) +
  geom_smooth(aes(x = df[, traits[1]],
                  y = df[, traits[8]],
                  alpha = 0.5)) +
  geom_point(aes(x = df[, traits[1]],
                 y = df[, traits[8]],
                 group = formation,
                 col = formation,
                 alpha = 0.5)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_y_continuous(name = traits[8]) +
  scale_color_manual(values = col.form)
#ln.o.side is not an issue really
summary(lm(df[, traits[7]] ~ df[, traits[1]]))
#slope = 1.092931; p-value < 2.2e-16; r2 = 0.8933

#### LOOK AT CHANGES OVER TIME AND BETWEEN FORMATIONS ------
## Add meta data
form.df <- form.meta[1:7,] #in same order as mean_by_formation

for(i in 1:nrow(form.df)){
  form.df$mean.age[i] <- mean(form.df$Start_age[i], form.df$End_age[i], na.rm = TRUE)
}

form.df$age.range <- ""
for(i in 1:nrow(form.df)){
  form.df$age.range[i] <- form.df$Start_age[i] - form.df$End_age[i]
}
form.df$age.range <- as.numeric(form.df$age.range)

reg.mean_by_formation.meta <- merge(reg.mean_by_formation, form.df,
                                by.x = "formation",
                                by.y = "formationCode")

## overall differnce in zooid height
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "NKLS"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "SHCSBSB"]
#-0.1514114 (decrease in length)
reg.mean_by_formation$avg.ow.m[reg.mean_by_formation$formation == "NKLS"] - reg.mean_by_formation$avg.ow.m[reg.mean_by_formation$formation == "SHCSBSB"]
#-0.1479244 (decrease in width)

##between formations
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "NKLS"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "NKBS"]
#0.1162435 (1.123269 diff)
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "NKBS"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Tewkesbury"]
#-0.08700402 (0.9166734 diff)
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Tewkesbury"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Waipuru"]
#0.06027717 (1.062131 diff)
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Waipuru"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Upper Kai-Iwi"]
#0.01439653 (1.014501 diff)
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Upper Kai-Iwi"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Tainui"]
#-0.2706406 (0.7628906 diff)
reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "Tainui"] - reg.mean_by_formation$avg.zh[reg.mean_by_formation$formation == "SHCSBSB"]
#0.01531607 (1.015434 diff)

## how is sd a function of sample size (number of zooids and number of colonies)?
#plot sd per colony by zooid no
ggplot(reg.mean_by_formation_colony) + 
  geom_point(aes(x = n.zooid, y = sd.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Number of Zooids per Colony") +
  scale_y_continuous(expression(sd~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

#plot sd per formation by colony no
ggplot(reg.mean_by_formation) + 
  geom_point(aes(x = num.col, y = sd.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Number of Colonies per Colony") +
  scale_y_continuous(expression(sd~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

anova(lm(reg.mean_by_formation$sd.zh ~ reg.mean_by_formation$num.col + reg.mean_by_formation$num.zooid))

#use formation means
#mean_by_formation

ggplot(data = reg.mean_by_formation.meta) +
  geom_point(aes(x = age.range, y = avg.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Age Range (Ma)") +
  scale_y_continuous(name = "Average Zooid Height (um)") +
  scale_color_manual(values = col.form)

ggplot(data = reg.mean_by_formation.meta) +
  geom_point(aes(x = mean.age, y = avg.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Age (Ma)") +
  scale_y_continuous(expression(Average~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)


df.form.meta <- merge(df, form.df,
                      by.x = "formation",
                      by.y = "formationCode")
                      
ggplot(data = df.form.meta) +
  geom_point(aes(x = mean.age, y = ln.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Age (Ma)") +
  scale_y_continuous(expression(ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

reg.form.meta <- merge(reg.colonies, form.df,
                      by.x = "formation",
                      by.y = "formationCode")

ggplot(data = reg.form.meta) +
  geom_point(aes(x = ln.cw.d, y = ln.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Cryptocyst Width at the Distal End (um)") +
  scale_y_continuous(expression(ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

box.ln.zh <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.zh",], 
       aes(x = formation, 
           y = measurementValue, 
           fill = formation)) +
  geom_boxplot() +
  scale_color_manual(values = col.form) +
  scale_fill_manual(values = col.form) +
  ggtitle("Boxplots of LN Zooid Heights") +
  scale_x_discrete(name = "Formation",
                   guide = guide_axis(angle = 45)) +
  ylab(expression(ln~Zooid~Height~(mu*m))) + 
  theme(text = element_text(size = 16),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggsave(box.ln.zh, 
       file = "./Results/boxplot.ln.zh.png", 
       width = 14, height = 10, units = "cm")