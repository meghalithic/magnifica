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
images.meta <- read.csv("./Data/meta.images.8Sept2023.csv",
                   header = TRUE,
                   sep = ",")

df.filter <- read.table("./Data/filteredImages.csv",
                        header = TRUE,
                        sep = ";")

form.meta <- read.csv("~/Documents/GitHub/bryozoa/stegino_metadata/newMetadata/formations.csv",
                      header = TRUE)

oxy.18 <- read.csv("Data/∂18O.csv",
                   header = TRUE)

locality.df <- read.csv("Data/All.NZ.Samples_EDM_31.07.2023_sheet1.csv",
                        header = TRUE)

#### EXPLORE DATA ----
nrow(df.filter) #1834
nrow(images.meta) #6443

#### REDUCE TO SELECTED IMAGES ####

df.filter$fileName.old <- c()
for(i in 1:nrow(df.filter)){
  df.filter$fileName.old[i] <- str_split(df.filter$path.tif[i], "/")[[1]][length(str_split(df.filter$path.tif[i], "/")[[1]])]
}

images.filter <- images.meta[images.meta$fileName.tif %in% df.filter$fileName.old,]
nrow(images.filter) #6438
length(unique(images.filter$fileName.tif)) #1394

images.filter$formation <- factor(images.filter$formation, 
                                  levels = c("NKLS", "NKBS", "Tewkesbury",
                                             "Waipuru", "Upper Kai-Iwi",
                                             "Tainui", "SHCSBSB")) 

#### CALCULATE DISTANCES ----
#measurements based off Voje et al. 2020 https://doi.org/10.5061/dryad.t4b8gthxm
#zh is similar to LZ
#cw.m is similar to WZ
#oh is similar to LO
#ow.m is similar to WO
#see "stegs_linear_24Mar2023.png" for linear measurements
#z = zooid
#h = height
#o = operculum
#l = left
#r = right
#w = width
#m = mid
#b = base
#mp = median process
#c = cryptocyst
#d = distal

#X = |(X1-X0)|
#Y = |(Y1-Y0)|
#D^2 = X^2 + Y^2

## Zooid height (maximum height at centerline): 4 to 12
zh.x <- abs((images.filter$X4-images.filter$X12))
zh.y <-  abs((images.filter$Y4-images.filter$Y12))
zh <- sqrt(((zh.x)^2 + (zh.y)^2))

## Operculum height of left side: 4 to 20
oh.l.x <- abs((images.filter$X4-images.filter$X20))
oh.l.y <-  abs((images.filter$Y4-images.filter$Y20))
oh.l <- sqrt(((oh.l.x)^2 + (oh.l.y)^2))

## Operculum height of right side: 4 to 21
oh.r.x <- abs((images.filter$X4-images.filter$X20))
oh.r.y <-  abs((images.filter$Y4-images.filter$Y20))
oh.r <- sqrt(((oh.r.x)^2 + (oh.r.y)^2))

## Operculum mid-width (maximum width at centerline): 19 to 0
ow.m.x <- abs((images.filter$X19-images.filter$X0))
ow.m.y <-  abs((images.filter$Y19-images.filter$Y0))
ow.m <- sqrt(((ow.m.x)^2 + (ow.m.y)^2))

## Operculum base width: 21 to 20
ow.b.x <- abs((images.filter$X21-images.filter$X20))
ow.b.y <-  abs((images.filter$Y21-images.filter$Y20))
ow.b <- sqrt(((ow.b.x)^2 + (ow.b.y)^2))

## Operculum side length of right side: 21 to 18
o.side.r.x <- abs((images.filter$X21-images.filter$X18))
o.side.r.y <- abs((images.filter$Y21-images.filter$Y18))
o.side.r <- sqrt(((o.side.r.x)^2 + (o.side.r.y)^2))

## Operculum side length of left side: 20 to 17
o.side.l.x <- abs((images.filter$X20-images.filter$X17))
o.side.l.y <- abs((images.filter$Y20-images.filter$Y17))
o.side.l <- sqrt(((o.side.l.x)^2 + (o.side.l.y)^2))

## Operculum average side length:
o.side.avg <- rowMeans(cbind(o.side.r, o.side.l))

p.o.side.rl <- ggplot() +
  geom_point(aes(x = o.side.l, y = o.side.r)) +
  ggtitle("Operculum Side Length, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum Length Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Length Left Side (pixels)")

ggsave(p.o.side.rl, file = "./Results/operculum_length.png", width = 14, height = 10, units = "cm")

summary(lm(o.side.r ~ o.side.l)) 
# slope = 0.912352; p-value < 2.2e-16; r2 = 0.9494

## Operculum height
oh <- (.5/ow.b)*sqrt(ow.b+oh.r+oh.l)

## Median process base width: 5 to 6
mpw.b.x <- abs((images.filter$X5-images.filter$X6))
mpw.b.y <-  abs((images.filter$Y5-images.filter$Y6))
mpw.b <- sqrt(((mpw.b.x)^2 + (mpw.b.y)^2))

## Cryptocyst mid-width: 10 to 11
cw.m.x <- abs((images.filter$X10-images.filter$X11))
cw.m.y <-  abs((images.filter$Y10-images.filter$Y11))
cw.m <- sqrt(((cw.m.x)^2 + (cw.m.y)^2))

## Cryptocyst base width: 9 to 1
cw.b.x <- abs((images.filter$X9-images.filter$X1))
cw.b.y <-  abs((images.filter$Y9-images.filter$Y1))
cw.b <- sqrt(((cw.b.x)^2 + (cw.b.y)^2))

## Cryptocyst distal width: 8 to 7
cw.d.x <- abs((images.filter$X8-images.filter$X7))
cw.d.y <-  abs((images.filter$Y8-images.filter$Y7))
cw.d <- sqrt(((cw.d.x)^2 + (cw.d.y)^2))

## Cryptocyst side length of right side: 1 to 7
c.side.r.x <- abs((images.filter$X1-images.filter$X7))
c.side.r.y <- abs((images.filter$Y1-images.filter$Y7))
c.side.r <- sqrt(((c.side.r.x)^2 + (c.side.r.y)^2))

## Cryptocyst side length of left side: 9 to 8
c.side.l.x <- abs((images.filter$X9-images.filter$X8))
c.side.l.y <- abs((images.filter$Y9-images.filter$Y8))
c.side.l <- sqrt(((c.side.l.x)^2 + (c.side.l.y)^2))

## Cryptocyst average side length:
c.side.avg <- rowMeans(cbind(c.side.r, c.side.l))

p.c.side.rl <- ggplot() +
  geom_point(aes(x = c.side.l, y = c.side.l)) +
  ggtitle("Cryptocyst Side Length, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Cryptocyst Length Right Side (pixels)") +
  scale_x_continuous(name = "Cryptocyst Length Left Side (pixels)")

ggsave(p.c.side.rl, file = "./Results/cryptocyst_length.png", width = 14, height = 10, units = "cm")

summary(lm(c.side.r ~ c.side.l)) 
# slope = 0.867776; p-value: < 2.2e-16; r2 = 0.7858

## right v left side of operculum
p.oh.rl <- ggplot() +
  geom_point(aes(x = oh.l, y = oh.r)) +
  ggtitle("Operculum Height, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum Height Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Height Left Side (pixels)")

oh.model <- lmodel2(formula = oh.r ~ oh.l,
                    range.x = "relative", 
                    range.y = "relative")
oh.model$regression.results #slope = 1, no asymmetry

#ggsave(p.oh.rl, file = "./Results/operculum_height.png", width = 14, height = 10, units = "cm")

##### MAKE TABLE & SCALE CORRECT -----
#For 30x it is 0.606 pixels per um(micrometer)

traits.df <- data.frame(boxID = images.filter$box_id,
                        magnification = images.filter$Mag,
                        imageName = images.filter$image,
                        specimenNR = images.filter$specimenNR.tif,
                        formation = images.filter$formation,
                        zh = zh/.606, #z = zooid; h = height
                        oh = oh/.606, #o = operculum
                        ow.m = ow.m/.606, #w = width; m = mid
                        ow.b = ow.b/.606, #b = base
                        o.side = o.side.avg/.606,
                        mpw.b = mpw.b/.606, #pt = polypide tube
                        cw.m = cw.m/.606, #c = cryptocyst
                        cw.b = cw.b/.606,
                        cw.d = cw.d/.606, #d = distal
                        c.side = c.side.avg/.606,
                        zh.zw = zh/cw.d, #similar to LZ/WZ
                        oh.ow = oh/ow.m) # similar to LO/WO 

##### LN TRANSFORM -----
traits.df$ln.zh <- log(traits.df$zh)
traits.df$ln.mpw.b <- log(traits.df$mpw.b)
traits.df$ln.cw.m <- log(traits.df$cw.m)
traits.df$ln.cw.d <- log(traits.df$cw.d)
traits.df$ln.ow.m <- log(traits.df$ow.m)
traits.df$ln.oh <- log(traits.df$oh)
traits.df$ln.o.side <- log(traits.df$o.side)
traits.df$ln.c.side <- log(traits.df$c.side)

#write.csv(traits.df,
#          "./Results/traits_8Sept2023.csv",
#          row.names = FALSE)

traits = names(traits.df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                             "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

##### SUMMARY STATS & FILTER -----
# 5 zooid per colony minimum
traits.df$zooid.id <- paste0(traits.df$boxID, "_", traits.df$image)
colnames(traits.df)[colnames(traits.df) == 'specimenNR'] <- 'colony.id'

tr.mean_by_formation_colony = traits.df %>% #use this going forward
  dplyr::group_by(formation, colony.id) %>%
  dplyr::summarize(n.zooid = length(unique(zooid.id)),
            avg.zh = mean(ln.zh, na.rm = T),
            sd.zh = sd(ln.zh, na.rm = T),
            avg.mpw.b = mean(ln.mpw.b, na.rm = T),
            sd.mpw.b = sd(ln.mpw.b, na.rm = T),
            avg.cw.m = mean(ln.cw.m, na.rm = T),
            sd.cw.m = sd(ln.cw.m, na.rm = T),
            avg.cw.d = mean(ln.cw.d, na.rm = T),
            sd.cw.d = sd(ln.cw.d, na.rm = T),
            avg.ow.m = mean(ln.ow.m, na.rm = T),
            sd.ow.m = sd(ln.ow.m, na.rm = T),
            avg.oh = mean(ln.oh, na.rm = T),
            sd.oh = sd(ln.oh, na.rm = T),
            avg.o.side = mean(ln.o.side, na.rm = T),
            sd.o.side = sd(ln.o.side, na.rm = T),
            avg.c.side = mean(ln.c.side, na.rm = T),
            sd.c.side = sd(ln.c.side, na.rm = T)) %>%
  as.data.frame()
nrow(tr.mean_by_formation_colony) #731

keep <- tr.mean_by_formation_colony$colony.id[tr.mean_by_formation_colony$n.zooid >= 5]
length(keep) #572 colonies

df <- traits.df[traits.df$colony.id %in% keep,]
nrow(df) #5971

tr.mean_by_formation_colony.keep <- tr.mean_by_formation_colony[tr.mean_by_formation_colony$n.zooid >= 5,]

tr.mean_by_formation = df %>%
  dplyr::group_by(formation) %>%
  dplyr::summarize(num.col = length(unique(colony.id)),
                   num.zooid = length(unique(zooid.id)),
                   avg.zooid = ceiling(num.zooid/num.col), #round up to nearest integer
                   avg.zh = mean(ln.zh, na.rm = T),
                   sd.zh = sd(ln.zh, na.rm = T),
                   avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                   sd.mpw.b = sd(ln.mpw.b, na.rm = T),
                   avg.cw.m = mean(ln.cw.m, na.rm = T),
                   sd.cw.m = sd(ln.cw.m, na.rm = T),
                   avg.cw.d = mean(ln.cw.d, na.rm = T),
                   sd.cw.d = sd(ln.cw.d, na.rm = T),
                   avg.ow.m = mean(ln.ow.m, na.rm = T),
                   sd.ow.m = sd(ln.ow.m, na.rm = T),
                   avg.oh = mean(ln.oh, na.rm = T),
                   sd.oh = sd(ln.oh, na.rm = T),
                   avg.o.side = mean(ln.o.side, na.rm = T),
                   sd.o.side = sd(ln.o.side, na.rm = T),
                   avg.c.side = mean(ln.c.side, na.rm = T),
                   sd.c.side = sd(ln.c.side, na.rm = T)) %>%
  as.data.frame()

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

#colonies with little humps: MORE THAN HALF THE SIZE
#NKBS --> 17% size of biggest
#Waipuru --> 24% size of biggest
#Upper Kai-Iwi --> 24% size of biggest
#NKLS are also small, but don't see hump....
#for a size difference of 1228.705, would need a temp diff of 204.7842C??? Doesn't make sense
#based that calculation off of temp.R in Dropbox/ROCKS-PARADOX/Bryozoans

###### DISTRIBUTION ------
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

ggsave(p.dist, file = "./Results/trait_distribution_8Sept2023.png", width = 14, height = 10, units = "cm")

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
  scale_x_continuous(name = "LN Zooid Height")

#ggsave(p.zh, file = "./Results/zooid_height_8Sept2023.png", width = 14, height = 10, units = "cm")

ggplot(df) +
  geom_histogram(aes(x = ln.zh)) +
  ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Frequency") +
  scale_x_continuous(name = "LN Zooid Height")

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
  scale_x_continuous(name = "LN Zooid Height")

#ggsave(p.zh.form, file = "./Results/zooid_height_by_formation_8Sept2023.png", width = 14, height = 10, units = "cm")

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
  scale_x_continuous(name = "LN Zooid height")

#ggsave(p.ow.zh, file = "./Results/ow.zh.scaling.png", width = 14, height = 10, units = "cm")

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

##### WRITE OUT TWO DATASETS ----
write.csv(small.colonies,
          "./Results/small.colonies.traits_8Sept2023.csv",
          row.names = FALSE)

write.csv(reg.colonies,
          "./Results/colonies.traits_8Sept2023.csv",
          row.names = FALSE)

range(reg.colonies$zh)
range(small.colonies$zh)
mean(reg.colonies$zh)
median(reg.colonies$zh)
sd(reg.colonies$zh)

range(reg.colonies$zh[reg.colonies$formation == "NKLS"]) #450.4467 1701.3353
mean(reg.colonies$zh[reg.colonies$formation == "NKLS"]) #803.359
sd(reg.colonies$zh[reg.colonies$formation == "NKLS"]) #101.0176

range(reg.colonies$zh[reg.colonies$formation == "NKBS"]) #494.0033 1495.0723
range(small.colonies$zh[small.colonies$formation == "NKBS"])
mean(reg.colonies$zh[reg.colonies$formation == "NKBS"]) #778.014
sd(reg.colonies$zh[reg.colonies$formation == "NKBS"]) #93.14028
mean(small.colonies$zh[small.colonies$formation == "NKBS"]) #392.4272; diff of 385.5868
median(reg.colonies$zh[reg.colonies$formation == "NKBS"])
median(small.colonies$zh[small.colonies$formation == "NKBS"])

range(reg.colonies$zh[reg.colonies$formation == "Tewkesbury"]) #330.5689 1668.3389
mean(reg.colonies$zh[reg.colonies$formation == "Tewkesbury"]) #780.6346
sd(reg.colonies$zh[reg.colonies$formation == "Tewkesbury"]) #99.51407

range(reg.colonies$zh[reg.colonies$formation == "Waipuru"]) #544.617 1361.006
range(small.colonies$zh[small.colonies$formation == "Waipuru"])
mean(reg.colonies$zh[reg.colonies$formation == "Waipuru"]) #801.5489
sd(reg.colonies$zh[reg.colonies$formation == "Waipuru"]) #95.42124
mean(small.colonies$zh[small.colonies$formation == "Waipuru"]) #379.6889; diff of 421.86
median(reg.colonies$zh[reg.colonies$formation == "Waipuru"])
median(small.colonies$zh[small.colonies$formation == "Waipuru"])

range(reg.colonies$zh[reg.colonies$formation == "Upper Kai-Iwi"]) #561.1774 1362.1250
range(small.colonies$zh[small.colonies$formation == "Upper Kai-Iwi"])
mean(reg.colonies$zh[reg.colonies$formation == "Upper Kai-Iwi"]) #885.2919
sd(reg.colonies$zh[reg.colonies$formation == "Upper Kai-Iwi"]) #118.5598
mean(small.colonies$zh[small.colonies$formation == "Upper Kai-Iwi"]) #392.8787; diff of 492.4132
median(reg.colonies$zh[reg.colonies$formation == "Upper Kai-Iwi"])
median(small.colonies$zh[small.colonies$formation == "Upper Kai-Iwi"])

range(reg.colonies$zh[reg.colonies$formation == "Tainui"]) #592.3886 1203.1921
mean(reg.colonies$zh[reg.colonies$formation == "Tainui"]) #947.1062
sd(reg.colonies$zh[reg.colonies$formation == "Tainui"]) #95.33073

range(reg.colonies$zh[reg.colonies$formation == "SHCSBSB"]) #418.482 1333.873
mean(reg.colonies$zh[reg.colonies$formation == "SHCSBSB"]) #932.8167
sd(reg.colonies$zh[reg.colonies$formation == "SHCSBSB"]) #95.79008

####### SAMPLING NUMBERS ------

length(unique(reg.form.meta$imageName))
length(unique(reg.form.meta$colony.id))
length(unique(reg.form.meta$zooid.id))

sum.reg <- reg.form.meta %>% 
  dplyr::group_by(formation) %>%
  dplyr::summarize(n.zooid = length(unique(zooid.id)),
            n.colony = length(unique(colony.id)),
            n.image = length(unique(imageName)),
            avg.zooid.colony = ceiling(n.zooid/n.colony)) #round up to nearest integer
col.mean <- reg.form.meta %>%
  dplyr::group_by(colony.id) %>%
  dplyr::summarise(n.zooid = length(unique(zooid.id)))
ceiling(mean(col.mean$n.zooid))

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

tr.mean_by_formation.meta <- merge(tr.mean_by_formation, form.df,
                                by.x = "formation",
                                by.y = "formationCode")

## overall differnce in zooid height
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "NKLS"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "SHCSBSB"]
#-0.1514114 (decrease in length)
tr.mean_by_formation$avg.ow.m[tr.mean_by_formation$formation == "NKLS"] - tr.mean_by_formation$avg.ow.m[tr.mean_by_formation$formation == "SHCSBSB"]
#-0.1479244 (decrease in width)

##between formations
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "NKLS"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "NKBS"]
#0.1162435 (1.123269 diff)
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "NKBS"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Tewkesbury"]
#-0.08700402 (0.9166734 diff)
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Tewkesbury"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Waipuru"]
#0.06027717 (1.062131 diff)
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Waipuru"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Upper Kai-Iwi"]
#0.01439653 (1.014501 diff)
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Upper Kai-Iwi"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Tainui"]
#-0.2706406 (0.7628906 diff)
tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "Tainui"] - tr.mean_by_formation$avg.zh[tr.mean_by_formation$formation == "SHCSBSB"]
#0.01531607 (1.015434 diff)

## how is sd a function of sample size (number of zooids and number of colonies)?
#plot sd per colony by zooid no
ggplot(tr.mean_by_formation_colony) + 
  geom_point(aes(x = n.zooid, y = sd.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Number of Zooids per Colony") +
  scale_y_continuous(name = "Zooid Height SD") +
  scale_color_manual(values = col.form)

#plot sd per formation by colony no
ggplot(tr.mean_by_formation) + 
  geom_point(aes(x = num.col, y = sd.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Number of Colonies per Colony") +
  scale_y_continuous(name = "Zooid Height SD") +
  scale_color_manual(values = col.form)

anova(lm(tr.mean_by_formation$sd.zh ~ tr.mean_by_formation$num.col + tr.mean_by_formation$num.zooid))

#use formation means
#mean_by_formation

ggplot(data = tr.mean_by_formation.meta) +
  geom_point(aes(x = age.range, y = avg.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Age Range (Ma)") +
  scale_y_continuous(name = "Average Zooid Height (um)") +
  scale_color_manual(values = col.form)

ggplot(data = tr.mean_by_formation.meta) +
  geom_point(aes(x = mean.age, y = avg.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Age (Ma)") +
  scale_y_continuous(name = "Average Zooid Height (um)") +
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
  scale_y_continuous(name = "Zooid Height (um)") +
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
  scale_y_continuous(name = "Zooid Height (um)") +
  scale_color_manual(values = col.form)

box.ln.zh <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.zh",], 
       aes(x = formation, 
           y = measurementValue, 
           fill = formation)) +
  geom_boxplot() +
  scale_color_manual(values = col.form) +
  scale_fill_manual(values = col.form) +
  ggtitle("Boxplots of LN Zooid Heights") +
  xlab("Formation") +
  ylab("LN Zooid Height (um)")

ggsave(box.ln.zh, 
       file = "./Results/boxplot.ln.zh.png", 
       width = 14, height = 10, units = "cm")

#### LOOK AT TRENDS RELATIVE TO LOCALITY ----
#use colony means
#mean_by_formation_colony
list.imageName <- str_split(df$imageName, fixed("_"))
df$specimenNum <- c()
for(i in 1:length(list.imageName)){
  df$specimenNum[i] <- list.imageName[[i]][1]
}
df$specimenNum
#remove all leading 0s
df$specimenNum <- str_remove(df$specimenNum, "^0+")

length(unique(df$specimenNum)) #566

locality.df$SAMPLE_ID 
length(unique(locality.df$SAMPLE_ID)) #531
nrow(locality.df) #538
# need to remove dupes
dupe.ids <- locality.df$SAMPLE_ID[duplicated(locality.df$SAMPLE_ID)]
#"115B"  "9.44"  "118"   "119"   "161"   "9.148" "9.95" 
# don't have 115B, 9.44, 9.148, 9.95 in df
# do have 118, 119, 161 in df
dupe.rid.ids <- c("115B", "9.148", "9.95", "9.44")
dupe.keep.ids <- c("118", "119", "161")
locality.df[locality.df$SAMPLE_ID %in% dupe.keep.ids,]
#118 from WABO I and II, different formations, can probably keep NKBS and not Landguard (esp since KV collected at NKBS)
#119 from WABO I and II, different formations, can probably keep NKBS and not Landguard (esp since KV collected at NKBS)
#161 both from WABO II, different formations, can probably keep Upper Castlecliff rather than Denby

locality.df.trim <- locality.df[-c(170, 171, 241),] #eliminate based on idnex
locality.df.trim[locality.df.trim$SAMPLE_ID %in% dupe.keep.ids,]
locality.df.trim <-  locality.df.trim[!(locality.df.trim$SAMPLE_ID %in% dupe.rid.ids),]
locality.df.trim[locality.df.trim$SAMPLE_ID %in% dupe.rid.ids,] #empty, few

loc.df.trim <- locality.df.trim %>%
  dplyr::select("Expedition", "SAMPLE_ID", "Formation_name",
         "GPS.lat", "GPS.long", "Physical.Description")

df.loc <- merge(df, loc.df.trim,
                by.x = "specimenNum",
                by.y = "SAMPLE_ID",
                all.x = TRUE, all.y = FALSE)
nrow(df.loc) #5971, same number as df

missing.loc <- c(unique(df.loc$specimenNum[is.na(df.loc$GPS.lat)]), unique(df.loc$specimenNum[df.loc$GPS.lat == ""]))
length(missing.loc) #357
#300: whiterock limestone, but we have it as NKBS (from WABO III)
length(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc]) #3510

length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "NKLS"])) #65
length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "NKBS"])) #87
length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "Tewkesbury"])) #106
length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "Waipuru"])) #11
length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "Upper Kai-Iwi"])) #18
length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "Tainui"])) #19
length(unique(df.loc$specimenNum[df.loc$specimenNum %in% missing.loc &
                                   df.loc$formation == "SHCSBSB"])) #50

length(unique(df.loc$specimenNum[df.loc$GPS.lat != "" & 
                                   !is.na(df.loc$GPS.lat)])) #210 for which have info; less than half

## are the small ones in certain localities?
df.loc$size <- ""
df.loc$size[df.loc$ln.zh <= 6.25] <- "small"

nkbs.loc <- unique(df.loc$GPS.lat[df.loc$formation == "NKBS"])
nkbs.sm.loc <- unique(df.loc$GPS.lat[df.loc$formation == "NKBS" &
                                      df.loc$size == "small"]) #18 localities with small zooids
nkbs.reg.loc <- unique(df.loc$GPS.lat[df.loc$formation == "NKBS" &
                                       df.loc$size != "small"]) #18 localities with small zooids
setdiff(nkbs.sm.loc, nkbs.loc) #no diff, so all the small localities are also in the regular localities
length(nkbs.sm.loc) #18
length(nkbs.loc) #125
length(nkbs.reg.loc) #113
length(setdiff(nkbs.loc, nkbs.sm.loc)) #107 localities without small zooids
length(setdiff(nkbs.sm.loc, nkbs.reg.loc))

wai.loc <- unique(df.loc$GPS.lat[df.loc$formation == "Waipuru"])
wai.sm.loc <- unique(df.loc$GPS.lat[df.loc$formation == "Waipuru" &
                                       df.loc$size == "small"]) #1 localities with small zooids
wai.reg.loc <- unique(df.loc$GPS.lat[df.loc$formation == "Waipuru" &
                                      df.loc$size != "small"]) #1 localities with small zooids
setdiff(wai.sm.loc, wai.loc) #no diff, so all the small localities are also in the regular localities
length(wai.loc) #6
length(wai.reg.loc)
length(setdiff(wai.loc, wai.sm.loc)) #5 localities without small zooids


uki.loc <- unique(df.loc$GPS.lat[df.loc$formation == "Upper Kai-Iwi"])
uki.sm.loc <- unique(df.loc$GPS.lat[df.loc$formation == "Upper Kai-Iwi" &
                                       df.loc$size == "small"]) #4 localities with small zooids
uki.reg.loc <- unique(df.loc$GPS.lat[df.loc$formation == "Upper Kai-Iwi" &
                                      df.loc$size != "small"]) #4 localities with small zooids
setdiff(uki.sm.loc, uki.loc) #no diff, so all the small localities are also in the regular localities
length(uki.loc) #7
length(uki.reg.loc) #4
length(setdiff(uki.loc, uki.sm.loc)) #3 localities without small zooids


#### LOOK AT TRENDS RELATIVE TO SUBSTRATE ----


#### LOOK AT TRENDS RELATIVE TO TEMPERATURE ----
#use formation means
#mean_by_formation

df.form <- merge(df, form.meta,
                 by.x = "formation",
                 by.y = "formationCode",
                 all.x = TRUE,
                 all.y = FALSE)

bottom = as.numeric(df.form$Isotope_Stage_Start)
top = as.numeric(df.form$Isotope_Stage_End)
df.form$med.O18 <- c()
df.form$sd.med.O18 <- c()
df.form$n.O18 <- c()
for (i in 1:nrow(df.form)){
  temp = oxy.18$d18O[which(oxy.18$Time <= bottom[i] & oxy.18$Time >= top[i])]
  df.form$med.O18[i] = median(temp)
  df.form$sd.med.O18[i] = sd(temp)
  df.form$n.O18[i] <- length(temp)
}

unique(df.form$med.O18)

ggplot(df.form) +
  geom_point(aes(x = med.O18, y = ln.zh)) + 
  theme(text = element_text(size = 16)) +
  scale_x_continuous(name = "mean Delta O18") +
  scale_y_continuous(name = "LN Zooid Height")
#NO PATTERN

ggplot(df.form) +
  geom_point(aes(x = sd.med.O18, y = ln.zh)) + 
  theme(text = element_text(size = 16)) +
  scale_x_continuous(name = "sd Delta O18") +
  scale_y_continuous(name = "LN Zooid Height")
#NO PATTERN

# zooids scale to 0 to -6 with temp, getting smaller when warmer

#### OLD ----
noPath <- output$id

list.trim <- gsub(list,
                  pattern = "/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/from lab computer/",
                  replacement = "")

##EXPLORE DATA

##images from shared "JPG" folder

nrow(output) #16924

colnames(output)
#id
#box_id is a combination of box_top, box_left, box_width, and box_height
#series of coordinates and landmark number [e.g., (X0, Y0)]
#see "stegs_landmarks.png" for landmarks

## ID & SPECIMEN.NR
id.parse <- str_split(output$id,
                      pattern = "/")

id.only <- unlist(lapply(id.parse, function (x) x[9]))
id.imageName <- str_remove(id.only, "[^BSE]*$")
output$imageName <- id.imageName

imageName.parse <- str_split(id.imageName, fixed("_"))
specimen.NR <- c()
mag <- c()
for(i in 1:length(imageName.parse)){
  specimen.NR[i] <- paste0(imageName.parse[[i]][1], imageName.parse[[i]][2])
  mag[i] <- imageName.parse[[i]][5]
}
output$SPECIMEN.NR <- specimen.NR
nrow(output[!duplicated(output$SPECIMEN.NR),]) #777

output$mag <- mag
output$mag[output$mag == "x30.BSE"] <- "x30"
output$mag[output$mag == "x35.BSE"] <- "x35"

#names(output)[names(output)=="id"] <- "path"
#output$id <- id.only

## MATCH WITH METADATA
##there are duplicated records, want to extract unique values
#upon first inspection, it seems dupes differ by NR.OF.PICS

meta.trim <- bryo.meta[!duplicated(bryo.meta$SPECIMEN.NR), 
                       c("DATE", "SPECIMEN.NR", "FORMATION")]
nrow(bryo.meta) #880
nrow(meta.trim) #778; 102 rows removed

output.meta <- merge(output, meta.trim,
                     by = "SPECIMEN.NR",
                     all.x = TRUE, all.y = FALSE)

nrow(output) #1624
nrow(output.meta) #1624

## TWO CLUSTERS
## are the two clusters driven by formation? - NO
ggplot(data = traits.df) +
  geom_smooth(aes(x = zh, y = ow.m), method = "lm") +
  geom_point(aes(x = zh, y = ow.m,
                 group = formation,
                 col = formation)) + #two clusters
  ggtitle("Scaling of operculum mid-width with zooid height, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum mid-width (pixels)") +
  scale_x_continuous(name = "Zooid height (pixels)")

## are they driven by differences in magnification? - NO
ggplot(data = traits.df) +
  geom_smooth(aes(x = zh, y = ow.m), method = "lm") +
  geom_point(aes(x = zh, y = ow.m,
                 group = magnification,
                 col = magnification,
                 alpha = .3)) + #two clusters
  ggtitle("Scaling of operculum mid-width with zooid height, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum mid-width (pixels)") +
  scale_x_continuous(name = "Zooid height (pixels)")

#### DIFFERENCES IN TRAIT MEANS ----
traits
formations <- c("NKLS", "NKBS", "Tewkesbury", 
                "Waipuru", "Upper Kai-Iwi", 
                "Tainui", "SHCSBSB")
## t tests
#ln zh
t.test(traits.df[traits.df$formation == formations[1], traits[1]],
       traits.df[traits.df$formation == formations[2], traits[1]])
t.test(traits.df[traits.df$formation == formations[2], traits[1]],
       traits.df[traits.df$formation == formations[3], traits[1]])
t.test(traits.df[traits.df$formation == formations[3], traits[1]],
       traits.df[traits.df$formation == formations[4], traits[1]])
t.test(traits.df[traits.df$formation == formations[4], traits[1]],
       traits.df[traits.df$formation == formations[5], traits[1]])
t.test(traits.df[traits.df$formation == formations[5], traits[1]],
       traits.df[traits.df$formation == formations[6], traits[1]])
t.test(traits.df[traits.df$formation == formations[6], traits[1]],
       traits.df[traits.df$formation == formations[7], traits[1]])

#ln.mpw.b
t.test(traits.df[traits.df$formation == formations[1], traits[2]],
       traits.df[traits.df$formation == formations[2], traits[2]])
t.test(traits.df[traits.df$formation == formations[2], traits[2]],
       traits.df[traits.df$formation == formations[3], traits[2]])
t.test(traits.df[traits.df$formation == formations[3], traits[2]],
       traits.df[traits.df$formation == formations[4], traits[2]])
t.test(traits.df[traits.df$formation == formations[4], traits[2]],
       traits.df[traits.df$formation == formations[5], traits[2]])
t.test(traits.df[traits.df$formation == formations[5], traits[2]],
       traits.df[traits.df$formation == formations[6], traits[2]])
t.test(traits.df[traits.df$formation == formations[6], traits[2]],
       traits.df[traits.df$formation == formations[7], traits[2]])

#ln.cw.m
t.test(traits.df[traits.df$formation == formations[1], traits[3]],
       traits.df[traits.df$formation == formations[2], traits[3]])
t.test(traits.df[traits.df$formation == formations[2], traits[3]],
       traits.df[traits.df$formation == formations[3], traits[3]])
t.test(traits.df[traits.df$formation == formations[3], traits[3]],
       traits.df[traits.df$formation == formations[4], traits[3]])
t.test(traits.df[traits.df$formation == formations[4], traits[3]],
       traits.df[traits.df$formation == formations[5], traits[3]])
t.test(traits.df[traits.df$formation == formations[5], traits[3]],
       traits.df[traits.df$formation == formations[6], traits[3]])
t.test(traits.df[traits.df$formation == formations[6], traits[3]],
       traits.df[traits.df$formation == formations[7], traits[3]])

#ln.cw.d
t.test(traits.df[traits.df$formation == formations[1], traits[4]],
       traits.df[traits.df$formation == formations[2], traits[4]])
t.test(traits.df[traits.df$formation == formations[2], traits[4]],
       traits.df[traits.df$formation == formations[3], traits[4]])
t.test(traits.df[traits.df$formation == formations[3], traits[4]],
       traits.df[traits.df$formation == formations[4], traits[4]])
t.test(traits.df[traits.df$formation == formations[4], traits[4]],
       traits.df[traits.df$formation == formations[5], traits[4]])
t.test(traits.df[traits.df$formation == formations[5], traits[4]],
       traits.df[traits.df$formation == formations[6], traits[4]])
t.test(traits.df[traits.df$formation == formations[6], traits[4]],
       traits.df[traits.df$formation == formations[7], traits[4]])

#ln.ow.m
t.test(traits.df[traits.df$formation == formations[1], traits[5]],
       traits.df[traits.df$formation == formations[2], traits[5]])
t.test(traits.df[traits.df$formation == formations[2], traits[5]],
       traits.df[traits.df$formation == formations[3], traits[5]])
t.test(traits.df[traits.df$formation == formations[3], traits[5]],
       traits.df[traits.df$formation == formations[4], traits[5]])
t.test(traits.df[traits.df$formation == formations[4], traits[5]],
       traits.df[traits.df$formation == formations[5], traits[5]])
t.test(traits.df[traits.df$formation == formations[5], traits[5]],
       traits.df[traits.df$formation == formations[6], traits[5]])
t.test(traits.df[traits.df$formation == formations[6], traits[5]],
       traits.df[traits.df$formation == formations[7], traits[5]])

#ln.oh
t.test(traits.df[traits.df$formation == formations[1], traits[6]],
       traits.df[traits.df$formation == formations[2], traits[6]])
t.test(traits.df[traits.df$formation == formations[2], traits[6]],
       traits.df[traits.df$formation == formations[3], traits[6]])
t.test(traits.df[traits.df$formation == formations[3], traits[6]],
       traits.df[traits.df$formation == formations[4], traits[6]])
t.test(traits.df[traits.df$formation == formations[4], traits[6]],
       traits.df[traits.df$formation == formations[5], traits[6]])
t.test(traits.df[traits.df$formation == formations[5], traits[6]],
       traits.df[traits.df$formation == formations[6], traits[6]])
t.test(traits.df[traits.df$formation == formations[6], traits[6]],
       traits.df[traits.df$formation == formations[7], traits[6]])

#ln.c.side
t.test(traits.df[traits.df$formation == formations[1], traits[7]],
       traits.df[traits.df$formation == formations[2], traits[7]])
t.test(traits.df[traits.df$formation == formations[2], traits[7]],
       traits.df[traits.df$formation == formations[3], traits[7]])
t.test(traits.df[traits.df$formation == formations[3], traits[7]],
       traits.df[traits.df$formation == formations[4], traits[7]])
t.test(traits.df[traits.df$formation == formations[4], traits[7]],
       traits.df[traits.df$formation == formations[5], traits[7]])
t.test(traits.df[traits.df$formation == formations[5], traits[7]],
       traits.df[traits.df$formation == formations[6], traits[7]])
t.test(traits.df[traits.df$formation == formations[6], traits[7]],
       traits.df[traits.df$formation == formations[7], traits[7]])

#ln.o.side
t.test(traits.df[traits.df$formation == formations[1], traits[8]],
       traits.df[traits.df$formation == formations[2], traits[8]])
t.test(traits.df[traits.df$formation == formations[2], traits[8]],
       traits.df[traits.df$formation == formations[3], traits[8]])
t.test(traits.df[traits.df$formation == formations[3], traits[8]],
       traits.df[traits.df$formation == formations[4], traits[8]])
t.test(traits.df[traits.df$formation == formations[4], traits[8]],
       traits.df[traits.df$formation == formations[5], traits[8]])
t.test(traits.df[traits.df$formation == formations[5], traits[8]],
       traits.df[traits.df$formation == formations[6], traits[8]])
t.test(traits.df[traits.df$formation == formations[6], traits[8]],
       traits.df[traits.df$formation == formations[7], traits[8]])
