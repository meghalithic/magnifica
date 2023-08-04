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

#### LOAD DATA ----
images.meta <- read.csv("./Data/meta.images.4Aug2023.csv",
                   header = TRUE,
                   sep = ",")

df.filter <- read.table("./Data/filteredImages.csv",
                        header = TRUE,
                        sep = ";")

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

#ggsave(p.o.side.rl, file = "./Results/operculum_length.png", width = 14, height = 10, units = "cm")

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

#ggsave(p.c.side.rl, file = "./Results/cryptocyst_length.png", width = 14, height = 10, units = "cm")

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

#write.csv(traits.df,
#          "./Results/traits_4Aug2023.csv",
#          row.names = FALSE)

##### LN TRANSFORM -----
traits.df$ln.zh <- log(traits.df$zh)
traits.df$ln.mpw.b <- log(traits.df$mpw.b)
traits.df$ln.cw.m <- log(traits.df$cw.m)
traits.df$ln.cw.d <- log(traits.df$cw.d)
traits.df$ln.ow.m <- log(traits.df$ow.m)
traits.df$ln.oh <- log(traits.df$oh)
traits.df$ln.o.side <- log(traits.df$o.side)
traits.df$ln.c.side <- log(traits.df$c.side)

traits = names(traits.df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                             "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

##### FILTER -----
# 5 zooid per colony minimum
traits.df$zooid.id <- paste0(traits.df$boxID, "_", traits.df$image)
colnames(traits.df)[colnames(traits.df) == 'specimenNR'] <- 'colony.id'

mean_by_formation_colony = traits.df %>% #use this going forward
  group_by(formation, colony.id) %>%
  summarize(n.zooid = length(unique(zooid.id)),
            avg.zh = mean(ln.zh, na.rm = T),
            avg.mpw.b = mean(ln.mpw.b, na.rm = T),
            avg.cw.m = mean(ln.cw.m, na.rm = T),
            avg.cw.d = mean(ln.cw.d, na.rm = T),
            avg.ow.m = mean(ln.ow.m, na.rm = T),
            avg.oh = mean(ln.oh, na.rm = T),
            avg.o.side = mean(ln.o.side, na.rm = T),
            avg.c.side = mean(ln.c.side, na.rm = T)) %>%
  as.data.frame()
nrow(mean_by_formation_colony) #731

keep <- mean_by_formation_colony$colony.id[mean_by_formation_colony$n.zooid >= 5]
length(keep) #572 colonies

df <- traits.df[traits.df$colony.id %in% keep,]
nrow(df) #5971

##### ABOUT TRAITS -----

nrow(images.filter) #6438
nrow(df) #5971

traits.melt <- melt(data = df,
                    id.vars = c("boxID", "zooid.id","imageName", "colony.id", "formation", "magnification"),
                    variable.name = "measurementType",
                    value.name = "measurementValue")
length(unique(traits.melt$colony.id)) #572 unique colonies [previously 742]

traits.stats <- traits.melt %>%
  group_by(measurementType) %>%
  summarise(avg = mean(measurementValue))

traits.stats.form <- traits.melt %>%
  group_by(measurementType, formation) %>%
  summarise(avg = mean(measurementValue))

##### HISTOGRAM -----

p.dist <- ggplot(traits.melt) +
  geom_density(aes(x = log(measurementValue),
                   group = measurementType,
                   col = measurementType)) + #lots are bimodal
  ggtitle("Distribution of traits, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN trait measurement (pixels)")

ggsave(p.dist, file = "./Results/trait_distribution_4Aug2023.png", width = 14, height = 10, units = "cm")

###### BIMODALITY ------  
##explore bimodality, using zooid height as an example then see if it generalizes
p.ln.zh <- ggplot(df) +
  geom_density(aes(x = ln.zh)) +
  ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN Zooid Height")

#ggsave(p.zh, file = "./Results/zooid_height_4Aug2023.png", width = 14, height = 10, units = "cm")

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

#ggsave(p.zh.form, file = "./Results/zooid_height_by_formation_4Aug2023.png", width = 14, height = 10, units = "cm")

#only see bimodal in Waipurpu, NKBS, Upper Kai-Iwi

##driven by MAGNIFICATION?
#note: previously had different magnifications, this is not the case anymore
#not that I can tell...
p.zh.mag <- ggplot(df) +
  geom_density(aes(x = ln.zh,
                   group = magnification,
                   col = magnification)) + 
  ggtitle(paste0("Zooid height by magnification, N zooids = ", nrow(df), ", N colony = ", length(keep))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "LN Zooid Height")

#ggsave(p.zh.mag, file = "./Results/zooid_height_by_magnification.png", width = 14, height = 10, units = "cm")
##doesn't seem to affect it, still get it with mag == 30 only

##what does it look like if all magnification is the same?
#still get bimodal hump

nrow(df[df$magnification == "x30",]) #6815
length(unique(df$colony.id[df$magnification == "x30"])) #604

##SAME INDIVIDUALS & COLONIES?
## really convince self that these are the same individuals
# make data more manageable by reducing it to the three formations
# make data even more manageable by reducing it to the small hump that is seen
nk.wa.uki <- df[df$formation == "NKBS" |
                df$formation == "Waipuru" |
                df$formation == "Upper Kai-Iwi",]
sm.zoo <- nk.wa.uki[nk.wa.uki$ln.zh <= 6.25,]
nrow(sm.zoo) #543 (was 943) zooids
length(unique(sm.zoo$colony.id)) #36 (was 64) colonies
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
imgj.res_077_1 <- read.csv("./Data/ImgJ/077_CV_1_Results.csv", header = TRUE)
imgj.res_077_2 <- read.csv("./Data/ImgJ/077_CV_2_Results.csv", header = TRUE)
imgj.res_077_3 <- read.csv("./Data/ImgJ/077_CV_3_Results.csv", header = TRUE)
imgj_077 <- rbind(imgj.res_077_1, imgj.res_077_2, imgj.res_077_3)
## scale
imgj_077$scale.zh <- imgj_077$Length/.606
## log
imgj_077$ln.zh <- log(imgj_077$scale.zh)
imgj_077$ln.zh # all small!!

##273S
unique(sm.zoo$imageName[sm.zoo$colony.id == "273S"]) #from images 1, 2, 3
imgj.res_273_1 <- read.csv("./Data/ImgJ/273_S_1_Results.csv", header = TRUE)
imgj.res_273_2 <- read.csv("./Data/ImgJ/273_S_2_Results.csv", header = TRUE)
imgj.res_273_3 <- read.csv("./Data/ImgJ/273_S_3_Results.csv", header = TRUE)
imgj_273 <- rbind(imgj.res_273_1, imgj.res_273_2, imgj.res_273_3)
## scale
imgj_273$scale.zh <- imgj_273$Length/.606
## log
imgj_273$ln.zh <- log(imgj_273$scale.zh)
imgj_273$ln.zh #all small!

#### SMALL HUMP ----
# know that these are all related, use zh as test
## test where hump is...
# frequency by size bin (quarter ln bins)
df.bins <- df %>% 
  mutate(zh.bin = cut(ln.zh, breaks = seq(5.5, 7.5, .1))) %>%
  as.data.frame()

df.bin.f <- df.bins %>%
  group_by(zh.bin) %>%
  summarise(n = n()) %>%
  as.data.frame()
View(df.bin.f)
#(6.2,6.3]
#6.25 like I eyeballed
write.csv(df.bin.f,
          "./Results/zh.bin.frequency.csv",
          row.names = FALSE)

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
  group_by(colony.id) %>%
  summarise(n.zooid = length(zooid.id),
            n.sm.zooid = sum(sm),
            prop.sm = n.sm.zooid/n.zooid)
View(prop.sm)
write.csv(prop.sm,
          "./Results/proportion.small.colonies.csv",
          row.names = FALSE)

#goes from 100% to 20%; perhaps make 20% the cut off
rm.col <- prop.sm$colony.id[prop.sm$prop.sm == 1]
length(rm.col) #31 colonies

small.colonies <- df[df$colony.id %in% rm.col,]
reg.colonies <- df[!(df$colony.id %in% rm.col),]

#for those with 20%, see where the sizes are
sm.zooids <- prop.sm$colony.id[prop.sm$prop.sm < 1 &
                                 prop.sm$prop.sm > 0]
length(sm.zooids) #19

range(df$ln.zh[df$colony.id %in% sm.zooids])
#5.800815 7.439169
sort(df$ln.zh[df$colony.id %in% sm.zooids])
#only 11 below 6.25, so probably fine

##### WRITE OUT TWO DATASETS ----
write.csv(small.colonies,
          "./Results/small.colonies.traits.csv",
          row.names = FALSE)

write.csv(reg.colonies,
          "./Results/colonies.traits.csv",
          row.names = FALSE)

##### CORRELATIONS & ALLOMETRIES -----
## are these coming from the same individuals??
## ask KLV for other metadata for sites
## OH hump is on the other side

col.form = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
             "#00A9FF", "#C77CFF", "#FF61CC")


col.traits = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
               "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")

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
