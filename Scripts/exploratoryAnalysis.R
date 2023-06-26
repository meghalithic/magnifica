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
images.meta <- read.csv("./Data/meta.images.22Jun2023.csv",
                   header = TRUE,
                   sep = ",")

df.filter <- read.table("./Data/filteredImages.csv",
                        header = TRUE,
                        sep = ";")

#### EXPLORE DATA ----
nrow(df.filter) #1834
nrow(images.meta) #15783

#### REDUCE TO SELECTED IMAGES ####

df.filter$fileName.old <- c()
for(i in 1:nrow(df.filter)){
  df.filter$fileName.old[i] <- str_split(df.filter$path.tif[i], "/")[[1]][length(str_split(df.filter$path.tif[i], "/")[[1]])]
}

images.filter <- images.meta[images.meta$fileName.tif %in% df.filter$fileName.old,]
nrow(images.filter) #15773
length(unique(images.filter$fileName.tif)) #1464; only 370 images removed

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

plot(o.side.l, o.side.r)
summary(lm(o.side.r ~ o.side.l)) #basically 1

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

plot(c.side.l, c.side.r)
summary(lm(c.side.r ~ c.side.l)) #basically 1

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

write.csv(traits.df,
          "./Results/traits_26Jun2023.csv",
          row.names = FALSE)

##### ABOUT TRAITS -----

nrow(images.filter) #15773
nrow(traits.df) #15773

traits.melt <- melt(data = traits.df,
                    id.vars = c("boxID","imageName", "specimenNR", "formation", "magnification"),
                    variable.name = "measurementType",
                    value.name = "measurementValue")
length(unique(traits.melt$specimenNR)) #742 unique colonies

traits.stats <- traits.melt %>%
  group_by(measurementType) %>%
  summarise(avg = mean(measurementValue))

##### HISTOGRAM -----

p.dist <- ggplot(traits.melt) +
  geom_density(aes(x = log10(measurementValue),
                   group = measurementType,
                   col = measurementType)) + #lots are bimodal
  ggtitle("Distribution of traits, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "log10 trait measurement (pixels)")

#ggsave(p.dist, file = "./Results/trait_distribution_22Jun2023.png", width = 14, height = 10, units = "cm")

###### BIMODALITY ------  
##explore bimodality, using zooid height as an example then see if it generalizes
p.zh <- ggplot(traits.df) +
  geom_density(aes(x = zh)) +
  ggtitle("Zooid height, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Zooid Height (pixels)")

#ggsave(p.zh, file = "./Results/zooid_height.png", width = 14, height = 10, units = "cm")

##driven by formation?
#all but Punneki Limestone are bimodal
p.zh.form <- ggplot(traits.df) +
  geom_density(aes(x = zh,
                   group = formation,
                   col = formation)) + 
  ggtitle("Zooid height by formation, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Zooid Height (pixels)")

#ggsave(p.zh.form, file = "./Results/zooid_height_by_formation.png", width = 14, height = 10, units = "cm")

#only see bimodal in Waipurpu, NKBS, Upper Kai-Iwi

##seems like first hump ends around 300 pixels
sm.traits <- traits.df[traits.df$zh < 500,]
length(unique(sm.traits$specimenNR)) #95 images out of 891

write.csv(sm.traits,
          "./Results/small_bimodal_hump.csv",
          row.names = FALSE)

##driven by magnification?
#not that I can tell...
p.zh.mag <- ggplot(traits.df) +
  geom_density(aes(x = zh,
                   group = magnification,
                   col = magnification)) + 
  ggtitle("Zooid height by magnification, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Zooid Height (pixels)")

#ggsave(p.zh.mag, file = "./Results/zooid_height_by_magnification.png", width = 14, height = 10, units = "cm")
##doesn't seem to affect it, still get it with mag == 30 only

##what does it look like if all magnification is the same?
#still get bimodal hump

nrow(traits.df[traits.df$magnification == "x30",]) #15773
length(unique(traits.df$specimenNR[traits.df$magnification == "x30"])) #742

p.zh.mag.30 <- ggplot(traits.df[traits.df$magnification == "x30",]) +
  geom_density(aes(x = zh[traits.df$magnification == "x30"])) + 
  ggtitle("Zooid height by magnification x30, N zooids = 18752, N colony = 884") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Zooid Height (pixels)")

#ggsave(p.zh.mag.30, file = "./Results/zooid_height_magnification_x30.png", width = 14, height = 10, units = "cm")

##### SCALING -----

## right v left side of operculum
p.oh.rl <- ggplot(data = traits.df) +
  geom_point(aes(x = oh.l, y = oh.r)) +
  ggtitle("Operculum Height, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum Height Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Height Left Side (pixels)")

oh.model <- lmodel2(formula = oh.r ~ oh.l,
                    data = traits.df,
                    range.x = "relative", 
                    range.y = "relative")
oh.model$regression.results #slope = 1, no asymmetry

#ggsave(p.oh.rl, file = "./Results/operculum_height.png", width = 14, height = 10, units = "cm")

##against zooid height

p.ow.zh <- ggplot(data = traits.df) +
  geom_smooth(aes(x = zh, y = ow.m), method = "lm") +
  geom_point(aes(x = zh, y = ow.m)) + #two clusters
  ggtitle("Scaling of operculum mid-width with zooid height, N zooids = 18890, N colony = 891") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum mid-width (pixels)") +
  scale_x_continuous(name = "Zooid height (pixels)")

#ggsave(p.ow.zh, file = "./Results/ow.zh.scaling.png", width = 14, height = 10, units = "cm")

#look at 4 images:
#2 with the largest operculum width
#2 with the largest zooid height

slice_max(traits.df, n = 2, order_by = ow.m)
#ow.m = 1061.6501; boxID = 290_1913_523_687; imageName = 806_CV_2_10v_x30_BSE
#ow.m = 959.4425; boxID = 927_2047_495_737; imageName = 768_CC_2_15v_x30_BSE

slice_max(traits.df, n = 2, order_by = zh)
#zh = 1498.012; boxID = 301_1136_705_1252; imageName = 511_CV_1_15v_x30_BSE
#zh = 1363.376; boxID = 21_769_409_746; imageName = 084_CV_2_10v_x30_BSE

zh.owm.model <- lm(ow.m ~ zh,
                   data = traits.df)






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