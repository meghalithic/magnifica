## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)
require(ggplot2)
require(reshape2)
require(lmodel2)
require(tidyverse)

#### LOAD DATA ----
output <- read.csv("./Data/output.csv", header = TRUE)
AP_images <- read.csv("./Data/images_from_AP.csv", header = TRUE)

#### EXPLORE DATA ----

nrow(AP_images) #1654
imageName.parse_AP <- str_split(AP_images$imageName, fixed("_"))
specimen.NR_AP <- c()
for(i in 1:length(imageName.parse_AP)){
  specimen.NR_AP[i] <- paste0(imageName.parse_AP[[i]][1], imageName.parse_AP[[i]][2])
}
AP_images$SPECIMEN.NR <- specimen.NR_AP
nrow(AP_images[!duplicated(AP_images$SPECIMEN.NR),]) #779, two more than in bryo metadata

nrow(output) #16924

colnames(output)
#id
#box_id is a combination of box_top, box_left, box_width, and box_height
#series of coordinates and landmark number [e.g., (X0, Y0)]
#see "stegs_landmarks.png" for landmarks

##### ID & SPECIMEN.NR -----
id.parse <- str_split(output$id,
                      pattern = "/")

id.only <- unlist(lapply(id.parse, function (x) x[9]))
id.imageName <- str_remove(id.only, "[^BSE]*$")
output$imageName <- id.imageName

#names(output)[names(output)=="id"] <- "path"
#output$id <- id.only

#### MATCH WITH METADATA ----
bryo.meta <- read.csv("./Data/Imaged Steginoporella magnifica specimens.csv",
                      header = TRUE)

##there are duplicated records, want to extract unique values
#upon first inspection, it seems dupes differ by NR.OF.PICS

meta.trim <- bryo.meta[!duplicated(bryo.meta$SPECIMEN.NR), c("DATE", "SPECIMEN.NR", "FORMATION")]
nrow(bryo.meta) #880
nrow(meta.trim) #777; 103 rows removed

imageName.parse <- str_split(id.imageName, fixed("_"))
specimen.NR <- c()
for(i in 1:length(imageName.parse)){
  specimen.NR[i] <- paste0(imageName.parse[[i]][1], imageName.parse[[i]][2], "_",
                           imageName.parse[[i]][])
}
output$SPECIMEN.NR <- specimen.NR
nrow(output[!duplicated(output$SPECIMEN.NR),]) #777 #two less than in AP_images, same number as bryo metadata


output.meta <- merge(output, meta.trim,
                     by = "SPECIMEN.NR",
                     all.x = TRUE, all.y = FALSE)

nrow(output) #1624
nrow(output.meta) #1624

#### CALCULATE DISTANCES ----
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
zh.x <- abs((output$X4-output$X12))
zh.y <-  abs((output$Y4-output$Y12))
zh <- sqrt(((zh.x)^2 + (zh.y)^2))

## Operculum height of left side: 4 to 20
oh.l.x <- abs((output$X4-output$X20))
oh.l.y <-  abs((output$Y4-output$Y20))
oh.l <- sqrt(((oh.l.x)^2 + (oh.l.y)^2))

## Operculum height of right side: 4 to 21
oh.r.x <- abs((output$X4-output$X20))
oh.r.y <-  abs((output$Y4-output$Y20))
oh.r <- sqrt(((oh.r.x)^2 + (oh.r.y)^2))

## Operculum mid-width (maximum width at centerline): 19 to 0
ow.m.x <- abs((output$X19-output$X0))
ow.m.y <-  abs((output$Y19-output$Y0))
ow.m <- sqrt(((ow.m.x)^2 + (ow.m.y)^2))

## Operculum base width: 21 to 20
ow.b.x <- abs((output$X21-output$X20))
ow.b.y <-  abs((output$Y21-output$Y20))
ow.b <- sqrt(((ow.b.x)^2 + (ow.b.y)^2))

## Median process base width: 5 to 6
mpw.b.x <- abs((output$X5-output$X6))
mpw.b.y <-  abs((output$Y5-output$Y6))
mpw.b <- sqrt(((mpw.b.x)^2 + (mpw.b.y)^2))

## Cryptocyst mid-width: 10 to 11
cw.m.x <- abs((output$X10-output$X11))
cw.m.y <-  abs((output$Y10-output$Y11))
cw.m <- sqrt(((cw.m.x)^2 + (cw.m.y)^2))

## Cryptocyst base width: 9 to 1
cw.b.x <- abs((output$X9-output$X1))
cw.b.y <-  abs((output$Y9-output$Y1))
cw.b <- sqrt(((cw.b.x)^2 + (cw.b.y)^2))

## Cryptocyst distal width: 8 to 7
cw.d.x <- abs((output$X8-output$X7))
cw.d.y <-  abs((output$Y8-output$Y7))
cw.d <- sqrt(((cw.d.x)^2 + (cw.d.y)^2))

traits.df <- data.frame(boxID = output.meta$box_id,
                        imageName = output.meta$imageName,
                        specimenNR = output.meta$SPECIMEN.NR,
                        formation = output.meta$FORMATION,
                        dates = output.meta$DATE,
                        zh = zh, #z = zooid; h = height
                        oh.l = oh.l, #o = operculum; l = left
                        oh.r = oh.r, #r = right
                        ow.m = ow.m, #w = width; m = mid
                        ow.b = ow.b, #b = base
                        mpw.b = mpw.b, #pt = polypide tube
                        cw.m = cw.m, #c = cryptocyst
                        cw.b = cw.b,
                        cw.d = cw.d) #d = distal

##### HISTOGRAM -----

traits.melt <- melt(data = traits.df,
                    id.vars = c("imageName", "specimenNR", "formation", "dates"),
                    variable.name = "measurementType",
                    value.name = "measurementValue")

ggplot(traits.melt) +
  geom_density(aes(x = measurementValue,
                   group = measurementType,
                   col = measurementType)) + #lots are bimodal
  ggtitle("Distribution of traits, N = 16924") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Trait measurement (pixels)")
  
##explore bimodality, using zooid height as an example then see if it generalizes
p <- ggplot(traits.df) +
  geom_density(aes(x = zh)) +
  ggtitle("Zooid height, N = 16924") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Zooid Height (pixels)")

ggsave(p, file = "./Results/zooid_height.png", width = 14, height = 10, units = "cm")

p <- ggplot(traits.df) +
  geom_density(aes(x = zh,
                   group = formation,
                   col = formation)) + #all but Punneki Limestone are bimodal
  ggtitle("Zooid height by formation, N = 16924") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Zooid Height (pixels)")

ggsave(p, file = "./Results/zooid_height_by_formation.png", width = 14, height = 10, units = "cm")


##seems like first hump ends around 300 pixels
sm.traits <- traits.df[traits.df$zh < 300,]
length(unique(sm.traits$specimenNR)) #174 images out of 778

##### ABOUT TRAITS -----

nrow(output) #16924
nrow(traits.df) #16924

traits.stats <- traits.melt %>%
  group_by(measurementType) %>%
  summarise(avg = mean(measurementValue))

# # A tibble: 9 × 2
# measurementType   avg
# <fct>           <dbl>
# zh               480.
# oh.l             265.
# oh.r             265.
# ow.m             259.
# ow.b             211.
# mpw.b            107.
# cw.m             240.
# cw.b             139.
# cw.d             261.

##expect oh.l and oh.r to be similar - YES
##expect zh to be largest value - YES
##expect mpw.b to be smallest value - YES

##### SCALING -----

## right v left side of operculum
p <- ggplot(data = traits.df) +
  geom_point(aes(x = oh.l, y = oh.r)) +
  ggtitle("Operculum Height, N = 16924") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum Height Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Height Left Side (pixels)")

oh.model <- lmodel2(formula = oh.r ~ oh.l,
                    data = traits.df,
                    range.x = "relative", 
                    range.y = "relative")
oh.model$regression.results #slope = 1, no asymmetry

ggsave(p, file = "./Results/operculum_height.png", width = 14, height = 10, units = "cm")

##against zooid height

ggplot(data = traits.df) +
  geom_smooth(aes(x = zh, y = ow.m), method = "lm") +
  geom_point(aes(x = zh, y = ow.m)) + #two clusters
  ggtitle("Scaling of operculum mid-width with zooid height, N = 16924") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum mid-width (pixels)") +
  scale_x_continuous(name = "Zooid height (pixels)")

#look at 4 images:
#2 with the largest operculum width
#2 with the largest zooid height

slice_max(traits.df, n = 2, order_by = ow.m)
#777_CV_1_15v_x30_BSE, boxID = 462_1457_541_660
#015i_CV_2_15v_x30, boxID = 741_1414_385_571

slice_max(traits.df, n = 2, order_by = zh)
#665_CV_3_10v_x30_BSE, boxID = 962_1180_457_606
#458_CV_1_15v_x30_BSE, boxID = 942_1361_445_585

zh.owm.model <- lm(ow.m ~ zh,
                   data = traits.df)

## are the two clusters driven by formation? - NO
ggplot(data = traits.df) +
  geom_smooth(aes(x = zh, y = ow.m), method = "lm") +
  geom_point(aes(x = zh, y = ow.m,
                 group = formation,
                 col = formation)) + #two clusters
  ggtitle("Scaling of operculum mid-width with zooid height, N = 16924") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum mid-width (pixels)") +
  scale_x_continuous(name = "Zooid height (pixels)")

