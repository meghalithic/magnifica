## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) creates trait measurements

## The output of this code is a dataset of all the traits, traits.df, saved
# as traits_date.csv

#### LOAD DATA ----

source("./Scripts/0-env.R")

images.df <- read.csv("./Data/images.filtered_26Feb2024.csv", #images.merged_30Nov2023.csv,
                      header = TRUE, 
                      sep = ",")

#### EXPLORE DATA ----
nrow(images.df) #6649

#### MANIPULATE DATA ####

images.df$formation <- factor(images.df$formation, 
                              levels = c("NKLS", "NKBS", "Tewkesbury",
                                         "Upper Kai-Iwi", "Tainui",
                                         "SHCSBSB", "modern")) 

colnames(images.df)[colnames(images.df) == 'new.id'] <- 'zooid.id'
colnames(images.df)[colnames(images.df) == 'specimenNR'] <- 'colony.id'

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
zh.x <- abs((images.df$X4-images.df$X12))
zh.y <-  abs((images.df$Y4-images.df$Y12))
zh <- sqrt(((zh.x)^2 + (zh.y)^2))

## Operculum height of left side: 4 to 20
oh.l.x <- abs((images.df$X4-images.df$X20))
oh.l.y <-  abs((images.df$Y4-images.df$Y20))
oh.l <- sqrt(((oh.l.x)^2 + (oh.l.y)^2))

## Operculum height of right side: 4 to 21
oh.r.x <- abs((images.df$X4-images.df$X20))
oh.r.y <-  abs((images.df$Y4-images.df$Y20))
oh.r <- sqrt(((oh.r.x)^2 + (oh.r.y)^2))

## Operculum mid-width (maximum width at centerline): 19 to 0
ow.m.x <- abs((images.df$X19-images.df$X0))
ow.m.y <-  abs((images.df$Y19-images.df$Y0))
ow.m <- sqrt(((ow.m.x)^2 + (ow.m.y)^2))

## Operculum base width: 21 to 20
ow.b.x <- abs((images.df$X21-images.df$X20))
ow.b.y <-  abs((images.df$Y21-images.df$Y20))
ow.b <- sqrt(((ow.b.x)^2 + (ow.b.y)^2))

## Operculum side length of right side: 21 to 18
o.side.r.x <- abs((images.df$X21-images.df$X18))
o.side.r.y <- abs((images.df$Y21-images.df$Y18))
o.side.r <- sqrt(((o.side.r.x)^2 + (o.side.r.y)^2))

## Operculum side length of left side: 20 to 17
o.side.l.x <- abs((images.df$X20-images.df$X17))
o.side.l.y <- abs((images.df$Y20-images.df$Y17))
o.side.l <- sqrt(((o.side.l.x)^2 + (o.side.l.y)^2))

## Operculum average side length:
o.side.avg <- rowMeans(cbind(o.side.r, o.side.l))

p.o.side.rl <- ggplot() +
  geom_point(aes(x = o.side.l, y = o.side.r)) +
  ggtitle(paste0("Operculum Side Length, N zooids = ", nrow(images.df), ", N colony = ", length(unique(images.df$colony.id)))) +
  plot.theme +
  scale_y_continuous(name = "Operculum Length Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Length Left Side (pixels)")

ggsave(p.o.side.rl, 
       file = "./Results/operculum.length.w.modern.png", width = 14, height = 10, units = "cm")

summary(lm(o.side.r ~ o.side.l)) 
# slope = 0.91; p-value < 2.2e-16; r2 = 0.95

## Operculum height
oh <- (.5/ow.b)*sqrt(ow.b+oh.r+oh.l)

## Median process base width: 5 to 6
mpw.b.x <- abs((images.df$X5-images.df$X6))
mpw.b.y <-  abs((images.df$Y5-images.df$Y6))
mpw.b <- sqrt(((mpw.b.x)^2 + (mpw.b.y)^2))

## Cryptocyst mid-width: 10 to 11
cw.m.x <- abs((images.df$X10-images.df$X11))
cw.m.y <-  abs((images.df$Y10-images.df$Y11))
cw.m <- sqrt(((cw.m.x)^2 + (cw.m.y)^2))

## Cryptocyst base width: 9 to 1
cw.b.x <- abs((images.df$X9-images.df$X1))
cw.b.y <-  abs((images.df$Y9-images.df$Y1))
cw.b <- sqrt(((cw.b.x)^2 + (cw.b.y)^2))

## Cryptocyst distal width: 8 to 7
cw.d.x <- abs((images.df$X8-images.df$X7))
cw.d.y <-  abs((images.df$Y8-images.df$Y7))
cw.d <- sqrt(((cw.d.x)^2 + (cw.d.y)^2))

## Cryptocyst side length of right side: 1 to 7
c.side.r.x <- abs((images.df$X1-images.df$X7))
c.side.r.y <- abs((images.df$Y1-images.df$Y7))
c.side.r <- sqrt(((c.side.r.x)^2 + (c.side.r.y)^2))

## Cryptocyst side length of left side: 9 to 8
c.side.l.x <- abs((images.df$X9-images.df$X8))
c.side.l.y <- abs((images.df$Y9-images.df$Y8))
c.side.l <- sqrt(((c.side.l.x)^2 + (c.side.l.y)^2))

## Cryptocyst average side length:
c.side.avg <- rowMeans(cbind(c.side.r, c.side.l))

p.c.side.rl <- ggplot() +
  geom_point(aes(x = c.side.l, y = c.side.l)) +
  ggtitle(paste0("Cryptocyst Side Length, N zooids = ", nrow(images.df), ", N colony = ", length(unique(images.df$colony.id)))) +
  plot.theme +
  scale_y_continuous(name = "Cryptocyst Length Right Side (pixels)") +
  scale_x_continuous(name = "Cryptocyst Length Left Side (pixels)")

ggsave(p.c.side.rl, 
       file = "./Results/cryptocyst.length.w.modern.png", 
       width = 14, height = 10, units = "cm")

summary(lm(c.side.r ~ c.side.l)) 
# slope = 0.87; p-value: < 2.2e-16; r2 = 0.78

## right v left side of operculum
p.oh.rl <- ggplot() +
  geom_point(aes(x = oh.l, y = oh.r)) +
  ggtitle(paste0("Operculum Height, N zooids = ", nrow(images.df), ", N colony = ", length(unique(images.df$colony.id)))) +
  plot.theme +
  scale_y_continuous(name = "Operculum Height Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Height Left Side (pixels)")

ggsave(p.oh.rl, 
       file = "./Results/operculum.height.w.modern.png", 
       width = 14, height = 10, units = "cm")

summary(lm(oh.r ~ oh.l)) 
# slope = 1; p-value < 2.2e-16; r2 = 1; no asymmetry

## centroid size
#It is the square root of the sum of squared distances of all the 
#landmarks you collected in a structure from their centroid.

##### MAKE TABLE & SCALE CORRECT -----
#For 30x it is 0.606 pixels per um(micrometer)

traits.df <- data.frame(boxID = images.df$box_id,
                        image = images.df$image,
                        formation = images.df$formation,
                        colony.id = images.df$colony.id,
                        zooid.id = images.df$zooid.id,
                        zh = zh/.606, #z = zooid; h = height
                        oh = oh/.606, #o = operculum
                        ow.m = ow.m/.606, #w = width; m = mid
                        ow.b = ow.b/.606, #b = base
                        o.side = o.side.avg/.606,
                        mpw.b = mpw.b/.606, #pt = polypide tube
                        cw.m = cw.m/.606, #c = cryptocyst
                        #cw.b = cw.b/.606,
                        cw.d = cw.d/.606, #d = distal
                        c.side = c.side.avg/.606)
                        #zh.zw = zh/cw.d, #similar to LZ/WZ
                        #oh.ow = oh/ow.m,
                        #area = zh*cw.d) # similar to LO/WO 

##### LN TRANSFORM -----
traits.df$ln.zh <- log(traits.df$zh)
traits.df$ln.mpw.b <- log(traits.df$mpw.b)
traits.df$ln.cw.m <- log(traits.df$cw.m)
traits.df$ln.cw.d <- log(traits.df$cw.d)
traits.df$ln.ow.m <- log(traits.df$ow.m)
traits.df$ln.oh <- log(traits.df$oh)
traits.df$ln.o.side <- log(traits.df$o.side)
traits.df$ln.c.side <- log(traits.df$c.side)
#cw.b too variable??

##### MINIMUM 5 ZOOIDS PER COLONY -----
samp.zoo <- traits.df %>%
    dplyr::group_by(colony.id) %>%
    dplyr::summarize(n.zooid = length(unique(zooid.id)),
                     formation = formation[1]) %>%
    as.data.frame()
nrow(samp.zoo) #749 colonies total

too.few <- samp.zoo[samp.zoo$n.zooid < 5,]
nrow(too.few) #161 colonies to remove
#low samples for: Waipuru, Upper Kai-Iwi, Tainui, Modern
#see which ones are removed and if can't redo them
low.samp <- c("Upper Kai-Iwi", "Tainui", "modern")
too.few[too.few$formation %in% low.samp,]
#would add 8 Upper Kai-Iwi; 2 modern; and 8 Tainui
#i.e., not enough still

keep <- samp.zoo$colony.id[samp.zoo$n.zooid >= 5]
length(keep) #588 colonies

df <- traits.df[traits.df$colony.id %in% keep,]
nrow(df) #6178

#### WRITE OUT DATASET ----

write.csv(df,
          "./Results/traits_26Feb2024.csv",
          row.names = FALSE)

####compare
traits.old <- read.csv("./Results/traits_29Sept2023.csv")
traits.new <- read.csv("./Results/traits_8Dec2023.csv")
nrow(traits.new) #6178

traits.new.foss <- traits.new[traits.new$formation != "modern",]
nrow(traits.old) #5971
nrow(traits.new.foss) #5964 FEWER BECAUSE MISSING 1200CC

setdiff(traits.old$imageName, traits.new.foss$image)
#only diff is that Pukenni Limestone is in old version
#Pukenni Limestone, 1200CC, was assigned NKBS
setdiff(traits.old$zooid.id, traits.new.foss$zooid.id) #7, all 1200CC

setdiff(traits.new.foss$image, traits.old$imageName)

traits.old.trim <- traits.old[traits.old$specimenNR != "1200CC",]
setdiff(traits.new.foss$ln.zh, traits.old.trim$ln.zh)
#calculate mean by formation, then look at changes over time
#only the diff with NKBS should be affected, not ALL the other formations...
# and it is in next script (exploratory analysis)



