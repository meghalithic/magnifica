## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) creates trait measurements

## The output of this code is a dataset of all the traits, traits.df, saved
# as traits_date.csv

#### LOAD DATA ----

source("./Scripts/0-env.R")

images.meta <- read.csv("./Data/meta.images_15Nov2023.csv", #meta.images_29Sept2023.csv", 
                        header = TRUE,
                        sep = ",")
#output from outputMetadata.R

df.filter <- read.table("./Data/filteredImages.csv",
                        header = TRUE,
                        sep = ";")

#### EXPLORE DATA ----
nrow(df.filter) #1834
nrow(images.meta) #7056

#### MANIPULATE DATA ####

#extract only the image name, not the entire path; helps match df.filter and images.meta
#merge by images.meta$image
df.filter$fileName.old <- c()
for(i in 1:nrow(df.filter)){
  df.filter$fileName.old[i] <- str_split(df.filter$path.tif[i], "/")[[1]][length(str_split(df.filter$path.tif[i], "/")[[1]])]
}

images.filter.fossil <- images.meta[images.meta$fileName.tif %in% df.filter$fileName.old,]
nrow(images.filter.fossil) #6438
length(unique(images.filter$fileName.tif)) #1395

images.filter.modern <- images.meta[images.meta$formation == "modern",]
nrow(images.filter.modern) #613

images.filter <- as.data.frame(rbind(images.filter.fossil,
                                     images.filter.modern))
nrow(images.filter) #7051

images.filter$formation <- factor(images.filter$formation, 
                                  levels = c("NKLS", "NKBS", "Tewkesbury",
                                             "Waipuru", "Upper Kai-Iwi",
                                             "Tainui", "SHCSBSB", "modern")) 

##### CREATE IDs -----
images.filter$zooid.id <- paste0(images.filter$box_id, "_", images.filter$fileName)
colnames(images.filter)[colnames(images.filter) == 'newSpecimenNR'] <- 'colony.id'

length(unique(images.filter$colony.id[images.filter$formation == "modern"])) #19

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
  ggtitle(paste0("Operculum Side Length, N zooids = ", nrow(images.filter), ", N colony = ", length(unique(images.filter$colony.id)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Operculum Length Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Length Left Side (pixels)")

ggsave(p.o.side.rl, 
       file = "./Results/operculum.length.w.modern.png", width = 14, height = 10, units = "cm")

summary(lm(o.side.r ~ o.side.l)) 
# slope = 0.91873; p-value < 2.2e-16; r2 = 0.95

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
  ggtitle(paste0("Cryptocyst Side Length, N zooids = ", nrow(images.filter), ", N colony = ", length(unique(images.filter$colony.id)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(name = "Cryptocyst Length Right Side (pixels)") +
  scale_x_continuous(name = "Cryptocyst Length Left Side (pixels)")

ggsave(p.c.side.rl, 
       file = "./Results/cryptocyst.length.w.modern.png", 
       width = 14, height = 10, units = "cm")

summary(lm(c.side.r ~ c.side.l)) 
# slope = 0.86534; p-value: < 2.2e-16; r2 = 0.7802

## right v left side of operculum
p.oh.rl <- ggplot() +
  geom_point(aes(x = oh.l, y = oh.r)) +
  ggtitle(paste0("Operculum Height, N zooids = ", nrow(images.filter), ", N colony = ", length(unique(images.filter$colony.id)))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
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

traits.df <- data.frame(boxID = images.filter$box_id,
                        magnification = images.filter$Mag,
                        imageName = images.filter$image,
                        specimenNR = images.filter$specimenNR.tif,
                        formation = images.filter$formation,
                        colony.id = images.filter$colony.id,
                        zooid.id = images.filter$zooid.id,
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
                        oh.ow = oh/ow.m,
                        area = zh*cw.d) # similar to LO/WO 

##### LN TRANSFORM -----
traits.df$ln.zh <- log(traits.df$zh)
traits.df$ln.mpw.b <- log(traits.df$mpw.b)
traits.df$ln.cw.m <- log(traits.df$cw.m)
traits.df$ln.cw.d <- log(traits.df$cw.d)
traits.df$ln.ow.m <- log(traits.df$ow.m)
traits.df$ln.oh <- log(traits.df$oh)
traits.df$ln.o.side <- log(traits.df$o.side)
traits.df$ln.c.side <- log(traits.df$c.side)
traits.df$ln.area <- log(traits.df$area)

##### MINIMUM 5 ZOOIDS PER COLONY -----
samp.zoo <- traits.df %>%
    dplyr::group_by(colony.id) %>%
    dplyr::summarize(n.zooid = length(unique(zooid.id)))
nrow(samp.zoo) #750 colonies total

length(samp.zoo$colony.id[samp.zoo$n.zooid < 5]) #159 colonies to remove
keep <- samp.zoo$colony.id[samp.zoo$n.zooid >= 5]
length(keep) #591 colonies

df <- traits.df[traits.df$colony.id %in% keep,]
nrow(df) #6584

#### WRITE OUT DATASET ----

write.csv(df,
          "./Results/traits_15Nov2023.csv",
          row.names = FALSE)
