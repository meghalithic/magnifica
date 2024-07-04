## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) creates trait measurements

## The output of this code is a dataset of all the traits, traits.df, saved
# as traits_date.csv

#### LOAD DATA ----

source("./Scripts/0-env.R")

images.df <- read.csv("./Data/images.filtered_1Jul2024.csv", #images.merged_30Nov2023.csv,
                      header = TRUE, 
                      sep = ",")

#### EXPLORE DATA ----
nrow(images.df) #7412

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
images.df$zh.px <- len.py(images.df$X4, images.df$X12,
                       images.df$Y4, images.df$Y12)

## Operculum height of left side: 4 to 20
images.df$oh.l.px <- len.py(images.df$X4, images.df$X20,
                         images.df$Y4, images.df$Y20)

## Operculum height of right side: 4 to 21
images.df$oh.r.px <- len.py(images.df$X4, images.df$X21,
                         images.df$Y4, images.df$Y21)

## Operculum mid-width (maximum width at centerline): 19 to 0
images.df$ow.m.px <- len.py(images.df$X19, images.df$X0,
                         images.df$Y19, images.df$Y0)

## Operculum base width: 21 to 20
images.df$ow.b.px <- len.py(images.df$X21, images.df$X20,
                         images.df$Y21, images.df$Y20)

## Operculum side length of right side: 21 to 18
images.df$o.side.r.px <- len.py(images.df$X21, images.df$X18,
                             images.df$Y21, images.df$Y18)

## Operculum side length of left side: 20 to 17
images.df$o.side.l.px <- len.py(images.df$X20, images.df$X17,
                             images.df$Y20, images.df$Y17)

## Median process base width: 5 to 6
images.df$mpw.b.px <- len.py(images.df$X5, images.df$X6,
                          images.df$Y5, images.df$Y6)

## Cryptocyst mid-width: 10 to 11
images.df$cw.m.px <- len.py(images.df$X10, images.df$X11,
                         images.df$Y10, images.df$Y11)

## Cryptocyst base width: 9 to 1
#images.df$cw.b.px <- len.py(images.df$X9, images.df$X1,
#                         images.df$Y9, images.df$Y1)

## Cryptocyst distal width: 8 to 7
images.df$cw.d.px <- len.py(images.df$X8, images.df$X7,
                         images.df$Y8, images.df$Y7)

## Cryptocyst side length of right side: 1 to 7
images.df$c.side.r.px <- len.py(images.df$X1, images.df$X7,
                             images.df$Y1, images.df$Y7)

## Cryptocyst side length of left side: 9 to 8
images.df$c.side.l.px <- len.py(images.df$X9, images.df$X8,
                             images.df$Y9, images.df$Y8)

## centroid size
#It is the square root of the sum of squared distances of all the 
#landmarks you collected in a structure from their centroid.

##### MAKE TABLE & SCALE CORRECT -----
#For 30x it is 0.606 pixels per um(micrometer)
#all bleed are at 40x, which is 0.825 pixels per um
#the only bleed numbers included are: 686, 244, 242; form.no = BLEED
#recalculate them manually
images.df$zh <- images.df$zh.px/images.df$scale
images.df$oh.l <- images.df$oh.l.px/images.df$scale
images.df$oh.r <- images.df$oh.r.px/images.df$scale
images.df$ow.b <- images.df$ow.b.px/images.df$scale
images.df$ow.m <- images.df$ow.m.px/images.df$scale
images.df$mpw.b <- images.df$mpw.b.px/images.df$scale
images.df$cw.m <- images.df$cw.m.px/images.df$scale
images.df$cw.d <- images.df$cw.d.px/images.df$scale
images.df$c.side.r <- images.df$c.side.r.px/images.df$scale
images.df$c.side.l <- images.df$c.side.l.px/images.df$scale
images.df$o.side.r <- images.df$o.side.r.px/images.df$scale
images.df$o.side.l <- images.df$o.side.l.px/images.df$scale

##### CHECK FOR ASYMMETRY -----
p.o.side.rl <- ggplot(images.df) +
  geom_point(aes(x = o.side.l, y = o.side.r)) + #,
                 #group = formation, fill = formation, col = formation)) +
  ggtitle(paste0("Operculum Side Length, N zooids = ", nrow(images.df), ", N colony = ", length(unique(images.df$colony.id)))) +
  plot.theme +
  theme(legend.position = "top") + 
  scale_y_continuous(name = "Operculum Length Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Length Left Side (pixels)")
#modern doesn't look that different from everything else

ggsave(p.o.side.rl, 
       file = "./Results/operculum.length.png", width = 14, height = 10, units = "cm")

summary(lm(images.df$o.side.r ~ images.df$o.side.l)) 
# slope = 0.92; p-value < 2.2e-16; r2 = 0.95

p.c.side.rl <- ggplot(images.df) +
  geom_point(aes(x = c.side.l, y = c.side.l)) + #,
    #group = formation, fill = formation, col = formation)) +
  ggtitle(paste0("Cryptocyst Side Length, N zooids = ", nrow(images.df), ", N colony = ", length(unique(images.df$colony.id)))) +
  plot.theme +
  scale_y_continuous(name = "Cryptocyst Length Right Side (pixels)") +
  scale_x_continuous(name = "Cryptocyst Length Left Side (pixels)")

ggsave(p.c.side.rl, 
       file = "./Results/cryptocyst.length.png", 
       width = 14, height = 10, units = "cm")

summary(lm(images.df$c.side.r ~ images.df$c.side.l)) 
# slope = 0.86; p-value: < 2.2e-16; r2 = 0.78

## right v left side of operculum height
p.oh.rl <- ggplot(images.df) +
  geom_point(aes(x = oh.l, y = oh.r)) + #,
    #group = formation, fill = formation, col = formation)) +
  ggtitle(paste0("Operculum Height, N zooids = ", nrow(images.df), ", N colony = ", length(unique(images.df$colony.id)))) +
  plot.theme +
  scale_y_continuous(name = "Operculum Height Right Side (pixels)") +
  scale_x_continuous(name = "Operculum Height Left Side (pixels)")

ggsave(p.oh.rl, 
       file = "./Results/operculum.height.png", 
       width = 14, height = 10, units = "cm")

summary(lm(images.df$oh.r ~ images.df$oh.l)) 
# slope = 0.92; p-value < 2.2e-16; r2 = .89; no asymmetry

## Operculum average side length:
images.df$o.side <- rowMeans(cbind(images.df$o.side.l, images.df$o.side.r))

## Cryptocyst average side length:
images.df$c.side <- rowMeans(cbind(images.df$c.side.l, images.df$c.side.r))

## Operculum height
images.df$oh <- (.5/images.df$ow.b)*sqrt(images.df$ow.b+images.df$oh.r+images.df$oh.l)

##### TRIM TO TRAITS ONLY ----
traits.df <- images.df %>%
  dplyr::select(box_id, image, 
                source, Mag, 
                colony.id, zooid.id, 
                formation, locality,
                zh, oh, ow.m, ow.b, 
                mpw.b, cw.m, cw.d,
                o.side, c.side)

colnames(traits.df)[colnames(traits.df) == 'box_id'] <- 'boxID'

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

##### MINIMUM 3 ZOOIDS PER COLONY -----
samp.zoo <- traits.df %>%
    dplyr::group_by(colony.id) %>%
    dplyr::summarize(n.zooid = length(unique(zooid.id)),
                     formation = formation[1]) %>%
    as.data.frame()
nrow(samp.zoo) #822 colonies total

too.few <- samp.zoo[samp.zoo$n.zooid < 3,]
nrow(too.few) #60 colonies to remove
table(too.few$formation)
#see which ones are removed and if can't redo them

keep <- samp.zoo$colony.id[samp.zoo$n.zooid >= 3]
length(keep) #760 colonies

df <- traits.df[traits.df$colony.id %in% keep,]
nrow(df) #7309

df %>%
    group_by(formation) %>%
    summarise(n.col = length(unique(colony.id)))

df %>%
    group_by(formation) %>%
    summarise(n.image = length(unique(image)))
length(unique(df$image)) #1578

#### CHECK FOR EFFECT OF MAGNIFICATION ----
colnames(df)

table(df$Mag, df$formation) #modern is the only one with different magnifications
#the most are in 50

ggplot(df) +
    geom_density(aes(x = log(zh),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#shifted larger with larger mag

ggplot(df) +
    geom_density(aes(x = log(ow.b),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#slightly shifted larger with larger mag

ggplot(df) +
    geom_density(aes(x = log(ow.m),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#shifted larger with larger mag

ggplot(df) +
    geom_density(aes(x = log(mpw.b),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#shifted larger with larger mag

ggplot(df) +
    geom_density(aes(x = log(cw.m),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#amazingly not shifted

ggplot(df) +
    geom_density(aes(x = log(cw.d),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#amazingly not shifted

ggplot(df) +
    geom_density(aes(x = log(o.side),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#only slightly shifted

ggplot(df) +
    geom_density(aes(x = log(c.side),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#amazingly not shifted

ggplot(df) +
    geom_density(aes(x = log(oh),
                     group = Mag, col = Mag)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#mostly not shifted

summary(manova(cbind(log(zh), log(ow.b), log(ow.m), log(mpw.b), log(cw.m), log(cw.d), log(o.side), log(c.side), log(oh)) ~ Mag, 
               data = df[df$formation == "modern",]))
#significant

#trim to just 30
df.30 <- df[df$Mag == "x30",]
df.30 %>% 
    group_by(formation) %>%
    summarise(n.col = length(unique(colony.id))) #still 30 in modern! woo!

df <- df.30

df.30.mod <- df[df$formation == "modern",]

ggplot(df.30.mod) +
    geom_density(aes(x = log(zh),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#TS-3 weird, but otherwise ok

ggplot(df.30.mod) +
    geom_density(aes(x = log(ow.b),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#weirdness with M797, TAN1108/233

ggplot(df.30.mod) +
    geom_density(aes(x = log(ow.m),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#

ggplot(df.30.mod) +
    geom_density(aes(x = log(mpw.b),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#something weird with Z8662

ggplot(df.30.mod) +
    geom_density(aes(x = log(cw.m),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#some weirdness with M797 and TS-10

ggplot(df.30.mod) +
    geom_density(aes(x = log(cw.d),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#some weirdness with TAN1108/233, M797

ggplot(df.30.mod) +
    geom_density(aes(x = log(o.side),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#Z8662 off

ggplot(df.30.mod) +
    geom_density(aes(x = log(c.side),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#weirdness again with M797, a bit with KWB_Feb

ggplot(df.30.mod) +
    geom_density(aes(x = log(oh),
                     group = locality, col = locality)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    theme(legend.position = "right") + 
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))
#weirdness with TS-3

summary(manova(cbind(log(zh), log(ow.b), log(ow.m), log(mpw.b), log(cw.m), log(cw.d), log(o.side), log(c.side), log(oh)) ~ locality, 
               data = df.30.mod))
#significant

#what if remove z8862 and M797 and TAN1108/233?
rm.loc <- c("M797", "TAN1108/233", "Z8662")
summary(manova(cbind(log(zh), log(ow.b), log(ow.m), log(mpw.b), log(cw.m), log(cw.d), log(o.side), log(c.side), log(oh)) ~ locality, 
               data = df.30.mod[!(df.30.mod$locality %in% rm.loc),])) #still sig
#what if remove all niwa?
summary(manova(cbind(log(zh), log(ow.b), log(ow.m), log(mpw.b), log(cw.m), log(cw.d), log(o.side), log(c.side), log(oh)) ~ locality, 
               data = df.30.mod[df.30.mod$source != "NIWA",])) #still sig
#still sig

#### WRITE OUT DATASET ----

write.csv(df,
          "./Results/traits_1Jul2024.csv",
          row.names = FALSE)

#### OLD ----
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



