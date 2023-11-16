## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) make list of images for comparison
# 2a) compare measurements to those taken using imageJ
# 2b) compare measurements to those extracted from landmarks using imageJ

## The output of this code is:
# - csv of images to import into imageJ
# - comparison of measurements to imageJ both of linear measurements and from landmarks

#### LOAD DATA ----

source("./Scripts/0-env.R")

df <- read.csv("./Results/colonies.traits_15Nov2023.csv",
               header = TRUE)

#### RANDOMLY SELECT IMAGES ----
## split by formation
## 3 zooids per colony
## 3 colonies per formation

# just need 3 random colonies
col.samp <- c()
for(i in 1:length(unique(df$formation))){
    x <- sample(unique(df$colony.id[df$formation == unique(df$formation)[i]]), 3, replace = FALSE)
    col.samp <- c(col.samp, x)
}
length(col.samp) #24; 3*8

#trim dataset
df.samp <- df[df$colony.id %in% col.samp,]

#extract images and select those that know have 3 at least zooids
n.img <- df.samp %>%
    dplyr::group_by(imageName) %>%
    dplyr::summarise(n = n()) %>%
    as.data.frame()
too.few <- n.img$imageName[n.img$n < 3]
df.samp.trim <- df.samp[!(df.samp$imageName %in% too.few),]

#randomly select 3 box_ids from each image
# zoo.samp <- c()
# for(i in 1:length(unique(df.samp.trim$imageName))){
#     x <- sample(df$boxID[df$imageName == unique(df.samp.trim$imageName)[i]], 3, replace = FALSE)
#     zoo.samp <- c(zoo.samp, x)
# }
# length(zoo.samp) #120; 3*3*8 at least
# 
# df.samp.id <- df.samp.trim[df.samp.trim$boxID %in% zoo.samp,]

write.csv(df.samp.trim,
          "./Data/imageJ_sensitivity/image.subsample.csv",
          row.names = FALSE)
#will have to manually compare to output from Steginator-magnifica by matching top (y) and left (x) pixel coordinates
#will manually select which 3 zooids for each formation

#### COMPARE IMAGES ----

## LOAD IMAGEJ RESULTS

##### COMPARE LINEAR MEASUREMENTS -----

##### COMPARE MEASUREMENTS FROM LANDMARKS -----