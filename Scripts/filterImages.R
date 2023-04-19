## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) reduces image list to those with 30 magnification

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)
require(ggplot2)
require(reshape2)
require(lmodel2)
require(tidyverse)

#### LOAD DATA ----
bryo.images <- read.table("./Data/Steginoporella_magnifica_image_metadata_17Apr2023.csv",
                          header = TRUE,
                          sep = ";")
nrow(bryo.images) #1890

bryo.images.30 <- bryo.images[bryo.images$Magnification == 30,]
nrow(bryo.images.30) #1835

#write.csv(bryo.images.30,
#          "./Data/filteredImages.csv",
#          row.names = FALSE)
