## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code adds the output from txtFiles to the output from imageFiles

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)

#### LOAD DATA ----

df.list <- read.table("./Data/imageList.csv",
                      header = TRUE,
                      sep = ";")

files.df <- read.table("./Data/txt_metadata.csv",
                      header = TRUE,
                      sep = ";")

#### COMBINE DATA ----

##### CREATE SHARED FILE NAME -----

df.list$image[1]

files.df$image[1]

nrow(files.df)
nrow(df.list)

df.images <- df.list[df.list$ext == "tif",] 
nrow(df.images) #there's an extra image...

df.image.meta <- merge(df.images, files.df,
                       by = "image",
                       all.x = TRUE, all.y = TRUE)
nrow(df.image.meta) #1940...some mismatches...

#write.csv(df.image.meta,
#          "./Data/imageFilesMetadata.csv",
#          row.names = FALSE)

## LOOK AT IT MANUALLY ##
