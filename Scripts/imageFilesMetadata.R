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

df.list$fileName[1]

files.df$fileName[1] #for txt
files.df$ImageName[1] #for tif; ImageName is inside of txt metadata file

nrow(files.df) #1889
nrow(df.list) #3779

## make two: one for txt one for tif

df.images <- df.list[df.list$ext == "tif",] #1890
df.txt <- df.list[df.list$ext == "txt",] #1889

length(setdiff(df.images$fileName, files.df$ImageName)) #49
length(setdiff(files.df$ImageName, df.images$fileName)) #39

df.image.meta <- merge(df.images, files.df,
                      by.x = "fileName", by.y = "ImageName",
                      all.x = TRUE, all.y = TRUE) #1938

setdiff(df.txt$fileName, files.df$fileName)
setdiff(files.df$fileName, df.txt$fileName)
# no difference in txt files


#write.csv(df.image.meta,
#          "./Data/image_txt.csv",
#          row.names = FALSE)

## RECONCILE MANUALLY ##
