## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
## 1) extracts file names from the zip file of Steginoporella images
## 2) creates a csv file with the information parsed
## 3) compares the files to the bryozoan metadata file

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)
require(splitstackshape)
require(data.table)

#### EXTRACT FILE NAMES ----
list = list.files(path = "/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images",
                  full.names = TRUE,
                  recursive = TRUE)
## remove abstract book
list.rm <- list[!grepl("/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images/Stegs/ISME 14 ABSTRACT BOOK",
                       list)]

#### TXT ONLY ----

list.trim <- list.rm[!grepl("*.tif",
                            list.rm)]

length(list.trim) #1889

##unlist
txtPath <- unlist(list.trim)

#### READ TXT FILES ----

##practice with one file
f <- read.table("/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images/Sara/002_CV_1_15v_x30.txt",
                sep = "^",
                fileEncoding="UTF-16",
                skip = 1)

## now make two columns, using "=" as deliminator

ff <- cSplit(f, 'V1',
             sep="=",
             stripWhite = TRUE,
             type.convert = FALSE)

#seems Condition is multiple "="
condition <- str_split(ff[ff$V1_1 == "Condition",],
                       pattern = "\ ")

av <- c("AV", gsub(".0kV", "v",condition[[3]][1]))
mag <- c(condition[[3]][2], condition[[4]][1])
wd <- c(condition[[4]][2], condition[[5]][1])
lensMode <- c(condition[[5]][2], condition[[6]][1])
path <- c("path", "/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images/Sara/002_CV_1_15v_x30.txt")

cond.paste <- paste(ff$V1_2[ff$V1_1 == "Condition"], 
                    ff$V1_3[ff$V1_1 == "Condition"],
                    ff$V1_4[ff$V1_1 == "Condition"], 
                    ff$V1_5[ff$V1_1 == "Condition"],
                    ff$V1_6[ff$V1_1 == "Condition"], 
                    sep = " ")

ff2 <- ff
ff2$V1_2[ff2$V1_1 == "Condition"] <- cond.paste

ff3 <- ff2[,1:2]

ff4 <- rbind(path, as.data.frame(ff3), av, mag, wd, lensMode)

ff4$V1_1

## now for all!

files.df <- data.frame()

for(i in 1:length(txtPath)){
  f <- read.table(txtPath[i],
                  sep = "^",
                  fileEncoding = "UTF-16",
                  skip = 1)
  
  ## now make two columns, using "=" as deliminator
  
  ff <- cSplit(f, 'V1',
               sep = "=",
               stripWhite = TRUE,
               type.convert = FALSE)
  
  #seems Condition is multiple "="
  condition <- str_split(ff[ff$V1_1 == "Condition",],
                         pattern = "\ ")
  
  av <- c(condition[[2]][1],condition[[3]][1])
  mag <- c(condition[[3]][2], condition[[4]][1])
  wd <- c(condition[[4]][2], condition[[5]][1])
  lensMode <- c(condition[[5]][2], condition[[6]][1])
  path <- c("path", txtPath[i])
   
  cond.paste <- paste(ff$V1_2[ff$V1_1 == "Condition"], 
                      ff$V1_3[ff$V1_1 == "Condition"],
                      ff$V1_4[ff$V1_1 == "Condition"], 
                      ff$V1_5[ff$V1_1 == "Condition"],
                      ff$V1_6[ff$V1_1 == "Condition"], 
                      sep = " ")
  
  ff2 <- ff
  
  ff2$V1_2[ff2$V1_1 == "Condition"] <- cond.paste
  
  ff3 <- ff2[,1:2]
  
  ff4 <- rbind(path, as.data.frame(ff3), av, mag, wd, lensMode)

  names <- ff4$V1_1
  ff5 <- as.data.frame(t(ff4[,-1]))
  colnames(ff5) <- names
  
  files.df <- rbind(files.df, ff5)
  
}

nrow(files.df) #1889

files.df$fileName <- basename(files.df$path)
files.df$image <- str_extract(files.df$fileName, pattern = "[^.]+")


imageName.parse <- str_split(files.df$image,
                             pattern = "_")
files.df$specimenNR <- ""
files.df$number <- ""
files.df$colonyCurve <- ""
files.df$pictureNumber <- ""
files.df$AV <- ""
files.df$magnification <- ""
files.df$backscatter <- ""
for(i in 1:nrow(files.df)){
  files.df$number[i] <- imageName.parse[[i]][1]
  files.df$colonyCurve[i] <- imageName.parse[[i]][2]
  files.df$pictureNumber[i] <- imageName.parse[[i]][3]
  files.df$AV[i] <- imageName.parse[[i]][4]
  files.df$magnification[i] <- str_extract(imageName.parse[[i]][5],
                                  pattern = "[^.]+")
  files.df$backscatter[i] <- str_extract(imageName.parse[[i]][6],
                                pattern = "[^.]+")
  files.df$specimenNR[i] <- paste0(imageName.parse[[i]][1], imageName.parse[[i]][2])
}

#write.csv(files.df,
#          "./Data/txt_metadata.csv",
#          row.names = FALSE)
