# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated:

# The purpose of this script is to run geometric morphometrics analysis on
# Steginoporella magnifica images.
# We test three different shape traits: linear, ratio, landmarks

#### LOAD LIBRARIES ----
library(geomorph)

#### LOAD DATA ----
## datafile of output and metadata merged (from "outputMetadata.R")
df <- read.csv("./Data/meta.images.Jun2023.csv",
               header = TRUE,
               sep = ",",
               stringsAsFactors = FALSE)

#### RECONFIGURE DATA ----

# $formation, factor
# $meta, matrix that contains newFileName, box_id, and image
# $land, landmarks as an array dim 2 (x,y)
df$id <- paste0(df$image, "_", df$box_id)
## formation, factor
formation <- as.factor(df$formation) #7 groups

## meta
meta.image <- df$image
meta.box_id <- df$box_id
meta.newFileName <- df$newFileName
meta.id <- df$id
meta.formation <- df$formation

## landmarks, array dim 2
row.names <- c(0:22)
column.names <- c("x","y")
id <- df$id
land <- array(dim = c(23, 2, nrow(df)),
              dimnames = list(row.names, column.names, id))

for(i in 1:nrow(df)){
  land[1,1,i] <- df$X0[i]
  land[1,2,i] <- df$Y0[i]
  
  land[2,1,i] <- df$X1[i]
  land[2,2,i] <- df$Y1[i]
  
  land[3,1,i] <- df$X2[i]
  land[3,2,i] <- df$Y2[i]
  
  land[4,1,i] <- df$X3[i]
  land[4,2,i] <- df$Y3[i]
  
  land[5,1,i] <- df$X4[i]
  land[5,2,i] <- df$Y4[i]
  
  land[6,1,i] <- df$X5[i]
  land[6,2,i] <- df$Y5[i]
  
  land[7,1,i] <- df$X6[i]
  land[7,2,i] <- df$Y6[i]
  
  land[8,1,i] <- df$X7[i]
  land[8,2,i] <- df$Y7[i]
  
  land[9,1,i] <- df$X8[i]
  land[9,2,i] <- df$Y8[i]
  
  land[10,1,i] <- df$X9[i]
  land[10,2,i] <- df$Y9[i]
  
  land[11,1,i] <- df$X10[i]
  land[11,2,i] <- df$Y10[i]
  
  land[12,1,i] <- df$X11[i]
  land[12,2,i] <- df$Y11[i]
  
  land[13,1,i] <- df$X12[i]
  land[13,2,i] <- df$Y12[i]
  
  land[14,1,i] <- df$X13[i]
  land[14,2,i] <- df$Y13[i]
  
  land[15,1,i] <- df$X14[i]
  land[15,2,i] <- df$Y14[i]
  
  land[16,1,i] <- df$X15[i]
  land[16,2,i] <- df$Y15[i]
  
  land[17,1,i] <- df$X16[i]
  land[17,2,i] <- df$Y16[i]
  
  land[18,1,i] <- df$X17[i]
  land[18,2,i] <- df$Y17[i]
  
  land[19,1,i] <- df$X18[i]
  land[19,2,i] <- df$Y18[i]
  
  land[20,1,i] <- df$X19[i]
  land[20,2,i] <- df$Y19[i]
  
  land[21,1,i] <- df$X20[i]
  land[21,2,i] <- df$Y20[i]
  
  land[22,1,i] <- df$X21[i]
  land[22,2,i] <- df$Y21[i]
  
  land[23,1,i] <- df$X22[i]
  land[23,2,i] <- df$Y22[i]
}



## put together into a list
bryo <- list(formation = formation, 
             meta.image = meta.image, 
             meta.box_id = meta.box_id, 
             meta.newFileName = meta.newFileName, 
             meta.id = meta.id, 
             meta.formation = meta.formation, 
             land = land)

dim(bryo$land)
dim(two.d.array(bryo$land))

#### PLOT ----

plotAllSpecimens(bryo$land, mean=FALSE)

#### PROCRUSTES ANALYSIS ----
## pretty sure this corrects alignment...
bryo.gpa <- gpagen(bryo$land, print.progress = FALSE)
summary(bryo.gpa)

bryo.gpa$coords # a 3D array of Procrustes coordinates 
bryo.gpa$Csize # a vector of centroid sizes

#### WHAT ARE LINKS?? ----
#plotAllSpecimens(plethodon.gpa$coords,links=plethodon$links)


