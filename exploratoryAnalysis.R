## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)

#### LOAD DATA ----
output <- read.csv("output.csv", header = TRUE)

#### EXPLORE DATA ----

nrow(output) #16924

colnames(output)
#id
#box_id seems to be the coordinates for the box WHAT ORDER??
#unsure what box_left is
#unsure what box_top is
#series of coordinates and landmark number [e.g., (X0, Y0)]

## WHAT ARE THE LANDMARKS??

##### RENAME ID -----
id.parse <- str_split(output$id,
                      pattern = "/")

id.only <- unlist(lapply(id.parse, function (x) x[9]))


names(output)[names(output)=="id"] <- "path"
output$id <- id.only
