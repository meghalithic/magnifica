#### LOAD PACKAGES ----
#library(coin)
library(corrplot)
#library(data.table)
#require(dplyr)
library(evolqg) 
#library(evolvability)
require(ggplot2)
#library(graphics)
#library(grid)
#library(gridBase)
#library(gridExtra)
#require(lmodel2)
#library(magrittr)
library(MASS)
#library(matrixcalc)
library(MCMCglmm)
#library(nse)
require(reshape2)
#library(rgl)
#library(scales)
#library(scatterplot3d)
require(stringr)
#library(tibble)
require(tidyverse)

#### PLOT THEME ----
#formations and colors: 
#NKLS = #F8766D
#NKBS = #CD9600
#Twekesbury = #7CAE00
#Waipuru = #00BE67
#Upper Kai-Iwi = #00A9FF
#Tainui = #C77CFF
#SHCSBSB = #FF61CC

col.form = c("#F8766D", "#CD9600",  "#00BE67", #"#7CAE00"
             "#00A9FF", "#C77CFF", "#FF61CC", "#00BFC4") #YES


col.traits = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
               "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")

form.meta <- read.csv("~/Documents/GitHub/bryozoa/stegino_metadata/newMetadata/formations.csv", header = TRUE)

#this is downloaded from: http://www.lorraine-lisiecki.com/LR04_MISboundaries.txt
oxy.18 <- read.csv("Data/∂18O.csv",
                   header = TRUE)

bottom = as.numeric(form.meta$Isotope_Stage_Start)
top = as.numeric(form.meta$Isotope_Stage_End)
form.meta$med.O18 <- c()
form.meta$sd.med.O18 <- c()
form.meta$n.O18 <- c()
for (i in 1:nrow(form.meta)){
    temp = oxy.18$d18O[which(oxy.18$Time <= bottom[i] & oxy.18$Time >= top[i])]
    form.meta$med.O18[i] = median(temp)
    form.meta$sd.med.O18[i] = sd(temp)
    form.meta$n.O18[i] <- length(temp)
}
