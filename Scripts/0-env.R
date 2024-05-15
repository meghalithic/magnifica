#### LOAD PACKAGES ----
#library(coin)
library(corrplot)
#library(data.table)
#require(dplyr)
library(evolqg) 
library(evolvability)
require(ggplot2)
#library(graphics)
library(grid)
library(gridBase)
library(gridExtra)
#require(lmodel2)
#library(magrittr)
library(MASS)
library(matrixcalc)
library(MCMCglmm)
#library(nse)
require(reshape2)
#library(rgl)
#library(scales)
#library(scatterplot3d)
require(stringr)
#library(tibble)
require(tidyverse)

#### FUNCTIONS ----
source("./Scripts/norm.vector.funct.R")
source("./Scripts/magnitudeFunction.R")


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

col.traits.repo <- c("#F8766D", #zh
                     "#CD9600", #oh
                     "#B79F00", #ow.m
                     "#00BE67", #o.side
                     "#00C094", #mpw.b
                     "#619CFF", #cw.d
                     "#00BFC4", #cw.m
                     "#C77CFF") #c.side)

plot.theme <- theme(text = element_text(size = 16),
      legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      plot.background = element_rect(fill = 'transparent', color = NA))

#pch code
# P = triangle = 17
# G = square = 15
# Global g = diamond = 18
# Temp = circle = 16

#### FORMATIONS ----
formation_transition <- c("NKLS to NKBS", 
                          "NKBS to Tewkesbury",
                          "Tewkesbury to Upper Kai-Iwi",
                          "Upper Kai-Iwi to Tainui", 
                          "Tainui to SHCSBSB",
                          "SHCSBSB to modern")
formation_transition <- factor(formation_transition,
                               levels = c("NKLS to NKBS", 
                                          "NKBS to Tewkesbury",
                                          "Tewkesbury to Upper Kai-Iwi",
                                          "Upper Kai-Iwi to Tainui", 
                                          "Tainui to SHCSBSB",
                                          "SHCSBSB to modern"))

formation_list <- c("NKLS", "NKBS", "Tewkesbury", 
                    "Upper Kai-Iwi", "Tainui", "SHCSBSB", "modern")
formation_list <- factor(formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                               "Upper Kai-Iwi", "Tainui", "SHCSBSB", "modern"))

#### TEMPERATURE ----
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

#convert ot ˚C
form.meta$temp <- 16.5 - (4.3*form.meta$med.O18) + (0.14*(form.meta$med.O18^2))

#make factors
form.meta$formationCode <- factor(form.meta$formationCode, 
                                  levels = c("NKLS", "NKBS", "Tewkesbury",
                                             "Upper Kai-Iwi", "Tainui",
                                             "SHCSBSB", "modern"))


#write.csv(form.meta,
#          "./Results/formation.with.isotopes.csv",
#          row.names = FALSE)
