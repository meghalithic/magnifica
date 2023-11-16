#### LOAD PACKAGES ----
library(coin)
library(corrplot)
library(data.table)
require(dplyr)
library(evolqg)
library(evolvability)
require(ggplot2)
library(graphics)
library(grid)
library(gridBase)
library(gridExtra)
require(lmodel2)
library(magrittr)
library(MASS)
library(matrixcalc)
library(MCMCglmm)
library(nse)
library(RColorBrewer)
require(reshape2)
library(rgl)
library(scales)
library(scatterplot3d)
require(stringr)
library(tibble)
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

col.form = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
             "#00A9FF", "#C77CFF", "#FF61CC", "#00BFC4")


col.traits = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
               "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
