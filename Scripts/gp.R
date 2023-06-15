# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated:

# The purpose of this script is to create a P and G matrix
# These matrices are initially made based on the following linear measurements:
# zh = zooid height
# mpw.b = median process base width
# cw.m = cryptocyst mid-width
# cw.d = cryptocyst distal-width
# ow.m = operculum mid-width

#### LOAD LIBRARIES ----
library(ggplot2)
library(data.table)
library(evolqg)
library(dplyr)
library(nse)
library(grid)
library(gridBase)
library(gridExtra)
library(MCMCglmm)
library(corrplot)
library(MASS)
library(reshape2)
library(scatterplot3d)
library(rgl)
library(scales)
library(RColorBrewer)
library(coin)

#### LOAD DATA ----

df <- read.csv("./Results/traits.csv",
               header = TRUE, 
               sep = ",",
               stringsAsFactors = FALSE)
df$new.id <- paste0(df$boxID, "_", df$image)

#### MANIPULATE DATA ----

# Extract unique elements and trait names

colony_list <- unique(df$image)
length(colony_list) #1824

levels(df$formation)<-c("NKLS", "NKBS", "Tewkesbury", "Waipuru", "Upper Kai-Iwi", "SHCSBSB", "Tainui") #old to young
formation_list<-unique(df$formation)
length(formation_list)

traits = names(df[, c("zh", "mpw.b", "cw.m", "cw.d", "ow.m")])
df$zh <- log10(df$zh)
df$mpw.b <- log10(df$mpw.b)
df$cw.m <- log10(df$cw.m)
df$cw.d <- log10(df$cw.d)
df$ow.m <- log10(df$ow.m)

colNums <- match(c(traits,"image"),names(df))

df = as.data.frame(df)

#### PLOT TRAITS ----
Fig <- list ()
for (i in 1:length(traits)){
  Fig[[i]] = ggplot(data = df)+ 
    geom_density(aes(x = log10(df[,5+i]), 
                     group = formation,
                     col = formation)) + 
    theme(legend.position="none", text = element_text(size=20)) +
    scale_x_continuous(name = traits[i])
  
}

ml <- marrangeGrob(Fig, nrow = 5, ncol = 1)
ml

ggsave(ml, file = "./Results/trait.interest_distribution.png", width = 14, height = 10, units = "cm")


## most would be normal without small hump...

#### REDUCE TO TRAITS OF INTEREST ----
trt_lg_N = c("formation","specimenNR", traits)
dat_lg_N = df[intersect(colnames(df), trt_lg_N)]
head(dat_lg_N)

#### SUMMARY STATISTICS ----
mean_by_formation = dat_lg_N %>%
  group_by(formation) %>%
  summarize(avg.zh = mean(zh, na.rm = T),
            avg.mpw.b = mean(mpw.b, na.rm = T),
            avg.cw.m = mean(cw.m, na.rm = T),
            avg.cw.d = mean(cw.d, na.rm = T),
            avg.ow.m = mean(ow.m, na.rm = T))

mean_by_colony = dat_lg_N %>%
  group_by(specimenNR) %>%
  summarize(avg.zh = mean(zh, na.rm = T),
            avg.mpw.b = mean(mpw.b, na.rm = T),
            avg.cw.m = mean(cw.m, na.rm = T),
            avg.cw.d = mean(cw.d, na.rm = T),
            avg.ow.m = mean(ow.m, na.rm = T))  

means = dat_lg_N %>%
  summarize(avg.zh = mean(zh, na.rm = T),
            avg.mpw.b = mean(mpw.b, na.rm = T),
            avg.cw.m = mean(cw.m, na.rm = T),
            avg.cw.d = mean(cw.d, na.rm = T),
            avg.ow.m = mean(ow.m, na.rm = T))  


#### SCALE DATA ----
dat_lda=dat_lg_N
dat_lda[,3:7]=scale(dat_lda[,3:7],center=F,scale=T) #just traits
data_discriminant=dat_lda

data_plot=na.omit(data_discriminant[,1:7]) #all rows

r3 <- lda(formula = formation ~ ., 
          data = data_plot[,2:7], method='mle') #just traits + formation 

plda = predict(object = r3, # predictions
               newdata = data_plot)
prop.lda = r3$svd^2/sum(r3$svd^2) #proportion of the variance explained by each LD axis

dataset = data.frame(formation = (data_plot)[,"formation"],
                     lda = plda$x)

p1 <- ggplot(dataset) + 
  geom_point(aes(lda.LD1, lda.LD2, color = formation), 
             size = 1, alpha = .75) + 
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
p1
ggsave(p1, file = "./Results/trait_discriminant.png", width = 14, height = 10, units = "cm")

#### P MATRIX ----

#### G MATRIX ----
#Preparing the data
colony_form=split.data.frame(dat_lg_N,dat_lg_N$formation) #check number of colonies not zooids
sample_sizes=lapply(colony_form,function(x){dim(x)[1]})
enough=as.matrix(sample_sizes)>100
form_final=colony_form[names(enough[enough==T,])]
form_data=lapply(form_final,function(x) x[complete.cases(x),])

##### PRIORS -----
phen.var=lapply(form_data,function (x){ (cov(x[,3:7]))}) #traits
prior=lapply(phen.var, function (x){list(G=list(V=x/2,nu=2),
                                         R=list(V=x/4,nu=2))})
#rescale by 100 or 1000

##### MCMC -----
#Running the MCMC chain
model_G=list()
for (i in 1:length(form_data)){ #length 7 because 7 formations
  model_G[[i]]<-MCMCglmm(cbind(zh, ow.m, mpw.b, cw.m, cw.d)~trait-1,
                         random = ~us(trait):specimenNR, #the number of these determines # of Gs
                         rcov=~us(trait):units,
                         family=rep("gaussian",5), #5 for 5 traits
                         data=form_data[[i]],
                         nitt=1500000,thin=1000,burnin=500000,
                         prior=prior[[i]],verbose=TRUE)
  
}


