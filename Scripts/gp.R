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
prior=lapply(phen.var, function (x){list(G=list(G1=list(V=x/2,nu=2)),
                                         R=list(V=x/4,nu=2))})

##### MCMC -----
#Running the MCMC chain
model_G=list()
for (i in 1:length(form_data)){ #length 7 because 7 formations
  model_G[[i]]<-MCMCglmm(cbind(zh, ow.m, mpw.b, cw.m, cw.d)~trait-1,
                         random = ~us(trait):specimenNR, #the number of these determines # of Gs #+ us(trait):formation
                         rcov=~us(trait):units,
                         family=rep("gaussian",5), #5 for 5 traits
                         data=form_data[[i]],
                         nitt=1500000,thin=1000,burnin=500000,
                         prior=prior[[i]],verbose=TRUE)
  
}

save(model_G, file="./Results/New_g_matrices.RData")

load(file = "New_g_matrices.RData") #load the g matrices calculated above 

summary(model_G[[1]])
summary(model_G[[2]])
summary(model_G[[3]])
summary(model_G[[4]])
summary(model_G[[5]])
summary(model_G[[6]])
summary(model_G[[7]])

##plots to see where sampling from:
plot(model_G[[1]]$VCV) #catepillar!
plot(model_G[[2]]$VCV) #catepillar!
plot(model_G[[3]]$VCV) #catepillar!
plot(model_G[[4]]$VCV) #catepillar!
plot(model_G[[5]]$VCV) #catepillar!
plot(model_G[[6]]$VCV) #catepillar!
plot(model_G[[7]]$VCV) #catepillar!
#formations: "NKLS", "NKBS", "Tewkesbury", "Waipuru", "Upper Kai-Iwi", "SHCSBSB", "Tainui" #old to young
## ORDER HERE: NKBS, NKLS, SSHCSBSB, Tainui, Twekesbury, Upper Kai-Iwi, Waipuru (despite releveling :/)

#Retrieving G from posterior
model = model_G
data = (dat_lg_N)
ntraits = 5
Gmat = lapply(model, function (x) { 
  matrix(posterior.mode(x$VCV)[1:ntraits^2], ntraits, ntraits)})


#Standardizing G by the trait means 

mean_by_form = setDT(na.omit(data[,2:7]))[, lapply(.SD, mean, na.rm = F), by = .(formation)] #omit col 1
u_form=split(mean_by_form,mean_by_form$formation)
test=lapply(u_form, function (x){ data.matrix(x[,2:6])}) #omit col and form
test_std=lapply(test,function (x){(as.numeric(x))%*%t(as.numeric(x))})
G_std=list()
for (i in 1:length(Gmat)){
  G_std[[i]]=Gmat[[i]]/(test_std[names(form_data[i])][[1]])
}
G_std
names(G_std)=names(form_data[1:i])

##Genetic variance in traits and eigenvectors

#load(file="New_g_matrices.RData") #load the g matrices calculated above 

lapply(G_std, isSymmetric)  
std_variances = lapply(G_std, diag)
paste("Trait variances")
head(std_variances)

eig_variances=lapply(G_std, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(eig_variances)

eig_percent=lapply(eig_variances, function (x) {x/sum(x)})
eig_per_mat=do.call(rbind, eig_percent)
eig_per_mat=data.frame(eig_per_mat,rownames(eig_per_mat))
eig_per=melt(eig_per_mat)
#dev.off()
PC_dist = ggplot(eig_per,
                 aes(x = variable, y = value,
                     group = rownames.eig_per_mat.,
                     colour = rownames.eig_per_mat.)) +
  geom_line(aes(linetype = rownames.eig_per_mat.)) +
  geom_point() +
  xlab("Principal component rank") +
  ylab("%Variation in the PC")
PC_dist

#ggsave(PC_dist, file = "./Results/PC_dist_form.png", width = 14, height = 10, units = "cm")

##Controlling for noise
#Extend G

G_ext = lapply(G_std, function (x){ ExtendMatrix(x, ret.dim = 7)$ExtMat})
lapply(G_ext,isSymmetric)  
Ext_std_variances=lapply(G_ext,diag)
Ext_eig_variances=lapply(G_ext,function (x) {eigen(x)$values})


comp_mat=RandomSkewers(G_ext)
corr_mat=comp_mat$correlations+t(comp_mat$correlations) 
diag(corr_mat)=1
paste("Random Skewers similarity matrix")
corrplot.mixed(corr_mat,upper="number", lower="pie")




