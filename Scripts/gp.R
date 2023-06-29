# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 28 Jun 2023
print("update 'last updated' & set working directory!")

# The purpose of this script is to create a P and G matrix for 
# Steginoporella magnifica for each age (formation) it is found
# Similar to Voje et al. where we are looking at formations rather than species

# These matrices are initially made based on the following linear measurements:
# zh = zooid height
# mpw.b = median process base width
# cw.m = cryptocyst mid-width
# cw.d = cryptocyst distal-width
# ow.m = operculum mid-width

# specimenNR = colony ID
# new.ID = individual zooid ID
# formation relates to the age

# Voje et al. used a 4 zooid/colony minimum; perhaps we should too

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

df <- read.csv("./Results/traits_26Jun2023.csv",
               header = TRUE, 
               sep = ",",
               stringsAsFactors = FALSE)
df$zooid.id <- paste0(df$boxID, "_", df$image)
colnames(df)[colnames(df) == 'specimenNR'] <- 'colony.id'


#### MANIPULATE DATA ----

# Extract unique elements and trait names

zooid_list <- unique(df$zooid.id)
length(zooid_list) #15773

colony_list <- unique(df$colony.id)
length(colony_list) #742

# arrange formations from oldest to youngest
df$formation <- factor(df$formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                "Waipuru", "Upper Kai-Iwi", 
                                                "Tainui", "SHCSBSB")) 
formation_list <- unique(df$formation)
length(formation_list) #7

df$ln.zh <- log(df$zh)
df$ln.mpw.b <- log(df$mpw.b)
df$ln.cw.m <- log(df$cw.m)
df$ln.cw.d <- log(df$cw.d)
df$ln.ow.m <- log(df$ow.m)
df$ln.oh <- log(df$oh)
df$ln.o.side <- log(df$o.side)
df$ln.c.side <- log(df$c.side)

traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

df.trim <- df %>%
  dplyr::select(zooid.id, colony.id, formation, matches(traits))

colNums <- match(c(traits, "zooid.id"), names(df.trim))
#  4  6  7  8  9 10 12 13  1

df = as.data.frame(df.trim)

#### PLOT TRAITS ----
Fig <- list ()
for (i in 1:length(traits)){
  Fig[[i]] = ggplot(data = df)+ 
    geom_density(aes(x = df[,traits[i]], 
                     group = formation,
                     col = formation)) + 
    theme(legend.position = "none", text = element_text(size = 20)) +
    scale_x_continuous(name = traits[i])
  
}

ml <- marrangeGrob(Fig, nrow = 4, ncol = 2)
ml #something wrong, when do it individually, they look different, but when loop they look the same...

#ggsave(ml, file = "./Results/trait.interest_distribution_26June2023.png", 
#       width = 14, height = 10, units = "cm")

## most would be normal without small hump...

#### REDUCE TO TRAITS OF INTEREST ----
trt_lg_N = c("formation", "colony.id", "zooid.id", traits)
dat_lg_N = df[intersect(colnames(df), trt_lg_N)]
head(dat_lg_N)

#### SUMMARY STATISTICS ----

mean_by_formation_colony = dat_lg_N %>% #use this going forward
  group_by(formation, colony.id) %>%
  summarize(n.zooid = length(unique(zooid.id)),
            avg.zh = mean(ln.zh, na.rm = T),
            avg.mpw.b = mean(ln.mpw.b, na.rm = T),
            avg.cw.m = mean(ln.cw.m, na.rm = T),
            avg.cw.d = mean(ln.cw.d, na.rm = T),
            avg.ow.m = mean(ln.ow.m, na.rm = T),
            avg.oh = mean(ln.oh, na.rm = T),
            avg.o.side = mean(ln.o.side, na.rm = T),
            avg.c.side = mean(ln.c.side, na.rm = T)) %>%
  as.data.frame()

#means of means
mean_by_formation = mean_by_formation_colony %>%
  group_by(formation) %>%
  summarize(n.col = length(unique(colony.id)),
            avg.zh = mean(avg.zh, na.rm = T),
            avg.mpw.b = mean(avg.mpw.b, na.rm = T),
            avg.cw.m = mean(avg.cw.m, na.rm = T),
            avg.cw.d = mean(avg.cw.d, na.rm = T),
            avg.ow.m = mean(avg.ow.m, na.rm = T),
            avg.oh = mean(avg.oh, na.rm = T),
            avg.o.side = mean(avg.o.side, na.rm = T),
            avg.c.side = mean(avg.c.side, na.rm = T)) %>%
  as.data.frame()

colony_means = dat_lg_N %>%
  group_by(colony.id) %>%
  summarize(formation = formation[1],
            n.zooid = length(unique(zooid.id)),
            avg.zh = mean(ln.zh, na.rm = T),
            avg.mpw.b = mean(ln.mpw.b, na.rm = T),
            avg.cw.m = mean(ln.cw.m, na.rm = T),
            avg.cw.d = mean(ln.cw.d, na.rm = T),
            avg.ow.m = mean(ln.ow.m, na.rm = T),
            avg.oh = mean(ln.oh, na.rm = T),
            avg.o.side = mean(ln.o.side, na.rm = T),
            avg.c.side = mean(ln.c.side, na.rm = T)) %>%
  as.data.frame()

means = dat_lg_N %>%
  summarize(avg.zh = mean(zh, na.rm = T),
            avg.mpw.b = mean(mpw.b, na.rm = T),
            avg.cw.m = mean(cw.m, na.rm = T),
            avg.cw.d = mean(cw.d, na.rm = T),
            avg.ow.m = mean(ow.m, na.rm = T),
            avg.o.side = mean(o.side, na.rm = T),
            avg.c.side = mean(c.side, na.rm = T),
            avg.oh = mean(oh, na.rm = T)) %>%
  as.data.frame()

## DISCRIMINANT ANALYSIS --> not needed now
# dat_lda = mean_by_formation_colony #dat_lg_N
# dat_lda[, 4:11] = scale(dat_lda[, 4:11], center = F, scale = T) #just traits
# data_discriminant = dat_lda
# 
# data_plot = na.omit(data_discriminant[, 1:11]) #all rows
# 
# r3 <- lda(formula = formation ~ ., 
#           data = data_plot[, c(1, 4:11)], method = 'mle') #just traits + formation 
# 
# plda = predict(object = r3, # predictions
#                newdata = data_plot)
# #proportion of the variance explained by each LD axis:
# prop.lda = r3$svd^2/sum(r3$svd^2) 
# 
# dataset = data.frame(formation = (data_plot)[,"formation"],
#                      lda = plda$x)
# 
# p1 <- ggplot(dataset) + 
#   geom_point(aes(lda.LD1, lda.LD2, color = formation), 
#              size = 1, alpha = .75) + 
#   labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep = ""),
#        y = paste("LD2 (", percent(prop.lda[2]), ")", sep = ""))
# p1
# ggsave(p1, file = "./Results/trait_discriminant_27Jun2023.png", 
#        width = 14, height = 10, units = "cm")

#### P MATRIX ----
#check number of colonies NOT zooids:
col_form = split.data.frame(mean_by_formation_colony, mean_by_formation_colony$formation) #zooids per formation
#just to look; max 328, smallest 19
col_form.n = lapply(col_form, function(x){dim(x)[1]})
by_form = split.data.frame(mean_by_formation_colony, mean_by_formation_colony$formation)
p.form_data = lapply(by_form, function(x) x[complete.cases(x),])
p.cov = lapply(p.form_data, function (x){ (cov(x[, 4:11]))}) #traits per colony (not variation within colony)


p.model = p.cov
p.data = (colony_means)
ntraits = 8
Pmat = lapply(model, function (x) { 
  matrix(posterior.mode(x$VCV)[1:ntraits^2], ntraits, ntraits)})


###### P STD ------

#Standardizing G by the trait means 

mean_by_form.col = setDT(na.omit(p.data[, c(2, 4:11)]))[, lapply(.SD, mean, na.rm = F),
                                                        by = .(formation)] #traits + formation
p.u_form = split(mean_by_form.col, mean_by_form.col$formation)
p.test = lapply(p.u_form, function (x){ data.matrix(x[, 2:9])}) #only traits from mean_by_form
p.test_std = lapply(p.test, function (x){(as.numeric(x))%*%t(as.numeric(x))})
P_std = list()
for (i in 1:length(Pmat)){
  P_std[[i]] = Pmat[[i]]/(p.test_std[names(p.form_data[i])][[1]])
}
P_std
names(P_std)=names(p.form_data[1:i])

##Genetic variance in traits and eigen vectors

#load(file="New_g_matrices.RData") #load the g matrices calculated above 

lapply(P_std, isSymmetric)  #is.symmetric.matrix
p.std_variances = lapply(P_std, diag)
paste("Trait variances")
head(p.std_variances)

###### P EIGEN ------

p.eig_variances=lapply(P_std, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(p.eig_variances)

p.eig_percent = lapply(p.eig_variances, function (x) {x/sum(x)})
p.eig_per_mat = do.call(rbind, p.eig_percent)
p.eig_per_mat = data.frame(p.eig_per_mat, rownames(p.eig_per_mat))
p.eig_per = melt(p.eig_per_mat)
#dev.off()
P_PC_dist = ggplot(p.eig_per,
                   aes(x = variable, y = value,
                       group = rownames.p.eig_per_mat.,
                       colour = rownames.p.eig_per_mat.)) +
  geom_line(aes(linetype = rownames.p.eig_per_mat.)) +
  geom_point() +
  xlab("Principal component rank") +
  ylab("%Variation in the PC")
P_PC_dist #one negative; none above 1!

ggsave(P_PC_dist, file = "./Results/P_PC_dist_form.png", 
       width = 14, height = 10, units = "cm")

###### P NOISE ------
##Controlling for noise
#Extend G
P_ext = lapply(P_std, function (x){ ExtendMatrix(x, ret.dim = 7)$ExtMat}) #not 8 because last eigen value (#8) was negative
#ignore warning from above
lapply(P_ext, isSymmetric)  
P_Ext_std_variances = lapply(P_ext, diag)
P_Ext_eig_variances = lapply(P_ext, function (x) {eigen(x)$values})

#### G MATRIX ----
#keep at zooid level because correct for this later
zoo_form = split.data.frame(dat_lg_N, dat_lg_N$formation) #zooids per formation
#just to look; highest 7836, smallest 454
zoo_form.n = lapply(zoo_form, function(x){dim(x)[1]})
zooid_by_form = split.data.frame(dat_lg_N, dat_lg_N$formation)
zooid_form_data = lapply(zooid_by_form, function(x) x[complete.cases(x),])

##### PRIORS -----
phen.var = lapply(zooid_form_data, function (x){ (cov(x[, 4:11]))}) #traits of ALL; correct for colony later
prior = lapply(phen.var, function (x){list(G = list(G1 = list(V = x/2, nu = 2)),
                                           R = list(V = x/4, nu = 2))})


##### MCMC -----
#Running the MCMC chain
model_G = list()
for (i in 1:length(formation_list)){ #length 7 because 7 formations
  model_G[[i]] <- MCMCglmm(cbind(ln.zh, ln.mpw.b, ln.cw.m, ln.cw.d, #same order as in priors
                                 ln.ow.m, ln.c.side, ln.o.side, ln.oh) ~ trait-1,
                         #account for variation w/in colony:
                         random = ~us(trait):colony.id, #the number of these determines # of Gs #+ us(trait):formation
                         rcov = ~us(trait):units,
                         family = rep("gaussian", 8), #num of traits
                         data = zooid_form_data[[i]],
                         nitt = 1500000, thin = 1000, burnin = 500000,
                         prior = prior[[i]], verbose = TRUE)
  
}

data.list = list(model_G, dat_lg_N, zooid_form_data, mean_by_formation_colony)

save(data.list, file = "./Results/g_matrices_data_form.RData")

#load(file = "./Results/g_matrices_data_form.RData") #load the g matrices calculated above 

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
#formations from oldest to youngest: "NKLS", "NKBS", "Tewkesbury", "Waipuru", 
#                                    "Upper Kai-Iwi", "SHCSBSB", "Tainui"

###### G MATRIX ------
#Retrieving G from posterior
g.model = model_G
g.data = (dat_lg_N)
ntraits = 8
Gmat = lapply(g.model, function (x) { 
  matrix(posterior.mode(x$VCV)[1:ntraits^2], ntraits, ntraits)})


###### G STD ------

#Standardizing G by the trait means 

mean_by_form = setDT(na.omit(g.data[, 3:11]))[, lapply(.SD, mean, na.rm = F), 
                                              by = .(formation)] #traits + formation
g.u_form = split(mean_by_form, mean_by_form$formation)
g.test = lapply(g.u_form, function (x){ data.matrix(x[, 2:9])}) #only traits from mean_by_form
g.test_std = lapply(g.test, function (x){(as.numeric(x))%*%t(as.numeric(x))})
G_std = list()
for (i in 1:length(Gmat)){
  G_std[[i]] = Gmat[[i]]/(g.test_std[names(zooid_form_data[i])][[1]]) #was p.form_data
}
G_std
names(G_std) = names(form_data[1:i])

##Genetic variance in traits and eigen vectors

#load(file="New_g_matrices.RData") #load the g matrices calculated above 

lapply(G_std, isSymmetric)  #is.symmetric.matrix
g.std_variances = lapply(G_std, diag)
paste("Trait variances")
head(g.std_variances)

#require(matrixcalc)
#m.1 <- round(G_std[[1]], 10)
#is.symmetric.matrix(m.1)

#is.positive.definite(m.1) 
#no spot with zero variance; 
#so variance along certain directions are negative

###### G EIGEN ------

g.eig_variances = lapply(G_std, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(g.eig_variances)

g.eig_percent = lapply(g.eig_variances, function (x) {x/sum(x)})
g.eig_per_mat = do.call(rbind, g.eig_percent)
g.eig_per_mat = data.frame(g.eig_per_mat, rownames(g.eig_per_mat))
g.eig_per = melt(g.eig_per_mat)
#dev.off()
G_PC_dist = ggplot(g.eig_per,
                   aes(x = variable, y = value,
                       group = rownames.g.eig_per_mat.,
                       colour = rownames.g.eig_per_mat.)) +
  geom_line(aes(linetype = rownames.g.eig_per_mat.)) +
  geom_point() +
  xlab("Principal component rank") +
  ylab("%Variation in the PC")
G_PC_dist #one negative; none above 1!

ggsave(G_PC_dist, file = "./Results/G_PC_dist_form.png", 
       width = 14, height = 10, units = "cm")

#Note that some matrices have negative eigenvalues. 
#This can cause a lot of trouble in analyses involving inverted matrices.
#Solution from evolqg Marroig et al. 2012

###### G NOISE ------
##Controlling for noise
#Extend G
G_ext = lapply(G_std, function (x){ ExtendMatrix(x, ret.dim = 7)$ExtMat}) #not 8 because last eigen value (#8) was negative
#ignore warning from above
lapply(G_ext, isSymmetric)  
Ext_std_variances = lapply(G_ext, diag)
Ext_eig_variances = lapply(G_ext, function (x) {eigen(x)$values})
##need to create random cov.m for comparison
#cov.m <- RandomMatrix(8, 1, 1, 100) 
#G_list <- list(G_ext$NKLS, cov.m)
comp_mat = RandomSkewers(G_ext) #need at least
corr_mat = comp_mat$correlations + t(comp_mat$correlations) 
diag(corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(corr_mat,upper = "number", lower = "pie")

##### RAREFACTION -----
##Rarefaction - 
#How much difference in matrix structure is induced by sampling error?

#Rarefaction of G
##Rarefaction - How much difference in matrix structure is induced by sampling error?

# We do a rarefaction analyses and compare a 'true' matrix with an 
# under-sampled version of itself and see whether the dissimilarities observed 
# across species can be explained by sampling alone.

# Computing a pooled P across all individuals in all samples
#P_pooled<-var(cbind(indata$tpg, indata$ds2, indata$lpt, indata$ds3), na.rm=TRUE)
#G_test <- G_ext$NKLS
#trait_means_global<-c(mean(indata$tpg, na.rm=TRUE), 
#                      mean(indata$ds2, na.rm=TRUE), 
#                      mean(indata$lpt, na.rm=TRUE), 
#                      mean(indata$ds3, na.rm=TRUE))
#P_pooled_mean_std <- meanStdG(P_pooled, trait_means_global)

# We start by generating "individuals" (i.e., colonies) 
# in a population with VCV patterns corresponding to a known G. 
# We then calculate VCV matrices for those individuals and obtain their 
# similarity with the 'true' G.
repititions <- 200
#Defining the population size. Each population size is repeated 200 times. 
#Starting at 2 individuals and going up to 100, 
col_form.n # max 328, min 19
pop.size = sort(rep(c(2:500), repititions)) #make 500 so much bigger than max no. colonies
pop.dist <- lapply(pop.size, function (x){mvrnorm(n = x, 
                                                  rep(0,8), # mean 0 for 8 traits 
                                                  G_ext[[1]])}) # Simulate individual trait data based on the first VCOV matrix in our time series and different sample sizes. 
# Calculate VCOV matrices 
pop.cv <- lapply(pop.dist, cov) 
# Estimate selection gradients based on Random Skewers, comparing our true 
# G matrix and the undersampled pooled P matrix. 
RS_result = RandomSkewers(pop.cv, 
                          G_ext[[1]]) 
out_results <- cbind(RS_result$correlation, pop.size)

# Computing Random Skewers between pairs of actual G matrices in our data
#(n = 7 for 7 formations) 
# sub-setting so that we only analyze the traits and formation of interest
mean.df.sub <- mean_by_formation_colony[, c(1, 4:11)] 
#remove rows with NA so that the sample size gets correct when we plot 
#the sample size for the VCOV.
mean.df_complete_cases <- mean.df.sub[complete.cases(mean.df.sub), ] #remove rows with NA so that the sample size gets correct when we plot the sample size for the VCOV.

# Calculating the sample size for each time point
colony_samples = split.data.frame(mean.df_complete_cases, mean.df_complete_cases$formation) #Sample size
sample_sizes_G = lapply(colony_samples, function(x){dim(x)[1]})

comp_sampleN = matrix(0, 7, 7) #calculating the smallest sample size in the comparison among pairs of G 

for (i in 1:length(sample_sizes_G)){
  for (j in 1:length(sample_sizes_G)){
    comp_sampleN[i,j]=min(as.numeric(sample_sizes_G[i]),as.numeric(sample_sizes_G[j]))
  }
}

# Here, we do the actual Random Skewers test and 
# combined the results (correlations in the response to selection) 
# with the smallest sample size for the pairs of G that are investigated
comp_mat = RandomSkewers(G_ext)
melt_comp = melt(comp_mat$correlations[lower.tri(comp_mat$correlations)])
melt_samples = melt(comp_sampleN[lower.tri(comp_sampleN)])
obs_melt = cbind(melt_comp, melt_samples)
colnames(obs_melt) = c("RS","N")

#Plotting the result
plot(out_results[, 2], out_results[, 1], 
       xlab = "Sample size", ylab = "Similarity", 
       pch = 19, col = "grey")
points(obs_melt$N, obs_melt$RS, 
       col = "#00BFC4", pch = 19) #color by formation

#### CORR OF G FOR P ----

#formations and colors: 
#NKLS = #F8766D
#NKBS = #CD9600
#Twekesbury = #7CAE00
#Waipuru = #00BE67
#Upper Kai-Iwi = #00A9FF
#Tainui = #C77CFF
#SHCSBSB = #FF61CC
col.form = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00A9FF", "#C77CFF", "#FF61CC")

##### RANDOM SKEWERS OF P & G OF EACH FORMATION -----
#NKLS
NKLS_comp_mat = RandomSkewers(list(G_ext[[1]], P_ext[[1]])) #need at least
NKLS_corr_mat = NKLS_comp_mat$correlations + t(NKLS_comp_mat$correlations) 
diag(NKLS_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(NKLS_corr_mat, upper = "number", lower = "pie")

#NKBS
NKBS_comp_mat = RandomSkewers(list(G_ext[[2]], P_ext[[2]])) #need at least
NKBS_corr_mat =NKBS_comp_mat$correlations + t(NKBS_comp_mat$correlations) 
diag(NKBS_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(NKBS_corr_mat, upper = "number", lower = "pie")

#Tewkesbury
Tewkesbury_comp_mat = RandomSkewers(list(G_ext[[3]], P_ext[[3]])) #need at least
Tewkesbury_corr_mat = Tewkesbury_comp_mat$correlations + t(Tewkesbury_comp_mat$correlations) 
diag(Tewkesbury_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Tewkesbury_corr_mat, upper = "number", lower = "pie")

#Waipuru
Waipuru_comp_mat = RandomSkewers(list(G_ext[[4]], P_ext[[4]])) #need at least
Waipuru_corr_mat = Waipuru_comp_mat$correlations + t(Waipuru_comp_mat$correlations) 
diag(Waipuru_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Waipuru_corr_mat, upper = "number", lower = "pie")

#Upper Kai-Iwi
UKai_Iwi_comp_mat = RandomSkewers(list(G_ext[[5]], P_ext[[5]])) #need at least
UKai_Iwi_corr_mat = UKai_Iwi_comp_mat$correlations + t(UKai_Iwi_comp_mat$correlations) 
diag(UKai_Iwi_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(UKai_Iwi_corr_mat, upper = "number", lower = "pie")

#Tainui
Tainui_comp_mat = RandomSkewers(list(G_ext[[6]], P_ext[[6]])) #need at least
Tainui_corr_mat = Tainui_comp_mat$correlations + t(Tainui_comp_mat$correlations) 
diag(Tainui_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Tainui_corr_mat, upper = "number", lower = "pie")

#SHCSBSB
SHCSBSB_comp_mat = RandomSkewers(list(G_ext[[7]], P_ext[[7]])) #need at least
SHCSBSB_corr_mat = SHCSBSB_comp_mat$correlations + t(SHCSBSB_comp_mat$correlations) 
diag(SHCSBSB_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(SHCSBSB_corr_mat, upper = "number", lower = "pie")

##### COMPARISON OF PC1 AND PC2 OF P & G OF EACH FORMATION -----
g.std_variances
p.std_variances

plot(g.std_variances[[1]], p.std_variances[[1]],
     pch = 19, col = col.form[1],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "NKLS")
abline(0,1)
summary(lm(p.std_variances[[1]] ~ g.std_variances[[1]]))

plot(g.std_variances[[2]], p.std_variances[[2]],
     pch = 19, col = col.form[2],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "NKBS")
abline(0,1)
summary(lm(p.std_variances[[2]] ~ g.std_variances[[2]]))

plot(g.std_variances[[3]], p.std_variances[[3]],
     pch = 19, col = col.form[3],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tewkesbury")
abline(0,1)
summary(lm(p.std_variances[[3]] ~ g.std_variances[[3]]))

plot(g.std_variances[[4]], p.std_variances[[4]],
     pch = 19, col = col.form[4],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Waipuru")
abline(0,1)
summary(lm(p.std_variances[[4]] ~ g.std_variances[[4]]))

plot(g.std_variances[[5]], p.std_variances[[5]],
     pch = 19, col = col.form[5],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Upper Kai-Iwi")
abline(0,1)
summary(lm(p.std_variances[[5]] ~ g.std_variances[[5]]))

plot(g.std_variances[[6]], p.std_variances[[6]],
     pch = 19, col = col.form[6],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tainui")
abline(0,1)
summary(lm(p.std_variances[[6]] ~ g.std_variances[[6]]))

plot(g.std_variances[[7]], p.std_variances[[7]],
     pch = 19, col = col.form[7],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "SHCSBSB")
abline(0,1)
summary(lm(p.std_variances[[7]] ~ g.std_variances[[7]]))
