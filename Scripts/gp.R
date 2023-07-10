# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 5 Jul 2023
print("update 'last updated' & set working directory!")
# set working directory to repo "magnifica"

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

#### PLOT THEME ----
#formations and colors: 
#NKLS = #F8766D
#NKBS = #CD9600
#Twekesbury = #7CAE00
#Waipuru = #00BE67
#Upper Kai-Iwi = #00A9FF
#Tainui = #C77CFF
#SHCSBSB = #FF61CC
col.form = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00A9FF", "#C77CFF", "#FF61CC")

#### MANIPULATE DATA ----

##### CREATE ID -----
# Extract unique elements and trait names

zooid_list <- unique(df$zooid.id)
length(zooid_list) #15773

colony_list <- unique(df$colony.id)
length(colony_list) #742

##### FORMATIONS ----
# arrange formations from oldest to youngest
df$formation <- factor(df$formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                "Waipuru", "Upper Kai-Iwi", 
                                                "Tainui", "SHCSBSB")) 
formation_list <- unique(df$formation)
length(formation_list) #7

##### LN TRANSFORM -----
df$ln.zh <- log(df$zh)
df$ln.mpw.b <- log(df$mpw.b)
df$ln.cw.m <- log(df$cw.m)
df$ln.cw.d <- log(df$cw.d)
df$ln.ow.m <- log(df$ow.m)
df$ln.oh <- log(df$oh)
df$ln.o.side <- log(df$o.side)
df$ln.c.side <- log(df$c.side)

#same order as in df
names(df)
traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

##### TRIM DATASET ----
df.trim <- df %>%
  dplyr::select(zooid.id, colony.id, formation, matches(traits))

colNums <- match(c(traits, "zooid.id"), names(df.trim))
#  4  6  7  8  9 10 11  1 (i.e., 4:11 are traits of interest)

df = as.data.frame(df.trim)

#### PLOT TRAITS ----

##### DISTRIBUTIONS -----
p.zh = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[1]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[1]) +
  scale_color_manual(values = col.form)

p.mpw.b = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[2]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[2]) +
  scale_color_manual(values = col.form)

p.cw.m = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[3]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[3]) +
  scale_color_manual(values = col.form)

p.cw.d = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[4]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[4]) +
  scale_color_manual(values = col.form)

p.ow.m = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[5]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[5]) +
  scale_color_manual(values = col.form)

p.oh = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[6]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[6]) +
  scale_color_manual(values = col.form)

p.c.side = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[7]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[7]) +
  scale_color_manual(values = col.form)

p.o.side = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[8]], 
                   group = formation,
                   col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = traits[8]) +
  scale_color_manual(values = col.form)

Fig = list(p.zh, p.mpw.b, p.cw.m, p.cw.d, p.ow.m, p.oh, p.c.side, p.o.side)
ml <- marrangeGrob(Fig, nrow = 4, ncol = 2)
ml
ggsave(ml, file = "./Results/trait.interest_distribution.png", 
       width = 14, height = 10, units = "cm")

## most would be normal without small hump...

#### REDUCE TO TRAITS OF INTEREST ----
trt_lg_N = c("formation", "colony.id", "zooid.id", traits)
dat_lg_N = df[intersect(colnames(df), trt_lg_N)]
head(dat_lg_N) #traits in same order as df and traits

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
  summarize(avg.zh = mean(ln.zh, na.rm = T),
            avg.mpw.b = mean(ln.mpw.b, na.rm = T),
            avg.cw.m = mean(ln.cw.m, na.rm = T),
            avg.cw.d = mean(ln.cw.d, na.rm = T),
            avg.ow.m = mean(ln.ow.m, na.rm = T),
            avg.oh = mean(ln.oh, na.rm = T),
            avg.o.side = mean(ln.o.side, na.rm = T),
            avg.c.side = mean(ln.c.side, na.rm = T)) %>%
  as.data.frame()

#### DISCRIMINANT ANALYSIS ----- 
# not needed now
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


#### CHECK SAMPLE SIZES ----
## number of zooids per colony
range(mean_by_formation_colony$n.zooid)
#4 53

#check number of zooids NOT colonies:
# by colonies use mean_by_formation_colony
# by zooid us dat_lg_N
col_form = split.data.frame(mean_by_formation_colony,  #by colonies
                            mean_by_formation_colony$formation) #zooids per formation
#just to look; max 328, smallest 19
col_form.n = lapply(col_form, function(x){dim(x)[1]})

#### SPLIT BY FORMATION ----
## by zooids:
by_form = split.data.frame(dat_lg_N, 
                           dat_lg_N$formation)
#just to look; highest 7836, smallest 454
by_form.n = lapply(by_form, function(x){dim(x)[1]})
form_data = lapply(by_form, function(x) x[complete.cases(x),])

#### P MATRIX ----
p.cov = lapply(form_data, function (x){ (cov(x[, 4:11]))}) #traits per colony (not variation within colony)
#p.cov is the same and the phen.var

###### P STD ------

#Standardizing P by the trait means 
#p.data = colony_means

#mean_by_form.col = setDT(na.omit(p.data[, c(2, 4:11)]))[, lapply(.SD, mean, na.rm = F),
#                                                        by = .(formation)] #traits + formation
#p.u_form = split(mean_by_form.col, mean_by_form.col$formation)
#p.test = lapply(p.u_form, function (x){ data.matrix(x[, 2:9])}) #only traits from mean_by_form
#p.test_std = lapply(p.test, function (x){(as.numeric(x))%*%t(as.numeric(x))})
#P_std = list()
#for (i in 1:length(p.cov)){
#  P_std[[i]] = p.cov[[i]]/(p.test_std[names(p.form_data[i])][[1]])
#}
#P_std
#names(P_std)=names(p.form_data[1:i])


##### P VARIANCES ----
##Phenotypic variance in traits and eigen vectors
Pmat = p.cov

lapply(Pmat, isSymmetric)  #is.symmetric.matrix
p.variances = lapply(Pmat, diag)
paste("Trait variances")
head(p.variances)

###### P EIGEN ------

p.eig_variances = lapply(Pmat, function (x) {eigen(x)$values})
# lapply(Pmat, function (x) {eigen(x)})
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
P_PC_dist #none negative; none above 1!

ggsave(P_PC_dist, file = "./Results/P_PC_dist_form.png", 
       width = 14, height = 10, units = "cm") 
#NKBS, Waipuru, Upper Kai-Iwi 

###### P NOISE ------
##Controlling for noise
#Extend G
P_ext = lapply(Pmat, function (x){ ExtendMatrix(x, ret.dim = 6)$ExtMat}) #not 8 because last eigen value (#8) was negative
#ignore warning from above
lapply(P_ext, isSymmetric)  
P_Ext_std_variances = lapply(P_ext, diag)
P_Ext_eig_variances = lapply(P_ext, function (x) {eigen(x)$values})

p.comp_mat = RandomSkewers(P_ext) #need at least
p.corr_mat = p.comp_mat$correlations + t(p.comp_mat$correlations) 
diag(p.corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(p.corr_mat, upper = "number", lower = "pie")


#### G MATRIX ----

##### PRIORS -----
#same as p.cov
phen.var = lapply(form_data, function (x){ (cov(x[, 4:11]))}) #traits of ALL; correct for colony later
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
                         data = form_data[[i]],
                         nitt = 1500000, thin = 1000, burnin = 500000,
                         prior = prior[[i]], verbose = TRUE)
  
}

data.list = list(model_G, dat_lg_N, form_data, mean_by_formation_colony)

save(data.list, file = "./Results/g_matrices_data_form.RData")

#load(file="./Results/g_matrices_data_form.RData") #load the g matrices calculated above 
#model_G <- data.list[[1]]

##### CHECK MODELS -----
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


##### CHECK P IS BIGGER THAN G -----

###### PRIORS ------
# diagonals of p.cov to priors
diag(phen.var[[2]]) #all larger
diag(prior$NKBS$G$G1$V)
diag(prior$NKBS$R$V)

diag(phen.var[[4]])
diag(prior$Waipuru$G$G1$V)
diag(prior$Waipuru$R$V)


diag(phen.var[[5]])
diag(prior$`Upper Kai-Iwi`$G$G1$V)
diag(prior$`Upper Kai-Iwi`$R$V)

###### RETRIEVE THE P MATRIX FROM THE MCMC OBJECT ------
# p matrix for each formation (maybe colony?)
#vcov(model_G[[1]])
post.vcv <- posterior.mode(model_G[[1]]$VCV)
col.vcv <- post.vcv[1:64]
d.col.vcv <- col.vcv[c(1,10,19,28,37,46, 55, 64)]
unt.vcv <- post.vcv[65:128]
p.unt.vcv <- unt.vcv[c(1,10,19,28,37,46, 55, 64)]
diag(Pmat[[1]])
diag(Gmat[[1]])
#zh is same as col.id, bigger than units
#mpb.w is same as col.id, bigger than units
#cw.m is same as col.id, smaller than units
#cw.d is same as col.id, bigger than units
#ow.m is same as col.id, bigger than units
#c.side is same as col.id, bigger than units
#o.side is same as col.id, bigger than units
#oh is same as col.id, bigger than units

###### POSTERIOR G MATRIX ------
#Retrieving G from posterior
g.model = model_G
ntraits = 8
Gmat = lapply(g.model, function (x) { 
  matrix(posterior.mode(x$VCV)[1:ntraits^2], ntraits, ntraits)})
#label lists as formations
names(Gmat) = names(by_form) #formation_list
#traits in Gmat are in different order than Pmat based on VCV 

# why aren't traits labeled??
for (i in seq_along(Gmat)){
  colnames(Gmat[[i]]) <- c("ln.zh", "ln.mpb.w", "ln.cw.m", "ln.cw.d",
                           "ln.ow.m", "ln.c.side", "ln.o.side", "ln.oh")
}
for (i in seq_along(Gmat)){
  rownames(Gmat[[i]]) <- c("ln.zh", "ln.mpb.w", "ln.cw.m", "ln.cw.d",
                           "ln.ow.m", "ln.c.side", "ln.o.side", "ln.oh")
}

diag(Gmat[[2]])
diag(Pmat[[2]])

diag(Gmat[[4]])
diag(Pmat[[4]])

diag(Gmat[[5]])
diag(Pmat[[5]])

###### G STD ------

#Standardizing G by the trait means 
#g.data = (dat_lg_N)

#of colonies
#mean_by_form = setDT(na.omit(g.data[, 3:11]))[, lapply(.SD, mean, na.rm = F), 
  #                                            by = .(formation)] #traits + formation
#g.u_form = split(mean_by_form, mean_by_form$formation)
#g.test = lapply(g.u_form, function (x){ data.matrix(x[, 2:9])}) #only traits from mean_by_form
#g.test_std = lapply(g.test, function (x){(as.numeric(x))%*%t(as.numeric(x))})
#G_std = list()
#for (i in 1:length(Gmat)){
#  G_std[[i]] = Gmat[[i]]/(g.test_std[names(zooid_form_data[i])][[1]]) #was p.form_data
#}
#G_std
#names(G_std) = names(form_data[1:i])

##Genetic variance in traits and eigen vectors

##### G VARIANCES -----
lapply(Gmat, isSymmetric)  #is.symmetric.matrix
g.variances = lapply(Gmat, diag)
paste("Trait variances")
head(g.variances)

#require(matrixcalc)
#m.1 <- round(G_std[[1]], 10)
#is.symmetric.matrix(m.1)

#is.positive.definite(m.1) 
#no spot with zero variance; 
#so variance along certain directions are negative

###### G EIGEN ------

g.eig_variances = lapply(Gmat, function (x) {eigen(x)$values})
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
G_ext = lapply(Gmat, function (x){ ExtendMatrix(x, ret.dim = 6)$ExtMat}) #not 8 because last eigen value (#8) was negative
#ignore warning from above
lapply(G_ext, isSymmetric)  
Ext_std_variances = lapply(G_ext, diag)
Ext_eig_variances = lapply(G_ext, function (x) {eigen(x)$values})
##need to create random cov.m for comparison
#cov.m <- RandomMatrix(8, 1, 1, 100) 
#G_list <- list(G_ext)

g.comp_mat = RandomSkewers(G_ext) #need at least
g.corr_mat = g.comp_mat$correlations + t(g.comp_mat$correlations) 
diag(g.corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(g.corr_mat,upper = "number", lower = "pie")

#### CORR OF G FOR P ----

##### CORR OF P & G DIAGONALS -----
#Gmat
#Pmat

d.gmat.nkbs <- diag(Gmat[[1]])[1,2,3,4,5,8,6,7]
plot(diag(Gmat[[1]]), diag(Pmat[[1]]),
     pch = 19, col = col.form[1],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKLS",
     xlim = c(0, .05),
     ylim = c(0, .2))
abline(0, 1)
summary(lm(diag(Pmat[[1]]) ~ diag(Gmat[[1]])))

plot(diag(Gmat[[2]]), diag(Pmat[[2]]),
     pch = 19, col = col.form[2],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKBS")
abline(0, 1)
summary(lm(diag(Pmat[[2]]) ~ diag(Gmat[[2]])))

plot(diag(Gmat[[3]]), diag(Pmat[[3]]),
     pch = 19, col = col.form[3],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tewkesbury",
     xlim = c(0, 0.01),
     ylim = c(0, 0.05))
abline(0, 1)
summary(lm(diag(Pmat[[3]]) ~ diag(Gmat[[3]])))

plot(diag(Gmat[[4]]), diag(Pmat[[4]]),
     pch = 19, col = col.form[4],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Waipuru")
abline(0, 1)
summary(lm(diag(Pmat[[4]]) ~ diag(Gmat[[4]])))

plot(diag(Gmat[[5]]), diag(Pmat[[5]]),
     pch = 19, col = col.form[5],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonals",
     main = "Upper Kai-Iwi")
abline(0, 1)
summary(lm(diag(Pmat[[5]]) ~ diag(Gmat[[5]])))

plot(diag(Gmat[[6]]), diag(Pmat[[6]]),
     pch = 19, col = col.form[6],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tainui")
abline(0, 1)
summary(lm(diag(Pmat[[6]]) ~ diag(Gmat[[6]])))

plot(diag(Gmat[[7]]), diag(Pmat[[7]]),
     pch = 19, col = col.form[7],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1)
summary(lm(diag(Pmat[[7]]) ~ diag(Gmat[[7]])))

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
abline(0, 1)
summary(lm(p.std_variances[[1]] ~ g.std_variances[[1]]))

plot(g.std_variances[[2]], p.std_variances[[2]],
     pch = 19, col = col.form[2],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "NKBS")
abline(0, 1)
summary(lm(p.std_variances[[2]] ~ g.std_variances[[2]]))

plot(g.std_variances[[3]], p.std_variances[[3]],
     pch = 19, col = col.form[3],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tewkesbury")
abline(0, 1)
summary(lm(p.std_variances[[3]] ~ g.std_variances[[3]]))

plot(g.std_variances[[4]], p.std_variances[[4]],
     pch = 19, col = col.form[4],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Waipuru")
abline(0, 1)
summary(lm(p.std_variances[[4]] ~ g.std_variances[[4]]))

plot(g.std_variances[[5]], p.std_variances[[5]],
     pch = 19, col = col.form[5],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Upper Kai-Iwi")
abline(0, 1)
summary(lm(p.std_variances[[5]] ~ g.std_variances[[5]]))

plot(g.std_variances[[6]], p.std_variances[[6]],
     pch = 19, col = col.form[6],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tainui")
abline(0, 1)
summary(lm(p.std_variances[[6]] ~ g.std_variances[[6]]))

plot(g.std_variances[[7]], p.std_variances[[7]],
     pch = 19, col = col.form[7],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "SHCSBSB")
abline(0, 1)
summary(lm(p.std_variances[[7]] ~ g.std_variances[[7]]))

##### RAREFACTION -----
##Rarefaction - 
#How much difference in matrix structure is induced by sampling error?

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
