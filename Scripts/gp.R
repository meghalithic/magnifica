# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 11 Jul 2023
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
library(tibble)
library(magrittr)
library(evoTS)
library(paleoTS)

#### LOAD DATA ----

#df <- read.csv("./Results/traits_26Jun2023.csv",
#               header = TRUE, 
#               sep = ",",
#               stringsAsFactors = FALSE)

df <- read.csv("./Results/colonies.traits.csv",
               header = TRUE, 
               sep = ",",
               stringsAsFactors = FALSE)

sm.df <- read.csv("./Results/small.colonies.traits.csv",
                  header = TRUE, 
                  sep = ",",
                  stringsAsFactors = FALSE)

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
             "#00A9FF", "#C77CFF", "#FF61CC")


col.traits = c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", 
               "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")

#### MANIPULATE DATA ----

##### CREATE ID -----

zooid_list <- unique(df$zooid.id)
length(zooid_list) #6274 (was 15773)

colony_list <- unique(df$colony.id)
length(colony_list) #572 (was 742)

##### FORMATIONS ----
# arrange formations from oldest to youngest
df$formation <- factor(df$formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                "Waipuru", "Upper Kai-Iwi", 
                                                "Tainui", "SHCSBSB")) 
formation_list <- unique(df$formation)
length(formation_list) #7

##### TRAITS -----

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
ggsave(ml, file = "./Results/trait.interest_distribution_reg.png", 
       width = 14, height = 10, units = "cm")

## most would be normal without small hump...

##### NORMALITY TESTS -----
## shapiro test; but need to subsample to be within 5000
## if significant, then significantly different from normal (i.e., non-normal)

#ln.zh
sub.ln.zh <- sample(df[, traits[1]], 
                    5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.zh) #p-value < 2.2e-16

#ln.mpw.b
sub.ln.mpw.b <- sample(df[, traits[2]], 
                       5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.mpw.b) #p-value < 2.2e-16

#ln.cw.m
sub.ln.cw.m <- sample(df[, traits[3]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.cw.m) #p-value < 0.1432

#ln.cw.d
sub.ln.cw.d <- sample(df[, traits[4]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.cw.d) #p-value < 0.003805

#ln.ow.m
sub.ln.ow.m <- sample(df[, traits[5]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.ow.m) #p-value < 2.2e-16

#ln.oh
sub.ln.oh <- sample(df[, traits[6]], 
                    5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.oh) #p-value < 2.2e-16

#ln.c.side
sub.ln.c.side <- sample(df[, traits[7]], 
                        5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.c.side) #p-value < 0.0003926

#ln.o.side
sub.ln.o.side <- sample(df[, traits[8]], 
                        5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.o.side) #p-value < 2.2e-16

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
            n.zooid = sum(n.zooid),
            avg.zooid = mean(n.zooid),
            avg.zh = mean(avg.zh, na.rm = T),
            avg.mpw.b = mean(avg.mpw.b, na.rm = T),
            avg.cw.m = mean(avg.cw.m, na.rm = T),
            avg.cw.d = mean(avg.cw.d, na.rm = T),
            avg.ow.m = mean(avg.ow.m, na.rm = T),
            avg.oh = mean(avg.oh, na.rm = T),
            avg.o.side = mean(avg.o.side, na.rm = T),
            avg.c.side = mean(avg.c.side, na.rm = T)) %>%
  as.data.frame()
## Grabowski & Porto claim sampling of 60 per sp...

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

ggsave(P_PC_dist, file = "./Results/P_PC_dist_form_reg.png", 
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
                                 ln.ow.m, ln.oh, ln.c.side, ln.o.side) ~ trait-1,
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
formation_list #order of formations
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
## chekcing NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

# NKBS
diag(phen.var[[2]]) #all larger
diag(prior$NKBS$G$G1$V)
diag(prior$NKBS$R$V)

# Waipuru
diag(phen.var[[4]]) #all larger
diag(prior$Waipuru$G$G1$V)
diag(prior$Waipuru$R$V)

# Upper Kai-Iwi
diag(phen.var[[5]]) #all larger
diag(prior$`Upper Kai-Iwi`$G$G1$V)
diag(prior$`Upper Kai-Iwi`$R$V)

###### RETRIEVE THE P MATRIX FROM THE MCMC OBJECT ------
# p matrix for each formation (maybe colony?)
## chekcing NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

#vcov(model_G[[1]]) #does not work

# NKBS
post.vcv.nkbs <- posterior.mode(model_G[[2]]$VCV)
col.vcv.nkbs <- post.vcv.nkbs[1:64] #has negatives
d.col.vcv.nkbs <- col.vcv.nkbs[c(1, 10, 19, 28, 37, 46, 55, 64)] #g matrix diagonal
unt.vcv.nkbs <- post.vcv.nkbs[65:128] #has negatives
#p.unt.vcv.nkbs <- unt.vcv.nkbs[c(1,10,19,28,37,46, 55, 64)]
#p.m.nkbs <- p.col.vcv.nkbs + p.unt.vcv.nkbs
p.m.nkbs <- col.vcv.nkbs + unt.vcv.nkbs #summed p-matrix
d.p.m.nkbs <- p.m.nkbs[c(1, 10, 9, 28, 37, 46, 55, 64)]
d.p.m.nkbs #extracted p-matrix diagonal
d.col.vcv.nkbs #is this "g" matrix smaller than the sum p matrix?
# "g' matrix all smaller than 'p' matrix [whew]
diag(Pmat[[2]]) # all larger

# Waipuru
post.vcv.wp <- posterior.mode(model_G[[4]]$VCV)
col.vcv.wp <- post.vcv.wp[1:64] #has negatives
d.col.vcv.wp <- col.vcv.wp[c(1, 10, 19, 28, 37, 46, 55, 64)]
unt.vcv.wp <- post.vcv.wp[65:128] #has negatives
#p.unt.vcv.wp <- unt.vcv.wp[c(1,10,19,28,37,46, 55, 64)]
#p.m.wp <- p.col.vcv.wp + p.unt.vcv.wp
p.m.wp <- col.vcv.wp + unt.vcv.wp
d.p.m.wp <- p.m.wp[c(1, 10, 19, 28, 37, 46, 55, 64)]
d.p.m.wp
d.col.vcv.wp #this 'g' matrix are smaller
diag(Pmat[[4]]) # cw.m; ow.m; oh; c.side; o.side are smaller than P matrix; all larger than G matrix 

# Upper Kai-Iwi
post.vcv.uki <- posterior.mode(model_G[[5]]$VCV)
col.vcv.uki <- post.vcv.uki[1:64] #has negatives
d.col.vcv.uki <- col.vcv.uki[c(1, 10, 19, 28, 37, 46, 55, 64)]
unt.vcv.uki <- post.vcv.uki[65:128] #has negatives
#p.unt.vcv.uki <- unt.vcv.uki[c(1,10,19,28,37,46, 55, 64)]
#p.m.uki <- p.col.vcv.uki + p.unt.vcv.uki
p.m.uki <- col.vcv.uki + unt.vcv.uki
d.p.m.uki <- p.m.uki[c(1, 10, 19, 28, 37, 46, 55, 64)]
d.p.m.uki
d.col.vcv.uki #this 'g' matrix are smaller
diag(Pmat[[5]]) # zh, mpw.b, cw.m, cw.d, ow.m, oh, c.side, o.side are smaller; all larger than G matrix 

###### POSTERIOR G MATRIX ------
#Retrieving G from posterior
g.model = model_G
ntraits = 8
Gmat = lapply(g.model, function (x) { 
  matrix(posterior.mode(x$VCV)[1:ntraits^2], ntraits, ntraits)})
#label lists as formations
names(Gmat) = names(by_form) #formation_list or form_data
#traits in Gmat are in different order than Pmat based on VCV 

# why aren't traits labeled??
for (i in seq_along(Gmat)){
  colnames(Gmat[[i]]) <- traits
}
for (i in seq_along(Gmat)){
  rownames(Gmat[[i]]) <- traits
}

## chekcing NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

diag(Gmat[[2]])
diag(Pmat[[2]]) # all larger

diag(Gmat[[4]])
diag(Pmat[[4]]) # all larger

diag(Gmat[[5]])
diag(Pmat[[5]]) # all larger

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

#### CORR OF G & P ----

##### CORR OF P & G DIAGONALS -----
#Gmat
#Pmat

###### REORDER TRAITS -----
formation_list

###### PLOT DIAGONALS -----

plot(diag(Gmat[[1]]), diag(Pmat[[1]]),
     pch = 19, col = col.form[1],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKLS",
     xlim = c(0, .05),
     ylim = c(0, .2))
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[1]]) ~ diag(Gmat[[1]]))) #p-value: 0.0008405 ***

plot(diag(Gmat[[2]]), diag(Pmat[[2]]),
     pch = 19, col = col.form[2],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKBS")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[2]]) ~ diag(Gmat[[2]]))) #p-value: 0.0002011 ***

plot(diag(Gmat[[3]]), diag(Pmat[[3]]),
     pch = 19, col = col.form[3],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tewkesbury",
     xlim = c(0, 0.01),
     ylim = c(0, 0.05))
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[3]]) ~ diag(Gmat[[3]]))) #p-value: 0.00396 **

plot(diag(Gmat[[4]]), diag(Pmat[[4]]),
     pch = 19, col = col.form[4],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Waipuru")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[4]]) ~ diag(Gmat[[4]]))) #p-value: 6.33e-05 ***

plot(diag(Gmat[[5]]), diag(Pmat[[5]]),
     pch = 19, col = col.form[5],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonals",
     main = "Upper Kai-Iwi")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[5]]) ~ diag(Gmat[[5]]))) #p-value: 2.999e-06 ***

plot(diag(Gmat[[6]]), diag(Pmat[[6]]),
     pch = 19, col = col.form[6],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tainui")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[6]]) ~ diag(Gmat[[6]]))) #p-value: 0.002156 **

plot(diag(Gmat[[7]]), diag(Pmat[[7]]),
     pch = 19, col = col.form[7],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[7]]) ~ diag(Gmat[[7]]))) #p-value: 2.477e-05 ***

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
g.eig_variances
p.eig_variances

plot(g.eig_variances[[1]], p.eig_variances[[1]],
     pch = 19, col = col.form[1],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "NKLS")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[1]] ~ g.eig_variances[[1]])) #p-value: 1.13e-05 ***

plot(g.eig_variances[[2]], p.eig_variances[[2]],
     pch = 19, col = col.form[2],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "NKBS")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[2]] ~ g.eig_variances[[2]])) #p-value: 1.073e-06 ***

plot(g.eig_variances[[3]], p.eig_variances[[3]],
     pch = 19, col = col.form[3],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tewkesbury")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[3]] ~ g.eig_variances[[3]])) #p-value: 3.526e-10 ***

plot(g.eig_variances[[4]], p.eig_variances[[4]],
     pch = 19, col = col.form[4],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Waipuru")
abline(0, 1) # looks okay...
summary(lm(p.eig_variances[[4]] ~ g.eig_variances[[4]])) #p-value: 1.533e-06 ***

plot(g.eig_variances[[5]], p.eig_variances[[5]],
     pch = 19, col = col.form[5],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Upper Kai-Iwi")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[5]] ~ g.eig_variances[[5]])) #p-value: 7.121e-10 ***

plot(g.eig_variances[[6]], p.eig_variances[[6]],
     pch = 19, col = col.form[6],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tainui")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[6]] ~ g.eig_variances[[6]])) #p-value: 4.171e-07 ***

plot(g.eig_variances[[7]], p.eig_variances[[7]],
     pch = 19, col = col.form[7],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "SHCSBSB")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[7]] ~ g.eig_variances[[7]])) #p-value: 3.279e-05 ***

##### PLOT P, E, AND G VARIATION ----
# E = units
# P = cov or units + colony id
# G = estimate or colony id

g.eig_vectors <- lapply(Gmat, function (x) {eigen(x)$vectors})
paste("Eigenvalue vectors")
head(g.eig_vectors)

r.names <- c(rep(formations[1], 8),
             rep(formations[2], 8),
             rep(formations[3], 8),
             rep(formations[4], 8),
             rep(formations[5], 8),
             rep(formations[6], 8),
             rep(formations[7], 8))

g.eig_vect_mat = do.call(rbind, g.eig_vectors)
g.eig_vect_mat = data.frame(g.eig_vect_mat,
                            rownames(r.names))
g.eig_vect = melt(g.eig_vect_mat)
ggplot(g.eig_vect,
       aes(x = variable, y = value,
           group = rep.rownames.g.eig_per_mat...8.,
           colour = rep.rownames.g.eig_per_mat...8.)) +
  geom_line(aes(linetype = rep.rownames.g.eig_per_mat...8.)) +
  geom_point() +
  xlab("Principal component vector") +
  ylab("%Variation in the PC")
G_PC_dist #one negative; none above 1!

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

#### EVOLVABILITY ----
##Is morphological evolution occurring along directions of high evolvability?
#Another question we might have is whether most changes are occurring along 
#dimensions of higher evolvability /conditional evolvability than we would 
#expect at random. To test that hypothesis, we can simply 
#calculate **e** and **ce** for the vectors of change and compare it to a 
#null distribution of **e** and **ce** derived from 10,000 random vectors. 

div_form_std_mean = mean_by_formation %>%
  magrittr::set_rownames(.$formation) %>%
  dplyr::select(-formation, -n.col) 
  
Bet_form = apply(div_form_std_mean, 1, Normalize)
colSums(Bet_form^2)#making sure they are normalized (should sum up to 1)

# differences standardized by ancestral mean

#split_data = split.data.frame(mean_by_formation, mean_by_formation$formation) #already age sorted
among_std_change = -diff(as.matrix(mean_by_formation[, 3:10]))

res_abs = lapply(std_change, abs)
deltazmi_vect_trait = unlist(res_abs)
deltazmi.norm <- lapply(std_change, function(x){colSums(t(x)^2)^0.5})
deltazmi_vect = unlist(deltazmi.norm)/(8^0.5)

Fig_deltaz.ana <- ggplot(as.data.frame(deltazmi_vect_trait), 
                         aes(x = deltazmi_vect_trait)) + 
  geom_histogram(color = "darkblue", fill = "lightblue") +
  xlab("DeltaZ (%mean)")

grid.arrange(Fig_deltaz.ana, 
             nrow = 1, ncol = 1)

#within formation
split_data = split.data.frame(mean_by_formation_colony, mean_by_formation_colony$formation)
within_std_change = lapply(split_data, function (y){
  -diff(as.matrix(y[, 4:11]))/as.matrix(y[-1, 4:11])
  })

#w.in_form = lapply(within_std_change, function(x) {sum(x != "NA")/8}) #8 is number of

w_form = lapply(within_std_change, function(x){apply(x, 1, Normalize)}) 
lapply(w_form, function (x){colSums(x^2)})

#Random correlations
num.traits <- 8
iterations = 400
beta.mat <- array(rnorm(num.traits * iterations), c(num.traits,iterations))
beta.mat <- apply(beta.mat, 2, Normalize)
bet_corr = t(beta.mat) %*% beta.mat
paste("95%CI for vector correlations between random vectors")
rand_corr = quantile(as.vector(abs(bet_corr)), 0.95)
rand_corr

#Deltaz vectors + Random Vectors
beta.mat.among = Bet_form
beta.mat.within = w_form

Gselect = Reduce("+", G_ext) / length(G_ext) #mean G to make drift tests during cladogenesis

evolv_within = list()
for (i in 1:5){
  evolv_within[[i]] = diag(t(beta.mat.within[[i]]) %*% G_ext[[i]] %*% beta.mat.within[[i]])
} 
evolv_among = diag(t(beta.mat.among) %*% Gselect %*% beta.mat.among)

evolv_mat_rand = diag(t(beta.mat) %*% Gselect %*% beta.mat)

evolv_final = list(unlist(evolv_within), evolv_among, evolv_mat_rand)
names(evolv_final) = c("within", "among", "random")
evolv_melted = melt(evolv_final)
evolv_ggp = ggplot(data = evolv_melted, aes(x = L1, y = value, fill = L1)) + 
  geom_boxplot(alpha = 0.6) +
  ylab("Evolvability") +
  xlab(NULL) + 
  theme(legend.position = 'none')

c_evolv_within = list()
for (i in 1:5){
  c_evolv_within[[i]]=(1/diag(t(beta.mat.within[[i]]) %*% 
                                solve(G_ext[[i]],
                                      beta.mat.within[[i]])))
} 

c_evolv_among = (1/diag(t(beta.mat.among) %*% 
                          solve(Gselect, 
                                beta.mat.among)))

c_evolv_mat_rand = (1/diag(t(beta.mat) %*% 
                             solve(Gselect, 
                                   beta.mat)))

c_evolv_final = list(unlist(c_evolv_within),
                     c_evolv_among,
                     c_evolv_mat_rand)
names(c_evolv_final) = c("within", "among", "random")
c_evolv_melted = melt(c_evolv_final)
cevol_ggp = ggplot(data = c_evolv_melted, 
                   aes(x = L1, y = value, fill = L1)) +
  geom_boxplot(alpha = 0.6) +
  ylab("Conditional Evolvability") +
  xlab(NULL) + 
  theme(legend.position = 'none')

grid.arrange(evolv_ggp, cevol_ggp, 
             nrow = 1, ncol = 2)


#### TIME SERIES ----
#### CREATE VARIABLES ----
## add time
df$ages <- c()
df$ages[df$formation == "NKLS"] <- 2
df$ages[df$formation == "NKBS"] <- 2.185
df$ages[df$formation == "Tewkesbury"] <- 1.4
df$ages[df$formation == "Waipuru"] <- 1.4
df$ages[df$formation == "Upper Kai-Iwi"] <- 0.65
df$ages[df$formation == "Tainui"] <- 0.4
df$ages[df$formation == "SHCSBSB"] <- 0.415

##### nsamp -----
time <- as.numeric(unique(sort(df$ages, decreasing = TRUE))) #order oldest to youngest
nsamp <- length(time)

##### nvar -----
nvar <- ncol(df[,grepl("^ln", colnames(df))])

##### tt -----
## from stegino_metadata/newMetadata/formations.csv

tt <- c()
for(i in 1:length(time)){
  tt[i] <- time[1]-time[i] #oldest is 0
}

##### M MATRIX & S COVARIANCE -----
M <- array(dim = c(nsamp, nvar))
#V <- array(dim = c(nsamp, nvar))
S <- list()
for(i in 1:nsamp){
  M[i,] <- colMeans(df[df$ages == time[i], grepl("^ln", colnames(df))], na.rm = TRUE)
  #V[i,] <- diag(var(df[df$ages == time[i], grepl("^ln", colnames(df))], na.rm = TRUE))
  S[[i]] <- cov(df[df$ages == time[i], grepl("^ln", colnames(df))])
}

colnames(M) <- colnames(df[,grepl("^ln", colnames(df))])
nn <- unname(table(df$ages)[as.character(time)]) #samples in order oldest to youngest
time.units <- "Myr"
sex <- "unknown"
reference <- "ROCKS-PARADOX"
taxon = "Steginoporella magnifica"

magnifica.timeSeries <- list(nsamp = nsamp,
                             nvar = nvar, 
                             M = M, 
                             S = S, 
                             nn = nn, 
                             tt = tt,
                             time.units = time.units, 
                             taxon = taxon, 
                             sex = sex, 
                             reference = reference)

source("check_TS.R") #from MacroEvolvability folder

check_TS(magnifica.timeSeries)

save(magnifica.timeSeries, file = "magnifica.timeSeries.RData")

##### GO THROUGH MODELS ----
d.var <- lapply(S, function(x){
  diag(x)
})

d.var[[1]] #all traits, first time
d.var.unlist <- unlist(d.var)
unique(names(d.var.unlist))
d.var.unlist[names(d.var.unlist) == "ln.zh"] #one trait, all times
M[1,] #all traits, first time
M[, 1] #one trait, all times

#one trait at a time?
#eventually turn into a list with names as traits
magnifica.zh <- as.paleoTS(mm = M[, 1], 
                           vv = d.var.unlist[names(d.var.unlist) == "ln.zh"], #diagonals of cov are variances 
                           tt = tt, 
                           nn = nn)

magnifica.mpw.b <- as.paleoTS(mm = M[, 2], 
                              vv = d.var.unlist[names(d.var.unlist) == "ln.mpw.b"], #diagonals of cov are variances 
                              tt = tt, 
                              nn = nn)

magnifica.cw.m <- as.paleoTS(mm = M[, 3], 
                             vv = d.var.unlist[names(d.var.unlist) == "ln.cw.m"], #diagonals of cov are variances 
                             tt = tt, 
                             nn = nn)

magnifica.cw.d <- as.paleoTS(mm = M[, 4], 
                             vv = d.var.unlist[names(d.var.unlist) == "ln.cw.d"], #diagonals of cov are variances 
                             tt = tt, 
                             nn = nn)

magnifica.ow.m <- as.paleoTS(mm = M[, 5], 
                             vv = d.var.unlist[names(d.var.unlist) == "ln.ow.m"], #diagonals of cov are variances 
                             tt = tt, 
                             nn = nn)

magnifica.oh <- as.paleoTS(mm = M[, 6], 
                           vv = d.var.unlist[names(d.var.unlist) == "ln.oh"], #diagonals of cov are variances 
                           tt = tt, 
                           nn = nn)

magnifica.c.side <- as.paleoTS(mm = M[, 7], 
                               vv = d.var.unlist[names(d.var.unlist) == "ln.c.side"], #diagonals of cov are variances 
                               tt = tt, 
                               nn = nn)

magnifica.o.side <- as.paleoTS(mm = M[, 8], 
                               vv = d.var.unlist[names(d.var.unlist) == "ln.o.side"], #diagonals of cov are variances 
                               tt = tt, 
                               nn = nn)

#visualize models
plotevoTS(magnifica.zh)
plotevoTS(magnifica.mpw.b)
plotevoTS(magnifica.cw.m)
plotevoTS(magnifica.cw.d)
plotevoTS(magnifica.ow.m)
plotevoTS(magnifica.oh) #inverse!
plotevoTS(magnifica.c.side)
plotevoTS(magnifica.o.side)

#fit models
fit.all.univariate(magnifica.zh, pool = FALSE) #URW
fit.all.univariate(magnifica.mpw.b, pool = FALSE) #URW
fit.all.univariate(magnifica.cw.m, pool = FALSE) #URW
fit.all.univariate(magnifica.cw.d, pool = FALSE) #URW
fit.all.univariate(magnifica.ow.m, pool = FALSE) #URW
fit.all.univariate(magnifica.oh, pool = FALSE) #URW
fit.all.univariate(magnifica.c.side, pool = FALSE) #URW
fit.all.univariate(magnifica.o.side, pool = FALSE) #URW

magnifica.mv <- make.multivar.evoTS(magnifica.zh,
                                    magnifica.mpw.b,
                                    magnifica.cw.m,
                                    magnifica.cw.d,
                                    magnifica.ow.m,
                                    magnifica.oh,
                                    magnifica.c.side,
                                    magnifica.o.side)
plotevoTS.multivariate(magnifica.mv,
                       y_min = -4, y_max = 8,
                       x.label = "Relative Time",
                       col = col.traits,
                       pch = c(20, 20, 20, 20,
                               20, 20, 20, 20))

fit.multivariate.URW(magnifica.mv, R = "diag", r = "fixed") #AICc = -175.344 WINNER; same as univariate
fit.multivariate.URW(magnifica.mv, R = "symmetric", r = "fixed") #AICc = -244.7925
fit.multivariate.URW(magnifica.mv, R = "symmetric", r = "accel") #AICc = -249.9289
fit.multivariate.URW(magnifica.mv, R = "symmetric", r = "decel") #AICc = -244.7943

fit.multivariate.OU(magnifica.mv, A.matrix = "diag", R.matrix = "symmetric") #AICc = -281.3557
fit.multivariate.OU(magnifica.mv, A.matrix = "diag", R.matrix = "diag") #AICc = -194.1703
fit.multivariate.OU(magnifica.mv, A.matrix = "upper.tri", R.matrix = "symmetric") #AICc = -298.1776
fit.multivariate.OU(magnifica.mv, A.matrix = "full", R.matrix = "symmetric") #AICc = -288.2023

