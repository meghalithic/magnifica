# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 29 Aug 2023
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
library(coin)
library(corrplot)
library(data.table)
library(dplyr)
library(evolqg)
library(evolvability)
library(ggplot2)
library(grid)
library(gridBase)
library(gridExtra)
library(magrittr)
library(MASS)
library(MCMCglmm)
library(nse)
library(RColorBrewer)
library(reshape2)
library(rgl)
library(scales)
library(scatterplot3d)
library(tibble)
library(matrixcalc)

#### LOAD DATA ----

#df <- read.csv("./Results/traits_26Jun2023.csv",
#               header = TRUE, 
#               sep = ",",
#               stringsAsFactors = FALSE)

df <- read.csv("./Results/colonies.traits.csv",
               header = TRUE, 
               sep = ",",
               stringsAsFactors = FALSE)
#already have small zooid removed and at least 5 zooids per colony

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

zooid_list <- unique(df$zooid.id)
length(zooid_list) #5480 (was 15773)

colony_list <- unique(df$colony.id)
length(colony_list) #541 (was 742)

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
shapiro.test(sub.ln.mpw.b) #p-value = 1.391e-07

#ln.cw.m
sub.ln.cw.m <- sample(df[, traits[3]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.cw.m) #p-value = 0.02664

#ln.cw.d
sub.ln.cw.d <- sample(df[, traits[4]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.cw.d) #p-value = 1.302e-05

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
shapiro.test(sub.ln.c.side) #p-value = 0.001001

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
  dplyr::group_by(formation, colony.id) %>%
  dplyr::summarize(n.zooid = length(zooid.id),
            avg.zh = mean(ln.zh, na.rm = T),
            sd.zh = sd(ln.zh, na.rm = T),
            avg.mpw.b = mean(ln.mpw.b, na.rm = T),
            sd.mpw.b = sd(ln.mpw.b, na.rm = T),
            avg.cw.m = mean(ln.cw.m, na.rm = T),
            sd.cw.m = sd(ln.cw.m, na.rm = T),
            avg.cw.d = mean(ln.cw.d, na.rm = T),
            sd.cw.d = sd(ln.cw.d, na.rm = T),
            avg.ow.m = mean(ln.ow.m, na.rm = T),
            sd.ow.m = sd(ln.ow.m, na.rm = T),
            avg.oh = mean(ln.oh, na.rm = T),
            sd.oh = sd(ln.oh, na.rm = T),
            avg.o.side = mean(ln.o.side, na.rm = T),
            sd.o.side = sd(ln.o.side, na.rm = T),
            avg.c.side = mean(ln.c.side, na.rm = T),
            sd.c.side = sd(ln.c.side, na.rm = T)) %>%
  as.data.frame()
min(mean_by_formation_colony$n.zooid) #5

mean(mean_by_formation_colony$sd.zh) #0.09225405
range(mean_by_formation_colony$sd.zh) #0.02640565 0.57287376

mean(mean_by_formation_colony$sd.mpw.b) #0.1217314
range(mean_by_formation_colony$sd.mpw.b) #0.03438424 0.36673432

mean(mean_by_formation_colony$sd.cw.m) #0.1490322
range(mean_by_formation_colony$sd.cw.m) #0.01840786 0.36565782

mean(mean_by_formation_colony$sd.cw.d) #0.1103786
range(mean_by_formation_colony$sd.cw.d) #0.0219327 0.3552131

mean(mean_by_formation_colony$sd.ow.m) #0.07753377
range(mean_by_formation_colony$sd.ow.m) #0.01233698 0.31205832

mean(mean_by_formation_colony$sd.oh) #0.075358
range(mean_by_formation_colony$sd.oh) #0.0139018 0.2877568

mean(mean_by_formation_colony$sd.o.side) #0.08611018
range(mean_by_formation_colony$sd.o.side) #0.01430516 0.58861649

mean(mean_by_formation_colony$sd.c.side) #0.132877
range(mean_by_formation_colony$sd.c.side) #0.0353629 0.4793441

#means of means
mean_by_formation = mean_by_formation_colony %>%
  dplyr::group_by(formation) %>%
  dplyr::summarize(num.col = length(unique(colony.id)),
            num.zooid = sum(n.zooid),
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
write.csv(mean_by_formation,
          "./Results/n.for.formations.csv",
          row.names = FALSE)

colony_means = dat_lg_N %>%
  dplyr::group_by(colony.id) %>%
  dplyr::summarize(formation = formation[1],
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
  dplyr::summarize(avg.zh = mean(ln.zh, na.rm = T),
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
#5 24

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

save(data.list, file = "./Results/g_matrices_data_form_reg.RData")

load(file="./Results/g_matrices_data_form_reg.RData") #load the g matrices calculated above 
model_G <- data.list[[1]]
dat_lg_N <- data.list[[2]]
form_data <- data.list[[3]]
mean_by_formation_colony <- data.list[[4]]

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

# diagonals of p.cov to priors
## chekcing NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

# PRIORS
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

#RETRIEVE THE P MATRIX FROM THE MCMC OBJECT
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

ggsave(G_PC_dist, file = "./Results/G_PC_dist_form_reg.png", 
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

#Gmat
#Pmat

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
summary(lm(diag(Pmat[[1]]) ~ diag(Gmat[[1]]))) #p-value: 0.01675 *

plot(diag(Gmat[[2]]), diag(Pmat[[2]]),
     pch = 19, col = col.form[2],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKBS")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[2]]) ~ diag(Gmat[[2]]))) #p-value: 0.004617 **

plot(diag(Gmat[[3]]), diag(Pmat[[3]]),
     pch = 19, col = col.form[3],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tewkesbury",
     xlim = c(0, 0.01),
     ylim = c(0, 0.05))
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[3]]) ~ diag(Gmat[[3]]))) #p-value: 0.0003241 ***

plot(diag(Gmat[[4]]), diag(Pmat[[4]]),
     pch = 19, col = col.form[4],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Waipuru")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[4]]) ~ diag(Gmat[[4]]))) #p-value: 0.01508 *

plot(diag(Gmat[[5]]), diag(Pmat[[5]]),
     pch = 19, col = col.form[5],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonals",
     main = "Upper Kai-Iwi")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[5]]) ~ diag(Gmat[[5]]))) #p-value: 0.001084 **

plot(diag(Gmat[[6]]), diag(Pmat[[6]]),
     pch = 19, col = col.form[6],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tainui")
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[6]]) ~ diag(Gmat[[6]]))) #p-value: 0.0003324 ***

plot(diag(Gmat[[7]]), diag(Pmat[[7]]),
     pch = 19, col = col.form[7],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1) # looks good!
summary(lm(diag(Pmat[[7]]) ~ diag(Gmat[[7]]))) #p-value: 0.01981 *

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
summary(lm(p.eig_variances[[1]] ~ g.eig_variances[[1]])) #p-value: 2.727e-08 ***

plot(g.eig_variances[[2]], p.eig_variances[[2]],
     pch = 19, col = col.form[2],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "NKBS")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[2]] ~ g.eig_variances[[2]])) #p-value: 1.287e-08 ***

plot(g.eig_variances[[3]], p.eig_variances[[3]],
     pch = 19, col = col.form[3],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tewkesbury")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[3]] ~ g.eig_variances[[3]])) #p-value: 1.477e-08 ***

plot(g.eig_variances[[4]], p.eig_variances[[4]],
     pch = 19, col = col.form[4],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Waipuru")
abline(0, 1) # looks okay...
summary(lm(p.eig_variances[[4]] ~ g.eig_variances[[4]])) #p-value: 4.228e-06 ***

plot(g.eig_variances[[5]], p.eig_variances[[5]],
     pch = 19, col = col.form[5],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Upper Kai-Iwi")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[5]] ~ g.eig_variances[[5]])) #p-value: 1.639e-09 ***

plot(g.eig_variances[[6]], p.eig_variances[[6]],
     pch = 19, col = col.form[6],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "Tainui")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[6]] ~ g.eig_variances[[6]])) #p-value: 1.461e-07 ***

plot(g.eig_variances[[7]], p.eig_variances[[7]],
     pch = 19, col = col.form[7],
     xlab = "G standardized variances",
     ylab = "P standardized variances",
     main = "SHCSBSB")
abline(0, 1) # looks good!
summary(lm(p.eig_variances[[7]] ~ g.eig_variances[[7]])) #p-value: 4.839e-08 ***

##### PLOT P, E, AND G VARIATION ----
# E = units
# P = cov or units + colony id
# G = estimate or colony id

g.eig_vectors <- lapply(Gmat, function (x) {eigen(x)$vectors})
paste("Eigenvalue vectors")
head(g.eig_vectors)

formations <- c("NKLS", "NKBS", "Tewkesbury", 
                "Waipuru", "Upper Kai-Iwi", 
                "Tainui", "SHCSBSB")

r.names <- c(rep(formations[1], 8),
             rep(formations[2], 8),
             rep(formations[3], 8),
             rep(formations[4], 8),
             rep(formations[5], 8),
             rep(formations[6], 8),
             rep(formations[7], 8))

g.eig_vect_mat = do.call(rbind, g.eig_vectors)
g.eig_vect_mat = data.frame(g.eig_vect_mat,
                            rownames = r.names)
g.eig_vect = melt(g.eig_vect_mat)
ggplot(g.eig_vect,
       aes(x = variable, y = value,
           group = rownames,
           colour = rownames)) +
  geom_line(aes(linetype = rownames)) +
  geom_point() +
  xlab("Principal component vector") +
  ylab("%Variation in the PC")

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

#### COMPARE Gs ----

# Function to convert vector to norm length 
f.normalize_vector <- function(vector) {
  norm_length <- sqrt(sum(vector^2))
  normalized_vector <- vector / norm_length
  return(normalized_vector)
}

#trait means by time
mean_by_formation
#order of formations:
# "NKLS"  "NKBS" "Tewkesbury"  "Waipuru" "Upper Kai-Iwi" "Tainui" "SHCSBSB"

### Calculate the vector that defines the observed divergence between sample/formation 1 an 2

NKLS <- as.numeric(mean_by_formation[1, 5:12]) # A vector containing trait means from sample/formation 1 
NKBS <- as.numeric(mean_by_formation[2, 5:12]) # A vector containing trait means from sample/formation 2 
tewk <- as.numeric(mean_by_formation[3, 5:12]) # A vector containing trait means from sample/formation 3
wai <- as.numeric(mean_by_formation[4, 5:12]) # A vector containing trait means from sample/formation 4
uki <- as.numeric(mean_by_formation[5, 5:12]) # A vector containing trait means from sample/formation 5 
tai <- as.numeric(mean_by_formation[6, 5:12]) # A vector containing trait means from sample/formation 6 
SHCSBSB <- as.numeric(mean_by_formation[7, 5:12]) # A vector containing trait means from sample/formation 7

#second - first
#t1 = NKBS - NKLS
#t2 = tewk - NKBS
#t3 = wai - tewk
#t4 = uki - wai
#t5 = tai - uki
#t6 = SHCSBSB - tai

## really need to learn how to name things in functions...
evolved_difference_unit_length_t1 <- f.normalize_vector(NKBS - NKLS)
evolved_difference_unit_length_t2 <- f.normalize_vector(tewk - NKBS)
evolved_difference_unit_length_t3 <- f.normalize_vector(wai - tewk)
evolved_difference_unit_length_t4 <- f.normalize_vector(uki - wai)
evolved_difference_unit_length_t5 <- f.normalize_vector(tai - uki)
evolved_difference_unit_length_t6 <- f.normalize_vector(SHCSBSB - tai)

G_matrix_NKLS = Gmat[[1]] # The G matrix estimated for sample/formation 1
G_matrix_NKBS = Gmat[[2]] # The G matrix estimated for sample/formation 2
G_matrix_tewk = Gmat[[3]] # The G matrix estimated for sample/formation 3
G_matrix_wai = Gmat[[4]] # The G matrix estimated for sample/formation 4
G_matrix_uki = Gmat[[5]] # The G matrix estimated for sample/formation 5
G_matrix_tai = Gmat[[6]] # The G matrix estimated for sample/formation 6
G_matrix_SHCSBSB = Gmat[[7]] # The G matrix estimated for sample/formation 7

### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
observed_evolvability_in_direction_of_change_t1 <- t(evolved_difference_unit_length_t1)%*%as.matrix(G_matrix_NKLS)%*%evolved_difference_unit_length_t1
observed_evolvability_in_direction_of_change_t2 <- t(evolved_difference_unit_length_t2)%*%as.matrix(G_matrix_NKBS)%*%evolved_difference_unit_length_t2
observed_evolvability_in_direction_of_change_t3 <- t(evolved_difference_unit_length_t3)%*%as.matrix(G_matrix_tewk)%*%evolved_difference_unit_length_t3
observed_evolvability_in_direction_of_change_t4 <- t(evolved_difference_unit_length_t4)%*%as.matrix(G_matrix_wai)%*%evolved_difference_unit_length_t4
observed_evolvability_in_direction_of_change_t5 <- t(evolved_difference_unit_length_t5)%*%as.matrix(G_matrix_uki)%*%evolved_difference_unit_length_t5
observed_evolvability_in_direction_of_change_t6 <- t(evolved_difference_unit_length_t6)%*%as.matrix(G_matrix_tai)%*%evolved_difference_unit_length_t6

### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
observed_conditional_evolvability_in_direction_of_change_t1 <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(G_matrix_NKLS))%*%evolved_difference_unit_length_t1)
observed_conditional_evolvability_in_direction_of_change_t2 <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(G_matrix_NKBS))%*%evolved_difference_unit_length_t2)
observed_conditional_evolvability_in_direction_of_change_t3 <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(G_matrix_tewk))%*%evolved_difference_unit_length_t3)
observed_conditional_evolvability_in_direction_of_change_t4 <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(G_matrix_wai))%*%evolved_difference_unit_length_t4)
observed_conditional_evolvability_in_direction_of_change_t5 <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(G_matrix_uki))%*%evolved_difference_unit_length_t5)
observed_conditional_evolvability_in_direction_of_change_t6 <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(G_matrix_tai))%*%evolved_difference_unit_length_t6)

### Generate 10,000 selection gradients in random directions in the n-dimensional space
n_dimensions <- 8 # number of traits in G matrix
Beta <- randomBeta(10000, n_dimensions)

## check positive definite
# round to 10 decimals to make it symmetric
# to make it positive definite, use extended matrix to fill in
is.symmetric.matrix(Beta)
is.positive.definite(Beta)

G_mat_NKLS <- round(as.matrix(G_matrix_NKLS), 10)
is.symmetric.matrix(G_mat_NKLS) #TRUE
is.positive.definite(G_mat_NKLS) #FALSE

G_mat_NKBS <- round(as.matrix(G_matrix_NKBS), 10)
is.symmetric.matrix(G_mat_NKBS) #TRUE
is.positive.definite(G_mat_NKBS) #FALSE

G_mat_tewk <- round(as.matrix(G_matrix_tewk), 10)
is.symmetric.matrix(G_mat_tewk) #TRUE
is.positive.definite(G_mat_tewk) #FALSE

G_mat_wai <- round(as.matrix(G_matrix_wai), 10)
is.symmetric.matrix(G_mat_wai) #TRUE
is.positive.definite(G_mat_wai) #FALSE

G_mat_uki <- round(as.matrix(G_matrix_uki), 10)
is.symmetric.matrix(G_mat_uki) #TRUE
is.positive.definite(G_mat_uki) #FALSE

G_mat_tai <- round(as.matrix(G_matrix_tai), 10)
is.symmetric.matrix(G_mat_tai) #TRUE
is.positive.definite(G_mat_tai) #FALSE

G_mat_SHCSBSB <- round(as.matrix(G_matrix_SHCSBSB), 10)
is.symmetric.matrix(G_mat_SHCSBSB) #TRUE
is.positive.definite(G_mat_SHCSBSB) #FALSE

#outputs e, r, c, a, i
#e = evolvability
#r = respondability
#c = conditional evolvability
#a = autonomy of each selection gradient
#i = integration
#Beta = matrix of selection gradients
#e and c are calculating variances of means; should not be negative
#conditional must be equal to or smaller than e; often much small

# Compute the mean, minimum and maximum evolvability (e_mean, e_min, e_max) for a G matrix based on 10,000 random selection gradients
X_t1 <- evolvabilityBeta(as.matrix(G_matrix_NKLS), Beta)
sumX_t1 <- summary(X_t1) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t2 <- evolvabilityBeta(as.matrix(G_matrix_NKBS), Beta)
sumX_t2 <- summary(X_t2) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t3 <- evolvabilityBeta(as.matrix(G_matrix_tewk), Beta)
sumX_t3 <- summary(X_t3) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t4 <- evolvabilityBeta(as.matrix(G_matrix_wai), Beta)
sumX_t4 <- summary(X_t4) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t5 <- evolvabilityBeta(as.matrix(G_matrix_uki), Beta)
sumX_t5 <- summary(X_t5) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t6 <- evolvabilityBeta(as.matrix(G_matrix_tai), Beta)
sumX_t6 <- summary(X_t6) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t7 <- evolvabilityBeta(as.matrix(G_matrix_SHCSBSB), Beta)
sumX_t7 <- summary(X_t7) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_sum <- data.frame(c.mean = c(sumX_t1$Averages[[3]], sumX_t2$Averages[[3]], sumX_t3$Averages[[3]],
                               sumX_t4$Averages[[3]], sumX_t5$Averages[[3]], sumX_t6$Averages[[3]], 
                               sumX_t7$Averages[[3]]),
                    c.min = c(sumX_t1$Minimum[[3]], sumX_t2$Minimum[[3]], sumX_t3$Minimum[[3]],
                              sumX_t4$Minimum[[3]], sumX_t5$Minimum[[3]], sumX_t6$Minimum[[3]],
                              sumX_t7$Minimum[[3]]),
                    c.max = c(sumX_t1$Maximum[[3]], sumX_t2$Maximum[[3]], sumX_t3$Maximum[[3]],
                              sumX_t4$Maximum[[3]], sumX_t5$Maximum[[3]], sumX_t6$Maximum[[3]],
                              sumX_t7$Maximum[[3]]),
                    e.mean = c(sumX_t1$Averages[[1]], sumX_t2$Averages[[1]], sumX_t3$Averages[[1]],
                               sumX_t4$Averages[[1]], sumX_t5$Averages[[1]], sumX_t6$Averages[[1]],
                               sumX_t7$Averages[[1]]),
                    e.min = c(sumX_t1$Minimum[[1]], sumX_t2$Minimum[[1]], sumX_t3$Minimum[[1]],
                              sumX_t4$Minimum[[1]], sumX_t5$Minimum[[1]], sumX_t6$Minimum[[1]],
                              sumX_t7$Minimum[[1]]),
                    e.max = c(sumX_t1$Maximum[[1]], sumX_t2$Maximum[[1]], sumX_t3$Maximum[[1]],
                              sumX_t4$Maximum[[1]], sumX_t5$Maximum[[1]], sumX_t6$Maximum[[1]],
                              sumX_t7$Maximum[[1]]),
                    observed_e = c(observed_evolvability_in_direction_of_change_t1,
                                   observed_evolvability_in_direction_of_change_t2,
                                   observed_evolvability_in_direction_of_change_t3,
                                   observed_evolvability_in_direction_of_change_t4,
                                   observed_evolvability_in_direction_of_change_t5,
                                   observed_evolvability_in_direction_of_change_t6,
                                   ""),
                    observed_c = c(observed_conditional_evolvability_in_direction_of_change_t1,
                                   observed_conditional_evolvability_in_direction_of_change_t2,
                                   observed_conditional_evolvability_in_direction_of_change_t3,
                                   observed_conditional_evolvability_in_direction_of_change_t4,
                                   observed_conditional_evolvability_in_direction_of_change_t5,
                                   observed_conditional_evolvability_in_direction_of_change_t6,
                                   ""),
                    row.names = formation_list)

write.csv(X_sum,
          "./Results/evolvability_summary.csv")

# By comparing the evolvabilities you estimated in the direction of change (lines 9 and 12) with the average evolvabilities calculated by running line 20, you get a sense of whether evolution happened in directions with above or below average evolvability.  

### Proportion of variance in n-dimensional trait space that is explained by PC1 (i.e., the first eigenvector)
#eigen(as.matrix(G_matrix_1))$values[1]/sum(eigen(as.matrix(G_matrix_1))$values)
eigen(as.matrix(G_matrix_NKLS))$values[1]/sum(eigen(as.matrix(G_matrix_NKLS))$values) #0.4676502
eigen(as.matrix(G_matrix_NKBS))$values[1]/sum(eigen(as.matrix(G_matrix_NKBS))$values) #0.5220323
eigen(as.matrix(G_matrix_tewk))$values[1]/sum(eigen(as.matrix(G_matrix_tewk))$values) #0.5024787
eigen(as.matrix(G_matrix_wai))$values[1]/sum(eigen(as.matrix(G_matrix_wai))$values) #0.400385
eigen(as.matrix(G_matrix_uki))$values[1]/sum(eigen(as.matrix(G_matrix_uki))$values) #0.7530483
eigen(as.matrix(G_matrix_tai))$values[1]/sum(eigen(as.matrix(G_matrix_tai))$values) #0.4238095
eigen(as.matrix(G_matrix_SHCSBSB))$values[1]/sum(eigen(as.matrix(G_matrix_SHCSBSB))$values) #0.5006772

### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 
Gmax_NKLS <- eigen(G_matrix_NKLS)$vectors[,1]
Gmax_NKBS <- eigen(G_matrix_NKBS)$vectors[,1]
Gmax_tewk <- eigen(G_matrix_tewk)$vectors[,1]
Gmax_wai <- eigen(G_matrix_wai)$vectors[,1]
Gmax_uki <- eigen(G_matrix_uki)$vectors[,1]
Gmax_tai <- eigen(G_matrix_tai)$vectors[,1]
Gmax_SHCSBSB <- eigen(G_matrix_SHCSBSB)$vectors[,1]

# Put Gmax to norm length
Gmax_NKLS_norm <- f.normalize_vector(Gmax_NKLS)
Gmax_NKBS_norm <- f.normalize_vector(Gmax_NKBS)
Gmax_tewk_norm <- f.normalize_vector(Gmax_tewk)
Gmax_wai_norm <- f.normalize_vector(Gmax_wai)
Gmax_uki_norm <- f.normalize_vector(Gmax_uki)
Gmax_tai_norm <- f.normalize_vector(Gmax_tai)
Gmax_SHCSBSB_norm <- f.normalize_vector(Gmax_SHCSBSB)

##Compute angles
# Calculate the dot product of the unit vectors
dot_product.Gmax_NKLS_NKBS <- sum(Gmax_NKLS_norm * Gmax_NKBS_norm)
# Calculate the angle in radians
angle_radians.Gmax_NKLS_NKBS <- acos(dot_product.Gmax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_NKBS <- angle_radians.Gmax_NKLS_NKBS * (180 / pi)
#5.25

dot_product.Gmax_NKBS_tewk <- sum(Gmax_NKBS_norm * Gmax_tewk_norm)
# Calculate the angle in radians
angle_radians.Gmax_NKBS_tewk <- acos(dot_product.Gmax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_tewk <- angle_radians.Gmax_NKBS_tewk * (180 / pi)
#5.47

dot_product.Gmax_tewk_wai <- sum(Gmax_tewk_norm * Gmax_wai_norm)
# Calculate the angle in radians
angle_radians.Gmax_tewk_wai <- acos(dot_product.Gmax_tewk_wai)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_wai <- angle_radians.Gmax_tewk_wai * (180 / pi)
#23.74

dot_product.Gmax_wai_uki <- sum(Gmax_wai_norm * Gmax_uki_norm)
# Calculate the angle in radians
angle_radians.Gmax_wai_uki <- acos(dot_product.Gmax_wai_uki)
# Convert the angle to degrees
angle_degrees.Gmax_wai_uki <- angle_radians.Gmax_wai_uki * (180 / pi)
#27.52

dot_product.Gmax_uki_tai <- sum(Gmax_uki_norm * Gmax_tai_norm)
# Calculate the angle in radians
angle_radians.Gmax_uki_tai <- acos(dot_product.Gmax_uki_tai)
# Convert the angle to degrees
angle_degrees.Gmax_uki_tai <- angle_radians.Gmax_uki_tai * (180 / pi)
#32.69

dot_product.Gmax_tai_SHCSBSB <- sum(Gmax_tai_norm * Gmax_SHCSBSB_norm)
# Calculate the angle in radians
angle_radians.Gmax_tai_SHCSBSB <- acos(dot_product.Gmax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Gmax_tai_SHCSBSB <- angle_radians.Gmax_tai_SHCSBSB * (180 / pi)
#23.03

### See if change is in direction of G max
# Calculate the dot product of the unit vectors
dot_product.Gmax_NKLS_max <- sum(Gmax_NKLS_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Gmax_NKLS_max <- acos(dot_product.Gmax_NKLS_max)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_max <- angle_radians.Gmax_NKLS_max * (180 / pi)
#73.14529

# Calculate the dot product of the unit vectors
dot_product.Gmax_NKBS_max <- sum(Gmax_NKBS_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Gmax_NKBS_max <- acos(dot_product.Gmax_NKBS_max)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_max <- angle_radians.Gmax_NKBS_max * (180 / pi)
#36.81499

# Calculate the dot product of the unit vectors
dot_product.Gmax_tewk_max <- sum(Gmax_tewk_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Gmax_tewk_max <- acos(dot_product.Gmax_tewk_max)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_max <- angle_radians.Gmax_tewk_max * (180 / pi)
#112.5576

# Calculate the dot product of the unit vectors
dot_product.Gmax_wai_max <- sum(Gmax_wai_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Gmax_wai_max <- acos(dot_product.Gmax_wai_max)
# Convert the angle to degrees
angle_degrees.Gmax_wai_max <- angle_radians.Gmax_wai_max * (180 / pi)
#153.9498

# Calculate the dot product of the unit vectors
dot_product.Gmax_uki_max <- sum(Gmax_uki_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Gmax_uki_max <- acos(dot_product.Gmax_uki_max)
# Convert the angle to degrees
angle_degrees.Gmax_uki_max <- angle_radians.Gmax_uki_max * (180 / pi)
#156.753

# Calculate the dot product of the unit vectors
dot_product.Gmax_tai_max <- sum(Gmax_tai_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Gmax_tai_max <- acos(dot_product.Gmax_tai_max)
# Convert the angle to degrees
angle_degrees.Gmax_tai_max <- angle_radians.Gmax_tai_max * (180 / pi)
#53.85211

#### GLOBAL G ----

##### PRIORS -----
#save(dat_lg_N, file = "./Results/dat_lg_N.RData")

#dat_lg_N.com = dat_lg_N[complete.cases(dat_lg_N),] #didn't fix anything
phen.var.glob = cov(dat_lg_N[, 4:11]) #traits of ALL; correct for colony and formation later
prior.glob = list(G = list(G1 = list(V = phen.var.glob/2, nu = 10), #V same as individual G matrices; nu is different
                           G2 = list(V = phen.var.glob/2, nu = 10)), #additional V, made same as above
                  R = list(V = phen.var.glob/4, nu = 5)) #V same as individual G matrices

##### MCMC -----
#Running the MCMC chain

model_Global <- MCMCglmm(cbind(ln.zh, ln.mpw.b, ln.cw.m, ln.cw.d, #same order as in priors
                               ln.ow.m, ln.oh, ln.c.side, ln.o.side) ~ trait-1,
                         #account for variation w/in colony:
                         random = ~us(trait):colony.id + us(trait):formation, #the number of these determines # of Gs #+ us(trait):formation
                         rcov = ~us(trait):units,
                         family = rep("gaussian", 8), #num of traits
                         data = dat_lg_N,
                         nitt = 1500000, thin = 1000, burnin = 500000,
                         prior = prior.glob, verbose = TRUE)

save(model_Global, file = "./Results/global_matrices_data_form_reg.RData")

#load(file="./Results/global_matrices_data_form_reg.RData") #load the g matrices calculated above 
#model_Global <- data.list[[1]]

##### CHECK MODELS -----
formation_list #order of formations
summary(model_Global)

##plots to see where sampling from:
plot(model_Global$VCV) #catepillar!

###### POSTERIOR G MATRIX ------
#Retrieving G from posterior
glob.model = model_Global
ntraits = 8
glob.Gmat = matrix(posterior.mode(glob.model$VCV)[1:ntraits^2], ntraits, ntraits)

##### G VARIANCES -----
isSymmetric(glob.Gmat)  #is.symmetric.matrix
glob.variances = diag(glob.Gmat)
paste("Trait variances")
head(glob.variances)

###### G EIGEN ------
glob.eig_variances = eigen(glob.Gmat)$values
paste("Eigenvalue variances")
head(glob.eig_variances)

glob.eig_percent = glob.eig_variances/sum(glob.eig_variances) 
glob.eig_per_mat = data.frame(glob.eig_percent)
#glob.eig_per = melt(glob.eig_per_mat)
glob.eig_per_mat$PC <- 1:nrow(glob.eig_per_mat)
#dev.off()
Glob_PC_dist = ggplot(glob.eig_per_mat,
                   aes(x = PC, y = glob.eig_percent)) +
  geom_line() +
  geom_point() +
  xlab("Principal component rank") +
  ylab("%Variation in the PC")
Glob_PC_dist #one negative; none above 1!

ggsave(Glob_PC_dist, file = "./Results/Global_G_PC_dist_form_reg.png", 
       width = 14, height = 10, units = "cm")

#Note that some matrices have negative eigenvalues. 
#This can cause a lot of trouble in analyses involving inverted matrices.
#Solution from evolqg Marroig et al. 2012

###### G NOISE ------
##Controlling for noise
#Extend G
Glob_ext = ExtendMatrix(glob.Gmat, ret.dim = 6)$ExtMat #not 8 because last eigen value (#8) was negative
#ignore warning from above
isSymmetric(Glob_ext)
Glob_Ext_std_variances = diag(Glob_ext) 
Glob_Ext_eig_variances = eigen(Glob_ext)$values

#### COMPARE TO GLOBAL G ----

### Calculate the vector that defines the observed divergence between sample/formation 1 an 2
global <- as.numeric(means)
#NKLS
#NKBS
#tewk
#wai
#uki
#tai
#SHCSBSB

evolved_difference_unit_length_NKLS_glob <- f.normalize_vector(NKLS - global)
evolved_difference_unit_length_NKBS_glob <- f.normalize_vector(NKBS - global)
evolved_difference_unit_length_tewk_glob <- f.normalize_vector(tewk - global)
evolved_difference_unit_length_wai_glob <- f.normalize_vector(wai - global)
evolved_difference_unit_length_uki_glob <- f.normalize_vector(uki - global)
evolved_difference_unit_length_tai_glob <- f.normalize_vector(tai - global)
evolved_difference_unit_length_SCHSBSB_glob <- f.normalize_vector(SHCSBSB-global)

G_matrix_glob = glob.Gmat 
#G_matrix_NKLS
#G_matrix_NKBS
#G_matrix_tewk
#G_matrix_wai
#G_matrix_uki
#G_matrix_tai
#G_matrix_SHCSBSB

### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
observed_evolvability_in_direction_of_change_NKLS <- t(evolved_difference_unit_length_NKLS_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_NKLS_glob
observed_evolvability_in_direction_of_change_NKBS <- t(evolved_difference_unit_length_NKBS_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_NKBS_glob
observed_evolvability_in_direction_of_change_tewk <- t(evolved_difference_unit_length_tewk_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_tewk_glob
observed_evolvability_in_direction_of_change_wai <- t(evolved_difference_unit_length_wai_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_wai_glob
observed_evolvability_in_direction_of_change_uki <- t(evolved_difference_unit_length_uki_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_uki_glob
observed_evolvability_in_direction_of_change_tai <- t(evolved_difference_unit_length_tai_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_tai_glob
observed_evolvability_in_direction_of_change_SHCSBSB <- t(evolved_difference_unit_length_SCHSBSB_glob)%*%as.matrix(G_matrix_glob)%*%evolved_difference_unit_length_SCHSBSB_glob

### The conditional evolvability in the direction of divergence
observed_conditional_evolvability_in_direction_of_change_NKLS <- 1/(t(evolved_difference_unit_length_NKLS_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_NKLS_glob)
observed_conditional_evolvability_in_direction_of_change_NKBS <- 1/(t(evolved_difference_unit_length_NKBS_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_NKBS_glob)
observed_conditional_evolvability_in_direction_of_change_tewk <- 1/(t(evolved_difference_unit_length_tewk_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_tewk_glob)
observed_conditional_evolvability_in_direction_of_change_wai <- 1/(t(evolved_difference_unit_length_wai_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_wai_glob)
observed_conditional_evolvability_in_direction_of_change_uki <- 1/(t(evolved_difference_unit_length_uki_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_uki_glob)
observed_conditional_evolvability_in_direction_of_change_tai <- 1/(t(evolved_difference_unit_length_tai_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_tai_glob)
observed_conditional_evolvability_in_direction_of_change_SCHSBSB <- 1/(t(evolved_difference_unit_length_SCHSBSB_glob)%*%solve(as.matrix(G_matrix_glob))%*%evolved_difference_unit_length_SCHSBSB_glob)

### Generate 10,000 selection gradients in random directions in the n-dimensional space
#n_dimensions <- 8 # number of traits in G matrix
#Beta <- randomBeta(10000, n_dimensions)

# Compute the mean, minimum and maximum evolvability (e_mean, e_min, e_max) for a G matrix based on 10,000 random selection gradients
X_glob <- evolvabilityBeta(as.matrix(G_matrix_glob), Beta)
summary(X_glob) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix
#e_mean 0.00664551
#c_mean -0.00010400
#e_min 1.237089e-04
#c_min -2.428818e-01
#e_max 0.02531620
#c_max 0.33979277

#X_t1 
#X_t2
#X_t3 
#X_t4 
#X_t5 
#X_t6 
#X_t7 

# By comparing the evolvabilities you estimated in the direction of change (lines 9 and 12) with the average evolvabilities calculated by running line 20, you get a sense of whether evolution happened in directions with above or below average evolvability.  

### Proportion of variance in n-dimensional trait space that is explained by PC1 (i.e., the first eigenvector)
eigen(as.matrix(G_matrix_glob))$values[1]/sum(eigen(as.matrix(G_matrix_glob))$values)
#0.5229979

### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 

Gmax_glob <- eigen(G_matrix_glob)$vectors[,1]

#Gmax_NKLS
#Gmax_NKBS
#Gmax_tewk
#Gmax_wai
#Gmax_uki
#Gmax_tai
#Gmax_SCHSBSB

# Put Gmax to norm length
Gmax_glob_norm <- f.normalize_vector(Gmax_glob)

#Gmax_NKLS_norm
#Gmax_NKBS_norm
#Gmax_tewk_norm
#Gmax_wai_norm
#Gmax_uki_norm
#Gmax_tai_norm
#Gmax_SCHSBSB_norm

##Compute angles

## NKLS
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_NKLS <- sum(Gmax_glob_norm * Gmax_NKLS_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_NKLS <- acos(dot_product.Gmax_glob_NKLS)
# Convert the angle to degrees
angle_degrees.Gmax_glob_NKLS <- angle_radians.Gmax_glob_NKLS * (180 / pi)
#7.953553

## NKBS
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_NKBS <- sum(Gmax_glob_norm * Gmax_NKBS_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_NKBS <- acos(dot_product.Gmax_glob_NKBS)
# Convert the angle to degrees
angle_degrees.Gmax_glob_NKBS <- angle_radians.Gmax_glob_NKBS * (180 / pi)
#3.790778

## Tewksbury
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_tewk <- sum(Gmax_glob_norm * Gmax_tewk_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_tewk <- acos(dot_product.Gmax_glob_tewk)
# Convert the angle to degrees
angle_degrees.Gmax_glob_tewk <- angle_radians.Gmax_glob_tewk * (180 / pi)
#5.831842

## Waipuru
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_wai <- sum(Gmax_glob_norm * Gmax_wai_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_wai <- acos(dot_product.Gmax_glob_wai)
# Convert the angle to degrees
angle_degrees.Gmax_glob_wai <- angle_radians.Gmax_glob_wai * (180 / pi)
#20.70651

## Upper Kai-Iwi
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_uki <- sum(Gmax_glob_norm * Gmax_uki_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_uki <- acos(dot_product.Gmax_glob_uki)
# Convert the angle to degrees
angle_degrees.Gmax_glob_uki <- angle_radians.Gmax_glob_uki * (180 / pi)
#17.58378

## Tainui
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_tai <- sum(Gmax_glob_norm * Gmax_tai_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_tai <- acos(dot_product.Gmax_glob_tai)
# Convert the angle to degrees
angle_degrees.Gmax_glob_tai <- angle_radians.Gmax_glob_tai * (180 / pi)
#22.45352

## SHCSBSB
# Calculate the dot product of the unit vectors
dot_product.Gmax_glob_SHCSBSB <- sum(Gmax_glob_norm * Gmax_SHCSBSB_norm)
# Calculate the angle in radians
angle_radians.Gmax_glob_SHCSBSB <- acos(dot_product.Gmax_glob_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Gmax_glob_SHCSBSB <- angle_radians.Gmax_glob_SHCSBSB * (180 / pi)
#10.81873

#### NOW INCLUDE SM ZOOIDS ----

df.sm.trim <- sm.df %>%
  dplyr::select(zooid.id, colony.id, formation, matches(traits))

sm.df = as.data.frame(df.sm.trim)

##### COMBINE -----

all.df <- rbind(df, sm.df)

##### REMAKE THINGS -----
zooid_list.all <- unique(all.df$zooid.id)
length(zooid_list.all) #5971

colony_list.all <- unique(all.df$colony.id)
length(colony_list.all) #572

# arrange formations from oldest to youngest
all.df$formation <- factor(all.df$formation, 
                           levels = c("NKLS", "NKBS", "Tewkesbury", 
                                      "Waipuru", "Upper Kai-Iwi", 
                                      "Tainui", "SHCSBSB")) 
