# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 29 Aug 2023

# The purpose of this script is to create a P and G matrix for 
# Steginoporella magnifica for each age (formation) it is found
# Similar to Voje et al. where we are looking at formations rather than species

# These matrices are initially made based on the following linear measurements:
# zh = zooid height
# mpw.b = median process base width
# cw.m = cryptocyst mid-width
# cw.d = cryptocyst distal-width
# ow.m = operculum mid-width
# we added:
# c.side = cryptocyst side length
# o.side = operculum side length

# specimenNR = colony ID
# new.ID = individual zooid ID
# formation relates to the age

# Voje et al. used a 4 zooid/colony minimum; perhaps we should too

# the outputs are:
# - order of formations (formation_list saved in data.list)
# - order of traits (traits saved in data.list)
# - plot of distribution of all the traits (saved as plot ml)
# - results of normality tests for all the trait distributions (saved as csv p.vals.df)
# - a trimmed data set (df) and dataset for downstream analyses (dat_lg_N saved in data.list)
# - sample sizes of colonies per formation (col_form.n) and zooids per formation (by_form.n) saved in data.list
# - P matrices for each formation (P_ext saved in data.list)
# - MCMC runs to estimates G (model_G saved in data.list) and G matrices for each formation (G_ext saved in data.list)
# - correlation plots of P and G for each formation (saved as graphs p.corr_mat, g.corr_mat)
# - rarefaction plot (saved as plot)


#### LOAD DATA ----

source("./Scripts/0-env.R")

df <- read.csv("./Results/colonies.traits_8Dec2023.csv",
               header = TRUE, 
               sep = ",",
               stringsAsFactors = FALSE)
#already have small zooid removed and at least 5 zooids per colony
#output from exploratoryAnalysis.R

form.meta <- read.csv("~/Documents/GitHub/bryozoa/stegino_metadata/newMetadata/formations.csv",
                      header = TRUE,
                      sep = ",",
                      stringsAsFactors = FALSE)

load(file = "./Results/sum.data.list.w.modern.RData") #load the g matrices calculated above 
mean_by_formation <- sum.data.list[[1]]
mean_by_formation_colony <- sum.data.list[[2]]

#### MANIPULATE DATA ----

zooid_list <- unique(df$zooid.id)
length(zooid_list) #5694 (was 15773 without modern)

colony_list <- unique(df$colony.id)
length(colony_list) #558 (was 742 without modern)

# arrange formations from oldest to youngest
df$formation <- factor(df$formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                "Waipuru", "Upper Kai-Iwi", 
                                                "Tainui", "SHCSBSB", "modern")) 
formation_list <- unique(df$formation)
length(formation_list) #8

#same order as in df
names(df)
traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])
length(traits) #8

##### TRIM DATASET ----
df.trim <- df %>%
  dplyr::select(zooid.id, colony.id, formation, matches(traits))

colNums <- match(c(traits, "zooid.id"), names(df.trim))
#  4  5 6  7  8  9 10 11  1 (i.e., 4:11 are traits of interest)

df = as.data.frame(df.trim)

#### PLOT TRAITS ----

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

ggsave(ml, 
       file = "./Results/trait.distribution.w.modern.png", 
       width = 14, height = 10, units = "cm")

#### NORMALITY TESTS ----
## shapiro test; but need to subsample to be within 5000
## if significant, then significantly different from normal (i.e., non-normal)
p.vals <- c()

#ln.zh
sub.ln.zh <- sample(df[, traits[1]], 
                    5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.zh)
p.vals[1] <- shapiro.test(sub.ln.zh)$p.value

#ln.mpw.b
sub.ln.mpw.b <- sample(df[, traits[2]], 
                       5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.mpw.b)
p.vals[2] <- shapiro.test(sub.ln.mpw.b)$p.value

#ln.cw.m
sub.ln.cw.m <- sample(df[, traits[3]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.cw.m) 
p.vals[3] <- shapiro.test(sub.ln.cw.m)$p.value

#ln.cw.d
sub.ln.cw.d <- sample(df[, traits[4]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.cw.d) 
p.vals[4] <- shapiro.test(sub.ln.cw.d)$p.value

#ln.ow.m
sub.ln.ow.m <- sample(df[, traits[5]], 
                      5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.ow.m) 
p.vals[5] <- shapiro.test(sub.ln.ow.m)$p.value

#ln.oh
sub.ln.oh <- sample(df[, traits[6]], 
                    5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.oh) 
p.vals[6] <- shapiro.test(sub.ln.oh)$p.value

#ln.c.side
sub.ln.c.side <- sample(df[, traits[7]], 
                        5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.c.side)
p.vals[7] <- shapiro.test(sub.ln.c.side)$p.value

#ln.o.side
sub.ln.o.side <- sample(df[, traits[8]], 
                        5000, replace = FALSE, prob = NULL)
shapiro.test(sub.ln.o.side) 
p.vals[8] <- shapiro.test(sub.ln.o.side)$p.value

p.vals.df <- as.data.frame(cbind(traits, p.vals))

write.csv(p.vals.df,
          "./Results/normality.test.w.modern.csv",
          row.names = FALSE)

#by formation
df.norm <- df %>%
    dplyr::group_by(formation) %>%
    dplyr::summarise(p.val.zh = shapiro.test(ln.zh)$p.value,
                     p.val.mpw.b = shapiro.test(ln.mpw.b)$p.value,
                     p.val.cw.m = shapiro.test(ln.cw.m)$p.value,
                     p.val.cw.d = shapiro.test(ln.cw.d)$p.value,
                     p.val.ow.m = shapiro.test(ln.ow.m)$p.value,
                     p.val.oh = shapiro.test(ln.oh)$p.value,
                     p.val.c.side = shapiro.test(ln.c.side)$p.value,
                     p.val.o.side = shapiro.test(ln.o.side)$p.value) %>%
    as.data.frame()
write.csv(df.norm,
          "./Results/normality.test.by.formation.w.modern.csv",
          row.names = FALSE)

#### REDUCE TO TRAITS OF INTEREST ----
trt_lg_N = c("formation", "colony.id", "zooid.id", traits)
dat_lg_N = df[intersect(colnames(df), trt_lg_N)]
head(dat_lg_N) #traits in same order as df and traits

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

#make sampling matrix
form <- levels(formation_list)
col.samp <- c(col_form.n[[1]], col_form.n[[2]], col_form.n[[3]],
              col_form.n[[4]], col_form.n[[5]], col_form.n[[6]],
              col_form.n[[7]], col_form.n[[8]])
zoo.samp <- c(by_form.n[[1]], by_form.n[[2]], by_form.n[[3]],
              by_form.n[[4]], by_form.n[[5]], by_form.n[[6]],
              by_form.n[[7]], by_form.n[[8]])
samp <- cbind(form, col.samp, zoo.samp)
write.csv(samp,
          "./Results/sampling.per.formation.w.moder.csv",
          row.names = FALSE)

#### P MATRIX ----
p.cov = lapply(form_data, function (x){ (cov(x[, 4:11]))}) #traits per colony (not variation within colony)
#p.cov is the same and the phen.var

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
P_PC_dist #none negative; none above 1; dim 8 close to 0

ggsave(P_PC_dist, 
       file = "./Results/P.PC.dist.w.modern.png", 
       width = 14, height = 10, units = "cm") 

###### P NOISE ------
##Controlling for noise
#Extend G
P_ext = lapply(Pmat, function (x){ ExtendMatrix(x, ret.dim = 5)$ExtMat}) #to match the G matrix
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

save(model_G,
     file = "./Results/model_G.w.modern.RData")

load(file = "./Results/model_G.w.modern.RData") #load the g matrices calculated above 

##### CHECK MODELS -----
formation_list #order of formations
summary(model_G[[1]])
summary(model_G[[2]])
summary(model_G[[3]])
summary(model_G[[4]])
summary(model_G[[5]])
summary(model_G[[6]])
summary(model_G[[7]])
summary(model_G[[8]])

##plots to see where sampling from:
plot(model_G[[1]]$VCV) #catepillar!
plot(model_G[[2]]$VCV) #catepillar!
plot(model_G[[3]]$VCV) #catepillar!
plot(model_G[[4]]$VCV) #catepillar!; high skew
plot(model_G[[5]]$VCV) #catepillar!
plot(model_G[[6]]$VCV) #catepillar!; high skew
plot(model_G[[7]]$VCV) #catepillar!
plot(model_G[[8]]$VCV) #catepillar!
#formations from oldest to youngest: "NKLS", "NKBS", "Tewkesbury", "Waipuru", 
#                                    "Upper Kai-Iwi", "Tainui", "SHCSBSB", "modern"

##### CHECK P IS BIGGER THAN G -----

# diagonals of p.cov to priors
## chekcing NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

# PRIORS
# NKLS
diag(phen.var[[1]]) #all larger
diag(prior$NKBS$G$G1$V) < diag(phen.var[[1]]) #all larger
diag(prior$NKBS$R$V) < diag(phen.var[[1]]) #all larger

# NKBS
diag(phen.var[[2]]) #all larger
diag(prior$NKBS$G$G1$V) < diag(phen.var[[2]]) #all larger
diag(prior$NKBS$R$V) < diag(phen.var[[2]]) #all larger

# Tewkesbury
diag(phen.var[[3]]) #all larger
diag(prior$NKBS$G$G1$V) < diag(phen.var[[3]]) #all larger
diag(prior$NKBS$R$V) < diag(phen.var[[3]]) #all larger

# Waipuru
diag(phen.var[[4]]) #all larger
diag(prior$Waipuru$G$G1$V) < diag(phen.var[[4]])
diag(prior$Waipuru$R$V) < diag(phen.var[[4]])

# Upper Kai-Iwi
diag(phen.var[[5]]) #all larger
diag(prior$`Upper Kai-Iwi`$G$G1$V) < diag(phen.var[[5]])
diag(prior$`Upper Kai-Iwi`$R$V) < diag(phen.var[[5]])

# Tainui
diag(phen.var[[6]]) #all larger
diag(prior$`Upper Kai-Iwi`$G$G1$V) < diag(phen.var[[6]])
diag(prior$`Upper Kai-Iwi`$R$V) < diag(phen.var[[6]])

# SHCSBSB
diag(phen.var[[7]]) #all larger
diag(prior$`Upper Kai-Iwi`$G$G1$V) < diag(phen.var[[7]])
diag(prior$`Upper Kai-Iwi`$R$V) < diag(phen.var[[7]])

# modern
diag(phen.var[[8]]) #all larger
diag(prior$NKBS$G$G1$V) < diag(phen.var[[8]]) #all larger
diag(prior$NKBS$R$V) < diag(phen.var[[8]]) #all larger

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
d.col.vcv.nkbs < diag(Pmat[[2]]) # all larger

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
d.col.vcv.wp < diag(Pmat[[4]]) # cw.m; ow.m; oh; c.side; o.side are smaller than P matrix; all larger than G matrix 

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
d.col.vcv.uki < diag(Pmat[[5]]) # zh, mpw.b, cw.m, cw.d, ow.m, oh, c.side, o.side are smaller; all larger than G matrix 

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

diag(Gmat[[2]]) < diag(Pmat[[2]]) # all larger

diag(Gmat[[4]]) < diag(Pmat[[4]]) # all larger

diag(Gmat[[5]]) < diag(Pmat[[5]]) # all larger

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
G_PC_dist #Tainui negative; none above 1; Upper Kai-Iwi looks FUNKY
## Waipuru 5th dim is wild; may only ret 4 dim?

ggsave(G_PC_dist, 
       file = "./Results/G.PC.dist.w.modern.png", 
       width = 14, height = 10, units = "cm")

#Note that some matrices have negative eigenvalues. 
#This can cause a lot of trouble in analyses involving inverted matrices.
#Solution from evolqg Marroig et al. 2012

###### G NOISE ------
##Controlling for noise
#Extend G
G_ext = lapply(Gmat, function (x){ ExtendMatrix(x, ret.dim = 5)$ExtMat}) #not 8 because last eigen value (#6) was negative
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

#### SAVE ALL OUTPUTS ----
data.list = list(P_ext, G_ext,
                 df, dat_lg_N, 
                 mean_by_formation, mean_by_formation_colony, 
                 form_data, by_form.n, col_form.n,
                 means)
save(data.list,
     file = "./Results/data.list.w.modern.RData")

load(file = "./Results/data.list.w.modern.RData") #load the g matrices calculated above 


#### CHECK P > G ----

## PLOT DIAGONALS
plot(diag(Gmat[[1]]), diag(Pmat[[1]]),
     pch = 19, col = col.form[1],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKLS",
     xlim = c(0, .05),
     ylim = c(0, .2))
abline(0, 1) 
summary(lm(diag(Pmat[[1]]) ~ diag(Gmat[[1]]))) 
#slope = 2.4; r2 = 0.65; p is sig 0.009

plot(diag(Gmat[[2]]), diag(Pmat[[2]]),
     pch = 19, col = col.form[2],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKBS")
abline(0, 1)
summary(lm(diag(Pmat[[2]]) ~ diag(Gmat[[2]])))
#slope = 3.4; r2 = 0.8; p is sig 0.003

plot(diag(Gmat[[3]]), diag(Pmat[[3]]),
     pch = 19, col = col.form[3],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tewkesbury",
     xlim = c(0, 0.01),
     ylim = c(0, 0.05))
abline(0, 1) 
summary(lm(diag(Pmat[[3]]) ~ diag(Gmat[[3]])))
#slope = 3.5; r2=0.8; p is sig 0.002

plot(diag(Gmat[[4]]), diag(Pmat[[4]]),
     pch = 19, col = col.form[4],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Waipuru")
abline(0, 1)
summary(lm(diag(Pmat[[4]]) ~ diag(Gmat[[4]]))) 
#slope = 1.7; r2 = 0.7; p is sig 0.009

plot(diag(Gmat[[5]]), diag(Pmat[[5]]),
     pch = 19, col = col.form[5],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonals",
     main = "Upper Kai-Iwi")
abline(0, 1) 
summary(lm(diag(Pmat[[5]]) ~ diag(Gmat[[5]])))
#slope = 1.5; r2 = 0.96; p is sig <0.001

plot(diag(Gmat[[6]]), diag(Pmat[[6]]),
     pch = 19, col = col.form[6],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tainui")
abline(0, 1) 
summary(lm(diag(Pmat[[6]]) ~ diag(Gmat[[6]]))) 
#slope = 1.9; r2 = 0.8; p is sig 0.001

plot(diag(Gmat[[7]]), diag(Pmat[[7]]),
     pch = 19, col = col.form[7],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1) 
summary(lm(diag(Pmat[[7]]) ~ diag(Gmat[[7]]))) 
#slope = 2.4; r2 = 0.5; p = 0.026

plot(diag(Gmat[[8]]), diag(Pmat[[8]]),
     pch = 19, col = col.form[8],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1) 
summary(lm(diag(Pmat[[8]]) ~ diag(Gmat[[8]]))) 
#slope = 1.3; r2 = 0.8; p = 0.003

#### RAREFACTION ----
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

comp_sampleN = matrix(0, 8, 8) #calculating the smallest sample size in the comparison among pairs of G 

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
       col = "#00BFC4", pch = 19)

p.rare <- ggplot() +
    geom_point(aes(out_results[, 2], out_results[, 1]),
               pch = 19, col = "grey", size = 2) +
    geom_point(aes(obs_melt$N, obs_melt$RS),
               col = "#00BFC4", pch = 19, size = 2) + 
    #geom_point(aes(obs_melt$N[4], obs_melt$RS[4]), #using to check where the points are
    #           col = "red", pch = 19, size = 2)
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA)) +
    scale_x_continuous(name = "Sample size") +
    scale_y_continuous(name = "Similarity") 

ggsave(p.rare, 
       file = "./Results/rarefaction.w.modern.png", 
       width = 14, height = 10, units = "cm")

## find which ones are outside of the gray
#low sample size, smallest similarity
which.min(obs_melt$RS)
obs_melt[25,] #0.6630677
comp_mat$correlations #.663 is corr between modern and Upper Kai-Iwi
#low sample size, next lowest similarity
sort(obs_melt$RS, decreasing = TRUE)
obs_melt[obs_melt$RS <= 0.7407344,] #index 18
comp_mat$correlations #0.74 is corr between modern and Tewkesbury

## COLOR MODERN

#### GLOBAL G ----

##### PRIORS -----

#dat_lg_N.com = dat_lg_N[complete.cases(dat_lg_N),] #didn't fix anything
phen.var.glob = cov(dat_lg_N[, 4:11]) #traits of ALL; correct for colony and formation later
prior.glob = list(G = list(G1 = list(V = phen.var.glob/2, nu = 2)), #nu = 10 #V same as individual G matrices; nu is different
                  R = list(V = phen.var.glob/4, nu = 2)) #nu = 5 #V same as individual G matrices

##### MCMC -----
#Running the MCMC chain

model_Global <- MCMCglmm(cbind(ln.zh, ln.mpw.b, ln.cw.m, ln.cw.d, #same order as in priors
                               ln.ow.m, ln.oh, ln.c.side, ln.o.side) ~ trait + trait:formation - 1, #get rid of "trait" alone?
                         #account for variation w/in colony:
                         random = ~us(trait):colony.id, #the number of these determines # of Gs
                         rcov = ~us(trait):units,
                         family = rep("gaussian", 8), #num of traits
                         data = dat_lg_N, #need to have formation as factor in data matrix
                         nitt = 1500000, thin = 1000, burnin = 500000,
                         prior = prior.glob, verbose = TRUE)

save(model_Global, 
     file = "./Results/global_matrix.w.modern.RData")

load(file="./Results/global_matrix.w.modern.RData") #load the g matrices calculated above 
model_Global

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

ggsave(Glob_PC_dist, 
       file = "./Results/GlobalG.PC.dist.reg.png", 
       width = 14, height = 10, units = "cm")

#Note that some matrices have negative eigenvalues. 
#This can cause a lot of trouble in analyses involving inverted matrices.
#Solution from evolqg Marroig et al. 2012

###### G NOISE ------
##Controlling for noise
#Extend G
Glob_ext = ExtendMatrix(glob.Gmat, ret.dim = 5)$ExtMat #not 8 because last eigen value (#8) was negative
Glob_ext <- as.matrix(Glob_ext)
#ignore warning from above
isSymmetric(Glob_ext)
Glob_Ext_std_variances = diag(Glob_ext) 
Glob_Ext_eig_variances = eigen(Glob_ext)$values
##need to create random cov.m for comparison
cov.m <- RandomMatrix(8, 1, 1, 100) 
G_list <- list(Glob_ext, as.matrix(cov.m))

glob.comp_mat = RandomSkewers(G_list) #need at least
glob.corr_mat = glob.comp_mat$correlations + t(glob.comp_mat$correlations) 
diag(glob.corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(glob.corr_mat,upper = "number", lower = "pie")

save(Glob_ext, 
     file = "./Results/global_ext.w.modern.RData")

load(file="./Results/global_ext.w.modern.RData") #load the g matrices calculated above 
Glob_ext

