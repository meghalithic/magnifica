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

df <- read.csv("./Results/colonies.traits_1Jul2024.csv",
               header = TRUE, 
               sep = ",",
               stringsAsFactors = FALSE)
#already have small zooid removed and at least 5 zooids per colony
#output from exploratoryAnalysis.R

load(file = "./Results/sum.data.list.RData") #load the g matrices calculated above 
mean_by_formation <- sum.data.list[[1]]
mean_by_formation_colony <- sum.data.list[[2]]
means <- sum.data.list[[3]]

mean(mean_by_formation_colony$n.zooid)

#### MANIPULATE DATA ----

zooid_list <- unique(df$zooid.id)
length(zooid_list) #6264 (6658 for 3 zoo)

colony_list <- unique(df$colony.id)
length(colony_list) #599 (711 for 3 zoo)

# arrange formations from oldest to youngest
df$formation <- factor(df$formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                "Upper Kai-Iwi",  "Tainui", 
                                                "SHCSBSB", "modern")) 
formation_list <- unique(df$formation)
length(formation_list) #7

#same order as in df
names(df)
traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])
length(traits) #8

##### TRIM DATASET ----
df.trim <- df %>%
  dplyr::select(zooid.id, colony.id, locality, formation, matches(traits))

colNums <- match(c(traits, "zooid.id"), names(df.trim))
#  4  5 6  7  8  9 10 11  1 (i.e., 4:11 are traits of interest)

df = as.data.frame(df.trim)

#### PLOT TRAITS ----

p.zh = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[1]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[1]) +
  scale_color_manual(values = col.form)

p.mpw.b = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[2]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[2]) +
  scale_color_manual(values = col.form)

p.cw.m = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[3]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[3]) +
  scale_color_manual(values = col.form)

p.cw.d = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[4]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[4]) +
  scale_color_manual(values = col.form)

p.ow.m = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[5]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[5]) +
  scale_color_manual(values = col.form)

p.oh = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[6]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[6]) +
  scale_color_manual(values = col.form)

p.c.side = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[7]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[7]) +
  scale_color_manual(values = col.form)

p.o.side = ggplot(data = df) + 
  geom_density(aes(x = df[, traits[8]], 
                   group = formation,
                   col = formation)) + 
  plot.theme +
  scale_x_continuous(name = traits[8]) +
  scale_color_manual(values = col.form)

Fig = list(p.zh, p.mpw.b, p.cw.m, p.cw.d, p.ow.m, p.oh, p.c.side, p.o.side)
ml <- marrangeGrob(Fig, nrow = 4, ncol = 2)
ml
#see a really big shift between Tewkesbury and earlier and Upper Kai-Iwi and after

ggsave(ml, 
       file = "./Results/trait.distribution.png", 
       width = 14, height = 20, units = "cm")

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
          "./Results/normality.test.csv",
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
          "./Results/normality.test.by.formation.csv",
          row.names = FALSE)

#### REDUCE TO TRAITS OF INTEREST ----
trt_lg_N = c("formation", "colony.id", "zooid.id", "locality", traits)
dat_lg_N = df[intersect(colnames(df), trt_lg_N)]
head(dat_lg_N) #traits in same order as df and traits

#### CHECK SAMPLE SIZES ----
## number of zooids per colony
range(mean_by_formation_colony$n.zooid)
#5 24 (max 35 for 3 zoo)

#average number of zooids:
mean_by_formation
mean(mean_by_formation$avg.zooid)
sum(mean_by_formation$num.zooid)
sum(mean_by_formation$num.col)

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
              col_form.n[[7]])
zoo.samp <- c(by_form.n[[1]], by_form.n[[2]], by_form.n[[3]],
              by_form.n[[4]], by_form.n[[5]], by_form.n[[6]],
              by_form.n[[7]])
samp <- as.data.frame(cbind(form, col.samp, zoo.samp))
write.csv(samp,
          "./Results/sampling.per.formation.csv",
          row.names = FALSE)

#### P MATRIX ----
p.cov = lapply(form_data, function (x){ (cov(x[, 5:12]))}) #traits per colony (not variation within colony)
#p.cov is the same and the phen.var
#pearson method is the default

##### P VARIANCES ----
##Phenotypic variance in traits and eigen vectors
Pmat = p.cov

lapply(Pmat, isSymmetric)  #is.symmetric.matrix
p.variances = lapply(Pmat, diag)
paste("Trait variances")
head(p.variances)

##### P EIGEN -----

p.eig_variances = lapply(Pmat, function (x) {eigen(x)$values})
# lapply(Pmat, function (x) {eigen(x)})
paste("Eigenvalue variances")
head(p.eig_variances)

p.eig_percent = lapply(p.eig_variances, function (x) {x/sum(x)})
p.eig_per_mat = do.call(rbind, p.eig_percent)
p.eig_per_mat = data.frame(p.eig_per_mat, rownames(p.eig_per_mat))
p.eig_per = melt(p.eig_per_mat)
p.eig_per$rownames.p.eig_per_mat. <- factor(p.eig_per$rownames.p.eig_per_mat., 
                                            levels = c("NKLS", "NKBS", "Tewkesbury",
                                                       "Upper Kai-Iwi", "Tainui",
                                                       "SHCSBSB", "modern"))
    
#dev.off()
P_PC_dist = ggplot(p.eig_per,
                   aes(x = variable, y = value,
                       group = rownames.p.eig_per_mat.,
                       colour = rownames.p.eig_per_mat.)) +
    geom_line(aes(linetype = rownames.p.eig_per_mat.)) +
    geom_point() +
    scale_y_continuous("%Variation in the PC",
                       limits = c(-.02, 0.7)) +
    plot.theme + 
    theme_linedraw() +
    scale_x_discrete("Principal component rank",
                     labels = c("PC1", "PC2", "PC3", "PC4",
                                "PC5", "PC6", "PC7", "PC8")) +
  scale_color_manual(values = col.form)

                      
P_PC_dist #none negative; none above 1; dim 8 close to 0; could keep 7

ggsave(P_PC_dist, 
       file = "./Results/P.PC.dist.png", 
       width = 14, height = 10, units = "cm") 

##### PC LOADINGS -----
#p.eig_variances values = Pmax; loadings 1-8 because of dimensions = PC1-PC8
#what does PC1 correlate with?

#https://www.statology.org/principal-components-analysis-in-r/
P_pc <- lapply(Pmat, function (x) {prcomp(x, scale = TRUE)})

P_pc_NKLS <- P_pc$NKLS
P_pc_NKBS <- P_pc$NKBS
P_pc_tewk <- P_pc$Tewkesbury
P_pc_uki <- P_pc$`Upper Kai-Iwi`
P_pc_tai <- P_pc$Tainui
P_pc_SHCSBSB <- P_pc$SHCSBSB
P_pc_mod <- P_pc$modern

biplot(P_pc_NKLS, scale = 0)
biplot(P_pc_NKBS, scale = 0)
biplot(P_pc_tewk, scale = 0)
biplot(P_pc_uki, scale = 0)
biplot(P_pc_tai, scale = 0)
biplot(P_pc_SHCSBSB, scale = 0)
biplot(P_pc_mod, scale = 0)
#all kinda look the same, with modern looking most different;
#where PC1 is length on left and width on right
#and PC2 is operculum height v everything else

##### P NOISE -----
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
phen.var = lapply(form_data, function (x){ (cov(x[, 5:12]))}) #traits of ALL; correct for colony later
prior = lapply(phen.var, function (x){list(G = list(G1 = list(V = x/2, nu = 2)),
                                           R = list(V = x/4, nu = 2))})

form.mod <- form_data$modern
phen.var.mod = cov(form.mod[, 5:12]) #traits of ALL; correct for colony later
prior.mod = list(G = list(G1 = list(V = phen.var.mod/2, nu = 2),
                           G2 = list(V = phen.var.mod/4,nu = 2)),
                  R = list(V = phen.var.mod/4, nu = 2))

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
                         prior = prior[[i]], verbose = TRUE,
                         pr = TRUE)
}

##need locality data for modern first
model_G.mod <- MCMCglmm(cbind(ln.zh, ln.mpw.b, ln.cw.m, ln.cw.d, #same order as in priors
                              ln.ow.m, ln.oh, ln.c.side, ln.o.side) ~ trait-1, #get rid of "trait" alone?
                        #account for variation w/in colony:
                        random = ~us(trait):colony.id + us(trait):locality, #the number of these determines # of Gs
                        rcov = ~us(trait):units,
                        family = rep("gaussian", 8), #num of traits
                        data = form.mod, #need to have formation as factor in data matrix
                        nitt = 1500000, thin = 1000, burnin = 500000,
                        prior = prior.mod, verbose = TRUE,
                        pr = TRUE)

save(model_G,
     file = "./Results/model_G.RData")

save(model_G.mod,
     file = "./Results/model_G.mod.locality.RData")

model_G_all <- model_G
model_G_all[[7]] <- model_G.mod
save(model_G_all,
     file = "./Results/model_G_all.RData")

load(file = "./Results/model_G_all.RData") #load the g matrices calculated above 


##### CHECK MODELS -----
formation_list #order of formations
summary(model_G_all[[1]])
summary(model_G_all[[2]])
summary(model_G_all[[3]])
summary(model_G_all[[4]])
summary(model_G_all[[5]])
summary(model_G_all[[6]])
summary(model_G_all[[7]])

##plots to see where sampling from:
plot(model_G_all[[1]]$VCV) #catepillar!
plot(model_G_all[[2]]$VCV) #catepillar!
plot(model_G_all[[3]]$VCV) #catepillar!
plot(model_G_all[[4]]$VCV) #catepillar!; high skew
plot(model_G_all[[5]]$VCV) #catepillar!
plot(model_G_all[[6]]$VCV) #catepillar!; high skew
plot(model_G_all[[7]]$VCV) #catepillar!
#formations from oldest to youngest: "NKLS", "NKBS", "Tewkesbury", 
#                                    "Upper Kai-Iwi", "Tainui", 
#                                    "SHCSBSB", "modern"

##### CHECK P IS BIGGER THAN G -----

# diagonals of p.cov to priors
## chekcing NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

# PRIORS
# NKLS
diag(phen.var[[1]]) #all larger
diag(prior$NKLS$G$G1$V) < diag(phen.var[[1]]) #all larger
diag(prior$NKLS$R$V) < diag(phen.var[[1]]) #all larger

# NKBS
diag(phen.var[[2]]) #all larger
diag(prior$NKBS$G$G1$V) < diag(phen.var[[2]]) #all larger
diag(prior$NKBS$R$V) < diag(phen.var[[2]]) #all larger

# Tewkesbury
diag(phen.var[[3]]) #all larger
diag(prior$Tewkesbury$G$G1$V) < diag(phen.var[[3]]) #all larger
diag(prior$Tewkesbury$R$V) < diag(phen.var[[3]]) #all larger

# Upper Kai-Iwi
diag(phen.var[[4]]) #all larger
diag(prior$`Upper Kai-Iwi`$G$G1$V) < diag(phen.var[[4]])
diag(prior$`Upper Kai-Iwi`$R$V) < diag(phen.var[[4]])

# Tainui
diag(phen.var[[5]]) #all larger
diag(prior$Tainui$G$G1$V) < diag(phen.var[[5]])
diag(prior$Tainui$R$V) < diag(phen.var[[5]])

# SHCSBSB
diag(phen.var[[6]]) #all larger
diag(prior$SHCSBSB$G$G1$V) < diag(phen.var[[6]])
diag(prior$SHCSBSB$R$V) < diag(phen.var[[6]])

# modern
diag(phen.var.mod) #all larger
diag(prior.mod$G$G1$V) < diag(phen.var.mod)
diag(prior.mod$R$V) < diag(phen.var.mod)

#RETRIEVE THE P MATRIX FROM THE MCMC OBJECT
# p matrix for each formation (maybe colony?)
## checking NKBS, Waipuru, Upper Kai-Iwi because they are being wonky

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

# Upper Kai-Iwi
post.vcv.uki <- posterior.mode(model_G[[4]]$VCV)
col.vcv.uki <- post.vcv.uki[1:64] #has negatives
d.col.vcv.uki <- col.vcv.uki[c(1, 10, 19, 28, 37, 46, 55, 64)]
unt.vcv.uki <- post.vcv.uki[65:128] #has negatives
#p.unt.vcv.uki <- unt.vcv.uki[c(1,10,19,28,37,46, 55, 64)]
#p.m.uki <- p.col.vcv.uki + p.unt.vcv.uki
p.m.uki <- col.vcv.uki + unt.vcv.uki
d.p.m.uki <- p.m.uki[c(1, 10, 19, 28, 37, 46, 55, 64)]
d.p.m.uki
d.col.vcv.uki #this 'g' matrix are smaller
d.col.vcv.uki < diag(Pmat[[4]]) # zh, mpw.b, cw.m, cw.d, ow.m, oh, c.side, o.side are smaller; all larger than G matrix 

## on all
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

##### POSTERIOR G MATRIX -----

#Retrieving G from posterior
g.model = model_G_all

#identify which columns are holding the posteriors
#(varX, covXY, covYX, varY)
colnames(g.model[[1]]$VCV)
#interested in within colony vcv; so want ~us(trait):colonyid
#have 8 traits, so total number of combos is 8^2 = 64

## on one
summary(g.model[[1]])
#e.g.,                                  post.mean   l-95% CI    u-95% CI    eff.samp
#traitln.zh:traitln.zh.colony.id        0.0065616   0.0046532   0.0091206   1000.0
#just the within colony
#gives summaries of the mean of the variances
summary(g.model[[1]]$VCV)$statistics[1:64,]
#gives mean, sd, se, and quantiles 5% CI (2.5 on either side)
summary(g.model[[1]]$VCV)$quantiles[1:64,]

mean(as.data.frame(g.model[[1]]$VCV[,1:64])$`traitln.zh:traitln.zh.colony.id`)
quantile(as.data.frame(g.model[[1]]$VCV[,1:64])$`traitln.zh:traitln.zh.colony.id`, 
         probs = c(0.05, 0.95))
#these are slightly different from the summary model...why?
#because I think summary does an HPD interval actually
HPDinterval(g.model[[1]]$VCV[,1:64]) #yep, these match
#great, so can recreate stuff



vcv.g.1 <- g.model[[1]]$VCV[,1:64]
#organized as vcv columns (e.g., traitln.zh:traitln.zh.colony.id)
#with each iteration down the size (e.g. 1000)
g.model[[1]]$VCV[1:1000] #gets 1000 variances of the first vcv column 
#traitln.zh:traitln.zh.colony.id 
#                    0.005992251
mean(g.model[[1]]$VCV[1:1000]) #equals what's given in the summary model

posterior.mode(g.model[[1]]$VCV[,1:64]) #get raw mean variances; similar to what is in the summary model
posterior.mode(posterior.cor(g.model[[1]]$VCV[,1:64])) #transforms it to correlation, which ranges from -1 to 1
HPDinterval(posterior.cor(g.model[[1]]$VCV[,1:64])) #get a lower and upper

#from chatgpt
#compare pair-wise
vcv.g.1 <- g.model[[1]]$VCV[,1:64] #1000 variances per variable
vcv.g.2 <- g.model[[2]]$VCV[,1:64]
vcv.g.3 <- g.model[[3]]$VCV[,1:64]
vcv.g.4 <- g.model[[4]]$VCV[,1:64]
vcv.g.5 <- g.model[[5]]$VCV[,1:64]
vcv.g.6 <- g.model[[6]]$VCV[,1:64]
vcv.g.7 <- g.model[[7]]$VCV[,1:64]
# Calculate differences
#distribution of how the covariances vary; i.e., are the differences 0 or not?
diffs.1.2 <- vcv.g.2-vcv.g.1 #differences between all 1000 variances for each variable 
diffs.2.3 <- vcv.g.3-vcv.g.2 #differences between all 1000 variances for each variable 
diffs.3.4 <- vcv.g.4-vcv.g.3 #differences between all 1000 variances for each variable 
diffs.4.5 <- vcv.g.5-vcv.g.4 #differences between all 1000 variances for each variable 
diffs.5.6 <- vcv.g.6-vcv.g.5 #differences between all 1000 variances for each variable 
diffs.6.7 <- vcv.g.7-vcv.g.6 #differences between all 1000 variances for each variable 
# Calculate 95% credible intervals for differences
# since want to check if it overlaps 0 (i.e., no difference, want to get CI around this)
ci_lower.1.2 <- apply(diffs.1.2, 2, quantile, probs = 0.025)
ci_upper.1.2 <- apply(diffs.1.2, 2, quantile, probs = 0.975)

ci_lower.2.3 <- apply(diffs.2.3, 2, quantile, probs = 0.025)
ci_upper.2.3 <- apply(diffs.2.3, 2, quantile, probs = 0.975)

ci_lower.3.4 <- apply(diffs.3.4, 2, quantile, probs = 0.025)
ci_upper.3.4 <- apply(diffs.3.4, 2, quantile, probs = 0.975)

ci_lower.4.5 <- apply(diffs.4.5, 2, quantile, probs = 0.025)
ci_upper.4.5 <- apply(diffs.4.5, 2, quantile, probs = 0.975)

ci_lower.5.6 <- apply(diffs.5.6, 2, quantile, probs = 0.025)
ci_upper.5.6 <- apply(diffs.5.6, 2, quantile, probs = 0.975)

ci_lower.6.7 <- apply(diffs.6.7, 2, quantile, probs = 0.025)
ci_upper.6.7 <- apply(diffs.6.7, 2, quantile, probs = 0.975)

# Check if the intervals contain zero (i.e., no difference)
#so it is asking if the lower is above zero or the upper is below (i.e., does not contain zero)
#false means no difference,
#because it means the lower is not greater than 0 and the upper is not lower than 0
different.1.2 <- (ci_lower.1.2 > 0 | ci_upper.1.2 < 0)
unique(different.1.2) 

different.2.3 <- (ci_lower.1.2 > 0 | ci_upper.1.2 < 0)
different.3.4 <- (ci_lower.1.2 > 0 | ci_upper.1.2 < 0)
different.4.5 <- (ci_lower.1.2 > 0 | ci_upper.1.2 < 0)
different.5.6 <- (ci_lower.1.2 > 0 | ci_upper.1.2 < 0)
different.6.7 <- (ci_lower.1.2 > 0 | ci_upper.1.2 < 0)

unique(different.2.3) #no diff
unique(different.3.4) #no diff
unique(different.4.5) #no diff
unique(different.5.6) #no diff
unique(different.6.7) #no diff

# Print results
if (any(different.1.2)) {
    print("The covariance matrices are significantly different.")
} else {
    print("No significant difference between covariance matrices.")
}

#using hdp:
# Calculate HPD intervals for each element of the difference
hpd_intervals.1.2 <- apply(diffs.1.2, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95))  # 95% HPD interval

hpd_intervals.2.3 <- apply(diffs.2.3, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95))  # 95% HPD interval
hpd_intervals.3.4 <- apply(diffs.3.4, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95))  # 95% HPD interval
hpd_intervals.4.5 <- apply(diffs.4.5, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95))  # 95% HPD interval
hpd_intervals.5.6 <- apply(diffs.5.6, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95))  # 95% HPD interval
hpd_intervals.6.7 <- apply(diffs.6.7, 2, function(x) HPDinterval(as.mcmc(x), prob = 0.95))  # 95% HPD interval

# Check if zero is included in the HPD intervals
#true means no difference
#here it is asking if it overlaps 0, so true means that it does
contains_zero.1.2 <- apply(hpd_intervals.1.2, 2, function(x) (x[1] <= 0 & x[2] >= 0))
contains_zero.2.3 <- apply(hpd_intervals.2.3, 2, function(x) (x[1] <= 0 & x[2] >= 0))
contains_zero.3.4 <- apply(hpd_intervals.3.4, 2, function(x) (x[1] <= 0 & x[2] >= 0))
contains_zero.4.5 <- apply(hpd_intervals.4.5, 2, function(x) (x[1] <= 0 & x[2] >= 0))
contains_zero.5.6 <- apply(hpd_intervals.5.6, 2, function(x) (x[1] <= 0 & x[2] >= 0))
contains_zero.6.7 <- apply(hpd_intervals.6.7, 2, function(x) (x[1] <= 0 & x[2] >= 0))
unique(contains_zero.1.2)

# Print results
if (any(!contains_zero.1.2)) {
    print("The covariance matrices are significantly different.")
} else {
    print("No significant difference between covariance matrices.")
}

#let's see if I get the same as above if I compare all the pairwise CI from each variable
#extract l-95% CI and u-95% CI for each
#see if for each variable, if the mean variance of time 1 is within the CI of time 1
vcv.g1 <- as.data.frame(g.model[[1]]$VCV[,1:64])
vcv.g2 <- as.data.frame(g.model[[2]]$VCV[,1:64])

#get all the means
mean.g1 <- apply(vcv.g1, 2, mean) %>% as.data.frame()
colnames(mean.g1) <- "mean"
mean.g1$variable <- rownames(mean.g1)

mean.g2 <- apply(vcv.g2, 2, mean) %>% as.data.frame()
colnames(mean.g2) <- "mean"
mean.g2$variable <- rownames(mean.g2)

#get all the ci
l.ci.g1 <- apply(vcv.g1, 2, quantile, probs = 0.05) %>% as.data.frame()
colnames(l.ci.g1) <- "lower.ci"
l.ci.g1$variable <- rownames(l.ci.g1)
u.ci.g1 <- apply(vcv.g1, 2, quantile, probs = 0.95) %>% as.data.frame()
colnames(u.ci.g1) <- "upper.ci"
u.ci.g1$variable <- rownames(u.ci.g1)
ci.g1 <- merge(l.ci.g1, u.ci.g1, by = "variable")

l.ci.g2 <- apply(vcv.g2, 2, quantile, probs = 0.05) %>% as.data.frame()
colnames(l.ci.g2) <- "lower.ci"
l.ci.g2$variable <- rownames(l.ci.g2)
u.ci.g2 <- apply(vcv.g2, 2, quantile, probs = 0.95) %>% as.data.frame()
colnames(u.ci.g2) <- "upper.ci"
u.ci.g2$variable <- rownames(u.ci.g2)
ci.g2 <- merge(l.ci.g2, u.ci.g2, by = "variable")

#get all the hpd intervals
var.g1 <- rownames(as.data.frame(HPDinterval(g.model[[1]]$VCV[,1:64])))
l.hpd.g1 <- as.data.frame(HPDinterval(g.model[[1]]$VCV[,1:64]))$lower
u.hpd.g1 <- as.data.frame(HPDinterval(g.model[[1]]$VCV[,1:64]))$upper
hdp.g1 <- as.data.frame(cbind(var.g1, l.hpd.g1, u.hpd.g1))
colnames(hdp.g1)[1] <- "variable"

var.g2 <- rownames(as.data.frame(HPDinterval(g.model[[2]]$VCV[,1:64])))
l.hpd.g2 <- as.data.frame(HPDinterval(g.model[[2]]$VCV[,1:64]))$lower
u.hpd.g2 <- as.data.frame(HPDinterval(g.model[[2]]$VCV[,1:64]))$upper
hdp.g2 <- as.data.frame(cbind(var.g2, l.hpd.g2, u.hpd.g2))
colnames(hdp.g2)[1] <- "variable"

#bring together
stats.g1 <- merge(mean.g1, ci.g1, by = "variable")
stats.g1 <- merge(stats.g1, hdp.g1, by = "variable")

stats.g2 <- merge(mean.g2, ci.g2, by = "variable")
stats.g2 <- merge(stats.g2, hdp.g2, by = "variable")

stats.g1$l.hpd.g1 <- as.numeric(stats.g1$l.hpd.g1)
stats.g1$u.hpd.g1 <- as.numeric(stats.g1$u.hpd.g1)
stats.g2$l.hpd.g2 <- as.numeric(stats.g2$l.hpd.g2)
stats.g2$u.hpd.g2 <- as.numeric(stats.g2$u.hpd.g2)

#see if means overlap ci and hpd
stats.g1$lower.ci <= stats.g2$mean & stats.g1$upper.ci >= stats.g2$mean

stats.g1$l.hpd.g1 <= stats.g2$mean & stats.g1$u.hpd.g1 >= stats.g2$mean

#see if ci and hpd overlap
stats.g1$lower.ci <= stats.g2$lower.ci & stats.g1$upper.ci >= stats.g2$lower.ci |
    stats.g1$lower.ci <= stats.g2$upper.ci & stats.g1$upper.ci >= stats.g2$upper.ci 

stats.g1$l.hpd.g1 <= stats.g2$l.hpd.g2 & stats.g1$u.hpd.g1 >= stats.g2$l.hpd.g2 |
    stats.g1$l.hpd.g1 <= stats.g2$u.hpd.g2 & stats.g1$u.hpd.g1 >= stats.g2$u.hpd.g2

#seems to be so

##### REIMANN DISTANCE -----
# look at magnitude of beta (x axis) as a function of reimann distance
# this is the distance (subtraction) between two matrices

dist.NKLS_NKBS <- RiemannDist(G_ext_NKBS, G_ext_NKLS)
dist.NKBS_tewk <- RiemannDist(G_ext_NKLS, G_ext_tewk)
dist.tewk_uki <- RiemannDist(G_ext_tewk, G_ext_uki)
dist.uki_tai <- RiemannDist(G_ext_uki, G_ext_tai)
dist.tai_SHCSBSB <- RiemannDist(G_ext_tai, G_ext_SHCSBSB)
dist.SHCSBSB_mod <- RiemannDist(G_ext_SHCSBSB, G_ext_mod)

dist.gmat <- c(dist.NKLS_NKBS, dist.NKBS_tewk,
               dist.tewk_uki, dist.uki_tai, 
               dist.tai_SHCSBSB, dist.SHCSBSB_mod, "")

# are these explained by sampling error alone?


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

##### G EIGEN -----
g.eig_variances = lapply(Gmat, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(g.eig_variances)

g.eig_percent = lapply(g.eig_variances, function (x) {x/sum(x)})
g.eig_per_mat = do.call(rbind, g.eig_percent)
g.eig_per_mat = data.frame(g.eig_per_mat, rownames(g.eig_per_mat))
g.eig_per = melt(g.eig_per_mat)
g.eig_per$rownames.g.eig_per_mat. <- factor(g.eig_per$rownames.g.eig_per_mat., 
                                            levels = c("NKLS", "NKBS", "Tewkesbury",
                                                       "Upper Kai-Iwi", "Tainui",
                                                       "SHCSBSB", "modern"))

#dev.off()
G_PC_dist = ggplot(g.eig_per,
                   aes(x = variable, y = value,
                       group = rownames.g.eig_per_mat.,
                       colour = rownames.g.eig_per_mat.)) +
    geom_line(aes(linetype = rownames.g.eig_per_mat.)) +
    geom_point() +
    plot.theme + 
    theme_linedraw() +
    scale_x_discrete("Principal component rank",
                     labels = c("PC1", "PC2", "PC3", "PC4",
                                "PC5", "PC6", "PC7", "PC8")) +
    scale_y_continuous("%Variation in the PC",
                       limits = c(-.02, 0.7)) +
    scale_color_manual(values = col.form)
G_PC_dist 
#none above 1
#modern negative at 6...dim 5
#upper kai iwi still weird
#have negative values because have negative variances,
#meaning that the PC value does not do a good job explaining the variance

ggsave(G_PC_dist, 
       file = "./Results/G.PC.dist.png", 
       width = 14, height = 10, units = "cm")

#Note that some matrices have negative eigenvalues. 
#This can cause a lot of trouble in analyses involving inverted matrices.
#Solution from evolqg Marroig et al. 2012

##### PC LOADINGS -----
#p.eig_variances values = Pmax; loadings 1-8 because of dimensions = PC1-PC8
#what does PC1 correlate with?

#https://www.statology.org/principal-components-analysis-in-r/
G_pc <- lapply(Gmat, function (x) {prcomp(x, scale = TRUE)})
G_pc_NKLS <- G_pc$NKLS
G_pc_NKBS <- G_pc$NKBS
G_pc_tewk <- G_pc$Tewkesbury
G_pc_uki <- G_pc$`Upper Kai-Iwi`
G_pc_tai <- G_pc$Tainui
G_pc_SHCSBSB <- G_pc$SHCSBSB
G_pc_mod <- G_pc$modern

biplot(G_pc_NKLS, scale = 0)
biplot(G_pc_NKBS, scale = 0)
biplot(G_pc_tewk, scale = 0)
biplot(G_pc_uki, scale = 0)
biplot(G_pc_tai, scale = 0)
biplot(G_pc_SHCSBSB, scale = 0)
biplot(G_pc_mod, scale = 0)
#all kinda look the same as the P matrix with modern and uki looking most different, 
#where PC1 is length on left and width on right
#and PC2 is operculum height v everything else

##### G NOISE -----
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
data.list = list(Pmat, Gmat,
                 p.eig_variances, g.eig_variances,
                 P_ext, G_ext,
                 df, dat_lg_N,
                 form_data, by_form.n, col_form.n)
save(data.list,
     file = "./Results/data.list.RData")

load(file = "./Results/data.list.RData") #load the g matrices calculated above 


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

plot(diag(Gmat[[2]]), diag(Pmat[[2]]),
     pch = 19, col = col.form[2],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "NKBS")
abline(0, 1)

plot(diag(Gmat[[3]]), diag(Pmat[[3]]),
     pch = 19, col = col.form[3],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tewkesbury",
     xlim = c(0, 0.01),
     ylim = c(0, 0.05))
abline(0, 1) 

plot(diag(Gmat[[4]]), diag(Pmat[[4]]),
     pch = 19, col = col.form[4],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonals",
     main = "Upper Kai-Iwi")
abline(0, 1) 

plot(diag(Gmat[[5]]), diag(Pmat[[5]]),
     pch = 19, col = col.form[5],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "Tainui")
abline(0, 1) 

plot(diag(Gmat[[6]]), diag(Pmat[[6]]),
     pch = 19, col = col.form[6],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1) 

plot(diag(Gmat[[7]]), diag(Pmat[[7]]),
     pch = 19, col = col.form[7],
     xlab = "G non-standardized diagonal",
     ylab = "P non-standardized diagonal",
     main = "SHCSBSB",
     xlim = c(0, 0.01),
     ylim = c(0,0.05))
abline(0, 1) 

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
mean.df.sub <- mean_by_formation_colony[, c(1, 4, 6, 8, 10, 12, 14, 16, 18)] 
#remove rows with NA so that the sample size gets correct when we plot 
#the sample size for the VCOV.
mean.df_complete_cases <- mean.df.sub[complete.cases(mean.df.sub), ] #remove rows with NA so that the sample size gets correct when we plot the sample size for the VCOV.

# Calculating the sample size for each time point
colony_samples = split.data.frame(mean.df_complete_cases, mean.df_complete_cases$formation) #Sample size
sample_sizes_G = lapply(colony_samples, function(x){dim(x)[1]})

comp_sampleN = matrix(0, 7, 7) #calculating the smallest sample size in the comparison among pairs of G for length of formation

for (i in 1:length(sample_sizes_G)){
    for (j in 1:length(sample_sizes_G)){
        comp_sampleN[i,j] = min(as.numeric(sample_sizes_G[i]), as.numeric(sample_sizes_G[j]))
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
     pch = 19, col = "grey", cex = .5)
points(obs_melt$N, obs_melt$RS, 
       col = "#00BFC4", pch = 19, cex = .5)

## none are outside the gray
# previously the modern sites were, but it seems that was due to locality
# and now that has been corrected for!

p.rare <- ggplot() +
    geom_point(aes(out_results[, 2], out_results[, 1]),
               pch = 19, col = "grey", size = 2) +
    geom_point(aes(obs_melt$N, obs_melt$RS),
               col = "#00BFC4", pch = 19, size = 2) + 
    plot.theme +
    scale_x_continuous(name = "Sample size") +
    scale_y_continuous(name = "Similarity") 

ggsave(p.rare, 
       file = "./Results/rarefaction.png", 
       width = 14, height = 10, units = "cm")

#### RAREFACTION FOR P ----
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
                                                  P_ext[[1]])}) # Simulate individual trait data based on the first VCOV matrix in our time series and different sample sizes. 
# Calculate VCOV matrices 
pop.cv <- lapply(pop.dist, cov) 
# Estimate selection gradients based on Random Skewers, comparing our true 
# G matrix and the undersampled pooled P matrix. 
RS_result = RandomSkewers(pop.cv, 
                          P_ext[[1]]) 
out_results <- cbind(RS_result$correlation, pop.size)

# Computing Random Skewers between pairs of actual G matrices in our data
#(n = 7 for 7 formations) 
# sub-setting so that we only analyze the traits and formation of interest
mean.df.sub <- mean_by_formation_colony[, c(1, 4, 6, 8, 10, 12, 14, 16, 18)] 
#remove rows with NA so that the sample size gets correct when we plot 
#the sample size for the VCOV.
mean.df_complete_cases <- mean.df.sub[complete.cases(mean.df.sub), ] #remove rows with NA so that the sample size gets correct when we plot the sample size for the VCOV.

# Calculating the sample size for each time point
colony_samples = split.data.frame(mean.df_complete_cases, mean.df_complete_cases$formation) #Sample size
sample_sizes_P = lapply(colony_samples, function(x){dim(x)[1]})

comp_sampleN = matrix(0, 7, 7) #calculating the smallest sample size in the comparison among pairs of G for length of formation

for (i in 1:length(sample_sizes_P)){
    for (j in 1:length(sample_sizes_P)){
        comp_sampleN[i,j] = min(as.numeric(sample_sizes_P[i]), as.numeric(sample_sizes_P[j]))
    }
}

# Here, we do the actual Random Skewers test and 
# combined the results (correlations in the response to selection) 
# with the smallest sample size for the pairs of G that are investigated
comp_mat = RandomSkewers(P_ext)
melt_comp = melt(comp_mat$correlations[lower.tri(comp_mat$correlations)])
melt_samples = melt(comp_sampleN[lower.tri(comp_sampleN)])
obs_melt = cbind(melt_comp, melt_samples)
colnames(obs_melt) = c("RS","N")

#Plotting the result
plot(out_results[, 2], out_results[, 1], 
     xlab = "Sample size", ylab = "Similarity", 
     pch = 19, col = "grey", cex = .5)
points(obs_melt$N, obs_melt$RS, 
       col = "#00BFC4", pch = 19, cex = .5)

## none are outside the gray
# previously the modern sites were, but it seems that was due to locality
# and now that has been corrected for!

p.rare <- ggplot() +
    geom_point(aes(out_results[, 2], out_results[, 1]),
               pch = 17, col = "grey", size = 2) +
    geom_point(aes(obs_melt$N, obs_melt$RS),
               col = "#00BFC4", pch = 17, size = 2) + 
    plot.theme +
    scale_x_continuous(name = "Sample size") +
    scale_y_continuous(name = "Similarity") 

ggsave(p.rare, 
       file = "./Results/rarefaction_p.png", 
       width = 14, height = 10, units = "cm")

#### GLOBAL G ----

##### PRIORS -----

#dat_lg_N.com = dat_lg_N[complete.cases(dat_lg_N),] #didn't fix anything
phen.var.glob = cov(dat_lg_N[, 5:12]) #traits of ALL; correct for colony and formation later
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
     file = "./Results/global_matrix.RData")

load(file="./Results/global_matrix.RData") #load the g matrices calculated above 
model_Global

##### CHECK MODELS -----
formation_list #order of formations
summary(model_Global)

##plots to see where sampling from:
plot(model_Global$VCV) #catepillar!

###### POSTERIOR GLOBAL G MATRIX ------
#Retrieving GLOBAL G from posterior
glob.model = model_Global
ntraits = 8
glob.Gmat = matrix(posterior.mode(glob.model$VCV)[1:ntraits^2], ntraits, ntraits)

##### GLOBAL G VARIANCES -----
isSymmetric(glob.Gmat)  #is.symmetric.matrix
glob.variances = diag(glob.Gmat)
paste("Trait variances")
head(glob.variances)

###### GLOBAL G EIGEN ------
glob.eig_variances = eigen(glob.Gmat)$values
paste("Eigenvalue variances")
head(glob.eig_variances)

glob.eig_percent = glob.eig_variances/sum(glob.eig_variances) 
glob.eig_per_mat = data.frame(glob.eig_percent)
#glob.eig_per = melt(glob.eig_per_mat)
glob.eig_per_mat$PC <- 1:nrow(glob.eig_per_mat)
#Sdev.off()
#glob.eig_per_mat$PC <- as.factor(glob.eig_per_mat$PC)
Glob_PC_dist = ggplot(glob.eig_per_mat,
                   aes(x = PC, y = glob.eig_percent)) +
  geom_point()  +
    geom_line() +
    plot.theme + 
    theme_linedraw() +
    scale_x_continuous("Principal component rank",
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8),
                     labels = c("PC1", "PC2", "PC3", "PC4",
                                "PC5", "PC6", "PC7", "PC8")) +
    scale_y_continuous("%Variation in the PC",
                       limits = c(-.02, 0.7))
Glob_PC_dist #one negative; none above 1!

ggsave(Glob_PC_dist, 
       file = "./Results/GlobalG.PC.dist.png", 
       width = 14, height = 10, units = "cm")

#Note that some matrices have negative eigenvalues. 
#This can cause a lot of trouble in analyses involving inverted matrices.
#Solution from evolqg Marroig et al. 2012

###### GLOBAL G NOISE ------
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
     file = "./Results/global_ext.RData")

load(file="./Results/global_ext.RData") #load the g matrices calculated above 
Glob_ext

