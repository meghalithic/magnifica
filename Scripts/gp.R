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
form_data = lapply(by_form, function(x) x[complete.cases(x),])
p.cov = lapply(form_data, function (x){ (cov(x[, 4:11]))}) #traits per colony (not variation within colony)

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
model = model_G
data = (dat_lg_N)
ntraits = 8
Gmat = lapply(model, function (x) { 
  matrix(posterior.mode(x$VCV)[1:ntraits^2], ntraits, ntraits)})


###### G STD ------

#Standardizing G by the trait means 

mean_by_form = setDT(na.omit(data[, 3:11]))[, lapply(.SD, mean, na.rm = F), 
                                            by = .(formation)] #traits + formation
u_form = split(mean_by_form, mean_by_form$formation)
test = lapply(u_form, function (x){ data.matrix(x[, 2:9])}) #only traits from mean_by_form
test_std = lapply(test, function (x){(as.numeric(x))%*%t(as.numeric(x))})
G_std = list()
for (i in 1:length(Gmat)){
  G_std[[i]] = Gmat[[i]]/(test_std[names(form_data[i])][[1]])
}
G_std
names(G_std)=names(form_data[1:i]) #only one negative; formation 4 (Tainui)

##Genetic variance in traits and eigenvectors

#load(file="New_g_matrices.RData") #load the g matrices calculated above 

lapply(G_std, isSymmetric)  #is.symmetric.matrix
std_variances = lapply(G_std, diag)
paste("Trait variances")
head(std_variances)

#require(matrixcalc)
#m.1 <- round(G_std[[1]], 10)
#is.symmetric.matrix(m.1)

#is.positive.definite(m.1) 
#no spot with zero variance; 
#so variance along certain directions are negative

###### G EIGEN ------

eig_variances=lapply(G_std, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(eig_variances)

eig_percent = lapply(eig_variances, function (x) {x/sum(x)})
eig_per_mat = do.call(rbind, eig_percent)
eig_per_mat = data.frame(eig_per_mat, rownames(eig_per_mat))
eig_per = melt(eig_per_mat)
#dev.off()
PC_dist = ggplot(eig_per,
                 aes(x = variable, y = value,
                     group = rownames.eig_per_mat.,
                     colour = rownames.eig_per_mat.)) +
  geom_line(aes(linetype = rownames.eig_per_mat.)) +
  geom_point() +
  xlab("Principal component rank") +
  ylab("%Variation in the PC")
PC_dist #one negative; none above 1!

ggsave(PC_dist, file = "./Results/PC_dist_form.png", 
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
# combined the results (correlations in the respons to selection) 
# with the smallest samle size for the pairs of G that are investigated
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
       col = "purple", pch = 19)






# G-based neutrality test
##Cladogenetic
clado_vecs=function (ancestor,descendant, dataframe){
  anc_matrix=dataframe[names(dataframe)==ancestor][[1]]
  desc_matrix=dataframe[names(dataframe)==descendant][[1]]
  if (max(anc_matrix$age)<=max(desc_matrix$age)){
    print('The descendant species is older or, at least, as old as the ancestral')
  }
  else{
    desc_value=desc_matrix[nrow(desc_matrix)][,4:11]
    anc_values=anc_matrix[anc_matrix$age>max(desc_matrix$age[nrow(desc_matrix)])][,4:11]
    std_diff=-sweep(as.matrix(anc_values),2,
                    as.matrix(desc_value))/as.matrix(anc_values)
    std_diff}
}

clado_vec=list()
clado_vec[[1]]=clado_vecs('coll','auric',age_sorted_clado)
clado_vec[[2]]=clado_vecs('auric','nsp09',age_sorted_clado)
clado_vec[[3]]=clado_vecs('nsp09','nsp10',age_sorted_clado)
clado_vec[[4]]=clado_vecs('lacry','nsp03',age_sorted_clado)
clado_vec[[5]]=clado_vecs('nsp03','nsp04',age_sorted_clado)

div_sp_std_mean=rbind(clado_vec[[2]],clado_vec[[3]],
                      clado_vec[[4]],clado_vec[[5]])
clado_B=t(div_sp_std_mean)%*%div_sp_std_mean/8
Gselect=Reduce("+", G_ext) / length(G_ext)#mean G to make drift tests during cladogenesis
pc_G_clado=eigen(Gselect)$vectors
var_B_clado=diag(t(pc_G_clado)%*%clado_B%*%pc_G_clado)
var_G_clado=eigen(Gselect)$values
Fig_drift_clado=ggplot(data=data.frame(var_G_clado,var_B_clado),
                       aes(x=log(var_G_clado),y=log(var_B_clado)))+
  geom_point()+geom_smooth(method="lm")
fit_obj_clado=lm(log(var_B_clado)~log(var_G_clado))
conf_int_drift_clado=stats::confint(fit_obj_clado)
head(conf_int_drift_clado)

#Anagenetic lineages
drift_N=lapply(results_mean, function(x) {sum(x!="NA")/8})
drift_N=drift_N[drift_N>9]
ana_drift=subset(results_mean,names(results_mean)%in%names(drift_N))
dft=lapply(ana_drift,function (x){t(x)%*%x})
Fig_drift=list()
conf_int_drift=list()
for (i in 1:5){
  
  dft[[i]]=dft[[i]]/drift_N[[i]]
  pc_G=eigen(G_std[[i]])$vectors
  var_B=diag(t(pc_G)%*%dft[[i]]%*%pc_G)
  var_G=eigen(G_std[[i]])$values
  Fig_drift[[i]]=ggplot(data=data.frame(var_G,var_B), aes(x=log(var_G),
                                                          y=log(var_B)))+
    geom_point()+geom_smooth(method="lm")
  fit_obj=lm(log(var_B)~log(var_G))
  conf_int_drift[[i]]=stats::confint(fit_obj)
  
}
names(conf_int_drift)=names(G_std[-4])
head(conf_int_drift)

##Is morphological evolution occurring along directions of high evolvability?
Bet_sp=apply(div_sp_std_mean,1,Normalize)
colSums(Bet_sp^2)#making sure they are normalized (should sum up to 1)

#within species
results_sub=subset(results_mean,names(results_mean)%in%names(drift_N))
w_sp=lapply(results_sub, function(x){apply(x,1,Normalize)}) 
lapply(w_sp,function (x){colSums(x^2)})

#Random correlations
num.traits <- 8
iterations=400
beta.mat <- array(rnorm(num.traits * iterations), c(num.traits,iterations))
beta.mat <- apply(beta.mat, 2, Normalize)
bet_corr=t(beta.mat)%*%beta.mat
paste("95%CI for vector correlations between random vectors")
rand_corr=quantile(as.vector(abs(bet_corr)),0.95)
rand_corr

#Deltaz vectors + Random Vectors
beta.mat.clado=Bet_sp
beta.mat.ana=w_sp
evolv_ana=list()
for (i in 1:5){
  evolv_ana[[i]]=diag(t(beta.mat.ana[[i]]) %*% G_ext[[i]] %*% beta.mat.ana[[i]])
} 
evolv_clado=diag(t(beta.mat.clado) %*% Gselect %*% beta.mat.clado)

evolv_mat_rand=diag(t(beta.mat) %*% Gselect %*% beta.mat)

evolv_final=list(unlist(evolv_ana),evolv_clado,evolv_mat_rand)
names(evolv_final)=c("Anagenetic", "Cladogenetic", "Random")
evolv_melted=melt(evolv_final)
evolv_ggp=ggplot(data=evolv_melted, aes(x=L1, y=value, fill=L1))+ 
  geom_boxplot(alpha=0.6)+
  ylab("Evolvability")+xlab(NULL)+ theme(legend.position='none')

c_evolv_ana=list()
for (i in 1:5){
  c_evolv_ana[[i]]=(1/diag(t(beta.mat.ana[[i]]) %*% solve(G_ext[[i]],
                                                          beta.mat.ana[[i]])))
} 
c_evolv_clado=(1/diag(t(beta.mat.clado) %*% solve(Gselect, beta.mat.clado)))

c_evolv_mat_rand=(1/diag(t(beta.mat) %*% solve(Gselect, beta.mat)))

c_evolv_final=list(unlist(c_evolv_ana),c_evolv_clado,c_evolv_mat_rand)
names(c_evolv_final)=c("Anagenetic", "Cladogenetic", "Random")
c_evolv_melted=melt(c_evolv_final)
cevol_ggp=ggplot(data=c_evolv_melted, aes(x=L1, y=value, fill=L1))+ 
  geom_boxplot(alpha=0.6)+ylab("Conditional Evolvability") +
  xlab(NULL)+ theme(legend.position='none')

grid.arrange(evolv_ggp,cevol_ggp,nrow=1,ncol=2)

# Selection
#Comparable anagenetic vectors
ana_morph_vec=function (ancestor,descendant, dataframe){
  anc_matrix=dataframe[names(dataframe)==ancestor][[1]]
  desc_matrix=dataframe[names(dataframe)==descendant][[1]]
  anc_values=anc_matrix[anc_matrix$age>max(desc_matrix$age[nrow(desc_matrix)])][,4:11]
  ana_values=anc_matrix[anc_matrix$age<=max(desc_matrix$age[nrow(desc_matrix)])][,4:11]
  std_diff=-sweep(as.matrix(anc_values),2,
                  as.matrix(ana_values[nrow(ana_values)]))/as.matrix(anc_values)
  std_diff
}


ana_vec=list()
ana_vec[[1]]=ana_morph_vec('coll','auric',age_sorted_clado)
ana_vec[[2]]=ana_morph_vec('auric','nsp09',age_sorted_clado)
ana_vec[[3]]=ana_morph_vec('nsp09','nsp10',age_sorted_clado)
ana_vec[[4]]=ana_morph_vec('lacry','nsp03',age_sorted_clado)
ana_vec[[5]]=ana_morph_vec('nsp03','nsp04',age_sorted_clado)

w_sp_beta=list()
w_sp_beta[[1]]=ana_vec[[2]]%*% solve(G_ext$auric)
w_sp_beta[[2]]=ana_vec[[3]]%*% solve(G_ext$nsp09)
w_sp_beta[[3]]=ana_vec[[4]]%*% solve(G_ext$lacry)
w_sp_beta[[4]]=ana_vec[[5]]%*% solve(G_ext$lacry)

beta_ana=unlist(lapply(w_sp_beta, function (x){colSums(t(x^2))^0.5}))

b_sp_beta=list()
b_sp_beta[[1]]=div_sp_std_mean[1:2,]%*% solve(G_ext$auric)
b_sp_beta[[2]]=div_sp_std_mean[3:4,]%*% solve(G_ext$nsp09)
b_sp_beta[[3]]=div_sp_std_mean[5,]%*% solve(G_ext$lacry)
b_sp_beta[[4]]=div_sp_std_mean[6:8,]%*% solve(G_ext$lacry)


beta_clado= unlist(lapply(b_sp_beta, function (x){colSums(t(x^2))^0.5}))

ana_lineage=NULL
subset_results=results_mean[names(G_ext)]
for (i in 1:5){ana_lineage[i]=lapply(subset_results,function (x){
  as.matrix(x) %*% solve(G_ext[[i]])})}
beta_ana_lineage=unlist(lapply(ana_lineage, function (x){colSums(t(x^2))^0.5}))

anag_p=as.vector(beta_ana)
anag_p=cbind(rep(3,length(anag_p)),anag_p)

clado_p=as.vector(beta_clado)
clado_p=cbind(rep(1,length(clado_p)),clado_p)

anag_l=as.vector(beta_ana_lineage)
anag_l=cbind(rep(2,length(anag_l)),anag_l)

selection_grad=as.data.frame(rbind(anag_p,clado_p,anag_l))
names(selection_grad)=c("T1","T2")
select_fig=ggplot(as.data.frame(selection_grad), aes(x=as.factor(T1), y=T2, 
                                                     fill=as.factor(T1)))+geom_boxplot(alpha=0.6)+
  ylab("Selection gradient (unitless elasticity)")+xlab(NULL)+ 
  theme(legend.position='none')
select_fig


#Permutation tests
library(lmPerm)
#random vs clado - conditional evolvability
c_evolv_melted_cladoRand=c_evolv_melted[c_evolv_melted$L1=='Random'|
                                          c_evolv_melted$L1=='Cladogenetic',]
t=aovp(value ~ L1,data=c_evolv_melted_cladoRand)
summary(t)

#random vs ana - conditional evolvability
c_evolv_melted_anaRand=c_evolv_melted[c_evolv_melted$L1=='Random'|
                                        c_evolv_melted$L1=='Anagenetic',]
t1=aovp(value ~ L1,data=c_evolv_melted_anaRand)
summary(t1)

#clado vs ana - conditional evolvability
c_evolv_melted_anaClado=c_evolv_melted[c_evolv_melted$L1=='Cladogenetic'|
                                         c_evolv_melted$L1=='Anagenetic',]
t2=aovp(value ~ L1,data=c_evolv_melted_anaClado)
summary(t2)

#random vs clado - evolvability
evolv_melted_cladoRand=evolv_melted[evolv_melted$L1=='Random'|
                                      evolv_melted$L1=='Cladogenetic',]
t3=aovp(value ~ L1,data=evolv_melted_cladoRand)
summary(t3)

#random vs ana - evolvability
evolv_melted_anaRand=evolv_melted[evolv_melted$L1=='Random'|
                                    evolv_melted$L1=='Anagenetic',]
t4=aovp(value ~ L1,data=evolv_melted_anaRand)
summary(t4)

#clado vs ana - evolvability
evolv_melted_anaClado=evolv_melted[evolv_melted$L1=='Cladogenetic'|
                                     evolv_melted$L1=='Anagenetic',]
t5=aovp(value ~ L1,data=evolv_melted_anaClado)
summary(t5)

#clado vs ana -  rates

melt_df_clado_ana=melt_df[melt_df$variable=='Cladogenesis'|
                            melt_df$variable=='Anagenesis-Comparable',]
t6=aovp(value ~ variable,data=melt_df_clado_ana)
summary(t6)

melt_df_clado_ana2=melt_df[melt_df$variable=='Cladogenesis'|
                             melt_df$variable=='Anagenesis-Sequential',]
t7=aovp(value ~ variable,data=melt_df_clado_ana2)
summary(t7)

#clado vs ana -  selection

gradient_df_clado_ana=as.data.frame(selection_grad)
t8=aovp(T2 ~ T1,data=gradient_df_clado_ana)
summary(t8)