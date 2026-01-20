source("./Scripts/0-env.R")
load(file = "./Results/evol.diff.list.RData") #load the g matrices calculated above 

G_ext_NKLS <- evol.diff.list[[1]]
G_ext_NKBS <- evol.diff.list[[2]]
G_ext_tewk <- evol.diff.list[[3]]
G_ext_uki <- evol.diff.list[[4]]
G_ext_tai <- evol.diff.list[[5]]
G_ext_SHCSBSB <- evol.diff.list[[6]]
G_ext_mod <- evol.diff.list[[7]]
Gmax_NKLS_norm <- evol.diff.list[[8]]
Gmax_NKBS_norm <- evol.diff.list[[9]]
Gmax_tewk_norm <- evol.diff.list[[10]]
Gmax_uki_norm <- evol.diff.list[[11]]
Gmax_tai_norm <- evol.diff.list[[12]]
Gmax_SHCSBSB_norm <- evol.diff.list[[13]]
Gmax_mod_norm <- evol.diff.list[[14]]
evolved_difference_unit_length_t1 <- evol.diff.list[[15]]
evolved_difference_unit_length_t2 <- evol.diff.list[[16]]
evolved_difference_unit_length_t3 <- evol.diff.list[[17]]
evolved_difference_unit_length_t4 <- evol.diff.list[[18]]
evolved_difference_unit_length_t5 <- evol.diff.list[[19]]
evolved_difference_unit_length_t6 <- evol.diff.list[[20]]

load(file = "./Results/sum.data.list.RData") #load the g matrices calculated above 
mean_by_formation <- sum.data.list[[1]]
mean_by_formation_colony <- sum.data.list[[2]]
means <- sum.data.list[[3]]

load(file = "./Results/data.list.RData") #load the g matrices calculated above 
Pmat <- data.list[[1]]
Gmat <- data.list[[2]]
p.eig_variances <- data.list[[3]]
g.eig_variances <- data.list[[4]]
P_ext <- data.list[[5]]
G_ext <- data.list[[6]]
df <- data.list[[7]]
dat_lg_N <- data.list[[8]]
form_data <- data.list[[9]]
by_form.n <- data.list[[10]]
col_form.n <- data.list[[11]]

#check number of zooids NOT colonies:
# by colonies use mean_by_formation_colony
# by zooid us dat_lg_N
col_form = split.data.frame(mean_by_formation_colony,  #by colonies
                            mean_by_formation_colony$formation) #zooids per formation
#just to look; max 328, smallest 19
col_form.n = lapply(col_form, function(x){dim(x)[1]})

#### RAREFACTION USING REIMANN DISTANCES ----
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
col_form.n #24 to 299 colonies
pop.size = sort(rep(c(2:500), repititions)) #make 500 so much bigger than max no. colonies
pop.dist <- lapply(pop.size, function (x){mvrnorm(n = x, 
                                                  rep(0,8), # mean 0 for 8 traits 
                                                  G_ext[[1]])}) # Simulate individual trait data based on the first VCOV matrix in our time series and different sample sizes. 
# Calculate VCOV matrices 
pop.cv <- lapply(pop.dist, cov) 

# Look at distribution of reimann distances
Reimann_result = c()

for(i in 1:length(pop.cv)){
    Reimann_result[i] = RiemannDist(pop.cv[[i]], G_ext[[1]])
}
    
out_results <- cbind(Reimann_result, pop.size)

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

diff.beta.gmax.df <- cbind(diff.beta.gmax.df, dist.gmat)

p.dist.gmat.b <- ggplot(diff.beta.gmax.df,
                        aes(x = as.numeric(mag.beta), y = as.numeric(dist.gmat))) + 
    geom_point() +
    geom_smooth(method = "lm",
                color = "#990000") +
    plot.theme +
    scale_x_continuous(expression(Magnitude~beta)) +
    scale_y_continuous(expression(distance~between~G~matrices~and~beta))

ggsave(p.dist.gmat.b, 
       file = "./Results/dist.g.mat.to.mag.beta.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(dist.gmat) ~ as.numeric(mag.beta),
           data = diff.beta.gmax.df[1:6,])) 
#nonsig at p = 0.2003; no relationship
