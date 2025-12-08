# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 29 Aug 2023

## The purpose of this script is to:
# 1. look at changes in E over time (i.e., does E expand in similar directions
# over time?)
# 2. assess the influence of E on direction of phenotypic change
# 3. test if E influences G over time (i.e., role in plasticity in directing
# evolution; should not affect it since underlying genetics is the same)

## the outputs are:
# -  comparison of adjacent # matrices through time (table and graph)
# -  comparison of estimated evolvability from E and Global G to 
# observed evolvability (table and graph)
# - comparison of direction of ∆z compared to Emax (table and graph)

#### LOAD DATA ----
source("./Scripts/0-env.R")

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

load(file = "./Results/size.list.RData")
P_size <- size.list[[1]]
G_size <- size.list[[2]]

load(file = "./Results/evol.diff.list.RData") #load the g matrices calculated above 

G_ext_NKLS <- evol.diff.list[[1]]
G_ext_NKBS <- evol.diff.list[[2]]
G_ext_tewk <- evol.diff.list[[3]]
G_ext_uki <- evol.diff.list[[4]]
G_ext_tai <- evol.diff.list[[5]]
G_ext_SHCSBSB <- evol.diff.list[[6]]
G_ext_mod <- evol.diff.list[[7]]
Gmax_NKLS_norm <- evol.diff.list[[8]]
Gmax_NKBS <- evol.diff.list[[9]]
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

#### CALCULATE E ----
Emat <- list()

for(i in 1:length(formation_list)){
    Emat[[i]] <- as.matrix(Pmat[[i]]) - as.matrix(Gmat[[i]])
}

names(Emat) = names(by_form.n)

##### E_ext -----

lapply(Emat, isSymmetric)  #is.symmetric.matrix
e.variances = lapply(Emat, diag)
paste("Trait variances")
head(e.variances)

e.eig_variances = lapply(Emat, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(e.eig_variances)

e.eig_percent = lapply(e.eig_variances, function (x) {x/sum(x)})
e.eig_per_mat = do.call(rbind, e.eig_percent)
e.eig_per_mat = data.frame(e.eig_per_mat, rownames(e.eig_per_mat))
e.eig_per = melt(e.eig_per_mat)
#dev.off()
E_PC_dist = ggplot(e.eig_per,
                   aes(x = variable, y = value,
                       group = rownames.e.eig_per_mat.,
                       colour = rownames.e.eig_per_mat.)) +
    geom_line(aes(linetype = rownames.e.eig_per_mat.)) +
    geom_point() +
    xlab("Principal component rank") +
    ylab("%Variation in the PC")
E_PC_dist #Tainui negative; none above 1; Upper Kai-Iwi looks FUNKY
## Waipuru 5th dim is wild; may only ret 4 dim?

E_ext = lapply(Emat, function (x){ ExtendMatrix(x, ret.dim = 5)$ExtMat}) #not 8 because last eigen value (#6) was negative
#ignore warning from above
#make 5th dim 0 if neg
lapply(E_ext, isSymmetric)  

E_ext_NKLS = round(as.matrix(E_ext[[1]]), 6) # The E matrix estimated for sample/formation 1
E_ext_NKBS = round(as.matrix(E_ext[[2]]), 6) # The E matrix estimated for sample/formation 2
E_ext_tewk = round(as.matrix(E_ext[[3]]), 6) # The E matrix estimated for sample/formation 3
E_ext_uki = round(as.matrix(E_ext[[4]]), 6) # The E matrix estimated for sample/formation 5
E_ext_tai = round(as.matrix(E_ext[[5]]), 6) # The E matrix estimated for sample/formation 6
E_ext_SHCSBSB = round(as.matrix(E_ext[[6]]), 6) # The E matrix estimated for sample/formation 7
E_ext_mod = round(as.matrix(E_ext[[7]]), 6) # The E matrix estimated for sample/formation 7

is.symmetric.matrix(E_ext_NKLS)
is.positive.definite(E_ext_NKLS)

is.symmetric.matrix(E_ext_NKBS)
is.positive.definite(E_ext_NKBS)

is.symmetric.matrix(E_ext_tewk)
is.positive.definite(E_ext_tewk)

is.symmetric.matrix(E_ext_uki)
is.positive.definite(E_ext_uki)

is.symmetric.matrix(E_ext_tai)
is.positive.definite(E_ext_tai)

is.symmetric.matrix(E_ext_SHCSBSB)
is.positive.definite(E_ext_SHCSBSB)

is.symmetric.matrix(E_ext_mod)
is.positive.definite(E_ext_mod)

#### Emax ----
Emax_NKLS <- eigen(E_ext_NKLS)$vectors[,1]
Emax_NKBS <- eigen(E_ext_NKBS)$vectors[,1]
Emax_tewk <- eigen(E_ext_tewk)$vectors[,1]
Emax_uki <- eigen(E_ext_uki)$vectors[,1]
Emax_tai <- eigen(E_ext_tai)$vectors[,1]
Emax_SHCSBSB <- eigen(E_ext_SHCSBSB)$vectors[,1]
Emax_mod <- eigen(E_ext_mod)$vectors[,1]

# Put Gmax to norm length
Emax_NKLS_norm <- f.normalize_vector(Emax_NKLS)
Emax_NKBS_norm <- f.normalize_vector(Emax_NKBS)
Emax_tewk_norm <- f.normalize_vector(Emax_tewk)
Emax_uki_norm <- f.normalize_vector(Emax_uki)
Emax_tai_norm <- f.normalize_vector(Emax_tai)
Emax_SHCSBSB_norm <- f.normalize_vector(Emax_SHCSBSB)
Emax_mod_norm <- f.normalize_vector(Emax_mod)

###### OBSERVED EVOLVABILITY ------
### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
observed_evolvability_in_direction_of_change_t1_E <- t(evolved_difference_unit_length_t1)%*%as.matrix(E_ext_NKLS)%*%evolved_difference_unit_length_t1
observed_evolvability_in_direction_of_change_t2_E <- t(evolved_difference_unit_length_t2)%*%as.matrix(E_ext_NKBS)%*%evolved_difference_unit_length_t2
observed_evolvability_in_direction_of_change_t3_E <- t(evolved_difference_unit_length_t3)%*%as.matrix(E_ext_tewk)%*%evolved_difference_unit_length_t3
observed_evolvability_in_direction_of_change_t4_E <- t(evolved_difference_unit_length_t4)%*%as.matrix(E_ext_uki)%*%evolved_difference_unit_length_t4
observed_evolvability_in_direction_of_change_t5_E <- t(evolved_difference_unit_length_t5)%*%as.matrix(E_ext_tai)%*%evolved_difference_unit_length_t5
observed_evolvability_in_direction_of_change_t6_E <- t(evolved_difference_unit_length_t6)%*%as.matrix(E_ext_SHCSBSB)%*%evolved_difference_unit_length_t6

###### OBSERVED CONDITIONAL EVOLVABILITY ------
### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
observed_conditional_evolvability_in_direction_of_change_t1_E <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(E_ext_NKLS))%*%evolved_difference_unit_length_t1)
observed_conditional_evolvability_in_direction_of_change_t2_E <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(E_ext_NKBS))%*%evolved_difference_unit_length_t2)
observed_conditional_evolvability_in_direction_of_change_t3_E <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(E_ext_tewk))%*%evolved_difference_unit_length_t3)
observed_conditional_evolvability_in_direction_of_change_t4_E <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(E_ext_uki))%*%evolved_difference_unit_length_t4)
observed_conditional_evolvability_in_direction_of_change_t5_E <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(E_ext_tai))%*%evolved_difference_unit_length_t5)
observed_conditional_evolvability_in_direction_of_change_t6_E <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(E_ext_SHCSBSB))%*%evolved_difference_unit_length_t6)

### Generate 10,000 selection gradients in random directions in the n-dimensional space
n_dimensions <- 8 # number of traits in G matrix
Beta <- randomBeta(10000, n_dimensions)

#outputs e, r, c, a, i
#e = evolvability
#r = respondability
#c = conditional evolvability
#a = autonomy of each selection gradient
#i = integration
#Beta = matrix of selection gradients
#e and c are calculating variances of means; should not be negative
#conditional must be equal to or smaller than e; often much small

###### ESTIMATED CONDITIONAL EVOLVABILITY & EVOLVABILITY ------

# Compute the mean, minimum and maximum evolvability (e_mean, e_min, e_max) for a G matrix based on 10,000 random selection gradients
X_t1_E <- evolvabilityBeta(as.matrix(E_ext_NKLS), Beta)
sumX_t1_E <- summary(X_t1_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t2_E <- evolvabilityBeta(as.matrix(E_ext_NKBS), Beta)
sumX_t2_E <- summary(X_t2_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t3_E <- evolvabilityBeta(as.matrix(E_ext_tewk), Beta)
sumX_t3_E <- summary(X_t3_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t4_E <- evolvabilityBeta(as.matrix(E_ext_uki), Beta)
sumX_t4_E <- summary(X_t4_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t5_E <- evolvabilityBeta(as.matrix(E_ext_tai), Beta)
sumX_t5_E <- summary(X_t5_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t6_E <- evolvabilityBeta(as.matrix(E_ext_SHCSBSB), Beta)
sumX_t6_E <- summary(X_t6_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t7_E <- evolvabilityBeta(as.matrix(E_ext_mod), Beta)
sumX_t7_E <- summary(X_t7_E) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_sum_E <- data.frame(c.mean = c(sumX_t1_E$Averages[[3]], sumX_t2_E$Averages[[3]], sumX_t3_E$Averages[[3]],
                               sumX_t4_E$Averages[[3]], sumX_t5_E$Averages[[3]], sumX_t6_E$Averages[[3]], 
                               sumX_t7_E$Averages[[3]]),
                    c.min = c(sumX_t1_E$Minimum[[3]], sumX_t2_E$Minimum[[3]], sumX_t3_E$Minimum[[3]],
                              sumX_t4_E$Minimum[[3]], sumX_t5_E$Minimum[[3]], sumX_t6_E$Minimum[[3]],
                              sumX_t7_E$Minimum[[3]]),
                    c.max = c(sumX_t1_E$Maximum[[3]], sumX_t2_E$Maximum[[3]], sumX_t3_E$Maximum[[3]],
                              sumX_t4_E$Maximum[[3]], sumX_t5_E$Maximum[[3]], sumX_t6_E$Maximum[[3]],
                              sumX_t7_E$Maximum[[3]]),
                    e.mean = c(sumX_t1_E$Averages[[1]], sumX_t2_E$Averages[[1]], sumX_t3_E$Averages[[1]],
                               sumX_t4_E$Averages[[1]], sumX_t5_E$Averages[[1]], sumX_t6_E$Averages[[1]],
                               sumX_t7_E$Averages[[1]]),
                    e.min = c(sumX_t1_E$Minimum[[1]], sumX_t2_E$Minimum[[1]], sumX_t3_E$Minimum[[1]],
                              sumX_t4_E$Minimum[[1]], sumX_t5_E$Minimum[[1]], sumX_t6_E$Minimum[[1]],
                              sumX_t7_E$Minimum[[1]]),
                    e.max = c(sumX_t1_E$Maximum[[1]], sumX_t2_E$Maximum[[1]], sumX_t3_E$Maximum[[1]],
                              sumX_t4_E$Maximum[[1]], sumX_t5_E$Maximum[[1]], sumX_t6_E$Maximum[[1]],
                              sumX_t7_E$Maximum[[1]]),
                    observed_e = c(observed_evolvability_in_direction_of_change_t1_E,
                                   observed_evolvability_in_direction_of_change_t2_E,
                                   observed_evolvability_in_direction_of_change_t3_E,
                                   observed_evolvability_in_direction_of_change_t4_E,
                                   observed_evolvability_in_direction_of_change_t5_E,
                                   observed_evolvability_in_direction_of_change_t6_E,
                                   ""),
                    observed_c = c(observed_conditional_evolvability_in_direction_of_change_t1_E,
                                   observed_conditional_evolvability_in_direction_of_change_t2_E,
                                   observed_conditional_evolvability_in_direction_of_change_t3_E,
                                   observed_conditional_evolvability_in_direction_of_change_t4_E,
                                   observed_conditional_evolvability_in_direction_of_change_t5_E,
                                   observed_conditional_evolvability_in_direction_of_change_t6_E,
                                   ""),
                    row.names = levels(formation_list))
#NO NEGATIVE VALUES!
#SAME AS WHAT WE SEE FOR G

write.csv(X_sum_E,
          "./Results/evolvability.summary_E.csv")

## PLOT
X_sum_E$formation <- rownames(X_sum_E)
X_sum_E$formation <- factor(X_sum_E$formation, levels = c("NKLS", "NKBS",
                                                      "Tewkesbury", "Upper Kai-Iwi",
                                                      "Tainui", "SHCSBSB", "modern"))

X_sum_E$form.trans <- formation_transition
X_sum_E$form.trans <- factor(X_sum_E$form.trans,
                           levels = c("NKLS to NKBS", 
                                      "NKBS to Tewkesbury",
                                      "Tewkesbury to Upper Kai-Iwi",
                                      "Upper Kai-Iwi to Tainui", 
                                      "Tainui to SHCSBSB",
                                      "SHCSBSB to modern", ""))

X_sum_E.trim <- X_sum_E[1:6,]
e.evol <- ggplot(X_sum_E.trim, aes(x = form.trans)) +
    geom_boxplot(aes(ymin = e.min, 
                     lower = e.min,
                     middle = e.mean,
                     ymax = e.max,
                     upper = e.max,
                     fill = "gray"),
                 stat = "identity", fill = "gray") +
    geom_point(aes(y = as.numeric(observed_e),
                   color = "black"),
               size = 5, shape = 17, color = "black") +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Evolvability") +
    plot.theme
#transition to modern is different

ggsave(e.evol, 
       file = "./Results/evolvability_E.png", 
       width = 14, height = 10, units = "cm")

# By comparing the evolvabilities you estimated in the direction of change (lines 9 and 12) with the average evolvabilities calculated by running line 20, you get a sense of whether evolution happened in directions with above or below average evolvability.  

##### CHANGE IN EMAX BETWEEN FORMATIONS -----

### Proportion of variance in n-dimensional trait space that is explained by PC1 (i.e., the first eigenvector)
#eigen(as.matrix(G_matrix_1))$values[1]/sum(eigen(as.matrix(G_matrix_1))$values)
eigen(as.matrix(E_ext_NKLS))$values[1]/sum(eigen(as.matrix(E_ext_NKLS))$values) #0.4749052
eigen(as.matrix(E_ext_NKBS))$values[1]/sum(eigen(as.matrix(E_ext_NKBS))$values) #0.4746287
eigen(as.matrix(E_ext_tewk))$values[1]/sum(eigen(as.matrix(E_ext_tewk))$values) #0.4843018
eigen(as.matrix(E_ext_uki))$values[1]/sum(eigen(as.matrix(E_ext_uki))$values) #0.4309931
eigen(as.matrix(E_ext_tai))$values[1]/sum(eigen(as.matrix(E_ext_tai))$values) #0.451781
eigen(as.matrix(E_ext_SHCSBSB))$values[1]/sum(eigen(as.matrix(E_ext_SHCSBSB))$values) #0.4146446
eigen(as.matrix(E_ext_mod))$values[1]/sum(eigen(as.matrix(E_ext_mod))$values) #0.5180123

### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 
Emax_NKLS <- eigen(E_ext_NKLS)$vectors[,1]
Emax_NKBS <- eigen(E_ext_NKBS)$vectors[,1]
Emax_tewk <- eigen(E_ext_tewk)$vectors[,1]
Emax_uki <- eigen(E_ext_uki)$vectors[,1]
Emax_tai <- eigen(E_ext_tai)$vectors[,1]
Emax_SHCSBSB <- eigen(E_ext_SHCSBSB)$vectors[,1]
Emax_mod <- eigen(E_ext_mod)$vectors[,1]

# Put Gmax to norm length
Emax_NKLS_norm <- f.normalize_vector(Emax_NKLS)
Emax_NKBS_norm <- f.normalize_vector(Emax_NKBS)
Emax_tewk_norm <- f.normalize_vector(Emax_tewk)
Emax_uki_norm <- f.normalize_vector(Emax_uki)
Emax_tai_norm <- f.normalize_vector(Emax_tai)
Emax_SHCSBSB_norm <- f.normalize_vector(Emax_SHCSBSB)
Emax_mod_norm <- f.normalize_vector(Emax_mod)

# Calculate the dot product of the unit vectors; tells number 0 to 1
dot_product.Emax_NKLS_NKBS <- sum(Emax_NKLS_norm * Emax_NKBS_norm) #0.9980535
# Calculate the angle in radians
angle_radians.Emax_NKLS_NKBS <- acos(dot_product.Emax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Emax_NKLS_NKBS <- angle_radians.Emax_NKLS_NKBS * (180 / pi)
#3.575492

dot_product.Emax_NKBS_tewk <- sum(Emax_NKBS_norm * Emax_tewk_norm) #0.9928007
# Calculate the angle in radians
angle_radians.Emax_NKBS_tewk <- acos(dot_product.Emax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Emax_NKBS_tewk <- angle_radians.Emax_NKBS_tewk * (180 / pi)
#6.879275

dot_product.Emax_tewk_uki <- sum(Emax_tewk_norm * Emax_uki_norm) #0.9638897
# Calculate the angle in radians
angle_radians.Emax_tewk_uki <- acos(dot_product.Emax_tewk_uki)
# Convert the angle to degrees
angle_degrees.Emax_tewk_uki <- angle_radians.Emax_tewk_uki * (180 / pi)
#15.44433

dot_product.Emax_uki_tai <- sum(Emax_uki_norm * Emax_tai_norm) #0.9926203
# Calculate the angle in radians
angle_radians.Emax_uki_tai <- acos(dot_product.Emax_uki_tai)
# Convert the angle to degrees
angle_degrees.Emax_uki_tai <- angle_radians.Emax_uki_tai * (180 / pi)
#6.965062

dot_product.Emax_tai_SHCSBSB <- sum(Emax_tai_norm * Emax_SHCSBSB_norm) #0.9689412
# Calculate the angle in radians
angle_radians.Emax_tai_SHCSBSB <- acos(dot_product.Emax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Emax_tai_SHCSBSB <- angle_radians.Emax_tai_SHCSBSB * (180 / pi)
#14.31728

dot_product.Emax_SHCSBSB_mod <- sum(Emax_SHCSBSB_norm * Emax_mod_norm) #0.9568188
# Calculate the angle in radians
angle_radians.Emax_SHCSBSB_mod <- acos(dot_product.Emax_SHCSBSB_mod)
# Convert the angle to degrees
angle_degrees.Emax_SHCSBSB_mod <- angle_radians.Emax_SHCSBSB_mod * (180 / pi)
#16.89897

#SUPER SIMILAR THROUGH TIME

corr.diff_Es <- c(dot_product.Emax_NKLS_NKBS, dot_product.Emax_NKBS_tewk,
                  dot_product.Emax_tewk_uki, dot_product.Emax_uki_tai,
                  dot_product.Emax_tai_SHCSBSB, dot_product.Emax_SHCSBSB_mod)

angle_diff_Es <- c(angle_degrees.Emax_NKLS_NKBS, angle_degrees.Emax_NKBS_tewk,
                   angle_degrees.Emax_tewk_uki, angle_degrees.Emax_uki_tai, 
                   angle_degrees.Emax_tai_SHCSBSB, angle_degrees.Emax_SHCSBSB_mod)

diff_between_Es <- as.data.frame(cbind(angle_diff_Es, corr.diff_Es))
diff_between_Es$angle_diff_Es <- as.numeric(diff_between_Es$angle_diff_Es)
diff_between_Es$angle.diff_Es.time <- formation_transition[-7]
diff_between_Es$angle.diff_Es.time <- factor(diff_between_Es$angle.diff_Es.time,
                                             levels = c("NKLS to NKBS", 
                                                        "NKBS to Tewkesbury",
                                                        "Tewkesbury to Upper Kai-Iwi",
                                                        "Upper Kai-Iwi to Tainui", 
                                                        "Tainui to SHCSBSB",
                                                        "SHCSBSB to modern"))

for(i in 1:nrow(diff_between_Es)){
    if(isTRUE(diff_between_Es$angle_diff_Es[i] > 90)){
        diff_between_Es$angle_diff_Es[i] <- 180 - diff_between_Es$angle_diff_Es[i]
    }
    else{
        next
    }
}

write.csv(diff_between_Es,
          "./Results/differences.between.Es.csv",
          row.names = FALSE)

p.ang_e <- ggplot(diff_between_Es) +
    geom_point(aes(x = angle.diff_Es.time, y = angle_diff_Es),
               size = 5, shape = 15) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between E matrices", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.ang_e, 
       file = "./Results/angle.e.diff.png", 
       width = 20, height = 20, units = "cm")

diff_between_Es$corr.diff_Es <- as.numeric(diff_between_Es$corr.diff_Es)
p.dot.prod_e <- ggplot(diff_between_Es) +
    geom_point(aes(x = angle.diff_Es.time, y = abs(corr.diff_Es)),
               size = 5, shape = 15) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Dot prodcut between E matrices") + 
    plot.theme

ggsave(p.dot.prod_e, 
       file = "./Results/dot.prod.e.diff.png", 
       width = 20, height = 20, units = "cm")

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO EMAX -----

### See if change is in direction of E max
## use Emax of t1 and compare to ∆z
# Calculate the dot product of the unit vectors
dot_product.Emax_NKLS <- sum(Emax_NKLS_norm * evolved_difference_unit_length_t1) #-0.08420543
# Calculate the angle in radians
angle_radians.Emax_NKLS <- acos(dot_product.Emax_NKLS)
# Convert the angle to degrees
angle_degrees.Emax_NKLS <- angle_radians.Emax_NKLS * (180 / pi)
#94.83034; 85.16966

# Calculate the dot product of the unit vectors
dot_product.Emax_NKBS <- sum(Emax_NKBS_norm * evolved_difference_unit_length_t2) #0.6802636
# Calculate the angle in radians
angle_radians.Emax_NKBS <- acos(dot_product.Emax_NKBS)
# Convert the angle to degrees
angle_degrees.Emax_NKBS <- angle_radians.Emax_NKBS * (180 / pi)
#47.13575

# Calculate the dot product of the unit vectors
dot_product.Emax_tewk <- sum(Emax_tewk_norm * evolved_difference_unit_length_t3) #-0.7823083
# Calculate the angle in radians
angle_radians.Emax_tewk <- acos(dot_product.Emax_tewk)
# Convert the angle to degrees
angle_degrees.Emax_tewk <- angle_radians.Emax_tewk * (180 / pi)
#141.4724; 38.5276

# Calculate the dot product of the unit vectors
dot_product.Emax_uki <- sum(Emax_uki_norm * evolved_difference_unit_length_t4) #-0.6355665
# Calculate the angle in radians
angle_radians.Emax_uki <- acos(dot_product.Emax_uki)
# Convert the angle to degrees
angle_degrees.Emax_uki <- angle_radians.Emax_uki * (180 / pi)
#129.462; 50.538

# Calculate the dot product of the unit vectors
dot_product.Emax_tai <- sum(Emax_tai_norm * evolved_difference_unit_length_t5) #0.5323496
# Calculate the angle in radians
angle_radians.Emax_tai <- acos(dot_product.Emax_tai)
# Convert the angle to degrees
angle_degrees.Emax_tai <- angle_radians.Emax_tai * (180 / pi)
#57.83565

# Calculate the dot product of the unit vectors
dot_product.Emax_SHCSBSB <- sum(Emax_SHCSBSB_norm * evolved_difference_unit_length_t6) #0.1301264
# Calculate the angle in radians
angle_radians.Emax_SHCSBSB <- acos(dot_product.Emax_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Emax_SHCSBSB <- angle_radians.Emax_SHCSBSB * (180 / pi)
#82.5231

corr.diff_Emax_to_z <- c(dot_product.Emax_NKLS, dot_product.Emax_NKBS,
                         dot_product.Emax_tewk, dot_product.Emax_uki,
                         dot_product.Emax_tai, dot_product.Emax_SHCSBSB)

angle_diff_Emax_to_z <- c(angle_degrees.Emax_NKLS, angle_degrees.Emax_NKBS,
                          angle_degrees.Emax_tewk, angle_degrees.Emax_uki, 
                          angle_degrees.Emax_tai, angle_degrees.Emax_SHCSBSB)

diff_between_Emax_z <- as.data.frame(cbind(angle_diff_Emax_to_z, corr.diff_Emax_to_z))
diff_between_Emax_z$angle.diff.time <- formation_transition[-7]
diff_between_Emax_z$angle.diff.time <- factor(diff_between_Emax_z$angle.diff.time,
                                              levels = c("NKLS to NKBS", 
                                                         "NKBS to Tewkesbury",
                                                         "Tewkesbury to Upper Kai-Iwi",
                                                         "Upper Kai-Iwi to Tainui", 
                                                         "Tainui to SHCSBSB",
                                                         "SHCSBSB to modern"))
colnames(diff_between_Emax_z) <- c("angle_diff_Emax_to_z", "corr.diff_Emax_to_z", "angle.diff.time")
diff_between_Emax_z$angle_diff_Emax_to_z <- as.numeric(diff_between_Emax_z$angle_diff_Emax_to_z)



for(i in 1:nrow(diff_between_Emax_z)){
    if(isTRUE(diff_between_Emax_z$angle_diff_Emax_to_z[i] > 90)){
        diff_between_Emax_z$angle_diff_Emax_to_z[i] <- 180 - as.numeric(diff_between_Emax_z$angle_diff_Emax_to_z[i])
    }
    else{
        next
    }
}

write.csv(diff_between_Emax_z,
          "./Results/differences.between.Emax.z.csv",
          row.names = FALSE)

p.e.pheno <- ggplot(diff_between_Emax_z) +
    geom_point(aes(x = angle.diff.time, y = angle_diff_Emax_to_z),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between E and ∆z", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.e.pheno, 
       file = "./Results/angle.e.pheno.png", 
       width = 14, height = 12, units = "cm")

p.dot.prod_emax <- ggplot(data = diff_between_Emax_z) +
    geom_point(aes(x = angle.diff.time, y = corr.diff_Emax_to_z),
               size = 5, shape = 15) + 
    plot.theme +
    #scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Dot Product of ∆Z from Emax")

ggsave(p.dot.prod_emax, 
       file = "./Results/corr.emax.p.diff.png", 
       width = 14, height = 10, units = "cm")

#### AMOUNT OF VARIATION E EXPLAINS ----
E_size <- c(sum(diag(E_ext_NKBS)),
            sum(diag(E_ext_NKLS)),
            sum(diag(E_ext_tewk)),
            sum(diag(E_ext_uki)),
            sum(diag(E_ext_tai)),
            sum(diag(E_ext_SHCSBSB)),
            sum(diag(E_ext_mod)))

#% not E
G_size/P_size
mean(G_size/P_size) #on avg 36% is not E; so 64% is E

#what explains E?
#temperature is in form.meta
#time averaging is calculated here
form.meta$amt.time = as.numeric(form.meta$Start_age)-as.numeric(form.meta$End_age)

##environment
form.meta$tract = c("TST-MCS", #nkls
                    "TST-MCS", #nkbs
                    "TST-MCS", #tewk 
                    "", #wai
                    "TST-MCS", #uki
                    "TST-MCS", #tai
                    "RST", #schsbsbs
                    "") #mod
#maybe not that useful

#add in shell side; doesn't make sense because of matrixes, unless it was a variable in its estimation?
#how do I think about this since I also used locality as uncertainty?

#for now, let's look at temp and time averaging
##1. REIMANN distance
e.dist.NKLS_NKBS <- RiemannDist(E_ext_NKBS, E_ext_NKLS)
e.dist.NKBS_tewk <- RiemannDist(E_ext_NKLS, E_ext_tewk)
e.dist.tewk_uki <- RiemannDist(E_ext_tewk, E_ext_uki)
e.dist.uki_tai <- RiemannDist(E_ext_uki, E_ext_tai)
e.dist.tai_SHCSBSB <- RiemannDist(E_ext_tai, E_ext_SHCSBSB)
e.dist.SHCSBSB_mod <- RiemannDist(E_ext_SHCSBSB, E_ext_mod)

dist.Emat <- c(e.dist.NKLS_NKBS, e.dist.NKBS_tewk,
               e.dist.tewk_uki, e.dist.uki_tai, 
               e.dist.tai_SHCSBSB, e.dist.SHCSBSB_mod)

diff.temp <- c(form.meta$temp[form.meta$formationName == "Nukumaru Brown Sand"] - form.meta$temp[form.meta$formationName == "Nukumaru Limestone"],
               form.meta$temp[form.meta$formationName == "Tewkesbury"] - form.meta$temp[form.meta$formationName == "Nukumaru Brown Sand"],
               form.meta$temp[form.meta$formationName == "Upper Kai-Iwi Shell Bed"] - form.meta$temp[form.meta$formationName == "Tewkesbury"],
               form.meta$temp[form.meta$formationName == "Tainui Shell Bed"] - form.meta$temp[form.meta$formationName == "Upper Kai-Iwi Shell Bed"],
               form.meta$temp[form.meta$formationName == "Shakespeare Basal Cliff Sand"] - form.meta$temp[form.meta$formationName == "Tainui Shell Bed"],
               form.meta$temp[form.meta$formationName == "modern"] - form.meta$temp[form.meta$formationName == "Shakespeare Basal Cliff Sand"])

formation_transition <- c("NKLS to NKBS", 
                          "NKBS to Tewkesbury",
                          "Tewkesbury to Upper Kai-Iwi",
                          "Upper Kai-Iwi to Tainui", 
                          "Tainui to SHCSBSB",
                          "SHCSBSB to modern")

diff.temp.emax.df <- as.data.frame(cbind(diff.temp, dist.Emat, formation_transition))

p.dist.emat.t <- ggplot(diff.temp.emax.df,
                        aes(x = as.numeric(diff.temp), y = as.numeric(dist.Emat))) + 
    geom_point() +
    geom_smooth(method = "lm",
                color = "#990000") +
    plot.theme +
    scale_x_continuous(expression(Temperature~Difference)) +
    scale_y_continuous(expression(distance~between~E~matrices))

ggsave(p.dist.emat.t, 
       file = "./Results/dist.e.mat.to.temp.diff.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(dist.Emat) ~ as.numeric(diff.temp),
           data = diff.temp.emax.df)) 
#nonsig at p = 0.2528; no relationship

#look at PC1 of emat
e.pc1 <- c(e.eig_variances$NKLS[1],
           e.eig_variances$NKBS[1],
           e.eig_variances$Tewkesbury[1],
           e.eig_variances$`Upper Kai-Iwi`[1],
           e.eig_variances$Tainui[1],
           e.eig_variances$SHCSBSB[1],
           e.eig_variances$modern[1])

e.eig_per.1 <- c(e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "NKLS" & e.eig_per$variable == "X1"],
                 e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "NKBS" & e.eig_per$variable == "X1"],
                 e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "Tewkesbury" & e.eig_per$variable == "X1"],
                 e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "Upper Kai-Iwi" & e.eig_per$variable == "X1"],
                 e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "Tainui" & e.eig_per$variable == "X1"],
                 e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "SHCSBSB" & e.eig_per$variable == "X1"],
                 e.eig_per$value[e.eig_per$rownames.e.eig_per_mat. == "modern" & e.eig_per$variable == "X1"])

# look at PC1 over all PCs in Emat
forms <- c("NKLS", "NKBS", "Tewkesbury", "Upper Kai-Iwi", "Tainui", "SHCSBSB", "modern")
e.pc.time <- as.data.frame(cbind(forms, e.pc1, e.eig_per.1))

e.pc1.diff <- c(e.pc1[2]-e.pc1[1],
                e.pc1[3]-e.pc1[2],
                e.pc1[4]-e.pc1[3],
                e.pc1[5]-e.pc1[4],
                e.pc1[6]-e.pc1[5],
                e.pc1[7]-e.pc1[6])

e.eig.per.1.diff <- c(e.eig_per.1[2]-e.eig_per.1[1],
                      e.eig_per.1[3]-e.eig_per.1[2],
                      e.eig_per.1[4]-e.eig_per.1[3],
                      e.eig_per.1[5]-e.eig_per.1[4],
                      e.eig_per.1[6]-e.eig_per.1[5],
                      e.eig_per.1[7]-e.eig_per.1[6])

diff.temp.emax.df <- as.data.frame(cbind(diff.temp.emax.df, e.pc1.diff, e.eig.per.1.diff))

p.diff.epc1.t <- ggplot(diff.temp.emax.df,
                        aes(x = as.numeric(diff.temp), y = as.numeric(e.pc1.diff))) + 
    geom_point() +
    geom_smooth(method = "lm",
                color = "#990000") +
    plot.theme +
    scale_x_continuous(expression(Temperature~Difference)) +
    scale_y_continuous(expression(PC1~Difference))

ggsave(p.diff.epc1.t, 
       file = "./Results/diff.e.pc1.to.temp.diff.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(e.pc1.diff) ~ as.numeric(diff.temp),
           data = diff.temp.emax.df)) 
#nonsig at p = 0.2056; no relationship

p.diff.epc1.per.t <- ggplot(diff.temp.emax.df,
                        aes(x = as.numeric(diff.temp), y = as.numeric(e.eig.per.1.diff))) + 
    geom_point() +
    geom_smooth(method = "lm",
                color = "#990000") +
    plot.theme +
    scale_x_continuous(expression(Temperature~Difference)) +
    scale_y_continuous("%PC1 Difference")

ggsave(p.diff.epc1.per.t, 
       file = "./Results/diff.e.pc1.per.to.temp.diff.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(e.eig.per.1.diff) ~ as.numeric(diff.temp),
           data = diff.temp.emax.df)) 
#nonsig at p = 0.1885; no relationship



P_G_E_size <- as.data.frame(cbind(form.name, P_size, G_size, E_size))
P_G_E_size$form.name <- factor(P_G_E_size$form.name, levels = c("NKLS", "NKBS",
                                                                "Tewkesbury", 
                                                                "Upper Kai-Iwi",
                                                                "Tainui",
                                                                "SHCSBSB",
                                                                "modern"))

#proportion G to P
P_G_E_size$prop.e.p <- as.numeric(P_G_E_size$E_size)/as.numeric(P_G_E_size$P_size)
range(P_G_E_size$prop.e.p)
#52 to 71 %

P_G_E_size$prop.g.p <- as.numeric(P_G_E_size$G_size)/as.numeric(P_G_E_size$P_size)
range(P_G_E_size$prop.g.p)

p.pge.size <- ggplot(P_G_E_size) +
    geom_point(aes(y = as.numeric(P_size), x = form.name,
                   group = form.name, col = form.name),
               size = 5, shape = 17) + 
    geom_point(aes(y = as.numeric(G_size), x = form.name,
                   group = form.name, col = form.name),
               size = 5, shape = 15) + 
    geom_point(aes(y = as.numeric(E_size), x = form.name,
                   group = form.name, col = form.name),
               size = 5, shape = 16) + 
    plot.theme +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Sizes of P and E",
                       limits = c(0, .25)) +
    scale_color_manual(values = col.form) +
    ggtitle("Size (diagonals) of P (triangle), G (cirlce), and E (X)")
#scale_shape_manual(name = "Matrix",
#                   values = c("P matrix" = 17,
#                              "G matrix" = 4,
#                              "E matrix" = 16)) 

ggsave(p.pge.size, 
       file = "./Results/p.g.e.size.png", 
       width = 20, height = 20, units = "cm")

## P v G
p.pg.size <- ggplot(P_G_E_size) +
    geom_point(aes(y = as.numeric(P_size), x = as.numeric(G_size),
                   group = form.name, col = form.name),
               size = 5) + 
    plot.theme +
    scale_x_continuous(name = "G matrix size") +
    scale_y_continuous(name = "P matrix size") +
    scale_color_manual(values = col.form) +
    ggtitle("Size of P v G")

summary(lm(as.numeric(P_G_E_size$P_size) ~ as.numeric(P_G_E_size$G_size)))
#slope = 0.24089, p-value = 0.7833, r2 = 0

## look at only those with high sample size
keep <- c("NKLS", "NKBS", "Tewkesbury", "SHCSBSB")
pge.trim <- P_G_E_size[P_G_E_size$form.name %in% keep,]
col.form.trim <- c(col.form[1], col.form[2], col.form[3], col.form[7])

p.pg.size.trim <- ggplot(pge.trim) +
    geom_point(aes(y = as.numeric(P_size), x = as.numeric(G_size),
                   group = form.name, col = form.name),
               size = 5) + 
    plot.theme +
    scale_x_continuous(name = "G matrix size") + 
    scale_y_continuous(name = "P matrix size") +
    scale_color_manual(values = col.form.trim) +
    ggtitle("Size of P v G for formations with high sample sizes")

ggsave(p.pg.size.trim, 
       file = "./Results/p.g.size.high.samples.png", 
       width = 20, height = 20, units = "cm")


## P v E
p.pe.size <- ggplot(P_G_E_size) +
    geom_point(aes(y = as.numeric(P_size), x = as.numeric(E_size),
                   group = form.name, col = form.name),
               size = 5, shape = 16) + 
    plot.theme +
    scale_x_continuous(name = "E matrix size") +
    scale_y_continuous(name = "P matrix size") +
    scale_color_manual(values = col.form) +
    ggtitle("Size of P v E")

summary(lm(as.numeric(P_G_E_size$P_size) ~ as.numeric(P_G_E_size$E_size)))
#slope = 0.89859, p-value = 0.006109, r2 = 0.7666
## almost sig!!

ggsave(p.pe.size, 
       file = "./Results/p.e.size.png", 
       width = 20, height = 20, units = "cm")

## look at only those with high sample size
keep <- c("NKLS", "NKBS", "Tewkesbury", "SHCSBSB")
pge.trim <- P_G_E_size[P_G_E_size$form.name %in% keep,]
col.form.trim <- c(col.form[1], col.form[2], col.form[3], col.form[7])

p.pe.size.trim <- ggplot(pge.trim) +
    geom_point(aes(y = as.numeric(P_size), x = as.numeric(E_size),
                   group = form.name, col = form.name),
               size = 5) + 
    plot.theme +
    scale_x_continuous(name = "E matrix size") + 
    scale_y_continuous(name = "P matrix size") +
    scale_color_manual(values = col.form.trim) +
    ggtitle("Size of P v E for formations with high sample sizes")

#### PROJECTING G & E ON P ----
#gather variance-covariance matrices
#E_ext_NKLS, etc.
#project
## X*E*X^T, where X is the matrix of eigen vectors as colums
#get PC values for all new matrices

p.eig_variances
#transform to columns rather than rows
#p.eig.trans <- lapply(p.eig_variances, function (x) t(t(x)))
#p.eig.trans.t <- lapply(p.eig_variances, function (x) t(x))

as.matrix(p.eig_variances$NKLS)%*%as.matrix(Gmat$NKLS)%*%t(as.matrix(p.eig_variances$NKLS))

Gmat.proj <- list()

###### PC VALUES -----



lapply(Emat, isSymmetric)  #is.symmetric.matrix
e.variances = lapply(Emat, diag)
paste("Trait variances")
head(e.variances)

e.eig_variances = lapply(Emat, function (x) {eigen(x)$values})
paste("Eigenvalue variances")
head(e.eig_variances)

e.eig_percent = lapply(e.eig_variances, function (x) {x/sum(x)})
e.eig_per_mat = do.call(rbind, e.eig_percent)
e.eig_per_mat = data.frame(e.eig_per_mat, rownames(e.eig_per_mat))
e.eig_per = melt(e.eig_per_mat)
#dev.off()
E_PC_dist = ggplot(e.eig_per,
                   aes(x = variable, y = value,
                       group = rownames.e.eig_per_mat.,
                       colour = rownames.e.eig_per_mat.)) +
    geom_line(aes(linetype = rownames.e.eig_per_mat.)) +
    geom_point() +
    xlab("Principal component rank") +
    ylab("%Variation in the PC")
E_PC_dist #Tainui negative; none above 1; Upper Kai-Iwi looks FUNKY
## Waipuru 5th dim is wild; may only ret 4 dim?
