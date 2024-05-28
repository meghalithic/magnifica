# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 29 Aug 2023

## The purpose of this script is to:
# 1. compare correlation of P to G
# 2. compare G through time
# 3. test if G changes in direction of above-average evolvability 
#    (i.e., conditional evolvability and evolvability)
# 4. test if G changes in direction of Gmax (angle)
# 5. test if Pmax changes in direction of Gmax (angle)
# 6. see how G changes through time
# 7. see how G changes in relation to temperature

## the outputs are:
# -  comparison of estimated evolvability from G and Global G to 
# observed evolvability (table and graph)
# - comparison of adjacent G matrices through time (table and graph)
# - comparison of direction of ∆z compared to Gmax (table and graph)
# - comparison of ∆z to temperature (linear model; graph)
# - calculation of beta (table and graph)
# - looking at P for G by redoing evolvability using P and comparing 
# effect of Pmax on ∆z (table and graph)

#### LOAD DATA ----

source("./Scripts/0-env.R")

load(file = "./Results/sum.data.list.w.modern.RData") #load the g matrices calculated above 
mean_by_formation <- sum.data.list[[1]]
mean_by_formation_colony <- sum.data.list[[2]]
means <- sum.data.list[[3]]

load(file = "./Results/data.list.w.modern.no.wai.RData") #load the g matrices calculated above 

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

load(file="./Results/global_ext.w.modern.no.wai.RData") #load the g matrices calculated above 
Glob_ext

load(file = "./Results/model_G.w.modern.no.wai.RData") #load the g matrices calculated above 
g.model <- model_G

#### CORR OF G & P ----

#Gmat
#Pmat

## RANDOM SKEWERS OF P & G OF EACH FORMATION
#NKLS
NKLS_comp_mat = RandomSkewers(list(G_ext[[1]], P_ext[[1]])) #need at least
NKLS_corr_mat = NKLS_comp_mat$correlations + t(NKLS_comp_mat$correlations) 
diag(NKLS_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(NKLS_corr_mat, upper = "number", lower = "pie")

#NKBS
NKBS_comp_mat = RandomSkewers(list(G_ext[[2]], P_ext[[2]])) #need at least
NKBS_corr_mat = NKBS_comp_mat$correlations + t(NKBS_comp_mat$correlations) 
diag(NKBS_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(NKBS_corr_mat, upper = "number", lower = "pie")

#Tewkesbury
Tewkesbury_comp_mat = RandomSkewers(list(G_ext[[3]], P_ext[[3]])) #need at least
Tewkesbury_corr_mat = Tewkesbury_comp_mat$correlations + t(Tewkesbury_comp_mat$correlations) 
diag(Tewkesbury_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Tewkesbury_corr_mat, upper = "number", lower = "pie")

#Upper Kai-Iwi
UKai_Iwi_comp_mat = RandomSkewers(list(G_ext[[4]], P_ext[[4]])) #need at least
UKai_Iwi_corr_mat = UKai_Iwi_comp_mat$correlations + t(UKai_Iwi_comp_mat$correlations) 
diag(UKai_Iwi_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(UKai_Iwi_corr_mat, upper = "number", lower = "pie")

#Tainui
Tainui_comp_mat = RandomSkewers(list(G_ext[[5]], P_ext[[5]])) #need at least
Tainui_corr_mat = Tainui_comp_mat$correlations + t(Tainui_comp_mat$correlations) 
diag(Tainui_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Tainui_corr_mat, upper = "number", lower = "pie")

#SHCSBSB
SHCSBSB_comp_mat = RandomSkewers(list(G_ext[[6]], P_ext[[6]])) #need at least
SHCSBSB_corr_mat = SHCSBSB_comp_mat$correlations + t(SHCSBSB_comp_mat$correlations) 
diag(SHCSBSB_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(SHCSBSB_corr_mat, upper = "number", lower = "pie")

#modern
modern_comp_mat = RandomSkewers(list(G_ext[[7]], P_ext[[7]])) #need at least
modern_corr_mat = modern_comp_mat$correlations + t(modern_comp_mat$correlations) 
diag(modern_corr_mat) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(modern_corr_mat, upper = "number", lower = "pie")

corr.p.g <- c(NKLS_corr_mat[1,2], NKBS_corr_mat[1,2],
              Tewkesbury_corr_mat[1,2], UKai_Iwi_corr_mat[1,2], 
              Tainui_corr_mat[1,2], SHCSBSB_corr_mat[1,2],
              modern_corr_mat[1,2])

corr.p.g.form <- cbind(levels(formation_list), corr.p.g)
colnames(corr.p.g.form) <- c("formation", "correlation")
corr.p.g.form

write.csv(corr.p.g.form,
          "Results/correlation.p.g.csv",
          row.names = FALSE)

#### POSITIVE DEFINITE ----
## check positive definite
# round to 10 decimals to make it symmetric
# to make it positive definite, use extended matrix to fill in

G_ext_NKLS = round(as.matrix(G_ext[[1]]), 6) # The G matrix estimated for sample/formation 1
G_ext_NKBS = round(as.matrix(G_ext[[2]]), 6) # The G matrix estimated for sample/formation 2
G_ext_tewk = round(as.matrix(G_ext[[3]]), 6) # The G matrix estimated for sample/formation 3
G_ext_uki = round(as.matrix(G_ext[[4]]), 6) # The G matrix estimated for sample/formation 5
G_ext_tai = round(as.matrix(G_ext[[5]]), 6) # The G matrix estimated for sample/formation 6
G_ext_SHCSBSB = round(as.matrix(G_ext[[6]]), 6) # The G matrix estimated for sample/formation 7
G_ext_mod = round(as.matrix(G_ext[[7]]), 6) # The G matrix estimated for sample/formation 7

is.symmetric.matrix(G_ext_NKLS)
is.positive.definite(G_ext_NKLS)

is.symmetric.matrix(G_ext_NKBS)
is.positive.definite(G_ext_NKBS)

is.symmetric.matrix(G_ext_tewk)
is.positive.definite(G_ext_tewk)

is.symmetric.matrix(G_ext_uki)
is.positive.definite(G_ext_uki)

is.symmetric.matrix(G_ext_tai)
is.positive.definite(G_ext_tai)

is.symmetric.matrix(G_ext_SHCSBSB)
is.positive.definite(G_ext_SHCSBSB)

is.symmetric.matrix(G_ext_mod)
is.positive.definite(G_ext_mod)

#### EVOLVABILITY ----

#trait means by time
mean_by_formation
#order of formations:
# "NKLS"  "NKBS" "Tewkesbury" "Upper Kai-Iwi" "Tainui" "SHCSBSB"

### Calculate the vector that defines the observed divergence between sample/formation 1 an 2
# A vector containing trait means from sample/formation
NKLS <- as.numeric(mean_by_formation[mean_by_formation == "NKLS", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")]) 
NKBS <- as.numeric(mean_by_formation[mean_by_formation == "NKBS", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")])
tewk <- as.numeric(mean_by_formation[mean_by_formation == "Tewkesbury", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")])
uki <- as.numeric(mean_by_formation[mean_by_formation == "Upper Kai-Iwi", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")])
tai <- as.numeric(mean_by_formation[mean_by_formation == "Tainui", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")])
SHCSBSB <- as.numeric(mean_by_formation[mean_by_formation == "SHCSBSB", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")])
mod <- as.numeric(mean_by_formation[mean_by_formation == "modern", c("avg.zh", "avg.mpw.b", "avg.cw.m", "avg.cw.d", "avg.ow.m", "avg.oh", "avg.o.side", "avg.c.side")])

#second - first
#t1 = NKBS - NKLS
#t2 = tewk - NKBS
#t3 = uki - tewk
#t4 = tai - uki
#t5 = SHCSBSB - tai
#t6 = mod - SHCSBSB

## really need to learn how to name things in functions...
evolved_difference_unit_length_t1 <- f.normalize_vector(NKBS - NKLS)
evolved_difference_unit_length_t2 <- f.normalize_vector(tewk - NKBS)
evolved_difference_unit_length_t3 <- f.normalize_vector(uki - tewk)
evolved_difference_unit_length_t4 <- f.normalize_vector(tai - uki)
evolved_difference_unit_length_t5 <- f.normalize_vector(SHCSBSB - tai)
evolved_difference_unit_length_t6 <- f.normalize_vector(mod - SHCSBSB)

###### OBSERVED EVOLVABILITY ------
### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
observed_evolvability_in_direction_of_change_t1 <- t(evolved_difference_unit_length_t1)%*%as.matrix(G_ext_NKLS)%*%evolved_difference_unit_length_t1
observed_evolvability_in_direction_of_change_t2 <- t(evolved_difference_unit_length_t2)%*%as.matrix(G_ext_NKBS)%*%evolved_difference_unit_length_t2
observed_evolvability_in_direction_of_change_t3 <- t(evolved_difference_unit_length_t3)%*%as.matrix(G_ext_tewk)%*%evolved_difference_unit_length_t3
observed_evolvability_in_direction_of_change_t4 <- t(evolved_difference_unit_length_t4)%*%as.matrix(G_ext_uki)%*%evolved_difference_unit_length_t4
observed_evolvability_in_direction_of_change_t5 <- t(evolved_difference_unit_length_t5)%*%as.matrix(G_ext_tai)%*%evolved_difference_unit_length_t5
observed_evolvability_in_direction_of_change_t6 <- t(evolved_difference_unit_length_t6)%*%as.matrix(G_ext_SHCSBSB)%*%evolved_difference_unit_length_t6

###### OBSERVED CONDITIONAL EVOLVABILITY ------
### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
observed_conditional_evolvability_in_direction_of_change_t1 <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(G_ext_NKLS))%*%evolved_difference_unit_length_t1)
observed_conditional_evolvability_in_direction_of_change_t2 <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(G_ext_NKBS))%*%evolved_difference_unit_length_t2)
observed_conditional_evolvability_in_direction_of_change_t3 <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(G_ext_tewk))%*%evolved_difference_unit_length_t3)
observed_conditional_evolvability_in_direction_of_change_t4 <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(G_ext_uki))%*%evolved_difference_unit_length_t4)
observed_conditional_evolvability_in_direction_of_change_t5 <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(G_ext_tai))%*%evolved_difference_unit_length_t5)
observed_conditional_evolvability_in_direction_of_change_t6 <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(G_ext_SHCSBSB))%*%evolved_difference_unit_length_t6)

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
X_t1 <- evolvabilityBeta(as.matrix(G_ext_NKLS), Beta)
sumX_t1 <- summary(X_t1) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t2 <- evolvabilityBeta(as.matrix(G_ext_NKBS), Beta)
sumX_t2 <- summary(X_t2) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t3 <- evolvabilityBeta(as.matrix(G_ext_tewk), Beta)
sumX_t3 <- summary(X_t3) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t4 <- evolvabilityBeta(as.matrix(G_ext_uki), Beta)
sumX_t4 <- summary(X_t4) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t5 <- evolvabilityBeta(as.matrix(G_ext_tai), Beta)
sumX_t5 <- summary(X_t5) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t6 <- evolvabilityBeta(as.matrix(G_ext_SHCSBSB), Beta)
sumX_t6 <- summary(X_t6) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t7 <- evolvabilityBeta(as.matrix(G_ext_mod), Beta)
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
                    row.names = levels(formation_list))
#NO NEGATIVE VALUES!

write.csv(X_sum,
          "./Results/evolvability.summary.csv")

## PLOT
X_sum$formation <- rownames(X_sum)
X_sum$formation <- factor(X_sum$formation, levels = c("NKLS", "NKBS",
                                                      "Tewkesbury", "Upper Kai-Iwi",
                                                      "Tainui", "SHCSBSB", "modern"))

X_sum$form.trans <- formation_transition
X_sum$form.trans <- factor(X_sum$form.trans,
                           levels = c("NKLS to NKBS", 
                                      "NKBS to Tewkesbury",
                                      "Tewkesbury to Upper Kai-Iwi",
                                      "Upper Kai-Iwi to Tainui", 
                                      "Tainui to SHCSBSB",
                                      "SHCSBSB to modern", ""))

X_sum.trim <- X_sum[1:6,]
p.evol <- ggplot(X_sum.trim, aes(x = form.trans)) +
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
        
ggsave(p.evol, 
       file = "./Results/evolvability.png", 
       width = 14, height = 10, units = "cm")

# By comparing the evolvabilities you estimated in the direction of change (lines 9 and 12) with the average evolvabilities calculated by running line 20, you get a sense of whether evolution happened in directions with above or below average evolvability.  

##### CHANGE IN GMAX BETWEEN FORMATIONS -----

### Proportion of variance in n-dimensional trait space that is explained by PC1 (i.e., the first eigenvector)
#eigen(as.matrix(G_matrix_1))$values[1]/sum(eigen(as.matrix(G_matrix_1))$values)
eigen(as.matrix(G_ext_NKLS))$values[1]/sum(eigen(as.matrix(G_ext_NKLS))$values) #0.4731959
eigen(as.matrix(G_ext_NKBS))$values[1]/sum(eigen(as.matrix(G_ext_NKBS))$values) #0.4819061
eigen(as.matrix(G_ext_tewk))$values[1]/sum(eigen(as.matrix(G_ext_tewk))$values) #0.4578517
eigen(as.matrix(G_ext_uki))$values[1]/sum(eigen(as.matrix(G_ext_uki))$values) #0.5058219
eigen(as.matrix(G_ext_tai))$values[1]/sum(eigen(as.matrix(G_ext_tai))$values) #0.3844157
eigen(as.matrix(G_ext_SHCSBSB))$values[1]/sum(eigen(as.matrix(G_ext_SHCSBSB))$values) #0.4479779
eigen(as.matrix(G_ext_mod))$values[1]/sum(eigen(as.matrix(G_ext_mod))$values) #0.5337849

### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 
Gmax_NKLS <- eigen(G_ext_NKLS)$vectors[,1]
Gmax_NKBS <- eigen(G_ext_NKBS)$vectors[,1]
Gmax_tewk <- eigen(G_ext_tewk)$vectors[,1]
Gmax_uki <- eigen(G_ext_uki)$vectors[,1]
Gmax_tai <- eigen(G_ext_tai)$vectors[,1]
Gmax_SHCSBSB <- eigen(G_ext_SHCSBSB)$vectors[,1]
Gmax_mod <- eigen(G_ext_mod)$vectors[,1]

# Put Gmax to norm length
Gmax_NKLS_norm <- f.normalize_vector(Gmax_NKLS)
Gmax_NKBS_norm <- f.normalize_vector(Gmax_NKBS)
Gmax_tewk_norm <- f.normalize_vector(Gmax_tewk)
Gmax_uki_norm <- f.normalize_vector(Gmax_uki)
Gmax_tai_norm <- f.normalize_vector(Gmax_tai)
Gmax_SHCSBSB_norm <- f.normalize_vector(Gmax_SHCSBSB)
Gmax_mod_norm <- f.normalize_vector(Gmax_mod)

# Calculate the dot product of the unit vectors; tells number 0 to 1
dot_product.Gmax_NKLS_NKBS <- sum(Gmax_NKLS_norm * Gmax_NKBS_norm) #0.9823321
# Calculate the angle in radians
angle_radians.Gmax_NKLS_NKBS <- acos(dot_product.Gmax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_NKBS <- angle_radians.Gmax_NKLS_NKBS * (180 / pi)
#10.78629

dot_product.Gmax_NKBS_tewk <- sum(Gmax_NKBS_norm * Gmax_tewk_norm) #0.9930473
# Calculate the angle in radians
angle_radians.Gmax_NKBS_tewk <- acos(dot_product.Gmax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_tewk <- angle_radians.Gmax_NKBS_tewk * (180 / pi)
#6.760309

dot_product.Gmax_tewk_uki <- sum(Gmax_tewk_norm * Gmax_uki_norm) #0.9590283
# Calculate the angle in radians
angle_radians.Gmax_tewk_uki <- acos(dot_product.Gmax_tewk_uki)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_uki <- angle_radians.Gmax_tewk_uki * (180 / pi)
#16.45787

dot_product.Gmax_uki_tai <- sum(Gmax_uki_norm * Gmax_tai_norm) #0.9656406
# Calculate the angle in radians
angle_radians.Gmax_uki_tai <- acos(dot_product.Gmax_uki_tai)
# Convert the angle to degrees
angle_degrees.Gmax_uki_tai <- angle_radians.Gmax_uki_tai * (180 / pi)
#15.06301

dot_product.Gmax_tai_SHCSBSB <- sum(Gmax_tai_norm * Gmax_SHCSBSB_norm) #0.7542782
# Calculate the angle in radians
angle_radians.Gmax_tai_SHCSBSB <- acos(dot_product.Gmax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Gmax_tai_SHCSBSB <- angle_radians.Gmax_tai_SHCSBSB * (180 / pi)
#41.03766

dot_product.Gmax_SHCSBSB_mod <- sum(Gmax_SHCSBSB_norm * Gmax_mod_norm) #0.1992655
# Calculate the angle in radians
angle_radians.Gmax_SHCSBSB_mod <- acos(dot_product.Gmax_SHCSBSB_mod)
# Convert the angle to degrees
angle_degrees.Gmax_SHCSBSB_mod <- angle_radians.Gmax_SHCSBSB_mod * (180 / pi)
#78.50599

corr.diff_Gs <- c(dot_product.Gmax_NKLS_NKBS, dot_product.Gmax_NKBS_tewk,
                  dot_product.Gmax_tewk_uki, dot_product.Gmax_uki_tai,
                  dot_product.Gmax_tai_SHCSBSB, dot_product.Gmax_SHCSBSB_mod)

angle_diff_Gs <- c(angle_degrees.Gmax_NKLS_NKBS, angle_degrees.Gmax_NKBS_tewk,
                   angle_degrees.Gmax_tewk_uki,
                   angle_degrees.Gmax_uki_tai, angle_degrees.Gmax_tai_SHCSBSB,
                   angle_degrees.Gmax_SHCSBSB_mod)

diff_between_Gs <- as.data.frame(cbind(angle_diff_Gs, corr.diff_Gs))
diff_between_Gs$angle_diff_Gs <- as.numeric(diff_between_Gs$angle_diff_Gs)
diff_between_Gs$angle.diff_Gs.time <- formation_transition[-7]
diff_between_Gs$angle.diff_Gs.time <- factor(diff_between_Gs$angle.diff_Gs.time,
                             levels = c("NKLS to NKBS", 
                                        "NKBS to Tewkesbury",
                                        "Tewkesbury to Upper Kai-Iwi",
                                        "Upper Kai-Iwi to Tainui", 
                                        "Tainui to SHCSBSB",
                                        "SHCSBSB to modern"))

for(i in 1:nrow(diff_between_Gs)){
    if(isTRUE(diff_between_Gs$angle_diff_Gs[i] > 90)){
        diff_between_Gs$angle_diff_Gs[i] <- 180 - diff_between_Gs$angle_diff_Gs[i]
    }
    else{
        next
    }
}

write.csv(diff_between_Gs,
          "./Results/differences.between.Gs.csv",
          row.names = FALSE)

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO GMAX -----

### See if change is in direction of G max
## use Gmax of t1 and compare to ∆z
# Calculate the dot product of the unit vectors
dot_product.Gmax_NKLS <- sum(Gmax_NKLS_norm * evolved_difference_unit_length_t1) #0.1673605
# Calculate the angle in radians
angle_radians.Gmax_NKLS <- acos(dot_product.Gmax_NKLS)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS <- angle_radians.Gmax_NKLS * (180 / pi)
#80.36561

# Calculate the dot product of the unit vectors
dot_product.Gmax_NKBS <- sum(Gmax_NKBS_norm * evolved_difference_unit_length_t2) #0.6870414
# Calculate the angle in radians
angle_radians.Gmax_NKBS <- acos(dot_product.Gmax_NKBS)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS <- angle_radians.Gmax_NKBS * (180 / pi)
#46.60364

# Calculate the dot product of the unit vectors
dot_product.Gmax_tewk <- sum(Gmax_tewk_norm * evolved_difference_unit_length_t3) #-0.9251647
# Calculate the angle in radians
angle_radians.Gmax_tewk <- acos(dot_product.Gmax_tewk)
# Convert the angle to degrees
angle_degrees.Gmax_tewk <- angle_radians.Gmax_tewk * (180 / pi)
#157.6932; 22.3068

# Calculate the dot product of the unit vectors
dot_product.Gmax_uki <- sum(Gmax_uki_norm * evolved_difference_unit_length_t4) #-0.7776185
# Calculate the angle in radians
angle_radians.Gmax_uki <- acos(dot_product.Gmax_uki)
# Convert the angle to degrees
angle_degrees.Gmax_uki <- angle_radians.Gmax_uki * (180 / pi)
#141.043; 38.957

# Calculate the dot product of the unit vectors
dot_product.Gmax_tai <- sum(Gmax_tai_norm * evolved_difference_unit_length_t5) #0.7512872
# Calculate the angle in radians
angle_radians.Gmax_tai <- acos(dot_product.Gmax_tai)
# Convert the angle to degrees
angle_degrees.Gmax_tai <- angle_radians.Gmax_tai * (180 / pi)
#41.298

# Calculate the dot product of the unit vectors
dot_product.Gmax_SHCSBSB <- sum(Gmax_SHCSBSB_norm * evolved_difference_unit_length_t6) #0.01569083
# Calculate the angle in radians
angle_radians.Gmax_SHCSBSB <- acos(dot_product.Gmax_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Gmax_SHCSBSB <- angle_radians.Gmax_SHCSBSB * (180 / pi)
#89.10094

corr.diff_Gmax_to_z <- c(dot_product.Gmax_NKLS, dot_product.Gmax_NKBS,
                    dot_product.Gmax_tewk, dot_product.Gmax_uki,
                    dot_product.Gmax_tai, dot_product.Gmax_SHCSBSB)

angle_diff_Gmax_to_z <- c(angle_degrees.Gmax_NKLS, angle_degrees.Gmax_NKBS,
                          angle_degrees.Gmax_tewk,
                          angle_degrees.Gmax_uki, angle_degrees.Gmax_tai,
                          angle_degrees.Gmax_SHCSBSB)

diff_between_Gmax_z <- as.data.frame(cbind(angle_diff_Gmax_to_z, corr.diff_Gmax_to_z))
diff_between_Gmax_z$angle.diff.time <- formation_transition[-7]
diff_between_Gmax_z$angle.diff.time <- factor(diff_between_Gmax_z$angle.diff.time,
                                              levels = c("NKLS to NKBS", 
                                                         "NKBS to Tewkesbury",
                                                         "Tewkesbury to Upper Kai-Iwi",
                                                         "Upper Kai-Iwi to Tainui", 
                                                         "Tainui to SHCSBSB",
                                                         "SHCSBSB to modern"))
colnames(diff_between_Gmax_z) <- c("angle_diff_Gmax_to_z", "corr.diff_Gmax_to_z", "angle.diff.time")
diff_between_Gmax_z$angle_diff_Gmax_to_z <- as.numeric(diff_between_Gmax_z$angle_diff_Gmax_to_z)



for(i in 1:nrow(diff_between_Gmax_z)){
    if(isTRUE(diff_between_Gmax_z$angle_diff_Gmax_to_z[i] > 90)){
        diff_between_Gmax_z$angle_diff_Gmax_to_z[i] <- 180 - as.numeric(diff_between_Gmax_P$angle_diff_Gmax_to_z[i])
    }
    else{
        next
    }
}

write.csv(diff_between_Gmax_z,
          "./Results/differences.between.Gmax.z.csv",
          row.names = FALSE)

p.g.pheno <- ggplot(diff_between_Gmax_z) +
    geom_point(aes(x = angle.diff.time, y = angle_diff_Gmax_to_z),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between G matrices", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.g.pheno, 
       file = "./Results/angle.g.pheno.png", 
       width = 14, height = 12, units = "cm")

###### PC2-4 compared to Gmax -----
## do for 2-4 as well
G2_NKLS <- eigen(G_ext_NKLS)$vectors[,2]
G2_NKBS <- eigen(G_ext_NKBS)$vectors[,2]
G2_tewk <- eigen(G_ext_tewk)$vectors[,2]
G2_uki <- eigen(G_ext_uki)$vectors[,2]
G2_tai <- eigen(G_ext_tai)$vectors[,2]
G2_SHCSBSB <- eigen(G_ext_SHCSBSB)$vectors[,2]
G2_mod <- eigen(G_ext_mod)$vectors[,2]

G3_NKLS <- eigen(G_ext_NKLS)$vectors[,3]
G3_NKBS <- eigen(G_ext_NKBS)$vectors[,3]
G3_tewk <- eigen(G_ext_tewk)$vectors[,3]
G3_uki <- eigen(G_ext_uki)$vectors[,3]
G3_tai <- eigen(G_ext_tai)$vectors[,3]
G3_SHCSBSB <- eigen(G_ext_SHCSBSB)$vectors[,3]
G3_mod <- eigen(G_ext_mod)$vectors[,3]

G4_NKLS <- eigen(G_ext_NKLS)$vectors[,4]
G4_NKBS <- eigen(G_ext_NKBS)$vectors[,4]
G4_tewk <- eigen(G_ext_tewk)$vectors[,4]
G4_uki <- eigen(G_ext_uki)$vectors[,4]
G4_tai <- eigen(G_ext_tai)$vectors[,4]
G4_SHCSBSB <- eigen(G_ext_SHCSBSB)$vectors[,4]
G4_mod <- eigen(G_ext_mod)$vectors[,4]

## normalize
G2_NKLS_norm <- f.normalize_vector(G2_NKLS)
G2_NKBS_norm <- f.normalize_vector(G2_NKBS)
G2_tewk_norm <- f.normalize_vector(G2_tewk)
G2_uki_norm <- f.normalize_vector(G2_uki)
G2_tai_norm <- f.normalize_vector(G2_tai)
G2_SHCSBSB_norm <- f.normalize_vector(G2_SHCSBSB)
G2_mod_norm <- f.normalize_vector(G2_mod)

G3_NKLS_norm <- f.normalize_vector(G3_NKLS)
G3_NKBS_norm <- f.normalize_vector(G3_NKBS)
G3_tewk_norm <- f.normalize_vector(G3_tewk)
G3_uki_norm <- f.normalize_vector(G3_uki)
G3_tai_norm <- f.normalize_vector(G3_tai)
G3_SHCSBSB_norm <- f.normalize_vector(G3_SHCSBSB)
G3_mod_norm <- f.normalize_vector(G3_mod)

G4_NKLS_norm <- f.normalize_vector(G4_NKLS)
G4_NKBS_norm <- f.normalize_vector(G4_NKBS)
G4_tewk_norm <- f.normalize_vector(G4_tewk)
G4_uki_norm <- f.normalize_vector(G4_uki)
G4_tai_norm <- f.normalize_vector(G4_tai)
G4_SHCSBSB_norm <- f.normalize_vector(G4_SHCSBSB)
G4_mod_norm <- f.normalize_vector(G4_mod)

## for 2-4 as well:
dot_product.G2_NKLS <- sum(G2_NKLS_norm * evolved_difference_unit_length_t1) #-0.9551744
angle_radians.G2_NKLS <- acos(dot_product.G2_NKLS)
angle_degrees.G2_NKLS <- angle_radians.G2_NKLS * (180 / pi)
#162.7799; 17.2201

dot_product.G2_NKBS <- sum(G2_NKBS_norm * evolved_difference_unit_length_t2) #-0.1888881
angle_radians.G2_NKBS <- acos(dot_product.G2_NKBS)
angle_degrees.G2_NKBS <- angle_radians.G2_NKBS * (180 / pi)
#100.8879; 79.1121

dot_product.G2_tewk <- sum(G2_tewk_norm * evolved_difference_unit_length_t3) #-0.159511
angle_radians.G2_tewk <- acos(dot_product.G2_tewk)
angle_degrees.G2_tewk <- angle_radians.G2_tewk * (180 / pi)
#99.17851; 80.82149

dot_product.G2_uki <- sum(G2_uki_norm * evolved_difference_unit_length_t4) #-0.3811011
angle_radians.G2_uki <- acos(dot_product.G2_uki)
angle_degrees.G2_uki <- angle_radians.G2_uki * (180 / pi)
#112.4019; 67.5981

dot_product.G2_tai <- sum(G2_tai_norm * evolved_difference_unit_length_t5) #0.2801589
angle_radians.G2_tai <- acos(dot_product.G2_tai)
angle_degrees.G2_tai <- angle_radians.G2_tai * (180 / pi)
#73.73031

dot_product.G2_SHCSBSB <- sum(G2_SHCSBSB_norm * evolved_difference_unit_length_t6) #-0.5470812
angle_radians.G2_SHCSBSB <- acos(dot_product.G2_SHCSBSB)
angle_degrees.G2_SHCSBSB <- angle_radians.G2_SHCSBSB * (180 / pi)
#123.167; 56.833

dot_product.G3_NKLS <- sum(G3_NKLS_norm * evolved_difference_unit_length_t1) #0.06532423
angle_radians.G3_NKLS <- acos(dot_product.G3_NKLS)
angle_degrees.G3_NKLS <- angle_radians.G3_NKLS * (180 / pi)
#86.25453

dot_product.G3_NKBS <- sum(G3_NKBS_norm * evolved_difference_unit_length_t3) #-0.1528877
angle_radians.G3_NKBS <- acos(dot_product.G3_NKBS)
angle_degrees.G3_NKBS <- angle_radians.G3_NKBS * (180 / pi)
#98.79431; 81.20569

dot_product.G3_tewk <- sum(G3_tewk_norm * evolved_difference_unit_length_t3) #-0.2831824
angle_radians.G3_tewk <- acos(dot_product.G3_tewk)
angle_degrees.G3_tewk <- angle_radians.G3_tewk * (180 / pi)
#106.4502; 73.5498

dot_product.G3_uki <- sum(G3_uki_norm * evolved_difference_unit_length_t4) #-0.3356188
angle_radians.G3_uki <- acos(dot_product.G3_uki)
angle_degrees.G3_uki <- angle_radians.G3_uki * (180 / pi)
#109.6102; 70.3898

dot_product.G3_tai <- sum(G3_tai_norm * evolved_difference_unit_length_t5) #0.4712344
angle_radians.G3_tai <- acos(dot_product.G3_tai)
angle_degrees.G3_tai <- angle_radians.G3_tai * (180 / pi)
#61.88555

dot_product.G3_SHCSBSB <- sum(G3_SHCSBSB_norm * evolved_difference_unit_length_t6) #0.7098425
angle_radians.G3_SHCSBSB <- acos(dot_product.G3_SHCSBSB)
angle_degrees.G3_SHCSBSB <- angle_radians.G3_SHCSBSB * (180 / pi)
#44.7779

dot_product.G4_NKLS <- sum(G4_NKLS_norm * evolved_difference_unit_length_t1) #0.1822874
angle_radians.G4_NKLS <- acos(dot_product.G4_NKLS)
angle_degrees.G4_NKLS <- angle_radians.G4_NKLS * (180 / pi)
#79.49698

dot_product.G4_NKBS <- sum(G4_NKBS_norm * evolved_difference_unit_length_t4) #-0.2379446
angle_radians.G4_NKBS <- acos(dot_product.G4_NKBS)
angle_degrees.G4_NKBS <- angle_radians.G4_NKBS * (180 / pi)
#103.7653; 76.2347

dot_product.G4_tewk <- sum(G4_tewk_norm * evolved_difference_unit_length_t3) #0.06053595
angle_radians.G4_tewk <- acos(dot_product.G4_tewk)
angle_degrees.G4_tewk <- angle_radians.G4_tewk * (180 / pi)
#86.52942

dot_product.G4_uki <- sum(G4_uki_norm * evolved_difference_unit_length_t4) #-0.2529781
angle_radians.G4_uki <- acos(dot_product.G4_uki)
angle_degrees.G4_uki <- angle_radians.G4_uki * (180 / pi)
#104.6538; 75.3462

dot_product.G4_tai <- sum(G4_tai_norm * evolved_difference_unit_length_t5) #-0.2860278
angle_radians.G4_tai <- acos(dot_product.G4_tai)
angle_degrees.G4_tai <- angle_radians.G4_tai * (180 / pi)
#106.6203; 73.3797

dot_product.G4_SHCSBSB <- sum(G4_SHCSBSB_norm * evolved_difference_unit_length_t6) #-0.01766609
angle_radians.G4_SHCSBSB <- acos(dot_product.G4_SHCSBSB)
angle_degrees.G4_SHCSBSB <- angle_radians.G4_SHCSBSB * (180 / pi)
#91.01225; 88.98775

#### LOOK AT TRENDS AS A FUNCTION OF TIME -----
form.df <- form.meta[c(1:3, 5:8),] #in same order as mean_by_formation
mean_by_formation

for(i in 1:nrow(form.df)){
    form.df$mean.age[i] <- mean(form.df$Start_age[i], form.df$End_age[i], na.rm = TRUE)
}

form.df$age.range <- ""
for(i in 1:nrow(form.df)){
    form.df$age.range[i] <- form.df$Start_age[i] - form.df$End_age[i]
}
form.df$age.range <- as.numeric(form.df$age.range)

last_row <- c(NA, NA, NA)
diff_between_Gmax_z <- rbind(diff_between_Gmax_z, last_row)
diff_between_Gmax_z$formation <- formation_list
df.diff <- merge(diff_between_Gmax_z, form.df,
                 by.x = "formation",
                 by.y = "formationCode")
df.diff <- merge(df.diff, mean_by_formation,
                 by.x = "formation",
                 by.y = "formation")
df.diff$formation <- factor(df.diff$formation,
                            levels = c("NKLS", "NKBS",
                                       "Tewkesbury", "Upper Kai-Iwi", "Tainui",
                                       "SHCSBSB", "modern"))

ggplot(data = df.diff) +
    geom_point(aes(x = age.range, y = angle_diff_Gmax_to_z),
               shape = 15) + 
    plot.theme +
    scale_x_continuous(name = "Age Range (Ma)") +
    scale_y_continuous(name = "Angle away from Gmax")

df.diff$corr.diff_Gmax_to_z <- as.numeric(df.diff$corr.diff_Gmax_to_z)

p.dot.prod_gmax <- ggplot(data = df.diff) +
    geom_point(aes(x = formation, y = corr.diff_Gmax_to_z),
               size = 5, shape = 15) + 
    plot.theme +
    #scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Dot Product of Phenotypic change from Gmax")

ggsave(p.dot.prod_gmax, 
       file = "./Results/corr.gmax.p.diff.png", 
       width = 14, height = 10, units = "cm")

p.ang_g <- ggplot(diff_between_Gs) +
    geom_point(aes(x = angle.diff_Gs.time, y = angle_diff_Gs),
               size = 5, shape = 15) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between G matrices", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.ang_g, 
       file = "./Results/angle.g.diff.png", 
       width = 20, height = 20, units = "cm")

diff_between_Gs$corr.diff_Gs <- as.numeric(diff_between_Gs$corr.diff_Gs)
p.dot.prod_g <- ggplot(diff_between_Gs) +
    geom_point(aes(x = angle.diff_Gs.time, y = abs(corr.diff_Gs)),
               size = 5, shape = 15) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Dot prodcut between G matrices") + 
    plot.theme

ggsave(p.dot.prod_g, 
       file = "./Results/dot.prod.g.diff.png", 
       width = 20, height = 20, units = "cm")

#### CALCULATE BETA ----
#∆z = ß*G
#ß = ∆z/G
#∆z = vector of change
#G is matrix at t1

beta_t1 = evolved_difference_unit_length_t1%*%solve(G_ext_NKLS)
beta_t2 = evolved_difference_unit_length_t2%*%solve(G_ext_NKBS)
beta_t3 = evolved_difference_unit_length_t3%*%solve(G_ext_tewk)
beta_t4 = evolved_difference_unit_length_t4%*%solve(G_ext_uki)
beta_t5 = evolved_difference_unit_length_t5%*%solve(G_ext_tai)
beta_t6 = evolved_difference_unit_length_t6%*%solve(G_ext_SHCSBSB)

beta_t1_norm <- f.normalize_vector(beta_t1)
beta_t2_norm <- f.normalize_vector(beta_t2)
beta_t3_norm <- f.normalize_vector(beta_t3)
beta_t4_norm <- f.normalize_vector(beta_t4)
beta_t5_norm <- f.normalize_vector(beta_t5)
beta_t6_norm <- f.normalize_vector(beta_t6)

##### DOT PRODUCT BETA T1 TO GMAX T1 -----
# Calculate the dot product of the unit vectors; tells number 0 to 1
dot_product.beta_t1_Gmax_t1 <- sum(beta_t1_norm * Gmax_NKLS_norm) #0.04595859
dot_product.beta_t2_Gmax_t2 <- sum(beta_t2_norm * Gmax_NKBS_norm) #0.07747256
dot_product.beta_t3_Gmax_t3 <- sum(beta_t3_norm * Gmax_tewk_norm) #-0.2514228
dot_product.beta_t4_Gmax_t4 <- sum(beta_t4_norm * Gmax_uki_norm) #-0.1437962
dot_product.beta_t5_Gmax_t5 <- sum(beta_t5_norm * Gmax_tai_norm) #0.1747052
dot_product.beta_t6_Gmax_t6 <- sum(beta_t6_norm * Gmax_SHCSBSB_norm) #0.00134973

#convert to angle
angle_radians.beta_t1_Gmax_t1 <- acos(dot_product.beta_t1_Gmax_t1)
# Convert the angle to degrees
angle_degrees.beta_t1_Gmax_t1 <- angle_radians.beta_t1_Gmax_t1 * (180 / pi) #87.36584

angle_radians.beta_t2_Gmax_t2 <- acos(dot_product.beta_t2_Gmax_t2)
angle_degrees.beta_t2_Gmax_t2 <- angle_radians.beta_t2_Gmax_t2 * (180 / pi) #85.5567

angle_radians.beta_t3_Gmax_t3 <- acos(dot_product.beta_t3_Gmax_t3)
angle_degrees.beta_t3_Gmax_t3 <- angle_radians.beta_t3_Gmax_t3 * (180 / pi) #104.5617; 75.4383

angle_radians.beta_t4_Gmax_t4 <- acos(dot_product.beta_t4_Gmax_t4)
angle_degrees.beta_t4_Gmax_t4 <- angle_radians.beta_t4_Gmax_t4 * (180 / pi) #98.26757; 81.73243

angle_radians.beta_t5_Gmax_t5 <- acos(dot_product.beta_t5_Gmax_t5)
angle_degrees.beta_t5_Gmax_t5 <- angle_radians.beta_t5_Gmax_t5 * (180 / pi) #79.9385

angle_radians.beta_t6_Gmax_t6 <- acos(dot_product.beta_t6_Gmax_t6)
angle_degrees.beta_t6_Gmax_t6 <- angle_radians.beta_t6_Gmax_t6 * (180 / pi) #89.92267

##### DOT PRODUCT BETA T1 TO BETA T2 -----
# Calculate the dot product of the unit vectors; tells number 0 to 1
dot_product.beta_t1_beta_t2 <- sum(beta_t1_norm * beta_t2_norm) #0.2412896
dot_product.beta_t2_beta_t3 <- sum(beta_t2_norm * beta_t3_norm) #-0.4428297
dot_product.beta_t3_beta_t4 <- sum(beta_t3_norm * beta_t4_norm) #-0.6314299
dot_product.beta_t4_beta_t5 <- sum(beta_t4_norm * beta_t5_norm) #-0.9118744
dot_product.beta_t5_beta_t6 <- sum(beta_t5_norm * beta_t6_norm) #-0.8769132

#convert to angle
angle_radians.beta_t1_beta_t2 <- acos(dot_product.beta_t1_beta_t2)
# Convert the angle to degrees
angle_degrees.beta_t1_beta_t2 <- angle_radians.beta_t1_beta_t2 * (180 / pi) #76.03734

angle_radians.beta_t2_beta_t3 <- acos(dot_product.beta_t2_beta_t3)
angle_degrees.beta_t2_beta_t3 <- angle_radians.beta_t2_beta_t3 * (180 / pi) #116.2846; 63.7154

angle_radians.beta_t3_beta_t4 <- acos(dot_product.beta_t3_beta_t4)
angle_degrees.beta_t3_beta_t4 <- angle_radians.beta_t3_beta_t4 * (180 / pi) #129.1557; 50.8443

angle_radians.beta_t4_beta_t5 <- acos(dot_product.beta_t4_beta_t5)
angle_degrees.beta_t4_beta_t5 <- angle_radians.beta_t4_beta_t5 * (180 / pi) #155.7657; 24.2343

angle_radians.beta_t5_beta_t6 <- acos(dot_product.beta_t5_beta_t6)
angle_degrees.beta_t5_beta_t6 <- angle_radians.beta_t5_beta_t6 * (180 / pi) #151.2722; 28.7278

##### DOT PRODUCT BETA T1 TO GMAX T2 -----
# Calculate the dot product of the unit vectors; tells number 0 to 1
dot_product.beta_t1_Gmax_t2 <- sum(beta_t1_norm * Gmax_NKBS_norm) #-0.06278445
dot_product.beta_t2_Gmax_t3 <- sum(beta_t2_norm * Gmax_tewk_norm) #0.14208
dot_product.beta_t3_Gmax_t4 <- sum(beta_t3_norm * Gmax_uki_norm) #-0.1573955
dot_product.beta_t4_Gmax_t5 <- sum(beta_t4_norm * Gmax_tai_norm) #-0.1613059
dot_product.beta_t5_Gmax_t6 <- sum(beta_t5_norm * Gmax_SHCSBSB_norm) #0.09014606
dot_product.beta_t6_Gmax_t7 <- sum(beta_t6_norm * Gmax_mod_norm) #-0.06645504

#convert to angle
angle_radians.beta_t1_Gmax_t2 <- acos(dot_product.beta_t1_Gmax_t2)
# Convert the angle to degrees
angle_degrees.beta_t1_Gmax_t2 <- angle_radians.beta_t1_Gmax_t2 * (180 / pi) #93.59965; 86.40035

angle_radians.beta_t2_Gmax_t3 <- acos(dot_product.beta_t2_Gmax_t3)
# Convert the angle to degrees
angle_degrees.beta_t2_Gmax_t3 <- angle_radians.beta_t2_Gmax_t3 * (180 / pi) #81.83178

angle_radians.beta_t3_Gmax_t4 <- acos(dot_product.beta_t3_Gmax_t4)
# Convert the angle to degrees
angle_degrees.beta_t3_Gmax_t4 <- angle_radians.beta_t3_Gmax_t4 * (180 / pi) #99.05575; 80.94425

angle_radians.beta_t4_Gmax_t5 <- acos(dot_product.beta_t4_Gmax_t5)
# Convert the angle to degrees
angle_degrees.beta_t4_Gmax_t5 <- angle_radians.beta_t4_Gmax_t5 * (180 / pi) #99.2827; 80.7173

angle_radians.beta_t5_Gmax_t6 <- acos(dot_product.beta_t5_Gmax_t6)
# Convert the angle to degrees
angle_degrees.beta_t5_Gmax_t6 <- angle_radians.beta_t5_Gmax_t6 * (180 / pi) #84.82799

angle_radians.beta_t6_Gmax_t7 <- acos(dot_product.beta_t6_Gmax_t7)
# Convert the angle to degrees
angle_degrees.beta_t6_Gmax_t7 <- angle_radians.beta_t6_Gmax_t7 * (180 / pi) #93.8104; 86.1896

##### MAGNITUDE BETA TO DIFF IN DOT PROD -----
## calculate magnitude ß
#use magnitude function

mag.beta_t1 <- magnitude(beta_t1)
mag.beta_t2 <- magnitude(beta_t2)
mag.beta_t3 <- magnitude(beta_t3)
mag.beta_t4 <- magnitude(beta_t4)
mag.beta_t5 <- magnitude(beta_t5)
mag.beta_t6 <- magnitude(beta_t6)

mag.beta <- c(mag.beta_t1, mag.beta_t2,
              mag.beta_t3, mag.beta_t4,
              mag.beta_t5, mag.beta_t6, "")

## calculate difference in dot product
# dot prod Gmax t1 and Gmax t2 - dot product ß t1 and Gmax t2
diff.gmax.b.nkls.nkbs <- abs(dot_product.Gmax_NKLS_NKBS) - abs(dot_product.beta_t1_Gmax_t2)
diff.gmax.b.nkbs.tewk <- abs(dot_product.Gmax_NKBS_tewk) - abs(dot_product.beta_t2_Gmax_t3)
diff.gmax.b.tewk.uki <- abs(dot_product.Gmax_tewk_uki) - abs(dot_product.beta_t3_Gmax_t4)
diff.gmax.b.uki.tai <- abs(dot_product.Gmax_uki_tai) - abs(dot_product.beta_t4_Gmax_t5)
diff.gmax.b.tai.shcsbsb <- abs(dot_product.Gmax_tai_SHCSBSB) - abs(dot_product.beta_t5_Gmax_t6)
diff.gmax.b.shcsbsb.mod <- abs(dot_product.Gmax_SHCSBSB_mod) - abs(dot_product.beta_t6_Gmax_t7)

diff.gmax.b <- c(diff.gmax.b.nkls.nkbs, diff.gmax.b.nkbs.tewk,
                 diff.gmax.b.tewk.uki, diff.gmax.b.uki.tai,
                 diff.gmax.b.tai.shcsbsb, diff.gmax.b.shcsbsb.mod, "")

#angle
ang.diff.gmax.b.nkls.nkbs <- abs(angle_degrees.Gmax_NKLS_NKBS) - abs(angle_radians.beta_t1_Gmax_t2)
ang.diff.gmax.b.nkbs.tewk <- abs(angle_degrees.Gmax_NKBS_tewk) - abs(angle_radians.beta_t2_Gmax_t3)
ang.diff.gmax.b.tewk.uki <- abs(angle_degrees.Gmax_tewk_uki) - abs(angle_radians.beta_t3_Gmax_t4)
ang.diff.gmax.b.uki.tai <- abs(angle_degrees.Gmax_uki_tai) - abs(angle_radians.beta_t4_Gmax_t5)
ang.diff.gmax.b.tai.shcsbsb <- abs(angle_degrees.Gmax_tai_SHCSBSB) - abs(angle_radians.beta_t5_Gmax_t6)
ang.diff.gmax.b.shcsbsb.mod <- abs(angle_degrees.Gmax_SHCSBSB_mod) - abs(angle_radians.beta_t6_Gmax_t7)

ang.diff.gmax.b <- c(ang.diff.gmax.b.nkls.nkbs, ang.diff.gmax.b.nkbs.tewk,
                     ang.diff.gmax.b.tewk.uki, ang.diff.gmax.b.uki.tai,
                     ang.diff.gmax.b.tai.shcsbsb, ang.diff.gmax.b.shcsbsb.mod, "")

## make df

diff.beta.gmax.df <- as.data.frame(cbind(mag.beta, diff.gmax.b, formation_transition, ang.diff.gmax.b))
diff.beta.gmax.df$formation_transition <- factor(diff.beta.gmax.df$formation_transition,
                                                 levels = c("NKLS to NKBS", 
                                                            "NKBS to Tewkesbury",
                                                            "Tewkesbury to Upper Kai-Iwi",
                                                            "Upper Kai-Iwi to Tainui", 
                                                            "Tainui to SHCSBSB",
                                                            "SHCSBSB to modern", ""))


## plot
#expectations: positive, negative, or no relationship
#no relationship means no effect of ß on Gmax

p.diff.b.gmax <- ggplot(diff.beta.gmax.df,
       aes(x = as.numeric(mag.beta), y = as.numeric(diff.gmax.b))) + 
    geom_point() +
    geom_smooth(method = "lm") +
    plot.theme +
    scale_x_continuous(expression(Magnitude~beta)) +
    scale_y_continuous(expression(reduction~of~correlation~of~G[max]~and~beta))

ggsave(p.diff.b.gmax, 
       file = "./Results/reduction.gmax.beta.to.mag.beta.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(diff.gmax.b) ~ as.numeric(mag.beta),
           data = diff.beta.gmax.df)) #barely sig at 0.04184; slope at 0, r2 at 0.61

p.ang.diff.b.gmax <- ggplot(diff.beta.gmax.df,
                        aes(x = as.numeric(mag.beta), y = as.numeric(ang.diff.gmax.b))) + 
    geom_point() +
    geom_smooth(method = "lm") +
    plot.theme +
    scale_x_continuous(expression(Magnitude~beta)) +
    scale_y_continuous(expression(reduction~of~angle~between~G[max]~and~beta))

ggsave(p.ang.diff.b.gmax, 
       file = "./Results/reduction.ang.gmax.beta.to.mag.beta.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(ang.diff.gmax.b) ~ as.numeric(mag.beta),
           data = diff.beta.gmax.df)) #nonsig at 0.1064; slope barely over 0, r2 0.3985

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
    geom_smooth(method = "lm") +
    plot.theme +
    scale_x_continuous(expression(Magnitude~beta)) +
    scale_y_continuous(expression(distance~between~G~matrices~and~beta))

ggsave(p.dist.gmat.b, 
       file = "./Results/dist.g.mat.to.mag.beta.png", 
       width = 14, height = 10, units = "cm")

summary(lm(as.numeric(dist.gmat) ~ as.numeric(mag.beta),
           data = diff.beta.gmax.df)) #nonsig at p = 0.4008; no relationship

#### SUBSTITUTE P FOR G ----
P_ext_NKLS = round(as.matrix(P_ext[[1]]), 6) # The G matrix estimated for sample/formation 1
P_ext_NKBS = round(as.matrix(P_ext[[2]]), 6) # The G matrix estimated for sample/formation 2
P_ext_tewk = round(as.matrix(P_ext[[3]]), 6) # The G matrix estimated for sample/formation 3
P_ext_uki = round(as.matrix(P_ext[[4]]), 6) # The G matrix estimated for sample/formation 5
P_ext_tai = round(as.matrix(P_ext[[5]]), 6) # The G matrix estimated for sample/formation 6
P_ext_SHCSBSB = round(as.matrix(P_ext[[6]]), 6) # The G matrix estimated for sample/formation 7
P_ext_mod = round(as.matrix(P_ext[[7]]), 6) # The G matrix estimated for sample/formation 7

is.symmetric.matrix(P_ext_NKLS)
is.positive.definite(P_ext_NKLS)

is.symmetric.matrix(P_ext_NKBS)
is.positive.definite(P_ext_NKBS)

is.symmetric.matrix(P_ext_tewk)
is.positive.definite(P_ext_tewk)

is.symmetric.matrix(P_ext_uki)
is.positive.definite(P_ext_uki)

is.symmetric.matrix(P_ext_tai)
is.positive.definite(P_ext_tai)

is.symmetric.matrix(P_ext_SHCSBSB)
is.positive.definite(P_ext_SHCSBSB)

is.symmetric.matrix(P_ext_mod)
is.positive.definite(P_ext_mod)

#trait means by time
mean_by_formation
#order of formations:
# "NKLS"  "NKBS" "Tewkesbury" "Upper Kai-Iwi" "Tainui" "SHCSBSB"

### Calculate the vector that defines the observed divergence between sample/formation 1 an 2
#use from section above
#NKLS <- as.numeric(mean_by_formation[1, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 1 
#NKBS <- as.numeric(mean_by_formation[2, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 2 
#tewk <- as.numeric(mean_by_formation[3, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 3
#uki <- as.numeric(mean_by_formation[5, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 5 
#tai <- as.numeric(mean_by_formation[6, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 6 
#SHCSBSB <- as.numeric(mean_by_formation[7, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 7
#mod <- as.numeric(mean_by_formation[8, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 7

#second - first
#t1 = NKBS - NKLS
#t2 = tewk - NKBS
#t3 = uki - tewk
#t4 = tai - uki
#t5 = SHCSBSB - tai
#t6 = mod - SHCSBSB

## really need to learn how to name things in functions...
#use from section above
#evolved_difference_unit_length_t1 <- f.normalize_vector(NKBS - NKLS)
#evolved_difference_unit_length_t2 <- f.normalize_vector(tewk - NKBS)
#evolved_difference_unit_length_t3 <- f.normalize_vector(uki - tewk)
#evolved_difference_unit_length_t4 <- f.normalize_vector(tai - uki)
#evolved_difference_unit_length_t5 <- f.normalize_vector(SHCSBSB - tai)
#evolved_difference_unit_length_t6 <- f.normalize_vector(mod - SHCSBSB)

###### OBSERVED EVOLVABILITY ------
### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
p_observed_evolvability_in_direction_of_change_t1 <- t(evolved_difference_unit_length_t1)%*%as.matrix(P_ext_NKLS)%*%evolved_difference_unit_length_t1
p_observed_evolvability_in_direction_of_change_t2 <- t(evolved_difference_unit_length_t2)%*%as.matrix(P_ext_NKBS)%*%evolved_difference_unit_length_t2
p_observed_evolvability_in_direction_of_change_t3 <- t(evolved_difference_unit_length_t3)%*%as.matrix(P_ext_tewk)%*%evolved_difference_unit_length_t3
p_observed_evolvability_in_direction_of_change_t4 <- t(evolved_difference_unit_length_t4)%*%as.matrix(P_ext_uki)%*%evolved_difference_unit_length_t4
p_observed_evolvability_in_direction_of_change_t5 <- t(evolved_difference_unit_length_t5)%*%as.matrix(P_ext_tai)%*%evolved_difference_unit_length_t5
p_observed_evolvability_in_direction_of_change_t6 <- t(evolved_difference_unit_length_t6)%*%as.matrix(P_ext_SHCSBSB)%*%evolved_difference_unit_length_t6

###### OBSERVED CONDITIONAL EVOLVABILITY ------
### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
p_observed_conditional_evolvability_in_direction_of_change_t1 <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(P_ext_NKLS))%*%evolved_difference_unit_length_t1)
p_observed_conditional_evolvability_in_direction_of_change_t2 <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(P_ext_NKBS))%*%evolved_difference_unit_length_t2)
p_observed_conditional_evolvability_in_direction_of_change_t3 <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(P_ext_tewk))%*%evolved_difference_unit_length_t3)
p_observed_conditional_evolvability_in_direction_of_change_t4 <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(P_ext_uki))%*%evolved_difference_unit_length_t4)
p_observed_conditional_evolvability_in_direction_of_change_t5 <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(P_ext_tai))%*%evolved_difference_unit_length_t5)
p_observed_conditional_evolvability_in_direction_of_change_t6 <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(P_ext_SHCSBSB))%*%evolved_difference_unit_length_t6)

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
p_X_t1 <- evolvabilityBeta(as.matrix(P_ext_NKLS), Beta)
p_sumX_t1 <- summary(p_X_t1) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_t2 <- evolvabilityBeta(as.matrix(P_ext_NKBS), Beta)
p_sumX_t2 <- summary(p_X_t2) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_t3 <- evolvabilityBeta(as.matrix(P_ext_tewk), Beta)
p_sumX_t3 <- summary(p_X_t3) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_t4 <- evolvabilityBeta(as.matrix(P_ext_uki), Beta)
p_sumX_t4 <- summary(p_X_t4) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_t5 <- evolvabilityBeta(as.matrix(P_ext_tai), Beta)
p_sumX_t5 <- summary(p_X_t5) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_t6 <- evolvabilityBeta(as.matrix(P_ext_SHCSBSB), Beta)
p_sumX_t6 <- summary(p_X_t6) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_t7 <- evolvabilityBeta(as.matrix(P_ext_mod), Beta)
p_sumX_t7 <- summary(p_X_t7) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

p_X_sum <- data.frame(p_c.mean = c(p_sumX_t1$Averages[[3]], p_sumX_t2$Averages[[3]], p_sumX_t3$Averages[[3]],
                                   p_sumX_t4$Averages[[3]], p_sumX_t5$Averages[[3]], p_sumX_t6$Averages[[3]], 
                                   p_sumX_t7$Averages[[3]]),
                      p_c.min = c(p_sumX_t1$Minimum[[3]], p_sumX_t2$Minimum[[3]], p_sumX_t3$Minimum[[3]],
                                  p_sumX_t4$Minimum[[3]], p_sumX_t5$Minimum[[3]], p_sumX_t6$Minimum[[3]],
                                  p_sumX_t7$Minimum[[3]]),
                      p_c.max = c(p_sumX_t1$Maximum[[3]], p_sumX_t2$Maximum[[3]], p_sumX_t3$Maximum[[3]],
                                  p_sumX_t4$Maximum[[3]], p_sumX_t5$Maximum[[3]], p_sumX_t6$Maximum[[3]],
                                  p_sumX_t7$Maximum[[3]]),
                      p_e.mean = c(p_sumX_t1$Averages[[1]], p_sumX_t2$Averages[[1]], p_sumX_t3$Averages[[1]],
                                   p_sumX_t4$Averages[[1]], p_sumX_t5$Averages[[1]], p_sumX_t6$Averages[[1]],
                                   p_sumX_t7$Averages[[1]]),
                      p_e.min = c(p_sumX_t1$Minimum[[1]], p_sumX_t2$Minimum[[1]], p_sumX_t3$Minimum[[1]],
                                  p_sumX_t4$Minimum[[1]], p_sumX_t5$Minimum[[1]], p_sumX_t6$Minimum[[1]],
                                  p_sumX_t7$Minimum[[1]]),
                      p_e.max = c(p_sumX_t1$Maximum[[1]], p_sumX_t2$Maximum[[1]], p_sumX_t3$Maximum[[1]],
                                  p_sumX_t4$Maximum[[1]], p_sumX_t5$Maximum[[1]], p_sumX_t6$Maximum[[1]],
                                  p_sumX_t7$Maximum[[1]]),
                      p_observed_e = c(p_observed_evolvability_in_direction_of_change_t1,
                                       p_observed_evolvability_in_direction_of_change_t2,
                                       p_observed_evolvability_in_direction_of_change_t3,
                                       p_observed_evolvability_in_direction_of_change_t4,
                                       p_observed_evolvability_in_direction_of_change_t5,
                                       p_observed_evolvability_in_direction_of_change_t6,
                                       ""),
                      p_observed_c = c(p_observed_conditional_evolvability_in_direction_of_change_t1,
                                       p_observed_conditional_evolvability_in_direction_of_change_t2,
                                       p_observed_conditional_evolvability_in_direction_of_change_t3,
                                       p_observed_conditional_evolvability_in_direction_of_change_t4,
                                       p_observed_conditional_evolvability_in_direction_of_change_t5,
                                       p_observed_conditional_evolvability_in_direction_of_change_t6,
                                       ""),
                      row.names = levels(formation_list))
#NO NEGATIVE VALUES!

write.csv(p_X_sum,
          "./Results/p_evolvability.summary.csv")

## PLOT
p_X_sum$formation <- rownames(p_X_sum)
p_X_sum$formation <- factor(p_X_sum$formation, levels = c("NKLS", "NKBS",
                                                          "Tewkesbury", "Upper Kai-Iwi", 
                                                          "Tainui", "SHCSBSB", "modern"))

p_X_sum$form.trans <- formation_transition
p_X_sum$form.trans <- factor(p_X_sum$form.trans,
                             levels = c("NKLS to NKBS", 
                                        "NKBS to Tewkesbury",
                                        "Tewkesbury to Upper Kai-Iwi",
                                        "Upper Kai-Iwi to Tainui", 
                                        "Tainui to SHCSBSB",
                                        "SHCSBSB to modern", ""))


p_X_sum.trim <- p_X_sum[1:6,]
p.p_evol <- ggplot(p_X_sum.trim, aes(x = form.trans)) +
    geom_boxplot(aes(ymin = p_e.min, 
                     lower = p_e.min,
                     middle = p_e.mean,
                     ymax = p_e.max,
                     upper = p_e.max,
                     fill = "gray"),
                 stat = "identity", fill = "gray") +
    geom_point(aes(y = as.numeric(p_observed_e),
                   color = "black"),
               size = 5, shape = 17, color = "black") +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Evolvability") +
    plot.theme

ggsave(p.p_evol, 
       file = "./Results/p_evolvability.png", 
       width = 14, height = 10, units = "cm")

# By comparing the evolvabilities you estimated in the direction of change (lines 9 and 12) with the average evolvabilities calculated by running line 20, you get a sense of whether evolution happened in directions with above or below average evolvability.  

##### ANGLE CHANGE IN P BETWEEN FORMATIONS -----

### Proportion of variance in n-dimensional trait space that is explained by PC1 (i.e., the first eigenvector)
#eigen(as.matrix(G_matrix_1))$values[1]/sum(eigen(as.matrix(G_matrix_1))$values)
eigen(as.matrix(P_ext_NKLS))$values[1]/sum(eigen(as.matrix(P_ext_NKLS))$values) #0.4448337; G: 0.4731959
eigen(as.matrix(P_ext_NKBS))$values[1]/sum(eigen(as.matrix(P_ext_NKBS))$values) #0.4611266; G: 0.4819061
eigen(as.matrix(P_ext_tewk))$values[1]/sum(eigen(as.matrix(P_ext_tewk))$values) #0.4759648; G: 0.4578517
eigen(as.matrix(P_ext_uki))$values[1]/sum(eigen(as.matrix(P_ext_uki))$values) #0.4968204; G: 0.5058219
eigen(as.matrix(P_ext_tai))$values[1]/sum(eigen(as.matrix(P_ext_tai))$values) #0.4023307; G: 0.3844157
eigen(as.matrix(P_ext_SHCSBSB))$values[1]/sum(eigen(as.matrix(P_ext_SHCSBSB))$values) #0.4271263; G: 0.4479779
eigen(as.matrix(P_ext_mod))$values[1]/sum(eigen(as.matrix(P_ext_mod))$values) #0.4686239; G: 0.5337849

### How much is the direction of Pmax (i.e., the direction first ) varying between different G-matrices? 
Pmax_NKLS <- eigen(P_ext_NKLS)$vectors[,1]
Pmax_NKBS <- eigen(P_ext_NKBS)$vectors[,1]
Pmax_tewk <- eigen(P_ext_tewk)$vectors[,1]
Pmax_uki <- eigen(P_ext_uki)$vectors[,1]
Pmax_tai <- eigen(P_ext_tai)$vectors[,1]
Pmax_SHCSBSB <- eigen(P_ext_SHCSBSB)$vectors[,1]
Pmax_mod <- eigen(P_ext_mod)$vectors[,1]

# Put Pmax to norm length
Pmax_NKLS_norm <- f.normalize_vector(Pmax_NKLS)
Pmax_NKBS_norm <- f.normalize_vector(Pmax_NKBS)
Pmax_tewk_norm <- f.normalize_vector(Pmax_tewk)
Pmax_uki_norm <- f.normalize_vector(Pmax_uki)
Pmax_tai_norm <- f.normalize_vector(Pmax_tai)
Pmax_SHCSBSB_norm <- f.normalize_vector(Pmax_SHCSBSB)
Pmax_mod_norm <- f.normalize_vector(Pmax_mod)

# Calculate the dot product of the unit vectors
dot_product.Pmax_NKLS_NKBS <- sum(Pmax_NKLS_norm * Pmax_NKBS_norm) #0.9934551
# Calculate the angle in radians
angle_radians.Pmax_NKLS_NKBS <- acos(dot_product.Pmax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Pmax_NKLS_NKBS <- angle_radians.Pmax_NKLS_NKBS * (180 / pi)
#6.558851

dot_product.Pmax_NKBS_tewk <- sum(Pmax_NKBS_norm * Pmax_tewk_norm) #0.9970584
# Calculate the angle in radians
angle_radians.Pmax_NKBS_tewk <- acos(dot_product.Pmax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Pmax_NKBS_tewk <- angle_radians.Pmax_NKBS_tewk * (180 / pi)
#4.395792

dot_product.Pmax_tewk_uki <- sum(Pmax_tewk_norm * Pmax_uki_norm) #0.9622356
# Calculate the angle in radians
angle_radians.Pmax_tewk_uki <- acos(dot_product.Pmax_tewk_uki)
# Convert the angle to degrees
angle_degrees.Pmax_tewk_uki <- angle_radians.Pmax_tewk_uki * (180 / pi)
#15.79629

dot_product.Pmax_uki_tai <- sum(Pmax_uki_norm * Pmax_tai_norm) #0.9729122
# Calculate the angle in radians
angle_radians.Pmax_uki_tai <- acos(dot_product.Pmax_uki_tai)
# Convert the angle to degrees
angle_degrees.Pmax_uki_tai <- angle_radians.Pmax_uki_tai * (180 / pi)
#13.36625

dot_product.Pmax_tai_SHCSBSB <- sum(Pmax_tai_norm * Pmax_SHCSBSB_norm) #0.9848919
# Calculate the angle in radians
angle_radians.Pmax_tai_SHCSBSB <- acos(dot_product.Pmax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Pmax_tai_SHCSBSB <- angle_radians.Pmax_tai_SHCSBSB * (180 / pi)
#9.972182

dot_product.Pmax_SHCSBSB_mod <- sum(Pmax_SHCSBSB_norm * Pmax_mod_norm) #0.8828773
# Calculate the angle in radians
angle_radians.Pmax_SHCSBSB_mod <- acos(dot_product.Pmax_SHCSBSB_mod)
# Convert the angle to degrees
angle_degrees.Pmax_SHCSBSB_mod <- angle_radians.Pmax_SHCSBSB_mod * (180 / pi)
#28.00858

corr.diff_Ps <- c(dot_product.Pmax_NKLS_NKBS, dot_product.Pmax_NKBS_tewk,
                  dot_product.Pmax_tewk_uki, dot_product.Pmax_uki_tai,
                  dot_product.Pmax_tai_SHCSBSB, dot_product.Pmax_SHCSBSB_mod, "")

angle_diff_Ps <- c(angle_degrees.Pmax_NKLS_NKBS, angle_degrees.Pmax_NKBS_tewk,
                   angle_degrees.Pmax_tewk_uki,
                   angle_degrees.Pmax_uki_tai, angle_degrees.Pmax_tai_SHCSBSB,
                   angle_degrees.Pmax_SHCSBSB_mod, "")

diff_between_Ps <- as.data.frame(cbind(corr.diff_Ps, angle_diff_Ps))
diff_between_Ps$angle_diff_Ps <- as.numeric(diff_between_Ps$angle_diff_Ps)
diff_between_Ps$angle.diff_Ps.time <- formation_transition
diff_between_Ps$angle.diff_Ps.time <- factor(diff_between_Ps$angle.diff_Ps.time,
                                             levels = c("NKLS to NKBS", 
                                                        "NKBS to Tewkesbury",
                                                        "Tewkesbury to Upper Kai-Iwi",
                                                        "Upper Kai-Iwi to Tainui", 
                                                        "Tainui to SHCSBSB",
                                                        "SHCSBSB to modern", ""))

for(i in 1:nrow(diff_between_Ps)){
    if(isTRUE(diff_between_Ps$angle_diff_Ps[i] > 90)){
        diff_between_Ps$angle_diff_Ps[i] <- 180 - diff_between_Ps$angle_diff_Ps[i]
    }
    else{
        next
    }
}

diff_between_Ps$angle.diff_Ps.time <- factor(diff_between_Ps$angle.diff_Ps.time,
                                                   levels = c("NKLS to NKBS", 
                                                              "NKBS to Tewkesbury",
                                                              "Tewkesbury to Upper Kai-Iwi",
                                                              "Upper Kai-Iwi to Tainui", 
                                                              "Tainui to SHCSBSB",
                                                              "SHCSBSB to modern", ""))

write.csv(diff_between_Ps,
          "./Results/differences.between.Ps.csv",
          row.names = FALSE)

diff_between_Ps$corr.diff_Ps <- as.numeric(diff_between_Ps$corr.diff_Ps )
p.dot.prod_p <- ggplot(diff_between_Ps) +
    geom_point(aes(x = angle.diff_Ps.time, y = abs(corr.diff_Ps)),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Dot prodcut between P matrices") + 
    plot.theme

ggsave(p.dot.prod_p, 
       file = "./Results/dot.prod.p.diff.png", 
       width = 20, height = 20, units = "cm")

p.ang_p <- ggplot(diff_between_Ps) +
    geom_point(aes(x = angle.diff_Ps.time, y = angle_diff_Ps),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle between P matrices") + 
    plot.theme

ggsave(p.ang_p, 
       file = "./Results/ang.p.diff.png", 
       width = 20, height = 20, units = "cm")

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO PMAX -----

### See if change is in direction of P max
## use Pmax of t1 and compare to ∆z
# Calculate the dot product of the unit vectors
dot_product.Pmax_NKLS <- sum(Pmax_NKLS_norm * evolved_difference_unit_length_t1) #-0.3122769
# Calculate the angle in radians
angle_radians.Pmax_NKLS <- acos(dot_product.Pmax_NKLS)
# Convert the angle to degrees
angle_degrees.Pmax_NKLS <- angle_radians.Pmax_NKLS * (180 / pi)
#108.1965; 71.8035

# Calculate the dot product of the unit vectors
dot_product.Pmax_NKBS <- sum(Pmax_NKBS_norm * evolved_difference_unit_length_t2) #0.7131075
# Calculate the angle in radians
angle_radians.Pmax_NKBS <- acos(dot_product.Pmax_NKBS)
# Convert the angle to degrees
angle_degrees.Pmax_NKBS <- angle_radians.Pmax_NKBS * (180 / pi)
#44.51169

# Calculate the dot product of the unit vectors
dot_product.Pmax_tewk <- sum(Pmax_tewk_norm * evolved_difference_unit_length_t3) #-0.8358853
# Calculate the angle in radians
angle_radians.Pmax_tewk <- acos(dot_product.Pmax_tewk)
# Convert the angle to degrees
angle_degrees.Pmax_tewk <- angle_radians.Pmax_tewk * (180 / pi)
#146.7081; 33.2919

# Calculate the dot product of the unit vectors
dot_product.Pmax_uki <- sum(Pmax_uki_norm * evolved_difference_unit_length_t4) #-0.6393634
# Calculate the angle in radians
angle_radians.Pmax_uki <- acos(dot_product.Pmax_uki)
# Convert the angle to degrees
angle_degrees.Pmax_uki <- angle_radians.Pmax_uki * (180 / pi)
#129.7444; 50.2556

# Calculate the dot product of the unit vectors
dot_product.Pmax_tai <- sum(Pmax_tai_norm * evolved_difference_unit_length_t5) #0.5922247
# Calculate the angle in radians
angle_radians.Pmax_tai <- acos(dot_product.Pmax_tai)
# Convert the angle to degrees
angle_degrees.Pmax_tai <- angle_radians.Pmax_tai * (180 / pi)
#53.68496

# Calculate the dot product of the unit vectors
dot_product.Pmax_SHCSBSB <- sum(Pmax_SHCSBSB_norm * evolved_difference_unit_length_t6) #-0.2165188
# Calculate the angle in radians
angle_radians.Pmax_SHCSBSB <- acos(dot_product.Pmax_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Pmax_SHCSBSB <- angle_radians.Pmax_SHCSBSB * (180 / pi)
#102.5046; 77.4954

corr.diff_Pmax_to_z <- c(dot_product.Pmax_NKLS, dot_product.Pmax_NKBS,
                         dot_product.Pmax_tewk, dot_product.Pmax_uki,
                         dot_product.Pmax_tai, dot_product.Pmax_SHCSBSB, "")

angle_diff_Pmax_to_z <- c(angle_degrees.Pmax_NKLS, angle_degrees.Pmax_NKBS,
                          angle_degrees.Pmax_tewk, angle_degrees.Pmax_uki, 
                          angle_degrees.Pmax_tai, angle_degrees.Pmax_SHCSBSB, "")
diff_between_Pmax_z <- as.data.frame(cbind(levels(formation_list), angle_diff_Pmax_to_z, corr.diff_Pmax_to_z))
colnames(diff_between_Pmax_z) <- c("formation", "angle_diff_Pmax_to_z", "corr.diff_Pmax_to_z")
diff_between_Pmax_z$angle_diff_Pmax_to_z <- as.numeric(diff_between_Pmax_z$angle_diff_Pmax_to_z)
diff_between_Pmax_z$corr.diff_Pmax_to_z <- as.numeric(diff_between_Pmax_z$corr.diff_Pmax_to_z)

for(i in 1:nrow(diff_between_Pmax_z)){
    if(isTRUE(diff_between_Pmax_z$angle_diff_Pmax_to_z[i] > 90)){
        diff_between_Pmax_z$angle_diff_Pmax_to_z[i] <- 180 - as.numeric(diff_between_Pmax_z$angle_diff_Pmax_to_z[i])
    }
    else{
        next
    }
}

write.csv(diff_between_Pmax_z,
          "./Results/differences.between.Pmax.z.csv",
          row.names = FALSE)

#### LOOK AT TRENDS AS A FUNCTION OF TIME -----

df.diff.p <- merge(diff_between_Pmax_z, form.df,
                 by.x = "formation",
                 by.y = "formationCode")
df.diff.p <- merge(df.diff.p, mean_by_formation,
                 by.x = "formation",
                 by.y = "formation")
df.diff.p$formation <- factor(df.diff.p$formation,
                            levels = c("NKLS", "NKBS",
                                       "Tewkesbury", "Upper Kai-Iwi",
                                       "Tainui", "SHCSBSB", "modern"))

ggplot(data = df.diff.p) +
    geom_point(aes(x = age.range, y = angle_diff_Pmax_to_z)) + 
    plot.theme +
    scale_x_continuous(name = "Age Range (Ma)") +
    scale_y_continuous(name = "Angle away from Gmax")

p.ang_pmax <- ggplot(df.diff.p) +
    geom_point(aes(x = formation, y = angle_diff_Pmax_to_z),
               size = 5, shape = 17) + 
    plot.theme +
    #scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle of Phenotypic change from Pmax",
                       lim = c(0, 90))

ggsave(p.ang_pmax, 
       file = "./Results/angle.pmax.z.diff.png", 
       width = 14, height = 10, units = "cm")


p.ang_p <- ggplot(diff_between_Ps) +
    geom_point(aes(x = angle.diff_Ps.time, y = angle_diff_Ps),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between P matrices", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.ang_p, 
       file = "./Results/angle.p.diff.png", 
       width = 20, height = 20, units = "cm")

p.dot.prod_p <- ggplot(diff_between_Ps) +
    geom_point(aes(x = angle.diff_Ps.time, y = as.numeric(corr.diff_Ps)),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between P matrices") + 
    plot.theme

#### COMPARE P AND G -----
##### SIZE OF P AND G ACROSS TIME AND ANGLE -----
P_size <- c(sum(diag(P_ext_NKBS)),
            sum(diag(P_ext_NKLS)),
            sum(diag(P_ext_tewk)),
            sum(diag(P_ext_uki)),
            sum(diag(P_ext_tai)),
            sum(diag(P_ext_SHCSBSB)),
            sum(diag(P_ext_mod)))

G_size <- c(sum(diag(G_ext_NKBS)),
            sum(diag(G_ext_NKLS)),
            sum(diag(G_ext_tewk)),
            sum(diag(G_ext_uki)),
            sum(diag(G_ext_tai)),
            sum(diag(G_ext_SHCSBSB)),
            sum(diag(G_ext_mod)))

form.name <- as.character(form.df$formationCode)

P_G_size <- as.data.frame(cbind(form.name, P_size, G_size))
P_G_size$form.name <- factor(P_G_size$form.name, levels = c("NKLS", "NKBS",
                                                            "Tewkesbury",
                                                            "Upper Kai-Iwi",
                                                            "Tainui",
                                                            "SHCSBSB",
                                                            "modern"))

#proportion G to P
P_G_size$prop.g.p <- as.numeric(P_G_size$G_size)/as.numeric(P_G_size$P_size)
range(P_G_size$prop.g.p)

ggplot(P_G_size) +
    geom_point(aes(y = as.numeric(P_size), x = form.name),
               size = 5, shape = 17,
               color = "black") + 
    geom_point(aes(y = as.numeric(G_size), x = form.name),
               size = 5, shape = 15, 
               color = "black") + 
    plot.theme +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Sizes of P and G",
                       limits = c(0, .25))

#add direction
#first PC value

P_dir <- c(p.eig_variances[[1]][1],
           p.eig_variances[[2]][1],
           p.eig_variances[[3]][1],
           p.eig_variances[[4]][1],
           p.eig_variances[[5]][1],
           p.eig_variances[[6]][1],
           p.eig_variances[[7]][1])

G_dir <- c(g.eig_variances[[1]][1],
           g.eig_variances[[2]][1],
           g.eig_variances[[3]][1],
           g.eig_variances[[4]][1],
           g.eig_variances[[5]][1],
           g.eig_variances[[6]][1],
           g.eig_variances[[7]][1])

P_G_dir <- as.data.frame(cbind(P_G_size, P_dir, G_dir))

ggplot(P_G_dir) +
    geom_point((aes(y = as.numeric(P_size), x = P_dir)),
               size = 5, shape = 17,
               color = col.form) + 
    geom_point((aes(y = as.numeric(G_size), x = G_dir)),
               size = 5, shape = 15,
               color = col.form) + 
    plot.theme +
    scale_x_continuous(name = "PC 1 of P and G",
                       limits = c(0, .15)) +
    scale_y_continuous(name = "Sizes of P and G",
                       limits = c(0, .25))

#### LOOK AT TRENDS AS A FUNCTION OF TEMPERATURE ----

df.form.pc <- merge(x = P_G_dir , form.meta,
                    by.x = "form.name",
                    by.y = "formationCode",
                    all.x = TRUE,
                    all.y = FALSE)

bottom = as.numeric(df.form.pc$Isotope_Stage_Start)
top = as.numeric(df.form.pc$Isotope_Stage_End)
df.form.pc$med.O18 <- c()
df.form.pc$sd.med.O18 <- c()
df.form.pc$n.O18 <- c()
for (i in 1:nrow(df.form.pc)){
    temp = oxy.18$d18O[which(oxy.18$Time <= bottom[i] & oxy.18$Time >= top[i])]
    df.form.pc$med.O18[i] = median(temp)
    df.form.pc$sd.med.O18[i] = sd(temp)
    df.form.pc$n.O18[i] <- length(temp)
}

p.pc1.temp <- ggplot(df.form.pc[-1,]) +
    geom_point(aes(x = med.O18,
                   y = as.numeric(P_dir))) +
    geom_smooth(aes(x = med.O18,
                    y = as.numeric(P_dir)),
                method = "lm") +
    plot.theme +
    scale_x_continuous(name = expression(Median~delta^18~O)) +
    scale_y_continuous(name = expression(PC1~of~P[mat])) +
    ggtitle("PC1 as a function of temperature without modern")

ggsave(p.pc1.temp, 
       file = "./Results/temp.pc.png", 
       width = 20, height = 20, units = "cm")

summary(lm(as.numeric(df.form.pc$P_dir[-1]) ~ df.form.pc$med.O18[-1]))
#slope = 0.01383, p-value = 0.4897, r2 = 0

#### GLOBAL G ----
##### CORR OF GLOBAL G TO P -----
#NKLS
NKLS_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[1]])) #need at least
NKLS_corr_mat.glob = NKLS_comp_mat.glob$correlations + t(NKLS_comp_mat.glob$correlations) 
diag(NKLS_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(NKLS_corr_mat.glob, upper = "number", lower = "pie")

#NKBS
NKBS_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[2]])) #need at least
NKBS_corr_mat.glob = NKBS_comp_mat.glob$correlations + t(NKBS_comp_mat.glob$correlations) 
diag(NKBS_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(NKBS_corr_mat.glob, upper = "number", lower = "pie")

#Tewkesbury
Tewkesbury_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[3]])) #need at least
Tewkesbury_corr_mat.glob = Tewkesbury_comp_mat.glob$correlations + t(Tewkesbury_comp_mat.glob$correlations) 
diag(Tewkesbury_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Tewkesbury_corr_mat.glob, upper = "number", lower = "pie")

#Upper Kai-Iwi
UKai_Iwi_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[4]])) #need at least
UKai_Iwi_corr_mat.glob = UKai_Iwi_comp_mat.glob$correlations + t(UKai_Iwi_comp_mat.glob$correlations) 
diag(UKai_Iwi_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(UKai_Iwi_corr_mat.glob, upper = "number", lower = "pie")

#Tainui
Tainui_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[5]])) #need at least
Tainui_corr_mat.glob = Tainui_comp_mat.glob$correlations + t(Tainui_comp_mat.glob$correlations) 
diag(Tainui_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(Tainui_corr_mat.glob, upper = "number", lower = "pie")

#SHCSBSB
SHCSBSB_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[6]])) #need at least
SHCSBSB_corr_mat.glob = SHCSBSB_comp_mat.glob$correlations + t(SHCSBSB_comp_mat.glob$correlations) 
diag(SHCSBSB_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(SHCSBSB_corr_mat.glob, upper = "number", lower = "pie")

#modern
modern_comp_mat.glob = RandomSkewers(list(Glob_ext, P_ext[[7]])) #need at least
modern_corr_mat.glob = modern_comp_mat.glob$correlations + t(modern_comp_mat.glob$correlations) 
diag(modern_corr_mat.glob) = 1
paste("Random Skewers similarity matrix")
corrplot.mixed(modern_corr_mat.glob, upper = "number", lower = "pie")

corr.p.glob <- c(NKLS_corr_mat.glob[1,2], NKBS_corr_mat.glob[1,2],
                 Tewkesbury_corr_mat.glob[1,2], UKai_Iwi_corr_mat.glob[1,2], 
                 Tainui_corr_mat.glob[1,2], SHCSBSB_corr_mat.glob[1,2],
                 modern_corr_mat.glob[1,2])

corr.p.glob.form <- cbind(levels(formation_list), corr.p.glob)
colnames(corr.p.glob.form) <- c("formation", "correlation")
corr.p.glob.form

write.csv(corr.p.glob.form,
          "Results/correlation.p.glob.w.modern.no.wai.csv",
          row.names = FALSE)

##### EVOLVABILITY -----

###### OBSERVED EVOLVABILITY ------
### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
observed_evolvability_in_direction_of_change_glob_t1 <- t(evolved_difference_unit_length_t1)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t1
observed_evolvability_in_direction_of_change_glob_t2 <- t(evolved_difference_unit_length_t2)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t2
observed_evolvability_in_direction_of_change_glob_t3 <- t(evolved_difference_unit_length_t3)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t3
observed_evolvability_in_direction_of_change_glob_t4 <- t(evolved_difference_unit_length_t4)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t4
observed_evolvability_in_direction_of_change_glob_t5 <- t(evolved_difference_unit_length_t5)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t5
observed_evolvability_in_direction_of_change_glob_t6 <- t(evolved_difference_unit_length_t6)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t6

###### OBSERVED CONDITIONAL EVOLVABILITY ------
### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
observed_conditional_evolvability_in_direction_of_change_glob_t1 <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t1)
observed_conditional_evolvability_in_direction_of_change_glob_t2 <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t2)
observed_conditional_evolvability_in_direction_of_change_glob_t3 <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t3)
observed_conditional_evolvability_in_direction_of_change_glob_t4 <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t4)
observed_conditional_evolvability_in_direction_of_change_glob_t5 <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t5)
observed_conditional_evolvability_in_direction_of_change_glob_t6 <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t6)

### Generate 10,000 selection gradients in random directions in the n-dimensional space
#n_dimensions <- 8 # number of traits in G matrix
#Beta <- randomBeta(10000, n_dimensions)

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
X_glob <- evolvabilityBeta(as.matrix(Glob_ext), Beta)
sumX_glob <- summary(X_glob) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix
sumX_glob

X_sum_glob <- data.frame(c.mean = c(sumX_glob$Averages[[3]], rep("", 6)),
                         c.min = c(sumX_glob$Minimum[[3]], rep("", 6)),
                         c.max = c(sumX_glob$Maximum[[3]], rep("", 6)),
                         e.mean = c(sumX_glob$Averages[[1]], rep("", 6)),
                         e.min = c(sumX_glob$Minimum[[1]], rep("", 6)),
                         e.max = c(sumX_glob$Maximum[[1]], rep("", 6)),
                         observed_e_glob = c(observed_evolvability_in_direction_of_change_glob_t1,
                                             observed_evolvability_in_direction_of_change_glob_t2,
                                             observed_evolvability_in_direction_of_change_glob_t3,
                                             observed_evolvability_in_direction_of_change_glob_t4,
                                             observed_evolvability_in_direction_of_change_glob_t5,
                                             observed_evolvability_in_direction_of_change_glob_t6,
                                             ""),
                         observed_c_glob = c(observed_conditional_evolvability_in_direction_of_change_glob_t1,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t2,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t3,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t4,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t5,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t6,
                                             ""),
                    row.names = levels(formation_list))

#NO NEGATIVE VALUES!

write.csv(X_sum_glob,
          "./Results/evolvability.global.summary.no.wai.csv")

## PLOT
X_sum_glob$formation <- rownames(X_sum_glob)
X_sum_glob$formation <- factor(X_sum_glob$formation, 
                               levels = c("NKLS", "NKBS",
                                          "Tewkesbury", 
                                          "Upper Kai-Iwi", "Tainui",
                                          "SHCSBSB", "modern"))

X_sum_glob$form.trans <- formation_transition
X_sum_glob$form.trans <- factor(X_sum_glob$form.trans,
                                levels = c("NKLS to NKBS", 
                                           "NKBS to Tewkesbury",
                                           "Tewkesbury to Upper Kai-Iwi",
                                           "Upper Kai-Iwi to Tainui", 
                                           "Tainui to SHCSBSB",
                                           "SHCSBSB to modern"))

X_sum_glob.trim <- X_sum_glob[1:6,]
p.evol_glob <- ggplot(X_sum_glob.trim, aes(x = form.trans)) +
    geom_hline(yintercept = as.numeric(X_sum_glob.trim$e.min[1]),
               color = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = as.numeric(X_sum_glob.trim$e.max[1]),
               color = "darkgray", linetype = "dashed") +
    geom_hline(yintercept = as.numeric(X_sum_glob.trim$e.mean[1])) +
    geom_point(aes(y = as.numeric(observed_e_glob)),
               size = 5, shape = 18) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Evolvability") +
    plot.theme

ggsave(p.evol_glob, 
       file = "./Results/globalG.evolvability.w.modern.no.wai.png", 
       width = 14, height = 10, units = "cm")

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO GLOBAL GMAX -----
Glob_ext_pos = round(as.matrix(Glob_ext), 6)
is.symmetric.matrix(Glob_ext_pos)
is.positive.definite(Glob_ext_pos)

#no differences in G matrices based on size
#if not due to sample size, then expect results to be similar to using individual Gs

### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 
Gmax_glob <- eigen(Glob_ext)$vectors[,1]

# Put Gmax to norm length
Gmax_glob_norm <- f.normalize_vector(Gmax_glob)

# Calculate the dot product of the unit vectors
dot_product.glob_Gmax_t1 <- sum(Gmax_glob_norm * evolved_difference_unit_length_t1) #0.12008
# Calculate the angle in radians
angle_radians.glob_Gmax_t1 <- acos(dot_product.glob_Gmax_t1)
# Convert the angle to degrees
angle_degrees.glob_Gmax_t1 <- angle_radians.glob_Gmax_t1 * (180 / pi)
#83.10328

# Calculate the dot product of the unit vectors
dot_product.glob_Gmax_t2 <- sum(Gmax_glob_norm * evolved_difference_unit_length_t2) #0.5914551
# Calculate the angle in radians
angle_radians.glob_Gmax_t2 <- acos(dot_product.glob_Gmax_t2)
# Convert the angle to degrees
angle_degrees.glob_Gmax_t2 <- angle_radians.glob_Gmax_t2 * (180 / pi)
#53.73966

# Calculate the dot product of the unit vectors
dot_product.glob_Gmax_t3 <- sum(Gmax_glob_norm * evolved_difference_unit_length_t3) #-0.9511204
# Calculate the angle in radians
angle_radians.glob_Gmax_t3 <- acos(dot_product.glob_Gmax_t3)
# Convert the angle to degrees
angle_degrees.glob_Gmax_t3 <- angle_radians.glob_Gmax_t3 * (180 / pi)
#162.0118; 17.9882

# Calculate the dot product of the unit vectors
dot_product.glob_Gmax_t4 <- sum(Gmax_glob_norm * evolved_difference_unit_length_t4) #-0.8357084
# Calculate the angle in radians
angle_radians.glob_Gmax_t4 <- acos(dot_product.glob_Gmax_t4)
# Convert the angle to degrees
angle_degrees.glob_Gmax_t4 <- angle_radians.glob_Gmax_t4 * (180 / pi)
#146.6897; 33.3103

# Calculate the dot product of the unit vectors
dot_product.glob_Gmax_t5 <- sum(Gmax_glob_norm * evolved_difference_unit_length_t5) #0.7010653
# Calculate the angle in radians
angle_radians.glob_Gmax_t5 <- acos(dot_product.glob_Gmax_t5)
# Convert the angle to degrees
angle_degrees.glob_Gmax_t5 <- angle_radians.glob_Gmax_t5 * (180 / pi)
#45.48746

# Calculate the dot product of the unit vectors
dot_product.glob_Gmax_t6 <- sum(Gmax_glob_norm * evolved_difference_unit_length_t6) #0.01431066
# Calculate the angle in radians
angle_radians.glob_Gmax_t6 <- acos(dot_product.glob_Gmax_t6)
# Convert the angle to degrees
angle_degrees.glob_Gmax_t6 <- angle_radians.glob_Gmax_t6 * (180 / pi)
#89.18003

corr.diff_Glob_max_to_z <- c(dot_product.Glob_max_NKLS, dot_product.Glob_max_NKBS,
                             dot_product.Glob_max_tewk, dot_product.Glob_max_uki,
                             dot_product.Glob_max_tai, dot_product.Glob_max_SHCSBSB, "")

angle_diff_Glob_max_to_z <- c(angle_degrees.Glob_max_NKLS, angle_degrees.Glob_max_NKBS,
                              angle_degrees.Glob_max_tewk,
                              angle_degrees.Glob_max_uki, angle_degrees.Glob_max_tai,
                              angle_degrees.Glob_max_SHCSBSB, "")
diff_between_Glob_max_z <- as.data.frame(cbind(levels(formation_list), corr.diff_Glob_max_to_z, angle_diff_Glob_max_to_z))
colnames(diff_between_Glob_max_z) <- c("formation", "corr.diff_Glob_max_to_P", "angle_diff_Glob_max_to_P")
diff_between_Glob_max_P$angle_diff_Glob_max_to_z <- as.numeric(diff_between_Glob_max_z$angle_diff_Glob_max_to_z)

for(i in 1:nrow(diff_between_Glob_max_z)){
    if(isTRUE(diff_between_Glob_max_z$angle_diff_Glob_max_to_z[i] > 90)){
        diff_between_Glob_max_z$angle_diff_Glob_max_to_z[i] <- 180 - as.numeric(diff_between_Glob_max_P$angle_diff_Glob_max_to_z[i])
    }
    else{
        next
    }
}

write.csv(diff_between_Glob_max_P,
          "./Results/differences.between.Glob_max.P.w.modern.csv",
          row.names = FALSE)

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

E_ext_NKLS = round(as.matrix(E_ext[[1]]), 6) # The G matrix estimated for sample/formation 1
E_ext_NKBS = round(as.matrix(E_ext[[2]]), 6) # The G matrix estimated for sample/formation 2
E_ext_tewk = round(as.matrix(E_ext[[3]]), 6) # The G matrix estimated for sample/formation 3
E_ext_uki = round(as.matrix(E_ext[[5]]), 6) # The G matrix estimated for sample/formation 5
E_ext_tai = round(as.matrix(E_ext[[6]]), 6) # The G matrix estimated for sample/formation 6
E_ext_SHCSBSB = round(as.matrix(E_ext[[7]]), 6) # The G matrix estimated for sample/formation 7
E_ext_mod = round(as.matrix(E_ext[[8]]), 6) # The G matrix estimated for sample/formation 7

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

#### AMOUNT OF VARIATION E EXPLAINS ----
E_size <- c(sum(diag(E_ext_NKBS)),
            sum(diag(E_ext_NKLS)),
            sum(diag(E_ext_tewk)),
            sum(diag(E_ext_uki)),
            sum(diag(E_ext_tai)),
            sum(diag(E_ext_SHCSBSB)),
            sum(diag(E_ext_mod)))

P_G_E_size <- as.data.frame(cbind(form.name, P_size, G_size, E_size))
P_G_E_size$form.name <- factor(P_E_size$form.name, levels = c("NKLS", "NKBS",
                                                            "Tewkesbury", 
                                                            "Upper Kai-Iwi",
                                                            "Tainui",
                                                            "SHCSBSB",
                                                            "modern"))

#proportion G to P
P_G_E_size$prop.e.p <- as.numeric(P_E_size$E_size)/as.numeric(P_E_size$P_size)
range(P_G_E_size$prop.e.p)
#38 to 68 %

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
#slope = -0.077, p-value = 0.93, r2 = 0

ggsave(p.pg.size, 
       file = "./Results/p.g.size.png", 
       width = 20, height = 20, units = "cm")

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
#slope = 0.70, p-value = 0.005655, r2 = 0.70

ggsave(p.pe.size, 
       file = "./Results/p.e.size.png", 
       width = 20, height = 20, units = "cm")

#### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO EMAX ----

### See if change is in direction of G max
## use Gmax of t1 and compare to ∆z
# Calculate the dot product of the unit vectors
dot_product.Emax_NKLS <- sum(Emax_NKLS_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Emax_NKLS <- acos(dot_product.Emax_NKLS)
# Convert the angle to degrees
angle_degrees.Emax_NKLS <- angle_radians.Emax_NKLS * (180 / pi)
#153.7483; Gmax is 80.90905

# Calculate the dot product of the unit vectors
dot_product.Emax_NKBS <- sum(Emax_NKBS_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Emax_NKBS <- acos(dot_product.Emax_NKBS)
# Convert the angle to degrees
angle_degrees.Emax_NKBS <- angle_radians.Emax_NKBS * (180 / pi)
#40.60758; Gmax is 43.14475

# Calculate the dot product of the unit vectors
dot_product.Emax_tewk <- sum(Emax_tewk_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Emax_tewk <- acos(dot_product.Emax_tewk)
# Convert the angle to degrees
angle_degrees.Emax_tewk <- angle_radians.Emax_tewk * (180 / pi)
#116.5479; Gmax is 122.1798

# Calculate the dot product of the unit vectors
dot_product.Emax_uki <- sum(Emax_uki_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Emax_uki <- acos(dot_product.Emax_uki)
# Convert the angle to degrees
angle_degrees.Emax_uki <- angle_radians.Emax_uki * (180 / pi)
#120.9149; Gmax is 145.7685

# Calculate the dot product of the unit vectors
dot_product.Emax_tai <- sum(Emax_tai_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Emax_tai <- acos(dot_product.Emax_tai)
# Convert the angle to degrees
angle_degrees.Emax_tai <- angle_radians.Emax_tai * (180 / pi)
#61.79681; Gmax is 40.81748

# Calculate the dot product of the unit vectors
dot_product.Emax_SHCSBSB <- sum(Emax_SHCSBSB_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Emax_SHCSBSB <- acos(dot_product.Emax_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Emax_SHCSBSB <- angle_radians.Emax_SHCSBSB * (180 / pi)
#100.5322; Gmax is 88.28694

angle_diff_Emax_to_z <- c(angle_degrees.Emax_NKLS, angle_degrees.Emax_NKBS,
                          angle_degrees.Emax_tewk, 
                          angle_degrees.Emax_uki, angle_degrees.Emax_tai,
                          angle_degrees.Emax_SHCSBSB, "")
angle_diff_between_Emax_z <- as.data.frame(cbind(levels(formation_list), angle_diff_Emax_to_z))
colnames(angle_diff_between_Emax_z) <- c("formation", "angle_diff_Emax_to_P")
angle_diff_between_Emax_P$angle_diff_Emax_to_z <- as.numeric(angle_diff_between_Emax_P$angle_diff_Emax_to_z)

for(i in 1:nrow(angle_diff_between_Emax_z)){
    if(isTRUE(angle_diff_between_Emax_z$angle_diff_Emax_to_z[i] > 90)){
        angle_diff_between_Emax_z$angle_diff_Emax_to_z[i] <- 180 - as.numeric(angle_diff_between_Emax_z$angle_diff_Emax_to_z[i])
    }
    else{
        next
    }
}

write.csv(angle_diff_between_Emax_z,
          "./Results/angle.differences.between.Emax.z.w.modern.csv",
          row.names = FALSE)
