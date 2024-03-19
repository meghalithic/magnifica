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
# -  

#### LOAD DATA ----

source("./Scripts/0-env.R")

load(file = "./Results/sum.data.list.w.modern.RData") #load the g matrices calculated above 
mean_by_formation <- sum.data.list[[1]]
mean_by_formation_colony <- sum.data.list[[2]]
means <- sum.data.list[[3]]

load(file = "./Results/data.list.w.modern.RData") #load the g matrices calculated above 

Pmat <- data.list[[1]]
Gmat <- data.list[[2]]
P_ext <- data.list[[3]]
G_ext <- data.list[[4]]
p.eig_variances <- data.list[[5]]
g.eig_variances <- data.list[[6]]
df <- data.list[[7]]
dat_lg_N <- data.list[[8]]
form_data <- data.list[[9]]
by_form.n <- data.list[[10]]
col_form.n <- data.list[[11]]

load(file="./Results/global_ext.w.modern.RData") #load the g matrices calculated above 
Glob_ext

load(file = "./Results/model_G.w.modern.RData") #load the g matrices calculated above 
g.model <- model_G

#### CORR OF G & P ----

#Gmat
#Pmat
df$formation <- factor(df$formation, levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                "Upper Kai-Iwi", "Tainui", "SHCSBSB", "modern")) 
formation_list <- unique(df$formation)
formation_list

##### G MATRICES THROUGH TIME -----
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
          "Results/correlation.p.g.w.modern.csv",
          row.names = FALSE)

##### GLOBAL G THROUGH TIME -----
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
          "Results/correlation.p.glob.w.modern.csv",
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
Glob_ext_pos = round(as.matrix(Glob_ext), 6)

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

is.symmetric.matrix(Glob_ext_pos)
is.positive.definite(Glob_ext_pos)

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

#### EVOLVABILITY ----

# Function to convert vector to norm length 
f.normalize_vector <- function(vector) {
    norm_length <- sqrt(sum(vector^2))
    normalized_vector <- vector / norm_length
    return(normalized_vector)
}

#trait means by time
mean_by_formation
#order of formations:
# "NKLS"  "NKBS" "Tewkesbury" "Upper Kai-Iwi" "Tainui" "SHCSBSB"

### Calculate the vector that defines the observed divergence between sample/formation 1 an 2

NKLS <- as.numeric(mean_by_formation[1, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 1 
NKBS <- as.numeric(mean_by_formation[2, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 2 
tewk <- as.numeric(mean_by_formation[3, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 3
uki <- as.numeric(mean_by_formation[5, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 5 
tai <- as.numeric(mean_by_formation[6, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 6 
SHCSBSB <- as.numeric(mean_by_formation[7, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 7
mod <- as.numeric(mean_by_formation[8, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 7

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
                              sumX_t7$Minimum[[3]],),
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

X_sum$form.trans <- c("NKLS to NKBS", 
                      "NKBS to Tewkesbury",
                      "Tewkesbury to Upper Kai-Iwi",
                      "Upper Kai-Iwi to Tainui", 
                      "Tainui to SHCSBSB",
                      "SHCSBSB to modern",
                      "")
X_sum$form.trans <- factor(X_sum$form.trans, levels = c("NKLS to NKBS", 
                                                        "NKBS to Tewkesbury",
                                                        "Tewkesbury to Upper Kai-Iwi",
                                                        "Upper Kai-Iwi to Tainui", 
                                                        "Tainui to SHCSBSB",
                                                        "SHCSBSB to modern",
                                                        ""))

X_sum.trim <- X_sum[1:7,]
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
               size = 5, shape = 15, color = "black") +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Evolvability") +
    plot.theme
        
ggsave(p.evol, 
       file = "./Results/evolvability.w.modern.png", 
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
dot_product.Gmax_NKLS_NKBS <- sum(Gmax_NKLS_norm * Gmax_NKBS_norm)
# Calculate the angle in radians
angle_radians.Gmax_NKLS_NKBS <- acos(dot_product.Gmax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_NKBS <- angle_radians.Gmax_NKLS_NKBS * (180 / pi)
#7.792976

dot_product.Gmax_NKBS_tewk <- sum(Gmax_NKBS_norm * Gmax_tewk_norm)
# Calculate the angle in radians
angle_radians.Gmax_NKBS_tewk <- acos(dot_product.Gmax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_tewk <- angle_radians.Gmax_NKBS_tewk * (180 / pi)
#6.538487

dot_product.Gmax_tewk_uki <- sum(Gmax_tewk_norm * Gmax_uki_norm)
# Calculate the angle in radians
angle_radians.Gmax_tewk_uki <- acos(dot_product.Gmax_tewk_uki)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_uki <- angle_radians.Gmax_tewk_uki * (180 / pi)
#97.61231

dot_product.Gmax_uki_tai <- sum(Gmax_uki_norm * Gmax_tai_norm)
# Calculate the angle in radians
angle_radians.Gmax_uki_tai <- acos(dot_product.Gmax_uki_tai)
# Convert the angle to degrees
angle_degrees.Gmax_uki_tai <- angle_radians.Gmax_uki_tai * (180 / pi)
#31.89257

dot_product.Gmax_tai_SHCSBSB <- sum(Gmax_tai_norm * Gmax_SHCSBSB_norm)
# Calculate the angle in radians
angle_radians.Gmax_tai_SHCSBSB <- acos(dot_product.Gmax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Gmax_tai_SHCSBSB <- angle_radians.Gmax_tai_SHCSBSB * (180 / pi)
#28.97957

dot_product.Gmax_SHCSBSB_mod <- sum(Gmax_SHCSBSB_norm * Gmax_mod_norm)
# Calculate the angle in radians
angle_radians.Gmax_SHCSBSB_mod <- acos(dot_product.Gmax_SHCSBSB_mod)
# Convert the angle to degrees
angle_degrees.Gmax_SHCSBSB_mod <- angle_radians.Gmax_SHCSBSB_mod * (180 / pi)
#64.5733

corr.diff_Gs <- c(dot_product.Gmax_NKLS_NKBS, dot_product.Gmax_NKBS_tewk,
                  dot_product.Gmax_tewk_uki, dot_product.Gmax_uki_tai,
                  dot_product.Gmax_tai_SHCSBSB, dot_product.Gmax_SHCSBSB_mod)

angle_diff_Gs <- c(angle_degrees.Gmax_NKLS_NKBS, angle_degrees.Gmax_NKBS_tewk,
                   angle_degrees.Gmax_tewk_uki,
                   angle_degrees.Gmax_uki_tai, angle_degrees.Gmax_tai_SHCSBSB,
                   angle_degrees.Gmax_SHCSBSB_mod)

angle.diff_Gs.time <- c("NKLS to NKBS", "NKBS to Tewkesbury",
                        "Tewkesbury to Upper Kai-Iwi",
                        "Upper Kai-Iwi to Tainui", "Tainui to SHCSBSB",
                        "SHCSBSB to modern")

diff_between_Gs <- as.data.frame(cbind(angle.diff_Gs.time, angle_diff_Gs, corr.diff_Gs))
diff_between_Gs$angle_diff_Gs <- as.numeric(diff_between_Gs$angle_diff_Gs)

for(i in 1:nrow(diff_between_Gs)){
    if(isTRUE(diff_between_Gs$angle_diff_Gs[i] > 90)){
        diff_between_Gs$angle_diff_Gs[i] <- 180 - diff_between_Gs$angle_diff_Gs[i]
    }
    else{
        next
    }
}

diff_between_Gs$angle.diff_Gs.time <- factor(diff_between_Gs$angle.diff_Gs.time,
                                             levels = c("NKLS to NKBS", 
                                                         "NKBS to Tewkesbury",
                                                         "Tewkesbury to Upper Kai-Iwi",
                                                         "Upper Kai-Iwi to Tainui", 
                                                         "Tainui to SHCSBSB",
                                                         "SHCSBSB to modern"))

write.csv(diff_between_Gs,
          "./Results/differences.between.Gs.w.modern.csv",
          row.names = FALSE)

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO GMAX -----

### See if change is in direction of G max
## use Gmax of t1 and compare to ∆z
# Calculate the dot product of the unit vectors
dot_product.Gmax_NKLS_max <- sum(Gmax_NKLS_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Gmax_NKLS_max <- acos(dot_product.Gmax_NKLS_max)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_max <- angle_radians.Gmax_NKLS_max * (180 / pi)
#80.90905

# Calculate the dot product of the unit vectors
dot_product.Gmax_NKBS_max <- sum(Gmax_NKBS_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Gmax_NKBS_max <- acos(dot_product.Gmax_NKBS_max)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_max <- angle_radians.Gmax_NKBS_max * (180 / pi)
#43.14475

# Calculate the dot product of the unit vectors
dot_product.Gmax_tewk_max <- sum(Gmax_tewk_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Gmax_tewk_max <- acos(dot_product.Gmax_tewk_max)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_max <- angle_radians.Gmax_tewk_max * (180 / pi)
#122.1798


# Calculate the dot product of the unit vectors
dot_product.Gmax_uki_max <- sum(Gmax_uki_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Gmax_uki_max <- acos(dot_product.Gmax_uki_max)
# Convert the angle to degrees
angle_degrees.Gmax_uki_max <- angle_radians.Gmax_uki_max * (180 / pi)
#145.7685

# Calculate the dot product of the unit vectors
dot_product.Gmax_tai_max <- sum(Gmax_tai_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Gmax_tai_max <- acos(dot_product.Gmax_tai_max)
# Convert the angle to degrees
angle_degrees.Gmax_tai_max <- angle_radians.Gmax_tai_max * (180 / pi)
#40.81748

# Calculate the dot product of the unit vectors
dot_product.Gmax_SHCSBSB_max <- sum(Gmax_SHCSBSB_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Gmax_SHCSBSB_max <- acos(dot_product.Gmax_SHCSBSB_max)
# Convert the angle to degrees
angle_degrees.Gmax_SHCSBSB_max <- angle_radians.Gmax_SHCSBSB_max * (180 / pi)
#88.28694

corr.diff_Gmax_to_P <- c(dot_product.Gmax_NKLS_max, dot_product.Gmax_NKBS_max,
                    dot_product.Gmax_tewk_max, dot_product.Gmax_uki_max,
                    dot_product.Gmax_tai_max, dot_product.Gmax_SHCSBSB_max, "")

angle_diff_Gmax_to_P <- c(angle_degrees.Gmax_NKLS_max, angle_degrees.Gmax_NKBS_max,
                          angle_degrees.Gmax_tewk_max,
                          angle_degrees.Gmax_uki_max, angle_degrees.Gmax_tai_max,
                          angle_degrees.Gmax_SHCSBSB_max, "")
diff_between_Gmax_P <- as.data.frame(cbind(levels(formation_list), corr.diff_Gmax_to_P, angle_diff_Gmax_to_P))
colnames(diff_between_Gmax_P) <- c("formation", "angle_diff_Gmax_to_P", "corr.diff_Gmax_to_P")
diff_between_Gmax_P$angle_diff_Gmax_to_P <- as.numeric(diff_between_Gmax_P$angle_diff_Gmax_to_P)

for(i in 1:nrow(diff_between_Gmax_P)){
    if(isTRUE(diff_between_Gmax_P$angle_diff_Gmax_to_P[i] > 90)){
        diff_between_Gmax_P$angle_diff_Gmax_to_P[i] <- 180 - as.numeric(diff_between_Gmax_P$angle_diff_Gmax_to_P[i])
    }
    else{
        next
    }
}

write.csv(diff_between_Gmax_P,
          "./Results/differences.between.Gmax.P.w.modern.csv",
          row.names = FALSE)

#### LOOK AT TRENDS AS A FUNCTION OF TIME -----
form.df <- form.meta[c(1:7, 12),] #in same order as mean_by_formation
mean_by_formation

for(i in 1:nrow(form.df)){
    form.df$mean.age[i] <- mean(form.df$Start_age[i], form.df$End_age[i], na.rm = TRUE)
}

form.df$age.range <- ""
for(i in 1:nrow(form.df)){
    form.df$age.range[i] <- form.df$Start_age[i] - form.df$End_age[i]
}
form.df$age.range <- as.numeric(form.df$age.range)

df.diff <- merge(angle_diff_between_Gmax_P, form.df,
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
    geom_point(aes(x = age.range, y = angle_diff_Gmax_to_P),
               shape = 15) + 
    plot.theme +
    scale_x_continuous(name = "Age Range (Ma)") +
    scale_y_continuous(name = "Angle away from Gmax")

p.ang_gmax <- ggplot(data = df.diff) +
    geom_point(aes(x = formation, y = angle_diff_Gmax_to_P),
               size = 5, shape = 15) + 
    plot.theme +
    #scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle of Phenotypic change from Gmax",
                       lim = c(0, 90))

ggsave(p.ang_gmax, 
       file = "./Results/angle.gmax.p.diff.w.modern.png", 
       width = 14, height = 10, units = "cm")

p.ang_g <- ggplot(angle_diff_between_Gs) +
    geom_point(aes(x = angle.diff_Gs.time, y = angle_diff_Gs),
               size = 5, shape = 15) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between G matrices", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.ang_g, 
       file = "./Results/angle.g.diff.w.modern.png", 
       width = 20, height = 20, units = "cm")

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

form.name <- form.df$formationCode
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
#31 to 57 %

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
           p.eig_variances[[7]][1],
           p.eig_variances[[8]][1])

G_dir <- c(g.eig_variances[[1]][1],
           g.eig_variances[[2]][1],
           g.eig_variances[[3]][1],
           g.eig_variances[[4]][1],
           g.eig_variances[[5]][1],
           g.eig_variances[[6]][1],
           g.eig_variances[[7]][1],
           g.eig_variances[[8]][1])

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

df.form.pc <- merge(x = PC1_P_G , form.meta,
                    by.x = "formation",
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
                   y = as.numeric(PC1_P))) +
    geom_smooth(aes(x = med.O18,
                    y = as.numeric(PC1_P)),
                method = "lm") +
    plot.theme +
    scale_x_continuous(name = expression(Median~delta^18~O)) +
    scale_y_continuous(name = expression(PC1~of~P[mat])) +
    ggtitle("PC1 as a function of temperature without modern")
    
ggsave(p.pc1.temp, 
       file = "./Results/temp.pc.without.modern.png", 
       width = 20, height = 20, units = "cm")

summary(lm(as.numeric(df.form.pc$PC1_P[-1]) ~ df.form.pc$med.O18[-1]))
#slope = 0.04, p-value = 0.53, r2 = 0


#### SUBSTITUTE P FOR G ----

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

p_X_sum$form.trans <- c("NKLS to NKBS", 
                        "NKBS to Tewkesbury",
                        "Tewkesbury to Upper Kai-Iwi",
                        "Upper Kai-Iwi to Tainui", 
                        "Tainui to SHCSBSB",
                        "SHCSBSB to modern",
                        "")
p_X_sum$form.trans <- factor(p_X_sum$form.trans, levels = c("NKLS to NKBS",
                                                            "NKBS to Tewkesbury",
                                                            "Tewkesbury to Upper Kai-Iwi",
                                                            "Upper Kai-Iwi to Tainui", 
                                                            "Tainui to SHCSBSB",
                                                            "SHCSBSB to modern",
                                                            ""))

p_X_sum.trim <- p_X_sum[1:7,]
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
       file = "./Results/p_evolvability.w.modern.png", 
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

# Put Gmax to norm length
Pmax_NKLS_norm <- f.normalize_vector(Pmax_NKLS)
Pmax_NKBS_norm <- f.normalize_vector(Pmax_NKBS)
Pmax_tewk_norm <- f.normalize_vector(Pmax_tewk)
Pmax_uki_norm <- f.normalize_vector(Pmax_uki)
Pmax_tai_norm <- f.normalize_vector(Pmax_tai)
Pmax_SHCSBSB_norm <- f.normalize_vector(Pmax_SHCSBSB)
Pmax_mod_norm <- f.normalize_vector(Pmax_mod)

# Calculate the dot product of the unit vectors
dot_product.Pmax_NKLS_NKBS <- sum(Pmax_NKLS_norm * Pmax_NKBS_norm)
# Calculate the angle in radians
angle_radians.Pmax_NKLS_NKBS <- acos(dot_product.Pmax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Pmax_NKLS_NKBS <- angle_radians.Pmax_NKLS_NKBS * (180 / pi)
#7.074267; G: 7.792976

dot_product.Pmax_NKBS_tewk <- sum(Pmax_NKBS_norm * Pmax_tewk_norm)
# Calculate the angle in radians
angle_radians.Pmax_NKBS_tewk <- acos(dot_product.Pmax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Pmax_NKBS_tewk <- angle_radians.Pmax_NKBS_tewk * (180 / pi)
#3.812412; G: 6.538487

dot_product.Pmax_tewk_uki <- sum(Pmax_tewk_norm * Pmax_uki_norm)
# Calculate the angle in radians
angle_radians.Pmax_tewk_uki <- acos(dot_product.Pmax_tewk_uki)
# Convert the angle to degrees
angle_degrees.Pmax_tewk_uki <- angle_radians.Pmax_tewk_uki * (180 / pi)
#28.94257; G: 97.61231

dot_product.Pmax_uki_tai <- sum(Pmax_uki_norm * Pmax_tai_norm)
# Calculate the angle in radians
angle_radians.Pmax_uki_tai <- acos(dot_product.Pmax_uki_tai)
# Convert the angle to degrees
angle_degrees.Pmax_uki_tai <- angle_radians.Pmax_uki_tai * (180 / pi)
#27.13413; G: 31.89257

dot_product.Pmax_tai_SHCSBSB <- sum(Pmax_tai_norm * Pmax_SHCSBSB_norm)
# Calculate the angle in radians
angle_radians.Pmax_tai_SHCSBSB <- acos(dot_product.Pmax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Pmax_tai_SHCSBSB <- angle_radians.Pmax_tai_SHCSBSB * (180 / pi)
#6.016497; G: 28.97957

dot_product.Pmax_SHCSBSB_mod <- sum(Pmax_SHCSBSB_norm * Pmax_mod_norm)
# Calculate the angle in radians
angle_radians.Pmax_SHCSBSB_mod <- acos(dot_product.Pmax_SHCSBSB_mod)
# Convert the angle to degrees
angle_degrees.Pmax_SHCSBSB_mod <- angle_radians.Pmax_SHCSBSB_mod * (180 / pi)
#23.43657; G: 64.5733

corr.diff_Ps <- c(dot_product.Pmax_NKLS_NKBS, dot_product.Pmax_NKBS_tewk,
                  dot_product.Pmax_tewk_uki, dot_product.Pmax_uki_tai,
                  dot_product.Pmax_tai_SHCSBSB, dot_product.Pmax_SHCSBSB_mod)

angle_diff_Ps <- c(angle_degrees.Pmax_NKLS_NKBS, angle_degrees.Pmax_NKBS_tewk,
                   angle_degrees.Pmax_tewk_uki,
                   angle_degrees.Pmax_uki_tai, angle_degrees.Pmax_tai_SHCSBSB,
                   angle_degrees.Pmax_SHCSBSB_mod)

angle.diff_Ps.time <- c("NKLS to NKBS", "NKBS to Tewkesbury",
                        "Tewkesbury to Upper Kai-Iwi",
                        "Upper Kai-Iwi to Tainui", "Tainui to SHCSBSB",
                        "SHCSBSB to modern")

diff_between_Ps <- as.data.frame(cbind(angle.diff_Ps.time, corr.diff_Ps, angle_diff_Ps))
diff_between_Ps$angle_diff_Ps <- as.numeric(diff_between_Ps$angle_diff_Ps)

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
                                                              "SHCSBSB to modern"))

write.csv(diff_between_Ps,
          "./Results/differences.between.Ps.w.modern.csv",
          row.names = FALSE)

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO GMAX -----

### See if change is in direction of G max
## use Gmax of t1 and compare to ∆z
# Calculate the dot product of the unit vectors
dot_product.Pmax_NKLS_max <- sum(Pmax_NKLS_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Pmax_NKLS_max <- acos(dot_product.Pmax_NKLS_max)
# Convert the angle to degrees
angle_degrees.Pmax_NKLS_max <- angle_radians.Pmax_NKLS_max * (180 / pi)
#103.5047; G: 80.90905

# Calculate the dot product of the unit vectors
dot_product.Pmax_NKBS_max <- sum(Pmax_NKBS_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Pmax_NKBS_max <- acos(dot_product.Pmax_NKBS_max)
# Convert the angle to degrees
angle_degrees.Pmax_NKBS_max <- angle_radians.Pmax_NKBS_max * (180 / pi)
#39.4028; G: 43.14475

# Calculate the dot product of the unit vectors
dot_product.Pmax_tewk_max <- sum(Pmax_tewk_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Pmax_tewk_max <- acos(dot_product.Pmax_tewk_max)
# Convert the angle to degrees
angle_degrees.Pmax_tewk_max <- angle_radians.Pmax_tewk_max * (180 / pi)
#118.6094; G: 122.1798


# Calculate the dot product of the unit vectors
dot_product.Pmax_uki_max <- sum(Pmax_uki_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Pmax_uki_max <- acos(dot_product.Pmax_uki_max)
# Convert the angle to degrees
angle_degrees.Pmax_uki_max <- angle_radians.Pmax_uki_max * (180 / pi)
#134.3997; G: 145.7685

# Calculate the dot product of the unit vectors
dot_product.Pmax_tai_max <- sum(Pmax_tai_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Pmax_tai_max <- acos(dot_product.Pmax_tai_max)
# Convert the angle to degrees
angle_degrees.Pmax_tai_max <- angle_radians.Pmax_tai_max * (180 / pi)
#54.31503; G: 40.81748

# Calculate the dot product of the unit vectors
dot_product.Pmax_SHCSBSB_max <- sum(Pmax_SHCSBSB_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Pmax_SHCSBSB_max <- acos(dot_product.Pmax_SHCSBSB_max)
# Convert the angle to degrees
angle_degrees.Pmax_SHCSBSB_max <- angle_radians.Pmax_SHCSBSB_max * (180 / pi)
#99.50977; G: 88.28694

corr.diff_Pmax_to_P <- c(dot_product.Pmax_NKLS_max, dot_product.Pmax_NKBS_max,
                         dot_product.Pmax_tewk_max, dot_product.Pmax_uki_max,
                         dot_product.Pmax_tai_max, dot_product.Pmax_SHCSBSB_max, "")

angle_diff_Pmax_to_P <- c(angle_degrees.Pmax_NKLS_max, angle_degrees.Pmax_NKBS_max,
                          angle_degrees.Pmax_tewk_max,
                          angle_degrees.Pmax_uki_max, angle_degrees.Pmax_tai_max,
                          angle_degrees.Pmax_SHCSBSB_max, "")
diff_between_Pmax_P <- as.data.frame(cbind(levels(formation_list), corr.diff_Pmax_to_P, angle_diff_Pmax_to_P))
colnames(diff_between_Pmax_P) <- c("formation", "angle_diff_Pmax_to_P", "corr.diff_Pmax_to_P")
diff_between_Pmax_P$diff_Pmax_to_P <- as.numeric(diff_between_Pmax_P$angle_diff_Pmax_to_P)

for(i in 1:nrow(diff_between_Pmax_P)){
    if(isTRUE(diff_between_Pmax_P$angle_diff_Pmax_to_P[i] > 90)){
        diff_between_Pmax_P$angle_diff_Pmax_to_P[i] <- 180 - as.numeric(diff_between_Pmax_P$angle_diff_Pmax_to_P[i])
    }
    else{
        next
    }
}

write.csv(diff_between_Pmax_P,
          "./Results/differences.between.Pmax.P.w.modern.csv",
          row.names = FALSE)

#### LOOK AT TRENDS AS A FUNCTION OF TIME -----
form.df <- form.meta[c(1:7, 12),] #in same order as mean_by_formation
mean_by_formation

for(i in 1:nrow(form.df)){
    form.df$mean.age[i] <- mean(form.df$Start_age[i], form.df$End_age[i], na.rm = TRUE)
}

form.df$age.range <- ""
for(i in 1:nrow(form.df)){
    form.df$age.range[i] <- form.df$Start_age[i] - form.df$End_age[i]
}
form.df$age.range <- as.numeric(form.df$age.range)

df.diff <- merge(angle_diff_between_Pmax_P, form.df,
                 by.x = "formation",
                 by.y = "formationCode")
df.diff <- merge(df.diff, mean_by_formation,
                 by.x = "formation",
                 by.y = "formation")
df.diff$formation <- factor(df.diff$formation,
                            levels = c("NKLS", "NKBS",
                                       "Tewkesbury", "Upper Kai-Iwi",
                                       "Tainui", "SHCSBSB", "modern"))

ggplot(data = df.diff) +
    geom_point(aes(x = age.range, y = angle_diff_Pmax_to_P)) + 
    plot.theme +
    scale_x_continuous(name = "Age Range (Ma)") +
    scale_y_continuous(name = "Angle away from Gmax")

p.ang_pmax <- ggplot(data = df.diff) +
    geom_point(aes(x = formation, y = angle_diff_Pmax_to_P),
               size = 5, shape = 17) + 
    plot.theme +
    #scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle of Phenotypic change from Pmax",
                       lim = c(0, 90))

ggsave(p.ang_pmax, 
       file = "./Results/angle.pmax.p.diff.w.modern.png", 
       width = 14, height = 10, units = "cm")


p.ang_p <- ggplot(angle_diff_between_Ps) +
    geom_point(aes(x = angle.diff_Ps.time, y = angle_diff_Ps),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between G matrices", 
                       lim = c(0, 90)) + 
    plot.theme

ggsave(p.ang_p, 
       file = "./Results/angle.p.diff.w.modern.png", 
       width = 20, height = 20, units = "cm")

#### USE GLOBAL G ----
#no differences in G matrices based on size
#if not due to sample size, then expect results to be similar to using individual Gs

###### OBSERVED EVOLVABILITY ------
### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
observed_evolvability_in_direction_of_change_glob_t1 <- t(evolved_difference_unit_length_t1)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t1
observed_evolvability_in_direction_of_change_glob_t2 <- t(evolved_difference_unit_length_t2)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t2
observed_evolvability_in_direction_of_change_glob_t3 <- t(evolved_difference_unit_length_t3)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t3
observed_evolvability_in_direction_of_change_glob_t4 <- t(evolved_difference_unit_length_t4)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t4
observed_evolvability_in_direction_of_change_glob_t5 <- t(evolved_difference_unit_length_t5)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t5
observed_evolvability_in_direction_of_change_glob_t6 <- t(evolved_difference_unit_length_t6)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t6
observed_evolvability_in_direction_of_change_glob_t7 <- t(evolved_difference_unit_length_t7)%*%as.matrix(Glob_ext)%*%evolved_difference_unit_length_t7

###### OBSERVED CONDITIONAL EVOLVABILITY ------
### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
observed_conditional_evolvability_in_direction_of_change_glob_t1 <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t1)
observed_conditional_evolvability_in_direction_of_change_glob_t2 <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t2)
observed_conditional_evolvability_in_direction_of_change_glob_t3 <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t3)
observed_conditional_evolvability_in_direction_of_change_glob_t4 <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t4)
observed_conditional_evolvability_in_direction_of_change_glob_t5 <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t5)
observed_conditional_evolvability_in_direction_of_change_glob_t6 <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t6)
observed_conditional_evolvability_in_direction_of_change_glob_t7 <- 1/(t(evolved_difference_unit_length_t7)%*%solve(as.matrix(Glob_ext))%*%evolved_difference_unit_length_t7)

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

X_sum_glob <- data.frame(c.mean = c(sumX_glob$Averages[[3]], rep("", 7)),
                         c.min = c(sumX_glob$Minimum[[3]], rep("", 7)),
                         c.max = c(sumX_glob$Maximum[[3]], rep("", 7)),
                         e.mean = c(sumX_glob$Averages[[1]], rep("", 7)),
                         e.min = c(sumX_glob$Minimum[[1]], rep("", 7)),
                         e.max = c(sumX_glob$Maximum[[1]], rep("", 7)),
                         observed_e_glob = c(observed_evolvability_in_direction_of_change_glob_t1,
                                             observed_evolvability_in_direction_of_change_glob_t2,
                                             observed_evolvability_in_direction_of_change_glob_t3,
                                             observed_evolvability_in_direction_of_change_glob_t4,
                                             observed_evolvability_in_direction_of_change_glob_t5,
                                             observed_evolvability_in_direction_of_change_glob_t6,
                                             observed_evolvability_in_direction_of_change_glob_t7,
                                             ""),
                         observed_c_glob = c(observed_conditional_evolvability_in_direction_of_change_glob_t1,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t2,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t3,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t4,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t5,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t6,
                                             observed_conditional_evolvability_in_direction_of_change_glob_t7,
                                             ""),
                    row.names = levels(formation_list))

#NO NEGATIVE VALUES!

write.csv(X_sum_glob,
          "./Results/evolvability.global.summary.csv")

## PLOT
X_sum_glob$formation <- rownames(X_sum_glob)
X_sum_glob$formation <- factor(X_sum_glob$formation, 
                               levels = c("NKLS", "NKBS",
                                          "Tewkesbury", 
                                          "Upper Kai-Iwi", "Tainui",
                                          "SHCSBSB", "modern"))

X_sum_glob$form.trans <- c("NKLS to NKBS",
                           "NKBS to Tewkesbury",
                           "Tewkesbury to Upper Kai-Iwi",
                           "Upper Kai-Iwi to Tainui", 
                           "Tainui to SHCSBSB",
                           "SHCSBSB to modern",
                           "")
X_sum_glob$form.trans <- factor(X_sum_glob$form.trans, 
                                levels = c("NKLS to NKBS", 
                                           "NKBS to Tewkesbury",
                                           "Tewkesbury to Upper Kai-Iwi",
                                           "Upper Kai-Iwi to Tainui", 
                                           "Tainui to SHCSBSB",
                                           "SHCSBSB to modern",
                                           ""))

X_sum_glob.trim <- X_sum_glob[1:7,]
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
       file = "./Results/globalG.evolvability.w.modern.png", 
       width = 14, height = 10, units = "cm")

##### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO GMAX -----
### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 
Gmax_glob <- eigen(Glob_ext)$vectors[,1]

# Put Gmax to norm length
Gmax_glob_norm <- f.normalize_vector(Gmax_glob)

### See if change is in direction of G max
# Calculate the dot product of the unit vectors
dot_product.Glob_max_NKLS <- sum(Gmax_glob_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Glob_max_NKLS <- acos(dot_product.Glob_max_NKLS)
# Convert the angle to degrees
angle_degrees.Glob_max_NKLS <- angle_radians.Glob_max_NKLS * (180 / pi)
#155.9968; was 80.90905

# Calculate the dot product of the unit vectors
dot_product.Glob_max_NKBS <- sum(Gmax_glob_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Glob_max_NKBS <- acos(dot_product.Glob_max_NKBS)
# Convert the angle to degrees
angle_degrees.Glob_max_NKBS <- angle_radians.Glob_max_NKBS * (180 / pi)
#43.31873; was 43.14475

# Calculate the dot product of the unit vectors
dot_product.Glob_max_tewk <- sum(Gmax_glob_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Glob_max_tewk <- acos(dot_product.Glob_max_tewk)
# Convert the angle to degrees
angle_degrees.Glob_max_tewk <- angle_radians.Glob_max_tewk * (180 / pi)
#126.4306; was 122.1798

# Calculate the dot product of the unit vectors
dot_product.Glob_max_uki <- sum(Gmax_glob_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Glob_max_uki <- acos(dot_product.Glob_max_uki)
# Convert the angle to degrees
angle_degrees.Glob_max_uki <- angle_radians.Glob_max_uki * (180 / pi)
#128.7347; was 145.7685

# Calculate the dot product of the unit vectors
dot_product.Glob_max_tai <- sum(Gmax_glob_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Glob_max_tai <- acos(dot_product.Glob_max_tai)
# Convert the angle to degrees
angle_degrees.Glob_max_tai <- angle_radians.Glob_max_tai * (180 / pi)
#49.55738; was 40.81748

# Calculate the dot product of the unit vectors
dot_product.Glob_max_SHCSBSB <- sum(Gmax_glob_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Glob_max_SHCSBSB <- acos(dot_product.Glob_max_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Glob_max_SHCSBSB <- angle_radians.Glob_max_SHCSBSB * (180 / pi)
#95.32893; was 88.28694

corr.diff_Glob_max_to_P <- c(dot_product.Glob_max_NKLS, dot_product.Glob_max_NKBS,
                             dot_product.Glob_max_tewk, dot_product.Glob_max_uki,
                             dot_product.Glob_max_tai, dot_product.Glob_max_SHCSBSB, "")

angle_diff_Glob_max_to_P <- c(angle_degrees.Glob_max_NKLS, angle_degrees.Glob_max_NKBS,
                              angle_degrees.Glob_max_tewk,
                              angle_degrees.Glob_max_uki, angle_degrees.Glob_max_tai,
                              angle_degrees.Glob_max_SHCSBSB, "")
diff_between_Glob_max_P <- as.data.frame(cbind(levels(formation_list), corr.diff_Glob_max_to_P, angle_diff_Glob_max_to_P))
colnames(diff_between_Glob_max_P) <- c("formation", "corr.diff_Glob_max_to_P", "angle_diff_Glob_max_to_P")
diff_between_Glob_max_P$angle_diff_Glob_max_to_P <- as.numeric(diff_between_Glob_max_P$angle_diff_Glob_max_to_P)

for(i in 1:nrow(diff_between_Glob_max_P)){
    if(isTRUE(diff_between_Glob_max_P$angle_diff_Glob_max_to_P[i] > 90)){
        diff_between_Glob_max_P$angle_diff_Glob_max_to_P[i] <- 180 - as.numeric(diff_between_Glob_max_P$angle_diff_Glob_max_to_P[i])
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
dot_product.Emax_NKLS_max <- sum(Emax_NKLS_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Emax_NKLS_max <- acos(dot_product.Emax_NKLS_max)
# Convert the angle to degrees
angle_degrees.Emax_NKLS_max <- angle_radians.Emax_NKLS_max * (180 / pi)
#153.7483; Gmax is 80.90905

# Calculate the dot product of the unit vectors
dot_product.Emax_NKBS_max <- sum(Emax_NKBS_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Emax_NKBS_max <- acos(dot_product.Emax_NKBS_max)
# Convert the angle to degrees
angle_degrees.Emax_NKBS_max <- angle_radians.Emax_NKBS_max * (180 / pi)
#40.60758; Gmax is 43.14475

# Calculate the dot product of the unit vectors
dot_product.Emax_tewk_max <- sum(Emax_tewk_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Emax_tewk_max <- acos(dot_product.Emax_tewk_max)
# Convert the angle to degrees
angle_degrees.Emax_tewk_max <- angle_radians.Emax_tewk_max * (180 / pi)
#116.5479; Gmax is 122.1798

# Calculate the dot product of the unit vectors
dot_product.Emax_uki_max <- sum(Emax_uki_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Emax_uki_max <- acos(dot_product.Emax_uki_max)
# Convert the angle to degrees
angle_degrees.Emax_uki_max <- angle_radians.Emax_uki_max * (180 / pi)
#120.9149; Gmax is 145.7685

# Calculate the dot product of the unit vectors
dot_product.Emax_tai_max <- sum(Emax_tai_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Emax_tai_max <- acos(dot_product.Emax_tai_max)
# Convert the angle to degrees
angle_degrees.Emax_tai_max <- angle_radians.Emax_tai_max * (180 / pi)
#61.79681; Gmax is 40.81748

# Calculate the dot product of the unit vectors
dot_product.Emax_SHCSBSB_max <- sum(Emax_SHCSBSB_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Emax_SHCSBSB_max <- acos(dot_product.Emax_SHCSBSB_max)
# Convert the angle to degrees
angle_degrees.Emax_SHCSBSB_max <- angle_radians.Emax_SHCSBSB_max * (180 / pi)
#100.5322; Gmax is 88.28694

angle_diff_Emax_to_P <- c(angle_degrees.Emax_NKLS_max, angle_degrees.Emax_NKBS_max,
                          angle_degrees.Emax_tewk_max, 
                          angle_degrees.Emax_uki_max, angle_degrees.Emax_tai_max,
                          angle_degrees.Emax_SHCSBSB_max, "")
angle_diff_between_Emax_P <- as.data.frame(cbind(levels(formation_list), angle_diff_Emax_to_P))
colnames(angle_diff_between_Emax_P) <- c("formation", "angle_diff_Emax_to_P")
angle_diff_between_Emax_P$angle_diff_Emax_to_P <- as.numeric(angle_diff_between_Emax_P$angle_diff_Emax_to_P)

for(i in 1:nrow(angle_diff_between_Emax_P)){
    if(isTRUE(angle_diff_between_Emax_P$angle_diff_Emax_to_P[i] > 90)){
        angle_diff_between_Emax_P$angle_diff_Emax_to_P[i] <- 180 - as.numeric(angle_diff_between_Emax_P$angle_diff_Emax_to_P[i])
    }
    else{
        next
    }
}

write.csv(angle_diff_between_Emax_P,
          "./Results/angle.differences.between.Emax.P.w.modern.csv",
          row.names = FALSE)
