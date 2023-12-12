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

load(file = "./Results/data.list.RData") #load the g matrices calculated above 

P_ext <- data.list[[1]]
G_ext <- data.list[[2]]
df <- data.list[[3]]
dat_lg_N <- data.list[[4]]
mean_by_formation <- data.list[[5]]
mean_by_formation_colony <- data.list[[6]]
form_data <- data.list[[7]]
by_form.n <- data.list[[8]]
col_form.n <- data.list[[9]]

#### CORR OF G & P ----

#Gmat
#Pmat

formation_list

## RANDOM SKEWERS OF P & G OF EACH FORMATION
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

corr.p.g <- c(NKLS_corr_mat[1,2], NKBS_corr_mat[1,2],
              Tewkesbury_corr_mat[1,2], Waipuru_corr_mat[1,2], 
              UKai_Iwi_corr_mat[1,2], Tainui_corr_mat[1,2], 
              SHCSBSB_corr_mat[1,2])

corr.p.g.form <- cbind(levels(formation_list), corr.p.g)
colnames(corr.p.g.form) <- c("formation", "correlation")
corr.p.g.form

write.csv(corr.p.g.form,
          "Results/correlation.p.g.csv",
          row.names = FALSE)


###### PC1 of P and PC1 of G through time -------
## NOTE THESE ARE NOT USING EXT MATRICES THAT ARE CORRECTED LIKE WE DO BELOW FOR
## EXAMINING CHANGE IN G THROUGH TIME

## PERCENT EXPLAINED BY PC1
PC1_P <- p.eig_per[p.eig_per$variable == "X1",]
PC1_G <- g.eig_per[g.eig_per$variable == "X1",]
PC1_P_G <- cbind(levels(formation_list), PC1_P$value, PC1_G$value)
colnames(PC1_P_G) <- c("formation", "PC1_P", "PC1_G")
PC1_P_G <- as.data.frame(PC1_P_G)

P_G_PC1 = ggplot(PC1_P_G) +
    geom_point(aes(x = as.numeric(PC1_G), y = as.numeric(PC1_P),
                   col = formation),
               size = 5, shape = 17) +
    geom_smooth(aes(x = as.numeric(PC1_G), y = as.numeric(PC1_P)),
                method = "lm") +
    xlab("% PC1 of P matrix") +
    ylab("% PC1 of G matrix") +
    scale_color_manual(values = col.form) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
P_G_PC1 #none negative; none above 1!
summary(lm(as.numeric(PC1_P_G$PC1_P) ~ as.numeric(PC1_P_G$PC1_G)))
#sig; slope = 0.51522; r2 = 0.9768

ggsave(P_G_PC1, 
       file = "./Results/PC1.P.G.png", 
       width = 14, height = 10, units = "cm")

## VALUES
p.eig_variances
g.eig_variances

p.eig_var_mat = do.call(rbind, p.eig_variances)
p.eig_var_mat = data.frame(p.eig_var_mat, rownames(p.eig_var_mat))
p.eig_var = melt(p.eig_var_mat)
colnames(p.eig_var) <- c("formation", "p.variable", "p.value")

g.eig_var_mat = do.call(rbind, g.eig_variances)
g.eig_var_mat = data.frame(g.eig_var_mat, rownames(g.eig_var_mat))
g.eig_var = melt(g.eig_var_mat)
colnames(g.eig_var) <- c("formation", "g.variable", "g.value")

p.eig_var.1 <- p.eig_var[p.eig_var$p.variable == "X1", ]
g.eig_var.1 <- g.eig_var[g.eig_var$g.variable == "X1", ]

p.g.eig_var.1 <- merge(p.eig_var.1,
                       g.eig_var.1,
                       by = "formation")

P_G_PC1.val = ggplot(p.g.eig_var.1,
                     aes(x = g.value, y = p.value,
                         group = formation,
                         colour = formation)) +
    geom_point() +
    xlab("P matrix PC1 value") +
    ylab("G matrix PC1 value")
P_G_PC1.val

#### CALCULATE E ----

##### PLOT P, E, AND G VARIATION ----
# E = units
# P = cov or units + colony id
# G = estimate or colony id

###### POSITIVE DEFINITE ------
## check positive definite
# round to 10 decimals to make it symmetric
# to make it positive definite, use extended matrix to fill in

G_ext_NKLS = round(as.matrix(G_ext[[1]]), 6) # The G matrix estimated for sample/formation 1
G_ext_NKBS = round(as.matrix(G_ext[[2]]), 6) # The G matrix estimated for sample/formation 2
G_ext_tewk = round(as.matrix(G_ext[[3]]), 6) # The G matrix estimated for sample/formation 3
G_ext_wai = round(as.matrix(G_ext[[4]]), 6) # The G matrix estimated for sample/formation 4
G_ext_uki = round(as.matrix(G_ext[[5]]), 6) # The G matrix estimated for sample/formation 5
G_ext_tai = round(as.matrix(G_ext[[6]]), 6) # The G matrix estimated for sample/formation 6
G_ext_SHCSBSB = round(as.matrix(G_ext[[7]]), 6) # The G matrix estimated for sample/formation 7

is.symmetric.matrix(G_ext_NKLS)
is.positive.definite(G_ext_NKLS)

is.symmetric.matrix(G_ext_NKBS)
is.positive.definite(G_ext_NKBS)

is.symmetric.matrix(G_ext_tewk)
is.positive.definite(G_ext_tewk)

is.symmetric.matrix(G_ext_wai)
is.positive.definite(G_ext_wai) #now TRUE

is.symmetric.matrix(G_ext_uki)
is.positive.definite(G_ext_uki)

is.symmetric.matrix(G_ext_tai)
is.positive.definite(G_ext_tai)

is.symmetric.matrix(G_ext_SHCSBSB)
is.positive.definite(G_ext_SHCSBSB)

P_ext_NKLS = round(as.matrix(P_ext[[1]]), 6) # The G matrix estimated for sample/formation 1
P_ext_NKBS = round(as.matrix(P_ext[[2]]), 6) # The G matrix estimated for sample/formation 2
P_ext_tewk = round(as.matrix(P_ext[[3]]), 6) # The G matrix estimated for sample/formation 3
P_ext_wai = round(as.matrix(P_ext[[4]]), 6) # The G matrix estimated for sample/formation 4
P_ext_uki = round(as.matrix(P_ext[[5]]), 6) # The G matrix estimated for sample/formation 5
P_ext_tai = round(as.matrix(P_ext[[6]]), 6) # The G matrix estimated for sample/formation 6
P_ext_SHCSBSB = round(as.matrix(P_ext[[7]]), 6) # The G matrix estimated for sample/formation 7

is.symmetric.matrix(P_ext_NKLS)
is.positive.definite(P_ext_NKLS)

is.symmetric.matrix(P_ext_NKBS)
is.positive.definite(P_ext_NKBS)

is.symmetric.matrix(P_ext_tewk)
is.positive.definite(P_ext_tewk)

is.symmetric.matrix(P_ext_wai)
is.positive.definite(P_ext_wai) #now TRUE

is.symmetric.matrix(P_ext_uki)
is.positive.definite(P_ext_uki)

is.symmetric.matrix(P_ext_tai)
is.positive.definite(P_ext_tai)

is.symmetric.matrix(P_ext_SHCSBSB)
is.positive.definite(P_ext_SHCSBSB)

#### COMPARE Gs THROUGH TIME ----

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

NKLS <- as.numeric(mean_by_formation[1, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 1 
NKBS <- as.numeric(mean_by_formation[2, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 2 
tewk <- as.numeric(mean_by_formation[3, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 3
wai <- as.numeric(mean_by_formation[4, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 4
uki <- as.numeric(mean_by_formation[5, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 5 
tai <- as.numeric(mean_by_formation[6, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 6 
SHCSBSB <- as.numeric(mean_by_formation[7, c(7, 11, 15, 19, 23, 27, 31, 35)]) # A vector containing trait means from sample/formation 7

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

###### OBSERVED EVOLVABILITY ------
### The evolvability in the direction of divergence from sample/formation 1 to sample/formation 2
#observed_evolvability_in_direction_of_change<-t(evolved_difference_unit_length)%*%as.matrix(G_matrix_1)%*%evolved_difference_unit_length
observed_evolvability_in_direction_of_change_t1 <- t(evolved_difference_unit_length_t1)%*%as.matrix(G_ext_NKLS)%*%evolved_difference_unit_length_t1
observed_evolvability_in_direction_of_change_t2 <- t(evolved_difference_unit_length_t2)%*%as.matrix(G_ext_NKBS)%*%evolved_difference_unit_length_t2
observed_evolvability_in_direction_of_change_t3 <- t(evolved_difference_unit_length_t3)%*%as.matrix(G_ext_tewk)%*%evolved_difference_unit_length_t3
observed_evolvability_in_direction_of_change_t4 <- t(evolved_difference_unit_length_t4)%*%as.matrix(G_ext_wai)%*%evolved_difference_unit_length_t4
observed_evolvability_in_direction_of_change_t5 <- t(evolved_difference_unit_length_t5)%*%as.matrix(G_ext_uki)%*%evolved_difference_unit_length_t5
observed_evolvability_in_direction_of_change_t6 <- t(evolved_difference_unit_length_t6)%*%as.matrix(G_ext_tai)%*%evolved_difference_unit_length_t6

###### OBSERVED CONDITIONAL EVOLVABILITY ------
### The conditional evolvability in the direction of divergence
#observed_conditional_evolvability_in_direction_of_change<-1/(t(evolved_difference_unit_length)%*%solve(as.matrix(G_matrix_1))%*%evolved_difference_unit_length)
observed_conditional_evolvability_in_direction_of_change_t1 <- 1/(t(evolved_difference_unit_length_t1)%*%solve(as.matrix(G_ext_NKLS))%*%evolved_difference_unit_length_t1)
observed_conditional_evolvability_in_direction_of_change_t2 <- 1/(t(evolved_difference_unit_length_t2)%*%solve(as.matrix(G_ext_NKBS))%*%evolved_difference_unit_length_t2)
observed_conditional_evolvability_in_direction_of_change_t3 <- 1/(t(evolved_difference_unit_length_t3)%*%solve(as.matrix(G_ext_tewk))%*%evolved_difference_unit_length_t3)
observed_conditional_evolvability_in_direction_of_change_t4 <- 1/(t(evolved_difference_unit_length_t4)%*%solve(as.matrix(G_ext_wai))%*%evolved_difference_unit_length_t4)
observed_conditional_evolvability_in_direction_of_change_t5 <- 1/(t(evolved_difference_unit_length_t5)%*%solve(as.matrix(G_ext_uki))%*%evolved_difference_unit_length_t5)
observed_conditional_evolvability_in_direction_of_change_t6 <- 1/(t(evolved_difference_unit_length_t6)%*%solve(as.matrix(G_ext_tai))%*%evolved_difference_unit_length_t6)

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

X_t4 <- evolvabilityBeta(as.matrix(G_ext_wai), Beta)
sumX_t4 <- summary(X_t4) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t5 <- evolvabilityBeta(as.matrix(G_ext_uki), Beta)
sumX_t5 <- summary(X_t5) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t6 <- evolvabilityBeta(as.matrix(G_ext_tai), Beta)
sumX_t6 <- summary(X_t6) #provides you with info on mean, minimum and maximum evolvability  (e_mean, e_min, e_max) and conditional evolvability  (c_mean, c_min, c_max) for a given G matrix

X_t7 <- evolvabilityBeta(as.matrix(G_ext_SHCSBSB), Beta)
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
                                                      "Tewkesbury", "Waipuru",
                                                      "Upper Kai-Iwi", "Tainui",
                                                      "SHCSBSB"))

X_sum$form.trans <- c("NKLS to NKBS", 
                      "NKBS to Tewkesbury",
                      "Tewkesbury to Waipuru", 
                      "Waipuru to Upper Kai-Iwi",
                      "Upper Kai-Iwi to Tainui", 
                      "Tainui to SHCSBSB",
                      "")
X_sum$form.trans <- factor(X_sum$form.trans, levels = c("NKLS to NKBS", 
                                                        "NKBS to Tewkesbury",
                                                        "Tewkesbury to Waipuru", 
                                                        "Waipuru to Upper Kai-Iwi",
                                                        "Upper Kai-Iwi to Tainui", 
                                                        "Tainui to SHCSBSB",
                                                        ""))

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
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA))
        
ggsave(p.evol, 
       file = "./Results/evolvability.png", 
       width = 14, height = 10, units = "cm")

# By comparing the evolvabilities you estimated in the direction of change (lines 9 and 12) with the average evolvabilities calculated by running line 20, you get a sense of whether evolution happened in directions with above or below average evolvability.  

#### ANGLE CHANGE IN G BETWEEN FORMATIONS ----

### Proportion of variance in n-dimensional trait space that is explained by PC1 (i.e., the first eigenvector)
#eigen(as.matrix(G_matrix_1))$values[1]/sum(eigen(as.matrix(G_matrix_1))$values)
eigen(as.matrix(G_ext_NKLS))$values[1]/sum(eigen(as.matrix(G_ext_NKLS))$values) #0.4319903
eigen(as.matrix(G_ext_NKBS))$values[1]/sum(eigen(as.matrix(G_ext_NKBS))$values) #0.4971828
eigen(as.matrix(G_ext_tewk))$values[1]/sum(eigen(as.matrix(G_ext_tewk))$values) #0.4953968
eigen(as.matrix(G_ext_wai))$values[1]/sum(eigen(as.matrix(G_ext_wai))$values) #0.4059567
eigen(as.matrix(G_ext_uki))$values[1]/sum(eigen(as.matrix(G_ext_uki))$values) #0.7550115
eigen(as.matrix(G_ext_tai))$values[1]/sum(eigen(as.matrix(G_ext_tai))$values) #0.397961
eigen(as.matrix(G_ext_SHCSBSB))$values[1]/sum(eigen(as.matrix(G_ext_SHCSBSB))$values) #0.4637151

### How much is the direction of Gmax (i.e., the direction first ) varying between different G-matrices? 
Gmax_NKLS <- eigen(G_ext_NKLS)$vectors[,1]
Gmax_NKBS <- eigen(G_ext_NKBS)$vectors[,1]
Gmax_tewk <- eigen(G_ext_tewk)$vectors[,1]
Gmax_wai <- eigen(G_ext_wai)$vectors[,1]
Gmax_uki <- eigen(G_ext_uki)$vectors[,1]
Gmax_tai <- eigen(G_ext_tai)$vectors[,1]
Gmax_SHCSBSB <- eigen(G_ext_SHCSBSB)$vectors[,1]

# Put Gmax to norm length
Gmax_NKLS_norm <- f.normalize_vector(Gmax_NKLS)
Gmax_NKBS_norm <- f.normalize_vector(Gmax_NKBS)
Gmax_tewk_norm <- f.normalize_vector(Gmax_tewk)
Gmax_wai_norm <- f.normalize_vector(Gmax_wai)
Gmax_uki_norm <- f.normalize_vector(Gmax_uki)
Gmax_tai_norm <- f.normalize_vector(Gmax_tai)
Gmax_SHCSBSB_norm <- f.normalize_vector(Gmax_SHCSBSB)

# Calculate the dot product of the unit vectors
dot_product.Gmax_NKLS_NKBS <- sum(Gmax_NKLS_norm * Gmax_NKBS_norm)
# Calculate the angle in radians
angle_radians.Gmax_NKLS_NKBS <- acos(dot_product.Gmax_NKLS_NKBS)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_NKBS <- angle_radians.Gmax_NKLS_NKBS * (180 / pi)
#5.907311

dot_product.Gmax_NKBS_tewk <- sum(Gmax_NKBS_norm * Gmax_tewk_norm)
# Calculate the angle in radians
angle_radians.Gmax_NKBS_tewk <- acos(dot_product.Gmax_NKBS_tewk)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_tewk <- angle_radians.Gmax_NKBS_tewk * (180 / pi)
#7.092774

dot_product.Gmax_tewk_wai <- sum(Gmax_tewk_norm * Gmax_wai_norm)
# Calculate the angle in radians
angle_radians.Gmax_tewk_wai <- acos(dot_product.Gmax_tewk_wai)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_wai <- angle_radians.Gmax_tewk_wai * (180 / pi)
#18.6635

dot_product.Gmax_wai_uki <- sum(Gmax_wai_norm * Gmax_uki_norm)
# Calculate the angle in radians
angle_radians.Gmax_wai_uki <- acos(dot_product.Gmax_wai_uki)
# Convert the angle to degrees
angle_degrees.Gmax_wai_uki <- angle_radians.Gmax_wai_uki * (180 / pi)
#27.97863

dot_product.Gmax_uki_tai <- sum(Gmax_uki_norm * Gmax_tai_norm)
# Calculate the angle in radians
angle_radians.Gmax_uki_tai <- acos(dot_product.Gmax_uki_tai)
# Convert the angle to degrees
angle_degrees.Gmax_uki_tai <- angle_radians.Gmax_uki_tai * (180 / pi)
#33.51788

dot_product.Gmax_tai_SHCSBSB <- sum(Gmax_tai_norm * Gmax_SHCSBSB_norm)
# Calculate the angle in radians
angle_radians.Gmax_tai_SHCSBSB <- acos(dot_product.Gmax_tai_SHCSBSB)
# Convert the angle to degrees
angle_degrees.Gmax_tai_SHCSBSB <- angle_radians.Gmax_tai_SHCSBSB * (180 / pi)
#24.61935

angle_diff_Gs <- c(angle_degrees.Gmax_NKLS_NKBS, angle_degrees.Gmax_NKBS_tewk,
                   angle_degrees.Gmax_tewk_wai, angle_degrees.Gmax_wai_uki,
                   angle_degrees.Gmax_uki_tai, angle_degrees.Gmax_tai_SHCSBSB)

angle.diff_Gs.time <- c("NKLS to NKBS", "NKBS to Tewkesbury",
                        "Tewkesbury to Waipuru", "Waipuru to Upper Kai-Iwi",
                        "Upper Kai-Iwi to Tainui", "Tainui to SHCSBSB")

angle_diff_between_Gs <- as.data.frame(cbind(angle.diff_Gs.time, angle_diff_Gs))
angle_diff_between_Gs$angle_diff_Gs <- as.numeric(angle_diff_between_Gs$angle_diff_Gs)

for(i in 1:nrow(angle_diff_between_Gs)){
    if(isTRUE(angle_diff_between_Gs$angle_diff_Gs[i] > 90)){
        angle_diff_between_Gs$angle_diff_Gs[i] <- 180 - angle_diff_between_Gs$angle_diff_Gs[i]
    }
    else{
        next
    }
}

angle_diff_between_Gs$angle.diff_Gs.time <- factor(angle_diff_between_Gs$angle.diff_Gs.time,
                                                   levels = c("NKLS to NKBS", 
                                                              "NKBS to Tewkesbury",
                                                              "Tewkesbury to Waipuru", 
                                                              "Waipuru to Upper Kai-Iwi",
                                                              "Upper Kai-Iwi to Tainui", 
                                                              "Tainui to SHCSBSB"))

write.csv(angle_diff_between_Gs,
          "./Results/angle.differences.between.Gs.csv",
          row.names = FALSE)

#### DIRECTION OF PHENOTYPIC CHANGE COMPARED TO GMAX ----

### See if change is in direction of G max
# Calculate the dot product of the unit vectors
dot_product.Gmax_NKLS_max <- sum(Gmax_NKLS_norm * evolved_difference_unit_length_t1)
# Calculate the angle in radians
angle_radians.Gmax_NKLS_max <- acos(dot_product.Gmax_NKLS_max)
# Convert the angle to degrees
angle_degrees.Gmax_NKLS_max <- angle_radians.Gmax_NKLS_max * (180 / pi)
#96.88325

# Calculate the dot product of the unit vectors
dot_product.Gmax_NKBS_max <- sum(Gmax_NKBS_norm * evolved_difference_unit_length_t2)
# Calculate the angle in radians
angle_radians.Gmax_NKBS_max <- acos(dot_product.Gmax_NKBS_max)
# Convert the angle to degrees
angle_degrees.Gmax_NKBS_max <- angle_radians.Gmax_NKBS_max * (180 / pi)
#158.9253

# Calculate the dot product of the unit vectors
dot_product.Gmax_tewk_max <- sum(Gmax_tewk_norm * evolved_difference_unit_length_t3)
# Calculate the angle in radians
angle_radians.Gmax_tewk_max <- acos(dot_product.Gmax_tewk_max)
# Convert the angle to degrees
angle_degrees.Gmax_tewk_max <- angle_radians.Gmax_tewk_max * (180 / pi)
#123.4847

# Calculate the dot product of the unit vectors
dot_product.Gmax_wai_max <- sum(Gmax_wai_norm * evolved_difference_unit_length_t4)
# Calculate the angle in radians
angle_radians.Gmax_wai_max <- acos(dot_product.Gmax_wai_max)
# Convert the angle to degrees
angle_degrees.Gmax_wai_max <- angle_radians.Gmax_wai_max * (180 / pi)
#25.64966

# Calculate the dot product of the unit vectors
dot_product.Gmax_uki_max <- sum(Gmax_uki_norm * evolved_difference_unit_length_t5)
# Calculate the angle in radians
angle_radians.Gmax_uki_max <- acos(dot_product.Gmax_uki_max)
# Convert the angle to degrees
angle_degrees.Gmax_uki_max <- angle_radians.Gmax_uki_max * (180 / pi)
#154.1939

# Calculate the dot product of the unit vectors
dot_product.Gmax_tai_max <- sum(Gmax_tai_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Gmax_tai_max <- acos(dot_product.Gmax_tai_max)
# Convert the angle to degrees
angle_degrees.Gmax_tai_max <- angle_radians.Gmax_tai_max * (180 / pi)
#34.34499

# Calculate the dot product of the unit vectors
dot_product.Gmax_SHCSBSB_max <- sum(Gmax_SHCSBSB_norm * evolved_difference_unit_length_t6)
# Calculate the angle in radians
angle_radians.Gmax_SHCSBSB_max <- acos(dot_product.Gmax_SHCSBSB_max)
# Convert the angle to degrees
angle_degrees.Gmax_SHCSBSB_max <- angle_radians.Gmax_SHCSBSB_max * (180 / pi)
#12.43953

angle_diff_Gmax_to_P <- c(angle_degrees.Gmax_NKLS_max, angle_degrees.Gmax_NKBS_max,
                          angle_degrees.Gmax_tewk_max, angle_degrees.Gmax_wai_max,
                          angle_degrees.Gmax_uki_max, angle_degrees.Gmax_tai_max,
                          angle_degrees.Gmax_SHCSBSB_max)
angle_diff_between_Gmax_P <- as.data.frame(cbind(levels(formation_list), angle_diff_Gmax_to_P))
colnames(angle_diff_between_Gmax_P) <- c("formation", "angle_diff_Gmax_to_P")
angle_diff_between_Gmax_P$angle_diff_Gmax_to_P <- as.numeric(angle_diff_between_Gmax_P$angle_diff_Gmax_to_P)

for(i in 1:nrow(angle_diff_between_Gmax_P)){
    if(isTRUE(angle_diff_between_Gmax_P$angle_diff_Gmax_to_P[i] > 90)){
        angle_diff_between_Gmax_P$angle_diff_Gmax_to_P[i] <- 180 - as.numeric(angle_diff_between_Gmax_P$angle_diff_Gmax_to_P[i])
    }
    else{
        next
    }
}

write.csv(angle_diff_between_Gmax_P,
          "./Results/angle.differences.between.Gmax.P.csv",
          row.names = FALSE)

#### LOOK AT TRENDS AS A FUNCTION OF TIME -----
form.df <- form.meta[1:7,] #in same order as mean_by_formation
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

ggplot(data = df.diff) +
    geom_point(aes(x = age.range, y = angle_diff_Gmax_to_P)) + 
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_x_continuous(name = "Age Range (Ma)") +
    scale_y_continuous(name = "Angle away from Gmax")

df.diff$formation <- factor(df.diff$formation,
                            levels = c("NKLS", "NKBS",
                                       "Tewkesbury", "Waipuru",
                                       "Upper Kai-Iwi", "Tainui",
                                       "SHCSBSB"))

p.ang_gmax <- ggplot(data = df.diff) +
    geom_point(aes(x = formation, y = angle_diff_Gmax_to_P),
               size = 5, shape = 17) + 
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA)) +
    #scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle of Phenotypic change from Gmax",
                       lim = c(0, 90))

ggsave(p.ang_gmax, 
       file = "./Results/angle.gmax.p.diff.png", 
       width = 14, height = 10, units = "cm")


ggplot(data = df.diff) +
    geom_point(aes(x = mean.age, y = avg.zh,
                   col = formationName)) + 
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
    scale_y_continuous(name = "Mean LN ZH") +
    scale_color_manual(values = col.form)

ggplot(data = df.diff) +
    geom_point(aes(x = age.range, y = avg.zh,
                   col = col.form)) + 
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_x_continuous(name = "Age Range (Ma)") +
    scale_y_continuous(name = "Mean LN ZH") +
    scale_color_manual(values = col.form) #seeingly no correlation

p.ang_g <- ggplot(angle_diff_between_Gs) +
    geom_point(aes(x = angle.diff_Gs.time, y = angle_diff_Gs),
               size = 5, shape = 17) +
    scale_x_discrete(name = "Formation Transition",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(name = "Angle difference between G matrices", 
                       lim = c(0, 90)) + 
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA))

ggsave(p.ang_g, 
       file = "./Results/angle.g.diff.png", 
       width = 20, height = 20, units = "cm")

##### SIZE OF P AND G ACROSS TIME AND ANGLE -----
P_size <- c(sum(diag(P_ext_NKBS)),
            sum(diag(P_ext_NKLS)),
            sum(diag(P_ext_tewk)),
            sum(diag(P_ext_wai)),
            sum(diag(P_ext_uki)),
            sum(diag(P_ext_tai)),
            sum(diag(P_ext_SHCSBSB)))

G_size <- c(sum(diag(G_ext_NKBS)),
            sum(diag(G_ext_NKLS)),
            sum(diag(G_ext_tewk)),
            sum(diag(G_ext_wai)),
            sum(diag(G_ext_uki)),
            sum(diag(G_ext_tai)),
            sum(diag(G_ext_SHCSBSB)))

form.name <- form.df$formationCode
P_G_size <- as.data.frame(cbind(form.name, P_size, G_size))
P_G_size$form.name <- factor(P_G_size$form.name, levels = c("NKLS", "NKBS",
                                                            "Tewkesbury", 
                                                            "Waipuru",
                                                            "Upper Kai-Iwi",
                                                            "Tainui",
                                                            "SHCSBSB"))

#proportion G to P
P_G_size$prop.g.p <- as.numeric(P_G_size$G_size)/as.numeric(P_G_size$P_size)
range(P_G_size$prop.g.p)
#30 to 65 %



ggplot(P_G_size) +
    geom_point(aes(y = as.numeric(P_size), x = form.name),
               size = 5, shape = 17,
               color = "black") + 
    geom_point(aes(y = as.numeric(G_size), x = form.name),
               size = 5,
               color = "black") + 
    theme(text = element_text(size = 16),
          #legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA)) +
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
               size = 5,
               color = col.form) + 
    theme(text = element_text(size = 16),
          #legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA)) +
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

p.temp.pc <- ggplot(df.form.pc) +
    geom_point(aes(x = med.O18, y = as.numeric(PC1_P)),
               col = col.form) + 
    geom_smooth(aes(x = med.O18, y = as.numeric(PC1_P)),
                method = "lm") +
    theme(text = element_text(size = 16)) +
    scale_x_continu ous(expression(mean~delta^18~O)) +
    scale_y_continuous(expression(PC1~of~P~matrix)) + 
    theme(text = element_text(size = 16),
          legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill='transparent', color=NA))

summary(lm(as.numeric(df.form.pc$PC1_P) ~ df.form.pc$med.O18))
#slope = 0.16682; r2 = 0.4162; p-value = 0.07002

ggsave(p.temp.pc, 
       file = "./Results/temp.pc.png", 
       width = 14, height = 10, units = "cm")

#NO PATTERN
summary(lm(df.form.pc$PC1_P ~ df.form.pc$med.O18))
#slope = 0.16682; p = 0.07002, R2 = 0.4162

