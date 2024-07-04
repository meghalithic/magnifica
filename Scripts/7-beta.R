# Meghan A. Balk
# meghan.balk@gmail.com
# initially created: Jun 2023
# last updated: 29 Aug 2023

## The purpose of this script is to:
# 1. calculate Beta
# 2. test directions of phenotypic change compared to Beta
# 3. test if Beta shapes G over time

## the outputs are:
# - comparison of Beta to Gmax at t1, t2, and between Betas at t1 and t2 
# (tables and graphs)

#### LOAD DATA ----

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

corr_b.t1_gmax.t1 <- c(dot_product.beta_t1_Gmax_t1, dot_product.beta_t2_Gmax_t2,
                       dot_product.beta_t3_Gmax_t3, dot_product.beta_t4_Gmax_t4,
                       dot_product.beta_t5_Gmax_t5, dot_product.beta_t6_Gmax_t6)
ang_b.t1_gmax.t1 <- c(angle_degrees.beta_t1_Gmax_t1, angle_degrees.beta_t2_Gmax_t2, 
                      angle_degrees.beta_t3_Gmax_t3, angle_degrees.beta_t4_Gmax_t4, 
                      angle_degrees.beta_t5_Gmax_t5, angle_degrees.beta_t6_Gmax_t6)
b.t1_gmax.t1 <- as.data.frame(cbind(corr_b.t1_gmax.t1, ang_b.t1_gmax.t1))
b.t1_gmax.t1$form <- formation_transition[1:6]

write.csv(b.t1_gmax.t1,
          "Results/b.t1_gmax.t1.csv",
          row.names = FALSE)

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

corr_b.t1_b.t2 <- c(dot_product.beta_t1_beta_t2, dot_product.beta_t2_beta_t3,
                    dot_product.beta_t3_beta_t4, dot_product.beta_t4_beta_t5,
                    dot_product.beta_t5_beta_t6)
ang_b.t1_b.t2 <- c(angle_degrees.beta_t1_beta_t2, angle_degrees.beta_t2_beta_t3, 
                   angle_degrees.beta_t3_beta_t4, angle_degrees.beta_t4_beta_t5, 
                   angle_degrees.beta_t5_beta_t6)
b.t1_b.t2 <- as.data.frame(cbind(corr_b.t1_b.t2, ang_b.t1_b.t2))
b.t1_b.t2$form_trans <- formation_transition[c(-6, -7)]

write.csv(b.t1_b.t2,
          "Results/b.t1_b.t2.csv",
          row.names = FALSE)

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

corr_b.t1_gmax.t2 <- c(dot_product.beta_t1_Gmax_t2, dot_product.beta_t2_Gmax_t3,
                       dot_product.beta_t3_Gmax_t4, dot_product.beta_t4_Gmax_t5,
                       dot_product.beta_t5_Gmax_t6, dot_product.beta_t6_Gmax_t7)
ang_b.t1_gmax.t2 <- c(angle_degrees.beta_t1_Gmax_t2, angle_degrees.beta_t2_Gmax_t3, 
                      angle_degrees.beta_t3_Gmax_t4, angle_degrees.beta_t4_Gmax_t5, 
                      angle_degrees.beta_t5_Gmax_t6, angle_degrees.beta_t6_Gmax_t7)
b.t1_gmax.t2 <- as.data.frame(cbind(corr_b.t1_gmax.t2, ang_b.t1_gmax.t2))
b.t1_gmax.t2$form_trans <- formation_transition[-7]

write.csv(b.t1_gmax.t2,
          "Results/b.t1_gmax.t2.csv",
          row.names = FALSE)

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
           data = diff.beta.gmax.df[1:6,])) 
#nonsig 0.2847; slope at 0, r2 at 0.09

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
           data = diff.beta.gmax.df[1:6,])) 
#nonsig at 0.08412; slope barely over 0, r2 0.4584

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
           data = diff.beta.gmax.df[1:6,])) 
#nonsig at p = 0.2003; no relationship