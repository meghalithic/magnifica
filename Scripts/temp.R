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

diff.df <- read.csv("./Results/diff.in.traits.csv",
                    header = TRUE,
                    sep = ",",
                    stringsAsFactors = FALSE)
diff.df <- diff.df[1:6,]

#### TEMPERATURE COMPARITONS ----
forms <- c("NKLS", "NKBS", "Tewkesbury",
           "Upper Kai-Iwi", "Tainui",
           "SHCSBSB", "modern")

form.meta.trim <- form.meta[form.meta$formationCode %in% forms,]

##### TEMPERATURE OVER TIME -----
p.deltaO <- ggplot(form.meta.trim) +
    geom_point(aes(formationCode, med.O18),
               shape = 16, size = 5) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    ylab(expression(delta^18~O)) + 
    xlab(label = "Formation") +
    plot.theme
#warmer is more negative

p.temp.form <- ggplot(form.meta.trim) +
    geom_point(aes(formationCode, temp),
               shape = 16, size = 5) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    ylab("Temperature (˚C)") + 
    xlab(label = "Formation") + 
    plot.theme

ggsave(p.temp.form, 
       file = "./Results/temp.form.png", 
       width = 14, height = 10, units = "cm")

form.meta.trim$mean.age <- (as.numeric(form.meta.trim$End_age) + as.numeric(form.meta.trim$Start_age))/2
p.temp.age <- ggplot(form.meta.trim) +
    geom_point(aes(mean.age, temp),
               shape = 16, size = 5) +
    ylab("Temperature (˚C)") + 
    xlab(label = "Age (mya)") + 
    plot.theme

ggsave(p.temp.age, 
       file = "./Results/temp.age.png", 
       width = 14, height = 10, units = "cm")

##### TEMPERATURE AND TRAIT CHANGES OVER TIME -----
mean_by_formation.meta <- merge(mean_by_formation,
                                form.meta,
                                by.x = "formation",
                                by.y = "formationCode")

#raw zh
ggplot(mean_by_formation.meta) +
    geom_point(aes(temp, exp(avg.zh)),
               shape = 16, size = 5) +
    xlab("Temperature (˚C)") + 
    ylab(expression(Zooid~Height~(mu*m))) + 
    plot.theme

#ln zh
summary(lm(mean_by_formation.meta$avg.zh ~ mean_by_formation.meta$temp))
#nonsig p = 0.9514; slope = 0.001625; r2 = 0

#raw zh
summary(lm(exp(mean_by_formation.meta$avg.zh) ~ mean_by_formation.meta$temp))
#nonsig p = 0.9711; slope = 0.8223; r2 = 0

summary(lm(mean_by_formation.meta$avg.mpw.b ~ mean_by_formation.meta$temp))
#nonsig
summary(lm(mean_by_formation.meta$avg.cw.m ~ mean_by_formation.meta$temp))
#nonsig
summary(lm(mean_by_formation.meta$avg.cw.d ~ mean_by_formation.meta$temp))
#nonsig
summary(lm(mean_by_formation.meta$avg.ow.m ~ mean_by_formation.meta$temp))
#nonsig
summary(lm(mean_by_formation.meta$avg.oh ~ mean_by_formation.meta$temp))
#nonsig
summary(lm(mean_by_formation.meta$avg.o.side ~ mean_by_formation.meta$temp))
#nonsig
summary(lm(mean_by_formation.meta$avg.c.side ~ mean_by_formation.meta$temp))
#nonsig

##### AMOUNT OF CHANGE COMPARED TO TEMPERATURE -----
temp.diff <- c((form.meta.trim$temp[form.meta.trim$formationCode == "NKLS"] - form.meta.trim$temp[form.meta.trim$formationCode == "NKBS"]),
               (form.meta.trim$temp[form.meta.trim$formationCode == "NKBS"] - form.meta.trim$temp[form.meta.trim$formationCode == "Tewkesbury"]),
               (form.meta.trim$temp[form.meta.trim$formationCode == "Tewkesbury"] - form.meta.trim$temp[form.meta.trim$formationCode == "Upper Kai-Iwi"]),
               (form.meta.trim$temp[form.meta.trim$formationCode == "Upper Kai-Iwi"] - form.meta.trim$temp[form.meta.trim$formationCode == "Tainui"]),
               (form.meta.trim$temp[form.meta.trim$formationCode == "Tainui"] - form.meta.trim$temp[form.meta.trim$formationCode == "SHCSBSB"]),
               (form.meta.trim$temp[form.meta.trim$formationCode == "SHCSBSB"] - form.meta.trim$temp[form.meta.trim$formationCode == "modern"]))
diff.df$temp.diff <- temp.diff

diff.df$diff.ln.zh <- as.numeric(diff.df$diff.ln.zh) 
diff.df$diff.ln.mpw.b <- as.numeric(diff.df$diff.ln.mpw.b)
diff.df$diff.ln.cw.m <- as.numeric(diff.df$diff.ln.cw.m)
diff.df$diff.ln.cw.d <- as.numeric(diff.df$diff.ln.cw.d)
diff.df$diff.ln.ow.m <- as.numeric(diff.df$diff.ln.ow.m)
diff.df$diff.ln.oh <- as.numeric(diff.df$diff.ln.oh)
diff.df$diff.ln.o.side <- as.numeric(diff.df$diff.ln.o.side)
diff.df$diff.ln.c.side <- as.numeric(diff.df$diff.ln.c.side)

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.zh),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Zooid~Height~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.zh ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.mpw.b),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Median~Process~Width~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.mpw.b ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.cw.m),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Cryptocyst~Width~at~Midline~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.cw.m ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.cw.d),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Cryptocyst~Width~at~Distal~end~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.cw.d ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.ow.m),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Operculum~Width~at~Midline~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.ow.m ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.oh),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Operculum~Height~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.oh ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.o.side),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Operculum~Side~Length~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.o.side ~ diff.df$temp.diff)) #nonsig

ggplot(diff.df) +
    geom_point(aes(temp.diff, diff.ln.c.side),
               shape = 16, size = 5) +
    xlab(expression(Delta~Temperature~(~degree~C))) + 
    ylab(expression(Delta~ln~Cryptocyst~Side~Length~(mu*m))) + 
    plot.theme
summary(lm(diff.df$diff.ln.c.side ~ diff.df$temp.diff)) #nonsig

### COMPARISON TO OTHER BRYOZOA ----
#FOR EVERY 1˚C, A CHANGE IN 0 TO 6.5 um IN LENGTH

# can also read in "temperature_zooid" in Dropbox/ROCKS-PARADOX/Bryozoans

#from Menon 1972
#ep = Electra pilosa
#cr = Conopenum reticulum
#units are um
#degrees in C
ep.temps <- c(6, 12, 18, 22)
ep.len <- c(686.34, 596.37, 586.35, 577)
ep.wid <- c(303.61, 311.81, 314.37, 314.30)
summary(lm(ep.len ~ ep.temps)) #slope: -6.429 
summary(lm(ep.wid ~ ep.temps)) #slope: 0.6621 

plot(ep.len ~ ep.temps)
plot(ep.wid ~ ep.temps)


cr.temps <- c(12, 18, 22)
cr.len <- c(558.31, 518.65, 500)
cr.wid <- c(315.26, 285.18, 314.24)
summary(lm(cr.len ~ cr.temps)) #slope: -5.8925 
summary(lm(cr.wid ~ cr.temps)) #slope: -0.4897 

plot(cr.len ~ cr.temps)
plot(cr.wid ~ cr.temps)


#from Silén & Harmelin
#Haplopoma sciaphium
#units are um
hs.temps <- c()
hs.temps[1] <- mean(c(12, 23))
hs.temps[2] <- mean(c(9, 12))
hs.temps[3] <- mean(c(1, 17))
hs.len <- c(440.2, 501.7, 566.8)
summary(lm(hs.len ~ hs.temps)) #slope = -12.990

plot(hs.len ~ hs.temps)

#Lombardi et al.
#Pentapora fascialis
#mart = mean annual range of temperature; C
#mat = mean annual temperature; C
#units are um
pf.loc <- c("Plymouth", "Lizard", "Grmac", "Tino Island",
            "Palau", "Praiano", "Cala Gonone", "Palmi",
            "Scoglitti")
pf.mart <- c(7.5, 5.5, 3.2, 8.9, 9.2, 10, 9.4, 8, 6.3)
pf.mat <- c(12, 11, 10.5, 17.2, 15.5, 17.4, 17.5, 16.6, 17.5)
pf.len <- c(824.41, 849.69, 753.16, 782.94, 797.14, 757.66, 766,38, 863.33)
pf.wid <- c(458.22, 415.77, 420.89, 409.91, 401.16, 429.22, 420.81, 427.32, 380.57)
summary(lm(pf.mat ~ pf.len)) #slope = -0.002522
summary(lm(pf.mat ~ pf.wid)) #slope = -0.05435
summary(lm(pf.mart ~ pf.len)) #slope = -0.0009694
summary(lm(pf.mart ~ pf.wid)) #slope = 0.009702

plot(pf.len ~ pf.mart)
plot(pf.wid ~ pf.mart)

plot(pf.len ~ pf.mat)
plot(pf.wid ~ pf.mat)

#### HOW MUCH TEMP CHANGE ----
#difference in small zooids to large zooids is 78.73 µm:
#y = mx + b
#y = zh
#figuring out x
#modern: 892.25 (log: 6.775514)
#NKLS: 795.29 (log: 6.681356)
#based on ep:
#starting temp
(795.2878-704.74)/-6.43 #-14.08208
#ending temp
(892.25-704.74)/-6.43 #-29.16174
#diff
-29.16174--14.08208 #-15.07966

#based on cr:
(795.29-627.79)/-5.89 #-28.43803
#ending temp
(892.25-627.79)/-5.89 #-44.89983
#diff
-44.89983--28.43803 #-16.4618

#based on hs:
(795.29-663.11)/-12.99 #-10.17552
#ending temp
(892.25-663.11)/-12.99 #-17.63972
#diff
-17.63972--10.17552 #-7.4642

#based on pf: NONSENSICAL
(795.29-16.82)/-0.003
#ending temp
(892.25-16.82)/-0.003

#based on stegino
(795.29-856.63)/0.82 #-74.80488
#ending temp
(892.25-856.63)/0.82 #43.43902
#diff
43.43902--74.80488 #118.2439


