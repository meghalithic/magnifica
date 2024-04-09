## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) creates trait measurements
# 2) investigates causes of measurement differences

## The output of this code is:
# - trimmed trait df saved as colonies.traits_date.csv
# - summary statistics of traits (mean_by_formation written, mean_by_formation_colony) saved as sum.data.list

#### LOAD DATA ----

source("./Scripts/0-env.R")

df <- read.csv("Results/traits_26Feb2024.csv", #30Nov2023,29Sept2023.csv",
               header = TRUE)
#output from traits.R

#### MANIPULATE DATA ----
traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

df$formation <- factor(df$formation, 
                       levels = c("NKLS", "NKBS", "Tewkesbury",
                                  "Upper Kai-Iwi", "Tainui",
                                  "SHCSBSB", "modern"))

##### DISTRIBUTION -----
p.ln.zh <- ggplot(df) +
    geom_density(aes(x = ln.zh)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

p.ln.zh.form <- ggplot(df) +
    geom_density(aes(x = ln.zh,
                     group = formation,
                     col = formation)) + #lots are bimodal
    ggtitle(paste0("Distribution of traits, N zooids = ", length(unique(df$zooid.id)),
                   ", N colony = ", length(unique(df$colony.id)))) +
    plot.theme +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

traits.melt <- melt(data = df,
                    id.vars = c("boxID", "zooid.id","image",
                                "colony.id", "formation"),
                    variable.name = "measurementType",
                    value.name = "measurementValue")
length(unique(traits.melt$colony.id)) #589 unique colonies [previously 742]

traits.stats <- traits.melt %>%
    dplyr::group_by(measurementType) %>%
    dplyr::summarise(avg = mean(measurementValue))

traits.stats.form <- traits.melt %>%
    dplyr::group_by(measurementType, formation) %>%
    dplyr::summarise(avg = mean(measurementValue))

traits.melt.trim <- traits.melt[traits.melt$measurementType %in% traits,]

p.dist <- ggplot(traits.melt.trim) +
    geom_density(aes(x = measurementValue,
                     group = measurementType,
                     col = measurementType)) + #lots are bimodal
    ggtitle(paste0("Distribution of traits, N zooids = ", length(unique(traits.melt$zooid.id)), ", N colony = ", length(unique(traits.melt$colony.id)))) +
    plot.theme +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~trait~(mu*m)))

##### REMOVE PUTATIVE CRYPTIC SPECIES -----
## bimodality in traits, driven by colonies with unusually small zooids

## Find cutoff for frequencies of size bins
# frequency by size bin (quarter ln bins)
df.bins <- df %>% 
  mutate(zh.bin = cut(ln.zh, breaks = seq(5.5, 7.5, .1))) %>%
  as.data.frame()

df.bin.f <- df.bins %>%
  dplyr::group_by(zh.bin) %>%
  dplyr::summarise(n = n()) %>%
  as.data.frame()
View(df.bin.f)
#(6.2,6.3]
#6.25 like I eyeballed if sort by bin size

sm.traits <- df[df$ln.zh < 6.25,]
sm.colonies <- unique(sm.traits$colony.id)
length(sm.colonies) #41; was 95 images out of 891

bins <- c("(5.5,5.6]", "(5.6,5.7]", "(5.7,5.8]",
          "(5.8,5.9]", "(5.9,6]", "(6,6.1]", 
          "(6.1,6.2]", "(6.2,6.3]")

df.bins$zh.bin <- as.character(df.bins$zh.bin)

df.bins$sm <- FALSE
df.bins$sm[df.bins$zh.bin %in% bins] <- TRUE

#look at proportions
prop.sm <- df.bins %>% 
  dplyr::group_by(colony.id) %>%
  dplyr::summarise(n.zooid = length(zooid.id),
                   n.sm.zooid = sum(sm),
                   prop.sm = n.sm.zooid/n.zooid)
View(prop.sm)

#goes from 100% to 20%; perhaps make 20% the cut off
rm.col <- prop.sm$colony.id[prop.sm$prop.sm == 1]
length(rm.col) #31 colonies

reg.colonies <- df[!(df$colony.id %in% rm.col),]
nrow(reg.colonies) #5687
length(unique(reg.colonies$colony.id)) #557

##### WRITE OUT DATASET ----
write.csv(reg.colonies,
          "./Results/colonies.traits_26Feb2024.csv",
          row.names = FALSE)


##### ABOUT TRAITS -----
df <- reg.colonies
str(df$formation)

mean_by_formation_colony = df %>% #use this going forward
    dplyr::group_by(formation, colony.id) %>%
    dplyr::summarize(n.zooid = length(zooid.id),
                     
                     avg.zh = mean(ln.zh, na.rm = T),
                     sd.zh = sd(ln.zh, na.rm = T),
                     
                     avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                     sd.mpw.b = sd(ln.mpw.b, na.rm = T),
                     
                     avg.cw.m = mean(ln.cw.m, na.rm = T),
                     sd.cw.m = sd(ln.cw.m, na.rm = T),
                     
                     avg.cw.d = mean(ln.cw.d, na.rm = T),
                     sd.cw.d = sd(ln.cw.d, na.rm = T),
                     
                     avg.ow.m = mean(ln.ow.m, na.rm = T),
                     sd.ow.m = sd(ln.ow.m, na.rm = T),
                     
                     avg.oh = mean(ln.oh, na.rm = T),
                     sd.oh = sd(ln.oh, na.rm = T),
                     
                     avg.o.side = mean(ln.o.side, na.rm = T),
                     sd.o.side = sd(ln.o.side, na.rm = T),
                     
                     avg.c.side = mean(ln.c.side, na.rm = T),
                     sd.c.side = sd(ln.c.side, na.rm = T)) %>%
    as.data.frame()
min(mean_by_formation_colony$n.zooid) #5

mean_by_formation = df %>%
    dplyr::group_by(formation) %>%
    dplyr::summarize(num.col = length(unique(colony.id)),
                     num.zooid = length(unique(zooid.id)),
                     avg.zooid = ceiling(num.zooid/num.col), #round up to nearest integer
                     
                     min.zh = min(ln.zh, na.rm = T),
                     max.zh = max(ln.zh, na.rm = T),
                     avg.zh = mean(ln.zh, na.rm = T),
                     sd.zh = sd(ln.zh, na.rm = T),
                     var.zh = var(ln.zh, na.rm = T),
                     
                     min.mpw.b = min(ln.mpw.b, na.rm = T),
                     max.mpw.b = max(ln.mpw.b, na.rm = T),
                     avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                     sd.mpw.b = sd(ln.mpw.b, na.rm = T),
                     var.mpw.b = var(ln.mpw.b, na.rm = T),
                     
                     min.cw.m = min(ln.cw.m, na.rm = T),
                     max.cw.m = max(ln.cw.m, na.rm = T),
                     avg.cw.m = mean(ln.cw.m, na.rm = T),
                     sd.cw.m = sd(ln.cw.m, na.rm = T),
                     var.cw.m = var(ln.cw.m, na.rm = T),
                     
                     min.cw.d = min(ln.cw.d, na.rm = T),
                     max.cw.d = max(ln.cw.d, na.rm = T),
                     avg.cw.d = mean(ln.cw.d, na.rm = T),
                     sd.cw.d = sd(ln.cw.d, na.rm = T),
                     var.cw.d = var(ln.cw.d, na.rm = T),
                     
                     min.ow.m = min(ln.ow.m, na.rm = T),
                     max.ow.m = max(ln.ow.m, na.rm = T),
                     avg.ow.m = mean(ln.ow.m, na.rm = T),
                     sd.ow.m = sd(ln.ow.m, na.rm = T),
                     var.ow.m = var(ln.ow.m, na.rm = T),
                     
                     min.oh = min(ln.oh, na.rm = T),
                     max.oh = max(ln.oh, na.rm = T),
                     avg.oh = mean(ln.oh, na.rm = T),
                     sd.oh = sd(ln.oh, na.rm = T),
                     var.oh = var(ln.oh, na.rm = T),
                     
                     min.o.side = min(ln.o.side, na.rm = T),
                     max.o.side = max(ln.o.side, na.rm = T),
                     avg.o.side = mean(ln.o.side, na.rm = T),
                     sd.o.side = sd(ln.o.side, na.rm = T),
                     var.o.side = var(ln.o.side, na.rm = T),
                     
                     min.c.side = min(ln.c.side, na.rm = T),
                     max.c.side = max(ln.c.side, na.rm = T),
                     avg.c.side = mean(ln.c.side, na.rm = T),
                     sd.c.side = sd(ln.c.side, na.rm = T),
                     var.c.side = var(ln.c.side, na.rm = T)) %>%
    as.data.frame()
## Grabowski & Porto claim sampling of 60 per sp...
#NKBS has heighest variance in traits because high sample size

colony_means = df %>%
    dplyr::group_by(colony.id) %>%
    dplyr::summarize(formation = formation[1],
                     n.zooid = length(unique(zooid.id)),
                     avg.zh = mean(ln.zh, na.rm = T),
                     avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                     avg.cw.m = mean(ln.cw.m, na.rm = T),
                     avg.cw.d = mean(ln.cw.d, na.rm = T),
                     avg.ow.m = mean(ln.ow.m, na.rm = T),
                     avg.oh = mean(ln.oh, na.rm = T),
                     avg.o.side = mean(ln.o.side, na.rm = T),
                     avg.c.side = mean(ln.c.side, na.rm = T)) %>%
    as.data.frame()

means = df %>%
    dplyr::summarize(avg.zh = mean(ln.zh, na.rm = T),
                     avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                     avg.cw.m = mean(ln.cw.m, na.rm = T),
                     avg.cw.d = mean(ln.cw.d, na.rm = T),
                     avg.ow.m = mean(ln.ow.m, na.rm = T),
                     avg.oh = mean(ln.oh, na.rm = T),
                     avg.o.side = mean(ln.o.side, na.rm = T),
                     avg.c.side = mean(ln.c.side, na.rm = T)) %>%
    as.data.frame()

sum.data.list = list(mean_by_formation, mean_by_formation_colony, means)
save(sum.data.list,
     file = "./Results/sum.data.list.w.modern.RData")

#### CORRELATIONS & ALLOMETRIES ----
## are these coming from the same individuals??
## ask KLV for other metadata for sites
## OH hump is on the other side

#from https://r-coder.com/correlation-plot-r/
panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    his <- hist(x, plot = FALSE)
    breaks <- his$breaks
    nB <- length(breaks)
    y <- his$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = rgb(0, 1, 1, alpha = 0.5), ...)
    # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- abs(cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

plot.df <- df %>%
    dplyr::select(ln.zh, ln.mpw.b, ln.cw.m, ln.cw.d, ln.ow.m, ln.oh, ln.o.side, ln.c.side)

pairs(plot.df,
      upper.panel = panel.cor,         # Disabling the upper panel
      diag.panel = panel.hist)

columns <- c("trait.y", "trait.x",
             "p-value", "slope", "r2", "df")
traitCorr <- data.frame(matrix(nrow = 1, ncol = length(columns)))
colnames(traitCorr) <- columns

trait.corr <- colnames(plot.df)
for(i in 1:length(plot.df)){
    for(j in 1:length(trait.corr)){
        model <- lm(plot.df[, trait.corr[j]] ~ plot.df[,i], 
                    na.action = na.exclude)
        sum.model <- summary(model)
        sub <- c(trait.corr[j],
                 colnames(plot.df[,i]),
                 sum.model$coefficients[8],
                 model$coefficients[[2]],
                 sum.model$adj.r.squared,
                 sum.model$df[2])
        
        traitCorr <- rbind(traitCorr, sub)
    }
}

write.csv(traitCorr,
          "./Results/trait.correlations.w.modern.csv",
          row.names = FALSE)

#### LOOK AT CHANGES OVER TIME AND BETWEEN FORMATIONS ------
## Add meta data
form.df <- form.meta[c(1:7,12),] #in same order as mean_by_formation

for(i in 1:nrow(form.df)){
  form.df$mean.age[i] <- mean(form.df$Start_age[i], form.df$End_age[i], na.rm = TRUE)
}

form.df$age.range <- ""
for(i in 1:nrow(form.df)){
  form.df$age.range[i] <- form.df$Start_age[i] - form.df$End_age[i]
}
form.df$age.range <- as.numeric(form.df$age.range)

mean_by_formation.meta <- merge(mean_by_formation, form.df,
                                by.x = "formation",
                                by.y = "formationCode")

## overall difference in zooid height
mean_by_formation$avg.zh[mean_by_formation$formation == "NKLS"] - mean_by_formation$avg.zh[mean_by_formation$formation == "modern"]
#-0.09371967 (decrease in length)
mean_by_formation$avg.ow.m[mean_by_formation$formation == "NKLS"] - mean_by_formation$avg.ow.m[mean_by_formation$formation == "modern"]
#-0.1627765 (decrease in width)

diff.ln.zh <- c()
diff.ln.mpw.b <- c()
diff.ln.cw.m <- c()
diff.ln.cw.d <- c()
diff.ln.ow.m <- c()
diff.ln.oh <- c()
diff.ln.o.side <- c()
diff.ln.c.side <- c()

diff.zh <- c()
diff.mpw.b <- c()
diff.cw.m <- c()
diff.cw.d <- c()
diff.ow.m <- c()
diff.oh <- c()
diff.o.side <- c()
diff.c.side <- c()

diff.zh.per <- c()
diff.mpw.b.per <- c()
diff.cw.m.per <- c()
diff.cw.d.per <- c()
diff.ow.m.per <- c()
diff.oh.per <- c()
diff.o.side.per <- c()
diff.c.side.per <- c()

for(i in 1:(nrow(mean_by_formation)-1)){
    diff.ln.zh[i] <- mean_by_formation$avg.zh[i+1] - mean_by_formation$avg.zh[i]
    diff.ln.mpw.b[i] <- mean_by_formation$avg.mpw.b[i+1] - mean_by_formation$avg.mpw.b[i]
    diff.ln.cw.m[i] <- mean_by_formation$avg.cw.m[i+1] - mean_by_formation$avg.cw.m[i]
    diff.ln.cw.d[i] <- mean_by_formation$avg.cw.d[i+1] - mean_by_formation$avg.cw.d[i]
    diff.ln.ow.m[i] <- mean_by_formation$avg.ow.m[i+1] - mean_by_formation$avg.ow.m[i]
    diff.ln.oh[i] <- mean_by_formation$avg.oh[i+1] - mean_by_formation$avg.oh[i]
    diff.ln.o.side[i] <- mean_by_formation$avg.o.side[i+1] - mean_by_formation$avg.o.side[i]
    diff.ln.c.side[i] <- mean_by_formation$avg.c.side[i+1] - mean_by_formation$avg.c.side[i]
    
    diff.zh[i] <- exp(mean_by_formation$avg.zh[i+1]) - exp(mean_by_formation$avg.zh[i])
    diff.mpw.b[i] <- exp(mean_by_formation$avg.mpw.b[i+1]) - exp(mean_by_formation$avg.mpw.b[i])
    diff.cw.m[i] <- exp(mean_by_formation$avg.cw.m[i+1]) - exp(mean_by_formation$avg.cw.m[i])
    diff.cw.d[i] <- exp(mean_by_formation$avg.cw.d[i+1]) - exp(mean_by_formation$avg.cw.d[i])
    diff.ow.m[i] <- exp(mean_by_formation$avg.ow.m[i+1]) - exp(mean_by_formation$avg.ow.m[i])
    diff.oh[i] <- exp(mean_by_formation$avg.oh[i+1]) - exp(mean_by_formation$avg.oh[i])
    diff.o.side[i] <- exp(mean_by_formation$avg.o.side[i+1]) - exp(mean_by_formation$avg.o.side[i])
    diff.c.side[i] <- exp(mean_by_formation$avg.c.side[i+1]) - exp(mean_by_formation$avg.c.side[i])
    
    diff.zh.per[i] <- (abs(exp(mean_by_formation$avg.zh[i+1]) - exp(mean_by_formation$avg.zh[i]))/exp(mean_by_formation$avg.zh[i]))*100
    diff.mpw.b.per[i] <- (abs(exp(mean_by_formation$avg.mpw.b[i+1]) - exp(mean_by_formation$avg.mpw.b[i]))/exp(mean_by_formation$avg.mpw.b[i]))*100
    diff.cw.m.per[i] <- (abs(exp(mean_by_formation$avg.cw.m[i+1]) - exp(mean_by_formation$avg.cw.m[i]))/exp(mean_by_formation$avg.cw.m[i]))*100
    diff.cw.d.per[i] <- (abs(exp(mean_by_formation$avg.cw.d[i+1]) - exp(mean_by_formation$avg.cw.d[i]))/exp(mean_by_formation$avg.cw.d[i]))*100
    diff.ow.m.per[i] <- (abs(exp(mean_by_formation$avg.ow.m[i+1]) - exp(mean_by_formation$avg.ow.m[i]))/exp(mean_by_formation$avg.ow.m[i]))*100
    diff.oh.per[i] <- (abs(exp(mean_by_formation$avg.oh[i+1]) - exp(mean_by_formation$avg.oh[i]))/exp(mean_by_formation$avg.oh[i]))*100
    diff.o.side.per[i] <- (abs(exp(mean_by_formation$avg.o.side[i+1]) - exp(mean_by_formation$avg.o.side[i]))/exp(mean_by_formation$avg.o.side[i]))*100
    diff.c.side.per[i] <- (abs(exp(mean_by_formation$avg.c.side[i+1]) - exp(mean_by_formation$avg.c.side[i]))/exp(mean_by_formation$avg.c.side[i]))*100
}


exp(mean_by_formation$avg.zh[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.zh[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.mpw.b[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.mpw.b[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.cw.m[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.cw.m[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.cw.d[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.cw.d[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.ow.m[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.ow.m[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.oh[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.oh[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.o.side[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.o.side[mean_by_formation$formation == "NKLS"])
exp(mean_by_formation$avg.c.side[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.c.side[mean_by_formation$formation == "NKLS"])

(abs(exp(mean_by_formation$avg.zh[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.zh[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.zh[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.mpw.b[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.mpw.b[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.mpw.b[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.cw.m[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.cw.m[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.cw.m[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.cw.d[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.cw.d[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.cw.d[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.ow.m[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.ow.m[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.ow.m[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.oh[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.oh[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.oh[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.o.side[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.o.side[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.o.side[mean_by_formation$formation == "NKLS"]))*100
(abs(exp(mean_by_formation$avg.c.side[mean_by_formation$formation == "modern"]) - exp(mean_by_formation$avg.c.side[mean_by_formation$formation == "NKLS"]))/exp(mean_by_formation$avg.c.side[mean_by_formation$formation == "NKLS"]))*100



diff.form <- c("NKLS to NKBS", "NKBS to Tewkesbury",
               "Tewkesbury to Upper Kai-Iwi", "Upper Kai-Iwi to Tainui",
               "Tainui to SHCSBSB", "SHCSBSB to modern")
diff.stats <- as.data.frame(cbind(diff.form, diff.ln.zh, diff.ln.mpw.b, 
                                  diff.ln.cw.m, diff.ln.cw.d, diff.ln.ow.m, 
                                  diff.ln.oh, diff.ln.o.side, diff.ln.c.side,
                                  
                                  diff.zh, diff.mpw.b, 
                                  diff.cw.m, diff.cw.d, diff.ow.m, 
                                  diff.oh, diff.o.side, diff.c.side,
                                  
                                  diff.zh.per, diff.mpw.b.per, 
                                  diff.cw.m.per, diff.cw.d.per, diff.ow.m.per, 
                                  diff.oh.per, diff.o.side.per, diff.c.side.per))

write.csv(diff.stats,
          "./Results/diff.in.traits.no.wai.csv",
          row.names = FALSE)

## how is sd a function of sample size (number of zooids and number of colonies)?
#plot sd per colony by zooid no
ggplot(mean_by_formation_colony) + 
  geom_point(aes(x = n.zooid, y = sd.zh,
                 col = formation)) + 
  plot.theme +
  scale_x_continuous(name = "Number of Zooids per Colony") +
  scale_y_continuous(expression(sd~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

summary(lm(mean_by_formation_colony$sd.zh ~ mean_by_formation_colony$n.zooid))
#sig but slope close to 0
summary(lm(mean_by_formation_colony$sd.mpw.b ~ mean_by_formation_colony$n.zooid))
#nonsig and slope close 0
summary(lm(mean_by_formation_colony$sd.cw.m ~ mean_by_formation_colony$n.zooid))
#nonsig and slop close to 0
summary(lm(mean_by_formation_colony$sd.cw.d ~ mean_by_formation_colony$n.zooid))
#nonsig and slope close to 0
summary(lm(mean_by_formation_colony$sd.ow.m ~ mean_by_formation_colony$n.zooid))
#nonsig and slope close to 0
summary(lm(mean_by_formation_colony$sd.oh ~ mean_by_formation_colony$n.zooid))
#nonsig and slope around 0
summary(lm(mean_by_formation_colony$sd.o.side ~ mean_by_formation_colony$n.zooid))
#nonsig and slope close to 0 
summary(lm(mean_by_formation_colony$sd.c.side ~ mean_by_formation_colony$n.zooid))
#sig and slope around 0

#plot sd per formation by colony no
ggplot(mean_by_formation) + 
  geom_point(aes(x = num.col, y = sd.zh,
                 col = formation)) + 
  plot.theme +
  scale_x_continuous(name = "Number of Colonies per Colony") +
  scale_y_continuous(expression(sd~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

anova(lm(mean_by_formation$sd.zh ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.mpw.b ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.cw.m ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.cw.d ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.ow.m ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.oh ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.o.side ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig
anova(lm(mean_by_formation$sd.c.side ~ mean_by_formation$num.col + mean_by_formation$num.zooid))
#nonsig

#use formation means
#mean_by_formation
ggplot(data = mean_by_formation.meta) +
  geom_point(aes(x = as.numeric(age.range), y = var.zh,
                 col = formation)) + 
  plot.theme +
  scale_x_continuous(name = "Age Range (Ma)") +
  scale_y_continuous(name = "Variation of Zooid Height (um)") +
  scale_color_manual(values = col.form)

summary(lm(mean_by_formation.meta$var.zh ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.mpw.b ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.cw.m ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.cw.d ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.ow.m ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.oh ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.o.side ~ mean_by_formation.meta$age.range))
#nonsig
summary(lm(mean_by_formation.meta$var.c.side ~ mean_by_formation.meta$age.range))
#nonsig

## figure goes with mean_by_formation table
df.form.meta <- merge(df, form.df,
                      by.x = "formation",
                      by.y = "formationCode")
          
traits.melt$formation <- factor(traits.melt$formation,
                                levels = c("NKLS", "NKBS", "Tewkesbury", 
                                           "Upper Kai-Iwi", "Tainui", "SHCSBSB", "modern")) 

box.ln.zh <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.zh",], 
       aes(x = formation, 
           y = measurementValue, 
           fill = formation)) +
  geom_boxplot() +
  scale_color_manual(values = col.form) +
  scale_fill_manual(values = col.form) +
  scale_x_discrete(name = "Formation",
                   guide = guide_axis(angle = 45)) +
  scale_y_continuous(expression(Zooid~Height~(mu*m)),
                     limits = c(5.5, 7.75),
                     breaks = c(log(250), log(400),
                                log(650), log(1000),
                                log(1450), log(2000)),
                     labels = c(250, 400,
                                650, 1000,
                                1450, 2000)) +
  plot.theme

ggsave(box.ln.zh, 
       file = "./Results/boxplot.ln.zh.w.modern.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.oh <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.oh",], 
                    aes(x = formation, 
                        y = measurementValue, 
                        fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Operculum~Height~(mu*m)),
                       limits = c(-3, -.5),
                       breaks = c(log(0.05), log(0.1),
                                  log(.15), log(.25),
                                  log(.5)),
                       labels = c(0.05, 0.1,
                                  0.15, 0.25,
                                  0.5)) +
    plot.theme

ggsave(box.ln.oh, 
       file = "./Results/boxplot.ln.oh.w.modern.new.y.axis.png", 
       width = 14, height = 10, units = "cm")


box.ln.mpw.b <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.mpw.b",], 
                    aes(x = formation, 
                        y = measurementValue, 
                        fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Median~Process~Width~(mu*m)),
                       limits = c(3.5, 6),
                       breaks = c(log(50), log(100),
                                  log(150), log(250),
                                  log(400)),
                       labels = c(50, 100,
                                  150, 250,
                                  400)) +
    plot.theme

ggsave(box.ln.mpw.b, 
       file = "./Results/boxplot.ln.mpw.b.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.mpw.b.simple <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.mpw.b",], 
       aes(x = formation, 
           y = measurementValue, 
           fill = formation)) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Median~Process~Width~(mu*m)),
                       limits = c(4.5, 6),
                       breaks = c(log(100),
                                  log(150), log(250),
                                  log(400)),
                       labels = c(100,
                                  150, 250,
                                  400)) +
    plot.theme
ggsave(box.mpw.b.simple, 
       file = "./Results/simple.boxplot.ln.mpw.b.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.cw.m <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.cw.m",], 
                       aes(x = formation, 
                           y = measurementValue, 
                           fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Cryptocyst~Width~at~Midline~(mu*m)),
                       limits = c(4, 7),
                       breaks = c(log(100), log(150),
                                  log(250), log(400),
                                  log(600), log(1000)),
                       labels = c(100, 150,
                                  250, 400,
                                  600, 1000)) +
    plot.theme

ggsave(box.ln.cw.m, 
       file = "./Results/boxplot.ln.cw.m.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.cw.m.simple <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.cw.m",], 
                           aes(x = formation, 
                               y = measurementValue, 
                               fill = formation)) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Cryptocyst~Width~at~Midline~(mu*m)),
                       limits = c(5, 7),
                       breaks = c(log(150),
                                  log(250), log(400),
                                  log(600), log(1000)),
                       labels = c(150,
                                  250, 400,
                                  600, 1000)) +
    plot.theme
ggsave(box.cw.m.simple, 
       file = "./Results/simple.boxplot.ln.cw.m.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.cw.d <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.cw.d",], 
                      aes(x = formation, 
                          y = measurementValue, 
                          fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Cryptocyst~Width~at~Distal~end~(mu*m)),
                       limits = c(4, 7),
                       breaks = c(log(100), log(150),
                                  log(250), log(400),
                                  log(600), log(1000)),
                       labels = c(100, 150,
                                  250, 400,
                                  600, 1000)) +
    plot.theme

ggsave(box.ln.cw.d, 
       file = "./Results/boxplot.ln.cw.d.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.ow.m <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.ow.m",], 
                      aes(x = formation, 
                          y = measurementValue, 
                          fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Operculum~Width~(mu*m)),
                       limits = c(4, 7),
                       breaks = c(log(100), log(150),
                                  log(250), log(400),
                                  log(600), log(1000)),
                       labels = c(100, 150,
                                  250, 400,
                                  600, 1000)) +
    plot.theme

ggsave(box.ln.ow.m, 
       file = "./Results/boxplot.ln.ow.m.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ow.m.simple <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.ow.m",], 
                          aes(x = formation, 
                              y = measurementValue, 
                              fill = formation)) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Operculum~Width~(mu*m)),
                       limits = c(5.5, 6.5),
                       breaks = c(log(250), log(300), 
                                  log(400),
                                  log(650)),
                       labels = c(250, 300,
                                  400, 650)) +
    plot.theme
ggsave(box.ow.m.simple, 
       file = "./Results/simple.boxplot.ln.ow.m.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.oh <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.oh",], 
                      aes(x = formation, 
                          y = measurementValue, 
                          fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Operculum~Height~(mu*m)),
                       limits = c(-3, -.75),
                       breaks = c(log(.05), log(.1),
                                  log(.2), log(.5)),
                       labels = c(0.05, 0.1,
                                  0.2, 0.5)) + 
    plot.theme

ggsave(box.ln.oh, 
       file = "./Results/boxplot.ln.oh.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.o.side <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.o.side",], 
                    aes(x = formation, 
                        y = measurementValue, 
                        fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Operculum~side~Length~(mu*m)),
                       limits = c(4, 7),
                       breaks = c(log(100), log(150),
                                  log(250), log(400),
                                  log(600), log(1000)),
                       labels = c(100, 150,
                                  250, 400,
                                  600, 1000)) +
    plot.theme

ggsave(box.ln.o.side, 
       file = "./Results/boxplot.ln.o.side.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.ln.c.side <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.c.side",], 
                        aes(x = formation, 
                            y = measurementValue, 
                            fill = formation)) +
    geom_boxplot() +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Cryptocyst~side~Length~(mu*m)),
                       limits = c(4, 7),
                       breaks = c(log(100), log(150),
                                  log(250), log(400),
                                  log(600), log(1000)),
                       labels = c(100, 150,
                                  250, 400,
                                  600, 1000)) +
    plot.theme

ggsave(box.ln.c.side, 
       file = "./Results/boxplot.ln.c.side.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")

box.c.side.simple <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.c.side",], 
                          aes(x = formation, 
                              y = measurementValue, 
                              fill = formation)) +
    geom_boxplot(outlier.shape = NA) +
    scale_color_manual(values = col.form) +
    scale_fill_manual(values = col.form) +
    scale_x_discrete(name = "Formation",
                     guide = guide_axis(angle = 45)) +
    scale_y_continuous(expression(Cryptocyst~side~Length~(mu*m)),
                       limits = c(5, 6.5),
                       breaks = c(log(150),
                                  log(250), log(400),
                                  log(650)),
                       labels = c(150,
                                  250, 400,
                                  650)) +
    plot.theme
ggsave(box.c.side.simple, 
       file = "./Results/simple.boxplot.ln.c.side.w.modern.no.wai.new.y.axis.png", 
       width = 14, height = 10, units = "cm")
