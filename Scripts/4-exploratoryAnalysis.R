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

df <- read.csv("Results/traits_8Dec2023.csv", #30Nov2023,29Sept2023.csv",
               header = TRUE)
#output from traits.R

form.meta <- read.csv("~/Documents/GitHub/bryozoa/stegino_metadata/newMetadata/formations.csv", header = TRUE)

#this is downloaded from: http://www.lorraine-lisiecki.com/LR04_MISboundaries.txt
oxy.18 <- read.csv("Data/∂18O.csv",
                   header = TRUE)

#this is from Emanuela and is in the stegino_metadata repository
locality.df <- read.csv("Data/All.NZ.Samples_EDM_31.07.2023_sheet1.csv",
                        header = TRUE)

#### MANIPULATE DATA ----
traits = names(df[, c("ln.zh", "ln.mpw.b", "ln.cw.m", "ln.cw.d", 
                      "ln.ow.m", "ln.oh", "ln.c.side", "ln.o.side")])

df$formation <- factor(df$formation, 
                       levels = c("NKLS", "NKBS", "Tewkesbury",
                                  "Waipuru", "Upper Kai-Iwi",
                                  "Tainui", "SHCSBSB", "modern"))

##### DISTRIBUTION -----
p.ln.zh <- ggplot(df) +
    geom_density(aes(x = ln.zh)) +
    ggtitle(paste0("Zooid height, N zooids = ", nrow(df), ", N colony = ", length(unique(df$colony.id)))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(expression(ln~Zooid~Height~(mu*m)))

p.ln.zh.form <- ggplot(df) +
    geom_density(aes(x = ln.zh,
                     group = formation,
                     col = formation)) + #lots are bimodal
    ggtitle(paste0("Distribution of traits, N zooids = ", length(unique(df$zooid.id)),
                   ", N colony = ", length(unique(df$colony.id)))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
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
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
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
#6.25 like I eyeballed

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
          "./Results/colonies.traits_8Dec2023.csv",
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
                     
                     min.mpw.b = min(ln.mpw.b, na.rm = T),
                     max.mpw.b = max(ln.mpw.b, na.rm = T),
                     avg.mpw.b = mean(ln.mpw.b, na.rm = T),
                     sd.mpw.b = sd(ln.mpw.b, na.rm = T),
                     
                     min.cw.m = min(ln.cw.m, na.rm = T),
                     max.cw.m = max(ln.cw.m, na.rm = T),
                     avg.cw.m = mean(ln.cw.m, na.rm = T),
                     sd.cw.m = sd(ln.cw.m, na.rm = T),
                     
                     min.cw.d = min(ln.cw.d, na.rm = T),
                     max.cw.d = max(ln.cw.d, na.rm = T),
                     avg.cw.d = mean(ln.cw.d, na.rm = T),
                     sd.cw.d = sd(ln.cw.d, na.rm = T),
                     
                     min.ow.m = min(ln.ow.m, na.rm = T),
                     max.ow.m = max(ln.ow.m, na.rm = T),
                     avg.ow.m = mean(ln.ow.m, na.rm = T),
                     sd.ow.m = sd(ln.ow.m, na.rm = T),
                     
                     min.oh = min(ln.oh, na.rm = T),
                     max.oh = max(ln.oh, na.rm = T),
                     avg.oh = mean(ln.oh, na.rm = T),
                     sd.oh = sd(ln.oh, na.rm = T),
                     
                     min.o.side = min(ln.o.side, na.rm = T),
                     max.o.side = max(ln.o.side, na.rm = T),
                     avg.o.side = mean(ln.o.side, na.rm = T),
                     sd.o.side = sd(ln.o.side, na.rm = T),
                     
                     min.c.side = min(ln.c.side, na.rm = T),
                     max.c.side = max(ln.c.side, na.rm = T),
                     avg.c.side = mean(ln.c.side, na.rm = T),
                     sd.c.side = sd(ln.c.side, na.rm = T)) %>%
    as.data.frame()
## Grabowski & Porto claim sampling of 60 per sp...

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

sum.data.list = list(mean_by_formation, mean_by_formation_colony)
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

##between formations
mean_by_formation$avg.zh[mean_by_formation$formation == "NKLS"] - mean_by_formation$avg.zh[mean_by_formation$formation == "NKBS"]
#0.03156643 (0.9298943 diff)
mean_by_formation$avg.zh[mean_by_formation$formation == "NKBS"] - mean_by_formation$avg.zh[mean_by_formation$formation == "Tewkesbury"]
#-0.002326976 (0.9946563 diff)
mean_by_formation$avg.zh[mean_by_formation$formation == "Tewkesbury"] - mean_by_formation$avg.zh[mean_by_formation$formation == "Waipuru"]
#-0.04762085 (0.8961468 diff)
mean_by_formation$avg.zh[mean_by_formation$formation == "Waipuru"] - mean_by_formation$avg.zh[mean_by_formation$formation == "Upper Kai-Iwi"]
#-0.1132389 (0.7704795 diff)
mean_by_formation$avg.zh[mean_by_formation$formation == "Upper Kai-Iwi"] - mean_by_formation$avg.zh[mean_by_formation$formation == "Tainui"]
#-0.03510717 (0.9223438 diff)
mean_by_formation$avg.zh[mean_by_formation$formation == "Tainui"] - mean_by_formation$avg.zh[mean_by_formation$formation == "SHCSBSB"]
#0.01531607 (1.015434 diff); THIS IS THE SAME
mean_by_formation$avg.zh[mean_by_formation$formation == "SHCSBSB"] - mean_by_formation$avg.zh[mean_by_formation$formation == "modern"]
#0.05769168 (1.142067 diff)

## how is sd a function of sample size (number of zooids and number of colonies)?
#plot sd per colony by zooid no
ggplot(mean_by_formation_colony) + 
  geom_point(aes(x = n.zooid, y = sd.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Number of Zooids per Colony") +
  scale_y_continuous(expression(sd~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

#plot sd per formation by colony no
ggplot(mean_by_formation) + 
  geom_point(aes(x = num.col, y = sd.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Number of Colonies per Colony") +
  scale_y_continuous(expression(sd~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

anova(lm(mean_by_formation$sd.zh ~ mean_by_formation$num.col + mean_by_formation$num.zooid))

#use formation means
#mean_by_formation

ggplot(data = mean_by_formation.meta) +
  geom_point(aes(x = as.numeric(age.range), y = avg.zh,
                 col = formation)) + 
  theme(text = element_text(size = 16),
        legend.position = "none") +
  scale_x_continuous(name = "Age Range (Ma)") +
  scale_y_continuous(name = "Average Zooid Height (um)") +
  scale_color_manual(values = col.form)

mean_by_formation.meta$formation <- factor(mean_by_formation.meta$formation, 
                                           levels = c("NKLS", "NKBS", "Tewkesbury", 
                                                      "Waipuru", "Upper Kai-Iwi", 
                                                      "Tainui", "SHCSBSB", "modern")) 
p.zh.age <- ggplot(data = mean_by_formation.meta) +
  geom_point(aes(x = as.numeric(mean.age), y = avg.zh,
                 col = formation),
             size = 5, shape = 17) + 
  theme(text = element_text(size = 16),
        #legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent', color=NA)) +
  scale_x_reverse(name = "Age (Ma)", limits = c(2.5, 0)) +
  scale_y_continuous(expression(Average~ln~Zooid~Height~(mu*m))) +
  scale_color_manual(values = col.form)

ggsave(p.zh.age, 
       file = "./Results/ln.zh.time.w.modern.png", 
       width = 14, height = 10, units = "cm")

df.form.meta <- merge(df, form.df,
                      by.x = "formation",
                      by.y = "formationCode")
          

traits.melt$formation <- factor(traits.melt$formation,
                                levels = c("NKLS", "NKBS", "Tewkesbury", 
                                           "Waipuru", "Upper Kai-Iwi", 
                                           "Tainui", "SHCSBSB", "modern")) 

box.ln.zh <- ggplot(data = traits.melt[traits.melt$measurementType == "ln.zh",], 
       aes(x = formation, 
           y = measurementValue, 
           fill = formation)) +
  geom_boxplot() +
  scale_color_manual(values = col.form) +
  scale_fill_manual(values = col.form) +
  ggtitle("Boxplots of LN Zooid Heights") +
  scale_x_discrete(name = "Formation",
                   guide = guide_axis(angle = 45)) +
  ylab(expression(ln~Zooid~Height~(mu*m))) + 
  theme(text = element_text(size = 16),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent', color=NA))

ggsave(box.ln.zh, 
       file = "./Results/boxplot.ln.zh.w.modern.png", 
       width = 14, height = 10, units = "cm")

#diff in variance across formations
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                 traits.melt$formation == "NKLS"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "NKBS"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "Tewkesbury"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "Waipuru"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "Upper Kai-Iwi"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "Tainui"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "SHCSBSB"])
var(traits.melt$measurementValue[traits.melt$measurementType == "ln.zh" &
                                     traits.melt$formation == "modern"])
#all very similar!! Except NKBS, which is quite high
