#these are extra analyses required by reviewers

#### EVOLUTIONARY MODELS ----
#include strict stasis
#mm = trait mean based on colony averages 
## (colony average will be the same as treating zooids as repeated measures of the same individual (treating colony as the real individual))
#vv = trait var based on average variance of colonies
#nn = sample size of total number of colonies
#tt = age based on formation

#goes from NKLS to modern

#number colonies per formation
trait.n = mean_by_formation$num.col

#reorganize formation order
mean_by_formation_colony$formation <- factor(mean_by_formation_colony$formation, 
                                             levels = c("NKLS", "NKBS", "Tewkesbury",
                                                        "Upper Kai-Iwi", "Tainui",
                                                        "SHCSBSB", "modern"))

#average of averages
mm.df = mean_by_formation_colony %>%
    dplyr::group_by(formation) %>%
    dplyr::summarise(mean.zl = mean(avg.zl),
                     mean.mpw.b = mean(avg.mpw.b),
                     mean.cw.m = mean(avg.cw.m),
                     mean.cw.d = mean(avg.cw.d),
                     mean.ow.m = mean(avg.ow.m),
                     mean.ol = mean(avg.ol),
                     mean.o.side = mean(avg.o.side),
                     mean.c.side = mean(avg.c.side)) %>%
    as.data.frame()

#average of variances
vv.df = mean_by_formation_colony %>%
    dplyr::group_by(formation) %>%
    dplyr::summarise(var.zl = mean(var.zl),
                     var.mpw.b = mean(var.mpw.b),
                     var.cw.m = mean(var.cw.m),
                     var.cw.d = mean(var.cw.d),
                     var.ow.m = mean(var.ow.m),
                     var.ol = mean(var.ol),
                     var.o.side = mean(var.o.side),
                     var.c.side = mean(var.c.side)) %>%
    as.data.frame()

#middle age in MYA
tt.df <- form.meta[-4, c(2, 13, 14)]
tt.df$Start_age <- as.numeric(tt.df$Start_age)
tt.df$End_age <- as.numeric(tt.df$End_age)
tt.df$Mid_age <- tt.df$Start_age - ((tt.df$Start_age-tt.df$End_age)/2)
trait.age <- tt.df$Mid_age

##ZL
pts.zl <- as.paleoTS(mm = mm.df$mean.zl,
                     vv = vv.df$var.zl,
                     nn = trait.n,
                     tt = trait.age,
                     oldest = "first")

zl.m1 <- fitSimple(pts.zl, model = "OU", method = "Joint", pool = FALSE)
zl.m2 <- fitSimple(pts.zl, model = "GRW", method = "Joint", pool = FALSE)
zl.m3 <- fitSimple(pts.zl, model = "URW", method = "Joint", pool = FALSE)
zl.m4 <- fitSimple(pts.zl, model = "Stasis", method = "Joint", pool = FALSE)
zl.m5 <- fitSimple(pts.zl, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(zl.m1, zl.m2, zl.m3, zl.m4, zl.m5)

## MPW.B
pts.mpw.b <- as.paleoTS(mm = mm.df$mean.mpw.b,
                        vv = vv.df$var.mpw.b,
                        nn = trait.n,
                        tt = trait.age,
                        oldest = "first")

mpw.b.m1 <- fitSimple(pts.mpw.b, model = "OU", method = "Joint", pool = FALSE)
mpw.b.m2 <- fitSimple(pts.mpw.b, model = "GRW", method = "Joint", pool = FALSE)
mpw.b.m3 <- fitSimple(pts.mpw.b, model = "URW", method = "Joint", pool = FALSE)
mpw.b.m4 <- fitSimple(pts.mpw.b, model = "Stasis", method = "Joint", pool = FALSE)
mpw.b.m5 <- fitSimple(pts.mpw.b, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(mpw.b.m1, mpw.b.m2, mpw.b.m3, mpw.b.m4, mpw.b.m5)

## CW.M
pts.cw.m <- as.paleoTS(mm = mm.df$mean.cw.m,
                       vv = vv.df$var.cw.m,
                       nn = trait.n,
                       tt = trait.age,
                       oldest = "first")

cw.m.m1 <- fitSimple(pts.cw.m, model = "OU", method = "Joint", pool = FALSE)
cw.m.m2 <- fitSimple(pts.cw.m, model = "GRW", method = "Joint", pool = FALSE)
cw.m.m3 <- fitSimple(pts.cw.m, model = "URW", method = "Joint", pool = FALSE)
cw.m.m4 <- fitSimple(pts.cw.m, model = "Stasis", method = "Joint", pool = FALSE)
cw.m.m5 <- fitSimple(pts.cw.m, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(cw.m.m1, cw.m.m2, cw.m.m3, cw.m.m4, cw.m.m5)

## CW.D
pts.cw.d <- as.paleoTS(mm = mm.df$mean.cw.d,
                       vv = vv.df$var.cw.d,
                       nn = trait.n,
                       tt = trait.age,
                       oldest = "first")

cw.d.m1 <- fitSimple(pts.cw.d, model = "OU", method = "Joint", pool = FALSE)
cw.d.m2 <- fitSimple(pts.cw.d, model = "GRW", method = "Joint", pool = FALSE)
cw.d.m3 <- fitSimple(pts.cw.d, model = "URW", method = "Joint", pool = FALSE)
cw.d.m4 <- fitSimple(pts.cw.d, model = "Stasis", method = "Joint", pool = FALSE)
cw.d.m5 <- fitSimple(pts.cw.d, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(cw.d.m1, cw.d.m2, cw.d.m3, cw.d.m4, cw.d.m5)

## OW.M
pts.ow.m <- as.paleoTS(mm = mm.df$mean.ow.m,
                       vv = vv.df$var.ow.m,
                       nn = trait.n,
                       tt = trait.age,
                       oldest = "first")

ow.m.m1 <- fitSimple(pts.ow.m, model = "OU", method = "Joint", pool = FALSE)
ow.m.m2 <- fitSimple(pts.ow.m, model = "GRW", method = "Joint", pool = FALSE)
ow.m.m3 <- fitSimple(pts.ow.m, model = "URW", method = "Joint", pool = FALSE)
ow.m.m4 <- fitSimple(pts.ow.m, model = "Stasis", method = "Joint", pool = FALSE)
ow.m.m5 <- fitSimple(pts.ow.m, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(ow.m.m1, ow.m.m2, ow.m.m3, ow.m.m4, ow.m.m5)

## OL
pts.ol <- as.paleoTS(mm = mm.df$mean.ol,
                     vv = vv.df$var.ol,
                     nn = trait.n,
                     tt = trait.age,
                     oldest = "first")

ol.m1 <- fitSimple(pts.ol, model = "OU", method = "Joint", pool = FALSE)
ol.m2 <- fitSimple(pts.ol, model = "GRW", method = "Joint", pool = FALSE)
ol.m3 <- fitSimple(pts.ol, model = "URW", method = "Joint", pool = FALSE)
ol.m4 <- fitSimple(pts.ol, model = "Stasis", method = "Joint", pool = FALSE)
ol.m5 <- fitSimple(pts.ol, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(ol.m1, ol.m2, ol.m3, ol.m4, ol.m5)

## O SIDE
pts.o.side <- as.paleoTS(mm = mm.df$mean.o.side,
                         vv = vv.df$var.o.side,
                         nn = trait.n,
                         tt = trait.age,
                         oldest = "first")

o.side.m1 <- fitSimple(pts.o.side, model = "OU", method = "Joint", pool = FALSE)
o.side.m2 <- fitSimple(pts.o.side, model = "GRW", method = "Joint", pool = FALSE)
o.side.m3 <- fitSimple(pts.o.side, model = "URW", method = "Joint", pool = FALSE)
o.side.m4 <- fitSimple(pts.o.side, model = "Stasis", method = "Joint", pool = FALSE)
o.side.m5 <- fitSimple(pts.o.side, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(o.side.m1, o.side.m2, o.side.m3, o.side.m4, o.side.m5)

## C SIDE
pts.c.side <- as.paleoTS(mm = mm.df$mean.c.side,
                         vv = vv.df$var.c.side,
                         nn = trait.n,
                         tt = trait.age,
                         oldest = "first")

c.side.m1 <- fitSimple(pts.c.side, model = "OU", method = "Joint", pool = FALSE)
c.side.m2 <- fitSimple(pts.c.side, model = "GRW", method = "Joint", pool = FALSE)
c.side.m3 <- fitSimple(pts.c.side, model = "URW", method = "Joint", pool = FALSE)
c.side.m4 <- fitSimple(pts.c.side, model = "Stasis", method = "Joint", pool = FALSE)
c.side.m5 <- fitSimple(pts.c.side, model = "StrictStasis", method = "Joint", pool = FALSE)
compareModels(c.side.m1, c.side.m2, c.side.m3, c.side.m4, c.side.m5)

#### HALDANES ----
#Rate (haldanes) = (log(X2) -log(xI))/alogx (t2 - tl)/((t2-t1)/tgen)
#alogx = pooled sd of the two times
#just zl
#assume generation time of 2 years
#second - first
#t1 = NKBS - NKLS
#t2 = tewk - NKBS
#t3 = uki - tewk
#t4 = tai - uki
#t5 = SHCSBSB - tai
#t6 = mod - SHCSBSB

#time interval
#mid-age / generation time
#mid-age = start age - (start age - end age/2)

#age is in Ka
form.meta$Start_age <- as.numeric(form.meta$Start_age)
form.meta$End_age <- as.numeric(form.meta$End_age)
form.meta$Mid_age <- form.meta$Start_age - ((form.meta$Start_age-form.meta$End_age)/2)

##pooled std:
sd(df$ln.zl[df$formation == "NKBS" | df$formation == "NKLS"], na.rm = T)
sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "NKBS"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "NKBS"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "NKLS"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "NKBS"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "NKBS"] + mean_by_formation$num.col[mean_by_formation$formation == "NKLS"]-2))
#very similar

t1_h <- (abs((form.meta$Mid_age[form.meta$formationCode == "NKBS"]-form.meta$Mid_age[form.meta$formationCode == "NKLS"])/2))*1000000
t2_h <- (abs((form.meta$Mid_age[form.meta$formationCode == "Tewkesbury"]-form.meta$Mid_age[form.meta$formationCode == "NKBS"])/2))*1000000
t3_h <- (abs((form.meta$Mid_age[form.meta$formationCode == "Upper Kai-Iwi"]-form.meta$Mid_age[form.meta$formationCode == "Tewkesbury"])/2))*1000000
t4_h <- (abs((form.meta$Mid_age[form.meta$formationCode == "Tainui"]-form.meta$Mid_age[form.meta$formationCode == "Upper Kai-Iwi"])/2))*1000000
t5_h <- (abs((form.meta$Mid_age[form.meta$formationCode == "SHCSBSB"]-form.meta$Mid_age[form.meta$formationCode == "Tainui"])/2))*1000000
t6_h <- (abs((form.meta$Mid_age[form.meta$formationCode == "modern"]-form.meta$Mid_age[form.meta$formationCode == "SHCSBSB"])/2))*1000000
tgen <- c(t1_h, t2_h, t3_h, t4_h, t5_h, t6_h)

zl_s_t1 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "NKBS"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "NKBS"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "NKLS"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "NKLS"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "NKBS"] + mean_by_formation$num.col[mean_by_formation$formation == "NKLS"]-2))
zl_s_t2 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "Tewkesbury"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "NKBS"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "NKBS"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "NKLS"] + mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"]-2))
zl_s_t3 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "Upper Kai-Iwi"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "Tewkesbury"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"] + mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"]-2))
zl_s_t4 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "Tainui"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "Tainui"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "Upper Kai-Iwi"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"] + mean_by_formation$num.col[mean_by_formation$formation == "Tainui"]-2))
zl_s_t5 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "SHCSBSB"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "Tainui"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "Tainui"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "Tainui"] + mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"]-2))
zl_s_t6 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "modern"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "modern"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"]-1)*(mean_by_formation$sd.zl[mean_by_formation$formation == "SHCSBSB"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"] + mean_by_formation$num.col[mean_by_formation$formation == "modern"]-2))

zl_t1 <- ((mean_by_formation$avg.zl[mean_by_formation$formation == "NKBS"] - mean_by_formation$avg.zl[mean_by_formation$formation == "NKLS"])/zl_s_t1)/t1_h
zl_t2 <- ((mean_by_formation$avg.zl[mean_by_formation$formation == "Tewkesbury"] - mean_by_formation$avg.zl[mean_by_formation$formation == "NKBS"])/zl_s_t2)/t2_h
zl_t3 <- ((mean_by_formation$avg.zl[mean_by_formation$formation == "Upper Kai-Iwi"] - mean_by_formation$avg.zl[mean_by_formation$formation == "Tewkesbury"])/zl_s_t3)/t3_h
zl_t4 <- ((mean_by_formation$avg.zl[mean_by_formation$formation == "Tainui"] - mean_by_formation$avg.zl[mean_by_formation$formation == "Upper Kai-Iwi"])/zl_s_t4)/t4_h
zl_t5 <- ((mean_by_formation$avg.zl[mean_by_formation$formation == "SHCSBSB"] - mean_by_formation$avg.zl[mean_by_formation$formation == "Tainui"])/zl_s_t5)/t5_h
zl_t6 <- ((mean_by_formation$avg.zl[mean_by_formation$formation == "modern"] - mean_by_formation$avg.zl[mean_by_formation$formation == "SHCSBSB"])/zl_s_t6)/t6_h

zl_haldanes <- c(zl_t1, zl_t2, zl_t3, zl_t4, zl_t5, zl_t6)

summary(glm(zl_haldanes ~ tgen))
summary(glm(log(abs(zl_haldanes)) ~ log(tgen)))
exp(summary(glm(log(abs(zl_haldanes)) ~ log(tgen)))$coefficient[1])

cw.m_s_t1 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "NKBS"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "NKBS"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "NKLS"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "NKLS"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "NKBS"] + mean_by_formation$num.col[mean_by_formation$formation == "NKLS"]-2))
cw.m_s_t2 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "Tewkesbury"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "NKBS"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "NKBS"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "NKLS"] + mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"]-2))
cw.m_s_t3 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "Upper Kai-Iwi"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "Tewkesbury"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "Tewkesbury"] + mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"]-2))
cw.m_s_t4 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "Tainui"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "Tainui"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "Upper Kai-Iwi"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "Upper Kai-Iwi"] + mean_by_formation$num.col[mean_by_formation$formation == "Tainui"]-2))
cw.m_s_t5 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "SHCSBSB"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "Tainui"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "Tainui"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "Tainui"] + mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"]-2))
cw.m_s_t6 <- sqrt(((mean_by_formation$num.col[mean_by_formation$formation == "modern"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "modern"])^2 + (mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"]-1)*(mean_by_formation$sd.cw.m[mean_by_formation$formation == "SHCSBSB"]^2))/(mean_by_formation$num.col[mean_by_formation$formation == "SHCSBSB"] + mean_by_formation$num.col[mean_by_formation$formation == "modern"]-2))

cw.m_t1 <- ((mean_by_formation$avg.cw.m[mean_by_formation$formation == "NKBS"] - mean_by_formation$avg.cw.m[mean_by_formation$formation == "NKLS"])/cw.m_s_t1)/t1_h
cw.m_t2 <- ((mean_by_formation$avg.cw.m[mean_by_formation$formation == "Tewkesbury"] - mean_by_formation$avg.cw.m[mean_by_formation$formation == "NKBS"])/cw.m_s_t2)/t2_h
cw.m_t3 <- ((mean_by_formation$avg.cw.m[mean_by_formation$formation == "Upper Kai-Iwi"] - mean_by_formation$avg.cw.m[mean_by_formation$formation == "Tewkesbury"])/cw.m_s_t3)/t3_h
cw.m_t4 <- ((mean_by_formation$avg.cw.m[mean_by_formation$formation == "Tainui"] - mean_by_formation$avg.cw.m[mean_by_formation$formation == "Upper Kai-Iwi"])/cw.m_s_t4)/t4_h
cw.m_t5 <- ((mean_by_formation$avg.cw.m[mean_by_formation$formation == "SHCSBSB"] - mean_by_formation$avg.cw.m[mean_by_formation$formation == "Tainui"])/cw.m_s_t5)/t5_h
cw.m_t6 <- ((mean_by_formation$avg.cw.m[mean_by_formation$formation == "modern"] - mean_by_formation$avg.cw.m[mean_by_formation$formation == "SHCSBSB"])/cw.m_s_t6)/t6_h

cw.m_haldanes <- c(cw.m_t1, cw.m_t2, cw.m_t3, cw.m_t4, cw.m_t5, cw.m_t6)

summary(glm(cw.m_haldanes ~ tgen))
summary(glm(log(abs(cw.m_haldanes)) ~ log(tgen)))
exp(-11.9298)#intercept
