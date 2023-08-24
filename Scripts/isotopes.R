
oxy.18 <- read.table("~/GitHub/bryozoa/magnifica/Data/∂18O.csv",
                     header = TRUE,
                     sep = ",")

f.ages <- read.table("~/GitHub/bryozoa/magnifica/Data/Age Wanganui.csv",
                     header = TRUE,
                     sep = ",")

formations = f.ages$Formation_name[c(2, 3, 4, 5, 12, 14, 15)]
bottom = as.numeric(f.ages$Isotope_Stage_Start[c(2, 3, 4, 5, 12, 14, 15)])
top = as.numeric(f.ages$Isotope_Stage_End[c(2, 3, 4, 5, 12, 14, 15)])
med.O18 = NULL
sd.med.O18 = NULL
n.O18 = NULL
for (i in 1:length(formations)){
  temp = oxy.18$d18O[which(oxy.18$Time <= bottom[i] & oxy.18$Time >= top[i])]
  med.O18[i] = median(temp)
  sd.med.O18[i] = sd(temp)
  n.O18[i] <- length(temp)
}

O18 = cbind(cbind(med.O18, sd.med.O18, n.O18), cbind(bottom, top))
O18 = as.data.frame(O18)
rownames(O18) = formations

temp.NKBS = oxy.18$d18O[which(oxy.18$Time <= bottom[2] & oxy.18$Time >= top[2])]
temp.wai = oxy.18$d18O[which(oxy.18$Time <= bottom[4] & oxy.18$Time >= top[4])]
temp.uki = oxy.18$d18O[which(oxy.18$Time <= bottom[5] & oxy.18$Time >= top[5])]
