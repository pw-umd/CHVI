library(psych);library(sp);library(sf);library(tigris)

`%notin%` <- Negate(`%in%`)
nys.shp <- st_read("shapefiles/nys_tract.shp")
nyc.shp <- st_read("shapefiles/nyc_tract.shp")

data.all <- read.csv("Data/data.all.v2.csv")

names(data.all)[3:46] <- c("Traffic","Truck_bus","PM2.5","Benzene","Wastewater","Remediation_site","Chemical_site","Oil_facility","Power_facility","Landfill",
                           "Waste_combustor","Metal_process","Inductrial_land","Vacant_home","Heat_projection","Flooding_risk","Non_vegetative","Agricultural_land","Hospital_time","income80",
                           "poverty100","No_bachelor","Unemployed","Single_parent","Hispanic","Black","Asian","Native","Limit_English","Redlining",
                           "Asthma","COPD","Heart_attack","Death","Low_birthweight","No_insurance","Disabled","Age65plus","Renter_home","Rental_cost",
                           "Energy_poverty","Mobile_home","Home1960","No_internet")

# sampling adequacy test
KMO(data.all[,3:21]) # risk group
KMO(data.all[,22:46]) # vulnerability group

# natural logarithm for the indicators that were highly skewed
data.all[,c(3:5,7:16,18,21,23,25:31,33:34,36,38:40,42:44,46)] <- data.all[,c(3:5,7:16,18,21,23,25:31,33:34,36,38:40,42:44,46)] + 0.0001 
data.all[,c(3:5,7:16,18,21,23,25:31,33:34,36,38:40,42:44,46)] <- lapply(data.all[,c(3:5,7:16,18,21,23,25:31,33:34,36,38:40,42:44,46)],log)

# normalization using z-score method
data.all[,3:46] <- scale(data.all[,3:46])
for (i in 3:46) { # PCA requires complete dataset so impute all missing as mean (0)
  data.all[,i] <- ifelse(is.na(data.all[,i]),0,data.all[,i])
}

# PCA
data.risk <- data.all[,2:21]
data.vuln <- data.all[,c(2,22:46)]

pca.risk <- prcomp(~ Traffic + Truck_bus + PM2.5 + Benzene + Wastewater + Remediation_site + Chemical_site + Oil_facility + Power_facility + Landfill + 
                     Waste_combustor + Metal_process + Inductrial_land + Vacant_home + Heat_projection + Flooding_risk + Non_vegetative + Agricultural_land + Hospital_time, data.risk)
screeplot(pca.risk,type = "l")
abline(h=1)
summary(pca.risk) # 5 PCs

pca.vuln <- prcomp(~ income80 + poverty100 + No_bachelor + Unemployed + Single_parent + Hispanic + Black + Asian + Native + Limit_English + Redlining + Asthma + COPD + Heart_attack + 
                     Death + Low_birthweight + No_insurance + Disabled + Age65plus + Renter_home + Rental_cost + Energy_poverty + Mobile_home + Home1960 + No_internet, data.vuln)
screeplot(pca.vuln,type = "l")
abline(h=1)
summary(pca.vuln) # 5 PCs

pca2.risk <- principal(data.risk[,2:20],nfactors=5,rotate="varimax")
pca2.vuln <- principal(data.vuln[,2:26],nfactors=5,rotate="varimax")
pca2.risk$loadings
pca2.vuln$loadings

score.final <- data.frame(matrix(nrow=nrow(data.all),ncol=5))
names(score.final) <- c("GEOID","score_risk","score_vuln","score_final","score_rank")
score.final$GEOID <- data.all$GEOID
for (c in c("risk","vuln")) {
  pca.scores.c <- data.frame(get(paste0("pca2.",c))$scores)
  # if (c=="risk") pca.scores.c[,3] <- pca.scores.c[,3]*(-1) # no need for direction correction
  n.RC <- ncol(pca.scores.c)
  
  for (j in 1:n.RC) {
    pca.scores.c[,j] <- ecdf(pca.scores.c[,j])(pca.scores.c[,j])*100
  }
  
  pca.scores.c$score <- rowSums(pca.scores.c)
  pca.scores.c <- cbind(get(paste0("data.",c))$GEOID,pca.scores.c)
  names(pca.scores.c) <- c("GEOID",paste0("RC",seq(1:n.RC),"_",c),paste0("score_",c))
  
  pca.scores.c[,paste0("score_",c)] <- ecdf(pca.scores.c[,paste0("score_",c)])(pca.scores.c[,paste0("score_",c)])*100
  
  if (c=="risk") score.final[,paste0("score_",c)] <- pca.scores.c[match(pca.scores.c$GEOID,score.final$GEOID),]$score_risk
  if (c=="vuln") score.final[,paste0("score_",c)] <- pca.scores.c[match(pca.scores.c$GEOID,score.final$GEOID),]$score_vuln
}
score.final$score_final <- score.final$score_risk*score.final$score_vuln
score.final$score_rank <- ecdf(score.final$score_final)(score.final$score_final)*100
score.final$dac <- ifelse(score.final$score_rank>=70,1,0)

score.final$GEOID <- as.character(score.final$GEOID)
main.pca.shp <- geo_join(nys.shp, score.final, "GEOID", "GEOID")

write.csv(score.final,"main.pca.csv")
st_write(main.pca.shp,"main.pca.shp",driver="ESRI Shapefile",layer="all")
main.pca.shp <- main.pca.shp[main.pca.shp$ALAND>0,]
st_write(main.pca.shp,"main.pca.land1.shp",driver="ESRI Shapefile",layer="all")
