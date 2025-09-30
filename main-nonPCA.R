library(sp);library(sf);library(tigris)

`%notin%` <- Negate(`%in%`)
nys.shp <- st_read("shapefiles/nys_tract.shp")
nyc.shp <- st_read("shapefiles/nyc_tract.shp")

data.all <- read.csv("Data/data.all.v2.csv")

names(data.all)[3:46] <- c("Traffic","Truck_bus","PM2.5","Benzene","Wastewater","Remediation_site","Chemical_site","Oil_facility","Power_facility","Landfill",
                           "Waste_combustor","Metal_process","Inductrial_land","Vacant_home","Heat_projection","Flooding_risk","Non_vegetative","Agricultural_land","Hospital_time","income80",
                           "poverty100","No_bachelor","Unemployed","Single_parent","Hispanic","Black","Asian","Native","Limit_English","Redlining",
                           "Asthma","COPD","Heart_attack","Death","Low_birthweight","No_insurance","Disabled","Age65plus","Renter_home","Rental_cost",
                           "Energy_poverty","Mobile_home","Home1960","No_internet")

for (i in c(3:46)) {
  if (i %in% c(8:17,20,24:31,38:41,44:46)) { # make zero raw values rank zero, and calculate percentiles for all other non-zero raw values 
    data.all.i <- data.all[,c(2,i)]
    data.all.ia <- data.all.i[data.all.i[,2]!=0,]
    data.all.ib <- data.all.i[data.all.i[,2]==0,]
    data.all.ia[,2] <- ecdf(data.all.ia[,2])(data.all.ia[,2])*100
    data.all.ic <- rbind(data.all.ia,data.all.ib)
    
    data.all <- merge(data.all,data.all.ic,by="GEOID",all.x=TRUE)
    data.all <- data.all[,c(2,1,3:ncol(data.all))]
  } else data.all[,i] <- ecdf(data.all[,i])(data.all[,i])*100
}
data.all <- data.all[,c("X","GEOID","Traffic","Truck_bus","PM2.5","Benzene","Wastewater","Remediation_site.y","Chemical_site.y","Oil_facility.y",
                        "Power_facility.y","Landfill.y","Waste_combustor.y","Metal_process.y","Inductrial_land.y","Vacant_home.y","Heat_projection.y","Flooding_risk","Non_vegetative","Agricultural_land.y",
                        "Hospital_time","income80","poverty100","No_bachelor.y","Unemployed.y","Single_parent.y","Hispanic.y","Black.y","Asian.y","Native.y",
                        "Limit_English.y","Redlining","Asthma","COPD","Heart_attack","Death","Low_birthweight","No_insurance.y","Disabled.y","Age65plus.y",
                        "Renter_home.y","Rental_cost","Energy_poverty","Mobile_home.y","Home1960.y","No_internet.y")]
names(data.all) <- gsub(".y", "", names(data.all))

data.all$score.a1 <- apply(data.all[,3:7],1,mean,na.rm=TRUE)
data.all$score.a2 <- apply(data.all[,8:16],1,mean,na.rm=TRUE)
data.all$score.a3 <- apply(data.all[,17:21],1,mean,na.rm=TRUE)
data.all$score.b1 <- apply(data.all[,22:26],1,mean,na.rm=TRUE)
data.all$score.b2 <- apply(data.all[,27:32],1,mean,na.rm=TRUE)
data.all$score.b3 <- apply(data.all[,33:40],1,mean,na.rm=TRUE)
data.all$score.b4 <- apply(data.all[,41:46],1,mean,na.rm=TRUE)
data.all$score_risk <- apply(data.all[,47:49],1,function(x) weighted.mean(x,c(1,1,2),na.rm=TRUE))
# data.all$score_risk <- apply(data.all[,47:48],1,mean,na.rm=TRUE) ## remove the group of climate change risks
data.all$score_vuln <- apply(data.all[,50:53],1,mean,na.rm=TRUE)
# data.all$score_vuln <- apply(data.all[,c(50,52:53)],1,mean,na.rm=TRUE) ## remove the group of race and language
data.all$score_final <- data.all$score_risk*data.all$score_vuln
data.all$score_rank <- ecdf(data.all$score_final)(data.all$score_final)*100
data.all$dac <- ifelse(data.all$score_rank>=70,1,0)

data.all$GEOID <- as.character(data.all$GEOID)
main.nonpca.shp <- geo_join(nys.shp, data.all[,c(2,54:58)], "GEOID", "GEOID")

write.csv(data.all[,c(2,54:58)],paste0("main.nonpca.csv"))
suppressWarnings(st_write(main.nonpca.shp,paste0("main.nonpca.shp"),driver="ESRI Shapefile",layer="all"))
main.nonpca.shp <- main.nonpca.shp[main.nonpca.shp$ALAND>0,]
suppressWarnings(st_write(main.nonpca.shp,paste0("main.nonpca.land.shp"),driver="ESRI Shapefile",layer="all"))
