library(sp);library(sf);library(tigris)

`%notin%` <- Negate(`%in%`)
nys.shp <- st_read("shapefiles/nys_tract.shp")
nyc.shp <- st_read("shapefiles/nyc_tract.shp")

data.all <- read.csv("Data/data.all.v2.csv")
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
data.all <- data.all[,c("X","GEOID","TRAFFIC","TRUCK","PM25","BENZENE","WASTEWATER","REME_AREA.y","CHEMICAL.y","OIL.y",
                        "POWER.y","LANDFILL.y","COMBUSTOR.y","METAL.y","INDUSTRIAL.y","VACANTHOME.y","HEAT90.y","FLOOD_RISK","NONVEG","AGRI_AREA.y",
                        "HOSPITAL","INCOME80","POVERTY100","NOBACHELOR.y","UNEMPLOYED.y","SGL_PARENT.y","HISPANIC.y","BLACK.y","ASIAN.y","NATIVE.y",
                        "ENGLISH.y","REDLINING","ASTHMA","COPD","MI","PREM_DEATH","LBW","NO_INS.y","DISABILITY.y","AGE65PLUS.y",
                        "RENTERHOME.y","RENT_COST","ENERGY","MOBILEHOME.y","HOME1960.y","NOINTERNET.y")]
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
