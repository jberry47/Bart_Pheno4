library(plyr)
library(ggplot2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
library(reshape2)
library(grid)
library(stringr)
library(RColorBrewer)
library(scales)
library(multcomp)
library(lme4)
library(lsmeans)
library(gridExtra)
library(ggdendro)
library(dendextend)
library(dendextendRcpp)
library(survminer)
library(ggridges)

setwd("/media/jberry/Extra Drive 1/Danforth/Sorghum/Pheno4")
img_to_barcode <- read.csv("SnapshotInfo.csv",header = T,stringsAsFactors = F)
img_to_barcode <- img_to_barcode[img_to_barcode$tiles != "",]
colnames(img_to_barcode)[3] <- "Barcodes"
img_to_barcode <- img_to_barcode[,c("id","Barcodes","timestamp","car.tag")]
sv_shapes <- read.table("shapes.txt",header = F,stringsAsFactors = F,sep = " ")
colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","oof","det")
sv_shapes$id <- substring(as.character(sapply(sv_shapes$meta,function(i) strsplit(i,"/")[[1]][2])),9)
sv_shapes$imgname <- as.character(sapply(sv_shapes$meta,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
sv_shapes <- join(sv_shapes,img_to_barcode,by="id")
assoc <- read.csv("pheno4_assoc.csv",header=T,stringsAsFactors = F)
sv_shapes <- join(sv_shapes,assoc,by="Barcodes")
sv_shapes <- sv_shapes[!is.na(sv_shapes$Drought),]
sv_shapes$timestamp <- as.POSIXct(strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S"))
beg <- min(sv_shapes$timestamp)
sv_shapes$DAP <- floor(as.numeric((sv_shapes$timestamp - beg)/60/60/24))+4
#sv_shapes$Drought <- ordered(sv_shapes$Drought, levels=c("AAA","ABB","ABA"))
#sv_shapes$Genotype <- "BTx623"
sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
sv_shapes <- sv_shapes[sv_shapes$Microbes != "Blank",]
sv_shapes$Genotype <- sapply(str_sub(sv_shapes$Barcodes,5,5),function(i)if(i==1){"BTx623"}else if(i==2){"San Chi San"}else{"Grassl"})
empties <- sv_shapes[sv_shapes$DAP == 20 & sv_shapes$area < 10,]
sv_shapes <- sv_shapes[!(sv_shapes$Barcodes %in% empties$Barcodes),]
sv_shapes$Microbes[sv_shapes$Microbes == "Mix(WW+WS)"] <- "Mix"
sv_shapes <- sv_shapes[(read.csv("pheno4_outliers_bool.csv",header = T,stringsAsFactors = F)$x),]
area_convert <- 13.2*3.7/46856
sv_shapes$blasted <- sv_shapes$car.tag %in% c(1304,1019)
tail(sv_shapes)
sv_shapes[sv_shapes$Genotype == "Grassl" & sv_shapes$DAP == 18 & sv_shapes$Drought == "WW_80" & sv_shapes$Microbes == "Control",]

write.csv(sv_shapes,"bart_pheno4_raw.csv",row.names = F,quote = F)

aggregate(data=empties,Barcodes~Microbes+Drought+Genotype,FUN=function(i)length(unique(i)))
length(unique(empties$Barcodes))

#*************************************************************************************************
# Watering data
#*************************************************************************************************
img_to_barcode <- read.csv("SnapshotInfo.csv",header = T,stringsAsFactors = F)
img_to_barcode <- img_to_barcode[img_to_barcode$tiles == "",]
img_to_barcode$timestamp <- as.POSIXct(strptime(img_to_barcode$timestamp,format = "%Y-%m-%d %H:%M:%S"))
colnames(img_to_barcode)[3] <- "Barcodes"
head(img_to_barcode)

ggplot(img_to_barcode,aes(timestamp,weight.before))+
    geom_point()


#*************************************************************************************************
# Outlier Detection - Done above, don't rerun
#*************************************************************************************************
library(gputools)
chooseGpu(1)
sv_shapes <- sv_shapes[order(sv_shapes$timestamp),]
cooksd <- cooks.distance(gpuGlm(data=sv_shapes,area~Microbes:Drought:Genotype:as.factor(DAP)))
sv_shapes1 <- sv_shapes[cooksd > 3*mean(cooksd),]
plot(cooksd)
abline(h=3*mean(cooksd),col="red")
head(sv_shapes1)

write.csv(cooksd < 3*mean(cooksd),"pheno4_outliers_bool.csv",row.names = F,quote = F)
write.csv(sv_shapes1,"pheno4_area_outliers.csv",row.names = F,quote = F)


#*************************************************************************************************
# Positional effects
#*************************************************************************************************
statefile <- NULL
for(my_dir in list.dirs("Pheno4_statefile")[-1]){
  for(i in list.files(my_dir)){
    temp <- read.table(paste0(substr(my_dir,0,24),"/",i),header = F,stringsAsFactors = F,sep = ";",skip = 1)
    num <- as.numeric(read.table(paste0(substr(my_dir,0,204),"/",i),header = F,stringsAsFactors = F,sep = ";",nrows = 1))
    meta <- substr(strsplit(i,"[.]")[[1]][3],5,7)
    statefile <- rbind(statefile,data.frame("Date"=substr(my_dir,18,24),"GH"=strsplit(meta,"_")[[1]][1],"Lane"=as.numeric(strsplit(meta,"_")[[1]][2]),"Position"=1:num,"Barcodes"=temp$V3,stringsAsFactors = F))
  }
}
head(statefile)
statefile$GH <- ordered(statefile$GH,levels=5:1)
statefile$timestamp <- as.POSIXlt(strptime(statefile$Date,format = "%m%d%y"))
statefile$Date <- paste0(lubridate::month(statefile$timestamp),lubridate::day(statefile$timestamp),lubridate::year(statefile$timestamp))
tail(statefile)
range(statefile$Lane[statefile$GH==1])

dat <- aggregate(data=sv_shapes,area~Barcodes+timestamp+Drought+Microbes+DAP,"mean")
dat$Date <- paste0(lubridate::month(dat$timestamp),lubridate::day(dat$timestamp),substr(lubridate::year(dat$timestamp),3,4))
dat <- join(statefile,dat,by=c("Barcodes","Date"))
dat <- na.omit(dat[,-6])
tail(dat)
gh13 <- unique(dat$Barcodes[dat$GH %in% c(1,3)])


#library(raster)
library(gstat)
library(sp)
#library(dismo)
#library(rgeos)

all_spat <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
for(day in 9:21){
  sub <- dat[dat$DAP == day & dat$Drought == "WS_25",]
  sub$area_n <- with(sub,(area-mean(area))/sd(area))
  sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
  coordinates(sub) <- ~Lane1+Position
  vgm1 <- variogram(area_n~1,sub)
  fit <- fit.variogram(vgm1,model=vgm(0.3,"Sph"))
  sub.grid <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
  coordinates(sub.grid) <- ~x+y
  kriged <- krige(area_n~1,sub,sub.grid,fit,maxdist=3)
  out <- setNames(data.frame(kriged)[,1:3],c("x","y",paste0("pred_",day)))
  all_spat <- join(all_spat,out,by=c("x","y"))  
}
all_spat$avg <- rowMeans(all_spat[,-(1:2)],na.rm = T)

p <- ggplot(all_spat,aes(x,y))+
  geom_tile(aes(fill=avg))+
  geom_vline(xintercept = c(8,16,24)+0.5)+
  scale_x_reverse()+
  scale_fill_gradient2(mid = "gray95",high="darkgreen",limits=c(-2.5,2.5))+
  theme_light()
p
ggsave("pheno4_wsOnly_allDAPavg_areaNorm_spatial.png",width=6.9,height=4.9,plot = p, dpi = 300)

#*************************************************************************************************
# Image variation 
#*************************************************************************************************
p <- ggplot(sv_shapes,aes(hour,det))+
  geom_jitter(width = 0.5)+
  scale_x_continuous(breaks = seq(from=0,to=24,by=2))+
  scale_y_continuous(limits=c(-200,2))+
  ylab("Correction Strength")+
  xlab("Time of Day (hr)")+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(legend.position='none')
p
ggsave("pheno4_image_variability.png",width=4.89,height=4.69,plot = p, dpi = 300)


#*************************************************************************************************
# Trends
#*************************************************************************************************
p <- ggplot(sv_shapes[sv_shapes$Microbes %in% c("Control","WW","WS","Mix"),],aes(DAP,area*area_convert))+
  facet_grid(Genotype~Drought)+
  geom_smooth(aes(color=Microbes))+
  ylab(~~Area~(cm^2))+
  scale_color_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno4_field_area_trends.png",width=6.69,height=6.29,plot = p, dpi = 300)

p <- ggplot(sv_shapes[sv_shapes$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(DAP,area*area_convert))+
  facet_grid(Genotype~Drought)+
  geom_smooth(aes(color=Microbes))+
  ylab(~~Area~(cm^2))+
  scale_color_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno4_dangl_area_trends.png",width=6.69,height=6.29,plot = p, dpi = 300)

p <- ggplot(sv_shapes[sv_shapes$Drought == "WS_25" & sv_shapes$Microbes %in% c("Control","SynCom A") & sv_shapes$Genotype == "Grassl",], aes(DAP,area))+
  facet_grid(~Microbes)+
  geom_smooth(aes(color=blasted,group=Barcodes),se=F)
p
head(sv_shapes)


aggregate(data=sv_shapes[sv_shapes$Drought == "WS_25" & sv_shapes$Microbes %in% c("Control","SynCom A") & sv_shapes$Genotype == "Grassl",],
  DAP ~ Barcodes,FUN = "max")
  


#*************************************************************************************************
# Boxplots
#*************************************************************************************************
p <- ggplot(sv_shapes[sv_shapes$DAP == 19 & sv_shapes$Genotype == "BTx623" &sv_shapes$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(Microbes,area))+
  facet_grid(Genotype~Drought)+
  geom_violin(aes(fill=Microbes,color=Microbes),alpha=.2)+
  geom_boxplot(width=.25)+
  ylab(~~Area~(cm^2))+
  xlab("")+
  scale_color_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  scale_fill_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position='none')
p
ggsave("pheno4_dangl_area_dap19_boxplots.png",width=4.8,height=6.29,plot = p, dpi = 300)


#*************************************************************************************************
# R^2 and Parital Correlations
#*************************************************************************************************
dat <- sv_shapes[sv_shapes$DAP==20,]

shapes <- c("area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd")
H2 <- c()
for(e in shapes){
  model <- lmer(eval(parse(text=e))~(1|Drought)+(1|Microbes)+(1|Genotype),data = dat)
  re<- VarCorr(model)
  res<-attr(VarCorr(model), "sc")^2
  des1.var <- as.numeric(attr(re[["Drought"]],"stddev"))^2
  des2.var <- as.numeric(attr(re[["Microbes"]],"stddev"))^2
  des3.var <- as.numeric(attr(re[["Genotype"]],"stddev"))^2
  tot.var<-sum(as.numeric(re),res)
  unexp <- 1-sum(as.numeric(re))/sum(as.numeric(re),res)
  h2 <- c((des1.var/tot.var),
    (des2.var/tot.var),
    (des3.var/tot.var),
    unexp)
  H2 <- rbind(H2,h2)
}
H2 <- data.frame(H2,row.names = shapes)
H2$Shape <- rownames(H2)
rownames(H2) <- NULL
colnames(H2) <- c("Drought","Microbes","Genotype","Unexplained","Shape")
H2$Shape <-  ordered(H2$Shape,levels=H2$Shape[order(H2$Unexplained)])
H2_melt <- melt(H2,id=c("Shape"))
H2_melt$variable <- ordered(H2_melt$variable,levels=c("Unexplained","Drought","Microbes","Genotype"))
head(H2_melt)
p <- ggplot(data=H2_melt,aes(Shape,value*100))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text = element_text(size = 14),
    axis.title.y= element_text(size = 18),
    axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
    plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
  theme(legend.position = "top")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno4_dap20_anova.png",width=6.17,height=4.54,plot = p, dpi = 300)


#*************************************************************************************************
# OOF Survival
#*************************************************************************************************
head(sv_shapes)
dat <- sv_shapes
dat$srv <- with(dat,Surv(time=DAP,event=oof))

mod1 <- summary(survfit(srv ~ Drought+Genotype+Microbes, data = dat, conf.type = "log-log"),time=min(dat$DAP):max(dat$DAP))
mod_df <- data.frame("DAP"=mod1$time,"strata"=as.character(mod1$strata),"surv"=mod1$surv,"low"=mod1$lower,"high"=mod1$upper,stringsAsFactors = F)
mod_df <- cbind(mod_df,setNames(data.frame(sapply(des,function(m){unlist(lapply(str_split(mod_df$strata,", "),function(i) trimws(str_split(i[str_detect(i,m)],"=")[[1]][2])))}),stringsAsFactors = F),des))

head(mod_df)
str(mod_df)
p <- ggplot(mod_df[mod_df$Drought == "WW_80",],aes(DAP,surv))+
  facet_grid(~Microbes)+
  #geom_ribbon(aes(ymin=low,ymax=high),fill="gray60")+
  geom_line(aes(color=Genotype))+
  ylab("Out Of Frame Risk")+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,.2))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
ggsave("pheno4_ww_oof_risk.png",width=9.95,height=3.9,plot = p, dpi = 300)


#*************************************************************************************************
# Fresh weight to area
#*************************************************************************************************
fw <- read.csv("pheno4_weights.csv",header = T,stringsAsFactors = F)
fw <- join(fw, sv_shapes[sv_shapes$DAP == 20,],by="Barcodes")
fw$Genotype
p <- ggplot(fw,aes(Fresh_weight,area*area_convert))+
  facet_grid(~Genotype)+
  geom_point(aes(color=Drought))+
  geom_smooth(aes(color=Drought),method = "lm",se=F)+
  ylab(~~Area~(cm^2))+
  xlab("Fresh Weight (g)")+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno4_weight-biomass.png",width=8.58,height=4.34,plot = p, dpi = 300)
head(fw)

coef(lm(data=fw,area~0+Genotype:Drought+Fresh_weight:Genotype:Drought))*area_convert


#*************************************************************************************************
# Specific t-tests
#*************************************************************************************************
dat <- sv_shapes[sv_shapes$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C") & sv_shapes$Genotype == "San Chi San" & sv_shapes$Drought == "WS_25",]
t.test(data=dat[dat$DAP == 18 & dat$Microbes %in% c("SynCom B","SynCom A"),],(area*area_convert)~Microbes)

head(dat)


#*************************************************************************************************
# Helper functions
#*************************************************************************************************
get_color <- function(file_name,start,stop){
  color_data <- read.table(file_name,header = F,stringsAsFactors = F,sep = " ")[,-257]
  color_data$id <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][2],"snapshot")[[1]][2]))
  color_data$imgname <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
  color_data <- join(color_data,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
  color_data <- join(color_data,assoc,by="Barcodes")
  color_data$timestamp <- strptime(color_data$timestamp,format = "%Y-%m-%d %H:%M:%S")
  color_data$DAP <- floor(as.numeric((color_data$timestamp - beg)/60/60/24))+2
  color_data[,start:stop] <- t(apply(color_data[,start:stop],1,function(i){i/(sum(i,na.rm = T)+1)}))*100
  color_data$hr <- as.POSIXlt(color_data$timestamp)$hour
  #color_data$status <- sapply(color_data$hr,function(i) if(i %in% 7:21){"Daytime"}else{"Nighttime"})
  return(color_data)
}

hist_avg <- function(data,start,stop){
  sub <- data
  test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
    data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
      data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
        colMeans(m[,start:stop],na.rm = T)
      }
      )))
    }
    )))
  })))
  return(test)
}

hist_sd <- function(data,start,stop){
  sub <- data
  test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
    data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
      data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
        apply(m[,start:stop],2,function(i){sd(i,na.rm = T)})
      }
      )))
    }
    )))
  })))
  return(test)
}


#*************************************************************************************************
# NIR Data
#*************************************************************************************************
nir <- get_color("nir2.txt",2,256)
nir$intensityAVG <- apply(nir[,3:255],1,function(i){sum((i/100)*(2:254),na.rm = T)})
nir$Microbes[nir$Microbes == "Mix(WW+WS)"] <- "Mix"
nir$Microbes <- ordered(nir$Microbes, levels=c("Control","WW","WS","Mix","SynCom A","SynCom B","SynCom C"))
nir$Genotype <- sapply(str_sub(nir$Barcodes,5,5),function(i)if(i==1){"BTx623"}else if(i==2){"San Chi San"}else{"Grassl"})
nir <- nir[!is.na(nir$Microbes),]
nir <- nir[!(nir$Barcodes %in% empties),]
outliers <- read.csv("pheno4_area_outliers.csv",header = T,stringsAsFactors = F)
outliers$camera_angle <- unlist(lapply(strsplit(outliers$meta,"_"),function(i) i[3]))
outliers$unique_id <- paste(outliers$Barcodes,outliers$DAP,outliers$camera_angle,sep="_")
nir$camera_angle <- unlist(lapply(strsplit(as.character(nir$V1),"_"),function(i) i[3]))
nir$unique_id <- paste(nir$Barcodes,nir$DAP,nir$camera_angle,sep="_")
nir <- nir[!(nir$unique_id %in% outliers$unique_id),]
cooksd_nir <- cooks.distance(lm(data=nir,intensityAVG~Microbes:Drought:Genotype:as.factor(DAP)))
plot(cooksd_nir)
abline(h=3*mean(cooksd_nir),col="red")
nir <- nir[cooksd < 3*mean(cooksd_nir),]

head(nir)

#heatmap broken out by microbes
test <- aggregate(data=nir[nir$intensityAVG != 0 & nir$DAP >4,],intensityAVG~Drought+Microbes+Genotype+DAP,FUN = function(i)mean(i,na.rm=T))
head(test)
p <- ggplot(test[test$DAP >5,],aes(as.factor(DAP),Microbes))+
  facet_grid(Genotype~Drought)+
  geom_tile(aes(fill=intensityAVG))+
  scale_fill_gradient2(limits=c(85,105),midpoint = 0.3+mean(nir$intensityAVG[nir$Drought == "WW_80" & nir$intensityAVG != 0 & nir$DAP == 15]),high ="gray10",low= "#56B1F7",mid = "#d7e4ef")+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title.x = element_blank())
p
ggsave("pheno4_nir_intensity_heatmap.png",width=8.75,height=6.29,plot = p, dpi = 300)


max(sv_shapes$DAP)
