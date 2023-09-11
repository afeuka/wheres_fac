### Title: Fitting poisson regressions to FAC attendance data
### Author: Abbey Feuka
### Date: 19052023
####################
library(tidyverse)
library(lubridate)
library(scales)

#data cleaning -----------------------------------
##dates -----------------------------
dat <- read.csv("FAC attendance data.csv")
dat <- subset(dat,Date!="")
dat$Date <- as.Date(dat$Date,format="%m/%d/%Y")
dat$month <- month(dat$Date)
dat$year <- year(dat$Date)
dat$julianDate <- as.POSIXlt(dat$Date)$yday 
dat$julianDate.txt <- paste0(dat$month,"-",day(dat$Date))

##covariate corellations -----------------------------
cor(cbind(dat$Dist2Wagar,dat$Dist2BikePath,dat$Dist2Maxline)) #highly correlated

dat.test <- subset(dat,!is.na(Temp4pm) & !is.na(ChancePrecip))
cor(cbind(dat.test$julianDate,dat.test$Temp4pm,dat.test$ChancePrecip))

dat <- dat %>% select(-c(Dist2BikePath,Dist2Maxline))

##scale covariates -----------------------------
dat$SpeakerType <- as.factor(dat$SpeakerType)
dat$VenueType <- as.factor(dat$VenueType)
dat$Dist2Wagar.sc <- scale(dat$Dist2Wagar)[,1]
dat$Temp4pm.sc <- scale(dat$Temp4pm)[,1]
dat$ChancePrecip.sc <- scale(dat$ChancePrecip)[,1]

dat <- dat %>% pivot_longer(cols=colnames(dat)[grepl("Grad",colnames(dat)) | 
                                                 grepl("FacPost",colnames(dat)) |
                                                 grepl("Other",colnames(dat))],
                            names_to=c("Attendee","Time"),
                            names_sep="_", values_to="count")
dat$Time <- factor(dat$Time,levels=c("5pm","530pm","6pm"))
dat$preCOVID <- ifelse(dat$Date<as.Date("03-01-2020",format="%m-%d-%Y"),1,0)

#recode speakers
dat$SpeakerType[dat$SpeakerType%in%c("Post Doc","Visiting Scientist")] <- "Research Scientist"

#analysis -----------------------------
##fit poisson regressions -----------------------------
source("poisAIC_fun.R") 

covs <- c("Null","VenueType","SpeakerType","Dist2Wagar.sc",
          "Temp4pm.sc","ChancePrecip.sc","preCOVID")
grad <- dat %>% filter(Attendee=="Grad" & Time=="530pm")
grad.tab <- pois.aic(grad,covs)

facpost <- dat %>% filter(Attendee=="FacPost" & Time=="530pm")
fp.tab <- pois.aic(facpost,covs)

other <- dat %>% filter(Attendee=="Other" & Time=="530pm")
other.tab <- pois.aic(other,covs)

# write.csv(grad.tab$AICtab,"grad_aic.csv")
# write.csv(fp.tab$AICtab,"facpost_aic.csv")
# write.csv(other.tab$AICtab,"other_aic.csv")

##coefficients from top models -----------------------------
gmod <- grad.tab$AICtab$mod.idx[1]
gmod <- grad.tab$mods[[gmod]]

gmod.coeffs <- 
  cbind.data.frame(mn=coefficients(gmod),lci=confint(gmod)[,1],uci=confint(gmod)[,2],
                   attendee=rep("Graduate Students",length(coefficients(gmod))))
gmod.coeffs$factor <- rownames(gmod.coeffs)

fpmod <- fp.tab$AICtab$mod.idx[1]
fpmod <- fp.tab$mods[[fpmod]]
fpmod.coeffs <- 
  cbind.data.frame(mn=coefficients(fpmod),lci=confint(fpmod)[,1],uci=confint(fpmod)[,2],
                   attendee=rep("Faculty/Postdocs",length(coefficients(fpmod))))
fpmod.coeffs$factor <- rownames(fpmod.coeffs)

omod <- other.tab$AICtab$mod.idx[1]
omod <- other.tab$mods[[omod]]
omod.coeffs <- 
  cbind.data.frame(mn=coefficients(omod),lci=confint(omod)[,1],uci=confint(omod)[,2],
                   attendee=rep("Non-FWCB",length(coefficients(omod))))
omod.coeffs$factor <- rownames(omod.coeffs)

##top model coefficients -----------------------------
topmod.coeffs <- rbind.data.frame(gmod.coeffs,fpmod.coeffs,omod.coeffs)
rownames(topmod.coeffs) <- NULL

topmod.coeffs %>% filter(factor=="ChancePrecip.sc")

###top model coefficient plot -----------------------------
ggplot(topmod.coeffs)+geom_point(aes(x=mn,y=factor),size=2.5)+
  geom_linerange(aes(y=factor,xmin=lci,xmax=uci),lwd=1.5) +
  geom_vline(xintercept=0,lty=2,col="darkgrey")+
  facet_wrap(.~attendee) + ylab("Predictor") + xlab("Coefficient for FAC attendance")+
  scale_y_discrete(labels=c("Chance of precipitation","Distance to Wagar",
                            "Graduate student seminar","No seminar speaker",
                            "Research scientist seminar","Temperature"))+
  theme(text = element_text(size=15))
# ggsave("top_mod_coeffs.png",device="png",height=5,width=8,unit="in")

#model predictions -----------------------------
##distance to wagar -----------------------------
###other -----------------------------
cols <- hue_pal()(3)
ndat <- 50
newdat <- data.frame('SpeakerType'=factor(rep("Grad Student",ndat)),
                     "Dist2Wagar.sc"=seq(min(dat$Dist2Wagar.sc),max(dat$Dist2Wagar.sc),length.out=ndat),
                     "Temp4pm.sc"=rep(0,ndat))
omod1 <- glm(count ~ 0+ SpeakerType + Dist2Wagar.sc + Temp4pm.sc,family="poisson",data=other)

omod.pred <- predict(omod1,newdat,se.fit = T)
dist.sc <- scale(dat$Dist2Wagar)
omod.pred <- data.frame(Dist2Wagar=newdat$Dist2Wagar.sc*attr(dist.sc, 'scaled:scale') + 
                          attr(dist.sc, 'scaled:center'),
                        pred=exp(omod.pred$fit),
                        se=exp(omod.pred$se.fit))
omod.pred$Attendee<- "Non-FWCB"

###graduate -----------------------------
gmod1 <- glm(count ~ 0+ SpeakerType + Dist2Wagar.sc + ChancePrecip.sc,family="poisson",data=grad)
newdat <- data.frame('SpeakerType'=factor(rep("Grad Student",ndat)),
                     "Dist2Wagar.sc"=seq(min(dat$Dist2Wagar.sc),max(dat$Dist2Wagar.sc),length.out=ndat),
                     "ChancePrecip.sc"=rep(0,ndat))
gmod.pred <- predict(gmod1,newdat,se.fit = T)
gmod.pred <- data.frame(Dist2Wagar=newdat$Dist2Wagar.sc*attr(dist.sc, 'scaled:scale') + 
                          attr(dist.sc, 'scaled:center'),
                        pred=exp(gmod.pred$fit),
                        se=exp(gmod.pred$se.fit))

gmod.pred$Attendee<- "Graduate Student"
dist.pred <- full_join(gmod.pred,omod.pred)

ggplot(data=dist.pred)+
  geom_line(aes(x=Dist2Wagar,y=pred,col=Attendee))+
  geom_ribbon(aes(x=Dist2Wagar,ymax=pred+se,ymin=pred-se,fill=Attendee),alpha=0.5)+
  ylab("Mean number of happy hour attendees")+
  xlab("Distance to Wagar (m)")+
  scale_color_manual(values=cols[c(2,3)])+
  scale_fill_manual(values=cols[c(2,3)])+
  theme(axis.title = element_text(size=15),
        axis.text=element_text(size=12),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

# ggsave(filename = "dist_eff.jpeg",device="jpeg",
#        width=7,height=6,units="in")

##temperature -----------------------------
###other -----------------------------
ndat <- 50
newdat <- data.frame('SpeakerType'=factor(rep("Grad Student",ndat)),
                     "Dist2Wagar.sc"= rep(0,ndat),
                     "Temp4pm.sc"=seq(min(dat$Temp4pm.sc,na.rm=T),max(dat$Temp4pm.sc,na.rm=T),length.out=ndat))
omod1 <- glm(count ~ 0+ SpeakerType + Dist2Wagar.sc + Temp4pm.sc,family="poisson",data=other)

omod.pred <- predict(omod1,newdat,se.fit = T)
temp.sc <- scale(dat$Temp4pm)
omod.pred <- data.frame(Temp4pm=newdat$Temp4pm.sc*attr(temp.sc, 'scaled:scale') + 
                          attr(temp.sc, 'scaled:center'),
                        pred=exp(omod.pred$fit),
                        se=exp(omod.pred$se.fit))
omod.pred$Attendee<- "Non-FWCB"

###graduate students -----------------------------
gmod1 <- glm(count ~ 0+ SpeakerType + Dist2Wagar.sc + Temp4pm.sc,family="poisson",data=grad)

gmod.pred <- predict(gmod1,newdat,se.fit = T)
gmod.pred <- data.frame(Temp4pm=newdat$Temp4pm.sc*attr(temp.sc, 'scaled:scale') + 
                          attr(temp.sc, 'scaled:center'),
                        pred=exp(gmod.pred$fit),
                        se=exp(gmod.pred$se.fit))
gmod.pred$Attendee<- "Gradaute Student"

###faculty -----------------------------
fmod1 <- glm(count ~ 0+ SpeakerType +Temp4pm.sc,family="poisson",data=facpost)
newdat <- data.frame('SpeakerType'=factor(rep("Grad Student",ndat)),
                     "Temp4pm.sc"=seq(min(dat$Temp4pm.sc,na.rm=T),max(dat$Temp4pm.sc,na.rm=T),length.out=ndat))

fmod.pred <- predict(fmod1,newdat,se.fit = T)
fmod.pred <- data.frame(Temp4pm=newdat$Temp4pm.sc*attr(temp.sc, 'scaled:scale') + 
                          attr(temp.sc, 'scaled:center'),
                        pred=exp(fmod.pred$fit),
                        se=exp(fmod.pred$se.fit))
fmod.pred$Attendee<- "Faculty/Postdoc"

temp.pred <- full_join(gmod.pred,omod.pred) %>% full_join(fmod.pred)

ggplot(data=temp.pred)+
  geom_line(aes(x=Temp4pm,y=pred,col=Attendee))+
  geom_ribbon(aes(x=Temp4pm,ymax=pred+se,ymin=pred-se,fill=Attendee),alpha=0.5)+
  ylab("Mean number of happy hour attendees")+
  xlab("Temperature at 1600hr (C)")+
  scale_color_manual(values=cols)+
  scale_fill_manual(values=cols)+
  theme(axis.title = element_text(size=15),
        axis.text=element_text(size=12),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.position = "bottom",
        legend.box="horizontal")

# ggsave(filename = "temp_eff.jpeg",device="jpeg",
#        width=7,height=6,units="in")

##chance of precipitation -----------------------------
###graduate -----------------------------
gmod1 <- glm(count ~ 0+ SpeakerType + Dist2Wagar.sc + ChancePrecip.sc,family="poisson",data=grad)
newdat <- data.frame('SpeakerType'=factor(rep("Grad Student",ndat)),
                     "Dist2Wagar.sc"=rep(0,ndat),
                     "ChancePrecip.sc"=seq(min(dat$ChancePrecip.sc,na.rm=T),max(dat$ChancePrecip.sc,na.rm=T),length.out=ndat))
gmod.pred <- predict(gmod1,newdat,se.fit = T)
precip.sc <- scale(dat$ChancePrecip)
gmod.pred <- data.frame(ChancePrecip=newdat$ChancePrecip.sc*attr(precip.sc, 'scaled:scale') + 
                          attr(precip.sc, 'scaled:center'),
                        pred=exp(gmod.pred$fit),
                        se=exp(gmod.pred$se.fit))

gmod.pred$Attendee<- "Graduate Student"

###faculty/postdocs -----------------------------
fmod1 <- glm(count ~ 0+ SpeakerType + ChancePrecip.sc,family="poisson",data=facpost)
newdat <- data.frame('SpeakerType'=factor(rep("Grad Student",ndat)),
                     "ChancePrecip.sc"=seq(min(dat$ChancePrecip.sc,na.rm=T),max(dat$ChancePrecip.sc,na.rm=T),length.out=ndat))
fmod.pred <- predict(fmod1,newdat,se.fit = T)
fmod.pred <- data.frame(ChancePrecip=newdat$ChancePrecip.sc*attr(precip.sc, 'scaled:scale') + 
                          attr(precip.sc, 'scaled:center'),
                        pred=exp(fmod.pred$fit),
                        se=exp(fmod.pred$se.fit))

fmod.pred$Attendee<- "Faculty/Postdoc"
precip.pred <- full_join(gmod.pred,fmod.pred)

ggplot(data=precip.pred)+
  geom_line(aes(x=ChancePrecip,y=pred,col=Attendee))+
  geom_ribbon(aes(x=ChancePrecip,ymax=pred+se,ymin=pred-se,fill=Attendee),alpha=0.5)+
  ylab("Mean number of happy hour attendees")+
  xlab("Chance of precipitation at 1600hr (%)")+
  scale_color_manual(values=cols[c(1,2)])+
  scale_fill_manual(values=cols[c(1,2)])+
  theme(axis.title = element_text(size=15),
        axis.text=element_text(size=12),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15))

# ggsave(filename = "precip_eff.jpeg",device="jpeg",
#        width=7,height=6,units="in")

