### Title: Function for fitting poissson GLMs to fac count data
### Author: Abbey Feuka
### Date: 25052023
### Notes: used in "fac_analysis.R"
######################

pois.aic <- function(dat,covs){
  gtable <- data.frame(Model=covs,
                       AIC=rep(NA,length(covs)),
                       Deviance=rep(NA,length(covs)))
  #single covariates
  g <- list()
  g[[1]] <- glm(count ~ 1 ,family="poisson", data=dat)
  gtable$AIC[1] <- g[[1]]$aic
  gtable$Deviance[1] <-g[[1]]$deviance
  gtable$mod.idx <- 1:nrow(gtable)
  for(i in 2:length(covs)){
    if(covs[i]%in%c("VenueType","SpeakerType")){
      g[[i]] <- glm(count ~ 0+ get(covs[i]) ,family="poisson", data=dat)
    } else {
      g[[i]] <- glm(count ~ get(covs[i]) ,family="poisson", data=dat)
    }
    gtable$AIC[i] <- g[[i]]$aic
    gtable$Deviance[i] <-g[[i]]$deviance
  }
  # gtable[order(gtable$AIC,decreasing=F),]
  poorfit <- gtable$Model[which(gtable$AIC>gtable$AIC[gtable$Model=="Null"])]
  
  #two covariates
  combos2 <-combn(covs[-(which(covs%in%c("Null",poorfit)))],2)
  combos2 <- sapply(1:ncol(combos2),function(i)paste0(combos2[1,i]," + ",combos2[2,i]))
  gtable<-rbind(gtable,data.frame(Model=combos2,
                                  AIC=rep(NA,length(combos2)),
                                  Deviance=rep(NA,length(combos2)),
                                  mod.idx=nrow(gtable)+1:length(combos2)))
  extra<- length(g)
  for(i in 1:length(combos2)){
    if(grepl("SpeakerType",combos2[i]) | grepl("VenueType",combos2[i])){
      g[[i+extra]] <- glm(paste("count", "~0+", combos2[i]),
                          family="poisson", data=dat)
    } else {
      g[[i+extra]] <- glm(paste("count", "~", combos2[i]),
                          family="poisson", data=dat)
    }
    
    gtable$AIC[i+extra] <- g[[i+extra]]$aic
    gtable$Deviance[i+extra] <-g[[i+extra]]$deviance
  }
  # gtable[order(gtable$AIC,decreasing=F),]
  
  #three covariates
  combos3 <-combn(covs[-(which(covs%in%c("Null",poorfit)))],3)
  combos3 <- sapply(1:ncol(combos3),function(i)paste0(combos3[1,i]," + ",combos3[2,i]," + ",
                                                      combos3[3,i]))
  gtable<-rbind(gtable,data.frame(Model=combos3,
                                  AIC=rep(NA,length(combos3)),
                                  Deviance=rep(NA,length(combos3)),
                                  mod.idx=nrow(gtable)+1:length(combos3)))
  extra<- length(g)
  for(i in 1:length(combos3)){
    if(grepl("SpeakerType",combos3[i]) | grepl("VenueType",combos3[i])){
    g[[i+extra]] <- glm(paste("count", "~0+", combos3[i]),
                        family="poisson", data=dat)
    } else {
      g[[i+extra]] <- glm(paste("count", "~", combos3[i]),
                          family="poisson", data=dat)
    }
    gtable$AIC[i+extra] <- g[[i+extra]]$aic
    gtable$Deviance[i+extra] <-g[[i+extra]]$deviance
  }
  gtable <- gtable[order(gtable$AIC,decreasing=F),]
  gtable$dAIC <- gtable$AIC-gtable$AIC[1]
  gtable$Weight <- exp(-0.5 * gtable$dAIC)
  
  list(AICtab = gtable, mods = g)
}