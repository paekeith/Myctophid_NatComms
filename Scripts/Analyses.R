# this script runs all of the analyses and produces all of the figures that are part of
# the manuscript. Last updated 16/12/2023

#loading packages
library(ggplot2)
library(tidyverse)
library(nlme)
library(MuMIn)
library(RColorBrewer)
library(dplyr)
library(effects)
library(ape)
library(sp)
library(gstat)

# library(car)
# library(stats)
# library(nlmeU)

## FISH DATA -----------------------------------------------------------------------------
data <- read.csv("Data/Processed/FINAL_FISH_DATA.csv")

summary(data)

#generating a weighted average prey mass in the stomach of each predator
data_weighted_average <- data %>%
  group_by(ID_unique) %>%
  dplyr::summarise(fish_mass=mean(weight_discovery_g),fish_mass_log10=mean(fish_mass_log10),SST=mean(SST),CHL=mean(CHL),temp_depth=mean(temp_depth),
            weight_prey_mass=weighted.mean(Prey.Mass..g.,Count),Lat = mean(Lat),Lon = mean(Lon)) %>%
  ungroup()

windows()

#now adding in the species codes, event, cruise again
data_weighted_average$Species <- NA
data_weighted_average$Cruise <- NA
data_weighted_average$Event <- NA

for(i in 1:nrow(data_weighted_average)){
  split <- strsplit(data_weighted_average$ID_unique[i],split="_")
  data_weighted_average$Cruise[i] <- unlist(split)[1]
  data_weighted_average$Event[i] <- unlist(split)[2]
  data_weighted_average$Species[i] <- unlist(split)[3]
}

data_weighted_average$PPMR_log10 <- log10(data_weighted_average$fish_mass/data_weighted_average$weight_prey_mass)

### Model selection ####
#### PPMR with SST & CHL-a ####
# First fitting most complex fixed effects structure with different random effects
# Models which have been commented out do not converge
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

mod_gls <- gls(PPMR_log10~SST*CHL,data=data_weighted_average,method="REML") 
mod_1 <- lme(PPMR_log10~SST*CHL,random = list(~1|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_2 <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise),data=data_weighted_average,method = "REML",control=lmec) 
mod_3 <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_4 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_5 <- lme(PPMR_log10~SST*CHL,random = list(~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence
mod_6 <- lme(PPMR_log10~SST*CHL,random = list(~1+CHL+SST|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence

mod_7 <- lme(PPMR_log10~SST*CHL,random = list(~1+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_8 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_9 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)

# mod_10 <- lme(PPMR_log10~SST*CHL,random = list(~1+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_11 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_12 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)

mod_13 <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_14 <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_15 <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST+CHL|SST),data=data_weighted_average,method = "REML",control=lmec)

# mod_16 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_17 <- lme(PPMR_log10~SST*CHL,random = list(~1+CHL|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_18 <- lme(PPMR_log10~SST*CHL,random = list(~1+SST|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_19 <- lme(PPMR_log10~SST*CHL,random = list(~1+CHL|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(mod_gls,mod_1,mod_2,mod_3,mod_4,mod_5,mod_6,mod_7,mod_8,mod_9,mod_13,mod_14,mod_15,mod_19))
anova_frame[order(anova_frame$AIC,decreasing = T),]

#best model is 13
summary(mod_13)

#checking model fit and variance
windows(record=T)
plot(mod_13)
qqnorm(resid(mod_13,type = "pearson"))
qqline(resid(mod_13,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(mod_13,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(mod_13,type="pearson"))

#Now adding variance weighting
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

lme_ppmr_ident_sp <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                         weights = varIdent(form= ~ 1 | Species),
                         data=data_weighted_average,method = "REML",control=lmec) 

lme_ppmr_ident_Cruise <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                             weights = varIdent(form= ~ 1 | Cruise),
                             data=data_weighted_average,method = "REML",control=lmec) 

lme_ppmr_fix_SST <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varFixed(~SST),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_ppmr_fix_chl <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varFixed(~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#exp_var
lme_ppmr_exp_SST <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varExp(form=~SST),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_ppmr_exp_CHL <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varExp(form=~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#const_var
# lme_ppmr_const_SST <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
#                           weights = varConstPower(form=~SST),
#                           data=data_weighted_average,method = "REML",control=lmec)

lme_ppmr_const_CHL <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                          weights = varConstPower(form=~CHL),
                          data=data_weighted_average,method = "REML",control=lmec)
#comb_var
lme_ppmr_comb_sp_Cruise <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                               weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                               data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- anova(mod_13,lme_ppmr_ident_sp,lme_ppmr_ident_Cruise,lme_ppmr_fix_SST,lme_ppmr_fix_chl,
                     lme_ppmr_exp_SST,lme_ppmr_exp_CHL,lme_ppmr_const_CHL,lme_ppmr_comb_sp_Cruise)

anova_frame[order(anova_frame$AIC,decreasing = T),]

#best model has combined variance structure
summary(lme_ppmr_comb_sp_Cruise)
qqnorm(resid(lme_ppmr_comb_sp_Cruise,type = "pearson"))
qqline(resid(lme_ppmr_comb_sp_Cruise,type = "pearson")) 
plot(lme_ppmr_comb_sp_Cruise)
plot(factor(data_weighted_average$Species),resid(lme_ppmr_comb_sp_Cruise,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_ppmr_comb_sp_Cruise,type="pearson"))

#now comparing different fixed effects structures with ML:
lme_PPMR_SST_int_CHL <- lme(PPMR_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                            weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                            data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR_SST_plus_CHL <- lme(PPMR_log10~SST+CHL,random = list(~1|Cruise,~1+SST|Species),
                             weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                             data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR_CHL <- lme(PPMR_log10~CHL,random = list(~1|Cruise,~1+SST|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR_SST <- lme(PPMR_log10~SST,random = list(~1|Cruise,~1+SST|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR <- lme(PPMR_log10~1,random = list(~1|Cruise,~1+SST|Species),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                data=data_weighted_average,method = "ML",control=lmec)

anova_frame <- anova(lme_PPMR_SST_int_CHL,lme_PPMR_SST_plus_CHL,lme_PPMR_CHL,lme_PPMR_SST,lme_PPMR)
anova_frame[order(anova_frame$AIC,decreasing = T),]

qqnorm(resid(lme_PPMR_SST,type = "pearson"))
qqline(resid(lme_PPMR_SST,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_PPMR_SST,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_PPMR_SST,type="pearson"))

#checking autocorrelation
temp_data = data.frame(error = residuals(lme_PPMR_SST), x = data_weighted_average$Lon, y = data_weighted_average$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")
dists <- as.matrix(dist(cbind(data_weighted_average$Lon,data_weighted_average$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(lme_PPMR_SST,type = "normalized"),dists) 

#now fitting the final mode with REML
lme_PPMR_SST <- lme(PPMR_log10~SST,random = list(~1|Cruise,~1+SST|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "REML",control=lmec)

summary(lme_PPMR_SST)
r.squaredGLMM(lme_PPMR_SST)

#plotting the model partial residuals
model_partial <- data.frame(effect("SST", lme_PPMR_SST, xlevels=list(SST=seq(min(data_weighted_average$SST),max(data_weighted_average$SST),length=100))))
data_weighted_average$resids <- resid(lme_PPMR_SST) + summary(lme_PPMR_SST)$tTable[1] + summary(lme_PPMR_SST)$tTable[2]*data_weighted_average$SST

data_weighted_average_abline <- as.data.frame(fixef(lme_PPMR_SST)[1])
data_weighted_average_abline[,2] <- fixef(lme_PPMR_SST)[2]
names(data_weighted_average_abline) <- c("Intercept","Slope")
data_weighted_average_abline$xmin <- min(data_weighted_average$SST)
data_weighted_average_abline$xmax <- max(data_weighted_average$SST)
data_weighted_average_abline$ymax <- (data_weighted_average_abline$Intercept)+((data_weighted_average_abline$xmin)*data_weighted_average_abline$Slope) #adjusting intercept to be correct for SST of -0.16
data_weighted_average_abline$ymin <- (data_weighted_average_abline$Intercept+(data_weighted_average_abline$xmax*data_weighted_average_abline$Slope)) #calculating the end point for the lines

ggplot(data_weighted_average,aes(x=SST,y=resids))+
  geom_point(size=3,alpha=0.3,pch=21,fill="grey",colour="black")+geom_segment(data=data_weighted_average_abline,aes(x=xmin,xend=xmax,y=ymax,yend=ymin),inherit.aes = FALSE,lwd=1.5,colour="black")+
  geom_ribbon(data=model_partial,aes(x=model_partial$SST,ymin=model_partial$lower,ymax=model_partial$upper),inherit.aes = FALSE,fill="darkgrey",alpha=0.6)+
  theme_bw()+labs(y="PPMR",x= "SST (°C)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),
                                                axis.title.y = element_text(size=30,vjust=1.5),legend.position="bottom",legend.text = element_text(size=15),legend.title = element_blank())+
  theme(legend.key = element_blank())+scale_y_continuous(breaks=seq(0,6,1),limits=c(0,5))+scale_x_continuous(breaks=seq(-2,8,2),limits=c(-2,8))+guides(colour=guide_legend(nrow=1,byrow=TRUE))

#### Fish mass with SST & CHL-a ####
# Models which have been commented out do not converge
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

mod_gls <- gls(fish_mass_log10~SST*CHL,data=data_weighted_average,method="REML") 
mod_1 <- lme(fish_mass_log10~SST*CHL,random = list(~1|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_2 <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise),data=data_weighted_average,method = "REML",control=lmec) 
mod_3 <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_4 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec) 
# mod_5 <- lme(fish_mass_log10~SST*CHL,random = list(~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_6 <- lme(fish_mass_log10~SST*CHL,random = list(~1+CHL+SST|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_7 <- lme(fish_mass_log10~SST*CHL,random = list(~1+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_8 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_9 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)

mod_10 <- lme(fish_mass_log10~SST*CHL,random = list(~1+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_11 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_12 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)

mod_13 <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_14 <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_15 <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST+CHL|SST),data=data_weighted_average,method = "REML",control=lmec)

# mod_16 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_17 <- lme(fish_mass_log10~SST*CHL,random = list(~1+CHL|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_18 <- lme(fish_mass_log10~SST*CHL,random = list(~1+SST|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_19 <- lme(fish_mass_log10~SST*CHL,random = list(~1+CHL|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(mod_gls,mod_1,mod_2,mod_3,mod_4,mod_6,mod_7,mod_8,mod_9,mod_10,mod_11,mod_12,mod_13,
                                   mod_14,mod_15,mod_17,mod_18,mod_19))

anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(mod_13)
qqnorm(resid(mod_13,type = "pearson"))
qqline(resid(mod_13,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(mod_13,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(mod_13,type="pearson"))

# adding variance structures
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)
lme_mass_ident_sp <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                         weights = varIdent(form= ~ 1 | Species),
                         data=data_weighted_average,method = "REML",control=lmec) 

lme_mass_ident_Cruise <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                             weights = varIdent(form= ~ 1 | Cruise),
                             data=data_weighted_average,method = "REML",control=lmec) 

lme_mass_fix_SST <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varFixed(~SST),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_mass_fix_chl <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varFixed(~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#exp_var
lme_mass_exp_SST <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varExp(form=~SST),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_mass_exp_CHL <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varExp(form=~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#const_var
# lme_mass_const_SST <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
#                           weights = varConstPower(form=~SST),
#                           data=data_weighted_average,method = "REML",control=lmec)

lme_mass_const_CHL <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                          weights = varConstPower(form=~CHL),
                          data=data_weighted_average,method = "REML",control=lmec)
#comb_var
lme_mass_comb_sp_Cruise <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                               weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                               data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- anova(mod_13,lme_mass_ident_sp,lme_mass_ident_Cruise,lme_mass_fix_SST,lme_mass_fix_chl,lme_mass_exp_SST,
                     lme_mass_exp_CHL,lme_mass_const_CHL,lme_mass_comb_sp_Cruise)
anova_frame[order(anova_frame$AIC,decreasing = T),]

#combined is best
plot(lme_mass_comb_sp_Cruise)
qqnorm(resid(lme_mass_comb_sp_Cruise,type = "pearson"))
qqline(resid(lme_mass_comb_sp_Cruise,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_mass_comb_sp_Cruise,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_mass_comb_sp_Cruise,type="pearson"))

#now comparing different fixed effects structures with ML:
lme_mass_SST_int_CHL <- lme(fish_mass_log10~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                            weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                            data=data_weighted_average,method = "ML",control=lmec)

lme_mass_SST_plus_CHL <- lme(fish_mass_log10~SST+CHL,random = list(~1|Cruise,~1+SST|Species),
                             weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                             data=data_weighted_average,method = "ML",control=lmec)

lme_mass_CHL <- lme(fish_mass_log10~CHL,random = list(~1|Cruise,~1+SST|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_mass_SST <- lme(fish_mass_log10~SST,random = list(~1|Cruise,~1+SST|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_mass <- lme(fish_mass_log10~1,random = list(~1|Cruise,~1+SST|Species),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                data=data_weighted_average,method = "ML",control=lmec)

anova_frame <- anova(lme_mass_SST_int_CHL,lme_mass_SST_plus_CHL,lme_mass_CHL,lme_mass_SST,lme_mass)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_mass_SST)
qqnorm(resid(lme_mass_SST,type = "pearson"))
qqline(resid(lme_mass_SST,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_mass_SST,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_mass_SST,type="pearson"))

#checking autocorrelation
temp_data = data.frame(error = residuals(lme_mass_SST), x = data_weighted_average$Lon, y = data_weighted_average$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")
dists <- as.matrix(dist(cbind(data_weighted_average$Lon,data_weighted_average$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(lme_mass_SST,type = "normalized"),dists) # suggests no autocorrelation in residuals

#now fitting model with REML
lme_mass_SST <- lme(fish_mass_log10~SST,random = list(~1|Cruise,~1+SST|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "REML",control=lmec)

summary(lme_mass_SST)
r.squaredGLMM(lme_mass_SST)

#plotting the model partial residuals
model_partial <- as.data.frame(effect("SST", lme_mass_SST, xlevels=list(SST=seq(min(data_weighted_average$SST),max(data_weighted_average$SST),length=100))))
data_weighted_average$resids <- resid(lme_mass_SST) + summary(lme_mass_SST)$tTable[1] + summary(lme_mass_SST)$tTable[2]*data_weighted_average$SST

data_abline <- as.data.frame(fixef(lme_mass_SST)[1])
data_abline[,2] <- fixef(lme_mass_SST)[2]
names(data_abline) <- c("Intercept","Slope")
data_abline$xmin <- min(data_weighted_average$SST)
data_abline$xmax <- max(data_weighted_average$SST)
data_abline$ymax <- (data_abline$Intercept)+((data_abline$xmin)*data_abline$Slope) #adjusting intercept to be correct for SST of -0.16
data_abline$ymin <- (data_abline$Intercept+(data_abline$xmax*data_abline$Slope)) #calculating the end point for the lines

ggplot(data_weighted_average,aes(x=SST,y=resids))+
  geom_point(size=3,alpha=0.3,pch=21,fill="grey",colour="black")+geom_segment(data=data_abline,aes(x=xmin,xend=xmax,y=ymax,yend=ymin),inherit.aes = FALSE,lwd=1.5,colour="black")+
  geom_ribbon(data=model_partial,aes(x=model_partial$SST,ymin=model_partial$lower,ymax=model_partial$upper),inherit.aes = FALSE,fill="darkgrey",alpha=0.6)+
  theme_bw()+labs(y="Predator mass",x= "SST (°C)")+theme(axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),axis.title.y = element_text(size=30,vjust=1.5),legend.position="bottom",legend.text = element_text(size=15),legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_blank())+scale_y_continuous(breaks=seq(-1,1,1),limits=c())+scale_x_continuous(breaks=seq(-2,8,2),limits=c(-2,8))+guides(colour=guide_legend(nrow=1,byrow=TRUE))

#### Prey mass with SST & CHL-a ####
# Models which have been commented out do not converge
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

mod_gls <- gls(log10(weight_prey_mass)~SST*CHL,data=data_weighted_average,method="REML") 
mod_1 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_2 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise),data=data_weighted_average,method = "REML",control=lmec) 
mod_3 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_4 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec) 
# mod_5 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence
# mod_6 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+CHL+SST|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence

mod_7 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_8 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST|Cruise),data=data_weighted_average,method = "REML",control=lmec)
# mod_9 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)

# mod_10 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_11 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_12 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)

mod_13 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_14 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_15 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST+CHL|SST),data=data_weighted_average,method = "REML",control=lmec)

# mod_16 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_17 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+CHL|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_18 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+SST|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_19 <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1+CHL|Cruise,~1+SST|Species),data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(mod_gls,mod_1,mod_2,mod_3,mod_4,mod_7,mod_8,mod_13,mod_14,mod_15))
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(mod_13)
qqnorm(resid(mod_13,type = "pearson"))
qqline(resid(mod_13,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(mod_13,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(mod_13,type="pearson"))

# adding variance weighting
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)
lme_ident_sp <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                    weights = varIdent(form= ~ 1 | Species),
                    data=data_weighted_average,method = "REML",control=lmec) 

lme_ident_Cruise <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                        weights = varIdent(form= ~ 1 | Cruise),
                        data=data_weighted_average,method = "REML",control=lmec) 

lme_fix_SST <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                   weights = varFixed(~SST),
                   data=data_weighted_average,method = "REML",control=lmec)

lme_fix_CHL <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                   weights = varFixed(~CHL),
                   data=data_weighted_average,method = "REML",control=lmec)

#exp_var
lme_exp_SST <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                   weights = varExp(form=~SST),
                   data=data_weighted_average,method = "REML",control=lmec)

lme_exp_CHL <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                   weights = varExp(form=~CHL),
                   data=data_weighted_average,method = "REML",control=lmec)
#const_var
lme_const_SST <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                     weights = varConstPower(form=~SST),
                     data=data_weighted_average,method = "REML",control=lmec)

lme_const_CHL <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                     weights = varConstPower(form=~CHL),
                     data=data_weighted_average,method = "REML",control=lmec)
#comb_var
lme_comb_sp_Cruise <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                          weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                          data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- anova(mod_13,lme_ident_sp,lme_ident_Cruise,lme_fix_SST,lme_fix_CHL,lme_exp_SST,lme_exp_CHL,lme_const_SST,lme_const_CHL,lme_comb_sp_Cruise)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_comb_sp_Cruise)
qqnorm(resid(lme_comb_sp_Cruise,type = "pearson"))
qqline(resid(lme_comb_sp_Cruise,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_comb_sp_Cruise,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_comb_sp_Cruise,type="pearson"))

#now comparing different fixed effects structures with ML:
lme_prey_mass_SST_int_CHL <- lme(log10(weight_prey_mass)~SST*CHL,random = list(~1|Cruise,~1+SST|Species),
                                 weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                                 data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass_SST_plus_CHL <- lme(log10(weight_prey_mass)~SST+CHL,random = list(~1|Cruise,~1+SST|Species),
                                  weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                                  data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass_CHL <- lme(log10(weight_prey_mass)~CHL,random = list(~1|Cruise,~1+SST|Species),
                         weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                         data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass_SST <- lme(log10(weight_prey_mass)~SST,random = list(~1|Cruise,~1+SST|Species),
                         weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                         data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass <- lme(log10(weight_prey_mass)~1,random = list(~1|Cruise,~1+SST|Species),
                     weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                     data=data_weighted_average,method = "ML",control=lmec)

anova_frame <- anova(lme_prey_mass_SST_int_CHL,lme_prey_mass_SST_plus_CHL,lme_prey_mass_CHL,lme_prey_mass_SST,lme_prey_mass)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_prey_mass)
qqnorm(resid(lme_prey_mass,type = "pearson"))
qqline(resid(lme_prey_mass,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_prey_mass,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_prey_mass,type="pearson"))

#checking autocorrelation
temp_data = data.frame(error = residuals(lme_prey_mass), x = data_weighted_average$Lon, y = data_weighted_average$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")
dists <- as.matrix(dist(cbind(data_weighted_average$Lon,data_weighted_average$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(lme_prey_mass,type = "normalized"),dists) # suggests no autocorrelation in residuals

#Fitting model with REML
lme_prey_mass <- lme(log10(weight_prey_mass)~1,random = list(~1|Cruise,~1+SST|Species),
                     weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                     data=data_weighted_average,method = "REML",control=lmec)

summary(lme_prey_mass)
r.squaredGLMM(lme_prey_mass)

#plotting the data
ggplot(data_weighted_average,aes(x=SST,y=log10(weight_prey_mass)))+
  geom_point(size=3,alpha=0.3,colour="black",fill="grey",pch=21)+
  labs(y="Mean prey mass",x= "SST (°C)")+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),axis.title.y = element_text(size=30,vjust=1.5),legend.position="bottom",legend.text = element_text(size=15),legend.title = element_blank())+
  theme(legend.key = element_blank())+scale_y_continuous(breaks=seq(-6,6,1),limits=c())+scale_x_continuous(breaks=seq(-2,8,2),limits=c(-2,8))+guides(colour=guide_legend(nrow=1,byrow=TRUE))


#### PPMR with temp_depth & CHL-a ####
# First fitting most complex fixed effects structure with different random effects
# Models which have been commented out do not converge
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

mod_gls <- gls(PPMR_log10~temp_depth*CHL,data=data_weighted_average,method="REML") 
mod_1 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_2 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise),data=data_weighted_average,method = "REML",control=lmec) 
mod_3 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_4 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_5 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence
mod_6 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+CHL+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence

mod_7 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_8 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise),data=data_weighted_average,method = "REML",control=lmec)
# mod_9 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)

# mod_10 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_11 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_12 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)

mod_13 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_14 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_15 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth+CHL|temp_depth),data=data_weighted_average,method = "REML",control=lmec)

# mod_16 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_17 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+CHL|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_18 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_19 <- lme(PPMR_log10~temp_depth*CHL,random = list(~1+CHL|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(mod_gls,mod_1,mod_2,mod_3,mod_4,mod_5,mod_6,mod_7,mod_8,mod_13,mod_14,mod_15,mod_19))
anova_frame[order(anova_frame$AIC,decreasing = T),]

#checking model fit and variance
plot(mod_13)
qqnorm(resid(mod_13,type = "pearson"))
qqline(resid(mod_13,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(mod_13,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(mod_13,type="pearson"))

#Now adding variance weighting
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

lme_ppmr_ident_sp <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                         weights = varIdent(form= ~ 1 | Species),
                         data=data_weighted_average,method = "REML",control=lmec) 

lme_ppmr_ident_Cruise <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                             weights = varIdent(form= ~ 1 | Cruise),
                             data=data_weighted_average,method = "REML",control=lmec) 

lme_ppmr_fix_temp_depth <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varFixed(~temp_depth),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_ppmr_fix_chl <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varFixed(~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#exp_var
lme_ppmr_exp_temp_depth <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varExp(form=~temp_depth),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_ppmr_exp_CHL <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varExp(form=~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#const_var
lme_ppmr_const_temp_depth <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                          weights = varConstPower(form=~temp_depth),
                          data=data_weighted_average,method = "REML",control=lmec)

lme_ppmr_const_CHL <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                          weights = varConstPower(form=~CHL),
                          data=data_weighted_average,method = "REML",control=lmec)
#comb_var
lme_ppmr_comb_sp_Cruise <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                               weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                               data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- anova(mod_13,lme_ppmr_ident_sp,lme_ppmr_ident_Cruise,lme_ppmr_fix_temp_depth,lme_ppmr_fix_chl,
                     lme_ppmr_exp_temp_depth,lme_ppmr_exp_CHL,lme_ppmr_const_temp_depth,lme_ppmr_const_CHL,lme_ppmr_comb_sp_Cruise)

anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_ppmr_comb_sp_Cruise)
qqnorm(resid(lme_ppmr_comb_sp_Cruise,type = "pearson"))
qqline(resid(lme_ppmr_comb_sp_Cruise,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_ppmr_comb_sp_Cruise,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_ppmr_comb_sp_Cruise,type="pearson"))

#now comparing different fixed effects structures with ML:
lme_PPMR_temp_depth_int_CHL <- lme(PPMR_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                            weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                            data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR_temp_depth_plus_CHL <- lme(PPMR_log10~temp_depth+CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                             weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                             data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR_CHL <- lme(PPMR_log10~CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR_temp_depth <- lme(PPMR_log10~temp_depth,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_PPMR <- lme(PPMR_log10~1,random = list(~1|Cruise,~1+temp_depth|Species),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                data=data_weighted_average,method = "ML",control=lmec)

anova_frame <- anova(lme_PPMR_temp_depth_int_CHL,lme_PPMR_temp_depth_plus_CHL,lme_PPMR_CHL,lme_PPMR_temp_depth,lme_PPMR)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_PPMR_temp_depth)
qqnorm(resid(lme_PPMR_temp_depth,type = "pearson"))
qqline(resid(lme_PPMR_temp_depth,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_PPMR_temp_depth,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_PPMR_temp_depth,type="pearson"))

#checking autocorrelation
temp_data = data.frame(error = residuals(lme_PPMR_temp_depth), x = data_weighted_average$Lon, y = data_weighted_average$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")
dists <- as.matrix(dist(cbind(data_weighted_average$Lon,data_weighted_average$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(lme_PPMR_temp_depth,type = "normalized"),dists) 

#now fitting the final model with REML
lme_PPMR_temp_depth <- lme(PPMR_log10~temp_depth,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "REML",control=lmec)

summary(lme_PPMR_temp_depth)
r.squaredGLMM(lme_PPMR_temp_depth)

#plotting the model partial residuals
model_partial <- data.frame(effect("temp_depth", lme_PPMR_temp_depth, xlevels=list(temp_depth=seq(min(data_weighted_average$temp_depth),max(data_weighted_average$temp_depth),length=100))))
data_weighted_average$resids <- resid(lme_PPMR_temp_depth) + summary(lme_PPMR_temp_depth)$tTable[1] + summary(lme_PPMR_temp_depth)$tTable[2]*data_weighted_average$temp_depth

data_weighted_average_abline <- as.data.frame(fixef(lme_PPMR_temp_depth)[1])
data_weighted_average_abline[,2] <- fixef(lme_PPMR_temp_depth)[2]
names(data_weighted_average_abline) <- c("Intercept","Slope")
data_weighted_average_abline$xmin <- min(data_weighted_average$temp_depth)
data_weighted_average_abline$xmax <- max(data_weighted_average$temp_depth)
data_weighted_average_abline$ymax <- (data_weighted_average_abline$Intercept)+((data_weighted_average_abline$xmin)*data_weighted_average_abline$Slope) #adjusting intercept to be correct for temp_depth of -0.16
data_weighted_average_abline$ymin <- (data_weighted_average_abline$Intercept+(data_weighted_average_abline$xmax*data_weighted_average_abline$Slope)) #calculating the end point for the lines

ggplot(data_weighted_average,aes(x=temp_depth,y=resids))+
  geom_point(size=3,alpha=0.3,pch=21,fill="grey",colour="black")+geom_segment(data=data_weighted_average_abline,aes(x=xmin,xend=xmax,y=ymax,yend=ymin),inherit.aes = FALSE,lwd=1.5,colour="black")+
  geom_ribbon(data=model_partial,aes(x=model_partial$temp_depth,ymin=model_partial$lower,ymax=model_partial$upper),inherit.aes = FALSE,fill="darkgrey",alpha=0.6)+
  theme_bw()+labs(y="PPMR",x= "temp_depth (°C)")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),
                                                axis.title.y = element_text(size=30,vjust=1.5),legend.position="bottom",legend.text = element_text(size=15),legend.title = element_blank())+
  theme(legend.key = element_blank())+scale_y_continuous(breaks=seq(0,6,1),limits=c(0,5))+scale_x_continuous(breaks=seq(-2,8,2),limits=c(-0.2,2.2))+guides(colour=guide_legend(nrow=1,byrow=TRUE))

#### Fish mass with temp_depth & CHL-a ####
# Models which have been commented out do not converge
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

mod_gls <- gls(fish_mass_log10~temp_depth*CHL,data=data_weighted_average,method="REML") 
mod_1 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_2 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise),data=data_weighted_average,method = "REML",control=lmec) 
mod_3 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_4 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec) 
# mod_5 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence
mod_6 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+CHL+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence

mod_7 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)
# mod_8 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise),data=data_weighted_average,method = "REML",control=lmec)
# mod_9 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)

# mod_10 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_11 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_12 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)

mod_13 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_14 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_15 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth+CHL|temp_depth),data=data_weighted_average,method = "REML",control=lmec)

mod_16 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_17 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+CHL|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_18 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_19 <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1+CHL|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(mod_gls,mod_1,mod_2,mod_3,mod_4,mod_6,mod_7,mod_13,
                                   mod_14,mod_16,mod_18))

anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(mod_13)
qqnorm(resid(mod_13,type = "pearson"))
qqline(resid(mod_13,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(mod_13,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(mod_13,type="pearson"))

# adding variance structures
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)
lme_mass_ident_sp <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                         weights = varIdent(form= ~ 1 | Species),
                         data=data_weighted_average,method = "REML",control=lmec) 

lme_mass_ident_Cruise <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                             weights = varIdent(form= ~ 1 | Cruise),
                             data=data_weighted_average,method = "REML",control=lmec) 

lme_mass_fix_temp_depth <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varFixed(~temp_depth),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_mass_fix_chl <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varFixed(~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#exp_var
lme_mass_exp_temp_depth <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varExp(form=~temp_depth),
                        data=data_weighted_average,method = "REML",control=lmec)

lme_mass_exp_CHL <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varExp(form=~CHL),
                        data=data_weighted_average,method = "REML",control=lmec)
#const_var
# lme_mass_const_temp_depth <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
#                           weights = varConstPower(form=~temp_depth),
#                           data=data_weighted_average,method = "REML",control=lmec)

lme_mass_const_CHL <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                          weights = varConstPower(form=~CHL),
                          data=data_weighted_average,method = "REML",control=lmec)
#comb_var
lme_mass_comb_sp_Cruise <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                               weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                               data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- anova(mod_13,lme_mass_ident_sp,lme_mass_ident_Cruise,lme_mass_fix_temp_depth,lme_mass_fix_chl,lme_mass_exp_temp_depth,
                     lme_mass_exp_CHL,lme_mass_const_CHL,lme_mass_comb_sp_Cruise)
anova_frame[order(anova_frame$AIC,decreasing = T),]

#combined is best
plot(lme_mass_comb_sp_Cruise)
qqnorm(resid(lme_mass_comb_sp_Cruise,type = "pearson"))
qqline(resid(lme_mass_comb_sp_Cruise,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_mass_comb_sp_Cruise,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_mass_comb_sp_Cruise,type="pearson"))

#now comparing different fixed effects structures with ML:
lme_mass_temp_depth_int_CHL <- lme(fish_mass_log10~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                            weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                            data=data_weighted_average,method = "ML",control=lmec)

lme_mass_temp_depth_plus_CHL <- lme(fish_mass_log10~temp_depth+CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                             weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                             data=data_weighted_average,method = "ML",control=lmec)

lme_mass_CHL <- lme(fish_mass_log10~CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_mass_temp_depth <- lme(fish_mass_log10~temp_depth,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "ML",control=lmec)

lme_mass <- lme(fish_mass_log10~1,random = list(~1|Cruise,~1+temp_depth|Species),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                data=data_weighted_average,method = "ML",control=lmec)

anova_frame <- anova(lme_mass_temp_depth_int_CHL,lme_mass_temp_depth_plus_CHL,lme_mass_CHL,lme_mass_temp_depth,lme_mass)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_mass_temp_depth)
qqnorm(resid(lme_mass_temp_depth,type = "pearson"))
qqline(resid(lme_mass_temp_depth,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_mass_temp_depth,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_mass_temp_depth,type="pearson"))

#checking autocorrelation
temp_data = data.frame(error = residuals(lme_mass_temp_depth), x = data_weighted_average$Lon, y = data_weighted_average$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")
dists <- as.matrix(dist(cbind(data_weighted_average$Lon,data_weighted_average$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(lme_mass_temp_depth,type = "normalized"),dists) # suggests no autocorrelation in residuals

#now fitting model with REML
lme_mass_temp_depth <- lme(fish_mass_log10~temp_depth,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                    data=data_weighted_average,method = "REML",control=lmec)

summary(lme_mass_temp_depth)
r.squaredGLMM(lme_mass_temp_depth)

#plotting the model partial residuals
model_partial <- as.data.frame(effect("temp_depth", lme_mass_temp_depth, xlevels=list(temp_depth=seq(min(data_weighted_average$temp_depth),max(data_weighted_average$temp_depth),length=100))))
data_weighted_average$resids <- resid(lme_mass_temp_depth) + summary(lme_mass_temp_depth)$tTable[1] + summary(lme_mass_temp_depth)$tTable[2]*data_weighted_average$temp_depth

data_abline <- as.data.frame(fixef(lme_mass_temp_depth)[1])
data_abline[,2] <- fixef(lme_mass_temp_depth)[2]
names(data_abline) <- c("Intercept","Slope")
data_abline$xmin <- min(data_weighted_average$temp_depth)
data_abline$xmax <- max(data_weighted_average$temp_depth)
data_abline$ymax <- (data_abline$Intercept)+((data_abline$xmin)*data_abline$Slope) #adjusting intercept to be correct for temp_depth of -0.16
data_abline$ymin <- (data_abline$Intercept+(data_abline$xmax*data_abline$Slope)) #calculating the end point for the lines

ggplot(data_weighted_average,aes(x=temp_depth,y=resids))+
  geom_point(size=3,alpha=0.3,pch=21,fill="grey",colour="black")+geom_segment(data=data_abline,aes(x=xmin,xend=xmax,y=ymax,yend=ymin),inherit.aes = FALSE,lwd=1.5,colour="black")+
  geom_ribbon(data=model_partial,aes(x=model_partial$temp_depth,ymin=model_partial$lower,ymax=model_partial$upper),inherit.aes = FALSE,fill="darkgrey",alpha=0.6)+
  theme_bw()+labs(y="Predator mass",x= "temp_depth (°C)")+theme(axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),axis.title.y = element_text(size=30,vjust=1.5),legend.position="bottom",legend.text = element_text(size=15),legend.title = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_blank())+scale_y_continuous(breaks=seq(-1,1,1),limits=c())+scale_x_continuous(breaks=seq(-2,8,2),limits=c(-0.2,3))+guides(colour=guide_legend(nrow=1,byrow=TRUE))

#### Prey mass with temp_depth & CHL-a ####
# Models which have been commented out do not converge
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

mod_gls <- gls(log10(weight_prey_mass)~temp_depth*CHL,data=data_weighted_average,method="REML") 
mod_1 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_2 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise),data=data_weighted_average,method = "REML",control=lmec) 
mod_3 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec) 

mod_4 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec) 
mod_5 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence
mod_6 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+CHL+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec) #no convergence

mod_7 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_8 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth|Cruise),data=data_weighted_average,method = "REML",control=lmec)
mod_9 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth+CHL|Cruise),data=data_weighted_average,method = "REML",control=lmec)

# mod_10 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_11 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_12 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth+CHL|Cruise,~1|Species),data=data_weighted_average,method = "REML",control=lmec)

mod_13 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_14 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
mod_15 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth+CHL|temp_depth),data=data_weighted_average,method = "REML",control=lmec)

# mod_16 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_17 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+CHL|Cruise,~1+CHL|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_18 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+temp_depth|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)
# mod_19 <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1+CHL|Cruise,~1+temp_depth|Species),data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(mod_gls,mod_1,mod_2,mod_3,mod_4,mod_5,mod_6,mod_7,mod_8,mod_9,mod_13,mod_14,mod_15))
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(mod_13)
qqnorm(resid(mod_13,type = "pearson"))
qqline(resid(mod_13,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(mod_13,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(mod_13,type="pearson"))

# adding variance weighting
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)
lme_ident_sp <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                    weights = varIdent(form= ~ 1 | Species),
                    data=data_weighted_average,method = "REML",control=lmec) 

lme_ident_Cruise <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                        weights = varIdent(form= ~ 1 | Cruise),
                        data=data_weighted_average,method = "REML",control=lmec) 

lme_fix_temp_depth <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                   weights = varFixed(~temp_depth),
                   data=data_weighted_average,method = "REML",control=lmec)

lme_fix_CHL <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                   weights = varFixed(~CHL),
                   data=data_weighted_average,method = "REML",control=lmec)

#exp_var
lme_exp_temp_depth <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                   weights = varExp(form=~temp_depth),
                   data=data_weighted_average,method = "REML",control=lmec)

lme_exp_CHL <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                   weights = varExp(form=~CHL),
                   data=data_weighted_average,method = "REML",control=lmec)
#const_var
lme_const_temp_depth <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                     weights = varConstPower(form=~temp_depth),
                     data=data_weighted_average,method = "REML",control=lmec)

lme_const_CHL <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                     weights = varConstPower(form=~CHL),
                     data=data_weighted_average,method = "REML",control=lmec)
#comb_var
lme_comb_sp_Cruise <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                          weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                          data=data_weighted_average,method = "REML",control=lmec)

anova_frame <- anova(mod_13,lme_ident_sp,lme_ident_Cruise,lme_fix_temp_depth,lme_fix_CHL,lme_exp_temp_depth,lme_exp_CHL,lme_const_temp_depth,lme_const_CHL,lme_comb_sp_Cruise)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_comb_sp_Cruise)
qqnorm(resid(lme_comb_sp_Cruise,type = "pearson"))
qqline(resid(lme_comb_sp_Cruise,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_comb_sp_Cruise,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_comb_sp_Cruise,type="pearson"))

#now comparing different fixed effects structures with ML:
lme_prey_mass_temp_depth_int_CHL <- lme(log10(weight_prey_mass)~temp_depth*CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                                 weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                                 data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass_temp_depth_plus_CHL <- lme(log10(weight_prey_mass)~temp_depth+CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                                  weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                                  data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass_CHL <- lme(log10(weight_prey_mass)~CHL,random = list(~1|Cruise,~1+temp_depth|Species),
                         weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                         data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass_temp_depth <- lme(log10(weight_prey_mass)~temp_depth,random = list(~1|Cruise,~1+temp_depth|Species),
                         weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                         data=data_weighted_average,method = "ML",control=lmec)

lme_prey_mass <- lme(log10(weight_prey_mass)~1,random = list(~1|Cruise,~1+temp_depth|Species),
                     weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                     data=data_weighted_average,method = "ML",control=lmec)

anova_frame <- anova(lme_prey_mass_temp_depth_int_CHL,lme_prey_mass_temp_depth_plus_CHL,lme_prey_mass_CHL,lme_prey_mass_temp_depth,lme_prey_mass)
anova_frame[order(anova_frame$AIC,decreasing = T),]

plot(lme_prey_mass)
qqnorm(resid(lme_prey_mass,type = "pearson"))
qqline(resid(lme_prey_mass,type = "pearson")) 
plot(factor(data_weighted_average$Species),resid(lme_prey_mass,type="pearson"))
plot(factor(data_weighted_average$Cruise),resid(lme_prey_mass,type="pearson"))

#checking autocorrelation
temp_data = data.frame(error = residuals(lme_prey_mass), x = data_weighted_average$Lon, y = data_weighted_average$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")
dists <- as.matrix(dist(cbind(data_weighted_average$Lon,data_weighted_average$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(lme_prey_mass,type = "normalized"),dists) # suggests no autocorrelation in residuals

r.squaredGLMM(lme_prey_mass)

#Fitting model with REML
lme_prey_mass <- lme(log10(weight_prey_mass)~1,random = list(~1|Cruise,~1+temp_depth|Species),
                     weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | Species)),
                     data=data_weighted_average,method = "REML",control=lmec)

summary(lme_prey_mass)
r.squaredGLMM(lme_prey_mass)

#plotting the data
ggplot(data_weighted_average,aes(x=temp_depth,y=log10(weight_prey_mass)))+
  geom_point(size=3,alpha=0.3,colour="black",fill="grey",pch=21)+
  labs(y="Mean prey mass",x= "temp_depth (°C)")+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(size=26,colour="black"),axis.title.x = element_text(size=30),axis.title.y = element_text(size=30,vjust=1.5),legend.position="bottom",legend.text = element_text(size=15),legend.title = element_blank())+
  theme(legend.key = element_blank())+scale_y_continuous(breaks=seq(-6,6,1),limits=c())+scale_x_continuous(breaks=seq(-2,8,2),limits=c(-2,8))+guides(colour=guide_legend(nrow=1,byrow=TRUE))


## DIETARY SIZE SELECTIVITY ##-----------------------------------------------------------

## firstly I need to calculate the size selectivity
rm(list=ls())
dev.off()
library("Hmisc")
library("dplyr")
library(dplyr)
library(tidyverse)
windows(record=T)
# increase rule of thumb kernel bandwidth by scaling factor to smooth out density
bw.scale.zoop = 3.0 
# bw.scale.fish = 2.0

# kernel function over sequence x, the observation xi and bandwidth h
kern = function(x,xi,h){
  z=(x-xi)/h
  # Gaussian kernel:
  out=dnorm(z, mean=0, sd=1)
  return(out)
}

data.zoop <- read.csv("Data/Processed/PLANKTON_DATA.csv")
data.fish <- read.csv("Data/Processed/FINAL_FISH_DATA.csv")
#excluding the rows which do not have a location group
data.fish <- subset(data.fish,data.fish$location_group!="NA")

n.env = max(data.zoop$location_group)
# delete extremely small values well outside range of fish diet
data.zoop <- subset(data.zoop,log10_mass>=(-4.5607)&log10_mass<=(-0.0195))
data.fish$fish_mass_log10 
data.fish$Prey.Mass..g. <- signif(data.fish$Prey.Mass..g.,3)

#subsetting the fish data to only retain relevant columns
names(data.fish)
data.fish <- data.fish%>%
  dplyr::select("Cruise","Event","Year","Species.Code","Prey","Count","Prey.Mass..g.","weight_discovery_g","Lat","Lon","SST","location_group","fish_mass_log10")
names(data.fish) <- c("Cruise","Event","Year","Fish_Species_Code","Prey","Prey_Count","Mean_Prey_Mass__g_","Fish_Mass_g_","Lat","Lon","SST","location_group","fish_mass_log10")

## PREFERENCE ESTIMATOR 
# bin fish by 0.05 and aggregate statistical fish
breaks <- round(seq(-1,1.7,by=0.05),3)
# specify interval/bin labels
tags <- seq(-0.975,1.675,by=0.05)
# bucketing values into bins
group_tags <- cut(data.fish$fish_mass_log10, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)

data.fish$log_10_fish_class_bin <- round(as.numeric(as.character(group_tags)),3)

data.fish$Statistical_Fish_ID_mass = paste0(data.fish$Fish_Species_Code,"_",
                                            data.fish$location_group,"_",
                                            data.fish$log_10_fish_class_bin,"g")
IDs = unique(data.fish$Statistical_Fish_ID_mass)
n.IDs = length(IDs)

# count total prey items per statistical fish
n.prey = rep(NA,n.IDs)
data.fish$n_prey_items = NA
for(i in 1:n.IDs){
  df = subset(data.fish,Statistical_Fish_ID_mass==IDs[i])
  n.prey[i] = nrow(df)
  data.fish$n_prey_items[data.fish$Statistical_Fish_ID_mass==IDs[i]] = n.prey[i]
}

# remove statistical fish with less than 6 differently sized prey items.
# this provides the right balance between sample size and power
data.fish <-  subset(data.fish,n_prey_items>=6)
IDs = unique(data.fish$Statistical_Fish_ID_mass)
n.IDs = length(IDs)

# code IDs again, that they run from 1:n.pred 
data.fish$ID_Pred_new <-  data.fish$Statistical_Fish_ID_mass |> 
  as.factor() |> 
  as.integer() 

n.total = nrow(data.fish)
n.pred = max(data.fish$ID_Pred_new)

# sequence for calculating densities
x.min = min(data.zoop$log10_mass)
x.max = max(data.zoop$log10_mass)
x.dens = seq(x.min,x.max,by=0.01)

# Running preference analysis:
preference_data <- as.data.frame(matrix(nrow=n.pred,ncol=7))
colnames(preference_data) <- c("Fish_no","zoo_mean","zoo_sd","fish_mean","fish_sd","pref_mean","pref_sd")
preference_data$Fish_no <- seq(1,n.pred,1)

for(i in 1:n.pred){
  # zoop density #############################################################
  # get location from fish subset i
  ID.loc = data.fish |> 
    subset(ID_Pred_new==i) |> 
    dplyr::select(location_group) |>
    slice(1) |>
    as.numeric()
  
  # location subset
  df = data.zoop |> subset(location_group==ID.loc)
  
  # weighted mean and sdev
  zoop.mean = weighted.mean(df$log10_mass, df$density_m2)
  zoop.sdev = sqrt(wtd.var(df$log10_mass, df$density_m2))
  
  preference_data$zoo_mean[i] <- zoop.mean
  preference_data$zoo_sd[i] <- zoop.sdev
  
  # initialize density vector
  y.dens.zoop = rep(0,length(x.dens))
  
  # kernel density 
  h.zoop = 1.06*zoop.sdev*nrow(df)^(-1/5) # Scott's rule bandwidth
  h.zoop = h.zoop*bw.scale.zoop
  h.zoop = max(h.zoop, 0.1)
  # add kernels
  for(j in 1:nrow(df)){
    y.dens.zoop = y.dens.zoop + df$density_m2[j] * kern(x=x.dens, xi=df$log10_mass[j], h=h.zoop)
  }
  # normalize that it integrates to 1
  y.dens.zoop = y.dens.zoop / sum(y.dens.zoop*(x.max-x.min)/length(x.dens) )
  
  # fish observations ##################################################
  
  df = data.fish |> subset(ID_Pred_new==i)
  
  # replicate observations from prey counts
  # prey.mass.log.rep = rep(log10(df$Mean.Prey.Mass..g.), df$Prey.Count) 
  prey.mass.log.rep = rep(log10(df$Mean_Prey_Mass__g_), df$Prey_Count) 
  
  fish.mean = mean(prey.mass.log.rep)
  fish.sdev = sd(prey.mass.log.rep)
  
  preference_data$fish_mean[i] <- fish.mean
  preference_data$fish_sd[i] <- fish.sdev
  # pref mean ############################################################
  
  # empirical pref mean
  zoop.pdf = rep(0, length(prey.mass.log.rep))
  for(j in 1:length(prey.mass.log.rep)){
    # which x.dens is closest
    zoop.pdf[j] = y.dens.zoop[ which.min(abs(prey.mass.log.rep[j]-x.dens)) ]
  }
  
  # weighted means and sdev. weights are inverse of environment density
  pref.mean = weighted.mean(prey.mass.log.rep, 1/zoop.pdf)
  pref.sdev = sqrt(wtd.var(prey.mass.log.rep, 1/zoop.pdf))
  
  preference_data$pref_mean[i] <- pref.mean
  preference_data$pref_sd[i] <- pref.sdev
}

#adding temperature values to these data
preference_data$sst_mean <- NA
preference_data$species <- NA
preference_data$fish_size_class <- NA

for(i in 1:nrow(preference_data)){
  index <- which(data.fish$ID_Pred_new==i)
  sst_subset <- data.fish$SST[index]
  # species_subset <- data.fish$Species.Code[index]
  species_subset <- data.fish$Fish_Species_Code[index]
  fish_size <- data.fish$log_10_fish_class_bin[index]
  
  preference_data$sst_mean[i] <- mean(sst_subset)
  preference_data$species[i] <- species_subset[1]
  preference_data$fish_size_class[i] <- fish_size[1]
  
}

write.csv(preference_data,"Data/Processed/diet_size_preferences.csv") #this is the main dataset to use for analyses
write.csv(data.fish,"Data/Processed/preference_data_fish_info.csv") #this is a record of the "statistical fish" 

### Selectivity models ####
#loading packages and data
library(ggplot2)
library(car)
library(dplyr)
library(nlme)
library(MuMIn)
library(gstat)
library(sp)
library(effects)
library(RColorBrewer)
library(car)

preferences <- read.csv("Data/Processed/diet_size_preferences.csv")
fish_data <- read.csv("Data/Processed/preference_data_fish_info.csv") #this is needed to associate the cruise and locations

preferences$Cruise <- NA
preferences$location_group <- NA
preferences$Lat <- NA
preferences$Lon <- NA

for(i in 1:nrow(preferences)){
  index <- which(fish_data$ID_Pred_new==i)
  cruise_subset <- fish_data$Cruise[index]
  location_subset <- fish_data$location_group[index]
  lat_subset <- fish_data$Lat[index]
  lon_subset <- fish_data$Lon[index]
  
  preferences$Cruise[i] <- cruise_subset[1] 
  preferences$location_group[i] <- location_subset[1] 
  preferences$Lat[i] <- mean(unique(lat_subset)) 
  preferences$Lon[i] <- mean(unique(lon_subset)) 
  
}

#calculating an average SST for each location 
SST <- preferences %>%
  group_by(location_group)%>%
  summarise(SST <- mean(sst_mean))
colnames(SST) <- c("Location","SST")

for(i in 1:nrow(SST)){
  index <- (which(preferences$location_group==SST$Location[i]))
  preferences$sst_mean[index] <- SST$SST[i]
}

## LMM of preference vs SST and size class 
lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

gls_mod <- gls(pref_mean~sst_mean*fish_size_class,data=preferences,method = "REML",control=lmec)
mod_1 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|species),data=preferences,method = "REML",control=lmec)
mod_2 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|location_group),data=preferences,method = "REML",control=lmec)
mod_3 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise),data=preferences,method = "REML",control=lmec)
mod_4 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|species,~1|location_group),data=preferences,method = "REML",control=lmec)
mod_5 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),data=preferences,method = "REML",control=lmec)
mod_6 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|location_group),data=preferences,method = "REML",control=lmec)
mod_7 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|location_group,~1|species),data=preferences,method = "REML",control=lmec)
mod_8 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1+sst_mean|species),data=preferences,method = "REML",control=lmec)
mod_9 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1+sst_mean|Cruise),data=preferences,method = "REML",control=lmec)
# mod_10 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|location_group,~1+sst_mean|species),data=preferences,method = "REML",control=lmec)
mod_11 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|location_group,~1+sst_mean|Cruise),data=preferences,method = "REML",control=lmec)
mod_12 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|location_group,~1+sst_mean|species),data=preferences,method = "REML",control=lmec)
# mod_13 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1+sst_mean|Cruise,~1|location_group,~1|species),data=preferences,method = "REML",control=lmec)
# mod_14 <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1+sst_mean|Cruise,~1|location_group,~1+sst_mean|species),data=preferences,method = "REML",control=lmec)

anova_frame <- as.data.frame(anova(gls_mod,mod_1,mod_2,mod_3,mod_4,mod_5,mod_6,mod_7,mod_8,mod_9,mod_11,mod_12))
anova_frame[order(anova_frame$AIC,decreasing = T),]

#best is mod 5
summary(mod_5)

windows()
plot(mod_5) #definite cone shape
qqnorm(resid(mod_5,type = "pearson"))
qqline(resid(mod_5,type = "pearson")) 
plot(factor(preferences$species),resid(mod_5,type="pearson")) #variation by species
plot(factor(preferences$Cruise),resid(mod_5,type="pearson")) #variation by Cruise
plot(factor(preferences$location_group),resid(mod_5,type="pearson")) #variation by location_group

E <- resid(mod_5)
coplot(E ~ sst_mean | factor(species), data = preferences) #spread against SST doesnt look too bad

#now running with weights
mod_ident_sp <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),
                      weights = varIdent(form= ~ 1 | species),
                      data=preferences,method = "REML",control=lmec) 

mod_ident_Cruise <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),
                        weights = varIdent(form= ~ 1 | Cruise),
                        data=preferences,method = "REML",control=lmec) 

mod_fix <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),
                 weights = varFixed(~sst_mean),
                 data=preferences,method = "REML",control=lmec) 

mod_exp <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),
                 weights = varExp(form = ~sst_mean),
                 data=preferences,method = "REML",control=lmec) 

mod_const <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),
                   weights = varConstPower(form=~sst_mean),
                   data=preferences,method = "REML",control=lmec) 

mod_comb_sp_Cruise <- lme(pref_mean~sst_mean*fish_size_class,random = list(~1|Cruise,~1|species),
                          weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
                          data=preferences,method = "REML",control=lmec)


anova_frame <- anova(mod_5,mod_ident_sp,mod_ident_Cruise,mod_fix,mod_exp,mod_const,
                     mod_comb_sp_Cruise)
anova_frame[order(anova_frame$AIC,decreasing = T),]


summary(mod_1_comb_sp_Cruise)

plot(mod_1_comb_sp_Cruise) #looks much better than mod_1_ident_Cruise
qqnorm(resid(mod_1_comb_sp_Cruise,type = "pearson"))
qqline(resid(mod_1_comb_sp_Cruise,type = "pearson")) #bottom tail is worse than in mod_1_ident_Cruise
plot(factor(preferences$species),resid(mod_1_comb_sp_Cruise,type="pearson")) #better than mod_1_ident_Cruise
plot(factor(preferences$Cruise),resid(mod_1_comb_sp_Cruise,type="pearson")) #similar to mod_1_ident_Cruise
plot(factor(preferences$location_group),resid(mod_1_comb_sp_Cruise,type="pearson")) 

E <- resid(mod_1_comb_sp_Cruise)
coplot(E ~sst_mean | factor(species), data = preferences) 

#comparing different fixed effects structures
mod_interact <- lme(pref_mean~sst_mean*fish_size_class,random = list(Cruise=~1,species=~1),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
                    data=preferences,method = "ML",control=lmec)
mod_plus <- lme(pref_mean~sst_mean+fish_size_class,random = list(Cruise=~1,species=~1),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
                data=preferences,method = "ML",control=lmec)
mod_SST <- lme(pref_mean~sst_mean,random = list(Cruise=~1,species=~1),
               weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
               data=preferences,method = "ML",control=lmec)
mod_size <- lme(pref_mean~fish_size_class,random = list(Cruise=~1,species=~1),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
                data=preferences,method = "ML",control=lmec)
mod_null <- lme(pref_mean~1,random = list(Cruise=~1,species=~1),
                weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
                data=preferences,method = "ML",control=lmec)

anova_frame <- anova(mod_interact,mod_plus,mod_SST,mod_size,mod_null)
anova_frame[order(anova_frame$AIC,decreasing = T),]

#now fitting model with REML to extract coefficients
mod_interact <- lme(pref_mean~sst_mean*fish_size_class,random = list(Cruise=~1,species=~1),
                    weights = varComb(varIdent(form =~ 1 | Cruise),varIdent(form =~ 1 | species)),
                    data=preferences,method = "REML",control=lmec)
summary(mod_interact)
r.squaredGLMM(mod_interact)

#checking autocorrelation
temp_data = data.frame(error = residuals(mod_interact), x = preferences$Lon, y = preferences$Lat)
coordinates(temp_data) <- c("x","y") 
bubble(temp_data, "error", col = c("black","grey"),
       main = "Residuals", xlab = "X-coordinates", ylab = "Y-coordinates")

dists <- as.matrix(dist(cbind(preferences$Lon,preferences$Lat))) #distance matrix
dists <- 1/dists
dists[is.infinite(dists)] <- 0   #Distance value is inf for repeated observations from the same station
Moran.I(resid(mod_interact,type = "normalized"),dists) # suggests no autocorrelation in residuals - 0.306

#investigating interactions
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(glmmTMB)
library(RColorBrewer)
library(nlme)

lmec = lmeControl(msMaxIter = 500, msMaxEval = 500)

windows(record=T)

brewer.pal(7,"Dark2")
ef1 <- as.data.frame(effect("sst_mean*fish_size_class", mod_interact, xlevels=list(fish_size_class=c(min(preferences$fish_size_class),max(preferences$fish_size_class)))))
colours <- c("#1B9E77","#D95F02")

svg("pref_interac_min_max.svg",width=8,height=8)
ggplot(ef1, aes(x=sst_mean, y=fit,group=fish_size_class))+
  geom_jitter(preferences,mapping=aes(y=pref_mean,x=sst_mean,colour=fish_size_class),size=3,alpha=0.8,width=0.1)+
  scale_color_gradient(low=colours[1],high=colours[2],breaks=c(min(preferences$fish_size_class),max(preferences$fish_size_class)),labels=c(-0.525,1.575))+
  scale_fill_manual(labels = c("0.4", "2.2"),values = colours)+geom_smooth(method="lm",size=1.2,alpha=0.8,aes(colour=fish_size_class))+
  geom_ribbon(aes(ymin=fit-se, ymax=fit+se, fill=as.factor(fish_size_class)),alpha=0.3,show.legend = F)+
  labs(x= "SST (°C)", y="Size preference (log10 g)", color="Size class (log10 g)") + theme_classic()+
  theme(axis.text = element_text(size=25,colour="black"),axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28,vjust=1.8),legend.position = c(0.6,-0.17),legend.direction="horizontal",
        legend.text = element_text(size=18),legend.title = element_text(size=20),legend.background = element_blank())+
  scale_y_continuous(limits=c(),breaks=seq(-3,0,1))+scale_x_continuous(limits=c(-2,8),breaks=seq(-2,8,2))+theme(plot.margin = margin(0.1,0.1,1.7,0.2, "cm"),legend.spacing.x = unit(1.0, 'cm'))
dev.off()
## ZOOPLANKTON SIZE DISTRIBUTION -------------------------------------
rm(list=ls())
dev.off()

library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(tidyverse)
library(nlme)
library(reldist)
library(ggridges)
library(viridis)

data.zoop <- read.csv("Data/Processed/PLANKTON_DATA.csv")

#subsetting the data to exclude prey sizes that are outwith the distribution of sizes commonly consumed by the fish.
# (i.e very small prey that are rarely found are assumed to be accidental ingestion)
data.zoop <- subset(data.zoop,data.zoop$log10_mass>=(-4)&data.zoop$log10_mass<=(-0.0195))

#normalising the densities to make plots comparable
data.zoop <- data.zoop%>%
  group_by(location_group)%>%
  mutate(norm_dens=density_m2/sum(density_m2),mean_SST=round(mean(unique(SST)),3))

#also making a grouping of sites into 1C bins
data.zoop$bin_1C <- NA
for(i in 1:nrow(data.zoop)){
  if(data.zoop$mean_SST[i]<=(-1)){
    data.zoop$bin_1C[i] <- -1
  } else if (data.zoop$mean_SST[i]>(-1) & data.zoop$mean_SST[i]<=(0)){
    data.zoop$bin_1C[i] <- 0
  } else if (data.zoop$mean_SST[i]>(0) & data.zoop$mean_SST[i]<=(1)){
    data.zoop$bin_1C[i] <- 1
  } else if (data.zoop$mean_SST[i]>(1) & data.zoop$mean_SST[i]<=(2)){
    data.zoop$bin_1C[i] <- 2
  } else if (data.zoop$mean_SST[i]>(2) & data.zoop$mean_SST[i]<=(3)){
    data.zoop$bin_1C[i] <- 3
  } else if (data.zoop$mean_SST[i]>(3) & data.zoop$mean_SST[i]<=(4)){
    data.zoop$bin_1C[i] <- 4
  } else if (data.zoop$mean_SST[i]>(4) & data.zoop$mean_SST[i]<=(5)){
    data.zoop$bin_1C[i] <- 5
  } else if (data.zoop$mean_SST[i]>(5) & data.zoop$mean_SST[i]<=(6)){
    data.zoop$bin_1C[i] <- 6
  } else if (data.zoop$mean_SST[i]>(6) & data.zoop$mean_SST[i]<=(7)){
    data.zoop$bin_1C[i] <- 7
  } else if (data.zoop$mean_SST[i]>(7) & data.zoop$mean_SST[i]<=(8)){
    data.zoop$bin_1C[i] <- 8
  }
}

#creating a column representing the middle temperature value for each bin
data.zoop$middle_bin <- data.zoop$bin_1C-0.5
#normalising the densities to sum to 1 witin each location, for comparability, and calculating the weighted average mass
data.zoop <- data.zoop%>%
  group_by(bin_1C)%>%
  mutate(norm_dens=density_m2/sum(density_m2),mean_SST=round(mean(SST),2),mean_mass=weighted.mean(log10_mass,density_m2))


#identifying the number of hauls per temperature bin
data.zoop%>%
  group_by(bin_1C)%>%
  dplyr::summarise(count=n_distinct(location_group))

sp_plot <- ggplot(data.zoop,aes(x=log10_mass,y=as.factor(middle_bin),weight=(density_m2)))+scale_fill_gradient(low = "blue",high = "red")+labs(y="Body mass (log10 g)",x="Mean SST (°C)")+
  geom_density_ridges_gradient(aes(height=..density..,weight=(density_m2),fill=mean_SST),bw= 0.08,scale=1.3,stat="density")+
  theme_classic()+  
  theme(axis.text = element_text(size=20,colour="black"),axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=25,vjust=1.8),legend.position = "none")+labs(x="Body mass (log10 g)",y="SST (°C)")+
  scale_x_continuous(expand = c(0,0),limits=c(-4,0))+scale_y_discrete(expand = c(0,0))+
  geom_segment(aes(x =-2.97, y = 1, xend = -2.97, yend = 1.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.11, y = 2, xend = -3.11, yend = 2.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.34, y = 3, xend = -3.34, yend = 3.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-2.94, y = 4, xend = -2.94, yend = 4.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.38, y = 5, xend = -3.38, yend = 5.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.07, y = 6, xend = -3.07, yend = 6.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.07, y = 7, xend = -3.07, yend = 7.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.26, y = 8, xend = -3.26, yend = 8.9),linewidth=1,linetype="dashed")+
  geom_segment(aes(x =-3.04, y = 9, xend = -3.04, yend = 9.9),linewidth=1,linetype="dashed")

sp_plot+ annotate("text", y = seq(1.4,9.9,1), x = rep(-0.1,9),label = c("(4)","(3)","(2)","(5)","(1)","(4)","(2)","(2)","(1)"),
                  color="black",size=4) 
