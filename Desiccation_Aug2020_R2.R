##This is the code for the Demes et al. 2020. Mutliple stressors drive convergent evolution
##of performance properties in marine macrophytes.

##Code by K Demes and S Starko
##Archived at http://github.com/sstarko/Demes_et_al_2020_Desiccation
##Last updated: Aug 19, 2020

setwd("./")

#Read in data file
x<-read.csv(file="Desiccation_data_Aug2020_R2.csv")

#Change columns to numeric where appropriate
x$Drying.time<-as.numeric(x$Drying.time)
x$RWC<-as.numeric(x$RWC)

#Load in packages necessary to run the following code
library(ggplot2)
library(gridExtra)
library(nlme)
library(lme4)
library(car)
library(ggpubr)
library(tidyverse)
library(plotrix)
library(AICcmodavg)
library(plyr)
library(dplyr)

######## first to separate out subsets of the data for each series of experiments/analyses
# the first subset removes all of the time out of water curve data to focus on the discrete treatment data
mainx<-subset(x, Treatment!="Desiccation.profile")

#the following subset removes the samples that were completely dried to 0% RWC or are terrestrial plant tissues and focuses on comparing wet (100% RWC) vs. partially dehydrated (75% RWC) samples
main1<-subset(mainx, Treatment!="Dry"&Habitat!="Terrestrial")

# the following subset retains the values that are either fresh (100% RWC) or completely dried (0% RWC). This will be used to compare against literature values for land plants to compare wet and dried seaweeds to land plants
dry<-subset(x, RWC=="0" | RWC=="100" | Habitat=="Terrestrial")

# the last subset is for the profile tests where samples were held out of water for varying amounts of time to test Relative Water Content (RWC) as a function of time
profiles<-subset(x, Treatment=="Desiccation.profile")


##############################################################################################################
##############################################################################################################
##############################################################################################################
# THE FUNCTION BELOW IS DESIGNED TO CREATED NEW DATAFRAMES WITH RESPONSE VARIABLES SUMMARISED BY FACTOR LEVELS
##############################################################################################################
##############################################################################################################
##############################################################################################################

##Load in summarySE function before beginning analyses and figures
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N     = length2(xx[,col], na.rm=na.rm),
                      mean  = mean   (xx[,col], na.rm=na.rm),
                      sd    = sd     (xx[,col], na.rm=na.rm),
                      median= median(xx[,col],  na.rm=na.rm),
                      min   = min   (xx[,col],  na.rm=na.rm),
                      max   = max   (xx[,col],  na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  datac$range <-datac$max - datac$min   # Calculates the range of sizes 
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


################
#Analyses
################



######################################################################
######################################################################
##### Comparing material properties when fresh #######################
######################################################################
######################################################################
######################################################################

# wet only for comparison among species in mechanical traits when tissues are fresh
wetONLY<-subset(x, Treatment=="Wet")

#Separate dataframes for each clade
Ala_wet<-subset(wetONLY, Clade=="Alariaceae")
Cos_wet<-subset(wetONLY, Clade=="Agaraceae")
Fuc_wet<-subset(wetONLY, Clade=="Fucales")
Giant_wet<-subset(wetONLY, Clade=="Giant Kelps")
Gig_wet<-subset(wetONLY, Clade=="Gigartinaceae")
Hedo_wet<-subset(wetONLY, Clade=="Hedophyllum")
Pyr_wet<-subset(wetONLY, Clade=="Pyropia")
Zos_wet<-subset(wetONLY, Clade=="Zosteraceae")

#Fit models comparing material properties when fresh for the different clades

#Alariaceae
lm(Breaking.Force~Habitat, data= Ala_wet) %>% anova()
lm(Strength~Habitat, data= Ala_wet) %>% anova()
lm(Modulus~Habitat, data= Ala_wet) %>% anova()
lm(Cross.sectional.area~Habitat, data=Ala_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Ala_wet) %>% anova()

#Agaraceae
lm(Breaking.Force~Habitat, data= Cos_wet) %>% anova()
lm(Strength~Habitat, data= Cos_wet) %>% anova()
lm(Modulus~Habitat, data= Cos_wet) %>% anova()
lm(Cross.sectional.area~Habitat, data= Cos_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Cos_wet) %>% anova()

#Fucales
lm(Breaking.Force~Habitat, data= Fuc_wet) %>% anova()
lm(Strength~Habitat, data= Fuc_wet) %>% anova()
lm(Modulus~Habitat, data= Fuc_wet) %>% anova()
lm(Cross.sectional.area~Habitat, data=Fuc_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Fuc_wet) %>% anova()

#Giant Kelps
lm(Breaking.Force~Habitat, data= Giant_wet) %>% anova()
lm(Strength~Habitat, data= Giant_wet) %>% anova()
lm(Modulus~Habitat, data= Giant_wet) %>% anova()
lm(Cross.sectional.area~Habitat, data=Giant_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Giant_wet) %>% anova()

#Gigartinaceae
lm(Breaking.Force~Habitat, data= Gig_wet)  %>% anova()
lm(Strength~Habitat, data= Gig_wet)  %>% anova()
lm(Modulus~Habitat, data= Gig_wet)  %>% anova()
lm(Cross.sectional.area~Habitat, data=Gig_wet)  %>% anova()
lm(Thickness_mm~Habitat, data = Gig_wet) %>% anova()

#Hedophyllum
lm(Breaking.Force~Habitat, data= Hedo_wet) %>% anova()
lm(Strength~Habitat, data= Hedo_wet) %>% anova()
lm(Modulus~Habitat, data= Hedo_wet) %>% anova()
lm(Modulus~Habitat, data= Hedo_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Hedo_wet) %>% anova()

#Pyropia
lm(Breaking.Force~Habitat, data= Pyr_wet) %>% anova()
lm(Strength~Habitat, data= Pyr_wet) %>% anova()
lm(Modulus~Habitat, data= Pyr_wet) %>% anova()
lm(Modulus~Habitat, data= Pyr_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Pyr_wet) %>% anova()

#Zosteraceae
lm(Breaking.Force~Habitat, data= Zos_wet) %>% anova()
lm(Strength~Habitat, data= Zos_wet) %>% anova()
lm(Modulus~Habitat, data=Zos_wet) %>% anova()
lm(Modulus~Habitat, data= Zos_wet) %>% anova()
lm(Thickness_mm~Habitat, data = Zos_wet) %>% anova()


######################################################
######################################################
######################################################
###### RWC as a function of time out of water ########
######################################################
######################################################
######################################################

#first to define the exponential curve
edc<- RWC~100*exp(-b*Drying.time)
#The intercept (a*exp(-b*x)) is fixed at 100 as we don't really need to estimate this

## Now to compete models: null model (RWC as a function of time only; no clade or habitat effects); h1 (differences by clade, but not habitat); h2 (where habitat response is constant across clades); and h3 (where habitat response can vary across clades) 
m_n<-nls(RWC~100*exp(-b*Drying.time), profiles, start = c(b=c(0.1)))
AICc(m_n)
#AICc = 1280.563
#AIC = 1380.482

m_h_1<- nlme(edc, ## Exp. decay curve, no differences by habitat - but differences by Clade
          data = profiles,
          fixed = b ~ 1,
          random = b ~1|Clade,
          start = c(b=c(0.1)))
AICc(m_h_1)
#AICc = 1151.7
#AIC = 1151.5


m_h_2<- nlme(edc, ## Exp. decay curve, differences by habitat - all Clades have same habitat response
             data = profiles,
             fixed = b ~ as.factor(Habitat),
             random = b ~ 1|Clade,
             start = c(b=c(0.1,0.1)))

AICc(m_h_2)
#AICc = 1074.785
#AIC = 1074.5

#Test main hypothesis
anova(m_h_2)
##Provides general support for the hypothesis that intertidal species tend to dry out more slowly than subtidal species 

m_h_3<- nlme(edc, ## Exp. decay curve, differences by habitat - Clades vary in habitat response
             data = profiles,
             fixed = b ~ as.factor(Habitat),
             random = b ~ as.factor(Habitat)|Clade,
             start = c(b=c(0.1,0.1)))

AICc(m_h_3)
#AICc = 1037.512
#AIC = 1036.925

##Suggests that there is no "Universal" effect of habitat, some clades have differences and others do not 


##Next we extract model predictions of the time required to achieve 75% RWC (25% desiccation) and compare subtidal and intertidal species using a paired t-test
#extract coefficients
m_h_3$coefficients
#Clade coefficients
clade.rates<- data.frame(Inter.tidal.rate=m_h_3$coefficients$fixed[1]+m_h_3$coefficients$random$Clade[,1],Sub.tidal.rate=m_h_3$coefficients$fixed[1]+m_h_3$coefficients$fixed[2]+m_h_3$coefficients$random$Clade[,1]+m_h_3$coefficients$random$Clade[,2])
#Calculate time to 75% RWC
clade.rates$t.75.RWC.Inter<- -log(75/100)/clade.rates$Inter.tidal.rate
clade.rates$t.75.RWC.Sub<- -log(75/100)/clade.rates$Sub.tidal.rate

#Paired t.test - do intertidal species generally take more time to dehydrate to 75%
t.test(clade.rates$t.75.RWC.Inter,clade.rates$t.75.RWC.Sub,paired=T,conf.level=0.95)
#T = 2.515, p = 0.04, mean = 24.7 min (1.48,48)

##Extract model predictions
p.data<- expand.grid(Drying.time=seq(0,210),Habitat=c("Intertidal","Subtidal"),Clade=clades)
p.data$predicted<- predict(m_h_3,newdata=p.data)

######################################################
######################################################
######################################################
#### Material properties and partial desiccation #####
######################################################
######################################################
######################################################

############## models comparing materials properties by treatment and habitat, nested within random effect clade 

#Set optimizer for complex models
ctrl <- lmeControl(opt='optim')

###Force per width
zforce <- lme(Force_per_width ~ 1  , random = ~ 1|Clade, data = main1)
zforce_1 <- lme(Force_per_width ~ Treatment  , random = ~ 1|Clade, data = main1)
zforce_2 <- lme(Force_per_width ~ Habitat  , random = ~ 1|Clade, data = main1)
zforce_3 <- lme(Force_per_width ~ Habitat*Treatment  , random = ~ 1|Clade, data = main1)
zforce_4 <- lme(Force_per_width ~ Habitat * Treatment, random = ~ Treatment+Habitat|Clade, data = main1)
zforce_5 <- lme(Force_per_width ~ Habitat * Treatment, random = ~ Treatment*Habitat|Clade, data = main1,control=ctrl)

lsmeans(zforce_3, pairwise~Habitat*Treatment, adjust="tukey")

AICc(zforce)
AICc(zforce_1)
AICc(zforce_2)
AICc(zforce_3)
AICc(zforce_4)
AICc(zforce_5)

#Model with differences based on habitat and treatment highly supported

##Extensibility
zext <- lme(Extensibility ~ 1  , random = ~ 1|Clade, data = main1)
zext_1 <- lme(Extensibility ~ Treatment  , random = ~ 1|Clade, data = main1)
zext_2 <- lme(Extensibility ~ Habitat  , random = ~ 1|Clade, data = main1)
zext_3 <- lme(Extensibility ~ Habitat*Treatment  , random = ~ 1|Clade, data = main1)
zext_4 <- lme(Extensibility ~ Habitat * Treatment, random = ~ Treatment+Habitat|Clade, data = main1,control=ctrl)
zext_5 <- lme(Extensibility ~ Habitat * Treatment, random = ~ Treatment*Habitat|Clade, data = main1,control=ctrl)

lsmeans(zext_3, pairwise~Habitat*Treatment, adjust="tukey")

zext %>% AICc()
zext_1 %>% AICc()
zext_2 %>% AICc()
zext_3 %>% AICc()
zext_4 %>% AICc()
zext_5 %>% AICc()


#Modulus
zmod <- lme(log10(Modulus) ~ 1  , random = ~ 1|Clade, data = main1)
zmod_1 <- lme(log10(Modulus) ~ Treatment  , random = ~ 1|Clade, data = main1)
zmod_2 <- lme(log10(Modulus) ~ Habitat  , random = ~ 1|Clade, data = main1)
zmod_3 <- lme(log10(Modulus) ~ Habitat*Treatment  , random = ~ 1|Clade, data = main1)
zmod_4 <- lme(log10(Modulus) ~ Habitat * Treatment, random = ~ Treatment+Habitat|Clade, data = main1,control=ctrl)
zmod_5 <- lme(log10(Modulus) ~ Habitat * Treatment, random = ~ Treatment*Habitat|Clade, data = main1,control=ctrl)

zmod %>% AICc() 
zmod_1 %>% AICc()
zmod_2 %>% AICc()
zmod_3 %>% AICc()
zmod_4 %>% AICc()
zmod_5 %>% AICc()

lsmeans(zmod_3, pairwise~Habitat*Treatment, adjust="tukey")


#Strength
zstr <- lme(log10(Strength) ~ 1  , random = ~ 1|Clade, data = main1)
zstr_1 <- lme(log10(Strength) ~ Treatment  , random = ~ 1|Clade, data = main1)
zstr_2 <- lme(log10(Strength) ~ Habitat  , random = ~ 1|Clade, data = main1)
zstr_3 <- lme(log10(Strength) ~ Habitat*Treatment  , random = ~ 1|Clade, data = main1)
zstr_4 <- lme(log10(Strength) ~ Habitat * Treatment, random = ~ Treatment+Habitat|Clade, data = main1,control=ctrl)
zstr_5 <- lme(log10(Strength) ~ Habitat * Treatment, random = ~ Treatment*Habitat|Clade, data = main1,control=ctrl)

zstr %>% AICc()
zstr_1 %>% AICc()
zstr_2 %>% AICc()
zstr_3 %>% AICc()
zstr_4 %>% AICc()
zstr_5 %>% AICc()

lsmeans(zstr_3, pairwise~Habitat*Treatment, adjust="tukey")


################
#Figures
################

##############################################################################################################
##############################################################################################################
##############################################################################################################
###### FIGURE 1 WAS MADE MANUALLY - Stylized cladogram showing relationships between taxa 
##############################################################################################################
##############################################################################################################
##############################################################################################################



##############################################################################################################
##############################################################################################################
##############################################################################################################
###### FIGURE 2 - DESICCATION PROFILE DATA
##############################################################################################################
##############################################################################################################
##############################################################################################################

# This will be a multi-panel figure with 8 panels (representing clade) comparing RWC~time for each subtidal vs. intertidal comparison
# Each panel will have a curve for the intertidal and subtidal species [black vs. grey] and semi-transparent data points on curves
##'profiles' is subset with profile data only
### Next to graph data for each clade

#Make ggplot object
RWC_profile<-ggplot(profiles, aes(x=Drying.time, y=RWC, colour=Habitat)) + 
  geom_point(size=5)+
  labs(x = "Time out of water [minutes]", y = "Relative Water Content [%]")+
  theme_minimal()+
  geom_smooth(method = "nls", formula=y ~ (100*exp(-b * x)),
              method.args = list(start = c(b = 0.01)), size=2, se = FALSE) + 
  scale_colour_grey()+
  geom_hline(yintercept = 75, lty=2) + 
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade)

#View plot
RWC_profile

#Figure exported to PDF as 5 x 7 image

##############################################################################################################
##############################################################################################################
##############################################################################################################
###### FIGURE 3, 4, SUPPLEMENTAL FIGURES X, Y - FRESH VS. 75% RWC IN INTERTIDAL VS SUBTIDAL RELATIVES
##############################################################################################################
##############################################################################################################
##############################################################################################################

## 'main1" dataset is the subset of the data with RWC = 100 and 75 (i.e., the fresh vs. semi-dried comparison)

#Make summaries of each of the properties of interest (breaking force Nmm-1, Extensibility, Strength, Modulus)
ext<-summarySE(main1, measurevar="Extensibility", groupvars=c("Clade", "RWC", "Habitat", "Species", "Treatment"))
mod<-summarySE(main1, measurevar="Modulus", groupvars=c("Clade","RWC","Habitat", "Species", "Treatment"))
strength<-summarySE(main1, measurevar="Strength", groupvars=c("Clade","RWC","Habitat","Species", "Treatment"))
break.force<-summarySE(main1, measurevar="Force_per_width", groupvars=c("Clade","RWC","Habitat", "Species", "Treatment"))

#Make a dataframe with summary of means and SE
figure3.data<- data.frame(Clade= ext$Clade, Habitat= ext$Habitat, Species= ext$Species, RWC= ext$RWC, Treatment=ext$Treatment, 
                          ext.mean = ext$Extensibility, ext.se = ext$se, ext.n=ext$N, 
                          mod.mean = mod$Modulus, mod.se = mod$se, mod.n=mod$N, 
                          break.force.mean = break.force$Force_per_width, break.force.se = break.force$se, break.force.n=break.force$N, 
                          strength.mean = strength$Strength, strength.se = strength$se, strength.n=strength$N)

########################################
########################################
#############FIGURE 3###################
########################################
########################################

force_fig<- ggplot(figure3.data, aes(x=RWC, y=break.force.mean, colour=Habitat)) + 
  geom_pointrange(aes(ymin=break.force.mean-break.force.se, ymax=break.force.mean+break.force.se), size=2) +
  geom_line(size=2) +
  geom_point(size=5)+
  scale_colour_grey()+
  labs(x = "Desiccation treatment (RWC)", y = "Breaking Force")+
  theme_minimal()+
  scale_x_reverse()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade, scale="free")

force_fig
#Export pdf at 7.5 x 10.5

########################################
########################################
#############FIGURE 4###################
########################################
########################################

ext_fig<- ggplot(figure3.data, aes(x=RWC, y=ext.mean, colour=Habitat)) + 
  geom_pointrange(aes(ymin=ext.mean-ext.se, ymax=ext.mean+ext.se), size=2) +
  geom_line(size=2) +
  geom_point(size=5)+
  scale_colour_grey()+
  labs(x = "Desiccation treatment (RWC)", y = "Extensibility")+
  theme_minimal()+
  scale_x_reverse()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade, scale="free")


ext_fig
#Export pdf at 7.5 x 10.5

#####################################################
#####################################################
#############SUPPLEMENTARY FIGURE ###################
#####################################################
#####################################################

stress_fig<- ggplot(figure3.data, aes(x=RWC, y=strength.mean, colour=Habitat)) + 
  geom_pointrange(aes(ymin=strength.mean-strength.se, ymax=strength.mean+strength.se), size=2) +
  geom_line(size=2) +
  geom_point(size=5)+
  scale_colour_grey()+
  labs(x = "Desiccation treatment (RWC)", y = "Breaking stress (MPa)")+
  theme_minimal()+
  scale_x_reverse()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade, scale="free")


stress_fig
#Export pdf at 7.5 x 10.5
#####################################################
#####################################################
#############SUPPLEMENTARY FIGURE ###################
#####################################################
#####################################################

stiff_fig<- ggplot(figure3.data, aes(x=RWC, y=mod.mean, colour=Habitat)) + 
  geom_pointrange(aes(ymin=mod.mean-mod.se, ymax=mod.mean+mod.se), size=2) +
  geom_line(size=2) +
  geom_point(size=5)+
  scale_colour_grey()+
  labs(x = "Desiccation treatment (RWC)", y = "Stiffness (MPa)")+
  theme_minimal()+
  scale_x_reverse()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade, scale="free")

stiff_fig

##############################################################################################################
##############################################################################################################
##############################################################################################################
###### FIGURE 5 - PLOTTING WET AND COMPLETELY DRY DATA WITH COMPARISON TO TERRESTRIAL PLANTS
##############################################################################################################
##############################################################################################################
##############################################################################################################

# for the species, I have dry data, I should get the mean +/- SE for wet and dry samples and plot in axes of strength vs. flexibility, adding other seaweeds and land plants
# land plant data is still to be populated

##getting summary stats for graph

#creating summary stats for each response variable
dry.ext<-summarySE(dry, measurevar="Extensibility", groupvars=c("Clade", "Habitat", "Species", "RWC"))
dry.mod<-summarySE(dry, measurevar="Modulus", groupvars=c("Clade","Habitat", "Species", "RWC"))
dry.break.force<-summarySE(dry, measurevar="Breaking.Force", groupvars=c("Clade","Habitat", "Species", "RWC"))
dry.strength<-summarySE(dry, measurevar="Strength", groupvars=c("Clade","Habitat","Species", "RWC"))

#creating new dataframe with summary data
dry.fig.data<- data.frame(Clade= dry.ext$Clade, Habitat= dry.ext$Habitat, Species= dry.ext$Species,RWC=dry.ext$RWC, 
                          ext.mean = dry.ext$Extensibility, ext.se = dry.ext$se, ext.n= dry.ext$N, 
                          mod.mean = dry.mod$Modulus, mod.se = dry.mod$se, mod.n= dry.mod$N, 
                          break.force.mean = dry.break.force$Breaking.Force, break.force.se = dry.break.force$se, break.force.n=dry.break.force$N, 
                          strength.mean = dry.strength$Strength, strength.se = dry.strength$se, strength.n=dry.strength$N)

#first to change RWC values to wet vs. dry
dry.fig.data$RWC <- as.character(dry.fig.data$RWC)
dry.fig.data$RWC[dry.fig.data$RWC == "0"] <- "DRY"
dry.fig.data$RWC[dry.fig.data$RWC == "100"] <- "WET"
dry.fig.data$RWC[dry.fig.data$RWC == "fibre"] <- "FIBRE"
dry.fig.data$RWC[dry.fig.data$RWC == "wood"] <- "WOOD"
dry.fig.data$RWC[dry.fig.data$RWC == "leaf"] <- "LEAF"


dry.fig.data$RWC<-factor(dry.fig.data$RWC, levels=c("WET", "DRY", "LEAF", "WOOD", "FIBRE"))

#now to create the figure

ggplot(data = dry.fig.data,aes(x = dry.strength$Strength, y = dry.mod$Modulus, colour= RWC, shape= Habitat)) + 
  geom_point(size=5) + 
  geom_errorbar(aes(ymin = dry.mod$Modulus - dry.mod$se, ymax = dry.mod$Modulus + dry.mod$se)) + 
  geom_errorbarh(aes(xmin = dry.strength$Strength - dry.strength$se,xmax = dry.strength$Strength + dry.strength$se))+
  geom_line(data = dry.fig.data, aes(x = dry.strength$Strength, y = dry.mod$Modulus, group = Species))+
  scale_color_manual(values=c("grey","black", "lightblue", "blue", "navy"))+
  labs(x = "Strength (MPa)", y = "Stiffness (MPa)")+
  theme_minimal()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  scale_x_log10()+
  scale_y_log10()



##############################################################################################################
##############################################################################################################
##############################################################################################################
###### SUPP FIGURE - DESICCATION PROFILE DATA
##############################################################################################################
##############################################################################################################
##############################################################################################################

# This will be a multi-panel figure with 8 panels (representing clade) comparing RWC~time for each subtidal vs. intertidal comparison
# Each panel will have a curve for the intertidal and subtidal species [black vs. grey] and semi-transparent data points on curves
##'profiles' is subset with profile data only
### Next to graph data for each clade

##Create columns describing water content of blades at different time points as well as initial water content
profiles$Total.water<-profiles$Tissue.Wet.Weight.g-profiles$Dry.Weight.g
profiles$WaterContent<-profiles$Desiccated.Weight.g-profiles$Dry.Weight.g
profiles$PercentWater<-(1 - profiles$Dry.Weight.g/profiles$Tissue.Wet.Weight.g)*100

##Supplementary Figure X
WaterContent_profile<-ggplot(profiles, aes(x=Drying.time, y=WaterContent, colour=Habitat)) + 
  geom_point(size=5)+
  labs(x = "Time out of water [minutes]", y = "Total water content [g]")+
  theme_minimal()+
  scale_colour_grey()+
  facet_wrap(~Clade, scale = "free")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  geom_smooth(method = "nls", formula=y ~ (a*exp(-b * x)),
  method.args = list(start = c(a = 1, b = 0.01)), size=2, se = FALSE) 

WaterContent_profile


#This figure (not included in paper) presents the mass of water loss as a function of time out of water
WaterLoss_profile<-ggplot(profiles, aes(x=Drying.time, y=Water.Lost.g, colour=Habitat)) + 
  geom_point(size=5)+
  labs(x = "Time out of water [minutes]", y = "Water lost [g]")+
  theme_minimal()+
  geom_smooth(method = "nls", formula=y ~ (a*(1-exp(-b * x))),
              method.args = list(start = c(b = 0.01, a = 0.5)), size=2, se = FALSE) + 
  scale_colour_grey()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade, scale = "free")

WaterLoss_profile


#Summarize dry weight and wet weight 
DryWeight<-summarySE(profiles, measurevar="Dry.Weight.g", groupvars=c("Clade","Habitat","Species"))
WetWeight<-summarySE(profiles, measurevar="Tissue.Wet.Weight.g", groupvars=c("Clade","Habitat","Species"))

#Combine dataframes with new column to distinguish whether it is wet or dry weight
colnames(WetWeight)[5]<-"Mean.weight"
WetWeight$Treatment<-1
colnames(DryWeight)[5]<-"Mean.weight"
DryWeight$Treatment<-0
WeightCompare<-rbind(WetWeight, DryWeight) 

WeightCompare %>% str()

Weight_fig<- ggplot(WeightCompare, aes(x=Treatment, y=Mean.weight, colour=Habitat)) + 
  geom_pointrange(aes(ymin=Mean.weight-se, ymax=Mean.weight+se), size=0.2) +
  geom_line(size=3) +
  geom_point(size=5)+
  scale_colour_grey()+
  labs(x = "Tissue state", y = "Tissue mass [g]")+
  theme_minimal()+
  scale_x_reverse(limits = c(1.2,-0.2))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))+
  facet_wrap(~Clade, scale="free")

Weight_fig


##Some additional analyses of tissue water content not included in the journal article

#Calculate total water summaries
TotalWater<-summarySE(profiles, measurevar="Total.water", groupvars=c("Clade","Habitat","Species"))

#Calculate percent water summaries
PercentWater<-summarySE(profiles, measurevar="PercentWater", groupvars=c("Clade","Habitat","Species"))



ggplot(TotalWater, aes(x=Habitat, y=Total.water, colour=Habitat)) + 
  geom_pointrange(aes(ymin=Total.water-se, ymax=Total.water+se), size=0.5) +
  theme_minimal()+
  geom_line(size=2) +
  geom_point(size=1)+
  scale_colour_grey()+
  labs(x = "Habitat", y = "Total Water (g)")+
  facet_wrap(~Clade, scale = "free")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))

lme(Total.water~Habitat, random = ~1|Clade, data=profiles) %>% summary() #Higher in intertidal
lme(Dry.Weight.g~Habitat, random = ~1|Clade, data=profiles) %>% summary() #Higher in intertidal
lme(PercentWater~Habitat, random = ~1|Clade, data=profiles) %>% summary() #Higher in subtidal


ggplot(DryWeight, aes(x=Habitat, y=Dry.Weight.g, colour=Habitat)) + 
    geom_pointrange(aes(ymin=Dry.Weight.g-se, ymax=Dry.Weight.g+se), size=0.5) +
    theme_minimal()+
    geom_line(size=2) +
    geom_point(size=1)+
    scale_colour_grey()+
    labs(x = "Habitat", y = "Dry Weight (g)")+
    facet_wrap(~Clade, scale = "free")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))

ggplot(PercentWater, aes(x=Habitat, y=PercentWater, colour=Habitat)) + 
  geom_pointrange(aes(ymin=PercentWater-se, ymax=PercentWater+se), size=0.5) +
  theme_minimal()+
  geom_line(size=2) +
  geom_point(size=1)+
  scale_colour_grey()+
  labs(x = "Habitat", y = "Percent Water (%)")+
  facet_wrap(~Clade, scale = "free")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face="bold"))+
  theme(strip.text.x = element_text(size = 16))+
  theme(legend.position = c(0.8, 0.2))+
  theme(legend.text=element_text(size=14))
  
