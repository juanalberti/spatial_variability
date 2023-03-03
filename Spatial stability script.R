#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                          spatial stability script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#used packages####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#use these commands if you are not sure you have the packages...
#install.packages("pacman")
#pacman::p_load(vegan, BiodiversityR, ggplot2, cowplot, gridExtra, nlme, lme4, 
#               tidyr, data.table,lmodel2,dplyr,broom,RColorBrewer, piecewiseSEM)
#otherwise..
library(vegan)
#library(BiodiversityR)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(nlme)
library("lme4")
library(tidyr)
library(data.table)
library(lmodel2)
library(dplyr)
library(broom)
library(RColorBrewer)
library("piecewiseSEM")
library(maptools)
library(mobr)
library("varhandle")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#used data####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#There are two data sets, one for the pre-treatment (observational_data.csv) data, and other 
#for the experimentally increased heterogeneity (increased heterogeneity.csv). In addition, pre-treatment data has a subset of sites
#in where soil conditions where measured (observational_data_soil.csv) and a version in where data for each site were sumarized (observational_data_site_climate.csv).

Pre_Data<- read.csv("observational_data.csv", header=TRUE, sep=",", 
                na.strings=c("NA", "NULL"), dec=".", strip.white=TRUE)
Pre_Data_site<- read.csv("observational_data_site_climate.csv", header=TRUE, sep=",", 
                         na.strings=c("NA", "NULL"), dec=".", strip.white=TRUE)
Pre_Data_soil<- read.csv("observational_data_soil.csv", header=TRUE, sep=",", 
                    na.strings=c("NA", "NULL"), dec=".", strip.white=TRUE)
Pre_Post_Data<- read.csv("increased heterogeneity.csv", header=TRUE, sep=",", 
                     na.strings="NA", dec=".", strip.white=TRUE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#mixed-effect models to evaluate the relation between different levels of diversity and spatial variability####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#variability vs alpha richness
lmer2<-lmer(log(variability)~richness + (1|site), data=Pre_Data)
lmer4<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2,lmer4)#this is reported as chi2 in the ms
#perform the same but with lme2 becasue is easier to extract coefficients to plot
ctrl <- lmeControl(opt='optim')
lme2<-lme(fixed=log(variability)~richness, random= ~1|site, 
          method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2)$fixed#reported in the ms

#variability VS gama richness
lmer2g<-lmer(log(variability)~gama + (1|site), data=Pre_Data)
lmer4g<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2g,lmer4g)
lme2_gama<-lme(fixed=log(variability)~gama, random= ~1|site, method = "REML", 
               control=ctrl, data=Pre_Data)
intervals(lme2_gama)$fixed

#variability VS beta diversity (multivariate presence/absence)
lmer2b<-lmer(log(variability)~mult_beta_pa + (1|site), data=Pre_Data)
lmer4b<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2b,lmer4b)
lme2_beta_pa<-lme(fixed=log(variability)~mult_beta_pa, random= ~1|site, control=ctrl, 
                  data=Pre_Data)
intervals(lme2_beta_pa)$fixed


###the two separate components of spatial variability (mean and sd) VS different levels of diversity####

#alpha richness
#mean
lmer2mean<-lmer(log(live_mass)~richness + (1|site), data=Pre_Data)
lmer4mean<-lmer(log(live_mass)~1 + (1|site), data=Pre_Data)
anova(lmer2mean,lmer4mean)
lme2mean<-lme(fixed=log(live_mass)~richness, random= ~1|site, 
              method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2mean)$fixed
#sd
lmer2sd<-lmer(log(sd_biomass)~richness + (1|site), data=Pre_Data)
lmer4sd<-lmer(log(sd_biomass)~1 + (1|site), data=Pre_Data)
anova(lmer2sd,lmer4sd)
lme2sd<-lme(fixed=log(sd_biomass)~richness, random= ~1|site, 
            method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2sd)$fixed

#gamma richness
#mean
lmer2meang<-lmer(log(live_mass)~gama + (1|site), data=Pre_Data)
lmer4meang<-lmer(log(live_mass)~1 + (1|site), data=Pre_Data)
anova(lmer2meang,lmer4meang)
lme2mean_g<-lme(fixed=log(live_mass)~gama, random= ~1|site, 
                method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2mean_g)$fixed
#sd
lmer2sdg<-lmer(log(sd_biomass)~gama + (1|site), data=Pre_Data)
lmer4sdg<-lmer(log(sd_biomass)~1 + (1|site), data=Pre_Data)
anova(lmer2sdg,lmer4sdg)
lme2sdg<-lme(fixed=log(sd_biomass)~gama, random= ~1|site, 
             method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2sdg)$fixed

#beta 
#mean
lmer2meanb<-lmer(log(live_mass)~mult_beta_pa + (1|site), data=Pre_Data)
lmer4meanb<-lmer(log(live_mass)~1 + (1|site), data=Pre_Data)
anova(lmer2meanb,lmer4meanb)
lme2mean_b<-lme(fixed=log(live_mass)~mult_beta_pa, random= ~1|site, 
                method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2mean_b)$fixed
#sd
lmer2sdb<-lmer(log(sd_biomass)~mult_beta_pa + (1|site), data=Pre_Data)
lmer4sdb<-lmer(log(sd_biomass)~1 + (1|site), data=Pre_Data)
anova(lmer2sdb,lmer4sdb)
lme2sdb<-lme(fixed=log(sd_biomass)~mult_beta_pa, random= ~1|site, 
             method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2sdb)$fixed

############################################################
#####species biomass spatial covariation VS richness#####
#alpha richness
lmer2_cov_r<-lmer(log(1-complementarity)~richness + (1|site), data=Pre_Data)
lmer4_cov_r<-lmer(log(1-complementarity)~1 + (1|site), data=Pre_Data)
anova(lmer2_cov_r,lmer4_cov_r)
lme2_cov_r=lme(fixed=log(1-complementarity)~richness, random= ~1|site, 
               control=ctrl, data=Pre_Data)
intervals(lme2_cov_r)$fixed
#gamma richness
lmer2_cov_g<-lmer(log(1-complementarity)~gama + (1|site), data=Pre_Data)
lmer4_cov_g<-lmer(log(1-complementarity)~1 + (1|site), data=Pre_Data)
anova(lmer2_cov_g,lmer4_cov_g)
lme2_cov_g=lme(fixed=log(1-complementarity)~gama, random= ~1|site, 
               control=ctrl, data=Pre_Data)
intervals(lme2_cov_g)$fixed
#beta
lmer2_cov_b<-lmer(log(1-complementarity)~mult_beta_pa + (1|site), data=Pre_Data)
lmer4_cov_b<-lmer(log(1-complementarity)~1 + (1|site), data=Pre_Data)
anova(lmer2_cov_b,lmer4_cov_b)
lme2_cov_b=lme(fixed=log(1-complementarity)~mult_beta_pa, random= ~1|site, 
               control=ctrl, data=Pre_Data)
intervals(lme2_cov_b)$fixed

#######################################################################
####variability vs covariation####
lmer2_cov_v<-lmer(log(variability)~log(1-complementarity) + (1|site), data=Pre_Data)
lmer4_cov_v<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2_cov_v,lmer4_cov_v)
lme2_syncrony=lme(fixed=log(variability)~log(1-complementarity), 
                  random= ~1|site, control=ctrl, data=Pre_Data)
intervals(lme2_syncrony)$fixed

######################################################################
####other measures of diversity (instead of richness)####

#variability vs alpha simpson
lmer2_sim<-lmer(log(variability)~simpson + (1|site), data=Pre_Data)
lmer4_sim<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2_sim,lmer4_sim)#this is reported as chi2 in the ms
lme2_sim<-lme(fixed=log(variability)~simpson, random= ~1|site, 
              method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2_sim)$fixed

#variability vs shannon
lmer2_shan<-lmer(log(variability)~shannon + (1|site), data=Pre_Data)
lmer4_shan<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2_shan,lmer4_shan)#this is reported as chi2 in the ms

lme2_shan<-lme(fixed=log(variability)~shannon, random= ~1|site, 
               method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2_shan)$fixed

#variability VS Spie
lmer2p<-lmer(log(variability)~Spie + (1|site), data=Pre_Data)
lmer4p<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2p,lmer4p)
lme2_pie<-lme(fixed=log(variability)~Spie, random= ~1|site, method = "REML", 
              control=ctrl, data=Pre_Data)
intervals(lme2_pie)$fixed

#variability vs gamma simpson
lmer2_simg<-lmer(log(variability)~simpson_g + (1|site), data=Pre_Data)
lmer4_simg<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2_simg,lmer4_simg)
lme2_simg<-lme(fixed=log(variability)~simpson_g, random= ~1|site, 
               method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2_simg)$fixed

#variability vs gamma shannon
lmer2_shang<-lmer(log(variability)~shannon_g + (1|site), data=Pre_Data)
lmer4_shang<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2_shang,lmer4_shang)#this is reported as chi2 in the ms
lme2_shang<-lme(fixed=log(variability)~shannon_g, random= ~1|site, 
                method = "REML", control=ctrl, data=Pre_Data)
intervals(lme2_shang)$fixed

###turnover beta
lmer2ab<-lmer(log(variability)~a_beta + (1|site), data=Pre_Data)
lmer4ab<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2ab,lmer4ab)
lme2_ab<-lme(fixed=log(variability)~a_beta, random= ~1|site, control=ctrl, data=Pre_Data)
intervals(lme2_ab)$fixed

###multiplicative beta
lmer2mb<-lmer(log(variability)~log(m_beta) + (1|site), data=Pre_Data)
lmer4mb<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2mb,lmer4mb)
lme2_mb<-lme(fixed=log(variability)~log(m_beta), random= ~1|site, control=ctrl, data=Pre_Data)
intervals(lme2_mb)$fixed

#multivariate abundance beta
lmer2_beta_a<-lmer(log(variability)~multi_beta + (1|site), data=Pre_Data)
lmer4_beta_a<-lmer(log(variability)~1 + (1|site), data=Pre_Data)
anova(lmer2_beta_a,lmer4_beta_a)
lme2_beta_a<-lme(fixed=log(variability)~multi_beta, random= ~1|site, control=ctrl, 
                 data=Pre_Data)
intervals (lme2_beta_a)$fixed

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#type II regression####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Here we used data sumarized by site (Pre_Data_site) 
type2<-lmodel2(log(variability) ~ richness, data=Pre_Data_site, nperm=999)
type2_gama<-lmodel2(log(variability) ~ gama, data=Pre_Data_site, nperm=999)
type2_beta<-lmodel2(log(variability) ~ mult_beta_pa, data=Pre_Data_site, nperm=999)
#to get the reported statistics
type2
type2_gama
type2_beta


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#piecewise Structural Equation Modeling (piecewiseSEM) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pre_Data$covariation <- 1-Pre_Data$complementarity
Pre_Data$variability_log<-log(Pre_Data$variability)
Pre_Data$biomass_log<-log(Pre_Data$live_mass)
Pre_Data$covariation_log<-log(Pre_Data$covariation)


#define individual models for response variables
variability_model<-lme(fixed=variability_log~gama + mult_beta_pa + richness +covariation_log
    + biomass_log, random= ~1|site, 
    method = "REML", control=ctrl, data=Pre_Data)
biomass_model<-lme(fixed=biomass_log~richness, random= ~1|site, 
    method = "REML", control=ctrl, data=Pre_Data)
covariation_model<-lme(fixed=covariation_log~richness + mult_beta_pa + gama, random= ~1|site, 
    method = "REML", control=ctrl, data=Pre_Data)

#initial model with bidirectional relations between levels of richness
model_full_sem<-psem(variability_model,
                     biomass_model,
                     covariation_model,
                     richness%~~%gama,
                     mult_beta_pa%~~%richness,
                     mult_beta_pa%~~%gama
)

# Get goodness-of-fit and AIC
summary(model_full_sem) #AIC=42; BIC=108

#remove richness->variability
variability_model1<-update(variability_model,.~.-richness)
model_1_sem<-psem(variability_model1,
                  biomass_model,
                  covariation_model,
                  richness%~~%gama,
                  mult_beta_pa%~~%richness,
                  mult_beta_pa%~~%gama
)
summary(model_1_sem) #AIC=38.2; BIC=99.4

#remove gamma->covariation
covariation_model1<-update(covariation_model,.~.-gama)
model_2_sem<-psem(variability_model1,
                  biomass_model,
                  covariation_model1,
                  richness%~~%gama,
                  mult_beta_pa%~~%richness,
                  mult_beta_pa%~~%gama
)
summary(model_2_sem) #AIC=39.1; BIC=96.8

#remove richness->biomass
model_3_sem<-psem(variability_model1,
                  covariation_model1,
                  richness%~~%gama,
                  mult_beta_pa%~~%richness,
                  mult_beta_pa%~~%gama
)
summary(model_3_sem) #AIC=32.8; BIC=76
#remove richness<->beta
model_4_sem<-psem(variability_model1,
                  covariation_model1,
                  richness%~~%gama,
                  mult_beta_pa%~~%gama
)
summary(model_4_sem) #AIC=32.8; BIC=76

final_model<-model_4_sem
coefs(final_model)


##################################################

#COMPLEMENTARY ANALYSIS FOR OBSERVATIONAL DATA####

#To remove the possible influence of key abiotic factors on the relationship between 
#different levels of biodiversity and the spatial variability of productivity, we used 
#a subset of bioclimatic variables. We performed a multiple regression of spatial variability 
#against these climatic variables, kept the residuals, and then modeled the relationship 
#between different levels of diversity and the obtained residuals, using type II regression 
##########residuals after regresion against all environmental variables####
xnam<-colnames(Pre_Data_site[3:25])
fmla_e <- as.formula(paste("log(variability) ~ ", paste(xnam, collapse= "+")))
lm_full_e<-lm(fmla_e, data=Pre_Data_site)
type2_env<-lmodel2(lm_full_e$residuals ~ richness, data=Pre_Data_site, nperm=999)
type2_env
type2_envg<-lmodel2(lm_full_e$residuals ~ gama, data=Pre_Data_site, nperm=999)
type2_envg
type2_envb<-lmodel2(lm_full_e$residuals ~ mult_beta_pa, data=Pre_Data_site, nperm=999)
type2_envb

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SEM adding variability in soil conditions#####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This analysis uses a subset of sites that measured soil conditions(Pre_Data_soil)

#define initial individual models 
variability_model_h<-lme(fixed=variability_log~gama + mult_beta_pa + richness +covariation_log
    + biomass_log+cv_soil, random= ~1|site, 
    method = "REML", control=ctrl, data=Pre_Data_soil)
biomass_model_h<-lme(fixed=biomass_log~richness, random= ~1|site, 
    method = "REML", control=ctrl, data=Pre_Data_soil)
covariation_model_h<-lme(fixed=covariation_log~richness + mult_beta_pa + gama, random= ~1|site, 
    method = "REML", control=ctrl, data=Pre_Data_soil)
beta_model_h<-lme(fixed=mult_beta_pa~cv_soil, random=~1|site, 
    method = "REML", control=ctrl, data=Pre_Data_soil)

#initial model
model_full_soil<-psem(variability_model_h,
                      biomass_model_h,
                      covariation_model_h,
                      beta_model_h,
                      richness%~~%gama,
                      mult_beta_pa%~~%richness,
                      mult_beta_pa%~~%gama
)
summary(model_full_soil)# 83   153

#remove biomass->variability
variability_model_h1<-update(variability_model_h, .~.-biomass_log)
model_1_soil<-psem(variability_model_h1,
                   biomass_model_h,
                   covariation_model_h,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%richness,
                   mult_beta_pa%~~%gama
                   )
summary(model_1_soil)# 67   134


#remove cv_soil->variability
variability_model_h2<-update(variability_model_h1, .~.-cv_soil)
model_2_soil<-psem(variability_model_h2,
                   biomass_model_h,
                   covariation_model_h,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%richness,
                   mult_beta_pa%~~%gama
)
summary(model_2_soil)# 65   129

#remove richness->biomass
model_3_soil<-psem(variability_model_h2,
                   covariation_model_h,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%richness,
                   mult_beta_pa%~~%gama
)
summary(model_3_soil)# 39.779   91.51


#remove richness->beta
model_4_soil<-psem(variability_model_h2,
                   covariation_model_h,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%gama
)
summary(model_4_soil) #39  91

#remove gama->variability
variability_model_h3<-update(variability_model_h2, .~.-gama)
model_5_soil<-psem(variability_model_h3,
                   covariation_model_h,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%gama
)
summary(model_5_soil) #39   88
anova(model_4_soil, model_5_soil)

#remove gama->covariation
covariation_model_h1<-update(covariation_model_h, .~.-gama)
model_6_soil<-psem(variability_model_h3,
                   covariation_model_h1,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%gama
)
summary(model_6_soil)#   35   81

#remove covariation~beta
covariation_model_h2<-update(covariation_model_h1, .~.-mult_beta_pa)
model_7_soil<-psem(variability_model_h3,
                   covariation_model_h2,
                   beta_model_h,
                   richness%~~%gama,
                   mult_beta_pa%~~%gama
)
summary(model_7_soil)

final_model_soil<-model_7_soil
##############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#After experimantally increased heterogeneity (post-treatment data)####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
post_h<-subset(Pre_Post_Data, Pre_Post_Data$year=="4")
pre_h<-subset(Pre_Post_Data, Pre_Post_Data$year=="0")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#test if slope changed (i.e. year*diversity significance) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#richness
lmer2_years=lmer(log(variability)~richness + year + richness:year + (1|site),  
                 data=Pre_Post_Data)
lmer3_years=lmer(log(variability)~richness + year + (1|site), data=Pre_Post_Data)
anova(lmer2_years, lmer3_years)
lme2_c=lme(fixed=log(variability)~richness, random= ~1|site, 
           control=ctrl, data=pre_h)
lme2_y4=lme(fixed=log(variability)~richness, random= ~1|site, 
            control=ctrl, data=post_h)
summary(lme2_c)
summary(lme2_y4)

#gamma
lmer2_g_years=lmer(log(variability)~gama + year + gama:year + (1|site),  
                   data=Pre_Post_Data)
lmer3_g_years=lmer(log(variability)~gama + year + (1|site), 
                   data=Pre_Post_Data)
anova(lmer2_g_years, lmer3_g_years)

lme2_g_c=lme(fixed=log(variability)~gama, random= ~1|site, control=ctrl, 
             data=pre_h)
lme2_g_y4=lme(fixed=log(variability)~gama, random= ~1|site, control=ctrl, 
              data=post_h)


#beta
lmer2_b_years=lmer(log(variability)~beta + year + beta:year + (1|site),  
                   data=Pre_Post_Data)
lmer3_b_years=lmer(log(variability)~beta + year + (1|site), 
                   data=Pre_Post_Data)
anova(lmer2_b_years, lmer3_b_years)
lme2_b_c=lme(fixed=log(variability)~beta, random= ~1|site, 
             control=ctrl, data=pre_h)
lme2_b_y4=lme(fixed=log(variability)~beta, random= ~1|site, 
              control=ctrl, data=post_h)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#SEM after increased heterogeneity ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#models for pre and post heterogeneity discriminated 
#pre heterogeneity:

variability_preh<-lme(fixed=variability_log~gama + beta + richness +covariation_log
    + biomass_log, random= ~1|site, 
    method = "REML", control=ctrl, data=pre_h)
biomass_preh<-lme(fixed=biomass_log~richness, random= ~1|site, 
    method = "REML", control=ctrl, data=pre_h)
covariation_preh<-lme(fixed=covariation_log~richness + beta + gama, random= ~1|site, 
    method = "REML", control=ctrl, data=pre_h)

model_full_sem_pre<-psem(variability_preh,
                         biomass_preh,
                         covariation_preh,
                         richness%~~%gama,
                         beta%~~%richness,
                         beta%~~%gama
                         )

# Get goodness-of-fit and AIC
summary(model_full_sem_pre)

#remove richness as predictor for variability 
variability_preh1<-update(variability_preh, .~.-richness)
model_1_sem_pre<-psem(variability_preh1,
                      biomass_preh,
                      covariation_preh,
                      richness%~~%gama,
                      beta%~~%richness,
                      beta%~~%gama
                      )


# Get goodness-of-fit and AIC
summary(model_1_sem_pre)

#remove richness as predictor for biomass

model_2_sem_pre<-psem(variability_preh1,
                       covariation_preh,
                       richness%~~%gama,
                       beta%~~%richness,
                       beta%~~%gama
                      )

# Get goodness-of-fit and AIC
summary(model_2_sem_pre)

#remove gama as predictor for covariation
covariation_preh1<-update(covariation_preh,.~.-gama)
model_3_sem_pre<-psem(variability_preh1,
                      covariation_preh1,
                      richness%~~%gama,
                      beta%~~%richness,
                      beta%~~%gama
                      )

# Get goodness-of-fit and AIC
summary(model_3_sem_pre)

#remove beta as predictor for covariation
covariation_preh2<-update(covariation_preh1,.~.-beta)
model_4_sem_pre<-psem(variability_preh1,
                      covariation_preh2,
                      richness%~~%gama,
                      beta%~~%richness,
                      beta%~~%gama
                      )

# Get goodness-of-fit and AIC
summary(model_4_sem_pre)

#remove biomass as predictors of variability
variability_preh2<-update(variability_preh1, .~.-biomass_log)
model_5_sem_pre<-psem(variability_preh2,
                      covariation_preh2,
                      richness%~~%gama,
                      beta%~~%richness,
                      beta%~~%gama
                      )

# Get goodness-of-fit and AIC
summary(model_5_sem_pre)

#remove beta<->richness
model_6_sem_pre<-psem(variability_preh2,
                      covariation_preh2,
                      richness%~~%gama,
                      beta%~~%gama
                      )

# Get goodness-of-fit and AIC
summary(model_6_sem_pre)

#model selection by AIC and BIC
summary(model_full_sem_pre)$IC
summary(model_1_sem_pre)$IC
summary(model_2_sem_pre)$IC
summary(model_3_sem_pre)$IC
summary(model_4_sem_pre)$IC
summary(model_5_sem_pre)$IC
summary(model_6_sem_pre)$IC
#Best model is model 6
coefs(model_6_sem_pre, pre_t, standardize = "scale")

#####################################
#post heterogeneity

variability_posth<-lme(fixed=variability_log~gama + beta + richness +covariation_log
                      + biomass_log, random= ~1|site, 
                      method = "REML", control=ctrl, data=post_h)
biomass_posth<-lme(fixed=biomass_log~richness, random= ~1|site, 
                  method = "REML", control=ctrl, data=post_h)
covariation_posth<-lme(fixed=covariation_log~richness + beta + gama, random= ~1|site, 
                      method = "REML", control=ctrl, data=post_h)

model_full_sem_post<-psem(variability_posth,
                          biomass_posth,
                          covariation_posth,
                          richness%~~%gama,
                          beta%~~%richness,
                          beta%~~%gama
                          )

# Get goodness-of-fit and AIC
summary(model_full_sem_post)

#remove richness->biomass

model_1_sem_post<-psem(variability_posth,
                       covariation_posth,
                       richness%~~%gama,
                       beta%~~%richness,
                       beta%~~%gama
                       )
# Get goodness-of-fit and AIC
summary(model_1_sem_post)

#remove richness->variability
variability_posth1<-update(variability_posth, .~.-richness)
model_2_sem_post<-psem(variability_posth1,
                       covariation_posth,
                       richness%~~%gama,
                       beta%~~%richness,
                       beta%~~%gama
                       )

# Get goodness-of-fit and AIC
summary(model_2_sem_post)

#remove gamma->covariation
covariation_posth1<-update(covariation_posth, .~.-gama)
model_3_sem_post<-psem(variability_posth1,
                       covariation_posth1,
                       richness%~~%gama,
                       beta%~~%richness,
                       beta%~~%gama
                       )

# Get goodness-of-fit and AIC
summary(model_3_sem_post)

#remove richness->covariation
covariation_posth2<-update(covariation_posth1, .~.-richness)
model_4_sem_post<-psem(variability_posth1,
                       covariation_posth2,
                       richness%~~%gama,
                       beta%~~%richness,
                       beta%~~%gama
                       )
summary(model_4_sem_post)

#remove beta<->richness
model_5_sem_post<-psem(variability_posth1,
                       covariation_posth2,
                       richness%~~%gama,
                       beta%~~%gama
                       )
summary(model_5_sem_post)


# Get goodness-of-fit and AIC
summary(model_full_sem_post)$IC
summary(model_1_sem_post)$IC
summary(model_2_sem_post)$IC
summary(model_3_sem_post)$IC
summary(model_4_sem_post)$IC
summary(model_5_sem_post)$IC

#multimodel to evaluate constrained and free parameters of paths

model_pre_post<-psem(
  lme(fixed=variability_log~gama + beta +covariation_log
      , random= ~1|site, 
      method = "REML", control=ctrl, data=Pre_Post_Data),
  lme(fixed=covariation_log~richness + beta, random= ~1|site, 
      method = "REML", control=ctrl, data=Pre_Post_Data),
  richness%~~%gama,
  beta%~~%gama
)
Pre_Post_Data$year<-as.factor(Pre_Post_Data$year)
multigroup <- multigroup(model_pre_post, group = "year")
#get the p values for the interactions with year 
multigroup$anovaInts
#get the coefficients for each group
multigroup$group.coefs
########################################################################

##########################PLOTS################################################

#fitted values for the fixed and random effects to plot the lines
to_plot<-Pre_Data
to_plot$F0 <- fitted(lme2, level = 0) 
to_plot$F1 <- fitted(lme2, level = 1) 
to_plot$F0_m <- fitted(lme2mean, level = 0) 
to_plot$F1_m <- fitted(lme2mean, level = 1)
to_plot$F0_mg <- fitted(lme2mean_g, level = 0) 
to_plot$F1_mg <- fitted(lme2mean_g, level = 1)
to_plot$F0_meb <- fitted(lme2mean_b, level = 0) 
to_plot$F1_meb <- fitted(lme2mean_b, level = 1)
to_plot$F0_sd <- fitted(lme2sd, level = 0) 
to_plot$F1_sd <- fitted(lme2sd, level = 1) 
to_plot$F0_sdg <- fitted(lme2sdg, level = 0) 
to_plot$F1_sdg <- fitted(lme2sdg, level = 1)
to_plot$F0_sdb <- fitted(lme2sdb, level = 0) 
to_plot$F1_sdb <- fitted(lme2sdb, level = 1)
to_plot$F0_sim <- fitted(lme2_sim, level = 0) 
to_plot$F1_sim <- fitted(lme2_sim, level = 1) 
to_plot$F0_shan <- fitted(lme2_shan, level = 0) 
to_plot$F1_shan <- fitted(lme2_shan, level = 1)
to_plot$F0_spie <- fitted(lme2_pie, level = 0)
to_plot$F1_spie <- fitted(lme2_pie, level = 1)
to_plot$F0_g <- fitted(lme2_gama, level = 0) 
to_plot$F1_g <- fitted(lme2_gama, level = 1) 
to_plot$F0_mb <- fitted(lme2_mb, level = 0) 
to_plot$F1_mb <- fitted(lme2_mb, level = 1)
to_plot$F0_bpa <- fitted(lme2_beta_pa, level = 0) 
to_plot$F1_bpa <- fitted(lme2_beta_pa, level = 1)
to_plot$F0_cov_r <- fitted(lme2_cov_r, level = 0) 
to_plot$F1_cov_r <- fitted(lme2_cov_r, level = 1) 
to_plot$F0_cov <- fitted(lme2_syncrony, level = 0) 
to_plot$F1_cov <- fitted(lme2_syncrony, level = 1) 
to_plot$F0_cov_g <- fitted(lme2_cov_g, level = 0) 
to_plot$F1_cov_g <- fitted(lme2_cov_g, level = 1)
to_plot$F0_cov_b <- fitted(lme2_cov_b, level = 0) 
to_plot$F1_cov_b <- fitted(lme2_cov_b, level = 1)
to_plot$F0_shang <- fitted(lme2_shang, level = 0) 
to_plot$F1_shang <- fitted(lme2_shang, level = 1)
to_plot$F0_simg <- fitted(lme2_simg, level = 0) 
to_plot$F1_simg <- fitted(lme2_simg, level = 1)
to_plot$F0_bab <- fitted(lme2_ab, level = 0) 
to_plot$F1_bab <- fitted(lme2_ab, level = 1)
to_plot$F0_mbab <- fitted(lme2_beta_a, level = 0) 
to_plot$F1_mbab <- fitted(lme2_beta_a, level = 1)

#create the shade to plot confidence intervals for typeII regresions  
#to plot alpha
shade_x<-seq(0,36)
shade_y1<-type2$confidence.intervals[2,2]+type2$confidence.intervals[2,5]*shade_x
shade_y2<-type2$confidence.intervals[2,3]+type2$confidence.intervals[2,4]*shade_x
shade<-as.data.frame(cbind(x=shade_x, y1=shade_y1, y2=shade_y2))
shade$ymax<-pmax(shade$y1, shade$y2)
shade$ymin<-pmin(shade$y1, shade$y2)
#to plot gama
shade_x_g<-seq(0,82)
shade_y1_g<-type2_gama$confidence.intervals[2,2]+type2_gama$confidence.intervals[2,5]*shade_x_g
shade_y2_g<-type2_gama$confidence.intervals[2,3]+type2_gama$confidence.intervals[2,4]*shade_x_g
shade_g<-as.data.frame(cbind(x=shade_x_g, y1=shade_y1_g, y2=shade_y2_g))
shade_g$ymax<-pmax(shade_g$y1, shade_g$y2)
shade_g$ymin<-pmin(shade_g$y1, shade_g$y2)
#to plot multivariate beta
shade_x_b<-seq(-0.1,0.45, by=0.01)
shade_y1_b<-type2_beta$confidence.intervals[1,2]+type2_beta$confidence.intervals[1,5]*shade_x_b
shade_y2_b<-type2_beta$confidence.intervals[1,3]+type2_beta$confidence.intervals[1,4]*shade_x_b
shade_b<-as.data.frame(cbind(x=shade_x_b, y1=shade_y1_b, y2=shade_y2_b))
shade_b$ymax<-pmax(shade_b$y1, shade_b$y2)
shade_b$ymin<-pmin(shade_b$y1, shade_b$y2)

#fitted values for the fixed and random effects to plot the lines (pre and post heterogeneity)


pre_h$F0_y0 <- fitted(lme2_c, level = 0) # to get the fitted values for the fixed effects (for plots)
pre_h$F1_y0 <- fitted(lme2_c, level = 1) # to get the fitted values for the random effects
post_h$F0_y4 <- fitted(lme2_y4, level = 0) # to get the fitted values for the fixed effects
post_h$F1_y4 <- fitted(lme2_y4, level = 1) # to get the fitted values for the random effects
pre_h$F0_y0g <- fitted(lme2_g_c, level = 0) # to get the fitted values for the fixed effects (for plots)
pre_h$F1_y0g <- fitted(lme2_g_c, level = 1) # to get the fitted values for the random effects
post_h$F0_y4g <- fitted(lme2_g_y4, level = 0) # to get the fitted values for the fixed effects
post_h$F1_y4g <- fitted(lme2_g_y4, level = 1) 
pre_h$F0_y0b <- fitted(lme2_b_c, level = 0) # to get the fitted values for the fixed effects (for plots)
pre_h$F1_y0b <- fitted(lme2_b_c, level = 1) # to get the fitted values for the random effects
post_h$F0_y4b <- fitted(lme2_b_y4, level = 0) # to get the fitted values for the fixed effects
post_h$F1_y4b <- fitted(lme2_b_y4, level = 1) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fig.3####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fig_3A<-ggplot(to_plot, aes(richness, log(variability)))+
  labs(x="alpha diversity")+
  geom_line(aes(y=F0),col="darkgrey", size=2) + 
  geom_line(aes(y=F1, col=site),alpha=0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))+
  theme(plot.title = element_text(size = 10))

Fig_3B<-ggplot(to_plot, aes(gama, log(variability)))+
  labs(x="gamma diversity")+
  geom_line(aes(y=F0_g),col="darkgrey", size=2) + 
  geom_line(aes(y=F1_g, col=site), alpha=0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))+
  theme(plot.title = element_text(size = 10))

Fig_3C<-ggplot(to_plot, aes(mult_beta_pa, log(variability)))+
  labs(x="beta diversity")+
  geom_line(aes(y=F0_bpa),col="darkgrey", size=2) + 
  geom_line(aes(y=F1_bpa, col=site), alpha=0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))+
  theme(plot.title = element_text(size = 10))

Fig_3D<-ggplot(to_plot, aes(richness, log(1-complementarity)))+
  labs(x="alpha diversity",y="species covariation (log scale)")+
  geom_line(aes(y=F0_cov_r),col="lightgrey", size=2.5) + 
  geom_line(aes(y=F1_cov_r, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_3E<-ggplot(to_plot, aes(gama, log(1-complementarity)))+
  labs(x="gamma diversity",y="species covariation (log scale)")+
  geom_line(aes(y=F0_cov_g),col="lightgrey", size=2.5) + 
  geom_line(aes(y=F1_cov_g, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_3F<-ggplot(to_plot, aes(mult_beta_pa, log(1-complementarity)))+
  labs(x="beta diversity",y="species covariation (log scale)")+
  geom_line(aes(y=F0_cov_b),col="lightgrey", size=2.5) + 
  geom_line(aes(y=F1_cov_b, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_3G<- ggplot(to_plot, aes(log(1-complementarity), log(variability)))+
  labs(x="species covariation (log scale)", 
       y="Spatial variability (log scale)")+
  geom_line(aes(y=F0_cov),col="lightgrey", size=2.5) + 
  geom_line(aes(y=F1_cov, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))
#~~~~~~~~~~~~~~~~~~~~~~~~
##Fig.4####
#~~~~~~~~~~~~~~~~~~~~~~~~
Fig_4A<-ggplot(post_h, aes(richness, log(variability)))+
  labs(x="Alpha diversity", y="Spatial variability \n of biomass (log scale)")+ 
  geom_line(alpha=0.3, aes(x=pre_h$richness, y=pre_h$F0_y0), 
            col="lightgrey",  size=1.6) + 
  geom_line(aes(y=F0_y4), col="lightgrey", size=1.6) + 
  geom_line(aes(y=F1_y4, col=site), alpha=0.6, size=0.75)+
  geom_point(size=3, alpha=0.6, aes(colour = site)) +
  geom_point(size=3, alpha=0.3, aes(x=pre_h$richness, y=log(pre_h$variability), 
                                    colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+
  scale_colour_discrete(drop=TRUE,
                        limits = levels(Biomass_NA$site))

Fig_4B<-ggplot(post_h, aes(gama, log(variability)))+
  labs(x="Gamma diversity", y="Spatial variability \n of biomass (log scale)")+ 
  geom_line(alpha=0.3, aes(x=pre_h$gama, y=pre_h$F0_y0g), 
            col="lightgrey",  size=1.6) + 
  geom_line(aes(y=F0_y4g), col="lightgrey", size=1.6) + 
  geom_line(aes(y=F1_y4g, col=site), size=1)+
  geom_point(size=3, alpha=0.6, aes(colour = site)) +
  geom_point(size=3, alpha=0.3, aes(x=pre_h$gama, y=log(pre_h$variability), 
                                    colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(Biomass_NA$site)) 

Fig_4C<-ggplot(post_h, aes(beta, log(variability)))+
  labs(x="Beta diversity", y="Spatial variability \n of biomass (log scale)")+ 
  geom_line(alpha=0.3, aes(x=pre_h$beta, y=pre_h$F0_y0b), 
            col="red",  size=1.6) + 
  geom_line(aes(y=F0_y4b), col="lightgrey", size=1.6) + 
  geom_line(aes(y=F1_y4b, col=site), size=1)+
  geom_point(size=3, alpha=0.6, aes(colour = site)) +
  geom_point(size=3, alpha=0.3, aes(x=pre_h$beta, y=log(pre_h$variability), 
                                             colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(Biomass_NA$site))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fig.S2####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fig_S2A<-ggplot(Pre_Data_site, aes(richness, log(variability))) +
  labs(x="alpha diversity")+
  scale_x_continuous(expand = c(0, 0))+
  geom_ribbon(data=shade,inherit.aes=F, aes(x=x, ymin=ymin, ymax=ymax, linetype=NA),
              fill="darkgrey", alpha=0.4)+
  geom_abline(size = 1.5, col="darkgrey", aes(slope=type2$regression.results$Slope[2], 
                                              intercept=type2$regression.results$Intercept[2], ))+
  geom_point(size=3,alpha=0.5, aes(colour = site)) +
  coord_cartesian(xlim = c(0, 36))+ 
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))+
  theme(plot.title = element_text(size = 10))

Fig_S2B<-ggplot(Pre_Data_site, aes(gama, log(variability))) +
  labs(x="gamma diversity")+
  scale_x_continuous(expand = c(0, 0))+
  geom_ribbon(data=shade_g,inherit.aes=F, aes(x=x, ymin=ymin, ymax=ymax, linetype=NA),
              fill="darkgrey", alpha=0.4)+
  geom_abline(size = 1.5, col="darkgrey", aes(slope=type2_gama$regression.results$Slope[2], 
                                              intercept=type2_gama$regression.results$Intercept[2], ))+
  geom_point(size=3,alpha=0.5, aes(colour = site)) +
  coord_cartesian(xlim = c(0, 82))+ 
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))+
  theme(plot.title = element_text(size = 10))

Fig_S2C<-ggplot(Pre_Data_site, aes(mult_beta_pa, log(variability))) +
  labs(x="beta diversity")+
  scale_x_continuous(expand = c(0, 0))+
  geom_ribbon(data=shade_b,inherit.aes=F, aes(x=x, ymin=ymin, ymax=ymax, linetype=NA),
              fill="darkgrey", alpha=0.4)+
  geom_abline(size = 1.5, col="darkgrey", aes(slope=type2_beta$regression.results$Slope[1], 
                                              intercept=type2_beta$regression.results$Intercept[1], ))+
  geom_point(size=3,alpha=0.5, aes(colour = site)) +
  coord_cartesian(xlim = c(-0.015, 0.45))+ 
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0,0.2,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))+
  theme(plot.title = element_text(size = 10))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fig.S3####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fig_S3A<-ggplot(to_plot, aes(simpson, log(variability)))+
  labs(x="alpha diversity (simpson)",
       y="spatial variability (log scale)")+
  geom_line(aes(y=F0_sim),col="lightgrey", size=2) + 
  geom_line(aes(y=F1_sim, col=site), alpha=0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_S3B<-ggplot(to_plot, aes(shannon, log(variability)))+
  labs(x="alpha diversity (shannon)",
       y="spatial variability (log scale)")+
  geom_line(aes(y=F0_shan),col="lightgrey", size=2)+ 
  geom_line(aes(y=F1_shan, col=site),alpha=0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_S3C<-ggplot(to_plot, aes(Spie, log(variability)))+
  labs(x="alpha diversity (sPIE)",
       y="spatial variability (log scale)")+
  geom_line(aes(y=F0_spie),col="lightgrey", size=2) + 
  geom_line(aes(y=F1_spie, col=site), alpha=0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0), "cm"))+ 
  scale_x_continuous(limits=c(0,25))
scale_colour_discrete(drop=TRUE,
                      limits = levels(to_plot$site))

Fig_S3D<-ggplot(to_plot, aes(log(m_beta), log(variability)))+
  labs(x="beta diversity (multiplicative)",
       y="spatial variability (log scale)")+
  geom_line(aes(y=F0_mb),col="lightgrey", size=2) + 
  geom_line(aes(y=F1_mb, col=site), alpha= 0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_S3E<-ggplot(to_plot, aes(a_beta, log(variability)))+
  labs(x="beta diversity (aditive)",
       y="spatial variability (log scale)")+
  geom_line(aes(y=F0_bab),col="lightgrey", size=2) + 
  geom_line(aes(y=F1_bab, col=site), alpha= 0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

  Fig_S3F<-ggplot(to_plot, aes(multi_beta, log(variability)))+
    labs(x="beta diversity (multivariate, abundance-based)",
         y="spatial variability (log scale)")+
  geom_line(aes(y=F0_mbab),col="lightgrey", size=2) + 
  geom_line(aes(y=F1_mbab, col=site), alpha= 0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

  Fig_S3G<-ggplot(to_plot, aes(simpson_g, log(variability)))+
    labs(x="gamma diversity (simpson)",
         y="spatial variability (log scale)")+
    geom_line(aes(y=F0_simg),col="lightgrey", size=2) + 
    geom_line(aes(y=F1_simg, col=site), alpha= 0.5, size=0.75)+
    geom_point(size=3, alpha=0.5, aes(colour = site)) +
    theme_cowplot()+
    theme(legend.position = "none")+
    theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
    theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
    scale_colour_discrete(drop=TRUE,
                          limits = levels(to_plot$site))
  Fig_S3H<-ggplot(to_plot, aes(shannon_g, log(variability)))+
  labs(x="gamma diversity (shannon)",
       y="spatial variability (log scale)")+
  geom_line(aes(y=F0_shang),col="lightgrey", size=2) + 
  geom_line(aes(y=F1_shang, col=site), alpha= 0.5, size=0.75)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fig.S4####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fig_S4A<-ggplot(to_plot, aes(richness, log(live_mass)))+
  labs(x="alpha diversity", y="spatial mean (log scale)")+
  geom_line(aes(y=F0_m),col="lightgrey", size=2.5)+ 
  geom_line(aes(y=F1_m, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_S4B<-ggplot(to_plot, aes(gama, log(live_mass)))+
  labs(x="gamma diversity", y="spatial mean (log scale)")+
  geom_line(aes(y=F0_mg),col="lightgrey", size=2.5)+ 
  geom_line(aes(y=F1_mg, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))
  
Fig_S4C<-ggplot(to_plot, aes(mult_beta_pa, log(live_mass)))+
  labs(x="beta diversity", y="spatial mean (log scale)")+
  geom_line(aes(y=F0_meb),col="lightgrey", size=2.5)+ 
  geom_line(aes(y=F1_meb, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_S4D<-ggplot(to_plot, aes(richness, log(sd_biomass)))+
  labs(x="alpha diversity", y="spatial SD (log scale)")+
  geom_line(aes(y=F0_sd),col="lightgrey", size=2.5)+ 
  geom_line(aes(y=F1_sd, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))


Fig_S4E<-ggplot(to_plot, aes(gama, log(sd_biomass)))+
  labs(x="gamma diversity", y="spatial SD (log scale)")+
  geom_line(aes(y=F0_sdg),col="lightgrey", size=2.5)+ 
  geom_line(aes(y=F1_sdg, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

Fig_S4F<-ggplot(to_plot, aes(mult_beta_pa, log(sd_biomass)))+
  labs(x="beta diversity", y="spatial SD (log scale)")+
  geom_line(aes(y=F0_sdb),col="lightgrey", size=2.5)+ 
  geom_line(aes(y=F1_sdb, col=site), size=1.5)+
  geom_point(size=3, alpha=0.5, aes(colour = site)) +
  theme_cowplot()+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0), "cm"))+ 
  scale_colour_discrete(drop=TRUE,
                        limits = levels(to_plot$site))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#log response ratios (to build Fig. S6) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lrr<-post_h
lrr$LRR_rich<-log(post_h$richness/ pre_h$richness)
lrr$LRR_biomass<-log(post_h$biomass/ pre_h$biomass)
lrr$LRR_variability<-log(post_h$variability/ pre_h$variability)
lrr$LRR_sd<-log(post_h$sd_biomass/ pre_h$sd_biomass)
lrr$LRR_g<-log(post_h$gama/ pre_h$gama)
lrr$LRR_b<-log(post_h$beta/ pre_h$beta)
lrr$LRR_cov<-log(post_h$covariation/ pre_h$covariation)
mean(pre_h$beta)
LRR_plot<-aggregate(lrr[,(ncol(lrr)-7):ncol(lrr)], 
                    by=list(lrr$site), mean)
m_rich<-mean(LRR_plot$LRR_rich)
m_biomass<-mean(LRR_plot$LRR_biomass)
m_variability<-mean(LRR_plot$LRR_variability)
m_sd<-mean(LRR_plot$LRR_sd)
m_gamma<-mean(LRR_plot$LRR_g)
m_beta<-mean(LRR_plot$LRR_b)
m_cov<-mean(LRR_plot$LRR_cov)

error_rich <- qt(0.975,df=length(LRR_plot$LRR_rich)-1)*
  sd(LRR_plot$LRR_rich)/sqrt(length(LRR_plot$LRR_rich))
inf_rich <-m_rich-error_rich
sup_rich <-m_rich+error_rich
error_biomass <- qt(0.975,df=length(LRR_plot$LRR_biomass)-1)*
  sd(LRR_plot$LRR_biomass)/sqrt(length(LRR_plot$LRR_biomass))
inf_biomass <-m_biomass-error_biomass
sup_biomass <-m_biomass+error_biomass
error_variability <- qt(0.975,df=length(LRR_plot$LRR_variability)-1)*
  sd(LRR_plot$LRR_variability)/sqrt(length(LRR_plot$LRR_variability))
inf_variability <-m_variability-error_variability
sup_variability <-m_variability+error_variability
error_sd <- qt(0.975,df=length(LRR_plot$LRR_sd)-1)*
  sd(LRR_plot$LRR_sd)/sqrt(length(LRR_plot$LRR_sd))
inf_sd <-m_sd-error_sd
sup_sd <-m_sd+error_sd
error_gamma <- qt(0.975,df=length(LRR_plot$LRR_g)-1)*
  sd(LRR_plot$LRR_g)/sqrt(length(LRR_plot$LRR_g))
inf_gamma <-m_gamma-error_gamma
sup_gamma <-m_gamma+error_gamma
error_beta <- qt(0.975,df=length(LRR_plot$LRR_b)-1)*
  sd(LRR_plot$LRR_b)/sqrt(length(LRR_plot$LRR_b))
inf_beta <-m_beta-error_beta
sup_beta <-m_beta+error_beta
error_cov <- qt(0.975,df=length(LRR_plot$LRR_cov)-1)*
  sd(LRR_plot$LRR_cov)/sqrt(length(LRR_plot$LRR_cov))
inf_cov <-m_cov-error_cov
sup_cov <-m_cov+error_cov


logresponces <-data.frame(variable=c("Alpha\ndiversity", "Beta\ndiversity","Gamma\ndiversity",
                                     "Species\n biomass\n spatial co-variation", "Biomass\nmean", "Biomass\nsd", 
                                     "Spatial\n variability\nof biomass"), 
                          LRR=c(m_rich, m_beta, m_gamma, m_cov, m_biomass, m_sd, m_variability), 
                          inf=c(inf_rich, inf_beta, inf_gamma, inf_cov, inf_biomass, inf_sd, inf_variability), 
                          sup=c(sup_rich, sup_beta, sup_gamma, sup_cov, sup_biomass, sup_sd, sup_variability)) 
logresponces$variable2 <- factor(logresponces$variable, as.character(logresponces$variable))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Fig.S6 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

se <- ggplot(logresponces, aes(x=variable2, y=LRR, ymin = inf, ymax=sup)) + 
  geom_point(size=6) + 
  geom_errorbar(width = 0.5) +
  labs(x=" ", y="Log response ratios") + 
  theme_cowplot() + 
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size="12"), axis.title.y=element_text(size="20"),
        axis.line = element_line(colour = "black"), panel.border= element_rect(fill = "NA"),
        panel.grid.major = element_blank(), panel.background = element_blank())

Fig_S6<-se+ geom_hline(yintercept=0)




