# ###
# title: "Pilot Study snail-parasite dynamics Tanzania"
# author: "NCS"
# date: 6 March 2024
# ###

####GATHERING AND CLEANING DATA####

#This script includes snail level datasets (phases2-3) and waterbody level
#datasets (phases 1-3). Data from both types is combined to create the datasets
#used for generalized linear mixed effects models (GLMMs)

#set working directory to where your data files are saved
#setwd()


##open dataframe of snail size and shedding data to bring in infection data 
#for phases 2 and 3
#phase 1 data will come from a different format dataset
P2P3_bulinus=read.csv("P2P3_measureshed_Upload.csv")

library(tidyverse)
library(glmmTMB)
library(effects)
#library(dplyr)

#create summary of snail, schistosome and non-schistosome infection data
#including infection fail data for binomial models
p2p3_summaries=P2P3_bulinus%>%
  group_by(District_Waterbody, Phase)%>%
  summarise(NonSchistoSnails=sum(NonSchistoStatus),
            SchistoSnails=sum(SchistoStatus),
            BulinusNumber=n(),
            SchistoFailures = BulinusNumber-SchistoSnails,
            NonSchistoFailures = BulinusNumber-NonSchistoSnails)

p2p3_summaries=P2P3_bulinus%>%
  select(District_Waterbody,Date)%>%
  distinct(District_Waterbody, .keep_all = TRUE)%>%
  mutate(Date=as.Date(Date, format = "%m/%d/%y"))%>%
  full_join(p2p3_summaries)
  


###Open data frames for phase 1-3 with additional waterbody character data
##Phase1##
#this does not have non-schistosome data
P1_WB_correct=read.csv("Phase1WaterbodyData2021_corrected.csv")

library(stringr)
#wrangling for correct data, date type and matching name format: 
P1_WB_correct$Date=(as.Date(P1_WB_correct$Date, format="%d-%b-%y"))
P1_WB_correct$number_snails_collected=as.numeric(P1_WB_correct$number_snails_collected)
colnames(P1_WB_correct)[colnames(P1_WB_correct)=="number_snails_collected"] = "BulinusNumber"
colnames(P1_WB_correct)[colnames(P1_WB_correct)=="Total..VE"] = "SchistoSnails"
colnames(P1_WB_correct)[colnames(P1_WB_correct)=="water_level"] = "Water_level"

#calculate infection failures for binomial model
P1_WB_correct= P1_WB_correct%>% 
  mutate(SchistoFailures = BulinusNumber-SchistoSnails)

P1_WB_correct$Phase="Phase1"  
P1_WB_correct=P1_WB_correct%>%
  select(District, Village, CorrectedWBName, Water_level, BulinusNumber, SchistoSnails,
         SchistoFailures, Phase, Date)

P1_WB_correct$CorrectedWBName=str_replace(P1_WB_correct$CorrectedWBName, " $", "")
P1_WB_correct$District=str_replace(P1_WB_correct$District, " $", "")
P1_WB_correct$Village=str_replace(P1_WB_correct$Village, " $", "")
P1_WB_correct= P1_WB_correct %>% 
  unite("District_Waterbody", District, CorrectedWBName, sep="_",remove=F)

P1_WB_correct=P1_WB_correct%>%
  mutate(Date=as.Date(Date, format = "%Y-%m-%d"))


##Phase2##
P2_WB_correct=read.csv("Phase2WaterbodyData2021.csv")

#wrangling for correct date, date type and matching name format: 
P2_WB_correct$Date=(as.Date(P2_WB_correct$Date, format="%d-%b-%y"))
colnames(P2_WB_correct)[colnames(P2_WB_correct)=="NumberSnailsCollected"] = "Bulinus_Collected"
colnames(P2_WB_correct)[colnames(P2_WB_correct)=="SchistoInfectedSnails"] = "Bulinus_Infected"
P2_WB_correct$Phase="Phase2"


##Phase3##
P3_WB_correct=read.csv("Phase3WaterbodyData2021.csv")

#wrangling for correct date, date type and matching name format: 
P3_WB_correct$Date=(as.Date(P3_WB_correct$Date, format="%d-%b-%y"))
colnames(P3_WB_correct)[colnames(P3_WB_correct)=="Number_of.Bulinus_collected"] = "Bulinus_Collected"
colnames(P3_WB_correct)[colnames(P3_WB_correct)=="Number.of.Infected.snails"] = "Bulinus_Infected"
colnames(P3_WB_correct)[colnames(P3_WB_correct)=="Long_dimension..M."] = "Long_dimension"
P3_WB_correct$Phase="Phase3"


#merging Phase2 and Phase3 waterbody data sets
common_cols = intersect(colnames(P3_WB_correct), colnames(P2_WB_correct))
two_phases=rbind(
  subset(P2_WB_correct, select = common_cols),
  subset(P3_WB_correct, select = common_cols))


#ensuring there are no unintentional spaces
two_phases$CorrectedWBName=str_replace(two_phases$CorrectedWBName, " $", "")
two_phases$District=str_replace(two_phases$District, " $", "")
two_phases$Village=str_replace(two_phases$Village, " $", "")
two_phases= two_phases %>% 
  unite("District_Waterbody", District, CorrectedWBName, sep="_",remove=F)

two_phases$District_Waterbody=as.factor(two_phases$District_Waterbody)



#####COMBINE THREE DATASETS####
###combining necessary columns from waterbody sheets and snail infection summary 
#to create data needed from phases 2 and 3
P2P3_assembly=left_join(two_phases[c("District", "Village", "CorrectedWBName", "Date", "Water_level", 
                                     #  "Bulinus_Collected", "Bulinus_Infected", 
                                     "District_Waterbody",
                                       "Phase", "Long_dimension")], 
                        p2p3_summaries[c("District_Waterbody","Phase","NonSchistoSnails",
                                         "SchistoSnails","BulinusNumber","SchistoFailures",
                                         "NonSchistoFailures")], 
                        by=c("District_Waterbody","Phase"))

#add 0 to cases of no snails and no infections
P2P3_assembly$SchistoSnails[is.na(P2P3_assembly$SchistoSnails)] <- 0
P2P3_assembly$NonSchistoSnails[is.na(P2P3_assembly$NonSchistoSnails)] <- 0
P2P3_assembly$SchistoFailures[is.na(P2P3_assembly$SchistoFailures)] <- 0
P2P3_assembly$NonSchistoFailures[is.na(P2P3_assembly$NonSchistoFailures)] <- 0
P2P3_assembly$BulinusNumber[is.na(P2P3_assembly$BulinusNumber)] <- 0

#This is then combined to phase 1 data for use in models
common_cols2 = intersect(colnames(P2P3_assembly), colnames(P1_WB_correct))
three_phases=rbind(
  subset(P2P3_assembly, select = common_cols2),
  subset(P1_WB_correct, select = common_cols2))


three_phases=left_join(three_phases,P2P3_assembly[c("District_Waterbody","Phase","Long_dimension",
                         "NonSchistoSnails","NonSchistoFailures")], by=c("District_Waterbody","Phase"))


#create a column for cattle use permission determined by village officials
#Kisimas are only for human use. this script accounts for potential typos and 
#data entry irregularities

three_phases= three_phases%>% 
  mutate(CattleUsePerm = case_when(str_detect(District_Waterbody, "_Kisima ") ~0,
                                   str_detect(District_Waterbody, "_kisima ") ~0,
                                   TRUE~1))


###add 0s for dimensions of dry waterbodies
three_phases=three_phases%>%
  mutate(Long_dimension=case_when(Water_level=="Dry"~paste(0),
                                  TRUE~paste0(Long_dimension)))

#ensure this is numerical and adds NAs to rows absent of data
three_phases$Long_dimension=as.numeric(three_phases$Long_dimension)

####GEOGRAPHIC COORDINATES####
#bring in waterbody coordinates
All_WB_Coords=read_csv("PilotData_FinalCoords.csv")
three_phases=left_join(three_phases, All_WB_Coords[c("District_Waterbody", "LatToUse",
                        "LongToUse")], by=c("District_Waterbody"))

# bring in validated permanence data
WB_permdata=read_csv("PhaseData_FinalPermanence.csv")
three_phases=left_join(three_phases, WB_permdata[c("District_Waterbody", 
                                  "Permanence")], by=c("District_Waterbody"))


####SCHOOL DISTANCE TO WATERBODIES####
#bring in coordinates for schools 
schools=read.csv("School_Coords.csv")

library(geosphere)

#the function below identifies the closest school to each waterbody within each village

#first we isolate all the unique villages:
k=unique(three_phases$Village)

Final_coords=All_WB_Coords

Final_coords=left_join(Final_coords, subset(three_phases, 
                                            !duplicated(District_Waterbody))[c("District_Waterbody", "Village")])

Final_coords=Final_coords[c("District_Waterbody","LatToUse" , "LongToUse" ,"Village")]

#then for each of those villages run the function
SchoolDist_func=function(k)
{
  x=as.data.frame((subset(Final_coords, Village==paste0(k)))[,1])
  mat=distm(subset(Final_coords, Village==paste0(k))[,c('LongToUse','LatToUse')], 
            subset(schools, Village==paste0(k))[,c('LongToUse','LatToUse')], fun=distVincentyEllipsoid)
  x$dist=apply(mat, MARGIN = 1, FUN=min, na.rm=T)
  x$NearestSchool=(subset(schools, Village==paste0(k)))$School_Name[apply(mat, 1, which.min)]
  colnames(x)[1] ="District_Waterbody"
  x
}

#coelesce the school distance data
SchoolDist_output=lapply(k, SchoolDist_func)
school_distance=do.call("rbind.data.frame", SchoolDist_output)

#and join it to the full dataset
three_phases=left_join(three_phases, school_distance)

three_phases$Long_dimension=as.numeric(three_phases$Long_dimension)

three_phases$CattleUsePerm=as.factor(three_phases$CattleUsePerm)
TWO=subset(three_phases, Phase!="Phase1")


#####GENERALIZED LINEAR MIXED EFFECTS MODELS; GLMMs####
#Random effects only models (variance partitioning)

#schistosome infection
Sc_one=glmmTMB(cbind(SchistoSnails,SchistoFailures)~1+ (1|District/Village/CorrectedWBName), 
            data=subset(three_phases, BulinusNumber>0), 
            family ="binomial", REML = T)
summary(Sc_one)

##proportion of total variance explained by random effects (i.e. conditional R2)
#install.packages('performance')
library(performance)
r2(Sc_one)

# refit, dropping District in REs due to negligible variance (3.575e-08)
Sc_one=glmmTMB(cbind(SchistoSnails,SchistoFailures)~1+ (1|Village/CorrectedWBName), 
               data=subset(three_phases, BulinusNumber>0), 
               family ="binomial", REML = T)
summary(Sc_one)


library(performance)
r2(Sc_one) 
#74.5% of total variance explained by random effects (village and waterbody)
# waterbody explains 46.6% of total variance (62.5% of conditional variance; (6 / (6 + 3.6)))
# village explains 27.9% of total variance (37.5% of conditional variance; (3.6 / (6 + 3.6)))


#non-schistosome infection
NS_one=glmmTMB(cbind(NonSchistoSnails,NonSchistoFailures)~(1|District/Village/CorrectedWBName), 
               data=subset(TWO, BulinusNumber>0), 
               family ="binomial", REML = T)

summary(NS_one)
r2(NS_one)

# refit, dropping district and village in REs due to negligible variances 
#(6.699e-13 and 9.103e-08, respectively)
NS_one=glmmTMB(cbind(NonSchistoSnails,NonSchistoFailures)~(1|CorrectedWBName), 
               data=subset(TWO, BulinusNumber>0), 
               family ="binomial", REML = T)
summary(NS_one)
r2(NS_one)
#46.5% of total variance explained by random effects (waterbody)

#snail presence
three_phases=three_phases%>%
  mutate(SnailPres=case_when(BulinusNumber>0~1,
                             TRUE~0))

SN_one=glmmTMB(SnailPres~(1|District/Village/CorrectedWBName), 
              data=three_phases, 
              family ="binomial", REML = T)
summary(SN_one)
r2(SN_one)
#12.9% of total variance explained by random effects (district, village and waterbody)
# waterbody explains 7.52% of total variance (58.3% of conditional variance; (0.285 / (0.285+0.081+0.123)))
# village explains 2.14% of total variance (16.6% of conditional variance; (0.081 / (0.285+0.081+0.123)))
# district explains 3.25% of total variance (25.2% of conditional variance; (0.123 / (0.285+0.081+0.123)))



######full GLMM for schistosome infection prevalence
SchistGLMM=glmmTMB(cbind(SchistoSnails,SchistoFailures)~Permanence+
                     log(dist)+ CattleUsePerm+Phase+
                     (1|Village/CorrectedWBName), 
                   data=three_phases, REML=T, family ="binomial")
summary(SchistGLMM) 
plot(allEffects(SchistGLMM))


#full GLMM for nonschistosome infection prevalence
NonSchistGLMM=glmmTMB(cbind(NonSchistoSnails,NonSchistoFailures)~Permanence+
                        log(Long_dimension)+
                        log(dist)+ 
                        CattleUsePerm+Phase+
                        (1|District/Village/CorrectedWBName), 
                      #data=data2, REML=T,family ="binomial")
                      data=subset(TWO, BulinusNumber>0), REML=T,family ="binomial")
summary(NonSchistGLMM) 
plot(allEffects(NonSchistGLMM))

#Figures for Phase related differences in infection of both groups
#identifying infection rates and confidence intervals
summary(predictorEffect("Phase", SchistGLMM))
summary(predictorEffect("Phase", NonSchistGLMM))

#compiling that into a data frame
phasedata=data.frame(x=c("Phase 1","Phase 2", "Phase 3","Phase 2", "Phase 3"),
                     y=c(0.0006373844, 0.0009384530, 0.0011194598, 0.02446072, 0.04604671),
                     sp=c("Schist", "Schist","Schist", "NS","NS"))


#graph infection data and include confidence intervals
#schistosome infection
ggplot(phasedata[1:3,],aes(x,y))+
  labs(x="Phase", y="Infection prevalence")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) + 
  geom_errorbar(aes(ymin=c(0.0002946002, 0.0004353506, 0.0005062413),
                    ymax=c(0.001378468, 0.002021778, 0.002473642)), width=0.2)+
  geom_point(col=c( "red2"), size=5) +ylim(0,0.01) 

#nonschistosome infection
ggplot(phasedata[4:5,],aes(x,y))+
  labs(x="Phase", y="Prevalence nonschistsosome")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) + 
  geom_errorbar(aes(ymin=c(0.01893165, 0.03513281 ),ymax=c(0.03155265, 0.06013967 )), width=0.2)+
  geom_point(col=c("darkgoldenrod2"), size=5) +ylim(0,0.065) 

#Figures for cattle permission related differences in infection of both groups
#identifying infection rates and confidence intervals
summary(predictorEffect("CattleUsePerm", SchistGLMM))
summary(predictorEffect("CattleUsePerm", NonSchistGLMM))

#compiling that into a data frame
Cattleuse_data=data.frame(x=c("Not Permitted", "Permitted","Not Permitted", "Permitted"),
                          y=c(0.0004111905, 0.0014577949 , 0.0335207, 0.0304141 ),
                          sp=c("Schist", "Schist", "NS","NS"))

#graph infection data and include confidence intervals
#schistosome infection
ggplot(Cattleuse_data[1:2,],aes(x,y))+
  labs(x="Cattle Use", y="Infection prevalence")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) + 
  geom_errorbar(aes(ymin=c(0.0001424258, 0.0006649575) ,ymax=c(0.001186526, 0.003192922)), width=0.2)+
  geom_point(col=c("red2"), size=5) +ylim(0,0.01) 

#nonschistosome infection
ggplot(Cattleuse_data[3:4,],aes(x,y))+
  labs(x="Cattle Use", y="Probability of infection")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) + 
  geom_errorbar(aes(ymin=c(0.02270131, 0.02265756 ),ymax=c(0.04923680, 0.04071538 )), width=0.2)+
  geom_point(col=c("darkgoldenrod2"), size=5) +ylim(0,0.065) #was "darkgoldenrod3"


#full GLMM for snail presence
snailGLMM=glmmTMB(SnailPres~Permanence+
                    log(dist)+ CattleUsePerm+Phase+
                    (1|District/Village/CorrectedWBName), 
                  data=three_phases, family ="binomial", REML=T)

summary(snailGLMM) 
plot(allEffects(snailGLMM))

#identifying snail presence probabilities and confidence intervals
summary(predictorEffect("Phase", snailGLMM))
summary(predictorEffect("CattleUsePerm", snailGLMM))

#compiling that into a dataframe
snailphasedata=data.frame(x=c("Phase 1","Phase 2", "Phase 3"),y=c(0.9894688, 0.8737589, 0.5500603 ))
snailcattledata=data.frame(x=c("Not Permitted", "Permitted"),y=c(0.9339328, 0.8968130))

#graph snail presence data and include confidence intervals
#snail presence for phase data
ggplot(snailphasedata,aes(x,y))+
  labs(x="Phase", y="Probability of snail presence")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) + 
  geom_errorbar(aes(ymin=c(0.9724114, 0.7861854, 0.4116515 ),ymax=c(0.9960232, 0.9287164, 0.6811303)), width=0.2)+
  geom_point(col=c("cornflowerblue"), size=5) +ylim(0,1)


#snail presence for cattle permission data
ggplot(snailcattledata,aes(x,y))+
  labs(x="Cattle Use", y="Probability of snail presence")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20)) + 
  geom_errorbar(aes(ymin=c(0.8723885, 0.8210243 ),ymax=c(0.9669210, 0.9427464 )), width=0.2)+
  geom_point(col=c("cornflowerblue"), size=5) +ylim(0,1)

####Data for table 1####

#Villages
#P1-3
three_phases%>%group_by(District)%>%
  summarize(n_distinct(Village))

three_phases%>%
  summarize(n_distinct(Village))

#P2-3
subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(n_distinct(Village))

subset(three_phases, Phase!="Phase1")%>%
  summarize(n_distinct(Village))


#Waterbodies
#P1-3
three_phases%>%group_by(District)%>%
  summarize(n_distinct(District_Waterbody))

District_WBs=three_phases%>%group_by(District)%>%
  summarize(WB_number=n_distinct(District_Waterbody))

three_phases%>%
  summarize(n_distinct(District_Waterbody))
  
#P2-3
subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(n_distinct(District_Waterbody))

p2p3_WBnumber=subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(WB_number=n_distinct(District_Waterbody))

subset(three_phases, Phase!="Phase1")%>%
  summarize(n_distinct(District_Waterbody))


#schools
schools%>%group_by(District)%>%
  summarize(n())

schools%>%
  summarize(n())

#Bulinus numbers
#p1-3
three_phases%>%group_by(District)%>%
  summarize(sum(BulinusNumber))

three_phases%>%
  summarize(sum(BulinusNumber))

#p2-3
subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(sum(BulinusNumber))

p2p3_snailtotals=subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(snails=sum(BulinusNumber))

subset(three_phases, Phase!="Phase1")%>%
  summarize(sum(BulinusNumber))


#Schistosome infection (P1-3)
#at waterbody level
three_phases%>%group_by(District_Waterbody)%>%
  arrange(desc(SchistoSnails)) %>%
  slice_head(n = 1)%>%
  filter(SchistoSnails > 0)%>%
  group_by(District)%>%
  summarize(n())

three_phases%>%group_by(District_Waterbody)%>%
  arrange(desc(SchistoSnails)) %>%
  slice_head(n = 1)%>%
  filter(SchistoSnails > 0)%>%
  group_by(District)%>%
  summarize(n=n())%>%left_join(District_WBs)%>%
  mutate(n/WB_number*100)

three_phases%>%group_by(District_Waterbody)%>%
  arrange(desc(SchistoSnails)) %>%
  slice_head(n = 1)%>%
  filter(SchistoSnails > 0)%>%
  group_by(District)%>%
  summarize(snailcount=n())%>%
  summarise(sum(snailcount),sum(snailcount)/467*100)


three_phases%>%group_by(District)%>%
  summarize(sum(SchistoSnails))

three_phases%>%group_by(District)%>%
  summarize(sum(SchistoSnails)/sum(BulinusNumber)*100)

three_phases%>%
  summarize(sum(SchistoSnails))

three_phases%>%
  summarize(sum(SchistoSnails)/sum(BulinusNumber)*100)

#Nonschistosome infection (P2-3)
#at waterbody level
subset(three_phases, Phase!="Phase1")%>%group_by(District_Waterbody)%>%
  arrange(desc(NonSchistoSnails)) %>%
  slice_head(n = 1)%>%
  filter(NonSchistoSnails > 0)%>%
  group_by(District)%>%
  summarize(n())

subset(three_phases, Phase!="Phase1")%>%group_by(District_Waterbody)%>%
  arrange(desc(NonSchistoSnails)) %>%
  slice_head(n = 1)%>%
  filter(NonSchistoSnails > 0)%>%
  group_by(District)%>%
  summarize(n=n())%>%left_join(p2p3_WBnumber)%>%
  mutate(n/WB_number*100)

subset(three_phases, Phase!="Phase1")%>%group_by(District_Waterbody)%>%
  arrange(desc(NonSchistoSnails)) %>%
  slice_head(n = 1)%>%
  filter(NonSchistoSnails > 0)%>%
  group_by(District)%>%
  summarize(snailcount=n())%>%
  summarise(sum(snailcount),sum(snailcount)/388*100)


subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(sum(NonSchistoSnails))

subset(three_phases, Phase!="Phase1")%>%group_by(District)%>%
  summarize(sum(NonSchistoSnails)/sum(BulinusNumber)*100)

subset(three_phases, Phase!="Phase1")%>%
  summarize(sum(NonSchistoSnails))

subset(three_phases, Phase!="Phase1")%>%
  summarize(sum(NonSchistoSnails)/sum(BulinusNumber)*100)

#Coinfected waterbodies
subset(three_phases, Phase!="Phase1")%>%
  mutate(coinfected=case_when(SchistoSnails>0 & NonSchistoSnails>0~1))%>%
  group_by(District_Waterbody)%>%
  arrange(desc(SchistoSnails)) %>%
  slice_head(n = 1)%>%
  group_by(District)%>%
  filter(coinfected > 0)%>%
  summarize(snailcount=n())%>%left_join(p2p3_WBnumber)%>%
  mutate(perc=snailcount/WB_number)

subset(three_phases, Phase!="Phase1")%>%
  mutate(coinfected=case_when(SchistoSnails>0 & NonSchistoSnails>0~1))%>%
  group_by(District_Waterbody)%>%
  arrange(desc(SchistoSnails)) %>%
  slice_head(n = 1)%>%
  group_by(District)%>%
  filter(coinfected > 0)%>%
  summarize(snailcount=n())%>%
  summarise(sum(snailcount), sum(snailcount)/388*100)
  
#Coinfected snails
P2P3_bulinus%>%
  mutate(coinfected=case_when(SchistoStatus>0 & NonSchistoStatus>0~1))%>%
  group_by(District)%>%
  filter(coinfected > 0)%>%
  summarize(snailcount=n())

P2P3_bulinus%>%
  mutate(coinfected=case_when(SchistoStatus>0 & NonSchistoStatus>0~1))%>%
  group_by(District)%>%
  filter(coinfected > 0)%>%
  summarize(snailcount=n())%>%left_join(p2p3_snailtotals)%>%
  mutate(perc=snailcount/snails*100)

P2P3_bulinus%>%
  mutate(coinfected=case_when(SchistoStatus>0 & NonSchistoStatus>0~1))%>%
  #group_by(District)%>%
  filter(coinfected > 0)%>%
  summarize(snailcount=n(), snailcount/25052*100)

