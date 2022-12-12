#----------------------------------------------
# PSYC 3010 A - Intermediate Research Methods
# Group Project - Team #2
# Reaction time study - Lawrence
#----------------------------------------------

# Set Working Directory
setwd("~/Desktop/PSYC3010/P3010GP_RFiles/Data")

#-------------------------------------------------------------------------------
# INSTALL PACKAGES/ LOAD LIBRARIES
#-------------------------------------------------------------------------------

install.packages("data.table")
install.packages("vcd")
install.packages("lmerTest")
# install.packages("devtools")
# install the stable version

library(dplyr)
library(ggplot2)
library(rstatix) # only used for get_summary_stats()
library(janitor) # only used for tabyl()
library(performance) # used for testing statistical assumptions
library(emmeans) # used for post-hoc t-tests
library(lsr) # used for effect sizes
library(stringr) # used for str_detect()
library(lme4)
require(lmerTest)
library(lattice)
library(skimr)
library(data.table)
library(performance) # remember to load the library after you install!
library(vcd)
library(tidyr)
library(MBESS)
library(flexplot)


#-------------------------------------------------------------------------------
# DATA CLEAN UP: STEP 1 COMBINE ALL .CSV FILES INTO ONE AND LOAD IT AS MASTER DATA
#-------------------------------------------------------------------------------

# Remove existing "combined.csv" and recombine all .csv files as a master dat.

file.remove("~/Desktop/PSYC3010/P3010GP_RFiles/Data/combined.csv") #This removes existing file, "combined.csv".
files <- list.files(pattern = ".csv")
temp <- lapply(files, fread, sep = ",")
dat <- rbindlist (temp, fill = TRUE )
write.csv(dat, file = "combined.csv", row.names = FALSE)

# Run below code only when you need to reload the combined.csv file as 
# your primary data = "dat"
dat = read.csv("combined.csv", header=TRUE, na.strings=c("", "NA", " "))

#-------------------------------------------------------------------------------
# DATA CLEAN UP: STEP 2 REMOVE PRACTICE ROWS (TWO COLUMNS TO LOOK AT)
#-------------------------------------------------------------------------------

# There are two columns that uniquely identify "practice rows" per participant. 
# Remove practice rows if "practrials.thisRepN" or "practrials.thisRepN" are "0"
dat <- subset(dat,is.na(practise_trials.thisRepN))
dat <- subset(dat,is.na(practrials.thisRepN))


# One important note: 
# If the variable in question has abbreviated word, such as "practrials"
# instead of "practise_trials..." it is indicative that  
# they are concerning a condition for "Lawrence."


# Confirm if the "0s" are removed from these variables
# If NAs are the only items you see in your list below,
# you are safe to move onto the next process.
tabyl(dat$practise_trials.thisRepN)
tabyl(dat$practrials.thisRepN)

# Total number of rows should be down to 3180 from 3330 at this time.
# (minus 150 rows = 15 files x 10 practice per participant)

# Confirm if all raw data files are combined and each participant has the 
# equal number of trials (n=212).
tabyl(dat$participant)

#-------------------------------------------------------------------------------
# DATA CLEAN UP: STEP 3: CONNECT TWO SEPERATE VARIABLES INTO ONE COLUMN 
#-------------------------------------------------------------------------------

#--- Step 3.1: Identify all variables relevant to our analysis ---

# Show all column names and identify all variables necessary.
colnames(dat)

#--- Step 3.2 ---
#    Combine separated variables 5) & 6), 7) & 8), 9) & 10), and 11) & 12)   
#    for ease of computation

# Combine two separated "correct" response keys into one column and call it "correct_keys"
# This column indicates the "correct answers" which you can compare against each response. 
dat$correct_keys <- coalesce(dat$Corr_Resp,dat$Corr_Response)


# Combine two separated response key data into one column and call it "response_keys"
# This one is the respondents' answers in one row.
dat$response_keys <- coalesce(dat$Presence_Resp.keys,dat$Presence_Response.keys)


# Confirm the "coalesce" function above by counting the number of Ms and Zs
# You can already see the discrepancies between the two lines.
tabyl(dat$correct_keys)
tabyl(dat$response_keys)


# Combine two response time variables into one column and call it "comb_RT"
# All participants' response time are in one line.
dat$combined_RT <- coalesce(dat$Presence_Resp.rt, dat$Presence_Response.rt)


#################change here: made this numeric
dat$comb_res_corr <-  coalesce(as.numeric(dat$Presence_Resp.corr),dat$Presence_Response.corr)
### NOTE: THIS VARIABLE IS NOT SUITABLE FOR OUR ANALYSIS DUE TO A FLAW IN THE DATA.

# Confirm the coalesce just to show you the problematic zeros...
dat$comb_res_corr

# What was missing from our data at this point was a proper identification of the conditions 
# (Laurence or Ours) under which each entry was set. Therefore...

# Command R to identify trial conditions between Lawrence and Ours by creating a 
# new column called, "WhichExperiment" and show "Lawrence" when participant's 
# "Presence_Resp.keys" says NA. Additionally, show "OurGroup" when 
# "Presence_Reponse.keys" says NA.

dat <- dat %>% 
  mutate(WhichExperiment = case_when(
    !is.na(Presence_Resp.corr) ~ "Lawrence",
    !is.na(Presence_Response.keys) ~ "OurGroup"))

# confirm the change
dat$WhichExperiment

#-------------------------------------------------------------------------------
# DATA CLEAN UP: STEP 4: COMPARE / CONTRAST & MAKE THEM INTO BINARY 
#-------------------------------------------------------------------------------

# In this step, the two variables are compared and a preliminary step is taken 
# to ascertain the proportion of each participant who answered correctly in 
# all or any of the conditions.

# Create a new binary column called Correct and display the number 1 if 
# 'correct_keys' and 'response_keys' match. Otherwise, display 0. 
dat = dat %>%
  mutate(Correct = case_when(
    correct_keys == response_keys ~ 1,
    TRUE ~ 0))

# Confirm the change
dat$Correct

# Filter out all NAs in the "correct_keys" by commanding R to 
# show results other than NAs included in the correct_keys.
dat = dat %>% filter(!is.na(correct_keys))

#-------------------------------------------------------------------------------
# DATA CLEAN UP: STEP 5: REMOVE PARTICIPANT WHEN ERROR >= 40% 
#-------------------------------------------------------------------------------

# Set a condition so that participants who exceeds 40% in the error rates are 
# identified and excluded from the primary data. 

# In this step, commanded R to sort by participant (ID), greate a new column call
# Percentage_Correct_Participant if their average of "Correct" variable is equal and 
# above 60% (inversely saying that people with 40% or error to be out of our study).

dat = dat %>%
  group_by(participant) %>%
  mutate(Percentage_Correct_Participant = mean(Correct)) %>%
  filter(Percentage_Correct_Participant >= 0.6)

# Confirm if the above screening worked properly
# Show only participant, n, and mean in the summary.

dat%>%
  group_by(participant)%>%
  get_summary_stats(Correct)%>%
  select(participant,n,mean)


#-------------------------------------------------------------------------------
# DATA CLEAN UP: STEP 3 REMOVE IF THEIR SEARCH ERROR EXCEEDS 10% ON EACH CONDITION
#-------------------------------------------------------------------------------

# Screen out the participants who made an error larger than 10% on any given conditions. 

dat = dat %>%
  group_by(participant, Set_Size, DistractorPresent, TargetPresent, WhichExperiment) %>%
  mutate(Percentage_Correct_Condition = mean(Correct)) %>%
  group_by(participant) %>%
  mutate(AllConditionsOkay = case_when(
    Percentage_Correct_Condition > 0.1 ~ 1,
    TRUE ~ 0),
    PercentageOkay = mean(AllConditionsOkay),
    ParticipantOkay = case_when(
      mean(AllConditionsOkay) == 1 ~ "ParticipantOkay",
      mean(AllConditionsOkay) < 1 ~ "ParticipantNotOkay"))

dat$AllConditionsOkay
dat$PercentageOkay

dat%>%
  group_by(participant) %>%
  tabyl(participant,AllConditionsOkay)


# Remove outliers by limiting response time to 5 seconds.
dat = dat %>% filter(combined_RT < 5)


#-------------------------------------------------------------------------------
# LMM: FITTING THE MODEL
#-------------------------------------------------------------------------------

#### Recap: ZOOM with Professor Joerges ####
#### Fitting our experiment into LMM framework ####

Model3 = lmer(combined_RT ~ DistractorPresent*WhichExperiment +
                (Set_Size + TargetPresent + as.factor(Line_Onset) | participant),data = dat)

# Model 3 has "Signularity" issue so remove "Onset" and fit as Model4
Model4 = lmer(combined_RT ~ DistractorPresent*WhichExperiment +
                (Set_Size + TargetPresent | participant),data = dat)

# Model4 continues to get "Singularity" so remove "TargetPresent" and fit as 
# fit as Model 5
Model5 = lmer(combined_RT ~ DistractorPresent*WhichExperiment +
                (Set_Size | participant),data = dat)


# Run a likelihood Ratioin test in anova() function to test if "onset" variable as
# random effect is causing the "SINGULARITY" issue...
anova(Model3,Model4,Model5)

table1 <- anova(Model3,Model4,Model5)
write.csv(table1, file = "table1.csv", row.names = T)


# Selected Model 4 for our analysis (Onset removed from random effect) 
summary(Model4)

Model4_CI = confint(Model4,method = "boot")
confint(Model4)

mod4_ranef <- ranef(Model4)

# Still does not provide "Confidence Intervals". Therefore, remove "Target_Present" 
# from the random effect, simplify the computation and re-fit the model4.
Model6 = lmer(Correct ~ DistractorPresent*WhichExperiment +
                (Set_Size + TargetPresent + as.factor(Line_Onset) | participant),data = dat)

Model7 = lmer(Correct ~ DistractorPresent*WhichExperiment +
                (Set_Size + TargetPresent | participant),data = dat)

Model8 = lmer(Correct ~ DistractorPresent*WhichExperiment +
                (Set_Size | participant),data = dat)

table2 <- anova(Model6,Model7,Model8)
write.csv(table2, file = "table2.csv", row.names = T)


# use Model7: get summary, CIs 
summary(Model7)

#-------------------------------------------------------------------------------
# LMM: DESCRIPTIVE STATISTICS
#-------------------------------------------------------------------------------
# Obrain mean values for RT and participants' Correct answers

dat=dat %>%
  group_by(participant, TargetPresent, 
           DistractorPresent, Set_Size, WhichExperiment)%>%
  mutate (MeanRT = mean(combined_RT, na.rm=T), 
          MeanAccuracy = mean(Correct*100, na.rm=T))

dat$MeanAccuracy
# Show mean score of "Correct" for each participant

# Show results in all conditions
dat%>% 
  get_summary_stats(MeanAccuracy)%>%
  print(n = 208)
 
# Combine results by participant (answers under 5 seconds)
dat%>%
  group_by(participant)%>%
  get_summary_stats(MeanAccuracy) 
  
# Show mean score of "MeanRT" for each participant
dat%>%
  group_by(participant)%>%
  get_summary_stats(MeanRT)

# Change variable type of DistroctorPresent (x axis) to factor for charts. 
dat$DistractorPresent <- as.factor(dat$DistractorPresent)

# Create a box-plot for RT (make sure to use all eligible participants' entries)
ggplot(dat,aes(x = DistractorPresent, y = combined_RT, color=WhichExperiment)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("Absent","Present")) +
  scale_color_discrete(labels = c("Lawrence","Our Study"))+
  ggtitle("Response Time") +
  labs(x="Distractor", y= "Response Time (sec)")


# Create a box-plot for accuracy (make sure to use all eligible participants' entries)
ggplot(dat,aes(x = DistractorPresent, y = MeanAccuracy, color=WhichExperiment)) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Absent","Present")) +
  ggtitle("Response Accuracy") +
  labs(x="Distractor", y= "Mean Accuracy (%)") +
  scale_color_discrete(labels = c("Lawrence","Our Study"))



