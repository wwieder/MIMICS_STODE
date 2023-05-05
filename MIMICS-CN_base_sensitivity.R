## Set working drive
rm(list = ls())
setwd("/Users/wwieder/Will/git_repos_local/MIMICS_STODE")

#Libraries
library(rootSolve)
library(boot)
library(dplyr)
library(purrr)
library(ggplot2)
library(Metrics)
library(deSolve)
library(tidyverse)

#bring in RXEQ function
source("CN_RXEQ.R")
source("calc_Tpars.R")


###########################################
# MIMICS single point function
###########################################
experiment = c('Control',
               'Clay=5','Clay=55',
               'CN=25,LIG=10',
               'CN=75,LIG=30')
for (e in 1:5)  {   # loop over exudation experiments
  data = read.csv("LTER_SITE_1.csv")
  df <- data[6,] #6=HAR, #14=LUQ
  
  Site = 'Temperate Forest'
  ANPP = df$ANPP 
  TSOI = df$MAT
  CLAY = df$CLAY2
  LIG = df$LIG
  CN = df$CN
  exud = 0.
  x=1
  ############################################################
  # select range of values to plot
  ############################################################
  CLAY=25
  CN=50
  LIG=20
  
  if (e==2) {CLAY=5}
  if (e==3) {CLAY=55}
  if (e==4) { CN=25
    LIG=10
  }
  if (e==5) {CN=75
    LIG=30
  }

  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud) #>> Same example site input as used for stode
  Tpars
  
  #----------initialize pools---------------
  LIT_1  <<- 1e-4
  LIT_2  <<- 1e-4
  MIC_1  <<- 1e-4
  MIC_2  <<- 1e-4
  SOM_1  <<- 1e-4
  SOM_2  <<- 1e-4
  SOM_3  <<- 1e-4
  
  LIT_1_N  <<- 1e-4
  LIT_2_N  <<- 1e-4
  MIC_1_N  <<- 1e-4
  MIC_2_N  <<- 1e-4
  SOM_1_N  <<- 1e-4
  SOM_2_N  <<- 1e-4
  SOM_3_N  <<- 1e-4
  DIN      <<- 1e-4
  LeachingLoss      <<- 1e-4
  
  pools = c('LIT_m',  'LIT_s',  'MIC_r',  'MIC_K',  'SOM_p',  'SOM_c',  'SOM_a',
            'LIT_m_N','LIT_s_N','MIC_r_N','MIC_K_N','SOM_p_N','SOM_c_N','SOM_a_N',
            'DIN')
  
  # Update MIMICS pools
  #---------------------------------------------
  
  Ty    <<- c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
              MIC_1 = MIC_1, MIC_2 = MIC_2, 
              SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3,
              LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N, 
              MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N, 
              SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
              DIN = DIN)
  
  test  <<- stode(y = Ty, time = 1e7, fun = CN_RXEQ, parms = Tpars, positive = TRUE)
  
  # save results to dataframe
  if (e == 1){
    test1 = test
    df_ss = data.frame(test) * MICROtoECO
    df_ss = data.frame(df_ss,pools)
    df_ss$experiment <- experiment[e]
  } else {
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[e]
    df_ss <- rbind(df_ss, df_exp)
  }    
  
}  # close experiment loop

#====================================
## Finished steady state calculations
## plot results
#====================================

df_ss$pools<- factor(df_ss$pools, levels = pools)
df_ss$experiment<- factor(df_ss$experiment, levels = experiment)

df_ss%>%
  filter(as.numeric(pools) <= 7) %>%
  ggplot(aes(experiment, y, fill = pools)) +
  geom_col() +
  labs(y = "Total C stocks (gC/m2)",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

df_ss %>%
  filter(as.numeric(pools) <= 14) %>% # exclude DIN
  mutate(CNPool = ifelse(grepl("_N", pools), "nitrogen", "carbon")) %>%
  group_by(experiment, CNPool) %>%
  summarize(TotPool = sum(y)) %>% ungroup() %>%
  pivot_wider(values_from = "TotPool", names_from = CNPool) %>%
  mutate(CNRatio = carbon/nitrogen) %>%
  rename(Experiment = experiment) %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  ggplot(aes(x=experiment, y=CNRatio,fill=experiment,show.legend = FALSE)) + 
  geom_bar(stat = "identity",show.legend = FALSE) + 
  coord_cartesian(ylim=c(6,12)) +
  labs(y = "Bulk C:N",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))




###########################################
# MIMICS single point function
###########################################
experiment = c('high clay, high quality',
               'control',
               'low clay, low quality')
for (e in 1:3)  {   # loop over experiments
  data = read.csv("LTER_SITE_1.csv")
  df <- data[6,] #6=HAR, #14=LUQ
  
  Site = 'Temperate Forest'
  ANPP = df$ANPP 
  TSOI = df$MAT
  CLAY = df$CLAY2
  LIG = df$LIG
  CN = df$CN
  exud = 0.
  x=1
  ############################################################
  # select range of values to plot
  ############################################################
  CLAY=25
  CN=50
  LIG=20
  
  if (e==1) {
    CLAY=55
    CN=25
    LIG=10
  }
  if (e==3) {
    CLAY=5
    CN=75
    LIG=30
  }

  Tpars = calc_Tpars(TSOI = TSOI, ANPP = ANPP, CLAY = CLAY, CN =CN, LIG = LIG,
                     x=x,exud=exud) #>> Same example site input as used for stode
  Tpars
  
  #----------initialize pools---------------
  LIT_1  <<- 1e-4
  LIT_2  <<- 1e-4
  MIC_1  <<- 1e-4
  MIC_2  <<- 1e-4
  SOM_1  <<- 1e-4
  SOM_2  <<- 1e-4
  SOM_3  <<- 1e-4
  
  LIT_1_N  <<- 1e-4
  LIT_2_N  <<- 1e-4
  MIC_1_N  <<- 1e-4
  MIC_2_N  <<- 1e-4
  SOM_1_N  <<- 1e-4
  SOM_2_N  <<- 1e-4
  SOM_3_N  <<- 1e-4
  DIN      <<- 1e-4
  LeachingLoss      <<- 1e-4
  
  pools = c('LIT_m',  'LIT_s',  'MIC_r',  'MIC_K',  'SOM_p',  'SOM_c',  'SOM_a',
            'LIT_m_N','LIT_s_N','MIC_r_N','MIC_K_N','SOM_p_N','SOM_c_N','SOM_a_N',
            'DIN')
  
  # Update MIMICS pools
  #---------------------------------------------
  
  Ty    <<- c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
              MIC_1 = MIC_1, MIC_2 = MIC_2, 
              SOM_1 = SOM_1, SOM_2 = SOM_2, SOM_3 = SOM_3,
              LIT_1_N = LIT_1_N, LIT_2_N = LIT_2_N, 
              MIC_1_N = MIC_1_N, MIC_2_N = MIC_2_N, 
              SOM_1_N = SOM_1_N, SOM_2_N = SOM_2_N, SOM_3_N = SOM_3_N,
              DIN = DIN)
  
  test  <<- stode(y = Ty, time = 1e7, fun = CN_RXEQ, parms = Tpars, positive = TRUE)
  
  # save results to dataframe
  if (e == 1){
    test1 = test
    df_ss = data.frame(test) * MICROtoECO
    df_ss = data.frame(df_ss,pools)
    df_ss$experiment <- experiment[e]
  } else {
    df_exp = data.frame(test) * MICROtoECO
    df_exp = data.frame(df_exp,pools)
    df_exp$experiment <- experiment[e]
    df_ss <- rbind(df_ss, df_exp)
  }    
  
}  # close experiment loop

#====================================
## Finished steady state calculations
## plot results
#====================================

df_ss$pools<- factor(df_ss$pools, levels = pools)
df_ss$experiment<- factor(df_ss$experiment, levels = experiment)

df_ss%>%
  filter(as.numeric(pools) <= 7) %>%
  ggplot(aes(experiment, y, fill = pools)) +
  geom_col() +
  labs(y = "Total C stocks (gC/m2)",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

df_ss %>%
  filter(as.numeric(pools) <= 14) %>% # exclude DIN
  mutate(CNPool = ifelse(grepl("_N", pools), "nitrogen", "carbon")) %>%
  group_by(experiment, CNPool) %>%
  summarize(TotPool = sum(y)) %>% ungroup() %>%
  pivot_wider(values_from = "TotPool", names_from = CNPool) %>%
  mutate(CNRatio = carbon/nitrogen) %>%
  rename(Experiment = experiment) %>%
  mutate(Experiment = factor(Experiment, levels = experiment)) %>% # fix duplicate names
  rename(experiment = Experiment) %>%
  ggplot(aes(x=experiment, y=CNRatio,fill=experiment,show.legend = FALSE)) + 
  geom_bar(stat = "identity",show.legend = FALSE) + 
  coord_cartesian(ylim=c(6,12)) +
  labs(y = "Bulk C:N",
       title = paste("Control (Clay=30,CN=50,LIG=20)",Site))

