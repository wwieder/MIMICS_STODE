#using Derek's litterbag code to run MSBio litter decomp simulations

library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)

source("functions/MIMICS_sim_litterbag.R")
source("functions/MIMICS_calc_steady_state_pools.R")
source("functions/calc_Tpars.R")
source("Parameters/MIMICS_parameters_sandbox_20231129.R")
source("functions/RXEQ.R")


##########
# MSBio runs
#############

#-------------------------------
#Using MSBio data
#-------------------------------
####
#load MSBio site and litter data and format to code structure
####

#load site data
MSBio <- read.csv("Example_simulations/Data/Site_annual_clim.csv")
#match input data strucutre
#AGNPP should be in gDW!! multiply by 2 here to remedy
#don't have gravimetric soil moisture, just volumetric, assuming a BD of 1g/cm3 makes them equivalent - could be bad assumption given this is BD of leaves
MSBio2 <- MSBio %>% mutate(SITE = Site, ANPP = AGNPP_sum*2, TSOI = TSOI_mean, CLAY = PCT_CLAY_mean, lig_N = LIG_N, GWC = H2OSOI_mean*100, W_SCALAR=W_SCALAR_mean) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG, C, N, CN, LIG_N, GWC, W_SCALAR) 
#daily data - change site name and MSBio2 row (1=BART, 8=SERC) to use different site daily input
BART_dailyinput <- read.csv("Example_simulations/Data//BART_clim.csv")

#check sum of daily fluxes
sum(BART_dailyinput$LITFALL)
plot(BART_dailyinput$LITFALL)
lines(BART_dailyinput$AGNPP)
2*sum(BART_dailyinput$AGNPP)
sum(BART_dailyinput$NPP)
sum(BART_dailyinput$LITFALL)
plot(BART_dailyinput$TSOI/10)
lines(BART_dailyinput$W_SCALAR)
lines(BART_dailyinput$H2OSOI)
lines(BART_dailyinput$TSOI/10)

BART_DI <- BART_dailyinput %>% mutate(DAY=X, ANPP = AGNPP*2, CLAY = rep(MSBio2[1,4], 366), LIG_N = rep(MSBio2[1,9], 366), GWC = H2OSOI*100) %>%
  select(DAY, ANPP, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 

# two test dataframes with daily and annual results
BART_DI2 <- BART_dailyinput %>% mutate(SITE= rep("BART", 366),DAY=X, ANPP = rep(sum(LITFALL)*2,366), LITFALL=LITFALL*2, CLAY = rep(MSBio2[1,4], 366), 
                                       LIG_N = rep(MSBio2[1,9], 366), GWC = H2OSOI*100) %>%
  select(SITE, DAY, ANPP, LITFALL, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 

BART_DI2_mean <- BART_DI2  %>% mutate( SITE= "BART", ANPP = mean(ANPP), TSOI = mean(TSOI), CLAY = mean(CLAY), 
                                       LIG_N= mean(LIG_N), GWC=mean(GWC), W_SCALAR=mean(W_SCALAR)) %>%
  select(SITE, ANPP, TSOI, CLAY, LIG_N, GWC, W_SCALAR) 
BART_DI2_mean = BART_DI2_mean[1,]

vmax_DI1 = exp(BART_DI2$TSOI * 0.063 + 5.47) * 0.000008 * BART_DI2$W_SCALAR *10
plot(vmax_DI1)
abline(h=mean(vmax_DI1),lt=2)
vmax_mean = exp(BART_DI2_mean$TSOI * 0.063 + 5.47) * 0.000008 * BART_DI2_mean$W_SCALAR *10
abline(h=vmax_mean,lw=2)
#Option 1: MSBio litter bags with just variation in NEON litter (not separated by species)
MSBio_BAGS <- read.csv("Example_simulations/Data/NEON_MSB_LitVars.csv")

#Option2: MSBio litter bags with leaf and litter chemistry combined
# MSBio_BAGS <- read.csv("NEON_MSB_LeafChem.csv")
# MSBio_BAGS <- MSBio_BAGS[,2:8]
# #rename to match input
# MSBio_BAGS2 <- MSBio_BAGS %>% mutate(Site = siteID, TYPE = taxonID, BAG_LIG = leaflig, BAG_N = leafN, BAG_CN = leafCN) %>%
#   select(Site, TYPE, BAG_LIG, BAG_N, BAG_CN)
# #add combined chem
# COMBO_BAGS <- data.frame(Site = c("BART", 'GRSM', 'HARV', 'LENO', 'MLBS', 'OSBS', 'SCBI', 'SERC', 'TALL', 'TREE', 'UNDE'),
#                          TYPE = rep("COMBO", 11),
#                          BAG_LIG = MSBio2$LIG,
#                          BAG_N = MSBio2$N,
#                          BAG_CN = MSBio2$CN)
# MSBio_BAGS3 <- rbind(MSBio_BAGS2, COMBO_BAGS)


### changed to fMET calculation in STODE script here!! Note that the two options are only somewhat related but less negatives in STODE equation

MSBio_BAGS$CALC_MET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * (MSBio_BAGS$BAG_LIG/MSBio_BAGS$BAG_N))
MSBio_BAGS$CALC_MET[MSBio_BAGS$CALC_MET <0] = 0 #setting negatives to zero - might want to reconsider this for future runs, all strucutral seems pretty unlikely
#MSBio_BAGS$CALC_N <- (1 / MSBio_BAGS$BAG_CN) / 2.5 * 100 #why calculating from CN and not N directly?
#MSBio_BAGS$CALC_MET2 <- 0.85 - 0.013 * MSBio_BAGS$BAG_LIG/MSBio_BAGS$CALC_N #calculate fMET

BAG_init_size <- 100
BAGS <- MSBio_BAGS %>% select(Site, TYPE, CALC_MET)
BAGS$BAG_LITm <- ((BAG_init_size * 1e3 / 1e4)/ depth) * BAGS$CALC_MET
BAGS$BAG_LITs <- ((BAG_init_size * 1e3 / 1e4)/ depth) * (1-BAGS$CALC_MET) #initial litter = 0.1 because of unit conversions here 


####
#run litterbag model 
####

#Individual litters decomposing too fast
#try N re-absorption
#try min and max of litter COMBO values
#compare to litter mass loss in observations

#SERC and all SERC litters (1) for 2 years
BAGS_SERC <- filter(BAGS, Site == "SERC" & TYPE == "mean")
BAGS_SERC <- BAGS_SERC[,2:5]
BAGS_out_SERC_SS <- BAGS_SERC %>% split(1:nrow(BAGS_SERC)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                           forcing_df=MSBio2[8,],
                                                                           #dailyInput = SERC_DI, 
                                                                           nspin_yrs=5,
                                                                           nspin_days=0,
                                                                           litadd_day=10,
                                                                           verbose=T)) %>% bind_rows()


#LENO and all LENO litters for 2 years
# BAGS_LENO <- filter(BAGS, Site == "LENO")
# BAGS_LENO <- BAGS_LENO[,2:5]
# MSBio_LENO_lsm <- MSBio2
# MSBio_LENO_lsm[4, 10] = 32 #from empirical data - LCI
# BAGS_out <- BAGS_LENO %>% split(1:nrow(BAGS_LENO)) %>% map(~ MIMICS_LITBAG(litBAG=.,
#                                                                            forcing_df=MSBio2[4,],
#                                                                            dailyInput = LENO_DI,
#                                                                            nspin_yrs=2,
#                                                                            nspin_days=0,
#                                                                            litadd_day=10,
#                                                                            verbose=T)) %>% bind_rows()

#UNDE and all UNDE litters (11) for 2 years
# BAGS_UNDE <- filter(BAGS, Site == "UNDE")
# BAGS_UNDE <- BAGS_UNDE[,2:5]
# BAGS_out <- BAGS_UNDE %>% split(1:nrow(BAGS_UNDE)) %>% map(~ MIMICS_LITBAG(litBAG=.,
#                                                                            forcing_df=MSBio2[11,],
#                                                                            nspin_yrs=2,
#                                                                            nspin_days=0,
#                                                                            litadd_day=10,
#                                                                            verbose=T)) %>% bind_rows()

#BART and all BART litters (1) for 2 years
BAGS_BART <- filter(BAGS, Site == "BART" & TYPE == "mean")
BAGS_BART <- BAGS_BART[,2:5]
BAGS_out_BART_SS <- BAGS_BART %>% split(1:nrow(BAGS_BART)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                           forcing_df=MSBio2[1,],
                                                                           #dailyInput = BART_DI,
                                                                           nspin_yrs=5,
                                                                           nspin_days=0,
                                                                           litadd_day=10,
                                                                           verbose=T)) %>% bind_rows()
ss_testM = MIMICS_SS(BART_DI2_mean[1,])
ss_testM[[2]]$VMAX
ss_testM[[2]]$KM
ss_testM[[2]]$I

ss_test = MIMICS_SS(BART_DI2)
ss_test[[2]]$VMAX
ss_test[[2]]$KM
ss_test[[2]]$I

ss_testM[[1]] - ss_test[[1]]
BAGS_out_BART_day1 <- BAGS_BART %>% split(1:nrow(BAGS_BART)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                                   forcing_df=BART_DI2,
                                                                                   dailyInput = BART_DI2,
                                                                                   nspin_yrs=5,
                                                                                   nspin_days=0,
                                                                                   litadd_day=10,
                                                                                   verbose=T)) %>% bind_rows()

# WRW tests for new code:
BAGS_BART2 <- filter(BAGS, Site == "BART" & TYPE == "mean")
BAGS_BART2 <- BAGS_BART2[,2:5]
BAGS_out_BART_daily <- BAGS_BART2 %>% split(1:nrow(BAGS_BART2)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                                   forcing_df=BART_DI2_mean[1,],
                                                                                   dailyInput = BART_DI2,
                                                                                   nspin_yrs=5,
                                                                                   nspin_days=0,
                                                                                   litadd_day=10,
                                                                                   verbose=T)) %>% bind_rows()

BAGS_out_BART_SS2 <- BAGS_BART2 %>% split(1:nrow(BAGS_BART2)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                                        forcing_df=BART_DI2_mean[1,],
                                                                                        #dailyInput = BART_DI2,
                                                                                        nspin_yrs=5,
                                                                                        nspin_days=0,
                                                                                        litadd_day=10,
                                                                                        verbose=T)) %>% bind_rows()
#all sites and all litters
# BAGS_mean <- filter(BAGS, TYPE == "mean")
# BAGS_input <- split(BAGS_mean, 1:nrow(BAGS_mean))
# forcing_input <- split(MSBio2, 1:nrow(MSBio2))
# BAGS_out_AllSites <- map2(forcing_input, BAGS_input, ~MIMICS_LITBAG(.x, .y, nspin_yrs=2, nspin_days=0, litadd_day=10, verbose=T)) %>% bind_rows()

####
#plot output
####

colorBlind7  <- c("#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #yellow (LENO), blue (SERC), green (UNDE)

#Formating observational data for comparing to field litter mass loss
Field_LML <- read.csv("Example_simulations/Data/Litter_decomp_all.csv")
#Add Species to group_by to get species-specific summary
LML_sum2 <- Field_LML %>% filter(site == 'BART') %>% group_by(time.point) %>% drop_na(percent.loss.litter) %>% summarize(mean.ML = mean(percent.loss.litter/2),
                                                                                                n = n(),
                                                                                                sd = sd(percent.loss.litter/2),
                                                                                                SE = sd/sqrt(n),
                                                                                                lci.ML = mean.ML - qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                                                uci.ML = mean.ML + qt(1 - ((1 - 0.95) / 2), n - 1) * SE,
                                                                                                doy = mean(days_elapsed))

###
#moisture function testing
###

#wide format MIMICS output for plotting
# BAGS_out_wide_fWmeth0 = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fWmeth1 = BAGS_out_fWm1 %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fWmeth2 = BAGS_out_fWm2 %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fWmeth3 = BAGS_out_fWm3 %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# 
# #wide format for fW effects on just Vmax or on both Vmax and tau
# BAGS_out_wide = BAGS_out_100y %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_fW.tau = BAGS_out_fWt %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# 
# #plot MIMICS output with different soil moisture and field litter mass loss together
# ggplot() +
#   #geom_line(data=BAGS_out_wide_fWmeth0, aes(y=(mean/0.1)*100, x=DAY, color ="fW=1"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fWmeth1, aes(y=(mean/0.1)*100, x=DAY, color ="CORPSE"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fWmeth2, aes(y=(mean/0.1)*100, x=DAY, color ="Calibrated"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fWmeth3, aes(y=(mean/0.1)*100, x=DAY, color ="W_SCALAR"), linewidth=2, alpha=0.5) +
#   geom_line(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, color ="Default"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_wide_fW.tau, aes(y=(mean/0.1)*100, x=DAY, color ="Tau Effects"), linewidth=2, alpha=0.5) +
#   #geom_ribbon(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, ymin = (lci/0.1)*100, ymax=(uci/0.1)*100), alpha = 0.3) +
#   #geom_point(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10), color = "#009E73", size = 3) +
#   #geom_errorbar(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10, ymin = (1-lci.ML)*100, ymax = (1-uci.ML)*100), width=0, color = "#009E73",linewidth=1) +
#   ylab("Litter Bag C Remaining (%)") +
#   xlab("Day") +
#   scale_color_manual(values = c("fW=1"="#E69F00", "CORPSE"="#56B4E9", "Calibrated"="#009E73", "W_SCALAR" = "#F0E442")) +
#   ggtitle("NPP option for turnover") +
#   theme_bw(base_size = 20)
# 
# 
# #plot MIMICS output with different soil moisture and microbial dynamics together
# ggplot() +
#   geom_line(data=BAGS_out_SERC_fW0, aes(y=MICr+MICk, x=DAY, color ="fW=1"), linewidth=2, alpha=0.3) +
#   geom_line(data=BAGS_out_SERC_fW1, aes(y=MICr+MICk, x=DAY, color ="CORPSE"), linewidth=2, alpha=0.3) +
#   geom_line(data=BAGS_out_SERC_fW2, aes(y=MICr+MICk, x=DAY, color ="Calibrated"), linewidth=2, alpha=0.3) +
#   geom_line(data=BAGS_out_SERC_fW3, aes(y=MICr+MICk, x=DAY, color ="W_SCALAR"), linewidth=2, alpha=0.3) +
#   #geom_line(data=BAGS_out, aes(y=MICr+MICk, x=DAY, color ="Default"), linewidth=2, alpha=0.5) +
#   #geom_line(data=BAGS_out_fWt, aes(y=MICr+MICk, x=DAY, color ="Tau Effects"), linewidth=2, alpha=0.5) +
#   ylab("microbial biomass") +
#   xlab("Day") +
#   #ylim(0,2) +
#   xlim(0,3650)+
#   scale_color_manual(values = c("fW=1"="#E69F00", "CORPSE"="#56B4E9", "Calibrated"="#009E73", "W_SCALAR" = "#F0E442")) +
#   ggtitle("beta option for turnover - SERC") +
#   theme_bw(base_size = 20)
# 
# #Mic biomass vs soil moisture
# BO_fWm0 <- BAGS_out_fWm0 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# BO_fWm1 <- BAGS_out_fWm1 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# BO_fWm2 <- BAGS_out_fWm2 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# BO_fWm3 <- BAGS_out_fWm3 %>% filter(DAY<366) %>% filter(Litter_Type=="mean") %>% left_join(SERC_DI, by="DAY")
# ggplot() +
#   geom_line(data=BO_fWm0, aes(y=MICr+MICk, x=GWC, color = "fW=1"), linewidth=2, alpha=0.5) +
#   geom_line(data=BO_fWm1, aes(y=MICr+MICk, x=GWC, color = "CORPSE"), linewidth=2, alpha=0.5) +
#   geom_line(data=BO_fWm2, aes(y=MICr+MICk, x=GWC, color = "Calibrated"), linewidth=2, alpha=0.5)+
#   geom_line(data=BO_fWm3, aes(y=MICr+MICk, x=GWC, color = "W_SCALAR"), linewidth=2, alpha=0.5)+
#   ylim(0,0.7) +
#   ylab("Microbial biomass") +
#   xlab("GWC") +
#   scale_linetype_manual(values = c("fW=1"="E69F00", "CORPSE"="#56B4E9", "Calibrated"="#009E73", "W_SCALAR" = "#F0E442"))


###
#Daily vs. STODE
####
ggplot() +
  geom_line(data=BAGS_out_BART_day1, aes(y=LITm+LITs, x=DAY/365, color = "Litter", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_day1, aes(y=MICr+MICk, x=DAY/365, color = "MIC", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_day1, aes(y=SOMa+SOMc+SOMp, x=DAY/365, color = "SOM", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=MICr+MICk, x=DAY/365, color = "MIC", linetype ="steady state"), linewidth=2, alpha=0.5,) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=LITm+LITs, x=DAY/365, color = "Litter", linetype ="steady state"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=SOMa+SOMc+SOMp, x=DAY/365, color = "SOM", linetype ="steady state"), linewidth=2, alpha=0.5) +
  scale_color_manual(values = c("MIC"="#E69F00", "Litter"="#56B4E9", "SOM"="#009E73", "W_SCALAR" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("Year") +
  ggtitle("W_SCALAR moisture, beta at BART") +
  theme_bw(base_size = 20)

# SERC_daily <- rbind(SERC_DI, SERC_DI)
# SERC_daily$DAY <- 1:732
# BAGS_daily <- inner_join(BAGS_out_SERC_lowANPP, SERC_daily, by='DAY')
# BAGS_out_4y <- filter(BAGS_out_BART, DAY <1461)
# BAGS_BART_sum <- BAGS_out_BART %>% mutate(YEAR = c(rep(1, 365), rep(2, 365), rep(3, 365), rep(4, 365), rep(5, 365))) %>% group_by(YEAR) %>%
#   summarise(daily_mean_mic = mean(MICr+MICk), daily_mean_lit = mean(LITm+LITs), daily_mean_som = mean(SOMa+SOMc+SOMp))
# BAGS_out_BART <- BAGS_out_BART %>% mutate(YEAR = c(rep(1, 365), rep(2, 365), rep(3, 365), rep(4, 365), rep(5, 365))) %>%
#   inner_join(BAGS_BART_sum)
ggplot() +
  geom_line(data=BAGS_out_BART_daily, aes(y=LITm+LITs, x=DAY/365, color = "Litter", linetype ="daily"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=MICr+MICk, x=DAY/365, color = "MIC", linetype ="daily"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=SOMa+SOMc+SOMp, x=DAY/365, color = "SOM", linetype ="daily"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=MICr+MICk, x=DAY/365, color = "MIC", linetype ="steady state"), linewidth=2, alpha=0.5,) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=LITm+LITs, x=DAY/365, color = "Litter", linetype ="steady state"), linewidth=2, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=SOMa+SOMc+SOMp, x=DAY/365, color = "SOM", linetype ="steady state"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_out_BART, aes(y=daily_mean_mic, x=DAY, color = "MIC", linetype ="daily mean"), linewidth=2, alpha=0.5,) +
  #geom_line(data=BAGS_out_BART, aes(y=daily_mean_lit, x=DAY, color = "Litter", linetype ="daily mean"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_out_BART, aes(y=daily_mean_som, x=DAY, color = "SOM", linetype ="daily mean"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_daily, aes(y=ANPP/2, x=DAY, color = "ANPP"), linewidth=2, alpha=0.5) +
  #geom_line(data=BAGS_daily, aes(y=W_SCALAR, x=DAY, color = "W_SCALAR"), linewidth=2, alpha=0.5) +
  #scale_y_log10() +
  scale_color_manual(values = c("MIC"="#E69F00", "Litter"="#56B4E9", "SOM"="#009E73", "W_SCALAR" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("Year") +
  ggtitle("W_SCALAR moisture, NPP at BART") +
  theme_bw(base_size = 20)


ggplot() +
  geom_line(data=BAGS_out_BART_daily, aes(y=LITm, x=DAY/365, color = "LITm", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=LITs, x=DAY/365, color = "LITs", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=MICr, x=DAY/365, color = "MICr", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=MICk, x=DAY/365, color = "MICk", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=LITm, x=DAY/365, color = "LITm", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=LITs, x=DAY/365, color = "LITs", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=MICr, x=DAY/365, color = "MICr", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=MICk, x=DAY/365, color = "MICk", linetype ="daily"), linewidth=2, alpha=0.9) +
  scale_y_log10() +
  scale_color_manual(values = c("MICr"="#E69F00", "MICk"="#56B4E9", "LITm"="#009E73", "LITs" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("year") +
  ggtitle("W_SCALAR moisture, NPP at BART") +
  theme_bw(base_size = 20)

ggplot() +
  geom_line(data=BAGS_out_BART_daily, aes(y=SOMa, x=DAY/365, color = "SOMa", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=SOMc, x=DAY/365, color = "SOMc", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_daily, aes(y=SOMp, x=DAY/365, color = "SOMp", linetype ="daily"), linewidth=1, alpha=0.5) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=SOMa, x=DAY/365, color = "SOMa", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=SOMc, x=DAY/365, color = "SOMc", linetype ="daily"), linewidth=2, alpha=0.9) +
  geom_line(data=BAGS_out_BART_SS2, aes(y=SOMp, x=DAY/365, color = "SOMp", linetype ="daily"), linewidth=2, alpha=0.9) +
  scale_color_manual(values = c("SOMa"="#E69F00", "SOMc"="#56B4E9", "SOMp"="#009E73", "LITs" = "#F0E442", "ANPP" = "#CC79A7")) +
  scale_linetype_manual(values = c("daily"=1, "steady state"=2, "daily mean"=3)) + 
  ylab("C pools") +
  xlab("Year") +
  ggtitle("W_SCALAR moisture, NPP at BART") +
  theme_bw(base_size = 20)

###
#within site testing
####


#can varying LQ get the same variability as at the sites?
# ggplot() +
#   geom_line(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY), linewidth=1, alpha=0.5, color = "#E69F00", linetype =1) +
#   geom_ribbon(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, ymin = (lci/0.1)*100, ymax=(uci/0.1)*100), alpha = 0.3) +
#   geom_point(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10), color = "#E69F00", size = 3) +
#   geom_errorbar(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10, ymin = (1-lci.ML)*100, ymax = (1-uci.ML)*100), width=0, color = "#E69F00",linewidth=1) +
#   ylab("Litter Bag C Remaining (%)") +
#   xlab("Day") +
#   labs(linetype="Parameter Set") +
#   ggtitle("LQ induced variablity - LENO") +
#   theme_bw(base_size = 20)
# 
# #can varying soil moisture get the same variability as at the sites?
# #wide format MIMICS output for plotting
# BAGS_out_wide = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_wet = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# BAGS_out_wide_dry = BAGS_out %>% mutate(LITBAG = LITBAGm+LITBAGs) %>% select(SITE, Litter_Type, DAY, LITBAG) %>% pivot_wider(names_from = Litter_Type, values_from = LITBAG)
# ggplot() +
#   geom_line(data=BAGS_out_wide_dry, aes(y=(mean/0.1)*100, x=DAY), linewidth=1.5, alpha=0.5, color = "#E69F00", linetype =1) +
#   #geom_line(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY), linewidth=1.5, alpha=0.5, color = "#E69F00", linetype =2) +
#   geom_line(data=BAGS_out_wide_mid, aes(y=(mean/0.1)*100, x=DAY), linewidth=1.5, alpha=0.5, color = "#E69F00", linetype =3) +
#   #geom_ribbon(data=BAGS_out_wide, aes(y=(mean/0.1)*100, x=DAY, ymin = (lci/0.1)*100, ymax=(uci/0.1)*100), alpha = 0.3) +
#   geom_point(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10), color = "#E69F00", size = 3) +
#   geom_errorbar(data=LML_sum2, aes(y=(1-mean.ML)*100, x=doy+10, ymin = (1-lci.ML)*100, ymax = (1-uci.ML)*100), width=0, color = "#E69F00",linewidth=1) +
#   ylab("Litter Bag C Remaining (%)") +
#   xlab("Day") +
#   labs(linetype="Parameter Set") +
#   ggtitle("fW=0.05 (solid); fW = 0.5 (dotted)- LENO") +
#   theme_bw(base_size = 20)

####
#plot output - comparing observed to MIMICS microbial community
####

#empirical microbe data
# MSBio_rK <- read.csv("MSBio_rK.csv")
# SERC_rK <- MSBio_rK %>% filter(site == "SERC" & time.point == 0) %>% mutate(rK = r/K) %>% mutate(CO = Copiotroph/Oligotroph)
# 
# #bringing data into one dataframe
# a <- data.frame(group = "model", value = (BAGS_out$MICr/BAGS_out$MICk))
# b <- data.frame(group = "obs_rK", value = SERC_rK$rK)
# c <- data.frame(group = "obs_CO", value = SERC_rK$CO)
# plot.data <- rbind(a,b,c)
# 
# #how does r:K change over time? very little...
# ggplot() + 
#   geom_line(data = BAGS_out, aes(x=DAY, y=(MICr/MICk)))
# #comparing r:K values in empirical and in MIMICS
# ggplot(plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot() +
#   xlab("Type of data") + ylab("r:K or C:O") + ggtitle("SERC Microbial Community Comparison") + 
#   theme_bw(base_size = 20) + theme(legend.position="none")
# #empirical data not directly related to LQ and related to moisture a little bit if you use r:K
# 
# #rwa
# df$SITE.rn = paste(df$SITE, df$run_num, sep = "")
# LIT_init <- df %>% filter(DAY == 10) %>% mutate(LITi = LITBAG_tot) %>% select(SITE.rn, LITi) 
# boxplot(LIT_init$LITi)
# df <- df %>% left_join(LIT_init, by = "SITE.rn") 
# df_730 <- df %>% filter(DAY==730) %>% mutate(LIT_PerLoss = ((LITi - LITBAG_tot)/LITi)*100) #(sample - recovered)/sample *100
# boxplot(df_730$LIT_PerLoss) #looks reasonable

