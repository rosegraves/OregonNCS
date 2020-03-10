# -------------------------------------------
#  R Scripts used to evaluate Natural Climate Solutions for the State of Oregon
#  See Graves et al. 2020 [[insert DOI link]]
#  
#
# -------------------------------------------
# Created by: Rose A. Graves
# The Nature Conservancy, Oregon
# rose.graves@tnc.org 
#  
#
# Modified from R Scripts created by David Marvin (dave@salo.ai) for use in Cameron et al. 2017 (https://doi.org/10.1073/pnas.1707811114)
#

#----------------------------------------------
# Load packages ####
library(tidyverse)
library(ggplot2)
library(ggforce)
library(grid)
library(truncnorm)

# Set the current working directory
setwd("C:/Users/rograves/Documents/TNC/Box Sync/ORFO Climate Science & Policy HET/ORFO CC Science/NCS Pathways Analysis/NCS_AnalysisFiles/OR_NCS_R_Code_PLEASE DO NOT DISTRIBUTE")
setwd('C:/NCSanalysis/') # Make sure this is correct for your computer!
getwd() # Show current working directory

# Source the functions needed for simulations
source('Oregon_NCS_simulation_functions.R') # 

# ----------------------------------------------

# Set up Scenario Simulation Parameters 

# Set years for scenarios
yr=seq(2020,2050,1) # Run simulation from 2020 to 2050, annual timestep
years = seq(1,31,1) # 
yrs=length(years)

# set number of monte carlo iterations
iterations = 1000

# set number of scenarios to run
imp_scenarios = 3  #
scenario.names <- c("Limited","Moderate","Ambitious") # Name scenarios

# Set discount rates for NCS activities
discount_rate_forest = 1 - 0.05 # 5% discount rate based on PNW forest research (Latta et al, Diaz et al) 
discount_rate_reforest = 1 - 0.075 # high mortality in riparian forest reforestation, 1.5x the discount rate for forests ~7.5% annual discount rate
discount_rate_other = 1 - 0.01 # 1% background mortality and reversal of other NCS pathways

# Before starting the next section, be sure that all input files and scenario definitons are correct #

# -----------------------------------------------
# Begin Natural Climate Solution Scenarios

# 1. Deferred Timber Harvest (IFM) -------------
# Load Timber Harvest Input Files
harvest_inputs <- read.csv("R_Inputs/timberharvest_inputs.csv") # input file containing timber harvest data
harvest_carbonstocks <- harvest_inputs[,c(2:3,6:7)] # carbon stock data, stock mean = CO2e per m3 of harvested timber
harvest_carbonstocks$activity<-as.factor(harvest_carbonstocks$activity) # harvest 'activity' refers to timber harvest by ownership and geography
harvest_carbonstocks$type<-as.factor(harvest_carbonstocks$type) # type: avoided conversion versus sequestration
harvest_carbonrates <- harvest_inputs[,c(2:3,4:5)] # carbon sequestration rates, only applies to even-age managed ownership
harvest_carbonrates$activity<-as.factor(harvest_carbonrates$activity)
harvest_carbonrates$type<-as.factor(harvest_carbonrates$type)

# Set up scenario implementation targets

# Timber harvest scenarios: narrative description
# All scenarios include a linear decrease to the target rate by 2030; Supplemental Table S2 provides values for each ownership
# Limited Implementation = decrease timber harvest by the historical variation (CV)
# Moderate Implementation = Gradual reduction of timber harvest to target reduction of 75% on most ownerships and 15% on State and Private Industrial ownerships by 2030. Allows 73% of current harvest volume (overall) from 2030 - 2050
# Ambitious Implementation = Gradual reduction of timber harvest to target reduction of 100% on most ownerships and ~20% on State and Private Industrial ownerships by 2030. Retains 60% of current harvest volume (overall) from 2030 - 2050

harvest_scenariotargets <- harvest_inputs[,c(2,8:12)] # timber harvest target reduction (%)

# build arrays to hold the annual implementation rates for scenario simulations 

harvest_imp <- array(NA,dim=c(length(harvest_inputs$activity),length(years),imp_scenarios), dimnames=list(harvest_inputs$activity,years,scenario.names))

# fill in annual implementation rates for scenario #1: limited implementation
for (ynum in 1:9){  # linear ramp up period of 10 years to reach the historical variation reduction
  harvest_imp[,ynum,1] <- ((harvest_scenariotargets$limited_target*harvest_scenariotargets$baseline)/10)*ynum
}
harvest_imp[,10:31,1] <- harvest_scenariotargets$limited_target*harvest_scenariotargets$baseline # years 10 - 31

# fill in annual implementation rates for scenario #2:moderate implementation
for (ynum in 1:10){ 
  harvest_imp[,ynum,2] <- ((harvest_scenariotargets$mod_target*harvest_scenariotargets$baseline)/10)*ynum
}
for (ynum in 11:31){
  harvest_imp[,ynum,2] <- harvest_scenariotargets$mod_target*harvest_scenariotargets$baseline
}

#  fill in annual implementation rates for scenario #3: ambitious implementation
for (ynum in 1:10){ 
  harvest_imp[,ynum,3] <- ((harvest_scenariotargets$amb_target*harvest_scenariotargets$baseline)/10)*ynum
}
for (ynum in 11:31){
  harvest_imp[,ynum,3] <- harvest_scenariotargets$amb_target*harvest_scenariotargets$baseline # maintain constant implementation rate for years 11 - 31
}

## RUN SCENARIOS FOR Timber Harvest 
# set up empty data frames to hold simulation results
harvest_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),dim(harvest_carbonstocks)[1]))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # ACTIVITY LOOP
  for (act in 1:dim(harvest_carbonstocks)[1]) {
    act_name = harvest_carbonstocks$activity[act]
    for (y in 1:length(years)){
      imp_rate = if (irate == 1) {
        imp_rate = harvest_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = harvest_imp[act,y,2]
      }else {
        imp_rate = harvest_imp[act,y,3]
      }
      
      d_rate = discount_rate_forest
      
      #########################################
      ##			 AVOIDED CONVERSION 	   ##
      #########################################
      act_stock_mean = harvest_carbonstocks[which(harvest_carbonstocks$activity %in% act_name), 3]
      act_stock_uc	 = harvest_carbonstocks[which(harvest_carbonstocks$activity %in% act_name), 4]
      act_seq_mean = harvest_carbonrates[which(harvest_carbonrates$activity %in% act_name), 3]
      act_seq_sd	 = harvest_carbonrates[which(harvest_carbonrates$activity %in% act_name), 4]
      
      imp_per_year = calc_area_mod(harvest_imp, d_rate, yrs)
      
      harvest_carbon_sum[irate,,y, act] = sapply(imp_per_year[y], function(x) {
        # ongoing sequestration due to lengthened rotation on evenage managed stands                                  
        (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
          # single event avoided emissions portion due to annual decrease in timber harvest
          (
            rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
          )
      })
      
    }
  }
  results_out[[irate]] = harvest_carbon_sum
  
  
  print(paste('finished harvest implementation scenario ', irate))
} 

## Summarize harvest Monte Carlo RESULTS to median and 95% CI for each year, each activity 
harvest_CI_High <- apply(harvest_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(harvest_CI_High) <- scenario.names
harvest_median <- apply(harvest_carbon_sum,c(1,3,4),median)/1e6
rownames(harvest_median) <- scenario.names
harvest_CI_Low <- apply(harvest_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(harvest_CI_Low) <- scenario.names

## reshape data
scen <- rep(c("Limited","Moderate", "Ambitious"), each = 31)
LimitedImp_harvest <- cbind(rowSums(cbind(harvest_median[1,1:31,1:14])),rowSums(cbind(harvest_CI_High[1,1:31,1:14])),rowSums(cbind(harvest_CI_Low[1,1:31,1:14])))
colnames(LimitedImp_harvest)=c("med","upper","lower")
Amb_harvest <- cbind(rowSums(cbind(harvest_median[3,1:31,1:14])),rowSums(cbind(harvest_CI_High[3,1:31,1:14])),rowSums(cbind(harvest_CI_Low[3,1:31,1:14])))
colnames(Amb_harvest)=c("med","upper","lower")
Mod_harvest <- cbind(rowSums(cbind(harvest_median[2,1:31,1:14])),rowSums(cbind(harvest_CI_High[2,1:31,1:14])),rowSums(cbind(harvest_CI_Low[2,1:31,1:14])))
colnames(Mod_harvest)=c("med","upper","lower")
harvest_by_year_all_ownerships <- cbind.data.frame(rbind(LimitedImp_harvest,Mod_harvest,Amb_harvest),scen)
harvest_by_year_all_ownerships$year <- rep(yr,3)
harvest_by_year_all_ownerships$actname <- rep("Timber Harvest",93)

# 2. Reforestation after wildfires (reforest) ####
# Reforestation INPUT FILES 
reforest_inputs <- read.csv("R_Inputs/reforestation_inputs.csv")
reforest_inputs$activity<-as.factor(reforest_inputs$activity)
reforest_inputs$activity <- droplevels(reforest_inputs$activity)
reforest_inputs$type<-as.factor(reforest_inputs$type)

# Set up scenario implementation targets

# Reforestation after wildfires scenario narrative description
# All scenarios include a linear decrease to the target rate by 2030
# Limited Implementation = Increased replanting rate by 100% over the first 10 years, steady after 2030
# Moderate Implementation = Gradually increased replanting rate by 100% in 2030 and  150% baseline by 2050
# Ambitious Implementation = Rapidly increased replanting rate by 100% in 2030, 300% in 2040, and 700% of baseline in 2050. In 2050, 63 to 84% of wildfire area replanted each year. 

reforest_target_2030 <- c(1.0, 1.0, 1.0)
reforest_target_2040 <- c(1.0, 1.25, 1.5)
reforest_target_2050 <- c(1.0, 3.0, 7.0)

# simulate fire area based on historical values for 30 year period
fire_area <- array(NA,dim=c(length(reforest_imprates$activity),1000,length(years)), dimnames=list(reforest_imprates$activity,iter,years))
for (i in 1:1000){
  fire_area[1,iter,1:31] <- calc_wildfire_area(reforest_imprates$fire_mean[1],reforest_imprates$fire_sd[1],years)
  fire_area[2,iter,1:31] <- calc_wildfire_area(reforest_imprates$fire_mean[2],reforest_imprates$fire_sd[2],years)
  fire_area[3,iter,1:31] <- calc_wildfire_area(reforest_imprates$fire_mean[3],reforest_imprates$fire_sd[3],years)
}

# build arrays to hold the annual implementation rates for scenario simulations

reforest_imp <- array(NA,dim=c(length(reforest_inputs$activity),length(years),imp_scenarios), dimnames=list(reforest_inputs$activity,years,scenario.names))

# fill in annual implementation rates for scenario #1: limited implementation
for (ynum in 1:9){  # ramp up period of 10 years to reach the historical variation (essentially a doubling by 2030)
  reforest_imp[,ynum,1] <- ((reforest_target_2030[1]*reforest_inputs$mean_reveg)/10)*ynum
}
reforest_imp[,10:31,1] <- (reforest_target_2030[1]*reforest_inputs$mean_reveg) # constant rate after 2030

# fill in annual implementation rates for scenario #2:moderate implementation scenario
for (ynum in 1:10){ # ramp up to 2030 target
  reforest_imp[,ynum,2] <- ((reforest_target_2030[2]*reforest_inputs$mean_reveg)/10)*ynum  
}
for (ynum in 11:20){  # ramp up to 2040 
  reforest_imp[,ynum,2] <- ((reforest_target_2040[2]*reforest_inputs$mean_reveg-reforest_target_2030[2]*reforest_inputs$mean_reveg)/10)*(ynum-10)+reforest_target_2030[2]*reforest_inputs$mean_reveg
}
for (ynum in 21:31){  # ramp up to 2050 target
  reforest_imp[,ynum,2] <- ((reforest_target_2050[2]*reforest_inputs$mean_reveg-reforest_target_2040[2]*reforest_inputs$mean_reveg)/10)*(ynum-20)+reforest_target_2040[2]*reforest_inputs$mean_reveg
}

# fill in annual implementation rates for scenario #3: ambitious implementation
for (ynum in 1:10){ # ramp up to 2030 target
  reforest_imp[,ynum,3] <- ((reforest_target_2030[3]*reforest_inputs$mean_reveg)/10)*ynum  
}
for (ynum in 11:20){  # ramp up to 2040 
  reforest_imp[,ynum,3] <- ((reforest_target_2040[3]*reforest_inputs$mean_reveg-reforest_target_2030[3]*reforest_inputs$mean_reveg)/10)*(ynum-10)+reforest_target_2030[3]*reforest_inputs$mean_reveg
}
for (ynum in 21:31){  # ramp up to 2050 target
  reforest_imp[,ynum,3] <- ((reforest_target_2050[3]*reforest_inputs$mean_reveg-reforest_target_2040[3]*reforest_inputs$mean_reveg)/10)*(ynum-20)+reforest_target_2040[3]*reforest_inputs$mean_reveg
}
reforest_imp[reforest_imp>1.0] <- 1.0

## implementation scenarios for reforestation 
# carbon data
model_carbon_reforest_low <- read.csv("R_Inputs/model_carbon_reforestation_low_v2.csv")
model_carbon_reforest_mod <- read.csv("R_Inputs/model_carbon_reforestation_mod_v2.csv")
model_carbon_reforest_hi <- read.csv("R_Inputs/model_carbon_reforestation_high_v2.csv")

# set up empty data frames or lists
reforest_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),length(reforest_inputs$activity)))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # ACTIVITY LOOP
  for (act in 1:length(reforest_inputs$activity)) {
    act_name = reforest_inputs$activity[act]
    # imp rate is calculated for each area in each scenario; this provide a matrix of iterations*years
    # such that for year 1, iteration 10, imp_rate == imp_rate[10,1]
    imp_rate = reforest_imp[act,1:31,irate]*fire_area[act,,]
    d_rate = discount_rate_forest
    
    # #########################################
    # ##			   SEQUESTRATION 		   ##
    # #########################################
    # get the correct C seq rates
    if (act_name == 'reforestation_low') {
      carbon_rates_reforestation = model_carbon_reforest_low
    }else if (act_name == 'reforestation_moderate') {
      carbon_rates_reforestation = model_carbon_reforest_mod
    }else {
      carbon_rates_reforestation = model_carbon_reforest_hi
    }
    
    # apply different seq rate and sd depending on age of cohort
    for (i in 1:iterations){
      # i=1
      i_imp = imp_rate[i,]
      # calculate cohort area
      cohort_area = calc_cohort_area_iter(i_imp,d_rate,length(years))
      cohort_carbon =	apply(cohort_area, 1, function(x) {
        y = x[!is.na(x)]
        act_mean = rev(carbon_rates_reforestation[1:length(y), 2])
        act_sd	 = rev(carbon_rates_reforestation[1:length(y), 3])
        #       
        out = sapply(1:length(y), function(p) {
          y[p] * rnorm(1, mean = act_mean[p], sd = act_sd[p])
        })
        return(out)
      })
      reforest_carbon_sum[irate,i, , act] = t(as.numeric(lapply(cohort_carbon,sum))*-1)
    }
  }
  results_out[[irate]] = reforest_carbon_sum
  
  print(paste('finished Reforest implementation scenario ', irate))
} 

## summarize and save Reforestation results 
## Summarize to median and 95% CI for each year, for low, moderate, and high productivity forests (see methods Graves et al. 2020)
reforest_CI_High <- apply(reforest_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(reforest_CI_High) <- scenario.names
reforest_median <- apply(reforest_carbon_sum,c(1,3,4),median)/1e6
rownames(reforest_median) <- scenario.names
reforest_CI_Low <- apply(reforest_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(reforest_CI_Low) <- scenario.names

## reshape data
# scen <- rep(c("Limited",   "Moderate" , "Ambitious"), each = 31)
LowImp_reforest<- cbind(rowSums(cbind(reforest_median[1,1:31,1:3])),rowSums(cbind(reforest_CI_High[1,1:31,1:3])),rowSums(cbind(reforest_CI_Low[1,1:31,1:3])))
colnames(LowImp_reforest)=c("med","upper","lower")
Amb_reforest <- cbind(rowSums(cbind(reforest_median[3,1:31,1:3])),rowSums(cbind(reforest_CI_High[3,1:31,1:3])),rowSums(cbind(reforest_CI_Low[3,1:31,1:3])))
colnames(Amb_reforest)=c("med","upper","lower")
Mod_reforest <- cbind(rowSums(cbind(reforest_median[2,1:31,1:3])),rowSums(cbind(reforest_CI_High[2,1:31,1:3])),rowSums(cbind(reforest_CI_Low[2,1:31,1:3])))
colnames(Mod_reforest)=c("med","upper","lower")
reforest_by_year_all_prodclass <- cbind.data.frame(rbind(LowImp_reforest,Mod_reforest,Amb_reforest),scen)
reforest_by_year_all_prodclass$year <- rep(yr,3)
reforest_by_year_all_prodclass$actname <- rep("Replant Wildfires",93)

# 3. Avoided conversion of forests (ACF) --------------------------------------
## load NCS input file
NCS_inputs = read.csv("R_Inputs/Oregon_NCS_all_other_model_inputs.csv") 
AC_forests_imprates = NCS_inputs[c(1:4),c(1:4)] # subest the forest conversion values
AC_forests_imprates$activity <-as.factor(AC_forests_imprates$activity)
AC_forests_imprates$activity <-droplevels(AC_forests_imprates$activity)
ACF_carbon_stocks = NCS_inputs[c(1:4),c(1:2,7:8)]
ACF_carbon_stocks$activity<-as.factor(ACF_carbon_stocks$activity)
ACF_carbon_stocks$activity <- droplevels(ACF_carbon_stocks$activity)
ACF_carbon_stocks$type<-as.factor(ACF_carbon_stocks$type)

ACF_carbon_rates = NCS_inputs[c(1:4),c(1:2,5:6)]
ACF_carbon_rates$activity<-as.factor(ACF_carbon_rates$activity)
ACF_carbon_rates$activity <- droplevels(ACF_carbon_rates$activity)
ACF_carbon_rates$type<-as.factor(ACF_carbon_rates$type)

# change in response to reviewer comments to reflect Woodbury et al. 2007 estimates of C loss in urban vs. rural forests (10/29/2019)
ACF_carbon_rates[1:2,3]<- ACF_carbon_rates[1:2,3]*0.84

## set implementation targets for the 3 scenarios
AC_forests_imprates # review baseline (annual ha converted) and historical variation (CV)

# Avoided Conversion of Forests narrative description
# All scenarios include a linear decrease to the target rate by 2030
# Limited Implementation = Reduce conversions by 10% by 2030, then stays constant
# Moderate Implementation = Reduce conversions by 50% by 2030, then stays constant
# Ambitious Implementation = Reduce conversions by 100% by 2030, then remain at zero

AC_forests_target2030 <- c(0.1,0.5,1.0) # percent reduction in conversion for each target decade in each of three scenarios
AC_forests_target2040 <- c(0.1,0.5,1.0)
AC_forests_target2050 <- c(0.1,0.5,1.0)

# build scenario implementation rates 

ACF_imp <- array(NA,dim=c(length(AC_forests_imprates$activity),length(years),imp_scenarios), dimnames=list(AC_forests_imprates$activity,years,scenario.names))
# limited implementation
for (ynum in 1:9){  # ramp up period of 10 years to reach the 2030 target
  ACF_imp[,ynum,1] <- ((AC_forests_imprates$baseline*AC_forests_target2030[1])/10)*ynum
}
ACF_imp[,10:31,1] <- AC_forests_imprates$baseline*AC_forests_target2050[1] # remain constant from 2030 - 2050

#moderate implementation
for (ynum in 1:9){  
  ACF_imp[,ynum,2] <- ((AC_forests_imprates$baseline*AC_forests_target2030[2])/10)*ynum
}
ACF_imp[,10:31,2] <- AC_forests_imprates$baseline*AC_forests_target2050[2] # remain constant from 2030 - 2050

# ambitious implementation
for (ynum in 1:9){  # ramp up period of 10 years to reach the 2030 target
  ACF_imp[,ynum,3] <- ((AC_forests_imprates$baseline*AC_forests_target2030[3])/10)*ynum
}
ACF_imp[,10:31,3] <- AC_forests_imprates$baseline*AC_forests_target2050[3] # remain constant from 2030 - 2050


## run scenarios of avoided forest conversion 
# set up empty data frames to hold simulation results
ACF_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),length(AC_forests_imprates$activity)))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # ACTIVITY LOOP
  for (act in 1:dim(ACF_carbon_stocks)[1]){
    act_name = ACF_carbon_stocks$activity[act]
    for (y in 1:length(years)){  ## get implementation rate for each year and activity
      imp_rate = if (irate == 1) {  ## IF NUMBER OF SCENARIOS CHANGES, DOUBLE CHECK THIS SECTION OF CODE!!
        imp_rate = ACF_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = ACF_imp[act,y,2]
      }else {
        imp_rate = ACF_imp[act,y,3]
      }
      
      d_rate = discount_rate_forest
      
      
      # #####################################
      ##			 AVOIDED CONVERSION 	   ##
      #########################################
      act_stock_mean = ACF_carbon_stocks[which(ACF_carbon_stocks$activity %in% act_name), 3]
      act_stock_uc	 = ACF_carbon_stocks[which(ACF_carbon_stocks$activity %in% act_name), 4]
      act_seq_mean = ACF_carbon_rates[act, 3]
      act_seq_sd	 = ACF_carbon_rates[act, 4]
      
      area_per_year = calc_area_mod(ACF_imp, d_rate, yrs)
      # y=1
      ACF_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
        # ongoing sequestration
        (rtruncnorm(iterations, a=0, b=Inf, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
          # single event AC portion
          (
            rtruncnorm(iterations, a=0,b=Inf,mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
          )
      })
    }
  }
  results_out[[irate]] = ACF_carbon_sum
  
  print(paste('finished implementation scenario ', irate))
}

## Summarize ACF RESULTS to median and 95% CI for each year, each activity 
ACF_CI_High <- apply(ACF_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(ACF_CI_High) <- scenario.names
ACF_median <- apply(ACF_carbon_sum,c(1,3,4),median)/1e6
rownames(ACF_median) <- scenario.names
ACF_CI_Low <- apply(ACF_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(ACF_CI_Low) <- scenario.names

## reshape data
LowImp_ACF <- cbind(rowSums(cbind(ACF_median[1,1:31,1:4])),rowSums(cbind(ACF_CI_High[1,1:31,1:4])),rowSums(cbind(ACF_CI_Low[1,1:31,1:4])))
colnames(LowImp_ACF)=c("med","upper","lower")
Amb_ACF <- cbind(rowSums(cbind(ACF_median[3,1:31,1:4])),rowSums(cbind(ACF_CI_High[3,1:31,1:4])),rowSums(cbind(ACF_CI_Low[3,1:31,1:4])))
colnames(Amb_ACF)=c("med","upper","lower")
Mod_ACF <- cbind(rowSums(cbind(ACF_median[2,1:31,1:4])),rowSums(cbind(ACF_CI_High[2,1:31,1:4])),rowSums(cbind(ACF_CI_Low[2,1:31,1:4])))
colnames(Mod_ACF)=c("med","upper","lower")
ACF_by_year_all <- cbind.data.frame(rbind(LowImp_ACF,Mod_ACF,Amb_ACF),scen)
ACF_by_year_all$year <- rep(yr,3)
ACF_by_year_all$actname <- rep("Forest- Avoided Conversion",93)


# 4. Reforestation of riparian areas (rref) -------------------------------------------------
rref_imprates = NCS_inputs[c(8:9),c(1:4)]
rref_imprates$activity <-as.factor(rref_imprates$activity)
rref_imprates$activity <-droplevels(rref_imprates$activity)
rref_carbon_stocks = NCS_inputs[c(8:9),c(1:2,7:8)]
rref_carbon_stocks$activity<-as.factor(rref_carbon_stocks$activity)
rref_carbon_stocks$activity <- droplevels(rref_carbon_stocks$activity)
rref_carbon_stocks$type<-as.factor(rref_carbon_stocks$type)

rref_carbon_rates = NCS_inputs[c(8:9),c(1:2,5:6)]
rref_carbon_rates$activity<-as.factor(rref_carbon_rates$activity)
rref_carbon_rates$activity <- droplevels(rref_carbon_rates$activity)
rref_carbon_rates$type<-as.factor(rref_carbon_rates$type)

# update 10/30/2019 to include maximum cumulative riparian reforestation area (ha) estimates
# max estimate methods, see Graves et al. 
rref_max_east = 76634.55
rref_max_west = 125780.4

rref_max=rbind(rref_max_east,rref_max_west)

# set implementation targets for the 3 scenarios
rref_imprates # review baseline (annual ha restored in east and western OR) and historical variation (CV)

# Riparian Reforestation narrative description
# All scenarios include a linear decrease to the target rate by 2030
# Limited Implementation = Increased riparian reforestation by 72% in eastern OR and 42% in western OR (based on historical variation); allows 10 years to reach target rates
# Moderate Implementation = Doubling (100% increase) in riparian reforestation rate by 2050, 250% increase by 2040, and 300% increase by 2050. 
# Ambitious Implementation = Rapidly increased riparian reforestation rates to reach the maximum area available by 2050. Maximum area estimated as 76,635 ha east of the Cascades and 125,780 ha west of the Cascades. East side riparian reforestation increases to 500% by 2030 and reaches the maximum area threshold in 2032. West side riparian reforestation increases by 1000% (10x) by 2030 and reaches maximum area threshold by 2044.

rref_target2030 <- cbind(c(0.72,1.0,5.0),c(0.42, 1.0,10.0)) # percent increase in restoration for each target decade in each of three scenarios
rref_target2040 <- cbind(c(0.72,2.5,5.0),c(0.42, 2.5,10.0))
rref_target2050 <- cbind(c(0.72,3.0,5.0),c(0.42, 3.0,10.0))

# build array to hold yearly implementation rates for scenarios
rref_imp <- array(NA,dim=c(length(rref_imprates$activity),length(years),imp_scenarios), dimnames=list(rref_imprates$activity,years,scenario.names))

# limited implementation: assume increase equivalent to historical variation to be reached by 2030
for (ynum in 1:9){  # ramp up period of 10 years to reach the historical variation
  rref_imp[1,ynum,1] <- ((rref_imprates$baseline[1]*rref_target2030[1,1])/10)*ynum
  rref_imp[2,ynum,1] <- ((rref_imprates$baseline[2]*rref_target2030[1,2])/10)*ynum
}
rref_imp[1,10:31,1] <- rref_imprates$baseline[1]*rref_target2030[1,1] # east riparian restoration
rref_imp[2,10:31,1] <- rref_imprates$baseline[2]*rref_target2030[1,2] # western riparian restoration

# moderate implementation
for (ynum in 1:10){ # first 10 years
  rref_imp[,ynum,2] <- ((rref_imprates$baseline*rref_target2030[2])/10)*ynum
}
for (ynum in 11:20){ # 2030 to 2040; starts at 2020 target and increases to reach 2040 target
  rref_imp[,ynum,2] <- rref_target2030[2]*rref_imprates$baseline + ((rref_imprates$baseline*rref_target2040[2] - rref_imprates$baseline*rref_target2030[2])/10)*(ynum-10)
}
for (ynum in 21:31){ # starts at 2040 target and increases to reach 2050 target
  rref_imp[,ynum,2] <- rref_target2040[2]*rref_imprates$baseline + ((rref_imprates$baseline*rref_target2050[2] - rref_imprates$baseline*rref_target2040[2])/11)*(ynum-20)
}

# ambitious implementation #
# need to include baseline implementation and maximum threshold
rref_MPtarget_east <- 57394.23 # assumes that baseline rate of restoration continues; this is cumulative max additional restoration potential
rref_MPtarget_west <- 118110.6 # assumes that baseline rate of restoration continues; this is cumulative max additional restoration potential
rref_maxtarget <- rbind(rref_MPtarget_east,rref_MPtarget_west)

#east side
for (ynum in 1:10){
  rref_imp[1,ynum,3] <- ((rref_imprates$baseline[1]*rref_target2030[3,1])/10)*ynum
}
for (ynum in 11:31){
  rref_imp[1,ynum,3] <- (rref_imprates$baseline[1]*rref_target2030[3,1])
}

#west side
for (ynum in 1:10){
  rref_imp[2,ynum,3] <- ((rref_imprates$baseline[2]*rref_target2030[3,2])/10)*ynum
}
for (ynum in 11:31){
  rref_imp[2,ynum,3] <- ((rref_imprates$baseline[2]*rref_target2030[3,2]))
}

## run riparian reforest scenarios
# set up empty data frames or lists
rref_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),2))

results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=3
  # ACTIVITY LOOP
  for (act in 1:2) {
    # act = 1
    act_name = rref_carbon_rates$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = rref_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = rref_imp[act,y,2]
      }else {
        imp_rate = rref_imp[act,y,3]
      }
      
      d_rate = discount_rate_reforest
      
      # #########################################
      # ##			   SEQUESTRATION 		   ##
      # #########################################
      seq_mean = rref_carbon_rates[act, 3]
      seq_sd	 = rref_carbon_rates[act, 4]
      #     
      area_per_year = calc_area_mod2(rref_imp, rref_maxtarget,d_rate, length(years))
      # 
      rref_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
        rnorm(iterations, mean = seq_mean, sd = seq_sd) * x * -1
      })
    }
  }
  results_out[[irate]] = rref_carbon_sum
  
  print(paste('finished Reforest implementation scenario ', irate))
} 

## summarize reforestation results 
rref_CI_High <- apply(rref_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(rref_CI_High) <- scenario.names
rref_median <- apply(rref_carbon_sum,c(1,3,4),median)/1e6
rownames(rref_median) <- scenario.names
rref_CI_Low <- apply(rref_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(rref_CI_Low) <- scenario.names

## reshape data
LowImp_rref <- cbind(rowSums(cbind(rref_median[1,1:31,1:2])),rowSums(cbind(rref_CI_High[1,1:31,1:2])),rowSums(cbind(rref_CI_Low[1,1:31,1:2])))
colnames(LowImp_rref)=c("med","upper","lower")
Amb_rref <- cbind(rowSums(cbind(rref_median[3,1:31,1:2])),rowSums(cbind(rref_CI_High[3,1:31,1:2])),rowSums(cbind(rref_CI_Low[3,1:31,1:2])))
colnames(Amb_rref)=c("med","upper","lower")
Mod_rref <- cbind(rowSums(cbind(rref_median[2,1:31,1:2])),rowSums(cbind(rref_CI_High[2,1:31,1:2])),rowSums(cbind(rref_CI_Low[2,1:31,1:2])))
colnames(Mod_rref)=c("med","upper","lower")
rref_by_year_all <- cbind.data.frame(rbind(LowImp_rref,Mod_rref,Amb_rref),scen)
rref_by_year_all$year <- rep(yr,3)
rref_by_year_all$actname <- rep("Riparian Reforestation",93)

# 5. Sage-steppe avoided conversion and restoration -----------------------------------------
sage_inputs <- NCS_inputs[c(5:6),]
sage_inputs$activity <-as.factor(sage_inputs$activity)
sage_inputs$activity <-droplevels(sage_inputs$activity)
sage_carbon_stocks = sage_inputs[,c(1:2,7:8)]
sage_carbon_stocks$type<-as.factor(sage_carbon_stocks$type)
sage_carbon_rates = sage_inputs[,c(1:2,5:6)]
## set scenario implementation rates 

sage_inputs <- sage_inputs %>%  # set maximum area estimates
  mutate(max=c(3869790,905970))

sage_imp <- array(NA,dim=c(length(sage_inputs$activity),length(years),imp_scenarios), dimnames=list(sage_inputs$activity,years,scenario.names))

# limited implementation: assume historical variation of 10%, to be reached by 2030
sage_imp[,10:31,1] <- sage_inputs$baseline*0.10 # historical variation
for (ynum in 1:9){  # ramp up period of 10 years to reach the historical variation
  sage_imp[,ynum,1] <- ((sage_inputs$baseline*0.10)/10)*ynum
}
# ambitious implementatio
for (ynum in 1:9) {
  sage_imp[1,ynum,3] <- sage_inputs$baseline[1]/10*ynum
  sage_imp[2,ynum,3] <- (sage_inputs$baseline[2]*3-sage_inputs$baseline[2])/10*ynum
}
sage_imp[1,10:31,3] <- sage_inputs$baseline[1]
sage_imp[2,10:31,3] <- sage_inputs$baseline[2]*3-sage_inputs$baseline[2]

# moderate implementation
for (ynum in 1:9) {
  sage_imp[1,ynum,2] <- sage_inputs$baseline[1]*0.1
  sage_imp[2,ynum,2] <- (sage_inputs$baseline[2])/10*ynum
}
for (ynum in 10:31) {
  sage_imp[1,ynum,2] <- sage_inputs$baseline[1]*0.2
  sage_imp[2,ynum,2] <- (sage_inputs$baseline[2])/20*(ynum-10)+sage_inputs$baseline[2]
}

# run implementation scenarios for sage 
sage_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),2))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # ACTIVITY LOOP
  for (act in 1:2) {
    act_name =sage_carbon_rates$activity[act]
    for (y in 1:length(years)){
      imp_rate = if (irate == 1) {
        imp_rate = sage_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = sage_imp[act,y,2]
      }else {
        imp_rate = sage_imp[act,y,3]
      }
      
      d_rate = discount_rate_other
      
      # #########################################
      # ##			   SEQUESTRATION 		   ##
      # #########################################
      if (sage_carbon_rates$type[act] == 'seq') {
        act_mean = sage_carbon_rates[act, 3]
        act_sd	 = sage_carbon_rates[act, 4]
        #     
        area_per_year = calc_area_mod(sage_imp, d_rate, length(years))
        #     
        sage_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          rnorm(iterations, mean = act_mean, sd = act_sd) * x * -1
        })
      }
      
      # #####################################
      ##			 AVOIDED CONVERSION 	   ##
      #########################################
      if (sage_carbon_rates$type[act] == 'AC') {
        act_stock_mean = sage_carbon_stocks[which(sage_carbon_stocks$activity %in% act_name), 3]
        act_stock_uc	 = sage_carbon_stocks[which(sage_carbon_stocks$activity %in% act_name), 4]
        act_seq_mean = sage_carbon_rates[act, 3]
        act_seq_sd	 = sage_carbon_rates[act, 4]
        
        area_per_year = calc_area_mod(sage_imp, d_rate, yrs)
        # y=1
        sage_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          # ongoing sequestration
          (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
            # single event AC portion
            (
              rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
            )
        })
      }
      ##### End Sage Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = sage_carbon_sum
  
  print(paste('finished sage implementation scenario ', irate))
} 

## summarize sage results
sage_CI_High <- apply(sage_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(sage_CI_High) <- scenario.names
sage_median <- apply(sage_carbon_sum,c(1,3,4),median)/1e6
rownames(sage_median) <- scenario.names
sage_CI_Low <- apply(sage_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(sage_CI_Low) <- scenario.names

## reshape data

LowImp_sage <- cbind(rowSums(cbind(sage_median[1,1:31,1:2])),rowSums(cbind(sage_CI_High[1,1:31,1:2])),rowSums(cbind(sage_CI_Low[1,1:31,1:2])))
colnames(LowImp_sage)=c("med","upper","lower")
Amb_sage <- cbind(rowSums(cbind(sage_median[3,1:31,1:2])),rowSums(cbind(sage_CI_High[3,1:31,1:2])),rowSums(cbind(sage_CI_Low[3,1:31,1:2])))
colnames(Amb_sage)=c("med","upper","lower")
Mod_sage <- cbind(rowSums(cbind(sage_median[2,1:31,1:2])),rowSums(cbind(sage_CI_High[2,1:31,1:2])),rowSums(cbind(sage_CI_Low[2,1:31,1:2])))
colnames(Mod_sage)=c("med","upper","lower")
sage_by_year_all <- cbind.data.frame(rbind(LowImp_sage,Mod_sage,Amb_sage),scen)
sage_by_year_all$year <- rep(yr,3)
sage_by_year_all$actname <- rep("Sagebrush-steppe pathways",93)


# 6. Restoration of tidal wetlands (tide) --------------------------------------
tide_inputs <- NCS_inputs[c(7),]
tide_inputs$activity <-as.factor(tide_inputs$activity)
tide_inputs$activity <-droplevels(tide_inputs$activity)

## set scenario implementation rates 
tide_inputs <- tide_inputs %>%
  mutate(max=5205)  # set max area

tide_imp <- array(NA,dim=c(length(tide_inputs$activity),length(years),imp_scenarios), dimnames=list(tide_inputs$activity,years,scenario.names))
# limited implementation: assume historical variation 
tide_imp[,10:31,1] <- 50 # increase by 50 ha in 2030 (>100% h.v.)
for (ynum in 1:9){  # ramp up period of 10 years to reach the historical variation
  tide_imp[,ynum,1] <- 5*ynum
}
# ambitious implementation
for (ynum in 1:5) {
  tide_imp[,ynum,3] <- (tide_inputs$baseline)/5*ynum # doubling = increase by baseline in 5 years
}
for (ynum in 6:15) {
  tide_imp[,ynum,3] <- (tide_inputs$baseline*4-tide_inputs$baseline*2)/10*(ynum-5)+tide_inputs$baseline # doubling again; target = 4*baseline-2*baseline
}

for (ynum in 16:26) {
  tide_imp[,ynum,3] <- (tide_inputs$baseline*8-tide_inputs$baseline*4)/15*(ynum-15)+tide_inputs$baseline*3
}

tide_imp[,27,3] <- 50  ### REACH MAXIMUM AREA OF 5200 HA OF HIGHLY SALINE TIDAL WETLANDS AT THIS POINT

for (ynum in 28:31) {
  tide_imp[,ynum,3] <- 0
}

# moderate implementation; take 10 years to incrase by 100 ha/year
for (ynum in 1:10) {
  tide_imp[,ynum,2] <- 100/10*ynum
}

tide_imp[,11:31,2] <- 100


# run implementation scenarios for tidal wetlands
tide_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),1))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # irate=1
  # ACTIVITY LOOP
  for (act in 1:length(tide_inputs$activity)) {
    # act = 1
    act_name =tide_inputs$activity[act]
    for (y in 1:length(years)){
      # y=1
      imp_rate = if (irate == 1) {
        imp_rate = tide_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = tide_imp[act,y,2]
      }else {
        imp_rate = tide_imp[act,y,3]
      }
      
      d_rate = discount_rate_other
      
      # #####################################
      ##			 AVOIDED CONVERSION 	   ##
      #########################################
      if (tide_inputs$type[act] == 'AC') {
        act_stock_mean = tide_inputs[which(tide_inputs$activity %in% act_name), 7]
        act_stock_uc	 = tide_inputs[which(tide_inputs$activity %in% act_name), 8]
        act_seq_mean = tide_inputs[act, 5]
        act_seq_sd	 = tide_inputs[act, 6]
        
        area_per_year = calc_area_mod(tide_imp, d_rate, yrs)
        # y=1
        tide_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          # ongoing sequestration
          (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
            # for tidal wetlands AC portion is a per ha per year rate
            (
              rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * x * -1
            )
        })
      }
      ##### End tide Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = tide_carbon_sum
  
  print(paste('finished tide implementation scenario ', irate))
} 

## summarize tide results 
tide_CI_High <- apply(tide_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(tide_CI_High) <- scenario.names
tide_median <- apply(tide_carbon_sum,c(1,3,4),median)/1e6
rownames(tide_median) <- scenario.names
tide_CI_Low <- apply(tide_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(tide_CI_Low) <- scenario.names

## reshape data
LowImp_tide <- cbind(rowSums(cbind(tide_median[1,1:31,])),rowSums(cbind(tide_CI_High[1,1:31,])),rowSums(cbind(tide_CI_Low[1,1:31,])))
colnames(LowImp_tide)=c("med","upper","lower")
Amb_tide <- cbind(rowSums(cbind(tide_median[3,1:31,])),rowSums(cbind(tide_CI_High[3,1:31,])),rowSums(cbind(tide_CI_Low[3,1:31,])))
colnames(Amb_tide)=c("med","upper","lower")
Mod_tide <- cbind(rowSums(cbind(tide_median[2,1:31,])),rowSums(cbind(tide_CI_High[2,1:31,])),rowSums(cbind(tide_CI_Low[2,1:31,])))
colnames(Mod_tide)=c("med","upper","lower")
tide_by_year_all <- cbind.data.frame(rbind(LowImp_tide,Mod_tide,Amb_tide),scen)
tide_by_year_all$year <- rep(yr,3)
tide_by_year_all$actname <- rep("Tidal Wetlands",93)


# 7. Grassland Avoided Conversion ----------------------------------
grass_inputs <- NCS_inputs[c(10),]
grass_inputs$activity <-as.factor(grass_inputs$activity)
grass_inputs$activity <-droplevels(grass_inputs$activity)

## set scenario implementation rates #

grass_imp <- array(NA,dim=c(length(grass_inputs$activity),length(years),imp_scenarios), dimnames=list(grass_inputs$activity,years,scenario.names))
# low implementation: assume historical variation 
grass_imp[,10:31,1] <- grass_inputs$baseline*0.10 # 10% of baseline
for (ynum in 1:9){  # ramp up period of 10 years to reach the historical variation
  grass_imp[,ynum,1] <- (grass_inputs$baseline*0.10)/10*ynum
}
# ambitious implementation: maximum potential 
for (ynum in 1:10) {
  grass_imp[,ynum,3] <- (grass_inputs$baseline)/10*ynum # decrease by equivalent of baseline amount by 2030
}
grass_imp[,11:31,3] <- grass_inputs$baseline # 100% of baseline (i.e. no conversion after 2030)

# moderate implementation - decreases conversion by 50% by 2030 and then keeps that
for (ynum in 1:10) {
  grass_imp[,ynum,2] <- ((grass_inputs$baseline)*0.5)/10*ynum 
}
for (ynum in 11:31) {
  grass_imp[,ynum,2] <- ((grass_inputs$baseline)*0.5)/21*(ynum-10) + grass_inputs$baseline*0.5
}


# run implementation scenarios for grass 
grass_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),1))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # ACTIVITY LOOP
  for (act in 1:length(grass_inputs$activity)) {
    act_name =grass_inputs$activity[act]
    for (y in 1:length(years)){
      imp_rate = if (irate == 1) {
        imp_rate = grass_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = grass_imp[act,y,2]
      }else {
        imp_rate = grass_imp[act,y,3]
      }
      
      d_rate = discount_rate_other
      
      # #####################################
      ##			 AVOIDED CONVERSION 	   ##
      #########################################
      if (grass_inputs$type[act] == 'AC') {
        act_stock_mean = grass_inputs[which(grass_inputs$activity %in% act_name), 7]
        act_stock_uc	 = grass_inputs[which(grass_inputs$activity %in% act_name), 8]
        act_seq_mean = grass_inputs[act, 5]
        act_seq_sd	 = grass_inputs[act, 6]
        
        area_per_year = calc_area_mod(grass_imp, d_rate, yrs)
        # y=1
        grass_carbon_sum[irate, ,y, act] = sapply(area_per_year[y], function(x) {
          # ongoing sequestration
          (rnorm(iterations, mean = act_seq_mean, sd = act_seq_sd) * x * -1) +
            # single event AC portion
            (
              rnorm(iterations, mean = act_stock_mean, sd = act_stock_uc) * imp_rate * d_rate * -1
            )
        })
      }
      ##### End grass Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = grass_carbon_sum
  
  print(paste('finished grass implementation scenario ', irate))
} 

## summarize grass results 
grass_CI_High <- apply(grass_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(grass_CI_High) <- scenario.names
grass_median <- apply(grass_carbon_sum,c(1,3,4),median)/1e6
rownames(grass_median) <- scenario.names
grass_CI_Low <- apply(grass_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(grass_CI_Low) <- scenario.names

## reshape data
LowImp_grass <- cbind(rowSums(cbind(grass_median[1,1:31,])),rowSums(cbind(grass_CI_High[1,1:31,])),rowSums(cbind(grass_CI_Low[1,1:31,])))
colnames(LowImp_grass)=c("med","upper","lower")
Amb_grass <- cbind(rowSums(cbind(grass_median[3,1:31,])),rowSums(cbind(grass_CI_High[3,1:31,])),rowSums(cbind(grass_CI_Low[3,1:31,])))
colnames(Amb_grass)=c("med","upper","lower")
Mod_grass <- cbind(rowSums(cbind(grass_median[2,1:31,])),rowSums(cbind(grass_CI_High[2,1:31,])),rowSums(cbind(grass_CI_Low[2,1:31,])))
colnames(Mod_grass)=c("med","upper","lower")
grass_by_year_all <- cbind.data.frame(rbind(LowImp_grass,Mod_grass,Amb_grass),scen)
grass_by_year_all$year <- rep(yr,3)
grass_by_year_all$actname <- rep("Grassland - AC",93)


# 8. Agricultural Pathways (ag) -------------------------------
ag_inputs <- NCS_inputs[c(11:13),]
ag_inputs$activity <-as.factor(ag_inputs$activity)
ag_inputs$activity <-droplevels(ag_inputs$activity)

ag_inputs$max <- c("1877500","916000",NA)

## set scenario implementation rates for agriculture 
ag_imp <- array(NA,dim=c(length(ag_inputs$activity),length(years),imp_scenarios), dimnames=list(ag_inputs$activity,years,scenario.names))

# assume historical variation 
for (ynum in 1:9){  # ramp up period of 10 years to reach the historical variation
  ag_imp[1,ynum,1] <- (ag_inputs$baseline[1]*ag_inputs$cv[1])/10*ynum
  ag_imp[2,ynum,1] <- (ag_inputs$baseline[2]*ag_inputs$cv[2])/10*ynum
  ag_imp[3,ynum,1] <- (ag_inputs$baseline[3]*(ag_inputs$cv[3]*0.4))/10*ynum # for N reduction, we assume 18% reduction on 40% of acres (18% is CV, 40% based on N.U.E.)
}
for (ynum in 10:31){ # historical variation
  ag_imp[1,ynum,1] <- (ag_inputs$baseline[1]*ag_inputs$cv[1])
  ag_imp[2,ynum,1] <- (ag_inputs$baseline[2]*ag_inputs$cv[2])
  ag_imp[3,ynum,1] <- (ag_inputs$baseline[3])*(ag_inputs$cv[3]*0.4)
}

# ambitious implementation potential
for (ynum in 1:31 ){ # linear increases in cover crop and no-till to reach max by 2050 (max is 50% of cropland and all tilled cropland)
  ag_imp[1,ynum,3] <- ((as.numeric(ag_inputs$max[1])*0.5)-ag_inputs$baseline[1])/30*ynum
  ag_imp[2,ynum,3] <- (as.numeric(ag_inputs$max[2])-ag_inputs$baseline[2])/30*ynum
}
for (ynum in 1:10) { # linear decrease to reach target of 40% on 40% of cropland by 2030
  ag_imp[3,ynum,3] <- ((ag_inputs$baseline[3]*(0.40*0.4))/10)*ynum
}
ag_imp[3,11:31,3] <- (ag_inputs$baseline[3]*(0.40*0.4))

# moderate implementation
for (ynum in 1:10){ # cover crops increase by 150% and no-till increases by 50%
  ag_imp[1,ynum,2] <- (ag_inputs$baseline[1]*1.5)/10*ynum
  ag_imp[2,ynum,2] <- (ag_inputs$baseline[2]*0.5)/10*ynum
}
for (ynum in 11:31) { # cover crops quadruple by 2050, no-till increase by 100%
  ag_imp[1,ynum,2] <- ((ag_inputs$baseline[1]*4-ag_inputs$baseline[1]*1.5)/20*(ynum-10))+(ag_inputs$baseline[1]*1.5)
  ag_imp[2,ynum,2] <- ag_inputs$baseline[2]
}

for (ynum in 1:10) { # linear decrease to reach target of 25% on 40% of acres (0.25*0.40) by 2030, decrease by 0.25 on 100% of acres by 2050
  ag_imp[3,ynum,2] <- ((ag_inputs$baseline[3]*(0.25*0.4))/10)*ynum
}
for (ynum in 11:31) {
  ag_imp[3,ynum,2] <- (((ag_inputs$baseline[3]*0.25*0.4-ag_inputs$baseline[3]*0.25*0.4)/21) *(ynum-10))+(ag_inputs$baseline[3]*0.25*0.40)
}

## run the ag implementation scenarios 
ag_carbon_sum = array(NA, dim=c(imp_scenarios,iterations,length(years),length(ag_inputs$activity)))
results_out = list()

# SCENARIO LOOP
for (irate in 1:imp_scenarios) {
  # ACTIVITY LOOP
  for (act in 1:length(ag_inputs$activity)) {
    act_name =ag_inputs$activity[act]
    for (y in 1:length(years)){
      imp_rate = if (irate == 1) {
        imp_rate = ag_imp[act,y,1]
      }else if (irate == 2) {
        imp_rate = ag_imp[act,y,2]
      }else {
        imp_rate = ag_imp[act,y,3]
      }
      
      d_rate = discount_rate_other
      
      # #########################################
      # ##			   SEQUESTRATION 		   ##
      # #########################################
      ###### Cover Crops ######
      if (act_name == 'cover_crop_use') {
        act_mean = ag_inputs[act, 5]
        act_sd	 = ag_inputs[act, 6]
        #     
        # area_per_year = calc_area_mod(ag_imp, d_rate, length(years))
        #
        ag_carbon_sum[irate, ,y, act] = rnorm(iterations, mean = act_mean, sd = act_sd) *imp_rate*d_rate* -1
      }
      
      ###### No-till Agriculture ####
      if (act_name == 'no_till_ag') {
        act_min = ag_inputs[act,5]
        act_max = ag_inputs[act,6]
        # area_per_year = calc_area_mod(ag_imp, d_rate, length(years))
        #
        ag_carbon_sum[irate, ,y,act] = runif(iterations,min=act_min, max=act_max)*imp_rate*d_rate* -1 # ongoing sequestration
      }
      
      # #####################################
      ##			 AVOIDED CONVERSION 	   ##
      #########################################
      if (act_name == 'N_management') {
        act_min = ag_inputs[act,7]
        act_max = ag_inputs[act,8]
        area_per_year = calc_area_mod(ag_imp, d_rate, length(years))
        #
        ag_carbon_sum[irate, ,y,act] = runif(iterations,min=act_min, max=act_max)*imp_rate * d_rate * -1 # avoided emissions
      }
      ##### End ag Simulations #####
    }  # close year loop
  }  # close activity loop
  results_out[[irate]] = ag_carbon_sum
  
  print(paste('finished ag implementation scenario ', irate))
} 

## Summarize agriculture results 
ag_CI_High <- apply(ag_carbon_sum,c(1,3,4),quantile,0.95)/1e6
rownames(ag_CI_High) <- scenario.names
ag_median <- apply(ag_carbon_sum,c(1,3,4),median)/1e6
rownames(ag_median) <- scenario.names
ag_CI_Low <- apply(ag_carbon_sum,c(1,3,4),quantile,0.05)/1e6
rownames(ag_CI_Low) <- scenario.names

# Ag results -- keep activity identity
CC_lowimp <- cbind(ag_median[1,1:31,1],ag_CI_High[1,1:31,1],ag_CI_Low[1,1:31,1])
CC_amb <- cbind(ag_median[3,1:31,1],ag_CI_High[3,1:31,1],ag_CI_Low[3,1:31,1])
CC_mod <- cbind(ag_median[2,1:31,1],ag_CI_High[2,1:31,1],ag_CI_Low[2,1:31,1])
CC_by_year <- cbind.data.frame(rbind(CC_lowimp,CC_mod,CC_amb),scen)
CC_by_year$year <- rep(yr,3)
CC_by_year$actname <- rep("Cover Crops",93)

NT_lowimp <- cbind(ag_median[1,1:31,2],ag_CI_High[1,1:31,2],ag_CI_Low[1,1:31,2])
NT_amb <- cbind(ag_median[3,1:31,2],ag_CI_High[3,1:31,2],ag_CI_Low[3,1:31,2])
NT_mod <- cbind(ag_median[2,1:31,2],ag_CI_High[2,1:31,2],ag_CI_Low[2,1:31,2])
NT_by_year <- cbind.data.frame(rbind(NT_lowimp,NT_mod,NT_amb),scen)
NT_by_year$year <- rep(yr,3)
NT_by_year$actname <- rep("No-Till",93)

Nmgmt_lowimp <- cbind(ag_median[1,1:31,3],ag_CI_High[1,1:31,3],ag_CI_Low[1,1:31,3])
Nmgmt_amb <- cbind(ag_median[3,1:31,3],ag_CI_High[3,1:31,3],ag_CI_Low[3,1:31,3])
Nmgmt_mod <- cbind(ag_median[2,1:31,3],ag_CI_High[2,1:31,3],ag_CI_Low[2,1:31,3])
Nmgmt_by_year <- cbind.data.frame(rbind(Nmgmt_lowimp,Nmgmt_mod,Nmgmt_amb),scen)
Nmgmt_by_year$year <- rep(yr,3)
Nmgmt_by_year$actname <- rep("N Mgmt",93)
ag_by_year_all_act <- rbind(CC_by_year, NT_by_year,Nmgmt_by_year)
colnames(ag_by_year_all_act)[1:3]=c("med","upper","lower")

# -----------------------------------------------
# Summarize all NCS Pathway Activities ----------------------------

 NCS_all <- rbind(IFM_by_year_all_ownerships,reforest_by_year_all_prodclass,rref_by_year_all,ACF_by_year_all,
                  tide_by_year_all, sage_by_year_all, grass_by_year_all, ag_by_year_all_act)
 unique(NCS_all$actname) 

 NCS_all$group <- ifelse(NCS_all$actname %in% c("Forest- Avoided Conversion","Grassland - AC",
                                                "Sagebrush-steppe pathways"),"Avoided Conversion",
                         ifelse(NCS_all$actname %in% c("Timber Harvest","Cover Crops","No-Till","N Mgmt"),
                                "Land Management", "Restoration"))  # rename pathways



 NCS_totals <- NCS_all %>%
   group_by(scen,year)%>%
   summarise(median=sum(med),
             upper=sum(upper),
             lower=sum(lower))

# calculate cumulative reductions
NCS_totals <- NCS_totals %>%
  mutate(median_cumsum = cumsum(median), upper_cumsum=cumsum(upper), lower_cumsum=cumsum(lower))

  
# calculate cumulative reduction by group
NCS_totals <- NCS_all %>%
  group_by(scen,year)%>%
  summarise(median=sum(med),
            upper=sum(upper),
            lower=sum(lower))

NCS_all_groups <- NCS_all %>%
  group_by(scen,year,group) %>% 
  summarise(median = sum(med), up = sum(upper), low=sum(lower))
NCS_all_groups$group <- as.factor(NCS_all_groups$group)

NCS_all_groups <- NCS_all_groups %>%
  group_by(scen,group) %>%
  mutate(median_cumsum=cumsum(median), upper_cumsum=cumsum(up), lower_cumsum=cumsum(low))

NCS_all_groups$scen <- factor(NCS_all_groups$scen,levels(NCS_all_groups$scen)[c(2,3,1)]) #reorder levels to be limited, moderate, ambitious
NCS_all_groups$group <- factor(NCS_all_groups$group,levels(NCS_all_groups$group)[c(2,3,1)]) #reorder levels to be limited, moderate, ambitious

# calculate cumulative reduction by activity
NCS_total_activity <- NCS_all %>%
  group_by(scen,actname) %>%
  mutate(median_cumsum=cumsum(med),upper_cumsum=cumsum(upper),lower_cumsum=cumsum(lower))

NCS_total_activity_2050 <- NCS_total_activity %>%
  filter(year=="2050") %>%
  group_by(scen) %>%
  mutate(total=sum(median_cumsum), prop=median_cumsum/total)

