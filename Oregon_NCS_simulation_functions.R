######################################
##									##
##									##
##		  Simulation Model		    ##
##		     Functions		        ##
##									##
##		 David Marvin (2017)
##     modified by Rose Graves (2019)
##     IN REVIEW
##     DO NOT DISTRIBUTE
######################################
  
  # calculate the total area of implementation by year, including the mortality rate	
	calc_area = function(imp_rate,discount,years){
		acres_in = imp_rate*discount
		acreage_out = acres_in
		for (i in 2:years){
			acreage_out[i] = (acreage_out[i-1] * discount) + acres_in
		}
		return(acreage_out)
	}
	
# calculate the total area of implementation by year, including mortality rate. Modified to allow for different implementation each year
	calc_area_mod = function(imp_rates,discount,years){
	  acres_in = imp_rates[act,1:years,irate]*discount
	  acreage_out = acres_in
	  for (i in 2:years){
	    acreage_out[i]=(acreage_out[i-1]*d_rate) + acres_in[i]
	  }
	  return(acreage_out)
	}
	
	# FOR RESTORATION PATHWAYS WITH A MAXIMUM: calculate the total area of implementation by year, including mortality rate. Modified to allow for different implementation each year and maximum area
	calc_area_mod2 = function(imp_rates,maximum, discount,years){
	  acres_in = imp_rates[act,1:years,irate]*discount
	  acreage_out = acres_in
	  for (i in 2:years){
	    acreage_out[i]=ifelse(((acreage_out[i-1]*d_rate) + acres_in[i]) > maximum[act], maximum[act], (acreage_out[i-1]*d_rate) + acres_in[i])
	  }
	  indx = which.max(acreage_out) # find which year maximum restoration is reached
	  acreage_out[indx:length(acreage_out)] = rep(acreage_out[indx]) # set all following years to the max value
	  return(acreage_out)
	}
	
# calculate the total area of implementation by year, including the mortality rate, and store cohort
	calc_cohort_area = function(imp_rate,discount,years){
		acreage_out = matrix(NA,years,years)
		cohort = 1:years
		for (c in cohort){
			yrs = seq(c,length(cohort),1) 
			acreage_out[yrs,c] = imp_rate * discount^(yrs-(c-1))
			}
		return(acreage_out)
	}

# calculate the total area of implementation by year, including mortality rate, and store cohort; modified to allow for changing imp-rates over years
	
	calc_cohort_area_iter = function(imp_rate,discount,years){
	  acreage_out = matrix(NA,years,years)
	  cohort = 1:years
	  for (c in cohort){
	    yrs = seq(c,length(cohort),1) 
	    acreage_out[yrs,c] = imp_rate[c] * discount^(yrs-(c-1))
	  }
	  return(acreage_out)
	}

	# calculate wildfire area
	calc_wildfire_area = function(fire_mean,fire_sd,years){
	  acreage_out =matrix(NA,nrow=length(years),ncol=1)
	  for(i in 1:length(years)){
	  acreage_out[i]=(rtruncnorm(1, mean = fire_mean, sd = fire_sd, a=0))
	  }
	  return(acreage_out)
	}
	
	#
	
	
	