library(stringr)

samples_to_parameters <- function(samples,parameters){
  for(j in names(samples)){
    if(str_detect(j,"vax\\d\\.(symp|beta|hosp|death)\\.effic")){###se é eficacia
      ex <- gregexpr("[0-9]+",j)
      vals <- as.numeric(regmatches(j,ex)[[1]])
      print(c(j,vals))
      if(length(vals) == 1){###aplica a todas as idades
        parameters[[as.character(j)]] <- as.numeric(samples[j])*rep(1,parameters$age.bins)
      }
      else if(length(vals) == 2){
        parameters[[str_extract(j,"vax\\d\\.(symp|beta|hosp|death)\\.effic\\.(P|A|C)")]][vals[2]] <- as.numeric(samples[j])
      }
      else if(length(vals) == 3){
        parameters[[str_extract(j,"vax\\d\\.(symp|beta|hosp|death)\\.effic\\.(P|A|C)")]][vals[2]:vals[3]] <- as.numeric(samples[j])
      }
    }
    else{
      parameters[[as.character(j)]] <- as.numeric(samples[j])
    }
  }
  print(parameters$vax1.symp.effic.C)
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.P"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                                 vax1.beta.effic.P,vax1.hosp.effic.P,
                                                                 vax1.death.effic.P,vax1.symp.effic.P))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$vax1.asymp.P = data.frame(agegroup = parameters$asymp$agegroup,vax1.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX1.EFFIC.CLIN))
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$vax1.ihr.P = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihr = parameters$ihr$ihr * (1-VAX1.EFFIC.SEVERE))
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$vax1.ihfr.P = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihfr = parameters$ihfr$ihfr * (1 - VAX1.EFFIC.DEATH))
    parameters$vax1.beta.P = (1- VAX1.EFFIC.BETA)
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.P"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                                 vax2.beta.effic.P,vax2.hosp.effic.P,
                                                                 vax2.death.effic.P,vax2.symp.effic.P))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    parameters$vax2.asymp.P = data.frame(agegroup = parameters$asymp$agegroup,vax2.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX2.EFFIC.CLIN))
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_w
    parameters$vax2.ihr.P = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihr = parameters$ihr$ihr * (1-VAX2.EFFIC.SEVERE))
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
    parameters$vax2.ihfr.P = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihfr = parameters$ihfr$ihfr * (1 - VAX2.EFFIC.DEATH))
    parameters$vax2.beta.P = (1- VAX2.EFFIC.BETA)
    
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.A"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                                 vax1.beta.effic.A,vax1.hosp.effic.A,
                                                                 vax1.death.effic.A,vax1.symp.effic.A))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$vax1.asymp.A = data.frame(agegroup = parameters$asymp$agegroup,vax1.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX1.EFFIC.CLIN))
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$vax1.ihr.A = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihr = parameters$ihr$ihr * (1-VAX1.EFFIC.SEVERE))
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$vax1.ihfr.A = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihfr = parameters$ihfr$ihfr * (1 - VAX1.EFFIC.DEATH))
    parameters$vax1.beta.A = (1- VAX1.EFFIC.BETA)
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.A"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                                 vax2.beta.effic.A,vax2.hosp.effic.A,
                                                                 vax2.death.effic.A,vax2.symp.effic.A))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    parameters$vax2.asymp.A = data.frame(agegroup = parameters$asymp$agegroup,vax2.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX2.EFFIC.CLIN))
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_w
    parameters$vax2.ihr.A = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihr = parameters$ihr$ihr * (1-VAX2.EFFIC.SEVERE))
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
    parameters$vax2.ihfr.A = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihfr = parameters$ihfr$ihfr * (1 - VAX2.EFFIC.DEATH))
    parameters$vax2.beta.A = (1- VAX2.EFFIC.BETA)
    
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.C"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                                 vax1.beta.effic.C,vax1.hosp.effic.C,
                                                                 vax1.death.effic.C,vax1.symp.effic.C))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$vax1.asymp.C = data.frame(agegroup = parameters$asymp$agegroup,vax1.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX1.EFFIC.CLIN))
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$vax1.ihr.C = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihr = parameters$ihr$ihr * (1-VAX1.EFFIC.SEVERE))
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$vax1.ihfr.C = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihfr = parameters$ihfr$ihfr * (1 - VAX1.EFFIC.DEATH))
    parameters$vax1.beta.C = (1- VAX1.EFFIC.BETA)
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.C"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                                 vax2.beta.effic.C,vax2.hosp.effic.C,
                                                                 vax2.death.effic.C,vax2.symp.effic.C))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    parameters$vax2.asymp.C = data.frame(agegroup = parameters$asymp$agegroup,vax2.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX2.EFFIC.CLIN))
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_w
    parameters$vax2.ihr.C = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihr = parameters$ihr$ihr * (1-VAX2.EFFIC.SEVERE))
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
    parameters$vax2.ihfr.C = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihfr = parameters$ihfr$ihfr * (1 - VAX2.EFFIC.DEATH))
    parameters$vax2.beta.C = (1- VAX2.EFFIC.BETA)
    
  }
  if(any(names(samples) == "prevalence")){
    ####reorganiza histórico
    PREV <- parameters$prevalence
    history <- parameters$history.P[,1:10] + parameters$history.P[,11:20]
    parameters$history.P <- cbind((1-PREV)*history,PREV*history)
    history <- parameters$history.A[,1:10] + parameters$history.A[,11:20]
    parameters$history.A <- cbind((1-PREV)*history,PREV*history)
    history <- parameters$history.C[,1:10] + parameters$history.C[,11:20]
    parameters$history.C <- cbind((1-PREV)*history,PREV*history)
    #######reorganiza condição inicial
    Y2 <- parameters$init.condition
    parameters$init.condition <- with(parameters,{
      Y <- Y2
      Y[Sindex] <- (1-PREV)*(Y2[Sindex]+Y2[Rindex])
      Y[Rindex] <- PREV*(Y2[Sindex]+Y2[Rindex])
      Y[SvPindex] <- (1-PREV)*(Y2[SvPindex]+Y2[RvPindex])
      Y[RvPindex] <- PREV*(Y2[SvPindex]+Y2[RvPindex])
      Y[SwPindex] <- (1-PREV)*(Y2[SwPindex]+Y2[RwPindex])
      Y[RwPindex] <- PREV*(Y2[SwPindex]+Y2[RwPindex])
      Y[SvAindex] <- (1-PREV)*(Y2[SvAindex]+Y2[RvAindex])
      Y[RvAindex] <- PREV*(Y2[SvAindex]+Y2[RvAindex])
      Y[SwAindex] <- (1-PREV)*(Y2[SwAindex]+Y2[RwAindex])
      Y[RwAindex] <- PREV*(Y2[SwAindex]+Y2[RwAindex])
      Y[SvCindex] <- (1-PREV)*(Y2[SvCindex]+Y2[RvCindex])
      Y[RvCindex] <- PREV*(Y2[SvCindex]+Y2[RvCindex])
      Y[SwCindex] <- (1-PREV)*(Y2[SwCindex]+Y2[RwCindex])
      Y[RwCindex] <- PREV*(Y2[SwCindex]+Y2[RwCindex])
      Y
    })
    
  }
  
  return(parameters)
}

