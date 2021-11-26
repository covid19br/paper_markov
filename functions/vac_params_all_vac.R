source("functions/epi_params_all_vac.R")
source("functions/opt_vax_rate_history.R")
source("functions/VAX_DISTR_RATE.R")
library(dplyr)

## functions ##

#### convert efficacy from studies to parameters
OBSERVED.EFFICACY.TO.PARAMETERS <- function(HOSP.PROP,MORT.PROP,ASYMP.PROP,SUSCEP.EFF,HOSP.EFF,MORT.EFF,SYMP.EFF){
  SUSCEP.EFF.PAR <- SUSCEP.EFF
  HOSP.EFF.PAR <- 1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)
  MORT.EFF.PAR <- 1 - (1-MORT.EFF)/(1-HOSP.EFF)
  SYMP.EFF.PAR <- 1 - ((1-SYMP.EFF)*(HOSP.PROP + (1 - HOSP.PROP)*(1-ASYMP.PROP))  - (1 - HOSP.EFF)*HOSP.PROP)/
    ((1-SUSCEP.EFF)*(1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)*HOSP.PROP)*(1-ASYMP.PROP))
  return(list(suscep = SUSCEP.EFF.PAR,hosp = HOSP.EFF.PAR,death = MORT.EFF.PAR,symp = SYMP.EFF.PAR))
}


## global parameters ##

POP.CITY.REL.FRAC = 1 #population of SP city compared to Brazil
MAX.VAC.RATE = 5e6 # Max number of vaccine applications per day
MAX.VAC.RATE2 <- MAX.VAC.RATE
#####PREVALENCIA ARTIFICIAL
PREV <- 0.5
# Fraction of people that take the first but not the second dose \theta
SECOND.VAX.LOSS.FRAC = 0.1
parameters$desist <- SECOND.VAX.LOSS.FRAC

# tabela de entrega de vacinas
rollout <- read.csv("./DATA/rollout_vacinas.csv")

data_base <- as.Date("2021-08-09")
rollout$date <- as.Date(rollout$date)
rollout <- rollout %>% filter(date >= data_base)
rollout[,2] <- rollout[,2]#######apagar o 10, não esquecer
date_range <- as.integer(diff(range(rollout$date)))

# time until end of vaccination schedule, in days
MAX.TIME.DAYS = min(300, date_range)

# tabelas de históricos
estados <- c("RO", "AC", "AM", "RR", "PA", "AP", "TO", "MA", "PI", "CE", "RN",
             "PB", "PE", "AL", "SE", "BA", "MG", "ES", "RJ", "SP", "PR", "SC",
             "RS", "MS", "MT", "GO", "DF")
hist_D1 <- list()
hist_D2 <- list()
hist_D1_acum <- list()
for (estado in estados) {
  hist_D1[[estado]] <- read.csv(paste0("DATA/historico/historico_D1_", estado, ".csv"))
  hist_D2[[estado]] <- read.csv(paste0("DATA/historico/historico_D2_", estado, ".csv"))
  hist_D1_acum[[estado]] <- hist_D1[[estado]] %>% 
    group_by(vacina, agegroup, tempo_doses) %>%
    summarize(n = sum(n)) %>%
    as.data.frame()
}
# agrupando BR todo
hist_D1 <- bind_rows(hist_D1, .id = "estado")
hist_D1 <- hist_D1 %>%
  group_by(vacina, agegroup, tempo_doses) %>%
  summarize(n = sum(n)) %>%
  as.data.frame()
hist_D1_acum <- bind_rows(hist_D1_acum, .id = "estado")
hist_D1_acum <- hist_D1_acum %>%
  group_by(vacina, agegroup) %>%
  summarize(dose1 = sum(n)) %>%
  as.data.frame()
hist_D2 <- bind_rows(hist_D2, .id = "estado")
hist_D2 <- hist_D2 %>%
  group_by(vacina, agegroup) %>%
  summarize(dose2 = sum(n)) %>%
  as.data.frame()
hist_D1D2 <- full_join(hist_D1_acum, hist_D2, by = c("vacina", "agegroup"))
hist_D1D2 <- hist_D1D2 %>%
  filter(!is.na(agegroup))
hist_D1D2$dose2[is.na(hist_D1D2$dose2)]<-0
## vaccine-specific parameters ##
# hist_D1$n<- hist_D1*POP.CITY.REL.FRAC
# hist_D1D2[,c("dose1","dose2")] <- hist_D1D2[,c("dose1","dose2")]*POP.CITY.REL.FRAC

list.parameters <- list()


### Pfizer ###

VAX.INITIAL.STORAGE.NUM = 0 # 6e6 # Number of READY vaccines on day one
#VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
VAX.PRODUCTION.RATE <- rollout$pfizer

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE

# Time window between first and second vaccines
VAX.WINDOW.DAYS = 56
#VAX.WINDOW.DAYS = 84 # 12 weeks

#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
parameters$vax2.beta.effic.P = 0.90*rep(1,age.bins)
parameters$vax2.symp.effic.P = 0.94*rep(1,age.bins)
parameters$vax2.hosp.effic.P = 0.87*rep(1,age.bins)
parameters$vax2.death.effic.P = 0.98*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
parameters$vax1.beta.effic.P   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.P
parameters$vax1.symp.effic.P   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.P
parameters$vax1.hosp.effic.P   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.P
parameters$vax1.death.effic.P  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.P

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                             vax1.beta.effic.P,vax1.hosp.effic.P,
                                                             vax1.death.effic.P,vax1.symp.effic.P))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                             vax2.beta.effic.P,vax2.hosp.effic.P,
                                                             vax2.death.effic.P,vax2.symp.effic.P))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$vax1.asymp.P = data.frame(agegroup = parameters$asymp$agegroup,vax1.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX1.EFFIC.CLIN))
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$vax1.ihr.P = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihr = parameters$ihr$ihr * (1-VAX1.EFFIC.SEVERE))
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$vax1.ihfr.P = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihfr = parameters$ihfr$ihfr * (1 - VAX1.EFFIC.DEATH))
parameters$vax1.beta.P = (1- VAX1.EFFIC.BETA)
parameters$vax.window.days.P <- VAX.WINDOW.DAYS

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_w
parameters$vax2.asymp.P = data.frame(agegroup = parameters$asymp$agegroup,vax2.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX2.EFFIC.CLIN))
# Fraction of severe cases in symptomatic cases (IHR) \sigma_w
parameters$vax2.ihr.P = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihr = parameters$ihr$ihr * (1-VAX2.EFFIC.SEVERE))
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
parameters$vax2.ihfr.P = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihfr = parameters$ihfr$ihfr * (1 - VAX2.EFFIC.DEATH))
parameters$vax2.beta.P = (1- VAX2.EFFIC.BETA)

######mexendo com historico
Y2 <- Y
dados_doses <- hist_D1D2 %>% filter(vacina == "Pfizer")
dados_historico <- hist_D1 %>% filter(vacina == "Pfizer")
dados_historico_2 <- dados_historico
dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
# dados_dose_1 <- rbind(c(1,0),c(2,0),dados_dose_1)

Y2[Sindex] <- Y2[Sindex] - dados_doses$dose2
Y2[SwPindex] <- Y2[SwPindex] + dados_doses$dose2
Y2[Sindex] <- Y2[Sindex] - dados_dose_1$n
Y2[SvPindex] <- Y2[SvPindex] + dados_dose_1$n
##########################
names(Y2) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
VAX.HISTORY <- as.numeric(hist$n)

optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
parameters$vac.rate.P <- vac.rate
parameters$vac.rate.v2.P <- vac.rate.v2
df <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
                             MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
df$Window <- "Pfizer"
plot_vac_schedule.3(df)

#############separa historico em faixa etaria###################################
dados_historico_2<-dados_historico_2[with(dados_historico_2,order(-tempo_doses,-agegroup)),]
age_hist <- matrix(0,nrow = MAX.TIME.DAYS+1,ncol = age.bins)
j <- 2
for(i in 1:nrow(dados_historico_2)){
  val <- dados_historico_2[i,]
  while(val$n > 0){
    if(sum(age_hist[j]) >= history[j]){ ###se já preencheu a linha
      if(j < MAX.TIME.DAYS+1){##e j não está no fim, incrementa
        j <- j+1
      }
      else{###se j está no fim, quebra o loop, colocar um warning aqui?
        break
      }
    }
    if(sum(age_hist[j,]) + val$n <= history[j]){###se cabe todas as pessoas
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val$n
      val$n <- 0
    }
    else{####se não cabe
      val2 <- history[j] - sum(age_hist[j,])
      # print(val2)
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val2
      val$n <- val$n-val2
      j <- j+1
    }
  }###fim while
}
#############separa historico usando prevalencia inicial########################
prev_hist <- cbind((1-PREV)*age_hist,PREV*age_hist)
parameters$history.P <- prev_hist
### AstraZeneca ###

VAX.INITIAL.STORAGE.NUM = 0 #6e6 # Number of READY vaccines on day one
#VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
VAX.PRODUCTION.RATE <- rollout$astrazeneca

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE2

# Time window between first and second vaccines
VAX.WINDOW.DAYS = 56
#VAX.WINDOW.DAYS = 84 # 12 weeks

#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
parameters$vax2.beta.effic.A = 0.6*rep(1,age.bins)
parameters$vax2.symp.effic.A = 0.81*rep(1,age.bins)
parameters$vax2.hosp.effic.A = 0.90*rep(1,age.bins)
parameters$vax2.death.effic.A = 0.95*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.4
parameters$vax1.beta.effic.A   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.A
parameters$vax1.symp.effic.A   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.A
parameters$vax1.hosp.effic.A   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.A
parameters$vax1.death.effic.A  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.A

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                             vax1.beta.effic.A,vax1.hosp.effic.A,
                                                             vax1.death.effic.A,vax1.symp.effic.A))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                             vax2.beta.effic.A,vax2.hosp.effic.A,
                                                             vax2.death.effic.A,vax2.symp.effic.A))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$vax1.asymp.A = data.frame(agegroup = parameters$asymp$agegroup,vax1.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX1.EFFIC.CLIN))
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$vax1.ihr.A = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihr = parameters$ihr$ihr * (1-VAX1.EFFIC.SEVERE))
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$vax1.ihfr.A = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihfr = parameters$ihfr$ihfr * (1 - VAX1.EFFIC.DEATH))
parameters$vax1.beta.A = (1- VAX1.EFFIC.BETA)

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_w
parameters$vax2.asymp.A = data.frame(agegroup = parameters$asymp$agegroup,vax2.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX2.EFFIC.CLIN))
# Fraction of severe cases in symptomatic cases (IHR) \sigma_w
parameters$vax2.ihr.A = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihr = parameters$ihr$ihr * (1-VAX2.EFFIC.SEVERE))
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
parameters$vax2.ihfr.A = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihfr = parameters$ihfr$ihfr * (1 - VAX2.EFFIC.DEATH))
parameters$vax2.beta.A = (1- VAX2.EFFIC.BETA)
parameters$vax.window.days.A <- VAX.WINDOW.DAYS
######mexendo com historico
dados_doses <- hist_D1D2 %>% filter(vacina == "AZ")
dados_historico <- hist_D1 %>% filter(vacina == "AZ")
dados_historico_2 <- dados_historico 
dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
# dados_dose_1 <- rbind(c(1,0),c(2,0),dados_dose_1)
dados_dose_1 <- dados_dose_1[1:10,]
Y2[Sindex] <- Y2[Sindex] - dados_doses$dose2
Y2[SwAindex] <- Y2[SwAindex] + dados_doses$dose2
Y2[Sindex] <- Y2[Sindex] - dados_dose_1$n
Y2[SvAindex] <- Y2[SvAindex] + dados_dose_1$n
######dados outras vacinas
# dados_doses <- hist_D1D2 %>% filter(vacina != "AZ") %>% group_by(agegroup) %>% summarize(n = sum(dose1+dose2))
# Y2[Sindex] <- Y2[Sindex] - dados_doses$n

hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
VAX.HISTORY <- as.numeric(hist$n)

optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
parameters$vac.rate.A <- vac.rate
parameters$vac.rate.v2.A <- vac.rate.v2
df.2 <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
                             MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
df.2$Window <- "AstraZeneca"
plot_vac_schedule.3(df.2)
#############separa historico em faixa etaria###################################
dados_historico_2<-dados_historico_2[with(dados_historico_2,order(-tempo_doses,-agegroup)),]
age_hist <- matrix(0,nrow = MAX.TIME.DAYS+1,ncol = age.bins)
j <- 2
for(i in 1:nrow(dados_historico_2)){
  val <- dados_historico_2[i,]
  while(val$n > 0){
    if(sum(age_hist[j]) >= history[j]){ ###se já preencheu a linha
      if(j < MAX.TIME.DAYS+1){##e j não está no fim, incrementa
        j <- j+1
      }
      else{###se j está no fim, quebra o loop, colocar um warning aqui?
        break
      }
    }
    if(sum(age_hist[j,]) + val$n <= history[j]){###se cabe todas as pessoas
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val$n
      val$n <- 0
    }
    else{####se não cabe
      val2 <- history[j] - sum(age_hist[j,])
      # print(val2)
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val2
      val$n <- val$n-val2
      j <- j+1
    }
  }###fim while
}
#############separa historico usando prevalencia inicial########################
prev_hist <- cbind((1-PREV)*age_hist,PREV*age_hist)
parameters$history.A <- prev_hist

### CoronaVac ###
VAX.INITIAL.STORAGE.NUM = 0 # 6e6 # Number of READY vaccines on day one
#VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
VAX.PRODUCTION.RATE <- rollout$coronavac

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE2

# Time window between first and second vaccines
VAX.WINDOW.DAYS = 28
#VAX.WINDOW.DAYS = 84 # 12 weeks

#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
parameters$vax2.beta.effic.C = 0.0*rep(1,age.bins)
parameters$vax2.symp.effic.C = 0.5*rep(1,age.bins)
parameters$vax2.hosp.effic.C = 0.83*rep(1,age.bins)
parameters$vax2.death.effic.C = 0.95*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
parameters$vax1.beta.effic.C   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.C
parameters$vax1.symp.effic.C   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.C
parameters$vax1.hosp.effic.C   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.C
parameters$vax1.death.effic.C  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.C

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                             vax1.beta.effic.C,vax1.hosp.effic.C,
                                                             vax1.death.effic.C,vax1.symp.effic.C))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr$ihr,ihfr$ihfr,asymp$asymp,
                                                             vax2.beta.effic.C,vax2.hosp.effic.C,
                                                             vax2.death.effic.C,vax2.symp.effic.C))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$vax1.asymp.C = data.frame(agegroup = parameters$asymp$agegroup,vax1.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX1.EFFIC.CLIN))
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$vax1.ihr.C = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihr = parameters$ihr$ihr * (1-VAX1.EFFIC.SEVERE))
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$vax1.ihfr.C = data.frame(agegroup = parameters$ihr$agegroup,vax1.ihfr = parameters$ihfr$ihfr * (1 - VAX1.EFFIC.DEATH))
parameters$vax1.beta.C = (1- VAX1.EFFIC.BETA)

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_w
parameters$vax2.asymp.C = data.frame(agegroup = parameters$asymp$agegroup,vax2.asymp = 1 - (1 - parameters$asymp$asymp) * (1-VAX2.EFFIC.CLIN))
# Fraction of severe cases in symptomatic cases (IHR) \sigma_w
parameters$vax2.ihr.C = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihr = parameters$ihr$ihr * (1-VAX2.EFFIC.SEVERE))
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
parameters$vax2.ihfr.C = data.frame(agegroup = parameters$ihr$agegroup,vax2.ihfr = parameters$ihfr$ihfr * (1 - VAX2.EFFIC.DEATH))
parameters$vax2.beta.C = (1- VAX2.EFFIC.BETA)
parameters$vax.window.days.C <- VAX.WINDOW.DAYS
######mexendo com historico
dados_doses <- hist_D1D2 %>% filter(vacina == "Coronavac")
dados_historico <- hist_D1 %>% filter(vacina == "Coronavac")
dados_historico_2 <- dados_historico
dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
# dados_dose_1 <- rbind(c(1,0),c(2,0),dados_dose_1)

Y2[Sindex] <- Y2[Sindex] - dados_doses$dose2
Y2[SwCindex] <- Y2[SwCindex] + dados_doses$dose2
Y2[Sindex] <- Y2[Sindex] - dados_dose_1$n
Y2[SvCindex] <- Y2[SvCindex] + dados_dose_1$n
# ######dados outras vacinas
# dados_doses <- hist_D1D2 %>% filter(vacina != "Coronavac") %>% group_by(agegroup) %>% summarize(n = sum(dose1+dose2))
# Y2[Sindex] <- Y2[Sindex] - dados_doses$n

hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
VAX.HISTORY <- as.numeric(hist$n)

optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
parameters$vac.rate.C <- vac.rate
parameters$vac.rate.v2.C <- vac.rate.v2
df.3 <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
                              MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
df.3$Window <- "CoronaVac"
plot_vac_schedule.3(df.3)
#############separa historico em faixa etaria###################################
dados_historico_2<-dados_historico_2[with(dados_historico_2,order(-tempo_doses,-agegroup)),]
age_hist <- matrix(0,nrow = MAX.TIME.DAYS+1,ncol = age.bins)
j <- 2
for(i in 1:nrow(dados_historico_2)){
  val <- dados_historico_2[i,]
  while(val$n > 0){
    if(sum(age_hist[j]) >= history[j]){ ###se já preencheu a linha
      if(j < MAX.TIME.DAYS+1){##e j não está no fim, incrementa
        j <- j+1
      }
      else{###se j está no fim, quebra o loop, colocar um warning aqui?
        break
      }
    }
    if(sum(age_hist[j,]) + val$n <= history[j]){###se cabe todas as pessoas
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val$n
      val$n <- 0
    }
    else{####se não cabe
      val2 <- history[j] - sum(age_hist[j,])
      # print(val2)
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val2
      val$n <- val$n-val2
      j <- j+1
    }
  }###fim while
}
#############separa historico usando prevalencia inicial########################
prev_hist <- cbind((1-PREV)*age_hist,PREV*age_hist)
parameters$history.C <- prev_hist

###############################################################################
Y2[Rindex] <- PREV*Y2[Sindex]
Y2[Sindex] <- (1-PREV)*Y2[Sindex]
Y2[RwPindex] <- PREV*Y2[SwPindex]
Y2[SwPindex] <- (1-PREV)*Y2[SwPindex]
Y2[RvPindex] <- PREV*Y2[SvPindex]
Y2[SvPindex] <- (1-PREV)*Y2[SvPindex]
Y2[RwAindex] <- PREV*Y2[SwAindex]
Y2[SwAindex] <- (1-PREV)*Y2[SwAindex]
Y2[RvAindex] <- PREV*Y2[SvAindex]
Y2[SvAindex] <- (1-PREV)*Y2[SvAindex]
Y2[RwCindex] <- PREV*Y2[SwCindex]
Y2[SwCindex] <- (1-PREV)*Y2[SwCindex]
Y2[RvCindex] <- PREV*Y2[SvCindex]
Y2[SvCindex] <- (1-PREV)*Y2[SvCindex]
parameters$init.condition <-  Y2
###############################################################################
parameters$vax.cov <- list(vax.cov = FALSE, cov = 1)

parameters$prevalence <- PREV

