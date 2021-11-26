source("functions/vac_params_all_vac.R")
source("functions/markov_all_vac.R")
source("functions/sensitivity_analysis_C1_2.R")
source("functions/samples_to_parameters.R")
library(stringr)
library(doParallel)

cbind.fill <- function(a,b){
  if(is_empty(a)){
    return(b)
  }
  else{
    return(cbind(a,b))
  }
}
make.unique.2 = function(x, sep='.'){
  ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}})
}


mat = Matrix(matrix(0, length(Y2), length(Y2)), sparse = TRUE)
rownames(mat) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
colnames(mat) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
nmax = MAX.TIME.DAYS # Tempo de passos
# all.res <- data.frame()
n.cores <- 5
cl <- makeCluster(n.cores,type = "FORK")
registerDoParallel(cl)
parameters$vax.cov <- list(vax.cov = TRUE, cov = 0.8)

all.days <- c()
for(p in c(0.0001,0.0025,0.005,0.01)){
  for (AZ.days in c(84,77,70,63,56)) {
    #copia para garantir
    parameters2 <- parameters
    #historico AZ
    ### AstraZeneca ###
    
    VAX.INITIAL.STORAGE.NUM = 0 #6e6 # Number of READY vaccines on day one
    #VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
    VAX.PRODUCTION.RATE <- 10*rollout$astrazeneca
    
    # above numbers are for whole country, correcting (roughly) by population
    VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
    VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
    MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE
    
    dados_doses <- hist_D1D2 %>% filter(vacina == "AZ")
    dados_historico <- hist_D1 %>% filter(vacina == "AZ")
    dados_historico_2 <- dados_historico 
    dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
    dados_dose_1 <- dados_dose_1[1:10,]
    hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
    VAX.HISTORY <- as.numeric(hist$n)
    #recomputa solução
    optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                     MAX.VAC.RATE, AZ.days, SECOND.VAX.LOSS.FRAC,
                                     MAX.TIME.DAYS,VAX.HISTORY)
    vac.rate <- optimal_solution$OPT.VAX.RATE
    history <- optimal_solution$HISTORY
    # print(history)
    vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
    parameters2$vac.rate.A <- vac.rate
    parameters2$vac.rate.v2.A <- vac.rate.v2
    parameters2$vax.window.days.A <- AZ.days
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
    # print("ok")
    #############separa historico usando prevalencia inicial########################
    prev_hist <- cbind((1-parameters2$prevalence)*age_hist,parameters2$prevalence*age_hist)
    parameters2$history.A <- prev_hist
    parameters2$p <- p
    print(c(AZ.days,p))
    all.res<- foreach(i = 1:nrow(samples),.combine = "rbind",.packages = "stringr") %dopar% {
      parameters3 <- samples_to_parameters(samples[i,],parameters2)
      mat2<- add_probabilities(mat,parameters3)
      res <- run_markov(mat2,nmax,parameters3)
      res <- group_vaccines.2(res)
      res <- sum(res[nrow(res),c("D","Dv","Dw")])
      res
    }
    
    all.days <- rbind(all.days,data.frame(deaths = all.res,days = AZ.days,p = p))
  }
}
stopCluster(cl)
rownames(all.days) <- NULL

##########saving the solution
# write.csv(all.days,"comparison_ideal_impacto.csv",row.names = FALSE)


###### if don't want to run everything
all.days <- read.csv("comparison_ideal_impacto.csv")

df.2 <- c()

for(p_ in c(0.0001,0.0025,0.005,0.01)){
  ddf <- all.days %>% filter(days == 84 & p == p_)
  for(days_ in c(77,70,63,56)){
    dfilter <- all.days %>% filter(days == days_ & p == p_)
    vals <- data.frame(difference = dfilter$deaths-ddf$deaths, days = days_,p= p_)
    df.2 <- rbind(df.2,vals)
  }
}
df.2$days <- as.factor(df.2$days/7)
df.2 <- df.2 %>% mutate(grau = case_when(
  p == 0.0001 ~ "Very Low",
  p == 0.0025 ~ "Low",
  p == 0.005 ~ "Medium",
  p == 0.01 ~ "High"                                
))
df.2$grau <- factor(df.2$grau, levels = c("Very Low","Low","Medium","High"))

ggplot(df.2,aes(x = grau,y = -difference,group = interaction(grau,as.factor(days)),fill = as.factor(days))) + 
  geom_boxplot()+
  theme_bw()+geom_hline(yintercept = 0) +
  labs(x="Probability of Infection",y="Reduction in deaths compared to\n12 weeks of interval between doses",fill = "Weeks") +
  theme(text = element_text(size = 16))


