source("functions/vac_params_all_vac.R")
source("functions/markov_all_vac.R")
source("functions/sensitivity_analysis_C1_2.R")
source("functions/samples_to_parameters.R")
library(stringr)
library(doParallel)
library(beepr)
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

all.days <- c()
for(p in c(0.0001,0.0025,0.005,0.01)){
  for (AZ.days in c(84,70,56)) {
    #copia para garantir
    parameters2 <- parameters
    #historico AZ
    ### AstraZeneca ###
    
    VAX.INITIAL.STORAGE.NUM = 0 #6e6 # Number of READY vaccines on day one
    #VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
    VAX.PRODUCTION.RATE <- rollout$astrazeneca
    
    # above numbers are for whole country, correcting (roughly) by population
    VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
    VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
    MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE
    
    dados_doses <- hist_D1D2 %>% filter(vacina == "AZ")
    dados_historico <- hist_D1 %>% filter(vacina == "AZ")
    dados_historico_2 <- dados_historico 
    dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
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
    for(cov in seq(100,10,-10)){
      if(cov == 100){
        parameters2$vax.cov <- list(vax.cov = FALSE, cov = 1)
      }
      else{
        parameters2$vax.cov <- list(vax.cov = TRUE, cov = 0.01*cov)
      }
      print(c(AZ.days,p,cov))
      all.res<- foreach(i = 1:nrow(samples),.combine = "rbind",.packages = "stringr") %dopar% {
        parameters3 <- samples_to_parameters(samples[i,],parameters2)
        mat2<- add_probabilities(mat,parameters3)
        res <- run_markov(mat2,nmax,parameters3)
        res <- group_vaccines.2(res)
        res <- sum(res[nrow(res),c("D","Dv","Dw")])
        res
      }
      all.days <- rbind(all.days,data.frame(deaths = all.res,days = AZ.days,p = p,cov = cov))
    }
  }
}
stopCluster(cl)
rownames(all.days) <- NULL
beep(3)

####saving the results
# write.csv(all.days,file = "comparison_coverage_cenario_2.csv")



all.days <- read.csv("comparison_coverage_cenario_2.csv")

df.2 <- c()

for(p_ in c(0.0001,0.0025,0.005,0.01)){
  # for(days_ in c(84,70,56)){
    ddf <- all.days %>% filter(cov == 100 & p == p_)
    
    for(cov_ in seq(100,10,-10)){
      dfilter <- all.days %>% filter(cov == cov_ & p == p_) 
      vals <- data.frame(difference = dfilter$deaths-ddf$deaths, cov = cov_, days = days_, p = p_)
      df.2 <- rbind(df.2,vals)
    }
  # }
}
df.2 <- df.2 %>% mutate(grau = case_when(
                                p == 0.0001 ~ "Very low",
                                p == 0.0025 ~ "Low",
                                p == 0.005 ~ "Medium",
                                p == 0.01 ~ "High"                                
))
df.2$grau <- factor(df.2$grau, levels = c("Very low","Low","Medium","High"))
ggplot(df.2,aes(x = cov, y = difference, fill = as.factor(days), group = interaction(as.factor(cov),as.factor(days_))))+
facet_wrap(~p,scales="free",labeller = labeller(p = label_both)) +
 geom_boxplot(aes(group = interaction(as.factor(cov),as.factor(days))))+
 theme_bw()+geom_hline(yintercept = 0) +
 labs(x="Cobertura mínima",y="Excesso de óbitos comparado a\n84 dias de intervalo e 100% de cobertura",fill = "Dias") +
 theme(text = element_text(size = 16))
ggplot(df.2,aes(x = cov, y = difference,  group = as.factor(cov)))+
  facet_wrap(~grau,scales="free") +
  geom_boxplot(aes(group = as.factor(cov)))+
  theme_bw()+geom_hline(yintercept = 0) +
  labs(x="Minimum vaccine coverage per age group",y="Excess of deaths compared to\n100% of coverage",fill = "Dias") +
  theme(text = element_text(size = 16))


####old code, ignore it
# all.days <- read_csv("comparison_coverage_cenario_2.csv")
# df.2 <- c()
# for(p_ in c(0.0001,0.0025,0.005,0.01)){
#   ddf <- all.days %>% filter(days == 84 & p == p_)
#   for(days_ in c(70,56)){
#     dfilter <- all.days %>% filter(days == days_ & p == p_)
#     vals <- data.frame(difference = dfilter$deaths-ddf$deaths, days = days_,p= p_)
#     print(mean(vals$difference))
#     print(quantile(vals$difference,c(0.025,0.5,0.975)))
#     df.2 <- rbind(df.2,vals)
#   }
# }
# df.2$days <- as.factor(df.2$days)
# ggplot(df.2,aes(x = difference, fill = days)) + facet_grid(rows = vars(p),cols = vars(days),scales="free") +
#   geom_histogram(position = "identity",bins = 50) + geom_vline(xintercept = 0)+
#   theme(legend.position = "none",text = element_text(size=16)) + labs(x = "Excesso de óbitos em relação à estratégia de\n84 dias de intervalo", y = "freq")

# 
# for(cov in seq(100,10,-10)){
#   if(cov == 100){
#     parameters$vax.cov <- list(vax.cov = FALSE, cov = 1)
#   }
#   else{
#     parameters$vax.cov <- list(vax.cov = TRUE, cov = 0.01*cov)
#   }
#   print(parameters$vax.cov$cov)
#   all.res<- foreach(i = 1:nrow(samples),.combine = "rbind",.packages = "stringr") %dopar% {
#     parameters2 <- samples_to_parameters(samples[i,],parameters)
#     mat2<- add_probabilities(mat,parameters2)
#     res <- run_markov(mat2,nmax,parameters2)
#     res <- group_vaccines.2(res)
#     res <- sum(res[nrow(res),c("D","Dv","Dw")])
#     res
#   }
#   print(tail(all.res))
#   df <- rbind(df,data.frame(deaths = all.res[,1],cov = cov))
# }
# stopCluster(cl)
# 
# rownames(df) <- NULL
############################################################


