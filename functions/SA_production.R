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
for(prod in c(1.0,1.25,1.50,1.75,2.00)){
  for(p in c(0.0001,0.0025,0.005,0.01)){
    for (AZ.days in c(84,77,70,63,56)) {
     #copia para garantir
      parameters2 <- parameters
      #historico AZ
      ### AstraZeneca ###
    
      VAX.INITIAL.STORAGE.NUM = 0 #6e6 # Number of READY vaccines on day one
      #VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
      VAX.PRODUCTION.RATE <- prod*rollout$astrazeneca
    
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
      print(c(AZ.days,p,prod))
      all.res<- foreach(i = 1:nrow(samples),.combine = "rbind",.packages = "stringr") %dopar% {
        parameters3 <- samples_to_parameters(samples[i,],parameters2)
        mat2<- add_probabilities(mat,parameters3)
        res <- run_markov(mat2,nmax,parameters3)
        res <- group_vaccines.2(res)
        res <- sum(res[nrow(res),c("D","Dv","Dw")])
        res
      }
    
      all.days <- rbind(all.days,data.frame(deaths = all.res,days = AZ.days,p = p,prod = prod))
    }
  }
}
stopCluster(cl)
rownames(all.days) <- NULL
#####write the solution######
# write.csv(all.days,"comparison_production_impacto.csv",row.names = FALSE)

#######reading the solution to avoid running everything######
all.days <- read.csv("comparison_production_impacto.csv")
df.2 <- c()
# for(prod_ in c(1,1.25,1.5,1.75,2.0)){
  for(p_ in c(0.0001,0.0025,0.005,0.01)){
    ddf <- all.days %>% filter(days == 84 & p == p_)
    for(days_ in c(77,70,63,56)){
      dfilter <- all.days %>% filter(days == days_ & p == p_)
      vals <- data.frame(difference = dfilter$deaths-ddf$deaths, days = days_,p= p_,prod = dfilter$prod)
      print(mean(vals$difference))
      print(quantile(vals$difference,c(0.025,0.5,0.975)))
      df.2 <- rbind(df.2,vals)
    }
  }
# }
df.2$days <- as.factor(df.2$days/7)
# ggplot(df.2,aes(x = difference, fill = days)) + facet_grid(rows = vars(p),cols = vars(days),scales="free") +
#   geom_histogram(position = "identity",bins = 50) + geom_vline(xintercept = 0)+
#   theme(legend.position = "none",text = element_text(size=16)) + labs(x = "Excesso de óbitos em relação à estratégia de\n84 dias de intervalo", y = "freq")
ggplot(df.2,aes(x = as.factor(p),y = difference,group = interaction(as.factor(p),as.factor(days)),fill = as.factor(days))) + 
  geom_boxplot()+
  theme_bw()+geom_hline(yintercept = 0) +
  labs(x="Probabilidade de infecção",y="Excesso de óbitos comparado a\n84 dias de intervalo",fill = "Semanas") +
  theme(text = element_text(size = 16))

# df <- all.days
df.2 <- df.2 %>% mutate(grau = case_when(
  p == 0.0001 ~ "Very Low",
  p == 0.0025 ~ "Low",
  p == 0.005 ~ "Medium",
  p == 0.01 ~ "High"                                
))
df.2$grau <- factor(df.2$grau, levels = c("Very Low","Low","Medium","High"))
df.2 <- df.2 %>% mutate(aumento = case_when(
  prod == 1 ~ "0%",
  prod == 1.25 ~ "25%",
  prod == 1.5 ~ "50%",
  prod == 1.75 ~ "75%",
  prod == 2.00 ~ "100%"
))
df.2$aumento <- factor(df.2$aumento, levels = c("0%","25%","50%","75%","100%"))
df.3 <- df.2 %>% filter(grau == "Alta")
ggplot(df.2,aes(x = aumento,y = -difference,group = interaction(as.factor(aumento),as.factor(days)),fill = as.factor(days))) + 
  facet_wrap(~grau,scales="free") + geom_boxplot()+
  theme_bw()+geom_hline(yintercept = 0) +
  labs(x="Increase in production",y="Reduction in number of deaths compared to\n12 weeks of interval between doses",fill = "Weeks") +
  theme(text = element_text(size = 16))
# for(p_ in c(0.0001,0.0025,0.005,0.01)){
#   for(days_ in c(70,56)){
#     dfilter <- df.2 %>% filter(days == days_ & p == p_)
#     print(c(p_,days_,mean(dfilter$difference)))
#     # print(quantile(dfilter$difference,c(0.025,0.5,0.975)))
#   }
# }

# names(all.res)<-make.unique.2(names(all.res),sep="")
# # epiclasses = c(c("S", "E", "A", "I", "H", "R", "D"),
# #                paste0(c("S", "E", "A", "I", "H", "R", "D"),"vA"),
# #                paste0(c("S", "E", "A", "I", "H", "R", "D"),"vP"),
# #                paste0(c("S", "E", "A", "I", "H", "R", "D"),"vC"),
# #                paste0(c("S", "E", "A", "I", "H", "R", "D"),"wA"),
# #                paste0(c("S", "E", "A", "I", "H", "R", "D"),"wP"),
# #                paste0(c("S", "E", "A", "I", "H", "R", "D"),"wC"))
# epiclasses = c(c("S", "E", "A", "I", "H", "R", "D"),
#                paste0(c("S", "E", "A", "I", "H", "R", "D"),"v"),
#                paste0(c("S", "E", "A", "I", "H", "R", "D"),"w"))
# 
# mean_values.2 <- c()
# for(classe in epiclasses){
#   classes <- paste0(classe,1:(nrow(samples)))
#   vals <- all.res %>% select(classes) 
#   aggregated <- data.frame(time = all.res[,"time1"],
#                            classe = classe,
#                            mean = apply(vals, 1, mean),
#                            sd = apply(vals, 1, sd),
#                            q2.5 = apply(vals, 1, quantile, 0.025),
#                            q50 = apply(vals, 1, quantile, 0.5),
#                            q97.5 = apply(vals, 1, quantile, 0.975))
#   mean_values.2 <- rbind(mean_values.2,aggregated)
# }
# mean_values.2$classe <- factor(mean_values.2$classe,levels = epiclasses)
# df <- mean_values.2 %>% filter(str_detect(classe,"D|Dv|Dw"))
# g <- ggplot(df,aes(x = time,y = q50,group = classe,fill = classe)) +
#   geom_line(aes(color = classe)) + geom_ribbon(aes(ymin = q2.5,ymax = q97.5),alpha = 0.2)+
#   scale_color_viridis_d() + scale_fill_viridis_d()+theme_minimal()
