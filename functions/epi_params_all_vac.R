require(Matrix)
require(plyr)
require(dplyr)
require(tidyverse)


############PREPARA DISTR ETARIA
age.bins <- 10 #numero de faixas etarias #########automatizar os breaks baseado no num de faixas
DistrEtaria <- as.data.frame(read.csv("DATA/2. DistrEtaria2020.csv"))
DistrEtaria <- DistrEtaria[2:nrow(DistrEtaria),]
colnames(DistrEtaria)[1] <- "agefloor"
popstruc <- data.frame(agefloor = as.numeric(DistrEtaria$agefloor),n = rowSums(DistrEtaria[,2:ncol(DistrEtaria)]))
popstruc0 <- popstruc
popstruc <- popstruc %>% mutate(agegroup = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,90,Inf),right = F,labels = F))
popstruc <- popstruc %>% select(c("agegroup","n")) %>% group_by(agegroup) %>% summarise(n = sum(n))
POP.TOTAL.NUM <- sum(popstruc$n)
##############PREPARA IHR
ihr <- data.frame(read.csv("DATA/14. ihr.csv"))
ihr[,2] <- as.numeric(ihr[,2])
ihr[,2]<-ihr[,2]/100   # csv data is in percentages
ihr <- ihr %>% mutate(agegroup = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,90,Inf),right = F,labels = F))
ihr$weighted_ihr <- ihr[,2]*popstruc0[,2]
ihr$age_distr <- popstruc0[,2]
ihr <- ihr %>% group_by(agegroup) %>%
  summarise_at(vars(weighted_ihr,age_distr),sum) %>%
  mutate(ihr = weighted_ihr/age_distr) %>% 
  select(c("agegroup","ihr"))
ihr <- as.data.frame(ihr)
###############PREPARA SINTOMATICOS
symp <- data.frame(read.csv("DATA/symptomatics_per_age.csv"))
symp[,2] <- as.numeric(symp[,2])
symp <- symp %>% mutate(agegroup = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,90,Inf),right = F,labels = F))
symp$weighted_symp <- symp[,2]*popstruc0[,2]
symp$age_distr <- popstruc0[,2]
symp <- symp %>% group_by(agegroup) %>%
  summarise_at(vars(weighted_symp,age_distr),sum) %>%
  mutate(symp = weighted_symp/age_distr) %>% 
  select(c("agegroup","symp"))
symp <- as.data.frame(symp)
##############PREPARA OBITOS
# ihfr_estados <- read.csv('DATA/ihfr.csv')
# DEATH.FRAC <- ihfr_estados %>%
#   filter(sg_uf == "SP") %>%
#   arrange(age_clas) %>%
#   .$ihfr.covid
ihfr_brasil <- read.csv("DATA/ihfrBrasil.csv")
###ajeitando na mao pq é mais facil
ihfr <- data.frame(agegroup = 1:10)
ihfr$ihfr <- ihfr_brasil$ihfr
# ihfr[1:2,2] <-DEATH.FRAC[1]
# ihfr[3:7,2] <-DEATH.FRAC[2]
# ihfr[7:10,2] <-DEATH.FRAC[3]
######################################
######################################
##### A === AstraZeneca
##### P === Pfizer
##### C === CoronaVac
#####################################
## Definir estados
states = c(c("S", "E", "A", "I", "H", "R", "D"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"vA"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"vP"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"vC"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"wA"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"wP"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"wC"))
for (v in 1:length(states)){
  eval(parse(text=paste0(states[v], 'index <- seq(', 1+age.bins*(v-1), ',', v*age.bins, ')')))
  eval(parse(text=paste0('init', states[v], ' <- 0*popstruc[,2]')))
}
initS <- popstruc[,2]
Y <- eval(parse(text=paste0('c(', paste(paste0('init', states), collapse=','), ')')))
Y <- unlist(Y)
Y <- as.numeric(Y)
names(Y) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
# _v: primeira dose
# _w: segunda dose

## Definir matriz de transição
mat = Matrix(matrix(0, length(Y), length(Y)), sparse = TRUE)
rownames(mat) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
colnames(mat) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))

parameters <- list()
parameters[['ihr']]<- ihr
parameters[["ihfr"]] <- ihfr
parameters[['symp']] <- symp
parameters[['asymp']] <- symp
parameters$asymp$symp <- 1-symp$symp
names(parameters$asymp) <- c("agegroup","asymp")
for (v in 1:length(states)){
  eval(parse(text=paste0("parameters[['",states[v], "index']] <- ",states[v],"index")))
}
nu <- read.csv("DATA/hosp_time_br.csv")
parameters[["p"]] <- 0.5
parameters[["gamma"]] <- 1.0/5.8
parameters[["nu"]] <- 1.0/nu$mean
parameters[["age.bins"]] <- age.bins

