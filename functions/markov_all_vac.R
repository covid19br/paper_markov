require(Matrix)
require(magrittr)

add_probabilities <- function(mat,params){
  mat2 <- mat
  mat2 <- within(c(params,mat2),{
    S_to_E <- p
    E_to_IAH <- 1-exp(-gamma)
    IA_to_R <- 1-exp(-nu)
    H_to_RD <- 1-exp(-nu) ###assumindo que é o mesmo valor provisoriamente
    ######################saindo de S############################
    mat2[cbind(Sindex,Eindex)] <- S_to_E  
    ######################saindo de E############################
    mat2[cbind(Eindex,Iindex)] <- (1-ihr[,2])*symp[,2]*E_to_IAH
    mat2[cbind(Eindex,Aindex)] <- (1-ihr[,2])*(1-symp[,2])*E_to_IAH
    mat2[cbind(Eindex,Hindex)] <- ihr[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(Iindex,Rindex)] <- IA_to_R
    mat2[cbind(Aindex,Rindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(Hindex,Rindex)] <- (1-ihfr[,2])*H_to_RD
    mat2[cbind(Hindex,Dindex)] <- ihfr[,2]*H_to_RD
    ########################PFIZER###############################
    #####################primeira dose###########################
    ######################saindo de S############################
    mat2[cbind(SvPindex,EvPindex)] <- vax1.beta.P*S_to_E  
    ######################saindo de E############################
    mat2[cbind(EvPindex,IvPindex)] <- (1-vax1.ihr.P[,2])*(1-vax1.asymp.P[,2])*E_to_IAH
    mat2[cbind(EvPindex,AvPindex)] <- (1-vax1.ihr.P[,2])*(vax1.asymp.P[,2])*E_to_IAH
    mat2[cbind(EvPindex,HvPindex)] <- vax1.ihr.P[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(IvPindex,RvPindex)] <- IA_to_R
    mat2[cbind(AvPindex,RvPindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(HvPindex,RvPindex)] <- (1-vax1.ihfr.P[,2])*H_to_RD
    mat2[cbind(HvPindex,DvPindex)] <- vax1.ihfr.P[,2]*H_to_RD
    
    #####################segunda dose###########################
    ######################saindo de S############################
    mat2[cbind(SwPindex,EwPindex)] <- vax2.beta.P*S_to_E  ####sem vacinas ainda
    ######################saindo de E############################
    mat2[cbind(EwPindex,IwPindex)] <- (1-vax2.ihr.P[,2])*(1-vax2.asymp.P[,2])*E_to_IAH
    mat2[cbind(EwPindex,AwPindex)] <- (1-vax2.ihr.P[,2])*vax2.asymp.P[,2]*E_to_IAH
    mat2[cbind(EwPindex,HwPindex)] <- vax2.ihr.P[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(IwPindex,RwPindex)] <- IA_to_R
    mat2[cbind(AwPindex,RwPindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(HwPindex,RwPindex)] <- (1-vax2.ihfr.P[,2])*H_to_RD
    mat2[cbind(HwPindex,DwPindex)] <- vax2.ihfr.P[,2]*H_to_RD
    ########################ASTRAZENECA##########################
    #####################primeira dose###########################
    ######################saindo de S############################
    mat2[cbind(SvAindex,EvAindex)] <- vax1.beta.A*S_to_E  
    ######################saindo de E############################
    mat2[cbind(EvAindex,IvAindex)] <- (1-vax1.ihr.A[,2])*(1-vax1.asymp.A[,2])*E_to_IAH
    mat2[cbind(EvAindex,AvAindex)] <- (1-vax1.ihr.A[,2])*(vax1.asymp.A[,2])*E_to_IAH
    mat2[cbind(EvAindex,HvAindex)] <- vax1.ihr.A[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(IvAindex,RvAindex)] <- IA_to_R
    mat2[cbind(AvAindex,RvAindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(HvAindex,RvAindex)] <- (1-vax1.ihfr.A[,2])*H_to_RD
    mat2[cbind(HvAindex,DvAindex)] <- vax1.ihfr.A[,2]*H_to_RD
    
    #####################segunda dose###########################
    ######################saindo de S############################
    mat2[cbind(SwAindex,EwAindex)] <- vax2.beta.A*S_to_E  
    ######################saindo de E############################
    mat2[cbind(EwAindex,IwAindex)] <- (1-vax2.ihr.A[,2])*(1-vax2.asymp.A[,2])*E_to_IAH
    mat2[cbind(EwAindex,AwAindex)] <- (1-vax2.ihr.A[,2])*vax2.asymp.A[,2]*E_to_IAH
    mat2[cbind(EwAindex,HwAindex)] <- vax2.ihr.A[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(IwAindex,RwAindex)] <- IA_to_R
    mat2[cbind(AwAindex,RwAindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(HwAindex,RwAindex)] <- (1-vax2.ihfr.A[,2])*H_to_RD
    mat2[cbind(HwAindex,DwAindex)] <- vax2.ihfr.A[,2]*H_to_RD
    
    ########################CORONAVAC############################
    #####################primeira dose###########################
    ######################saindo de S############################
    mat2[cbind(SvCindex,EvCindex)] <- vax1.beta.C*S_to_E  
    ######################saindo de E############################
    mat2[cbind(EvCindex,IvCindex)] <- (1-vax1.ihr.C[,2])*(1-vax1.asymp.C[,2])*E_to_IAH
    mat2[cbind(EvCindex,AvCindex)] <- (1-vax1.ihr.C[,2])*(vax1.asymp.C[,2])*E_to_IAH
    mat2[cbind(EvCindex,HvCindex)] <- vax1.ihr.C[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(IvCindex,RvCindex)] <- IA_to_R
    mat2[cbind(AvCindex,RvCindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(HvCindex,RvCindex)] <- (1-vax1.ihfr.C[,2])*H_to_RD
    mat2[cbind(HvCindex,DvCindex)] <- vax1.ihfr.C[,2]*H_to_RD
    
    #####################segunda dose###########################
    ######################saindo de S############################
    mat2[cbind(SwCindex,EwCindex)] <- vax2.beta.C*S_to_E 
    ######################saindo de E############################
    mat2[cbind(EwCindex,IwCindex)] <- (1-vax2.ihr.C[,2])*(1-vax2.asymp.C[,2])*E_to_IAH
    mat2[cbind(EwCindex,AwCindex)] <- (1-vax2.ihr.C[,2])*vax2.asymp.C[,2]*E_to_IAH
    mat2[cbind(EwCindex,HwCindex)] <- vax2.ihr.C[,2]*E_to_IAH
    #####################saindo de I e A#########################
    mat2[cbind(IwCindex,RwCindex)] <- IA_to_R
    mat2[cbind(AwCindex,RwCindex)] <- IA_to_R
    ####################saindo de H##############################
    mat2[cbind(HwCindex,RwCindex)] <- (1-vax2.ihfr.C[,2])*H_to_RD
    mat2[cbind(HwCindex,DwCindex)] <- vax2.ihfr.C[,2]*H_to_RD
    
    
    #############################################################
    mat2[cbind(Sindex,Sindex)] <- 1 -  rowSums(mat2[Sindex,])
    mat2[cbind(Eindex,Eindex)] <- 1 -  rowSums(mat2[Eindex,])
    mat2[cbind(Aindex,Aindex)] <- 1 -  rowSums(mat2[Aindex,])
    mat2[cbind(Iindex,Iindex)] <- 1 -  rowSums(mat2[Iindex,])
    mat2[cbind(Hindex,Hindex)] <- 1 -  rowSums(mat2[Hindex,])
    mat2[cbind(Rindex,Rindex)] <- 1 -  rowSums(mat2[Rindex,])
    mat2[cbind(Dindex,Dindex)] <- 1 -  rowSums(mat2[Dindex,])
    ##########################PFIZER############################
    mat2[cbind(SvPindex,SvPindex)] <- 1 -  rowSums(mat2[SvPindex,])
    mat2[cbind(EvPindex,EvPindex)] <- 1 -  rowSums(mat2[EvPindex,])
    mat2[cbind(AvPindex,AvPindex)] <- 1 -  rowSums(mat2[AvPindex,])
    mat2[cbind(IvPindex,IvPindex)] <- 1 -  rowSums(mat2[IvPindex,])
    mat2[cbind(HvPindex,HvPindex)] <- 1 -  rowSums(mat2[HvPindex,])
    mat2[cbind(RvPindex,RvPindex)] <- 1 -  rowSums(mat2[RvPindex,])
    mat2[cbind(DvPindex,DvPindex)] <- 1 -  rowSums(mat2[DvPindex,])
    mat2[cbind(SwPindex,SwPindex)] <- 1 -  rowSums(mat2[SwPindex,])
    mat2[cbind(EwPindex,EwPindex)] <- 1 -  rowSums(mat2[EwPindex,])
    mat2[cbind(AwPindex,AwPindex)] <- 1 -  rowSums(mat2[AwPindex,])
    mat2[cbind(IwPindex,IwPindex)] <- 1 -  rowSums(mat2[IwPindex,])
    mat2[cbind(HwPindex,HwPindex)] <- 1 -  rowSums(mat2[HwPindex,])
    mat2[cbind(RwPindex,RwPindex)] <- 1 -  rowSums(mat2[RwPindex,])
    mat2[cbind(DwPindex,DwPindex)] <- 1 -  rowSums(mat2[DwPindex,])
    ########################ASTRAZENECA#########################
    mat2[cbind(SvAindex,SvAindex)] <- 1 -  rowSums(mat2[SvAindex,])
    mat2[cbind(EvAindex,EvAindex)] <- 1 -  rowSums(mat2[EvAindex,])
    mat2[cbind(AvAindex,AvAindex)] <- 1 -  rowSums(mat2[AvAindex,])
    mat2[cbind(IvAindex,IvAindex)] <- 1 -  rowSums(mat2[IvAindex,])
    mat2[cbind(HvAindex,HvAindex)] <- 1 -  rowSums(mat2[HvAindex,])
    mat2[cbind(RvAindex,RvAindex)] <- 1 -  rowSums(mat2[RvAindex,])
    mat2[cbind(DvAindex,DvAindex)] <- 1 -  rowSums(mat2[DvAindex,])
    mat2[cbind(SwAindex,SwAindex)] <- 1 -  rowSums(mat2[SwAindex,])
    mat2[cbind(EwAindex,EwAindex)] <- 1 -  rowSums(mat2[EwAindex,])
    mat2[cbind(AwAindex,AwAindex)] <- 1 -  rowSums(mat2[AwAindex,])
    mat2[cbind(IwAindex,IwAindex)] <- 1 -  rowSums(mat2[IwAindex,])
    mat2[cbind(HwAindex,HwAindex)] <- 1 -  rowSums(mat2[HwAindex,])
    mat2[cbind(RwAindex,RwAindex)] <- 1 -  rowSums(mat2[RwAindex,])
    mat2[cbind(DwAindex,DwAindex)] <- 1 -  rowSums(mat2[DwAindex,])
    #######################CORONAVAC############################
    mat2[cbind(SvCindex,SvCindex)] <- 1 -  rowSums(mat2[SvCindex,])
    mat2[cbind(EvCindex,EvCindex)] <- 1 -  rowSums(mat2[EvCindex,])
    mat2[cbind(AvCindex,AvCindex)] <- 1 -  rowSums(mat2[AvCindex,])
    mat2[cbind(IvCindex,IvCindex)] <- 1 -  rowSums(mat2[IvCindex,])
    mat2[cbind(HvCindex,HvCindex)] <- 1 -  rowSums(mat2[HvCindex,])
    mat2[cbind(RvCindex,RvCindex)] <- 1 -  rowSums(mat2[RvCindex,])
    mat2[cbind(DvCindex,DvCindex)] <- 1 -  rowSums(mat2[DvCindex,])
    mat2[cbind(SwCindex,SwCindex)] <- 1 -  rowSums(mat2[SwCindex,])
    mat2[cbind(EwCindex,EwCindex)] <- 1 -  rowSums(mat2[EwCindex,])
    mat2[cbind(AwCindex,AwCindex)] <- 1 -  rowSums(mat2[AwCindex,])
    mat2[cbind(IwCindex,IwCindex)] <- 1 -  rowSums(mat2[IwCindex,])
    mat2[cbind(HwCindex,HwCindex)] <- 1 -  rowSums(mat2[HwCindex,])
    mat2[cbind(RwCindex,RwCindex)] <- 1 -  rowSums(mat2[RwCindex,])
    mat2[cbind(DwCindex,DwCindex)] <- 1 -  rowSums(mat2[DwCindex,])
    mat2
  })$mat2
}

## Rodar simulação
run_markov <- function(mat,nmax,parameters){
  Yn <- parameters$init.condition
  #####prepara linhas
  solution <- matrix(0,nrow = nmax+1,ncol = length(Y))
  vac_in_time.P <- matrix(0,nrow = nmax+1,ncol = 2*parameters$age.bins)
  vac_in_time.A <- matrix(0,nrow = nmax+1,ncol = 2*parameters$age.bins)
  vac_in_time.C <- matrix(0,nrow = nmax+1,ncol = 2*parameters$age.bins)
  
  solution[1,] <- Yn
  
  open.cov <- NULL
  solution <- within(c(parameters,solution),{
    transfer <- rep(0,length(Yn))
    for(i in 2:(nmax+1)){
      ###########
      Yn <- crossprod(mat,Yn)
      if(vax.cov$vax.cov == TRUE){
        for(age in age.bins:2){###### checa toda iteração, mais simples, faz for sobre as idades
          index <- c(SvPindex[age],EvPindex[age],IvPindex[age],AvPindex[age],HvPindex[age],RvPindex[age],DvPindex[age],
                     SwPindex[age],EwPindex[age],IwPindex[age],AwPindex[age],HwPindex[age],RwPindex[age],DwPindex[age],
                     SvAindex[age],EvAindex[age],IvAindex[age],AvAindex[age],HvAindex[age],RvAindex[age],DvAindex[age],
                     SwAindex[age],EwAindex[age],IwAindex[age],AwAindex[age],HwAindex[age],RwAindex[age],DwAindex[age],
                     SvCindex[age],EvCindex[age],IvCindex[age],AvCindex[age],HvCindex[age],RvCindex[age],DvCindex[age],
                     SwCindex[age],EwCindex[age],IwCindex[age],AwCindex[age],HwCindex[age],RwCindex[age],DwCindex[age])
          index.total <- c(index,Sindex[age],Eindex[age],Iindex[age],Aindex[age],Hindex[age],Rindex[age],Dindex[age])
          age.coverage<- sum(Yn[index])/sum(init.condition[index.total])
          if(age.coverage >= vax.cov$cov){
            open.cov <= age - 1
          }
          else{
            open.cov <- age
            break
          }
        }
      }
      ####### taxa de vacinacao primeira dose
      v1<- VAX.DISTR.RATE(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i],Yn[Sindex]+Yn[Rindex],open = open.cov)
      v1.P <- v1*vac.rate.P[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.A <- v1*vac.rate.A[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.C <- v1*vac.rate.C[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.S.P <- v1.P*Yn[Sindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.R.P <- v1.P*Yn[Rindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.S.A <- v1.A*Yn[Sindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.R.A <- v1.A*Yn[Rindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.S.C <- v1.C*Yn[Sindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.R.C <- v1.C*Yn[Rindex]/(Yn[Rindex]+Yn[Sindex]+1)
      vac_in_time.P[i,] <- c(v1.S.P,v1.R.P)
      vac_in_time.A[i,] <- c(v1.S.A,v1.R.A)
      vac_in_time.C[i,] <- c(v1.S.C,v1.R.C)
      ####### taxa de vacinacao segunda dose para Pfizer
      if(i > vax.window.days.P){
        v2.S.P <- (1-desist)*vac_in_time.P[i-vax.window.days.P,1:age.bins]
        v2.R.P <- (1-desist)*vac_in_time.P[i-vax.window.days.P,(age.bins+1):(2*age.bins)]
        ###
        v2.S.hist.P <- (1-desist)*history.P[i,1:age.bins]
        v2.R.hist.P <- (1-desist)*history.P[i,(age.bins+1):(2*age.bins)]
        ###
        v2.S2.P <- ((1-vax1.beta.P*p)**vax.window.days.P)*v2.S.P + ((1-vax1.beta.P*p)**(i-1))*v2.S.hist.P
        v2.R2.P <- v2.R.P+v2.R.hist.P+
          (1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**vax.window.days.P)*v2.S.P +
          (1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**(i-1))*v2.S.hist.P
        ###
        v2.S.P <- pmax(pmin(v2.S2.P,Yn[SvPindex]),0)
        v2.R.P <- pmax(pmin(v2.R2.P,Yn[RvPindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.P <- (1-desist)*history.P[i,1:age.bins]
        v2.R.P <- (1-desist)*history.P[i,(age.bins+1):(2*age.bins)]
        v2.S2.P <- pmax(pmin(((1-vax1.beta.P*p)**(i-1))*v2.S.P,Yn[SvPindex]),0)
        v2.R.P <- pmax(pmin(v2.R.P+(1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**(i-1))*v2.S.P,Yn[RvPindex]  ),0)
        v2.S.P <- v2.S2.P
      }
      if(i > vax.window.days.A){ #vacinacao segunda dose Astrazeneca
        v2.S.A <- (1-desist)*vac_in_time.A[i-vax.window.days.A,1:age.bins]
        v2.R.A <- (1-desist)*vac_in_time.A[i-vax.window.days.A,(age.bins+1):(2*age.bins)]
        ###
        v2.S.hist.A <- (1-desist)*history.A[i,1:age.bins]
        v2.R.hist.A <- (1-desist)*history.A[i,(age.bins+1):(2*age.bins)]
        ###
        v2.S2.A <- ((1-vax1.beta.A*p)**vax.window.days.A)*v2.S.A + ((1-vax1.beta.A*p)**(i-1))*v2.S.hist.A
        v2.R2.A <- v2.R.A+v2.R.hist.A+
          (1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**vax.window.days.A)*v2.S.A +
          (1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**(i-1))*v2.S.hist.A
        ###
        v2.S.A <- pmax(pmin(v2.S2.A,Yn[SvAindex]),0)
        v2.R.A <- pmax(pmin(v2.R2.A,Yn[RvAindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.A <- (1-desist)*history.A[i,1:age.bins]
        v2.R.A <- (1-desist)*history.A[i,(age.bins+1):(2*age.bins)]
        v2.S2.A <- pmax(pmin(((1-vax1.beta.A*p)**(i-1))*v2.S.A,Yn[SvAindex]),0)
        v2.R.A <- pmax(pmin(v2.R.A+(1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**(i-1))*v2.S.A,Yn[RvAindex]  ),0)
        v2.S.A <- v2.S2.A
      }
      if(i > vax.window.days.C){
        v2.S.C <- (1-desist)*vac_in_time.C[i-vax.window.days.C,1:age.bins]
        v2.R.C <- (1-desist)*vac_in_time.C[i-vax.window.days.C,(age.bins+1):(2*age.bins)]
        ###
        v2.S.hist.C <- (1-desist)*history.C[i,1:age.bins]
        v2.R.hist.C <- (1-desist)*history.C[i,(age.bins+1):(2*age.bins)]
        ###
        v2.S2.C <- ((1-vax1.beta.C*p)**vax.window.days.C)*v2.S.C + ((1-vax1.beta.C*p)**(i-1))*v2.S.hist.C
        v2.R2.C <- v2.R.C+v2.R.hist.C+
          (1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**vax.window.days.C)*v2.S.C +
          (1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**(i-1))*v2.S.hist.C
        ###
        v2.S.C <- pmax(pmin(v2.S2.C,Yn[SvCindex]),0)
        v2.R.C <- pmax(pmin(v2.R2.C,Yn[RvCindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.C <- (1-desist)*history.C[i,1:age.bins]
        v2.R.C <- (1-desist)*history.C[i,(age.bins+1):(2*age.bins)]
        v2.S2.C <- pmax(pmin(((1-vax1.beta.C*p)**(i-1))*v2.S.C,Yn[SvCindex]),0)
        v2.R.C <- pmax(pmin(v2.R.C+(1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**(i-1))*v2.S.C,Yn[RvCindex]  ),0)
        v2.S.C <- v2.S2.C
      }
      transfer[Sindex] <- -v1.S.P - v1.S.A - v1.S.C
      transfer[SvPindex] <- v1.S.P - v2.S.P
      transfer[SvAindex] <- v1.S.A - v2.S.A
      transfer[SvCindex] <- v1.S.C - v2.S.C
      transfer[Rindex] <- -v1.R.P - v1.R.A - v1.R.C
      transfer[RvPindex] <- v1.R.P-v2.R.P
      transfer[RvAindex] <- v1.R.A-v2.R.A
      transfer[RvCindex] <- v1.R.C-v2.R.C
      transfer[SwPindex] <- v2.S.P
      transfer[SwAindex] <- v2.S.A
      transfer[SwCindex] <- v2.S.C
      transfer[RwPindex] <- v2.R.P
      transfer[RwAindex] <- v2.R.A
      transfer[RwCindex] <- v2.R.C
      Yn <- Yn + transfer
      solution[i,] <- as.numeric(Yn)
    }
    solution <- cbind(time = 0:nmax,solution)
  })$solution
  colnames(solution) <- c("time",colnames(mat))
  solution <- as.data.frame(solution)
}

######apenas para solução densa AJEITAR PRA AGREGAR IDADE
plot_solution<- function(res,states = NULL,group.vac = FALSE){
  st <- c(c("S", "E", "A", "I", "H", "R", "D"),
          paste0(c("S", "E", "A", "I", "H", "R", "D"),"v"),
          paste0(c("S", "E", "A", "I", "H", "R", "D"),"w"))
  if(group.vac){
    df <- data.frame(time = res[,1],
                     S = rowSums(res[,Sindex+1]),
                     E = rowSums(res[,Eindex+1]),
                     A = rowSums(res[,Aindex+1]),
                     I = rowSums(res[,Iindex+1]),
                     H = rowSums(res[,Hindex+1]),
                     R = rowSums(res[,Rindex+1]),
                     D = rowSums(res[,Dindex+1]),
                     Sv = rowSums(res[,c(SvPindex,SvAindex,SvCindex)+1]),
                     Ev = rowSums(res[,c(EvPindex,EvAindex,EvCindex)+1]),
                     Av = rowSums(res[,c(AvPindex,AvAindex,AvCindex)+1]),
                     Iv = rowSums(res[,c(IvPindex,IvAindex,IvCindex)+1]),
                     Hv = rowSums(res[,c(HvPindex,HvAindex,HvCindex)+1]),
                     Rv = rowSums(res[,c(RvPindex,RvAindex,RvCindex)+1]),
                     Dv = rowSums(res[,c(DvPindex,DvAindex,DvCindex)+1]),
                     Sw = rowSums(res[,c(SwPindex,SwAindex,SwCindex)+1]),
                     Ew = rowSums(res[,c(EwPindex,EwAindex,EwCindex)+1]),
                     Aw = rowSums(res[,c(AwPindex,AwAindex,AwCindex)+1]),
                     Iw = rowSums(res[,c(IwPindex,IwAindex,IwCindex)+1]),
                     Hw = rowSums(res[,c(HwPindex,HwAindex,HwCindex)+1]),
                     Rw = rowSums(res[,c(RwPindex,RwAindex,RwCindex)+1]),
                     Dw = rowSums(res[,c(DwPindex,DwAindex,DwCindex)+1])
    )
  }
  else{
    df <- data.frame(time = res[,1],
                     S = rowSums(res[,Sindex+1]),
                     E = rowSums(res[,Eindex+1]),
                     A = rowSums(res[,Aindex+1]),
                     I = rowSums(res[,Iindex+1]),
                     H = rowSums(res[,Hindex+1]),
                     R = rowSums(res[,Rindex+1]),
                     D = rowSums(res[,Dindex+1]),
                     SvP = rowSums(res[,SvPindex+1]),
                     EvP = rowSums(res[,EvPindex+1]),
                     AvP = rowSums(res[,AvPindex+1]),
                     IvP = rowSums(res[,IvPindex+1]),
                     HvP = rowSums(res[,HvPindex+1]),
                     RvP = rowSums(res[,RvPindex+1]),
                     DvP = rowSums(res[,DvPindex+1]),
                     SwP = rowSums(res[,SwPindex+1]),
                     EwP = rowSums(res[,EwPindex+1]),
                     AwP = rowSums(res[,AwPindex+1]),
                     IwP = rowSums(res[,IwPindex+1]),
                     HwP = rowSums(res[,HwPindex+1]),
                     RwP = rowSums(res[,RwPindex+1]),
                     DwP = rowSums(res[,DwPindex+1]),
                     SvA = rowSums(res[,SvAindex+1]),
                     EvA = rowSums(res[,EvAindex+1]),
                     AvA = rowSums(res[,AvAindex+1]),
                     IvA = rowSums(res[,IvAindex+1]),
                     HvA = rowSums(res[,HvAindex+1]),
                     RvA = rowSums(res[,RvAindex+1]),
                     DvA = rowSums(res[,DvAindex+1]),
                     SwA = rowSums(res[,SwAindex+1]),
                     EwA = rowSums(res[,EwAindex+1]),
                     AwA = rowSums(res[,AwAindex+1]),
                     IwA = rowSums(res[,IwAindex+1]),
                     HwA = rowSums(res[,HwAindex+1]),
                     RwA = rowSums(res[,RwAindex+1]),
                     DwA = rowSums(res[,DwAindex+1]),
                     SvC = rowSums(res[,SvCindex+1]),
                     EvC = rowSums(res[,EvCindex+1]),
                     AvC = rowSums(res[,AvCindex+1]),
                     IvC = rowSums(res[,IvCindex+1]),
                     HvC = rowSums(res[,HvCindex+1]),
                     RvC = rowSums(res[,RvCindex+1]),
                     DvC = rowSums(res[,DvCindex+1]),
                     SwC = rowSums(res[,SwCindex+1]),
                     EwC = rowSums(res[,EwCindex+1]),
                     AwC = rowSums(res[,AwCindex+1]),
                     IwC = rowSums(res[,IwCindex+1]),
                     HwC = rowSums(res[,HwCindex+1]),
                     RwC = rowSums(res[,RwCindex+1]),
                     DwC = rowSums(res[,DwCindex+1])
    )
  }
  
  if(!is.null(states)){
    df <- df[,c("time",states)]
  }
  df <- df %>% gather("state", "value", -time)
  g <- ggplot(df, aes(x = time, y = value, col = state)) +
    geom_line(size = 1) + theme_minimal()
  g
}

pivot_age_classes <- function(res) {
  epiclasses <- c("S", "E", "A", "I", "H", "R", "D")
  vacclasses <- c("", paste0(rep(c("v", "w"), each=3), c("P", "A", "C")))
  all.res <- NULL
  for (epi in epiclasses) {
    for (vac in vacclasses) {
      classe <- paste0(epi, vac)
      res.classe <- res %>%
        select(time, matches(paste0("^", classe, "\\d+"))) %>%
        pivot_longer(
                     cols = !time,
                     names_to = "age_class",
                     names_prefix = classe,
                     names_transform = list(age_class = as.integer),
                     values_to = classe
        ) %>%
        as.data.frame()
      if (is.null(all.res)) {
        all.res <- res.classe
      } else {
        all.res <- cbind(all.res, res.classe[,3])
        colnames(all.res)[ncol(all.res)] <- classe
      }
    }
  }
  return(all.res)
}

dose_class_helper <- function(x) {
  ifelse(x == "", 0,
         ifelse(x=="v", 1,
                ifelse(x=="w", 2, NA
                )
         )
  )
}

pivot_vaccines <- function(res) {
  epiclasses <- c("S", "E", "A", "I", "H", "R", "D")
  all.res <- NULL
  for (epi in epiclasses) {
    res.epi <- res %>%
      select(time, age_class, starts_with(epi, ignore.case = FALSE)) %>%
      pivot_longer(
                   cols = starts_with(epi, ignore.case = FALSE),
                   names_to = "vaccine",
                   names_prefix = epi,
                   values_to = epi
      ) %>%
      as.data.frame()
    if (is.null(all.res)) {
      all.res <- res.epi
    } else {
      all.res <- cbind(all.res, res.epi[,epi])
      colnames(all.res)[ncol(all.res)] <- epi
    }
  }
  all.res$dose <- dose_class_helper(substring(all.res$vaccine, 1, 1))
  all.res$vaccine <- substring(all.res$vaccine, 2, 2)
  return(all.res)
}

pivot_all <- function(res)
  return(pivot_vaccines(pivot_age_classes(res)))

group_age_classes <- function(res){
  epiclasses = c(c("S", "E", "A", "I", "H", "R", "D"),
             paste0(c("S", "E", "A", "I", "H", "R", "D"),"vA"),
             paste0(c("S", "E", "A", "I", "H", "R", "D"),"vP"),
             paste0(c("S", "E", "A", "I", "H", "R", "D"),"vC"),
             paste0(c("S", "E", "A", "I", "H", "R", "D"),"wA"),
             paste0(c("S", "E", "A", "I", "H", "R", "D"),"wP"),
             paste0(c("S", "E", "A", "I", "H", "R", "D"),"wC"))
  all.res <- data.frame(matrix(0,nrow = nrow(res),ncol = length(epiclasses)+1))
  names(all.res) <- c("time",epiclasses)
  all.res[,"time"] <- res[,"time"]
  for (classe in epiclasses) {
      classes <- paste0(classe,1:10)
      all.res[,as.character(classe)] <- res %>% select(classes) %>% rowSums()
  }
  return(all.res)
}

group_vaccines <- function(res){
  res <- pivot_all(res)
  epiclasses = c("S", "E", "A", "I", "H", "R", "D")

  # all.res <- data.frame(matrix(0,nrow = nrow(res),ncol = length(epiclasses)+1))
  # names(all.res) <- c("time",epiclasses)
  all.res <- c()
  for (classe in epiclasses) {
    res2 <- res %>% select(time, age_class,vaccine,as.character(classe),dose) %>%
      group_by(time,as.character(classe),dose) %>% summarize(n = sum())
    all.res <- rbind(all.res,as.data.frame(res2))
  }
  names(all.res) <- c("time","classe","dose","n")
  return(all.res)
}

group_vaccines.2 <- function(res){
  res <- group_age_classes(res)
  no_vac <- c("S", "E", "A", "I", "H", "R", "D")
  epiclasses = c(paste0(c("S", "E", "A", "I", "H", "R", "D"),"v"),
                 paste0(c("S", "E", "A", "I", "H", "R", "D"),"w"))
  vac <- c("A","P","C")
  all.res <- data.frame(matrix(0,nrow = nrow(res),ncol = 3*length(no_vac)+1))
  names(all.res) <- c("time",no_vac,epiclasses)
  all.res[,no_vac] <- res[,no_vac]
  all.res[,"time"] <- res[,"time"]
  for (classe in epiclasses) {
    classes <- paste0(classe,vac)
    all.res[,as.character(classe)] <- res %>% select(classes) %>% rowSums()
  }
  return(all.res)
}

