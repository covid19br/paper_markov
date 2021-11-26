###### taxa vacinacao
VAX.DISTR.RATE <- function(V, POP, tol=100,open = NULL) {
  vac <- rep(0,length(POP))
  V2 <- V
  if(is.null(open)){
    i <- length(POP) ##### vacinando idosos primeiro
    while(i > 0) {
      if(POP[i] > tol & V2 > 0){
        vac[i] <- max(min(V2,POP[i]),0)
        V2 <- V2 - vac[i]
      }
      else{
        vac[i] <- 0
      }
      i <- i - 1
    }
  }
  else { #### para i >= open, passa vacina 
    if(V <= sum(POP[length(POP):open])){
      for(i in length(POP):open){
        vac[i] <- V*POP[i]/sum(POP[length(POP):open])
      }
    }
    else{ 
      vac[length(POP):open] <- POP[length(POP):open] ####vacina todos das outras faixas
      V2 <- V2 - sum(POP[length(POP):open])
      while(open > 0 && V2 > 0){### abre a próxima faixa e torce pra não bugar
        vac[open] <- max(min(V2,POP[open]),0)
        V2 <- V2 - vac[open]
        open <- open - 1
      }
    }
  }
  return(vac)
}