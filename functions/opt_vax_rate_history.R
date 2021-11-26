library(lpSolve)

#' opt_vax_rate.
#' 
#' Optimal vaccination rate function, 
#' its optimaze the problem of rollout vaccine solving a by non-linear programming solver.
#'
#' @param VAX.INITIAL.STORAGE.NUM # Number of READY vaccines on day one 'VAX.INITIAL.STORAGE.NUM'.
#' @param VAX.PRODUCTION.RATE # Number of vaccines produced per day 'VAX.PRODUCTION.RATE'.
#' @param MAX.VAC.RATE # Max number of vaccine applications per day 'MAX.VAC.RATE'.
#' @param VAX.WINDOW.DAYS # Time window between first and second vaccines 'VAX.WINDOWS.DAYS'. 
#' Its varies with the vaccine type.
#' @param SECOND.VAX.LOSS.FRAC # Fraction of people that take the first but not the second dose 'SECODNE.VAX.LOSS.FRAC'.
#' @param MAX.TIME.DAYS # time until end of vaccination schedule, in days 'MAX.TIME.DAYS'.
#' @param VAX.HISTORY #number of people with time i between doses, in days 'VAX.HISTORY'.
#' @import lpSolve
#'
#' @return
#' @export 
#' @examples
opt_vax_rate <- function(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                         MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                         MAX.TIME.DAYS,VAX.HISTORY = NULL) {
  times <- seq(0, MAX.TIME.DAYS)
  
  
  #####transforma vax.production.rate em vetor
  ####check vax.production.rate size is == max.time.days, if not, append with zeros
  # VAX.PRODUCTION.RATE <- VAX.PRODUCTION.RATE
  if (length(VAX.PRODUCTION.RATE) == 1){
    VAX.PRODUCTION.RATE <- c(0,rep(VAX.PRODUCTION.RATE,MAX.TIME.DAYS))
  }
  if(length(VAX.PRODUCTION.RATE) < length(times)){
    VAX.PRODUCTION.RATE <- c(0,VAX.PRODUCTION.RATE,rep(0,length(times)-length(VAX.PRODUCTION.RATE)+VAX.WINDOW.DAYS))
  }
  else if(length(VAX.PRODUCTION.RATE) >= length(times) & length(VAX.PRODUCTION.RATE) < length(times)+VAX.WINDOW.DAYS ){
    VAX.PRODUCTION.RATE <- c(0,VAX.PRODUCTION.RATE,rep(VAX.PRODUCTION.RATE[length(VAX.PRODUCTION.RATE)],length(times)+VAX.WINDOW.DAYS-length(VAX.PRODUCTION.RATE)))
  }
  
  if(!is.null(VAX.HISTORY)){
    if(length(VAX.HISTORY) > VAX.WINDOW.DAYS){ 
      #se tem pessoas com dose atrasada, acumula tudo no ultimo dia
      VAX.HISTORY[VAX.WINDOW.DAYS] <- sum(VAX.HISTORY[VAX.WINDOW.DAYS:length(VAX.HISTORY)])
      VAX.HISTORY <- VAX.HISTORY[1:VAX.WINDOW.DAYS]
    }
    ###reverte vetor para pessoas com dose recebida a mais tempo tomem antes
    VAX.HISTORY <- rev(VAX.HISTORY)
    VAX.HISTORY <- c(0,VAX.HISTORY)##append zero no começo
    VAX.HISTORY <- c(VAX.HISTORY,rep(0,length(times)-length(VAX.HISTORY)))###append zero no fim
  }
  ####garante que sum_i=1^j g_i <= V_0 + sum_i=1^j p_i
  ####como essa condicao é mais restritiva para feasibility, computamos primeiro
  prod.cum <- VAX.INITIAL.STORAGE.NUM+cumsum(VAX.PRODUCTION.RATE) #V_0 + sum_i=1^j p_i
  HISTORY.2 <- rep(0,length(VAX.HISTORY))
  sobra <- 0
  for(i in 2:length(VAX.HISTORY)){
    val.1 <- sum(HISTORY.2[1:i-1]) ### o que já foi usado
    val.2 <- val.1 + sobra + VAX.HISTORY[i] #### o que quer ser aplicado
    sobra <- max(0, val.2 - prod.cum[i]) ##quanto sobra se satisfizer a restricao
    HISTORY.2[i] <- val.2 - sobra - val.1 ###valor aplicado
  }
  # HISTORY.2 <- HISTORY.2[2:length(HISTORY.2)]
  HISTORY <- rep(0,length(HISTORY.2))
  
  ##### satisfazendo g_i <= v_max
  SURPLUS <- 0
    for (i in 1:length(HISTORY.2)) {
      val <- min(HISTORY.2[i]+SURPLUS,MAX.VAC.RATE)
      SURPLUS <- SURPLUS + HISTORY.2[i]-val
      HISTORY[i] <- val
    }
  # HISTORY <- HISTORY[2:length(HISTORY)]
  


  
  alpha = 1 - SECOND.VAX.LOSS.FRAC

  V.T = max(VAX.WINDOW.DAYS *alpha * MAX.VAC.RATE - sum(VAX.PRODUCTION.RATE[length(times):(length(times)+VAX.WINDOW.DAYS)]), 0)
  # print(all(VAX.INITIAL.STORAGE.NUM + cumsum(VAX.PRODUCTION.RATE-MAX.VAC.RATE)-
  #             (MAX.TIME.DAYS - VAX.WINDOW.DAYS) * alpha * MAX.VAC.RATE >= V.T ))
  # boundary-free case
  #####modificar para funcionar com history
  if(all(VAX.INITIAL.STORAGE.NUM + cumsum(VAX.PRODUCTION.RATE-MAX.VAC.RATE)-
         (MAX.TIME.DAYS - VAX.WINDOW.DAYS) * alpha * MAX.VAC.RATE >= V.T )){
    return((MAX.VAC.RATE - HISTORY)/(1+alpha))
  }
  
  # feasibility is trivial: V(T) > 0 if VAX.RATE == 0
  # solution using LP
  N <- length(times)
  M <- matrix(0, N, N)
  M[row(M)-col(M) >= 0] = 1
  M.2 <- matrix(0, N, N)
  M.2[seq(1 + VAX.WINDOW.DAYS, N),] <- M[seq(1, N - VAX.WINDOW.DAYS),]
  
  M.total <- - M - alpha*M.2
  
  # objective
  C.T <- colSums(M.total)
  ## inequality constraints (Ax <= b)
  A <- M.total
  # diagonal matrix offset by window
  offdiag <- matrix(0, nrow = N, ncol = N)
  offdiag[row(offdiag) - col(offdiag) == VAX.WINDOW.DAYS] = 1
  A <- rbind(- M.total,
             diag(nrow = N, ncol = N) + alpha * offdiag,
             M[N,]
  )
  b <- c(VAX.INITIAL.STORAGE.NUM + cumsum(VAX.PRODUCTION.RATE[1:length(times)])-cumsum(HISTORY),
         rep(MAX.VAC.RATE, N)-HISTORY,
         (VAX.INITIAL.STORAGE.NUM + sum(VAX.PRODUCTION.RATE[1:(length(times)+VAX.WINDOW.DAYS)]))/(1+alpha)
  )
  xopt <-  lp(direction="min",
              objective.in = C.T,
              const.mat = A,
              const.dir = rep("<=", length(b)),
              const.rhs = b)
  #print(paste("Vacinas restantes no tempo final:",
  #            sum(M.total[nrow(M.total),] * xopt$solution) + b[length(times)]))
  return(list(OPT.VAX.RATE = xopt$solution,VAX.PRODUCTION.RATE = VAX.PRODUCTION.RATE,HISTORY = HISTORY))
}

#' plot_vac_schedule.
#' 
#' A function to plot the rollout vaccine after the optimal solution.
#' 
#' @param OPT.VAX.RATE # Optimal vaccine rate from the `opt_vax_rate` 'OPT.VAX.RATE'.
#' @param VAX.INITIAL.STORAGE.NUM # Number of READY vaccines on day one 'VAX.INITIAL.STORAGE.NUM'.
#' @param VAX.PRODUCTION.RATE # Number of vaccines produced per day 'VAX.PRODUCTION.RATE'.
#' @param MAX.VAC.RATE # Max number of vaccine applications per day 'MAX.VAC.RATE'.
#' @param VAX.WINDOW.DAYS # Time window between first and second vaccines 'VAX.WINDOWS.DAYS'.
#' Its varies with the vaccine type.
#' @param SECOND.VAX.LOSS.FRAC # Fraction of people that take the first but not the second dose 'SECODNE.VAX.LOSS.FRAC'.
#' @param MAX.TIME.DAYS # time until end of vaccination schedule, in days 'MAX.TIME.DAYS'.
#' @import lpsolve
#' 
#' @return
#' @export
#' 
#' @examples
# plot_vac_schedule <- function(OPT.VAX.RATE, VAX.INITIAL.STORAGE.NUM,
#                               VAX.PRODUCTION.RATE, MAX.VAC.RATE,
#                               VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
#                               MAX.TIME.DAYS) {
#   times <- seq(1, MAX.TIME.DAYS)
#   alpha = 1 - SECOND.VAX.LOSS.FRAC
#   V <- VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * times -
#     cumsum(OPT.VAX.RATE) -
#     alpha * c(rep(0, VAX.WINDOW.DAYS),
#               cumsum(OPT.VAX.RATE[seq(1, length(times) - VAX.WINDOW.DAYS)]))
#   
#   par(mar = c(5, 4, 4, 4) + 0.3)
#   plot(times, OPT.VAX.RATE,
#        xlab="time (days)",
#        ylab="vaccination rate (doses/day)",
#        type='l',lwd=2.5)
#   lines(times, alpha * c(rep(0, VAX.WINDOW.DAYS),
#                          OPT.VAX.RATE[-seq(length(times)+1-VAX.WINDOW.DAYS, length(times))]),
#         type='l', lty=2,lwd = 2.5)
#   par(new = TRUE)
#   plot(times, V, type = "l", axes = FALSE, bty = "n", lty = 3, xlab = "", ylab = "", col="blue",lwd =2.5)
#   axis(side=4, at = pretty(round(range(V), -5)))
#   mtext("Vaccine doses stored", side=4, line=3, col="blue")
# }

#' interpolate.VAX.RATE.
#' 
#' A function to interpolate vaccine rate rollout give the vaccine rate
#'
#' @param OPT.VAX.RATE # Optimal vaccine rate from the `opt_vax_rate` 'OPT.VAX.RATE'.
#'
#' @return
#' @export
#'
#' @examples
interpolate.VAX.RATE <- function(OPT.VAX.RATE) {
  VAX.RATE <- function(t){
    if (t <= 0)
      return(0)
    tu <- ceiling(t)
    x <- tu - t
    s = 0.1
    if (x > 1-s && tu > 1)
      return((x)*OPT.VAX.RATE[tu] +
               (1-x)*(OPT.VAX.RATE[tu] + OPT.VAX.RATE[tu-1])/2)
    if (x < s)
      return((1-x/s)*OPT.VAX.RATE[tu] +
               (x/s)*(OPT.VAX.RATE[tu] + OPT.VAX.RATE[tu+1])/2)
    return(OPT.VAX.RATE[tu])
  }
}

generate_vac_schedule <- function( VAX.INITIAL.STORAGE.NUM,
                                   VAX.PRODUCTION.RATE, MAX.VAC.RATE,
                                   VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                   MAX.TIME.DAYS,VAX.HISTORY = NULL){
  
  OPT.VAX.RATE <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                               MAX.VAC.RATE, VAX.WINDOW.DAYS,
                               SECOND.VAX.LOSS.FRAC, MAX.TIME.DAYS,VAX.HISTORY)
  # print(OPT.VAX.RATE)
  # print(MAX.TIME.DAYS)
  HISTORY <- OPT.VAX.RATE$HISTORY
  
  OPT.VAX.RATE <- OPT.VAX.RATE$OPT.VAX.RATE
  times <- seq(0, MAX.TIME.DAYS)
  if(length(VAX.PRODUCTION.RATE) == 1){
    v2 <- cumsum(c(0,rep(VAX.PRODUCTION.RATE,MAX.TIME.DAYS)))
  }
  else if(length(VAX.PRODUCTION.RATE) > 1 & length(VAX.PRODUCTION.RATE) < MAX.TIME.DAYS){
    v2 <- cumsum(c(0,VAX.PRODUCTION.RATE,rep(0,length(times)-length(VAX.PRODUCTION.RATE))))
  }
  else{
    v2<- cumsum(c(0,VAX.PRODUCTION.RATE))
  }
  v2 <- v2[1:(MAX.TIME.DAYS+1)]
  # print(c(length(HISTORY),length(OPT.VAX.RATE),length(v2)))
  alpha = 1 - SECOND.VAX.LOSS.FRAC
  V <- VAX.INITIAL.STORAGE.NUM + v2 -
    cumsum(OPT.VAX.RATE) - cumsum(HISTORY) - 
    alpha * c(rep(0, VAX.WINDOW.DAYS),
              cumsum(OPT.VAX.RATE[seq(1, length(times) - VAX.WINDOW.DAYS)]))
  # V <- V[1:(length(V)-1)]
  V[abs(V) < 1e-5] <- 0
  V2 <- c(rep(0,VAX.WINDOW.DAYS),alpha*OPT.VAX.RATE)
  V2 <- V2[1:(MAX.TIME.DAYS+1)]
  df <- data.frame(V,OPT.VAX.RATE,V2,HISTORY)
  colnames(df) <- c("Stock","Vax1","Vax2","History")
  df$t <- as.numeric(row.names(df))
  return(df)
}

plot_vac_schedule.2 <- function (df){
  sc <- max(df$Stock)/max(df$Vax1)+1
  g <- ggplot(df, aes(x = t, y = Vax1/POP.TOTAL.NUM)) +facet_wrap(~Window,nrow = 1)+ geom_line(size = 1) +
    geom_line(aes(y = Vax2/POP.TOTAL.NUM), linetype = "longdash",size = 1,color = "#b1b1b1") +
    geom_line(aes(y = (Stock/sc)/POP.TOTAL.NUM),color = "#466496",linetype = "dotdash",size=1)+
    geom_line(aes(y = History/POP.TOTAL.NUM),color = "red",linetype = "dotdash",size=1)+
    scale_y_continuous("Vaccination rate (doses/pop-day)",sec.axis = sec_axis(~(.*sc),name = "Vaccine doses stored (doses/pop)",labels = scales::percent_format(accuracy=0.1)),labels = scales::percent_format(accuracy = 0.1))+
    scale_x_continuous("t (days)")+
    theme_bw()+
    theme(axis.title.y.right = element_text(color = "#466496"),text = element_text(size = 16))
  g
}

plot_vac_schedule.3 <- function (df){
  sc <- max(df$Stock)/max(df$Vax1)+1
  g <- ggplot(df, aes(x = as.Date.numeric(t,origin = "2021-08-08"), y = Vax1)) +facet_wrap(~Window,nrow = 1)+ geom_line(size = 1) +
    geom_line(aes(y = Vax2), linetype = "longdash",size = 1,color = "#b1b1b1") +
    geom_line(aes(y = (Stock/sc)),color = "#466496",linetype = "dotdash",size=1)+
    geom_line(aes(y = History),color = "red",linetype = "dotdash",size=1)+
    scale_y_continuous("Taxa de vacinação (doses/dia)",sec.axis = sec_axis(~(.*sc),name = "Doses estocadas",labels = scales::comma),labels = scales::comma)+
    labs(x = "Data")+
    theme_bw()+
    theme(axis.title.y.right = element_text(color = "#466496"),text = element_text(size = 16))
  g
}
############## to plot
# df <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
#                             MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
