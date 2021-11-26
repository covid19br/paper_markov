library(rlist)
library(LearnBayes)
params <- list()
p <- list(param = "p",distr = "unif",vals = c(0.005,0.02))
prevalence <- list(param = "prevalence", distr = "unif", vals = c(0.40,0.60))
#####pfizer


vax1.beta.effic.P <- list(param = "vax1.beta.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.53),list(p=0.975,x=0.66)))
vax2.beta.effic.P <- list(param = "vax2.beta.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.88),list(p=0.975,x=0.95)))

vax1.symp.effic.P <- list(param = "vax1.symp.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.56),list(p=0.975,x=0.71)))
vax2.symp.effic.P <- list(param = "vax2.symp.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.70),list(p=0.975,x=0.93)))

vax1.hosp.effic.P <- list(param = "vax1.hosp.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.75),list(p=0.975,x=0.88)))
vax2.hosp.effic.P <- list(param = "vax2.hosp.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.82),list(p=0.975,x=0.99)))

vax1.death.effic.P <- list(param = "vax1.death.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.75),list(p=0.975,x=0.88)))
vax2.death.effic.P <- list(param = "vax2.death.effic.P",distr = "beta",vals = beta.select(list(p=0.025,x=0.82),list(p=0.975,x=0.99)))

#####coronavac

### https://doi.org/10.1056/NEJMoa2107715
vax1.symp.effic.C <- list(param = "vax1.symp.effic.C",distr = "beta",vals = beta.select(list(p=0.025,x=0.14),list(p=0.975,x=0.18)))
vax2.symp.effic.C <- list(param = "vax2.symp.effic.C",distr = "beta",vals = beta.select(list(p=0.025,x=0.65),list(p=0.975,x=0.69)))

vax1.hosp.effic.C1 <- list(param = "vax1.hosp.effic.C1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.271),list(p=0.975,x=0.397)))
vax1.hosp.effic.C2 <- list(param = "vax1.hosp.effic.C7",distr = "beta",vals = beta.select(list(p=0.025,x=0.258),list(p=0.975,x=0.330)))
vax1.hosp.effic.C3 <- list(param = "vax1.hosp.effic.C8",distr = "beta",vals = beta.select(list(p=0.025,x=0.299),list(p=0.975,x=0.351)))
vax1.hosp.effic.C4 <- list(param = "vax1.hosp.effic.C9",distr = "beta",vals = beta.select(list(p=0.025,x=0.021),list(p=0.975,x=0.138)))
vax1.hosp.effic.C5 <- list(param = "vax1.hosp.effic.C10",distr = "unif",vals = c(0,0))

vax2.hosp.effic.C1 <- list(param = "vax2.hosp.effic.C1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.813),list(p=0.975,x=0.867)))
vax2.hosp.effic.C2 <- list(param = "vax2.hosp.effic.C7",distr = "beta",vals = beta.select(list(p=0.025,x=0.763),list(p=0.975,x=0.798)))
vax2.hosp.effic.C3 <- list(param = "vax2.hosp.effic.C8",distr = "beta",vals = beta.select(list(p=0.025,x=0.72),list(p=0.975,x=0.76)))
vax2.hosp.effic.C4 <- list(param = "vax2.hosp.effic.C9",distr = "beta",vals = beta.select(list(p=0.025,x=0.599),list(p=0.975,x=0.660)))
vax2.hosp.effic.C5 <- list(param = "vax2.hosp.effic.C10",distr = "beta",vals = beta.select(list(p=0.025,x=0.228),list(p=0.975,x=0.413)))

vax1.death.effic.C1 <- list(param = "vax1.death.effic.C1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.264),list(p=0.975,x=0.539)))
vax1.death.effic.C2 <- list(param = "vax1.death.effic.C7",distr = "beta",vals = beta.select(list(p=0.025,x=0.303),list(p=0.975,x=0.407)))
vax1.death.effic.C3 <- list(param = "vax1.death.effic.C8",distr = "beta",vals = beta.select(list(p=0.025,x=0.347),list(p=0.975,x=0.415)))
vax1.death.effic.C4 <- list(param = "vax1.death.effic.C9",distr = "beta",vals = beta.select(list(p=0.025,x=0.027),list(p=0.975,x=0.107)))
vax1.death.effic.C5 <- list(param = "vax1.death.effic.C10",distr = "unif",vals = c(0,0))

vax2.death.effic.C1 <- list(param = "vax2.death.effic.C1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.514),list(p=0.975,x=0.821)))
vax2.death.effic.C2 <- list(param = "vax2.death.effic.C7",distr = "beta",vals = beta.select(list(p=0.025,x=0.763),list(p=0.669,x=0.833)))
vax2.death.effic.C3 <- list(param = "vax2.death.effic.C8",distr = "beta",vals = beta.select(list(p=0.025,x=0.766),list(p=0.975,x=0.800)))
vax2.death.effic.C4 <- list(param = "vax2.death.effic.C9",distr = "beta",vals = beta.select(list(p=0.025,x=0.636),list(p=0.975,x=0.706)))
vax2.death.effic.C5 <- list(param = "vax2.death.effic.C10",distr = "beta",vals = beta.select(list(p=0.025,x=0.238),list(p=0.975,x=0.451)))

### https://doi.org/10.1101/2021.05.19.21257472 
# vax1.symp.effic.C2 <- list(param = "vax1.symp.effic.C8-10",distr = "betaneg",vals = beta.select(list(p=0.025,x=0.478),list(p=0.975,x=0.6165)))
# vax2.symp.effic.C2 <- list(param = "vax2.symp.effic.C8-10",distr = "beta",vals = beta.select(list(p=0.025,x=0.269),list(p=0.975,x=0.533)))
# 
# vax1.hosp.effic.C2 <- list(param = "vax1.hosp.effic.C8-10",distr = "betaneg",vals = beta.select(list(p=0.025,x=0.495),list(p=0.975,x=0.671)))
# vax2.hosp.effic.C2 <- list(param = "vax2.hosp.effic.C8-10",distr = "beta",vals = beta.select(list(p=0.025,x=0.442),list(p=0.975,x=0.698)))
# 
# vax1.death.effic.C2 <- list(param = "vax1.death.effic.C8-10",distr = "betaneg",vals = beta.select(list(p=0.025,x=0.5355),list(p=0.975,x=0.7485)))
# vax2.death.effic.C2 <- list(param = "vax2.death.effic.C8-10",distr = "beta",vals = beta.select(list(p=0.025,x=0.537),list(p=0.975,x=0.823)))

######
######astrazeneca

### https://doi.org/10.1016/S0140-6736(21)00432-3
vax1.beta.effic.A <- list(param = "vax1.beta.effic.A",distr = "beta",vals = beta.select(list(p=0.025,x=0.460),list(p=0.975,x=0.759)))
vax2.beta.effic.A <- list(param = "vax2.beta.effic.A",distr = "beta",vals = beta.select(list(p=0.025,x=0.358),list(p=0.975,x=0.750)))

vax1.symp.effic.A1 <- list(param = "vax1.symp.effic.A1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.27),list(p=0.975,x=0.66)))
vax2.symp.effic.A1 <- list(param = "vax2.symp.effic.A1-6",distr = "unif",vals = c(0.099,0.99))
vax1.symp.effic.A2 <- list(param = "vax1.symp.effic.A7-10",distr = "beta",vals = beta.select(list(p=0.025,x=0.264),list(p=0.975,x=0.397)))
vax2.symp.effic.A2 <- list(param = "vax2.symp.effic.A7-10",distr = "beta",vals = beta.select(list(p=0.025,x=0.692),list(p=0.975,x=0.842)))

vax1.hosp.effic.A1 <- list(param = "vax1.hosp.effic.A1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.62),list(p=0.975,x=0.66)))

vax1.hosp.effic.A2 <- list(param = "vax1.hosp.effic.A7",distr = "beta",vals = beta.select(list(p=0.025,x=0.424),list(p=0.975,x=0.474)))
vax1.hosp.effic.A3 <- list(param = "vax1.hosp.effic.A8",distr = "beta",vals = beta.select(list(p=0.025,x=0.252),list(p=0.975,x=0.398)))
vax1.hosp.effic.A4 <- list(param = "vax1.hosp.effic.A9",distr = "beta",vals = beta.select(list(p=0.025,x=0.280),list(p=0.975,x=0.374)))
vax1.hosp.effic.A5 <- list(param = "vax1.hosp.effic.A10",distr = "unif",vals = c(0,0))

vax2.hosp.effic.A1 <- list(param = "vax2.hosp.effic.A1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.898),list(p=0.975,x=0.966)))
vax2.hosp.effic.A2 <- list(param = "vax2.hosp.effic.A7",distr = "beta",vals = beta.select(list(p=0.025,x=0.843),list(p=0.975,x=0.956)))
vax2.hosp.effic.A3 <- list(param = "vax2.hosp.effic.A8",distr = "beta",vals = beta.select(list(p=0.025,x=0.846),list(p=0.975,x=0.912)))
vax2.hosp.effic.A4 <- list(param = "vax2.hosp.effic.A9",distr = "beta",vals = beta.select(list(p=0.025,x=0.849),list(p=0.975,x=0.887)))
vax2.hosp.effic.A5 <- list(param = "vax2.hosp.effic.A10",distr = "beta",vals = beta.select(list(p=0.025,x=0.354),list(p=0.975,x=0.685)))

vax1.death.effic.A1 <- list(param = "vax1.death.effic.A1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.618),list(p=0.975,x=0.676)))
vax1.death.effic.A2 <- list(param = "vax1.death.effic.A7",distr = "beta",vals = beta.select(list(p=0.025,x=0.410),list(p=0.975,x=0.494)))
vax1.death.effic.A3 <- list(param = "vax1.death.effic.A8",distr = "beta",vals = beta.select(list(p=0.025,x=0.269),list(p=0.975,x=0.458)))
vax1.death.effic.A4 <- list(param = "vax1.death.effic.A9",distr = "beta",vals = beta.select(list(p=0.025,x=0.322),list(p=0.975,x=0.434)))
vax1.death.effic.A5 <- list(param = "vax1.death.effic.A10",distr = "unif",vals = c(0,0))
vax2.death.effic.A1 <- list(param = "vax2.death.effic.A1-6",distr = "beta",vals = beta.select(list(p=0.025,x=0.721),list(p=0.975,x=0.984)))
vax2.death.effic.A2 <- list(param = "vax2.death.effic.A7",distr = "beta",vals = beta.select(list(p=0.025,x=0.718),list(p=0.669,x=0.962)))
vax2.death.effic.A3 <- list(param = "vax2.death.effic.A8",distr = "beta",vals = beta.select(list(p=0.025,x=0.881),list(p=0.975,x=0.953)))
vax2.death.effic.A4 <- list(param = "vax2.death.effic.A9",distr = "beta",vals = beta.select(list(p=0.025,x=0.891),list(p=0.975,x=0.929)))
vax2.death.effic.A5 <- list(param = "vax2.death.effic.A10",distr = "beta",vals = beta.select(list(p=0.025,x=0.514),list(p=0.975,x=0.821)))


params <- list(prevalence = prevalence,
               ## pfizer D1
               vax1.beta.effic.P = vax1.beta.effic.P,
               vax1.symp.effic.P = vax1.symp.effic.P,
               vax1.hosp.effic.P = vax1.hosp.effic.P,
               vax1.death.effic.P = vax1.death.effic.P,
               ## pfizer D2
               vax2.beta.effic.P = vax2.beta.effic.P,
               vax2.symp.effic.P = vax2.symp.effic.P,
               vax2.hosp.effic.P = vax2.hosp.effic.P,
               vax2.death.effic.P = vax2.death.effic.P,
               ## coronavac symp
               vax1.symp.effic.C = vax1.symp.effic.C,
               vax1.symp.effic.C = vax1.symp.effic.C,
               ## coronavac hosp D1
               "vax1.hosp.effic.C1-6" = vax1.hosp.effic.C1,
               "vax1.hosp.effic.C7" = vax1.hosp.effic.C2,
               "vax1.hosp.effic.C8" = vax1.hosp.effic.C3,
               "vax1.hosp.effic.C9" = vax1.hosp.effic.C4,
               "vax1.hosp.effic.C10" = vax1.hosp.effic.C5,
               ## coronavac hosp D2
               "vax2.hosp.effic.C1-6" = vax2.hosp.effic.C1,
               "vax2.hosp.effic.C7" = vax2.hosp.effic.C2,
               "vax2.hosp.effic.C8" = vax2.hosp.effic.C3,
               "vax2.hosp.effic.C9" = vax2.hosp.effic.C4,
               "vax2.hosp.effic.C10" = vax2.hosp.effic.C5,
               ## coronavac death D1
               "vax1.death.effic.C1-6" = vax1.death.effic.C1,
               "vax1.death.effic.C7" = vax1.death.effic.C2,
               "vax1.death.effic.C8" = vax1.death.effic.C3,
               "vax1.death.effic.C9" = vax1.death.effic.C4,
               "vax1.death.effic.C10" = vax1.death.effic.C5,
               ## coronavac death D2
               "vax2.death.effic.C1-6" = vax2.death.effic.C1,
               "vax2.death.effic.C7" = vax2.death.effic.C2,
               "vax2.death.effic.C8" = vax2.death.effic.C3,
               "vax2.death.effic.C9" = vax2.death.effic.C4,
               "vax2.death.effic.C10" = vax2.death.effic.C5,               
               ## AZ infec symp
               vax1.beta.effic.A = vax1.beta.effic.A,
               vax2.beta.effic.A = vax2.beta.effic.A,
               'vax1.symp.effic.A1-6' = vax1.symp.effic.A1,
               'vax1.symp.effic.A7-10' = vax1.symp.effic.A2,
               'vax2.symp.effic.A1-6' = vax2.symp.effic.A1,
               'vax2.symp.effic.A7-10' = vax2.symp.effic.A2,
               ## AZ hosp D1
               "vax1.hosp.effic.A1-6" = vax1.hosp.effic.A1,
               "vax1.hosp.effic.A7" = vax1.hosp.effic.A2,
               "vax1.hosp.effic.A8" = vax1.hosp.effic.A3,
               "vax1.hosp.effic.A9" = vax1.hosp.effic.A4,
               "vax1.hosp.effic.A10" = vax1.hosp.effic.A5,
               ## AZ hosp D2
               "vax2.hosp.effic.A1-6" = vax2.hosp.effic.A1,
               "vax2.hosp.effic.A7" = vax2.hosp.effic.A2,
               "vax2.hosp.effic.A8" = vax2.hosp.effic.A3,
               "vax2.hosp.effic.A9" = vax2.hosp.effic.A4,
               "vax2.hosp.effic.A10" = vax2.hosp.effic.A5,
               ## AZ death D1
               "vax1.death.effic.A1-6" = vax1.death.effic.A1,
               "vax1.death.effic.A7" = vax1.death.effic.A2,
               "vax1.death.effic.A8" = vax1.death.effic.A3,
               "vax1.death.effic.A9" = vax1.death.effic.A4,
               "vax1.death.effic.A10" = vax1.death.effic.A5,
               ## AZ death D2
               "vax2.death.effic.A1-6" = vax2.death.effic.A1,
               "vax2.death.effic.A7" = vax2.death.effic.A2,
               "vax2.death.effic.A8" = vax2.death.effic.A3,
               "vax2.death.effic.A9" = vax2.death.effic.A4,
               "vax2.death.effic.A10" = vax2.death.effic.A5
               )
# p = p, gamma = gamm,
n.samples <- 500
set.seed(42)
n.params <- length(names(params))
samples <- matrix(0,nrow = n.samples,ncol = n.params)
samples <-as.data.frame(samples)
colnames(samples) <- names(params)

for (param in params) {
  distr<- param$distr
  if(distr == "unif"){
    samples[,param$param] <- runif(n.samples,min = param$vals[1],param$vals[2])
  }
  if(distr == "normal"){
    samples[,param$param] <- rnorm(n.samples,mean = param$vals[1],sd = param$vals[2])
  }
  if(distr == "gamma"){
    samples[,param$param] <- rgamma(n.samples,shape = param$vals[1],rate = param$vals[2])
  }
  if(distr == "exp"){
    samples[,param$param] <- rexp(n.samples,rate = param$vals[1])
  }
  if(distr == "beta"){
    samples[,param$param] <- rbeta(n.samples,param$vals[1],param$vals[2])
  }
  if(distr == "betabinom"){
    samples[,param$param] <- 1-rbeta(n.samples,param$vals[1],param$vals[2])/rbeta(n.samples,param$vals[3],param$vals[4])
  }
  if(distr == "betaneg"){
    samples[,param$param] <- 2*rbeta(n.samples,param$vals[1],param$vals[2])-1
  }
}
