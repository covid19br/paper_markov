source("functions/vac_params_all_vac.R")
source("functions/markov_all_vac.R")
source("functions/sensitivity_analysis_C2.R")
source("functions/samples_to_parameters.R")
###############to run###########################
mat = Matrix(matrix(0, length(Y2), length(Y2)), sparse = TRUE)
rownames(mat) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
colnames(mat) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
parameters <-samples_to_parameters(samples[1,],parameters)
mat2<- add_probabilities(mat,parameters)
nmax = MAX.TIME.DAYS # Tempo de passos
res <- run_markov(mat2,nmax,parameters)
# res <- group_age_classes(res)
y <- group_age_classes(res)[nrow(res),]



plot_solution(res,states = c("D","DvP","DwP","DvA","DwA","DvC","DwC"),group.vac = F)
