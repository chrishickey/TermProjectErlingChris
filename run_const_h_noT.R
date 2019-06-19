source("./general.R")
library(ggplot2)
library(rstan)
#setwd("")

# MODEL REFERENCED IN FINAL REPORT
model_name = "const_h_noT"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 4

# get pars vector
paramList = c("RiskAversion","PainAvoidance","mu_p", "sigma_p", "mu_RiskAversion","mu_PainAvoidance","log_lik","PredictedResponse")
dataList = get_dataList()

output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)

BIC(output, dataList, 4)
# Run for BBC
PPC(output, dataList)
# Run for LOOIC score
LOOIC(output)
# Run for posreriors (click back through plots to see all (first plot will be blank due to no tau value))
compare_PA(output)

# ## traceplot (can run for traceplots if you want) 
# pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("RiskAversion"))
traceplot(output, pars=c("PainAvoidance"))
traceplot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))
#
#
# # posterior plots
# pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("RiskAversion"))
stan_plot(output, pars=c("PainAvoidance"))
stan_plot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))
stan_dens(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))

# Run for parameter recovery
parameter_recovery(model_name, paramList, iterations, warmups, 4)
