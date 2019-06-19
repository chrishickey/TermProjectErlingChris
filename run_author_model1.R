source("./general.R")
library(ggplot2)
library(rstan)

#setwd()

model_name = "author1"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 4

# get pars vector
paramList = c("beta_mu","theta", "ddb", "PointPosteriors", "PredictedResponse", "log_lik")
dataList = get_dataList2()


output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)

#parameters <- rstan::extract(output)


#BIC_2(output, dataList, 1)

# Run for LOOIC
LOOIC(output)

# Run for PPC
PPC_2(output, dataList, iterations-warmups)


## traceplot
pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("beta_mu"))
traceplot(output, pars=c("theta"))
traceplot(output, pars=c("ddb"))


# posterior plots
pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("beta_mu"))
stan_plot(output, pars=c("theta"))
stan_plot(output, pars=c("ddb"))

