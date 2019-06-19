library(boot)
source("./general.R")
library("pracma")
#setwd("/home/chris/ExperimentalPsychology/TermProjectErlingChris")
# RUN THIS FILE TO GENERATE SIMULATED DATA FOR PARAMETER RECOVERY

# Simulation parameters
seed <- 211196
num_subjs  <- 36 # number of subjects | just setting a number larger than subjects in trials

model_name = "author2"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 4
set.seed(seed)   # always set a seed number for this homework!


# get pars vector
paramList = c("beta_mu", "beta", "beta_sig", "theta", "ddb", "PointPosteriors", "PredictedResponse", "log_lik")
dataList = get_dataList2()


output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)



# Set seed


params = rstan::extract(output)
beta1 <- c()
beta2 <- c()
beta3 <- c()
beta4 <- c()
skip <- 0
for (i in 1:36){
  if (i > 13){
    skip <- 1
  }
  beta1[i + skip] <- mean(params$beta[,i,1])
  beta2[i + skip] <- mean(params$beta[,i,2])
  beta3[i + skip] <- mean(params$beta[,i,3])
  beta4[i + skip] <- mean(params$beta[,i,4])
}


# True parameters
simul_pars <- data.frame(
                         beta1 = beta1,
                         beta2 = beta2,
                         beta3 = beta3,
                         beta4 = beta4
                         )


# read the data file | because we'll just use the presented conditions cus they are equally devided
dat_1 = read.table("./data/AllSubjectsProcessed.tsv", header=T, sep="\t")

# For storing simulated choice data for all subjects
all_data = data.frame(Trial = dat_1$Trial,
                     RiskType = dat_1$RiskType,
                     RewardType = dat_1$RewardType,
                     ResponseType = dat_1$ResponseType,
                     Reward = dat_1$Reward,
                     Shock = dat_1$Shock,
                     SubjID = dat_1$SubjID)

preprocessedData <- get_dataList2()

responses = rep(0, nrow(dat_1))
for (i in 1:nrow(dat_1)) {
    # Individual-level (i.e. per subject) parameter values
    beta <- c(simul_pars$beta1[dat_1$SubjID[i]], simul_pars$beta2[dat_1$SubjID[i]], simul_pars$beta3[dat_1$SubjID[i]], simul_pars$beta4[dat_1$SubjID[i]])
    X <- unlist(preprocessedData$X[i])
    pGamble <- inv.logit(dot(X, beta))

    ResponseType <- rbinom(1,1,pGamble)

    # Append current subject with all subjects' data
    responses[i] = ResponseType
}

all_data$ResponseType = responses

# Write out data
write.table(all_data, file = "./simulators/simul_data_2_author2.txt", row.names = F, col.names = T, sep = "\t")
# write out parameters
write.table(simul_pars[-c(14),], file ="./simulators/simul_param_2_author2.txt", row.names=F, col.names = T, sep = "\t")

