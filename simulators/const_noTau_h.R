library(boot)
source("./general.R")
#setwd("/home/chris/ExperimentalPsychology/RA-pain-model") SET TO PROJECT ROOT
# RUN TO GENERATED DATA FOR PARAMETER RECOVERY

# Simulation parameters
seed <- 211196
num_subjs  <- 36 # number of subjects | just setting a number larger than subjects in trials

model_name = "const_h_noT"
print("running model")
print(model_name)


# Set seed
set.seed(seed)   # always set a seed number for this homework!

output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)

params = rstan::extract(output)
risk_aversion = c()
pain_avoidance_low <- c()
pain_avoidance_med <- c()
pain_avoidance_high <- c()
skip <- 0
for (i in 1:36){
  if (i > 13){
    skip <- 1
  }
  risk_aversion[i + skip] <- mean(params$RiskAversion[,i])
  pain_avoidance_low[i + skip] <- mean(params$PainAvoidance[,i,1])
  pain_avoidance_med[i + skip] <- mean(params$PainAvoidance[,i,2])
  pain_avoidance_high[i + skip] <- mean(params$PainAvoidance[,i,3])
  tau[i + skip] <- 12

}


# True parameters
simul_pars <- data.frame(
                         RiskAversion = risk_aversion,
                         PainAvoidance_low = pain_avoidance_low,
                         PainAvoidance_med = pain_avoidance_med,
                         PainAvoidance_high = pain_avoidance_high,
                         tau = tau
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

responses = rep(0, nrow(dat_1))
for (i in 1:nrow(dat_1)) {
    # Individual-level (i.e. per subject) parameter values
    RiskAversion <- simul_pars$RiskAversion[dat_1$SubjID[i]]
    PainAvoidance <- c(simul_pars$PainAvoidance_low[dat_1$SubjID[i]], simul_pars$PainAvoidance_med[dat_1$SubjID[i]], simul_pars$PainAvoidance_high[dat_1$SubjID[i]])
    tau <- simul_pars$tau[dat_1$SubjID[i]]

    evSafe  <- 0.01^RiskAversion

    evGamble <- (all_data$RewardType[i]*0.33)^RiskAversion - log(PainAvoidance[all_data$RiskType[i]] + 1)

    pGamble <- inv.logit(tau * (evGamble - evSafe))

    ResponseType <- rbinom(1,1,pGamble)

    # Append current subject with all subjects' data
    responses[i] = ResponseType
}

all_data$ResponseType = responses

# Write out data
write.table(all_data, file = "./simulators/simul_data_2_const_h_noT.txt", row.names = F, col.names = T, sep = "\t")
# write out parameters
write.table(simul_pars[-c(14),], file ="./simulators/simul_param_2_const_h_noT.txt", row.names=F, col.names = T, sep = "\t")
write.table(
            data.frame(
              RiskAversion = 0.5,
              PainAvoidance_low = 0.1,
              PainAvoidance_med = 0.5,
              PainAvoidance_high = 0.9,
              tau = 10
            ),
            file ="./simulators/simul_param_hyp_2_const_h_noT.txt", row.names=F, col.names = T, sep = "\t")
