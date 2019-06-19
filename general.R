# Programmed by Erling Ljunggren, cus that's important to know.

library(rstan)
library(loo)

library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
rstan_options(auto_write = TRUE)

# set seed
set.seed(21111996) # such eggs of easter ^^

sample_model <- function(model_name, dataList, paramList, iterations, warmups, chains, init="random"){
    # run!
    if (!file.exists(paste("./stanfits/", model_name, ".rds",sep=""))){
        print("fitting stan model")
        output = stan(paste("./models/", model_name, ".stan", sep=""),
              data = dataList, pars = paramList, init=init,
              iter = iterations, warmup=warmups, chains=chains, cores=chains)
        saveRDS(output, paste("./stanfits/", model_name, ".rds",sep=""))
    } else {
        print("recovering stan model")
        output = readRDS(paste("./stanfits/", model_name, ".rds",sep=""))
    }
    return(output)
}

get_dataList <- function(path="./data/AllSubjectsProcessed.tsv"){
    # read the data file
    dat_1 = read.table(path, header=T, sep="\t")

    allSubjs = unique(dat_1$SubjID)  # all subject IDs
    N = length(allSubjs)                   # number of subjects
    # T = table(dat_1)[1]     # number of trials per subject (=108)
    # T = nrow(dat_1)     # number of trials per subject (=108)
    T = 108
    numIter = 108           # max number of iterations to find global minimum values
    numPars = 3             # number of parameters


    # data frames for fill-in
    # -1 will be values for not used fields |  NaN
    Tsubj = array(0, c(N))
    RewardType = array(0, c(N,T))
    RiskType = array(0, c(N,T))
    ResponseType = array(0, c(N,T))
    Shock = array(0, c(N,T))

    #fill in with data
    for (n in 1:N){
        subjdat = subset(dat_1, SubjID == allSubjs[n])
        trials = nrow(subjdat)
        Tsubj[n] = trials
        RiskType[n,1:trials] = subjdat$RiskType
        RewardType[n,1:trials] = subjdat$RewardType
        ResponseType[n,1:trials] = subjdat$ResponseType
        Shock[n,1:trials] = subjdat$Shock
    }

    dataList <- list(
                     N       = N,
                     T       = T, #108
                     Tsubj   = Tsubj, # <= 108
                     RiskType = RiskType,
                     RewardType = RewardType,
                     ResponseType =  ResponseType,
                     Shock = Shock
    )
    return(dataList)
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

PPC <- function(output, dataList){
    print("running PPC")
    #pdf(paste("./plots/", model_name, "_prediction.pdf", sep=""))
    params = rstan::extract(output)

    # Set up the vectors
    true_gambling <- c("tSafe", "tGamble")
    prediction <- c("pSafe","pGamble")

    total_values <- c(0,0,0,0)

    for(n in 1:dataList$N){
                 #true: S, G  pred:
        pred_values <- c(0,0, # S
                         0,0) # G
        for(i in 1:dataList$Tsubj[n]){
            yPred <- getmode(params$PredictedResponse[,n,i])
            yTrue <- dataList$ResponseType[n,i]
            index = yPred*2 + yTrue + 1
            pred_values[index] = pred_values[index] + 1
        }
        total_values = total_values + pred_values

        # Create the data frame
        df <- expand.grid(true_gambling, prediction)
        df$value <- pred_values

        #Plot the Data
        g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
            theme(legend.position="none") + xlab("") + ylab("") + ggtitle(paste(model_name, "- Subject:",n)) +
            scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))
        print(g)
    }
    #Plot the Data for all subjects together
    df <- expand.grid(true_gambling, prediction)
    df$value <- total_values
    g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
        theme(legend.position="none") + xlab("") + ylab("") + ggtitle(paste(model_name,"- All Subjects")) +
        scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))
    print(g)
    print((df[1,3] + df[4,3] ) / sum(df[,3]))
    
}

getBIC <- function(ll, num_params, samples){
   (-2) * mean(ll) + log(samples) * num_params
}

BIC <- function(output, dataList, num_params){
    print("running BIC")
    #pdf(paste("./plots/", model_name, "_BIC.pdf", sep=""))
    params = rstan::extract(output)
    individ_BIC = rep(0, dataList$N)
    for(n in 1:dataList$N){
        individ_BIC[n] = getBIC(params$log_lik[,n], num_params, dataList$Tsubj[n])
    }
    df_individ_BIC = data.frame(BIC = individ_BIC, id = 1:dataList$N)
    g <- ggplot(data = df_individ_BIC, mapping = aes(x=id, y=BIC)) + geom_point() + ggtitle(paste(model_name, "- BIC per subject | total BIC:", sum(individ_BIC)))
    print(g)
}

BIC_2 <- function(output, dataList, num_params){
  print("running BIC")
  #pdf(paste("./plots/", model_name, "_BIC.pdf", sep=""))
  params = rstan::extract(output)
  individ_BIC = rep(0, dataList$n_subj)
  for(n in 1:dataList$n_subj){
    individ_BIC[n] = getBIC(params$log_lik[,n], num_params, length(dataList$ix[dataList$ix == n]))
  }
  df_individ_BIC = data.frame(BIC = individ_BIC, id = 1:dataList$n_subj)
  g <- ggplot(data = df_individ_BIC, mapping = aes(x=id, y=BIC)) + geom_point() + ggtitle(paste(model_name, "- BIC per subject | total BIC:", sum(individ_BIC)))
  print(g)
}

LOOIC <- function(output){
    print("running LOOIC")
    #pdf(paste("./plots/", model_name, "_LOOIC.pdf", sep=""), height=2, width=4.5)
    # Extract pointwise log-likelihood and compute LOO
    log_lik <- extract_log_lik(output, merge_chains = FALSE)

    # as of loo v2.0.0 we can optionally provide relative effective sample sizes
    # when calling loo, which allows for better estimates of the PSIS effective
    # sample sizes and Monte Carlo error
    r_eff <- relative_eff(exp(log_lik))

    looic <- loo::loo(log_lik, r_eff = r_eff, cores = 2)

    table <- tableGrob(looic$estimates)
    title <- textGrob(paste(model_name, "- LOOIC"), gp = gpar(fontsize = 20))
    padding <- unit(0.5,"line")
    table <- gtable_add_rows(
      table, heights = grobHeight(title) + padding, pos = 0
    )
    table <- gtable_add_grob(
      table, list(title),
      t = 1, l = 1, r = ncol(table)
    )
    grid.draw(table)
}


PPC_2 <- function(output, dataList, samples){
  print("running PPC")
  #pdf(paste("./plots/", model_name, "_prediction.pdf", sep=""))
  params = rstan::extract(output)

  # Set up the vectors
  true_gambling <- c("tSafe", "tGamble")
  prediction <- c("pSafe","pGamble")

  total_values <- c(0,0,0,0)

  for(n in 1:length(dataList$N)){
    #true: S, G  pred:
    pred_values <- c(0,0, # S
                     0,0) # G

    yPred <- getmode(params$PredictedResponse[,n])
    yTrue <- dataList$N[n]
    index = yPred*2 + yTrue + 1
    pred_values[index] = pred_values[index] + 1
    total_values = total_values + pred_values

    # Create the data frame
    df <- expand.grid(true_gambling, prediction)
    df$value <- pred_values

    #Plot the Data
    #g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
    #  theme(legend.position="none") + xlab("") + ylab("") + ggtitle(paste("Subject:",n)) +
    #  scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))
    #print(g)
  }
  #Plot the Data for all subjects together
  df <- expand.grid(true_gambling, prediction)
  df$value <- total_values
  g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
    theme(legend.position="none") + xlab("") + ylab("") + ggtitle(paste("All Subjects")) +
    scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))
  print(g)
  print((df[1,3] + df[4,3] ) / sum(df[,3]))
}

get_dataList2 <- function(path="./data/AllSubjectsProcessed.tsv"){
  dat_1 = read.table(path, header=T, sep="\t")
  allSubjs = unique(dat_1$SubjID)  # all subject IDs
  N = length(allSubjs)
  T = 108


  Reward = array(0, c(N,T))
  ZScore = array(0, c(N,T))
  Tsubj = array(0, c(N))
  Zscore = rep(0, length(dat_1$SubjID))
  # This is where you need to look tomorrow
  for (n in 1:N){
    subjdat = subset(dat_1, SubjID == allSubjs[n])
    trials = nrow(subjdat)
    Tsubj[n] = trials
    Reward[n,1:trials] = subjdat$Reward
    ZScore[n,1:trials] = (subjdat$Reward - mean(subjdat$Reward)) / sd(subjdat$Reward)
  }
  Zscore = array(ZScore)
  X = list()
  for (n in 1:length(dat_1$SubjID)){
    subjdat = subset(dat_1, SubjID == dat_1$SubjID[n])
    X = c(X, list(c(1, as.integer(dat_1$RiskType[n] == 2), as.integer(dat_1$RiskType[n] == 3), (dat_1$Reward[n] - mean(subjdat$Reward)) / sd(subjdat$Reward))))
  }

  for (n in 1:length(dat_1$SubjID)){
    if(dat_1$SubjID[n] > 14){
      dat_1$SubjID[n] = dat_1$SubjID[n] - 1
    }
  }


  dataList2 = get_dataList()

  dataList <- list(
    n_obs       = length(dat_1$SubjID),
    n_pred       = 4,
    n_subj   = N, # <= 108
    N = dat_1$ResponseType,
    X = X,
    ix = dat_1$SubjID

    #RiskType = dataList2$RiskType,
    #RewardType = dataList2$RewardType,
    #ResponseType =  dataList2$ResponseType,
    #Shock = dataList2$Shock,
    #Reward = Reward,
    #ZScore = ZScore

  )
  return(dataList)
}

compare_PA <- function(output){
    print("comparing PA")
    #png(paste("./plots/", model_name, "_comparison_PA.png", sep=""), width=900, height=900, res=150)
    parameters = extract(output)
    PainAvoidance_low = data.frame(parameters$PainAvoidance[,,1])
    PainAvoidance_low = reshape::melt(PainAvoidance_low)
    PainAvoidance_med = data.frame(parameters$PainAvoidance[,,2])
    PainAvoidance_med = reshape::melt(PainAvoidance_med)
    PainAvoidance_hig = data.frame(parameters$PainAvoidance[,,3])
    PainAvoidance_hig = reshape::melt(PainAvoidance_hig)
    low <- ggplot(PainAvoidance_low, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("low risk") + theme(legend.position="none") + scale_x_continuous(limits=c(0,4))
    med <- ggplot(PainAvoidance_med, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("medium risk") + theme(legend.position="none") + scale_x_continuous(limits=c(0,4))
    hig <- ggplot(PainAvoidance_hig, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("high risk") + theme(legend.position="none") + scale_x_continuous(limits=c(0,4))
    grid.arrange(low,med,hig, ncol=1)
}

parameter_recovery <- function(model_name, paramList, iterations, warmups, chains, init="random"){
    dataList = get_dataList(path=paste("./simulators/simul_data_2_", model_name, ".txt", sep=""))

    if (!file.exists(paste("./simulators/stanfits/", model_name, ".rds",sep=""))){
        print("fitting stan model")
        output = stan(paste("./models/", model_name, ".stan", sep=""),
              data = dataList, pars = paramList, init=init,
              iter = iterations, warmup=warmups, chains=chains, cores=chains)
        saveRDS(output, paste("./simulators/stanfits/", model_name, ".rds",sep=""))
    } else {
        print("recovering stan model")
        output = readRDS(paste("./simulators/stanfits/", model_name, ".rds",sep=""))
    }
    # extract Stan fit object (parameters)
    parameters <- rstan::extract(output)

    #pdf(paste("./plots/", model_name, "_recovery.pdf", sep=""))
    true_params = read.table(paste("./simulators/simul_param_2_", model_name, ".txt", sep=""), header=T, sep="\t")

    RA_mean = apply(parameters$RiskAversion, 2, mean)
    RA_sd = apply(parameters$RiskAversion, 2, sd)

    # print(RA_mean)
    # print(true_params$RiskAversion)
    #
    df = data.frame(TrueParameters=true_params$RiskAversion, SimulatedParameters=RA_mean, errY=RA_sd)
    l1_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=0, y=0, xend=2, yend=2, colour=1), show.legend=FALSE)
    print(l1_plot + ggtitle("(Risk Aversion)Rho Estimated Parameters vs Recovered Parameters"))

    PA_l_mean = apply(parameters$PainAvoidance[,,1], 2, mean)
    PA_l_sd = apply(parameters$PainAvoidance[,,1], 2, sd)
    df = data.frame(TrueParameters=true_params$PainAvoidance_low, SimulatedParameters=PA_l_mean, errY=PA_l_sd)
    l1_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=0, y=0, xend=1, yend=1, colour=1), show.legend=FALSE)
    print(l1_plot + ggtitle("(Low Pain Avoidance)Lambda1 Estimated Parameters vs Recovered Parameters"))

    PA_m_mean = apply(parameters$PainAvoidance[,,2], 2, mean)
    PA_m_sd = apply(parameters$PainAvoidance[,,2], 2, sd)
    df = data.frame(TrueParameters=true_params$PainAvoidance_med, SimulatedParameters=PA_m_mean, errY=PA_m_sd)
    l2_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=0, y=0, xend=3, yend=3, colour=1), show.legend=FALSE)
    print(l2_plot + ggtitle("(Medium Pain Avoidance)Lambda2 Estimated Parameters vs Recovered Parameters"))

    PA_h_mean = apply(parameters$PainAvoidance[,,3], 2, mean)
    PA_h_sd = apply(parameters$PainAvoidance[,,3], 2, sd)
   # plot(true_params$PainAvoidance_high, PA_h_mean,
  #       xlab="true params", ylab="recovered params")
    #arrows(true_params$PainAvoidance_high, PA_h_mean-PA_h_sd, true_params$PainAvoidance_high, PA_h_mean+PA_h_sd, length=0.05, angle=90, code=3)
    #lines(c(0,5), c(0,5))
    df = data.frame(TrueParameters=true_params$PainAvoidance_high, SimulatedParameters=PA_h_mean, errY=PA_h_sd)
    l3_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=0, y=0, xend=3, yend=3, colour=1), show.legend=FALSE)
    print(l3_plot + ggtitle("(High Pain Avoidance)Lambda3 Estimated Parameters vs Recovered Parameters"))
    
    tau_mean = apply(parameters$tau, 2, mean)
    tau_sd = apply(parameters$tau, 2, sd)
    # plot(true_params$PainAvoidance_high, PA_h_mean,
    #       xlab="true params", ylab="recovered params")
    #arrows(true_params$PainAvoidance_high, PA_h_mean-PA_h_sd, true_params$PainAvoidance_high, PA_h_mean+PA_h_sd, length=0.05, angle=90, code=3)
    #lines(c(0,5), c(0,5))
    df = data.frame(TrueParameters=true_params$tau, SimulatedParameters=tau_mean, errY=tau_sd)
    l3_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=5, y=5, xend=15, yend=15, colour=1), show.legend=FALSE)
    print(l3_plot + ggtitle("Tau Estimated Parameters vs Recovered Parameters"))


}


parameter_recovery2 <- function(model_name, paramList, iterations, warmups, chains, init="random"){
  dataList = get_dataList2(path=paste("./simulators/simul_data_2_", model_name, ".txt", sep=""))
  
  if (!file.exists(paste("./simulators/stanfits/", model_name, ".rds",sep=""))){
    print("fitting stan model")
    output = stan(paste("./models/", model_name, ".stan", sep=""),
                  data = dataList, pars = paramList, init=init,
                  iter = iterations, warmup=warmups, chains=chains, cores=chains)
    saveRDS(output, paste("./simulators/stanfits/", model_name, ".rds",sep=""))
  } else {
    print("recovering stan model")
    output = readRDS(paste("./simulators/stanfits/", model_name, ".rds",sep=""))
  }
  # extract Stan fit object (parameters)
  parameters <- rstan::extract(output)
  
  #pdf(paste("./plots/", model_name, "_recovery.pdf", sep=""))
  true_params = read.table(paste("./simulators/simul_param_2_", model_name, ".txt", sep=""), header=T, sep="\t")
  
  beta1_mean = apply(parameters$beta[,,1], 2, mean)
  beta1_sd = apply(parameters$beta[,,1], 2, sd)
  
  # print(RA_mean)
  # print(true_params$RiskAversion)
  #
  df = data.frame(TrueParameters=true_params$beta1, SimulatedParameters=beta1_mean, errY=beta1_sd)
  beta_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=-3, y=-3, xend=15, yend=15, colour=1), show.legend=FALSE)
  print(beta_plot + ggtitle("Beta1 Estimated Parameters vs Recovered Parameters"))
  
  beta2_mean = apply(parameters$beta[,,2], 2, mean)
  beta2_sd = apply(parameters$beta[,,2], 2, sd)
  df = data.frame(TrueParameters=true_params$beta2, SimulatedParameters=beta2_mean, errY=beta2_sd)
  beta_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=1, y=1, xend=-7, yend=-7, colour=1), show.legend=FALSE)
  print(beta_plot + ggtitle("Beta2 Estimated Parameters vs Recovered Parameters"))
  
  beta3_mean = apply(parameters$beta[,,3], 2, mean)
  beta3_sd = apply(parameters$beta[,,3], 2, sd)
  df = data.frame(TrueParameters=true_params$beta3, SimulatedParameters=beta3_mean, errY=beta3_sd)
  beta_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=0, y=0, xend=-10, yend=-10, colour=1), show.legend=FALSE)
  print(beta_plot + ggtitle("Beta3 Estimated Parameters vs Recovered Parameters"))
  
  beta4_mean = apply(parameters$beta[,,4], 2, mean)
  beta4_sd = apply(parameters$beta[,,4], 2, sd)
  df = data.frame(TrueParameters=true_params$beta4, SimulatedParameters=beta4_mean, errY=beta4_sd)
  beta_plot <- ggplot(data=df, aes(x=TrueParameters, y=SimulatedParameters, show.legend=F)) + geom_point() + geom_errorbar(aes(ymin=SimulatedParameters-errY, ymax=SimulatedParameters+errY)) + geom_segment(aes(x=0, y=0, xend=10, yend=10, colour=1), show.legend=FALSE)
  print(beta_plot + ggtitle("Beta4 Estimated Parameters vs Recovered Parameters"))
}
