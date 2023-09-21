install.packages("stats4")
install.packages("tidyverse")
install.packages("stargazer")
install.packages("readxl")
install.packages("dbplyr")
install.packages("xtable")


library(stats4)
library(tidyverse)
library(stargazer)
library(readxl)
library(dbplyr)
library(tidyverse)
library(xtable)

rm(list=ls())
setwd("/Users/mihai/Desktop/Caltech/SupernormalStimuli/sns data/")

################# prepare data ####################
data = read_excel("data sns w&c.xlsx",sheet = "R_testing_competence") # change here with "R_test_warmth" , "R_train_competence" and "R_test_competence"
## rename and clean data
data = data%>%rename(states = state)%>%
  group_by(user_id) %>%
  dplyr::mutate(ID = cur_group_id())

#data = data%>%rename(engage = correct, nengage = incorrect, trial = qn)

## count the number of a=1 and a=0 in each state for each individual
count_action = data%>%dplyr::select(ID, states, engage, nengage)%>%
  group_by(ID, states)%>%
  summarize(sum_e = sum(engage),
            sum_n = sum(nengage))

count_action = count_action%>%arrange(ID, states)
data = data%>%arrange(ID, trial)

nsubj = length(unique(data$ID))

engage = matrix(NA, nrow = nsubj, ncol = 3)
nengage = matrix(NA, nrow = nsubj, ncol = 3)
for (i in 1:nsubj){
  engage[i,1] = count_action$sum_e[(i-1)*3+1]
  engage[i,2] = count_action$sum_e[(i-1)*3+2]
  engage[i,3] = count_action$sum_e[(i-1)*3+3]
  nengage[i,1] = count_action$sum_n[(i-1)*3+1]
  nengage[i,2] = count_action$sum_n[(i-1)*3+2]
  nengage[i,3] = count_action$sum_n[(i-1)*3+3]  
}


ntrials = rowSums(engage)+rowSums(nengage)

#################### grid search to get a proper starting value ######################
# parameter for NH cost
kappa = seq(0.1,20,0.1)
# parameter for LLR cost
beta = seq(0.001,5,0.01)
# parameter for mutual info
rho = seq(0.1, 20,0.1)
nparam = max(length(kappa), length(beta), length(rho))


loglike_NH = rep(NA, length(kappa))
loglike_LLR = rep(NA, length(beta))
loglike_MI = rep(NA, length(rho))
initval_NH = rep(NA, nsubj)
initval_LLR = rep(NA, nsubj)
initval_MI = rep(NA, nsubj)

multinom=function(x1,x2){
  return(sum(log(1:max(1,(x1+x2))))-sum(log(c(1:max(1,x1),1:max(1,x2)))))
}

# incentive
r = c(-1,1,0.5);
# number of states
N = length(r);
# cost function
cost_llr = function(prob, k){
  summ = 0
  for (i in 1:N){
    for (j in 1:N){
      if (i != j){
        summ = summ+(k/abs(i-j)^2)*(prob[i]*log(prob[i]/prob[j])+(1-prob[i])*log((1-prob[i])/(1-prob[j])))
      }
    }
  }
  return(summ)
}

cost_nh = function(prob, k){
  summ = 0
  for (i in 1:(N-1)){
    summ = summ + prob[i]*log(prob[i]/(prob[i]+prob[i+1]))+prob[i+1]*log(prob[i+1]/(prob[i]+prob[i+1]))+(1-prob[i])*log((1-prob[i])/(2-prob[i]-prob[i+1]))+(1-prob[i+1])*log((1-prob[i+1])/(2-prob[i]-prob[i+1]))   
  }
  return(k*summ/N)
}

cost_mi = function(prob, k){
  summ = 0
  for (i in 1:N){
    summ = summ + (prob[i]/N)*log(prob[i]/sum(prob))+((1-prob[i])/N)*log((1-prob[i])/(N-sum(prob))) ### double-check if this is correct
  }
  return(k*summ)
}

# constraint matrix
conmat=rbind(diag(N),-diag(N))
convec=c(rep(0.000001,N), rep(-0.999999,N))

# grid search
for (indii in 1:nsubj){
  multinomconst = 0
  for (i in 1:N){
    multinomconst = multinomconst + multinom(engage[indii,i],nengage[indii,i])
  }
  for (j in 1:length(beta)){
    obj_llr = function(prob){
      return(-sum(r*prob)/N+cost_llr(prob, beta[j]))
    }
    optimout_llr=constrOptim(rep(1/20,N),obj_llr, grad = NULL, ui=conmat,ci=convec)
    loglike_LLR[j]= multinomconst + sum(engage[indii,]*log(optimout_llr$par))+sum(nengage[indii,]*(log(1-optimout_llr$par)))
  }
  for (j in 1:length(kappa)){
    obj_nh = function(prob){
      return(-sum(r*prob)/N+cost_nh(prob, kappa[j]))
    }
    obj_mi = function(prob){
      return(-sum(r*prob)/N+cost_mi(prob, rho[j]))
    }
    optimout_nh=constrOptim(rep(1/20,N),obj_nh, grad = NULL, ui=conmat,ci=convec)
    loglike_NH[j]= multinomconst + sum(engage[indii,]*log(optimout_nh$par))+sum(nengage[indii,]*(log(1-optimout_nh$par)))
    optimout_mi=constrOptim(rep(1/20,N),obj_mi, grad = NULL, ui=conmat,ci=convec)
    loglike_MI[j]= multinomconst + sum(engage[indii,]*log(optimout_mi$par))+sum(nengage[indii,]*(log(1-optimout_mi$par)))
  }
  initval_LLR[indii] = beta[which.max(loglike_LLR)]
  initval_NH[indii] = kappa[which.max(loglike_NH)]
  initval_MI[indii] = rho[which.max(loglike_MI)]
}

# check whether the initial value is at the boundary
initval_LLR == beta[1] 
initval_LLR == beta[length(beta)]
initval_NH == kappa[1] 
initval_NH == kappa[length(kappa)] 
initval_MI == rho[1] # 3 reaches the lower bound
initval_MI == rho[length(rho)]


## Run this section only if lower or upper bound are met and need to be re-estimated or drop subjects. For competence domain , so far we only obtain FALSE, so no need to run
# Alternatively re-estimate for either MI, NH or LLR. Below is an example for subjects 6, 12, 13 for NH
kappa.refined = seq(50,100,0.1)
loglike_NH.refined = rep(NA, 3)
initval_NH.refined = rep(NA, 3)
subj.refined = c(6,12,13)

for (indii in subj.refined){
  multinomconst = 0
  for (i in 1:N){
    multinomconst = multinomconst + multinom(engage[indii,i],nengage[indii,i])
  }
  
  for (j in 1:length(kappa.refined)){
    obj_nh = function(prob){
      return(-sum(r*prob)/N+cost_nh(prob, kappa.refined[j]))
    }
    optimout_nh=constrOptim(rep(1/20,N),obj_nh, grad = NULL, ui=conmat,ci=convec)
    loglike_NH.refined[j]= multinomconst + sum(engage[indii,]*log(optimout_nh$par))+sum(nengage[indii,]*(log(1-optimout_nh$par)))
  }
  initval_NH.refined[indii] = kappa.refined[which.max(loglike_NH.refined)]
}
## This ends the "refined" sub-section

# get lower and upper bound
lower_LLR = initval_LLR - 0.01
upper_LLR = initval_LLR + 0.01
lower_NH = initval_NH - 0.1
upper_NH = initval_NH + 0.1
lower_MI = initval_MI - 0.1
upper_MI = initval_MI + 0.1


#########################  estimation  ############################
# MLE
prob_llr = function(k){
  obj_llr = function(prob){
    return(-sum(r*prob)/N+cost_llr(prob, k))
  }
  optimout_llr=constrOptim(rep(1/20,N),obj_llr, grad = NULL, ui=conmat,ci=convec)
  return(optimout_llr$par)
}
prob_nh = function(k){
  obj_nh = function(prob){
    return(-sum(r*prob)/N+cost_nh(prob, k))
  }
  optimout_nh=constrOptim(rep(1/20,N),obj_nh, grad = NULL, ui=conmat,ci=convec)
  return(optimout_nh$par)
}
prob_mi = function(k){
  obj_mi = function(prob){
    return(-sum(r*prob)/N+cost_mi(prob, k))
  }
  optimout_mi=constrOptim(rep(1/20,N),obj_mi, grad = NULL, ui=conmat,ci=convec)
  return(optimout_mi$par)
}

par_llr = rep(NA, nsubj)
par_nh = rep(NA, nsubj)
par_mi = rep(NA, nsubj)
llout_llr = rep(NA, nsubj)
llout_nh = rep(NA, nsubj)
llout_mi = rep(NA, nsubj)

count_subj = 1:nsubj
#count_subj = count_subj[! count_subj %in% exclude]
for (indii in count_subj){
  multinomconst = 0
  for (i in 1:N){
    multinomconst = multinomconst + multinom(engage[indii,i],nengage[indii,i])
  }
  ll_llr = function(k){
    probvec = prob_llr(k)
    like = multinomconst + sum(engage[indii,]*log(probvec))+sum(nengage[indii,]*(log(1-probvec)))
    return(-like)
  }
  ll_nh = function(k){
    probvec = prob_nh(k)
    like = multinomconst + sum(engage[indii,]*log(probvec))+sum(nengage[indii,]*(log(1-probvec)))
    return(-like)
  }
  ll_mi = function(k){
    probvec = prob_mi(k)
    like = multinomconst + sum(engage[indii,]*log(probvec))+sum(nengage[indii,]*(log(1-probvec)))
    return(-like)
  }
  model_mi=mle(ll_mi,start=list(k=initval_MI[indii]), method="L-BFGS-B",lower=lower_MI[indii],upper=upper_MI[indii])
  par_mi[indii]=coef(model_mi)
  llout_mi[indii]=logLik(model_mi)[1]
  model_llr=mle(ll_llr,start=list(k=initval_LLR[indii]), method="L-BFGS-B",lower=lower_LLR[indii],upper=upper_LLR[indii])
  model_nh=mle(ll_nh,start=list(k=initval_NH[indii]), method="L-BFGS-B",lower=lower_NH[indii],upper=upper_NH[indii])
  par_llr[indii]=coef(model_llr)
  par_nh[indii]=coef(model_nh)
  llout_llr[indii]=logLik(model_llr)[1]
  llout_nh[indii]=logLik(model_nh)[1]
}

count_trial = rowSums(engage) + rowSums(nengage)
indioutframe = data.frame(matrix(nrow = nsubj, ncol = 7))
colnames(indioutframe) = c('ID', 'AICLLR', 'AICNH', 'AICMI', 'BICLLR', 'BICNH', 'BICMI')
indioutframe$AICLLR = llout_llr*(-2)+2
indioutframe$AICNH = llout_nh*(-2)+2
indioutframe$AICMI = llout_mi*(-2)+2
indioutframe$BICLLR = llout_llr*(-2)+log(count_trial)
indioutframe$BICNH = llout_nh*(-2)+log(count_trial)
indioutframe$BICMI = llout_mi*(-2)+log(count_trial)
indioutframe$ID = unique(data$ID)
indioutframe$user_id = unique(data$user_id)
AIC = indioutframe%>%dplyr::select(AICLLR,AICNH,AICMI)
AIC$minCol = base::apply(AIC, 1, which.min)
BIC = indioutframe%>%select(BICLLR,BICNH,BICMI)
BIC$minCol = base::apply(BIC, 1, which.min)
AIC = AIC%>%mutate(
  select_model.aic = case_when(
    minCol == 1 ~ 'LLR',
    minCol == 2 ~ 'NH',
    minCol == 3 ~ 'MI',
    minCol == "NA" ~ '',
  )
)
AIC <- AIC%>%dplyr::select(AICLLR, AICNH, AICMI, select_model.aic)
AIC <- AIC %>%
  mutate(select_model.aic = map_chr(select_model.aic, ~ toString(.)))
colnames(AIC) <- c("LLR", "NH", "MI", "Selected model")


BIC = BIC%>%mutate(
  select_model.bic = case_when(
    minCol == 1 ~ 'LLR',
    minCol == 2 ~ 'NH',
    minCol == 3 ~ 'MI'
  )
)
BIC <- BIC%>%dplyr::select(BICLLR, BICNH, BICMI, select_model.bic)
colnames(BIC) <- c("LLR", "NH", "MI", "Selected model")
table_aic <- xtable(AIC, digits = 4, caption = "AIC for different cost functions")
results_path = "/Users/mihai/Desktop/Caltech/SupernormalStimuli/warmth and competence/data Pavlovia"
html_output <- print(table_aic, type = "html", comment = FALSE, caption.placement = "top", printl.results = FALSE)
footnote <- "<p> Notes on missing values: unable to find a global maximum for at least one of the cost models</p>"
con <- file(file.path(results_path, "table_aic.html"), open = "wt")
cat(html_output, footnote, sep = "\n", file = con)
close(con)

# get subj-level choice probabilities
prob_llr_subj = matrix(NA, nsubj, N)
prob_nh_subj = matrix(NA, nsubj, N)
prob_mi_subj = matrix(NA, nsubj, N)
for (indii in count_subj){
  multinomconst = 0
  for (i in 1:N){
    multinomconst = multinomconst + multinom(engage[indii,i],nengage[indii,i])
  }
  obj_llr = function(prob){
    return(-sum(r*prob)/N+cost_llr(prob, par_llr[indii]))
  }
  optimout_llr=constrOptim(rep(1/20,N),obj_llr, grad = NULL, ui=conmat,ci=convec)
  obj_nh = function(prob){
    return(-sum(r*prob)/N+cost_nh(prob, par_nh[indii]))
  }
  obj_mi = function(prob){
    return(-sum(r*prob)/N+cost_mi(prob, par_mi[indii]))
  }
  optimout_nh=constrOptim(rep(1/20,N),obj_nh, grad = NULL, ui=conmat,ci=convec)
  optimout_mi=constrOptim(rep(1/20,N),obj_mi, grad = NULL, ui=conmat,ci=convec)
  prob_llr_subj[indii, ] = optimout_llr$par
  prob_nh_subj[indii, ] = optimout_nh$par
  prob_mi_subj[indii, ] = optimout_mi$par
}

# subj id for each best fit
best_model_llr = c(2)
best_model_nh = c(4,7,8,9,11,12,13,14)
best_model_mi = c(1,6,15)


# aggregate level
# grid search 
loglike_NH = rep(NA, length(kappa))
loglike_LLR = rep(NA, length(beta))
loglike_MI = rep(NA, length(rho))
multinomconst = 0
data.agg <- data%>%ungroup()%>%group_by(states)%>%summarise(sum_e = sum(engage), sum_n = sum(nengage), mean_e = mean(engage), mean_n = mean(nengage))
engage = data.agg$sum_e
nengage = data.agg$sum_n
for (i in 1:N){
  multinomconst = multinomconst + multinom(engage,nengage)
}
for (j in 1:length(beta)){
  obj_llr = function(prob){
    return(-sum(r*prob)/N+cost_llr(prob, beta[j]))
  }
  optimout_llr=constrOptim(rep(1/20,N),obj_llr, grad = NULL, ui=conmat,ci=convec)
  loglike_LLR[j]= multinomconst + sum(engage*log(optimout_llr$par))+sum(nengage*(log(1-optimout_llr$par)))
}
for (j in 1:length(kappa)){
  obj_nh = function(prob){
    return(-sum(r*prob)/N+cost_nh(prob, kappa[j]))
  }
  obj_mi = function(prob){
    return(-sum(r*prob)/N+cost_mi(prob, rho[j]))
  }
  optimout_nh=constrOptim(rep(1/20,N),obj_nh, grad = NULL, ui=conmat,ci=convec)
  loglike_NH[j]= multinomconst + sum(engage*log(optimout_nh$par))+sum(nengage*(log(1-optimout_nh$par)))
  optimout_mi=constrOptim(rep(1/20,N),obj_mi, grad = NULL, ui=conmat,ci=convec)
  loglike_MI[j]= multinomconst + sum(engage*log(optimout_mi$par))+sum(nengage*(log(1-optimout_mi$par)))
}
initval_LLR.agg = beta[which.max(loglike_LLR)]
initval_NH.agg = kappa[which.max(loglike_NH)]
initval_MI.agg = rho[which.max(loglike_MI)]

# estimation
multinomconst = 0
for (i in 1:N){
  multinomconst = multinomconst + multinom(engage,nengage)
}
ll_llr = function(k){
  probvec = prob_llr(k)
  like = multinomconst + sum(engage*log(probvec))+sum(nengage*(log(1-probvec)))
  return(-like)
}
ll_nh = function(k){
  probvec = prob_nh(k)
  like = multinomconst + sum(engage*log(probvec))+sum(nengage*(log(1-probvec)))
  return(-like)
}
ll_mi = function(k){
  probvec = prob_mi(k)
  like = multinomconst + sum(engage*log(probvec))+sum(nengage*(log(1-probvec)))
  return(-like)
}
model_mi=mle(ll_mi,start=list(k=initval_MI.agg), method="L-BFGS-B",lower=initval_MI.agg - 0.1,upper=initval_MI.agg+0.1)
par_mi.agg=coef(model_mi)
llout_mi.agg=logLik(model_mi)[1]
model_llr=mle(ll_llr,start=list(k=initval_LLR.agg), method="L-BFGS-B",lower=initval_LLR.agg - 0.01,upper=initval_LLR.agg + 0.01)
model_nh=mle(ll_nh,start=list(k=initval_NH.agg), method="L-BFGS-B",lower=initval_NH.agg-0.1,upper=initval_NH.agg+0.1)
par_llr.agg=coef(model_llr)
par_nh.agg=coef(model_nh)
llout_llr.agg=logLik(model_llr)[1]
llout_nh.agg=logLik(model_nh)[1]
# agg best fit: NH > LLR > MI

# get agg-level choice probabilities
for (indii in count_subj){
  multinomconst = 0
  for (i in 1:N){
    multinomconst = multinomconst + multinom(engage,nengage)
  }
  obj_llr = function(prob){
    return(-sum(r*prob)/N+cost_llr(prob, par_llr.agg))
  }
  optimout_llr=constrOptim(rep(1/20,N),obj_llr, grad = NULL, ui=conmat,ci=convec)
  obj_nh = function(prob){
    return(-sum(r*prob)/N+cost_nh(prob, par_nh.agg))
  }
  obj_mi = function(prob){
    return(-sum(r*prob)/N+cost_mi(prob, par_mi.agg))
  }
  optimout_nh=constrOptim(rep(1/20,N),obj_nh, grad = NULL, ui=conmat,ci=convec)
  optimout_mi=constrOptim(rep(1/20,N),obj_mi, grad = NULL, ui=conmat,ci=convec)
  prob_llr.agg = optimout_llr$par
  prob_nh.agg = optimout_nh$par
  prob_mi.agg = optimout_mi$par
} 




# Create a data frame from vectors
df <- data.frame(
  state = c(1, 2, 3),
  mean_e = data.agg$mean_e,
  prob_llr_agg = prob_llr.agg,
  prob_nh_agg = prob_nh.agg,
  prob_mi_agg = prob_mi.agg
)

ggplot(df, aes(x = state)) +
  # Add bar plot for mean_e
  geom_bar(aes(y = mean_e, fill = "Empirical"), stat = "identity") +
  # Add lines and points for the other vectors
  geom_line(aes(y = prob_llr_agg, color = "LLR"), size = 1) +
  geom_point(aes(y = prob_llr_agg, color = "LLR"), size = 3) +
  
  geom_line(aes(y = prob_nh_agg, color = "NH"), size = 1) +
  geom_point(aes(y = prob_nh_agg, color = "NH"), size = 3) +
  
  geom_line(aes(y = prob_mi_agg, color = "MI"), size = 1) +
  geom_point(aes(y = prob_mi_agg, color = "MI"), size = 3) +
  
  # Label axes
  labs(x = "State", y = "Choice Probability") +
  scale_fill_manual(values = c("Empirical" = "grey"),
                    name = "", 
                    breaks = "Empirical") +
  scale_color_manual(values = c("LLR" = "red", "NH" = "blue", "MI" = "green"),
                     name = "") +
  scale_x_continuous(breaks = df$state, labels = 1:3) +
  # Theme adjustments
  theme_minimal() +
  theme(legend.position = "right",
        axis.ticks = element_blank())




save.image("output_competence.Rdata")










