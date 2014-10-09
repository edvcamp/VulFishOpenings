### Simulation code to represent vulnerable pool dynamics from CR vul exchange for ed working spreadsheet
### Looks at a single lake for a year, assuming fixed biologial and fisheries parms
### Estimates the open sequence.  Note, can change the length of the sequence estimated
rm(list=ls())

### Parameter Section ###
unit_time <- 360        #days in a year, almost.  reduced to get an integer of 12 30-day months or 36 10-day periods
No <- 1000              #Initial abundance
M  <- 0.5/unit_time     #natural mortality
q  <- 0.05              #catchability
v1 <- 2/unit_time       #exchange rate to the vulnerable pool
v2 <- 3/unit_time       #exchange rate to the invulnerable pool
rel_prop <- 1           #proportion released
rel_surv <- 0.9         #release survival rate
recov <- 0.2            #rate of recovery from the refractory pool
pvul_recov <- 0.5       #Propotion of recovered fish that go to the vulnerable pool
surv <- exp(-M)         #survival
sd_obs <- .2            #standard deviation of observation error of fishers, such that greater values lead to greater divergence from perfect knowledge
cpue_base <- 0.5          #catch rate required for a satisfaction or utility of 1
val_power <- 2          #power function of how satisaction increases with catch rate
max_eff_wknd <- 10      #max effort possible on a weekend
max_eff_wkdy <- 10      #max effort possible on a weekday
vul_init <- v1/(v1+v2)*No     #initial number of vulnerable fish
max_cpue <- q*vul_init        #max cpue, don't think this is used
cpue_half <- 3                #inflection point for logistic effort
sdcpue <- 5                   #sigma of logistic effort

### Procedure Section ###
fun <- function(theta){
  open=exp(theta)/(1+exp(theta))         #theta has been logit transformed to 0, now back transformed to 0.5
  #vectors
  day <- seq(1,unit_time, by=1)
  vul=NULL; invul = NULL; refract = NULL; max_eff = NULL; pmax_eff = NULL;
  effort = NULL; catch = NULL; vul_surv = NULL; invul_surv = NULL; refract_surv = NULL;
  alive = NULL; eff_cum = NULL; cpue = NULL; value = NULL; obs_err =NULL;
  open_seq=NULL;
  
  #independent vectors
  max_eff = rep(c(max_eff_wknd, rep(max_eff_wkdy,5),max_eff_wknd),length.out=unit_time)   #gives the max effort for any given day, assuming potential difference between weekends and weekdays
  open1 = ifelse(open>0.5,1,0)             #0.5 and lower will become "0's", anything greater will become "1's"
  open_seq = rep(open1,unit_time/length(theta))                 #repeats the monthly pattern 12 times

  #initial values
  vul[1] = vul_init
  invul[1] = No-vul_init
  refract[1] = No-(vul[1]+invul[1])
  obs_err[1] =    1 #rnorm(1,1, sd_obs)     #set to 1 now, switch to rnorm for observation error and departure from perfect knowledge
  pmax_eff[1] = 1/(1+exp(-(q*(obs_err[1]*vul[1])-cpue_half)/sdcpue))
  effort[1] = max_eff[1]*pmax_eff[1]*open_seq[1]
  eff_cum[1] = effort[1]
  catch[1] = vul[1]*(1-exp(-q*effort[1]))
  vul_surv[1] = surv*(vul[1]-catch[1])
  invul_surv[1] = invul[1]*surv
  refract_surv[1] = (refract[1] + catch[1]*rel_prop*rel_surv)*surv
  alive = vul[1]+invul[1]+refract[1]
  cpue[1] = ifelse(effort[1]>0,catch[1]/effort[1],0)
  value[1] = effort[1]*(cpue[1]/cpue_base)^val_power
  
  #time dynamic values
  for(i in 2:length(day)){
      vul[i] = vul_surv[i-1] - v2*vul_surv[i-1] + v1*invul_surv[i-1] + recov*refract_surv[i-1]*pvul_recov
      invul[i] = invul_surv[i-1] - v1*invul_surv[i-1] + v2*vul_surv[i-1] + refract_surv[i-1]*recov*(1-pvul_recov)
      refract[i] = refract_surv[i-1]*(1-recov)
      obs_err[i] = 1  #rnorm(1,1,sd_obs)
      pmax_eff[i] = 1/(1+exp(-(q*(obs_err[i]*vul[i])-cpue_half)/sdcpue))
      effort[i] = max_eff[i]*pmax_eff[i]*open_seq[i]
      catch[i] = vul[i]*(1-exp(-q*effort[i]))
      vul_surv[i] = surv*(vul[i]-catch[i])
      invul_surv[i] = invul[i]*surv
      refract_surv[i] = (refract[i] + catch[i]*rel_prop*rel_surv)*surv
      alive[i] = vul[i]+invul[i]+refract[i]
      eff_cum[i] = eff_cum[i-1]+effort[i]                                       #modified from Carl's sheet, I think this makes more sense
      cpue[i] = ifelse(effort[i]>0,catch[i]/effort[i],0)
      value[i] = effort[i]*(cpue[i]/cpue_base)^val_power
  }
  avg_cpue <- mean(cpue[which(effort>0)])                                          #average cpue where effort>0,
  avg_cpue_x_100 <- (mean(cpue[which(effort>0)]))*100                              #*100 to make plotting easier
  tot_catch <- sum(catch);   tot_effort <- sum(effort);   tot_val <- sum(value);   #response metrics we might care about
  #metrics <- as.matrix(data.frame(avg_cpue_x_100, tot_catch, tot_effort, tot_val)) #dumb wrappers to make barplots easy
  -1*tot_val
}

ii=0.5                                         #initial values
theta=rep(log(ii/(1-ii)),unit_time/12)         #logit transformed starting values (turns .5 to 0)
fit=optim(theta, fun, method="SANN")                          #should estimate the open sequence
fit1 = optim(fit$par, fun, method="SANN")                     #re-estimate the open sequence, using initial estimates and new starting values
#fit2 = optim(fit1$par, fun, method="SANN")                   #can be repeated any number of times
#fit3 = optim(fit2$par, fun, method="SANN")
parms = fit1$par                               #label parameters
bt_parms <- exp(fit1$par)/(1+exp(fit1$par))    #back transform parameters to look at (remember, values >0.5 =1, else 0)
bt_parms


### Report Section ###
#This will take the estimated parameters from above, and feed them through the model to provide whatever response metrics are desired
fun1 <- function(theta){
  open=exp(theta)/(1+exp(theta))
  #vectors
  day <- seq(1,unit_time, by=1)
  vul=NULL; invul = NULL; refract = NULL; max_eff = NULL; pmax_eff = NULL;
  effort = NULL; catch = NULL; vul_surv = NULL; invul_surv = NULL; refract_surv = NULL;
  alive = NULL; eff_cum = NULL; cpue = NULL; value = NULL; obs_err =NULL;
  open_seq=NULL;

  #independent vectors
  max_eff = rep(c(max_eff_wknd, rep(max_eff_wkdy,5),max_eff_wknd),length.out=unit_time)   #gives the max effort for any given day, assuming potential difference between weekends and weekdays
  open1 = ifelse(open>0.5,1,0)
  open_seq = rep(open1,unit_time/length(theta))

  #initial values
  vul[1] = vul_init
  invul[1] = No-vul_init
  refract[1] = No-(vul[1]+invul[1])
  obs_err[1] =    1 #rnorm(1,1, sd_obs)     #set to 1 now, switch to rnorm for observation error and departure from perfect knowledge
  pmax_eff[1] = 1/(1+exp(-(q*(obs_err[1]*vul[1])-cpue_half)/sdcpue))
  effort[1] = max_eff[1]*pmax_eff[1]*open_seq[1]
  eff_cum[1] = effort[1]
  catch[1] = vul[1]*(1-exp(-q*effort[1]))
  vul_surv[1] = surv*(vul[1]-catch[1])
  invul_surv[1] = invul[1]*surv
  refract_surv[1] = (refract[1] + catch[1]*rel_prop*rel_surv)*surv
  alive = vul[1]+invul[1]+refract[1]
  cpue[1] = ifelse(effort[1]>0,catch[1]/effort[1],0)
  value[1] = effort[1]*(cpue[1]/cpue_base)^val_power

  #time dynamic values
  for(i in 2:length(day)){
      vul[i] = vul_surv[i-1] - v2*vul_surv[i-1] + v1*invul_surv[i-1] + recov*refract_surv[i-1]*pvul_recov
      invul[i] = invul_surv[i-1] - v1*invul_surv[i-1] + v2*vul_surv[i-1] + refract_surv[i-1]*recov*(1-pvul_recov)
      refract[i] = refract_surv[i-1]*(1-recov)
      obs_err[i] = 1  #rnorm(1,1,sd_obs)
      pmax_eff[i] = 1/(1+exp(-(q*(obs_err[i]*vul[i])-cpue_half)/sdcpue))
      effort[i] = max_eff[i]*pmax_eff[i]*open_seq[i]
      catch[i] = vul[i]*(1-exp(-q*effort[i]))
      vul_surv[i] = surv*(vul[i]-catch[i])
      invul_surv[i] = invul[i]*surv
      refract_surv[i] = (refract[i] + catch[i]*rel_prop*rel_surv)*surv
      alive[i] = vul[i]+invul[i]+refract[i]
      eff_cum[i] = eff_cum[i-1]+effort[i]                                       #modified from Carl's sheet, I think this makes more sense
      cpue[i] = ifelse(effort[i]>0,catch[i]/effort[i],0)
      value[i] = effort[i]*max(cpue[i]/cpue_base-1)^val_power
  }
  avg_cpue <- mean(cpue[which(effort>0)])                                          #average cpue where effort>0,
  avg_cpue_x_100 <- (mean(cpue[which(effort>0)]))*100                              #*100 to make plotting easier
  tot_catch <- sum(catch);   tot_effort <- sum(effort);   tot_val <- sum(value);   #response metrics we might care about
  metrics <- as.matrix(data.frame(avg_cpue_x_100, tot_catch, tot_effort, tot_val)) #dumb wrappers to make barplots easy

  list(day=day, vul=vul, invul=invul, refract=refract, pmax_eff=pmax_eff, effort=effort, catch=catch,
     vul_surv=vul_surv, invul_surv=invul_surv, refract_surv=refract_surv, alive=alive, eff_cum=eff_cum, cpue=cpue,
     value=value, avg_cpue_x_100=avg_cpue_x_100, tot_catch= tot_catch, tot_effort = tot_effort, tot_val=tot_val,
     metrics=metrics, open_seq=open_seq)

}

op1 = fun1(parms)
op1$open_seq
op1$tot_val