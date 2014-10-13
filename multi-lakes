### Simulation code to represent vulnerable pool dynamics from CR vul exchange for ed working spreadsheet
### Looks at a single lake for a year, assuming fixed biologial and fisheries parms, across different closure regimes

rm(list=ls())

### Parameter Section ###
lakereps <- 7               #number of lakes per distance class should equal or be a mulitple of the number of open sequences to exist simoltaneously
ndisclass <- 3             #number of distance classes
disclass <- c(10, 100, 250) #distance classes in kilometers     #c(10,100,250)
nlakes <- lakereps*ndisclass      #number of lakes total--currently set for the fewest number of lakes
lakedis <- rep(disclass, each=lakereps)    #vector of lake distances   rep(1,nlakes) #
maxdis <- 400                              #theoretical average max distance one would travel for a trip
grav_pow <- .2                             #power function of gravity model

unit_time <- 360        #days in a year
No <- rep(500,nlakes)  #rnorm(nlakes,1000,100)              #Initial abundance
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

cpue_base <- .5          #catch rate required for a satisfaction or utility of 1
val_power <- 1.5          #power function of how satisaction increases with catch rate

max_eff_wknd <- 10            #max effort possible on a weekend
max_eff_wkdy <- 10            #max effort possible on a weekday

vul_init <- v1/(v1+v2)*No       #initial number of vulnerable fish
max_cpue <- q*vul_init          #max cpue, don't think this is used
cpue_half <- 2                #inflection point for logistic effort
sdcpue <- 0.5                   #sigma of logistic effort

day <- seq(1,unit_time, by=1)

#open_seq=allsame_open

### Procedure Section ###
fun <- function(open_seq){
  #vectors
  dis_prop =  rep(1,nlakes) # ((maxdis-lakedis)/maxdis) # #   #
  cost = .2*lakedis
  vul=matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  invul = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  refract = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  max_eff = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  pmax_eff = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  effort = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  catch = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  vul_surv = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  invul_surv = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  refract_surv = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  alive = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  eff_cum = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  cpue = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  value = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  obs_err = matrix(NA,nrow=length(day), ncol=nlakes) #NULL;
  grav_wt = matrix(NA, nrow=length(day), ncol=nlakes)
  #open_seq= matrix(c(all_open, week_2, week_1),nrow=length(day), ncol=nlakes) #Note, the number of scenarios must be a factor of nlakes to work like this
  max_eff = matrix(rep(c(max_eff_wknd, rep(max_eff_wkdy,5),max_eff_wknd),length.out=unit_time),nrow=length(day), ncol=nlakes)
  #independent vectors
  #open_seq = week_2;
  #max_eff = rep(c(max_eff_wknd, rep(max_eff_wkdy,5),max_eff_wknd),length.out=unit_time)   #gives the max effort for any given day, assuming potential difference between weekends and weekdays

  #initial values
  vul[1,] = vul_init
  invul[1,] = No[1]-vul_init
  refract[1,] = No[1]-(vul[1,]+invul[1,])
  obs_err[1,] =    1 #rnorm(1,1, sd_obs)     #set to 1 now, switch to rnorm for observation error and departure from perfect knowledge
  #pmax_eff[1,] = 1/(1+exp(-(q*(obs_err[1,]* mean(vul[1,]))-cpue_half)/sdcpue))
  pmax_eff[1,] = 1/(1+exp(-(q*obs_err[1,]* mean(vul[1,]*dis_prop) -cpue_half)/sdcpue))
  grav_wt[1,] = ((q*vul[1,])/cost)^grav_pow
  #effort[1,] = ((sum(max_eff[1,])*pmax_eff[1,]*open_seq[1,])*(vul[1,]/sum(vul[1,])))*dis_prop^grav_pow
  effort[1,] = sum(max_eff[1,]) * pmax_eff[1,] * open_seq[1,] * grav_wt[1,] / sum(grav_wt[1,])
  eff_cum[1,] = effort[1,]
  catch[1,] = vul[1,]*(1-exp(-q*effort[1,]))
  vul_surv[1,] = surv*(vul[1,]-catch[1,])
  invul_surv[1,] = invul[1,]*surv
  refract_surv[1,] = (refract[1,] + catch[1,]*rel_prop*rel_surv)*surv
  alive[1,] = vul[1,]+invul[1,]+refract[1,]
  cpue[1,] = ifelse(effort[1,]>0,catch[1,]/effort[1,],0)
  value[1,] = effort[1,]*(cpue[1,]/cpue_base)^val_power
  #value[i] = effort[i]*max(0,cpue[i]/cpue_base-1)^val_power
  #time dynamic values
  for(i in 2:length(day)){
      vul[i,] = vul_surv[i-1,] - v2*vul_surv[i-1,] + v1*invul_surv[i-1,] + recov*refract_surv[i-1,]*pvul_recov
      invul[i,] = invul_surv[i-1,] - v1*invul_surv[i-1,] + v2*vul_surv[i-1,] + refract_surv[i-1,]*recov*(1-pvul_recov)
      refract[i,] = refract_surv[i-1,]*(1-recov)
      obs_err[i,] = 1  #rnorm(1,1,sd_obs)
      #pmax_eff[i,] = 1/(1+exp(-(q*(obs_err[i,]*mean(vul[i,]))-cpue_half)/sdcpue))
      pmax_eff[i,] = 1/(1+exp(-(q*obs_err[i,]* mean(vul[i,]*dis_prop)-cpue_half)/sdcpue))
      grav_wt[i,] = ((q*vul[i,])/cost)^grav_pow
      #effort[i,] = ((sum(max_eff[i,])*pmax_eff[i,]*open_seq[i,])*(vul[i,]/sum(vul[i,])))*dis_prop^grav_pow
      effort[i,] = sum(max_eff[i,]) * pmax_eff[i,] * open_seq[i,] * grav_wt[i,] / sum(grav_wt[i,])
      catch[i,] = vul[i,]*(1-exp(-q*effort[i,]))
      vul_surv[i,] = surv*(vul[i,]-catch[i,])
      invul_surv[i,] = invul[i,]*surv
      refract_surv[i,] = (refract[i,] + catch[i,]*rel_prop*rel_surv)*surv
      alive[i,] = vul[i,]+invul[i,]+refract[i,]
      eff_cum[i,] = eff_cum[i-1,]+effort[i,]                                       #modified from Carl's sheet, I think this makes more sense
      cpue[i,] = ifelse(effort[i,]>0,catch[i,]/effort[i,],0)
      value[i,] = effort[i,]*max(0,cpue[i,]/cpue_base-1)^val_power              #Using this interpretation, there is 0 sat_catch at cpue = 0.5
      #value[i] = effort[i]*max(0,cpue[i]/cpue_base-1)^val_power
      #value[i,] = effort[i,]*(cpue[i,]/cpue_base)^val_power
  }
  avg_cpue <- mean(cpue[which(effort>0)])                                          #average cpue where effort>0,
  avg_cpue_x_100 <- (mean(cpue[which(effort>0)]))*100                              #*100 to make plotting easier
  tot_catch <- sum(catch);   tot_effort <- sum(effort);   val_lake <- colSums(value);   #response metrics we might care about
  tot_val=sum(val_lake);

  list(day=day, vul=vul, invul=invul, refract=refract, pmax_eff=pmax_eff, effort=effort, catch=catch,
     vul_surv=vul_surv, invul_surv=invul_surv, refract_surv=refract_surv, alive=alive, eff_cum=eff_cum, cpue=cpue,
     value=value, avg_cpue_x_100=avg_cpue_x_100, tot_catch= tot_catch, tot_effort = tot_effort, val_lake=val_lake,
     tot_val=tot_val, cost=cost, grav_wt=grav_wt)
}

#effort sequences
all_open <- rep(1,unit_time)                              #open all the time
week_2 <- rep(c(1,rep(0,5),1),length.out=unit_time)       #open weekends or two days a week
week_1 <- rep(c(1,rep(0,6)),length.out=unit_time)         #Open one weekend day a week
twoweek_2 <- rep(c(1,rep(0,12),1),length.out=unit_time)   #Open two weekend days every two weeks
twoweek_1 <- rep(c(1,rep(0,13)),length.out=unit_time)     #open 1 weekend day every two weeks
month_1 <- rep(c(1,rep(0,29)),length.out=unit_time)       #open 1 weekend day a month
week_1_1 <- rep(c(1,0,0,0,0,0,0),length.out=unit_time)
week_1_2 <- rep(c(0,1,0,0,0,0,0),length.out=unit_time)
week_1_3 <- rep(c(0,0,1,0,0,0,0),length.out=unit_time)
week_1_4 <- rep(c(0,0,0,1,0,0,0),length.out=unit_time)
week_1_5 <- rep(c(0,0,0,0,1,0,0),length.out=unit_time)
week_1_6 <- rep(c(0,0,0,0,0,1,0),length.out=unit_time)
week_1_7 <- rep(c(0,0,0,0,0,0,1),length.out=unit_time)


#lake schedules
allsame_open = matrix(c(all_open),nrow=length(day), ncol=nlakes)
allsame_week_2 = matrix(c(week_2),nrow=length(day), ncol=nlakes)
allsame_week1  = matrix(c(week_1),nrow=length(day), ncol=nlakes)
allsame_2week2 = matrix(c(twoweek_2),nrow=length(day),ncol=nlakes)
allsame_2week1 = matrix(c(twoweek_1),nrow=length(day), ncol=nlakes)
allsame_month1 = matrix(c(month_1), nrow=length(day), ncol=nlakes)
rotate = matrix(c(week_1_1, week_1_2, week_1_3, week_1_4, week_1_5, week_1_6, week_1_7), nrow=length(day), ncol=nlakes)
#oneeach   = matrix(c(all_open, week_2, week_1, twoweek_1), nrow=length(day), ncol=nlakes)
#threeopen = matrix(c(all_open, all_open, all_open, twoweek_1), nrow=length(day), ncol=nlakes)

#output
op1 = fun(allsame_open)
op2 = fun(allsame_week_2)
op3 = fun(allsame_week1)
op4 = fun(allsame_2week2)
op5 = fun(allsame_2week1)
op6 = fun(allsame_month1)
op7 = fun(rotate)


#checking val
op1$tot_val
op2$tot_val
op3$tot_val
op4$tot_val
op5$tot_val
op6$tot_val
op7$tot_val


colSums(op1$effort)
par(mfrow=c(2,2))
plot(day, op1$pmax_eff[,1], ylim=c(0,1), ylab="pmax_eff", xlab="day")
plot(rowSums(op1$vul),op1$pmax_eff[,1], ylab="pmax_eff", xlab="#vul")
plot(day,rowSums(op1$effort), ylab="effort", xlab="day")
plot(rowSums(op1$vul), rowSums(op1$effort), ylab="effort", xlab="vul#", ylim=c(0,1.2*max(rowSums(op1$effort))))

yy=max(op1$tot_val, op2$tot_val, op3$tot_val, op4$tot_val, op5$tot_val, op6$tot_val, op7$tot_val)

#checking val
op1$tot_val/yy
op2$tot_val/yy
op3$tot_val/yy
op4$tot_val/yy
op5$tot_val/yy
op6$tot_val/yy
op7$tot_val/yy
