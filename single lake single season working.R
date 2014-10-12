### Simulation code to represent vulnerable pool dynamics from CR vul exchange for ed working spreadsheet
### Looks at a single lake for a year, assuming fixed biologial and fisheries parms, across different closure regimes

rm(list=ls())

### Parameter Section ###
unit_time <- 360        #days in a year
No <- 500              #Initial abundance
M  <- 0.5/unit_time     #natural mortality
q  <- 0.05              #catchability
v1 <- 2/unit_time       #exchange rate to the vulnerable pool
v2 <- 3/unit_time       #exchange rate to the invulnerable pool
rel_prop <- 1           #proportion released
rel_surv <- 0.9         #release survival rate
recov <- 0.2            #rate of recovery from the refractory pool
pvul_recov <- 0.5       #Propotion of recovered fish that go to the vulnerable pool
sd_obs <- .2            #standard deviation of observation error of fishers, such that greater values lead to greater divergence from perfect knowledge

cpue_base <- 0.5          #catch rate required for a satisfaction or utility of 1
val_power <- 1.5          #power function of how satisaction increases with catch rate
#scrap code for picturing alternative value functions
#y=seq(0,10, length=100)
#sat=(y/cpue_base)^val_power
#plot(y,sat, type="l", xlab="cpue",ylab="sat")

max_eff_wknd <- 10      #max effort possible on a weekend
max_eff_wkdy <- 10      #max effort possible on a weekday

cpue_half <- 2                #inflection point for logistic effort
sdcpue <- 0.5                   #sigma of logistic effort

theta<-list()
theta$No <- No
theta$M <- M
theta$q <- q
theta$v1 <- v1
theta$v2 <- v2
theta$rel_surv <- rel_surv
theta$recov <- recov
theta$pvul_recov <- pvul_recov
theta$cpue_base <- cpue_base          #catch rate required for a satisfaction or utility of 1
theta$val_power <- val_power          #power function of how satisaction increases with catch rate
theta$cpue_half <- cpue_half                #inflection point for logistic effort
theta$sdcpue <- sdcpue                   #sigma of logistic effort
theta$max_eff_wknd <- max_eff_wknd      #max effort possible on a weekend
theta$max_eff_wkdy <- max_eff_wkdy      #max effort possible on a weekday

#effort sequences
all_open <- rep(1,unit_time)                              #open all the time
week_2 <- rep(c(1,rep(0,5),1),length.out=unit_time)       #open weekends or two days a week
week_1 <- rep(c(1,rep(0,6)),length.out=unit_time)         #Open one weekend day a week
twoweek_2 <- rep(c(1,rep(0,12),1),length.out=unit_time)   #Open two weekend days every two weeks
twoweek_1 <- rep(c(1,rep(0,13)),length.out=unit_time)     #open 1 weekend day every two weeks
month_1 <- rep(c(1,rep(0,29)),length.out=unit_time)       #open 1 weekend day a month



### Procedure Section ###
fun <- function(theta,open_seq){
  No <- theta$No
  M <- theta$M
  q <- theta$q
  v1 <- theta$v1
  v2 <- theta$v2
  rel_surv <- theta$rel_surv
  recov <- theta$recov
  pvul_recov <- theta$pvul_recov
  cpue_base <- theta$cpue_base
  val_power <- theta$val_power
  cpue_half <- theta$cpue_half
  sdcpue <- theta$sdcpue
  max_eff_wknd <- theta$max_eff_wknd
  max_eff_wkdy <- theta$max_eff_wkdy
  surv <- exp(-M)         #survival
    
  #vectors
  day <- seq(1,unit_time, by=1)
  vul=NULL; invul = NULL; refract = NULL; max_eff = NULL; pmax_eff = NULL;
  effort = NULL; catch = NULL; vul_surv = NULL; invul_surv = NULL; refract_surv = NULL;
  alive = NULL; eff_cum = NULL; cpue = NULL; value = NULL; obs_err =NULL;
#  open_seq=NULL;
  
  #independent vectors
#  open_seq = week_2;
  vul_init <- v1/(v1+v2)*No     #initial number of vulnerable fish
  max_cpue <- q*vul_init        #max cpue, don't think this is used

  max_eff = rep(c(max_eff_wknd, rep(max_eff_wkdy,5),max_eff_wknd),length.out=unit_time)   #gives the max effort for any given day, assuming potential difference between weekends and weekdays
  
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
  value[1] = effort[1]*max(0,cpue[1]/cpue_base-1)^val_power
  
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
      value[i] = effort[i]*max(0,cpue[i]/cpue_base-1)^val_power
  }
  avg_cpue <- mean(cpue[which(effort>0)])                                          #average cpue where effort>0,
  avg_cpue_x_100 <- (mean(cpue[which(effort>0)]))*100                              #*100 to make plotting easier
  tot_catch <- sum(catch);   tot_effort <- sum(effort);   tot_val <- sum(value);   #response metrics we might care about
  metrics <- as.matrix(data.frame(avg_cpue_x_100, tot_catch, tot_effort, tot_val)) #dumb wrappers to make barplots easy

  list(day=day, vul=vul, invul=invul, refract=refract, pmax_eff=pmax_eff, effort=effort, catch=catch,  
     vul_surv=vul_surv, invul_surv=invul_surv, refract_surv=refract_surv, alive=alive, eff_cum=eff_cum, cpue=cpue,
     value=value, avg_cpue_x_100=avg_cpue_x_100, tot_catch= tot_catch, tot_effort = tot_effort, tot_val=tot_val,
     avg_cpue=avg_cpue, metrics=metrics)
}


### Report Section ###
#outputs of different closure schedules
op1 = fun(theta, open_seq=all_open)
op2 = fun(theta, open_seq=week_2)
op3 = fun(theta, open_seq=week_1)
op4 = fun(theta, open_seq=twoweek_2)
op5 = fun(theta, open_seq=twoweek_1)
op6 = fun(theta, open_seq=month_1)


#Plots of numbers of fish in different pools over time
par(mfrow=c(2,3), mai=c(.3,.2,.3,.2), omi=c(.6,.7,.1,.1), mgp=c(1,.4,0))
#plot1, all days open to fishing
plot(op1$day, op1$vul, type="l", col="black", lwd=2, ylim=c(0,800), xlim=c(0,150), xlab="", ylab="", xaxt="n", main="All Open")
lines(op1$day, op1$invul, type="l", col="blue", lwd=2)
lines(op1$day, op1$refract, type="l", col="green", lwd=2)
legend('topright', legend='Invulnerable', col='blue', bty="n", lty=1, lwd=2)
#plot 2, Open 2 days a week
plot(op2$day, op2$vul, type="l", col="black", lwd=2, ylim=c(0,800), xlim=c(0,150), xlab="", ylab="", main="2 days/week", yaxt="n", xaxt="n")
lines(op2$day, op2$invul, type="l", col="blue", lwd=2)
lines(op2$day, op2$refract, type="l", col="green", lwd=2)
legend('topright', legend='Vulnerable', col='black', bty="n", lty=1, lwd=2)
#plot 3, open 1 day a week
plot(op3$day, op3$vul, type="l", col="black", lwd=2, ylim=c(0,800), xlim=c(0,150), xlab="", ylab="", main="1 day/week", yaxt="n", xaxt="n")
lines(op3$day, op3$invul, type="l", col="blue", lwd=2)
lines(op3$day, op3$refract, type="l", col="green", lwd=2)
legend('topright', legend='Refractory', col='green', bty="n", lty=1, lwd=2)
#plot 4, open 2 consecutive days every two weeks
plot(op4$day, op4$vul, type="l", col="black", lwd=2, ylim=c(0,800), xlim=c(0,150), xlab="", ylab="", main="2 days/2 weeks")
lines(op4$day, op4$invul, type="l", col="blue", lwd=2)
lines(op4$day, op4$refract, type="l", col="green", lwd=2)
mtext("Numbers of fish", side=2, line=3, adj=2.2, cex=1.5, outer=F)
#plot 5, open 1 day every two weeks
plot(op5$day, op5$vul, type="l", col="black", lwd=2, ylim=c(0,800), xlim=c(0,150), xlab="", ylab="", main="1 day/2 weeks", yaxt="n")
lines(op5$day, op5$invul, type="l", col="blue", lwd=2)
lines(op5$day, op5$refract, type="l", col="green", lwd=2)
mtext("Days", side=1, line=3, cex=1.5, outer=F)
#plot 6, open 1 day every two weeks
plot(op6$day, op6$vul, type="l", col="black", lwd=2, ylim=c(0,800), xlim=c(0,150), xlab="", ylab="", main="1 day/month", yaxt="n")
lines(op6$day, op6$invul, type="l", col="blue", lwd=2)
lines(op6$day, op6$refract, type="l", col="green", lwd=2)


#plots of effort over time
par(mfrow=c(2,3), mai=c(.3,.2,.3,.2), omi=c(.6,.7,.1,.1), mgp=c(1,.4,0))
plot(op1$day, op1$effort, type="l", col="black", lwd=2, ylim=c(0,10), xlim=c(0,300), xlab="", ylab="", xaxt="n", main="All Open")
plot(op2$day, op2$effort, type="l", col="black", lwd=2, ylim=c(0,10), xlim=c(0,300), xlab="", ylab="", xaxt="n", yaxt ="n", main="2 days/week")
plot(op3$day, op3$effort, type="l", col="black", lwd=2, ylim=c(0,10), xlim=c(0,300), xlab="", ylab="", xaxt="n", yaxt ="n", main="1 day/week")
plot(op4$day, op4$effort, type="l", col="black", lwd=2, ylim=c(0,10), xlim=c(0,300), xlab="", ylab="",  main="2 days/2 weeks")
mtext("Effort (trips/day)", side=2, line=3, adj=2.2, cex=1.5, outer=F)
plot(op5$day, op5$effort, type="l", col="black", lwd=2, ylim=c(0,10), xlim=c(0,300), xlab="", ylab="",  yaxt ="n", main="1 day/2 weeks")
mtext("Days", side=1, line=3, cex=1.5, outer=F)
plot(op6$day, op6$effort, type="l", col="black", lwd=2, ylim=c(0,10), xlim=c(0,300), xlab="", ylab="",  yaxt ="n", main="1 day/month")


#Compare 1 day/wk with all open
"Compare.plots" <- function()
{
  par(mar=c(1,1,2,1),oma=c(4,4,1,4))
  layout(matrix(c(1:6),3,2,byrow=T))
  plot(op1$day, op1$vul, type="l", lty=1, lwd=2, ylim=c(0,400), xlim=c(0,150), xlab="", ylab="", xaxt="n", yaxt="n",main="All Open",cex.main=1.5)
  lines(op1$day, op1$invul, lty=3, lwd=2)
  lines(op1$day, op1$refract, lwd=2,col="grey")
  legend("center", legend=c("Invulnerable","Vulnerable","Refractory"), bty="n", lty=c(3,1,1),col=c(rep(1,2),"grey"), lwd=2,cex=1.3)
  axis(2,at=seq(0,400,200),labels=TRUE,cex.axis=1.3)
  mtext(expression(paste("Density (Fish "^{.}," ha"^{-1},")")),2,cex=1.2,line=2.5)
  plot(op3$day, op3$vul, type="l", lty=1, lwd=2, ylim=c(0,400), xlim=c(0,150), xlab="", ylab="", main="1 day/week", yaxt="n", xaxt="n",cex.main=1.5)
  lines(op3$day, op3$invul, lty=3, lwd=2)
  lines(op3$day, op3$refract, col="grey", lwd=2)
  plot(op1$day, op1$effort, type="l",lty=1,lwd=2,xlim=c(0,150),xlab="",ylab="",cex.axis=1.3,ylim=c(0,10),xaxt="n")
  lines(op1$day,op1$cpue,lty=2,lwd=2)
  legend("top",legend=c("Effort","CPUE"),bty="n",lty=1:2,lwd=2,cex=1.3)
  mtext("Effort (AD)",2,cex=1.2,line=2.5)
  x.tmp<-which(op3$effort>0)
  plot(op3$day[x.tmp], op3$effort[x.tmp], pch=19,lwd=2,xlim=c(0,150),xlab="",ylab="",yaxt="n",cex=0.8,ylim=c(0,10),xaxt="n")
  points(op3$day[x.tmp],op3$cpue[x.tmp],lty=2,xlim=c(0,150),cex=0.8)
  legend("center",legend=c("Effort","CPUE"),bty="n",pch=c(19,1),cex=1.3,pt.cex=1)
  axis(4,at=seq(0,10,2),seq(0,10,2),cex.axis=1.3)
  mtext(expression(paste("CPUE (Fish "^{.}," AD"^{-1},")")),4,cex=1.2,line=3)
  plot(op1$day,op1$value,type="l",lwd=2,xlab="",ylab="",xlim=c(0,150),ylim=c(0,600),yaxt="n",cex.axis=1.3)
  axis(2,at=seq(0,600,200),labels=TRUE,cex.axis=1.3)
  mtext("Satisfaction",2,cex=1.2,line=2.5)
  plot(op3$day[x.tmp],op3$value[x.tmp],pch=19,xlab="",ylab="",yaxt="n",cex=0.8,xlim=c(0,150),ylim=c(0,600),cex.axis=1.3)
  par(mfcol=c(1,1))
  mtext("Time (d)",1,cex=1.2,line=3.5)
}

#plots response metrics
"Response.metrics" <- function()
{
  par(mar=c(2,2.5,2,2.5),oma=c(2,2,1,0))
  layout(matrix(c(1:4),2,2,byrow=T))
  barplot(c(op1$tot_catch,op2$tot_catch,op3$tot_catch,op4$tot_catch,op5$tot_catch,op6$tot_catch),main="Total catch",names.arg=1:6)
  mtext("Total catch",2,line=3.5)
  mtext(expression(paste("(Fish "^{.}," ha"^{"-1 ."}," y"^{-1},")")),2,line=2.2,cex=0.9)
  barplot(c(op1$tot_effort,op2$tot_effort,op3$tot_effort,op4$tot_effort,op5$tot_effort,op6$tot_effort),main="Total effort",names.arg=1:6)
  mtext("Total effort",2,line=3.5)
  mtext(expression(paste("(AD "^{.}," ha"^{"-1 ."}," y"^{-1},")")),2,line=2.2,cex=0.9)
  legend("topright",legend=c("All open","2 days/week","1 day/week","2 days/2 weeks","1 day/2 weeks","1 day/month"),bty="n",lty=0,pch=49:54,pt.cex=1.3)#,title="Scenarios")
  barplot(c(op1$avg_cpue,op2$avg_cpue,op3$avg_cpue,op4$avg_cpue,op5$avg_cpue,op6$avg_cpue),main="Catch-per-unit effort",names.arg=1:6)
  mtext("Catch-per-unit-effort",2,line=3.5)
  mtext(expression(paste("(Fish "^{.}," d"^{-1},")")),2,line=2.2,cex=0.9)
  barplot(c(op1$tot_val,op2$tot_val,op3$tot_val,op4$tot_val,op5$tot_val,op6$tot_val),main="Value",names.arg=1:6)
  mtext("Value (AD)",2,line=2.2)
  par(mfcol=c(1,1))
  mtext("Scenario",1,font=2,cex=1,line=2.5)
}

"No.Exchange" <- function()
{
  op1 = fun(theta, open_seq=all_open)
  op3 = fun(theta, open_seq=week_1)
  
  resp <- matrix(nrow=4,ncol=4)
  resp[,3]<-c(op1$tot_catch,op1$tot_effort,op1$avg_cpue,op1$tot_val)
  resp[,4]<-c(op3$tot_catch,op3$tot_effort,op3$avg_cpue,op3$tot_val)
  theta2<-theta
  theta2$v1<-1
  theta2$v2<-0
  theta2$pvul_recov<-1
  theta2$recov<-0.2
  op1b = fun(theta2, open_seq=all_open)
  op3b = fun(theta2, open_seq=week_1)
  resp[,1]<-c(op1b$tot_catch,op1b$tot_effort,op1b$avg_cpue,op1b$tot_val)
  resp[,2]<-c(op3b$tot_catch,op3b$tot_effort,op3b$avg_cpue,op3b$tot_val)
  return(resp)
}

"Satis.exponent" <- function()
{
  op1 = fun(theta, open_seq=all_open)
  op3 = fun(theta, open_seq=week_1)
  
  resp <- matrix(nrow=4,ncol=4)
  resp[,3]<-c(op1$tot_catch,op1$tot_effort,op1$avg_cpue,op1$tot_val)
  resp[,4]<-c(op3$tot_catch,op3$tot_effort,op3$avg_cpue,op3$tot_val)
  theta2<-theta
  theta2$v1<-1
  theta2$v2<-0
  theta2$pvul_recov<-1
  theta2$recov<-0.2
  op1b = fun(theta2, open_seq=all_open)
  op3b = fun(theta2, open_seq=week_1)
  resp[,1]<-c(op1b$tot_catch,op1b$tot_effort,op1b$avg_cpue,op1b$tot_val)
  resp[,2]<-c(op3b$tot_catch,op3b$tot_effort,op3b$avg_cpue,op3b$tot_val)
  return(resp)
}

"Response.to.beta" <- function()
{
  names=list()
  out.names[[1]]<-c(expression(paste(beta,"=0.5")),expression(paste(beta,"=1.0")),expression(paste(beta,"=1.5")))
  out.names[[2]]<-c("Catch","Effort","CPUE","Value")
  val.tab<-matrix(nrows=3,ncol=4,dimnames=out.names)
  
  par(mar=c(2,2.5,2,2.5),oma=c(2,2,1,0))
  layout(matrix(c(1:3),1,3,byrow=T))
  legend("topright",legend=c("All open","2 days/week","1 day/week","2 days/2 weeks","1 day/2 weeks","1 day/month"),bty="n",lty=0,pch=49:54,pt.cex=1.3)#,title="Scenarios")
  barplot(c(op1$avg_cpue,op2$avg_cpue,op3$avg_cpue,op4$avg_cpue,op5$avg_cpue,op6$avg_cpue),main="Catch-per-unit effort",names.arg=1:6)
  mtext("Catch-per-unit-effort",2,line=3.5)
  mtext(expression(paste("(Fish "^{.}," d"^{-1},")")),2,line=2.2,cex=0.9)
  barplot(c(op1$tot_val,op2$tot_val,op3$tot_val,op4$tot_val,op5$tot_val,op6$tot_val),main="Value",names.arg=1:6)
  mtext("Value (AD)",2,line=2.2)
  par(mfcol=c(1,1))
  mtext("Scenario",1,font=2,cex=1,line=2.5)
}

"Sens.Analysis" <- function()
{
  l<-length(theta)
  elast<-array(dim=c(6,l,2))
  theta$val_pow=1.0
  for( i in 1:l)
  {
    op1 = fun(theta, open_seq=all_open)$tot_val
    op2 = fun(theta, open_seq=week_2)$tot_val
    op3 = fun(theta, open_seq=week_1)$tot_val
    op4 = fun(theta, open_seq=twoweek_2)$tot_val
    op5 = fun(theta, open_seq=twoweek_1)$tot_val
    op6 = fun(theta, open_seq=month_1)$tot_val
    theta2<-theta
    theta2[i]<-theta[[i]]*0.9
    elast[1,i,1] = (fun(theta2, open_seq=all_open)$tot_val-op1)/op1
    elast[2,i,1] = (fun(theta2, open_seq=week_2)$tot_val-op2)/op2
    elast[3,i,1] = (fun(theta2, open_seq=week_1)$tot_val-op3)/op3
    elast[4,i,1] = (fun(theta2, open_seq=twoweek_2)$tot_val-op4)/op4
    elast[5,i,1] = (fun(theta2, open_seq=twoweek_1)$tot_val-op5)/op5
    elast[6,i,1] = (fun(theta2, open_seq=month_1)$tot_val-op6)/op6
    theta2[i]<-theta[[i]]*1.1
    elast[1,i,2] = (fun(theta2, open_seq=all_open)$tot_val-op1)/op1
    elast[2,i,2] = (fun(theta2, open_seq=week_2)$tot_val-op2)/op2
    elast[3,i,2] = (fun(theta2, open_seq=week_1)$tot_val-op3)/op3
    elast[4,i,2] = (fun(theta2, open_seq=twoweek_2)$tot_val-op4)/op4
    elast[5,i,2] = (fun(theta2, open_seq=twoweek_1)$tot_val-op5)/op5
    elast[6,i,2] = (fun(theta2, open_seq=month_1)$tot_val-op6)/op6
  }
  theta$val_pow=1.5
  plot(1:(l-1),elast[3,1:(l-1),1],pch=6,ylab="",xlab="",xaxt="n",yaxt="n",ylim=c(-0.3,0.3))
  points(1:(l-1),elast[3,1:(l-1),2],pch=2)
  for(j in 1:(l-1))
    lines(rep(j,2),elast[3,j,1:2],lty=3)
  axis(2,at=seq(-0.3,0.45,0.15),labels=seq(-0.3,0.45,0.15))
  abline(h=0,lty=2,lwd=2)
  xlabs<-c(expression(paste(italic("N"[0]))),
           expression(italic("M")),
           expression(italic("q")),
           expression(paste(italic("v"[1]))),
           expression(paste(italic("v"[2]))),
           expression(paste(italic("S"[r]))),
           expression(paste(italic("v"[r]))),
           expression(paste(italic("p"[v]))),
           expression(paste(italic(C[0]))),
           expression(italic(beta)),
           expression(paste(italic("C"[50]))),
           expression(paste(italic("C"[sigma]))),
           expression(paste(italic("E"["WE"]))),
           expression(paste(italic("E"["WD"]))))
  axis(1,at=1:(l-1),labels=c(xlabs[1:(l-2)],expression(paste(italic("E"["max"])))))
  mtext("Elasticity",2,line=2.5,cex=1.5)
  mtext("Parameter",1,line=2.5,cex=1.5)
}

"Value.plot" <- function()
{
  plot(seq(0,5,0.001),(seq(0,5,0.001)/cpue_base-1)^1,t="l",ylim=c(0,5),xaxs="i",yaxs="i",xlab="",ylab="")
  lines(seq(0,5,0.001),(seq(0,5,0.001)/cpue_base-1)^1.5,lty=2)
  lines(seq(0,5,0.001),(seq(0,5,0.001)/cpue_base-1)^0.5,lty=3)
  arrows(x0=cpue_base,y0=2,x1=cpue_base,y1=0.5,code=2,length=0.1,lwd=1.5)
  text(cpue_base,2.25,labels=expression(italic(paste("C"[0]))),cex=1.4)
  legend("topleft",legend=c(expression(paste(beta,"=0.5")),
                            expression(paste(beta,"=1.0")),
                            expression(paste(beta,"=1.5"))),lty=c(3,1,2),bty="n",cex=1.4)
  mtext(expression(bold("Value (AD)")),2,cex=1.25,line=2.3)
  mtext(expression(bold(paste("CPUE (fish "^{.}," AD"^{-1},")"))),1,cex=1.25,line=2.75)
}

"Sens.to.E" <- function()
{
  theta2<-theta
  val.vec <- matrix(nrow=100,ncol=2)
  for(i in 1:100)
  {
    theta2$max_eff_wknd<-i
    theta2max_eff_wkdy<-i
    val.vec[i,1] = fun(theta2, open_seq=all_open)$tot_val
    val.vec[i,2]= fun(theta2, open_seq=week_1)$tot_val
  }
  plot(val.vec[,2],t="l",yaxt="n",ylab="",xlab="",ylim=c(1000,4000),lwd=2,cex.axis=1.2)
  lines(val.vec[,1],col="grey",lwd=2)
  axis(2,at=seq(1000,4000,500),labels=T,cex.axis=1.2)
  mtext("Maximum daily effort (AD)",1,cex=1.5,line=2.5)
  mtext("Value (AD)",2,cex=1.5,line=3)
  legend(x=45,y=3500,legend=c("All open","1 day/week"),col=c("grey","black"),lwd=2,bty="n",lty=1,cex=1.2)
}
