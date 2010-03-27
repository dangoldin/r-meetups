library(rjags)

TrainingFraction<- 1.0

MyCounter <- function(vector) {
  sum(!is.na(vector))
}

line.inits1 <- list(tau=1,beta=c(0,0,1,-0.02),beta3=c(1,1,1,1,1,1,1,1),beta4=c(-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,-0.02),tauBeta1=1,tauBeta2=1,tauBeta3=1,tauBeta4=1,beta3int=0,beta3slope=1,beta4int=0,beta4slope=1)
#list(tau=1,beta=c(0,0),beta3=c(1,1,1,1,1,1,1),beta4=c(-0.02,-0.02,-0.02,-0.02,-0.02,-0.02,-0.02),tauBeta2=1,tauBeta3=1,tauBeta4=1,beta3int=0,beta3slope=1,beta4int=0,beta4slope=1)

y <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsIgg.txt",header=F,sep="\t"))
time <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsTime.txt",header=F,sep="\t"))
timesincedose <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose.txt",header=F,sep="\t"))
timesincedose1 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose1.txt",header=F,sep="\t"))
timesincedose2 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose2.txt",header=F,sep="\t"))
timesincedose3 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose3.txt",header=F,sep="\t"))
timesincedose4 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose4.txt",header=F,sep="\t"))
timesincedose5 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose5.txt",header=F,sep="\t"))
timesincedose6 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose6.txt",header=F,sep="\t"))
timesincedose7 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose7.txt",header=F,sep="\t"))
timesincedose8 <- as.matrix(read.table("/Users/dbm/cdc/extractTimeHumanBugsDose8.txt",header=F,sep="\t"))

# remove last columns which is junk
y <- y[,-dim(y)[2]]
time <- time[,-dim(time)[2]]
timesincedose <- timesincedose[,-dim(timesincedose)[2]]
timesincedose1 <- timesincedose1[,-dim(timesincedose1)[2]]
timesincedose2 <- timesincedose2[,-dim(timesincedose2)[2]]
timesincedose3 <- timesincedose3[,-dim(timesincedose3)[2]]
timesincedose4 <- timesincedose4[,-dim(timesincedose4)[2]]
timesincedose5 <- timesincedose5[,-dim(timesincedose5)[2]]
timesincedose6 <- timesincedose6[,-dim(timesincedose6)[2]]
timesincedose7 <- timesincedose7[,-dim(timesincedose7)[2]]
timesincedose8 <- timesincedose8[,-dim(timesincedose8)[2]]



# this probably doesn't make sense - better to forecast within subject
RandomNumbers <- runif(dim(y)[1])
train <- (RandomNumbers <= TrainingFraction)
test <- (RandomNumbers > TrainingFraction)


# this computes the number of IgG measurements that each subject has
NumberOfMeasurements <- apply(y[train,],1,MyCounter)


# the next block of code zaps the last measurement for each subject 
# so that the measurement can be forecast
#ActualSecondLast <- ActualLast <- rep(NA,1303)
#for (i in 1:1303) {
#	ActualLast[i] <- y[i,NumberOfMeasurements[i]]
#    if (NumberOfMeasurements[i] > 1) {
#  	  ActualSecondLast[i] <- y[i,NumberOfMeasurements[i]-1]
#  	}
#	y[i,NumberOfMeasurements[i]] <- NA
#}



line.data <- list(N=dim(y[train,])[1],T=NumberOfMeasurements,time=time[train,],timesincedose=timesincedose[train,],timesincedose1=timesincedose1[train,],timesincedose2=timesincedose2[train,],timesincedose3=timesincedose3[train,],timesincedose4=timesincedose4[train,],timesincedose5=timesincedose5[train,],timesincedose6=timesincedose6[train,],timesincedose7=timesincedose7[train,],timesincedose8=timesincedose8[train,],y=y[train,],zeroes=c(0,0,0))
m <- jags.model("/Users/dbm/cdc/JAGS/line.bugRandomSlope", 
                data=line.data,
		inits=line.inits1,
                n.chains=1
               )
update(m,1000)   # burn-in

#dic <- dic.samples(m,n.iter=10000,type="popt")
#print(dic)

parameters <- c("beta","individual.beta2","beta3","beta4")

# the next setting for parameters allows one to track the elisaLast
#parameters <- NULL
#for (i in 1:1303) {
#   parameters <- c(parameters,paste("y[",i,",",NumberOfMeasurements[i],"]",sep=""))
#}

samples <- coda.samples(m, parameters, 10000)

xxx <- apply(samples[[1]],2,mean)    # these are in alpha order 1,10,100,1000,1001,...

# RMSE calculations
alphaNumbers <- read.table("/Users/dbm/cdc/alphaNumbers.txt")
sqrt( sum((ActualLast[alphaNumbers[,1]] - xxx)^2)/1303 )

sqrt( sum((ActualLast - ActualSecondLast)^2,na.rm=TRUE)/1303 )


#summary(samples)
#plot(samples[[1]])
beta1 <- xxx[1]
beta3 <- xxx[5:12]
beta4 <- xxx[13:20]
beta2 <- xxx[21:1323]



# computation of SSE
SSE <- 0
for (i in (1:1303)[test]) {
  temp <- beta1+(beta2[i]*log(time[i,]+1)) + (beta3[2]*exp(beta4[2]*timesincedose2[i,])) + (beta3[5]*exp(beta4[5]*timesincedose5[i,])) + (beta3[6]*exp(beta4[6]*timesincedose6[i,])) + (beta3[7]*exp(beta4[7]*timesincedose7[i,])) + (beta3[1]*exp(beta4[1]*timesincedose1[i,])) + (beta3[3]*exp(beta4[3]*timesincedose3[i,])) + (beta3[4]*exp(beta4[4]*timesincedose4[i,])) + (beta3[8]*exp(beta4[8]*timesincedose8[i,]))
	SSE <- SSE + sum((y[i,] - temp)^2,na.rm=TRUE)
}

# computation of lagged predictions
lag <- 52
PredictedLast <- rep(0,1303)
temp <- matrix(0,1303,17)
for (i in 1:1303) {
  temp[i,] <- beta1+(beta2[i]*log(lag+time[i,]+1)) + (beta3[1]*exp(beta4[1]*(lag+timesincedose1[i,]))) + (beta3[2]*exp(beta4[2]*(lag+timesincedose2[i,]))) + (beta3[3]*exp(beta4[3]*(lag+timesincedose3[i,]))) + (beta3[4]*exp(beta4[4]*(lag+timesincedose4[i,]))) + (beta3[5]*exp(beta4[5]*(lag+timesincedose5[i,]))) + (beta3[6]*exp(beta4[6]*(lag+timesincedose6[i,]))) + (beta3[7]*exp(beta4[7]*(lag+timesincedose7[i,]))) + (beta3[8]*exp(beta4[8]*(lag+timesincedose8[i,])))

  PredictedLast[i] <- temp[i,NumberOfMeasurements[i]]
}

# next few lines are when doing group 4 only
#beta2 <- xxx[21:81]
#beta3[2] <- 0
#beta3[5] <- 0
#beta3[6] <- 0
#beta3[7] <- 0

# x <- jags.samples(m,parameters,n.iter=1000)

pdf(paste("/Users/dbm/cdc/JAGS/tracesGroup4COP.pdf"),height=20,width=20)
par(mfrow=c(4,2))
par(mai=c(0.5,0.5,0.25,0.25))
for (i in sample(1:1303,8)) {
  plot(time[i,],y[i,],ylim=c(0,max(y[i,],na.rm=T)),type="l",xlab="",ylab="")
  par(new=TRUE)
  plot(time[i,],beta1+(beta2[i]*log(time[i,]+1))
#  + (beta3[2]*exp(beta4[2]*timesincedose2[i,])) 
#  + (beta3[5]*exp(beta4[5]*timesincedose5[i,])) 
#  + (beta3[6]*exp(beta4[6]*timesincedose6[i,])) 
#  + (beta3[7]*exp(beta4[7]*timesincedose7[i,])) 
  + (beta3[1]*exp(beta4[1]*timesincedose1[i,])) 
  + (beta3[3]*exp(beta4[3]*timesincedose3[i,])) 
  + (beta3[4]*exp(beta4[4]*timesincedose4[i,])) 
  + (beta3[8]*exp(beta4[8]*timesincedose8[i,])), 
  ylim=c(0,max(y[i,],na.rm=T)),type="l",col="red",xlab="time (weeks)",ylab="IgG")
       
  legend(20,1,c("actual","fitted"),col=c(1,2),lty=c(1,1))
}                      
dev.off()

    

# coda.samples            Generate posterior samples in mcmc.list format
# jags.model              Create a JAGS model object
# jags.module             Dynamically load JAGS modules
# jags.samples            Generate posterior samples
# print.mcarray           Objects for representing MCMC output
# update.jags             Functions for manipulating jags model objects

# as.mcmc.list
# coef(m)
