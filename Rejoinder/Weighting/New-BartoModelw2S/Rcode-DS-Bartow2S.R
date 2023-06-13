set.seed(12345)
expit <- function(x){return((1+exp(-x))^(-1))}


nrep <- 1000
nsamp<- 2000
nburn <- 1000
nthin <- 5
nits<-nburn+nthin*nsamp

nvec <- c(200,500,1000,2000)

ests.mat.new<-array(0,c(length(nvec),nrep,nsamp))

estsATE.new <- array(NA,c(length(nvec),nrep))
coverATE.new <- array(NA,c(length(nvec),nrep))

for(i.ssize in 1:length(nvec)){
  n <- nvec[i.ssize]
  e0<-proc.time()[1]
  for(irep in 1:nrep){
    # data
   	X <- rnorm(n)
	Z <- rbinom(n,1,expit(X))
	X1<- exp(X)/5
	Y <- rpois(n,Z + X1)

    	#MCMC on Bartolucci model with 2S weights

	prec0<-prec1<-0.01

	old.lam0<-old.lam1<-2
	ind<-Z==1
	Y0<-Y[!ind];S0<-sum(Y0)
	Y1<-Y[ind];S1<-sum(Y1)
	n1<-sum(ind)
	n0<-n-n1

	fit0<-glm(Z~X,family=binomial)
	gam.est<-coef(fit0)

	old.gam0<-gam.est[1]
	old.gam1<-gam.est[2]
	old.pi<-expit(old.gam0+old.gam1*X)
	old.w<-Z/old.pi+(1-Z)/(1-old.pi)

	S0w<-sum(old.w[!ind]*Y0)
	S1w<-sum(old.w[ind]*Y1)
	w0<-sum(old.w[!ind])
	w1<-sum(old.w[ind])
	

	old.like0<-log(old.lam0)*S0w-w0*old.lam0
	old.prior0<--0.5*prec0*(log(old.lam0))^2-log(old.lam0)
	old.like1<-log(old.lam1)*S1w-w1*old.lam1
	old.prior1<--0.5*prec1*(log(old.lam1))^2-log(old.lam1)


	ico<-0
	u0<-log(runif(nits))
	u1<-log(runif(nits))
	z0<-rnorm(nits)*0.1
	z1<-rnorm(nits)*0.1

	for(iter in 1:nits){

		#Update lambdas given gammas

		new.lam0<-abs(old.lam0+z0[iter])
		new.like0<-log(new.lam0)*S0w-w0*new.lam0
		new.prior0<--0.5*prec0*(log(old.lam0))^2-log(old.lam0)
		if(u0[iter]<new.like0+new.prior0-old.like0-old.prior0){
			old.lam0<-new.lam0
			old.prior0<-new.prior0
			old.like0<-new.like0
		}
	
		new.lam1<-abs(old.lam1+z1[iter])
		new.like1<-log(new.lam1)*S1w-w1*new.lam1
		new.prior1<--0.5*prec1*(log(old.lam1))^2-log(old.lam1)
		if(u1[iter]<new.like1+new.prior1-old.like1-old.prior1){
			old.lam1<-new.lam1
			old.prior1<-new.prior1
			old.like1<-new.like1
		}

		if(iter > nburn & iter %% nthin == 0){
			ico<-ico+1
			ests.mat.new[i.ssize,irep,ico]<-old.lam1-old.lam0
		}
	}
	e1<-proc.time()[1]
	if(irep %% 10 == 0){
		print(c(n,irep,as.numeric(e1-e0)))
	}
  }
  save.image("new.RData")
}

#quit(save="yes")

hist(apply(ests.mat.new[1,,],1,mean))

ests.ci<-apply(ests.mat.new,1:2,quantile,prob=c(0.025,0.975))

ests.cover<-(ests.ci[1,,]<1 & ests.ci[2,,] > 1)

apply(ests.cover,1,mean)

boxplot(t(apply(ests.mat.new,1:2,mean)))

boxplot(t(ests.mat.new[4,1:50,]))
