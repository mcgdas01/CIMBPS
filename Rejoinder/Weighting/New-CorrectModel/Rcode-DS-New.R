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

    # MCMC on correct model

	old.th<-c(1,1)
	prec0<-0.01
	old.eta<-old.th[1]*Z+old.th[2]*X1
	old.like<-sum(log(old.eta)*Y-old.eta)
	old.prior<--0.5*prec0*sum(old.th^2)
	ico<-0
	umat<-log(matrix(runif(nits*2),nrow=nits,ncol=2))
	zmat<-matrix(rnorm(nits*2),nrow=nits,ncol=2)
	for(iter in 1:nits){
		new.th<-old.th
		for(j in 1:2){			
			new.th[j]<-old.th[j]+zmat[iter,j]
			new.eta<-new.th[1]*Z+new.th[2]*X1
			if(min(new.eta)<0){
				new.like<--1e10				
			}else{
				new.like<-sum(log(new.eta)*Y-new.eta)
			}
			new.prior<--0.5*prec0*sum(new.th^2)

			if(umat[iter,j] < new.like+new.prior-old.like-old.prior){
				old.th[j]<-new.th[j]
				old.like<-new.like
				old.prior<-new.prior
			}
		}
		#print(old.th)
		if(iter > nburn & iter %% nthin == 0){
			ico<-ico+1
			ests.mat.new[i.ssize,irep,ico]<-old.th[1]
		}
	}
	e1<-proc.time()[1]
	if(irep %% 10 == 0){
		print(c(n,irep,as.numeric(e1-e0)))
	}
  }
  save.image("correct.RData")
}

#quit(save="yes")
