#Compute ATE by BB

library(fastglm)

set.seed(23498)
nvec <- c(200,500,1000,2000)
nreps<-1000
L<-2000
expit<-function(x){1/(1+exp(-X))}

est.mat<-array(0,c(length(nvec),nreps,L))

for(ni in 1:length(nvec)){
	n<-nvec[ni]
	e0<-proc.time()[1]
	for(irep in 1:nreps){
		X<- rnorm(n)
    	    	Z <- rbinom(n,1,expit(X))
          	Y <- rpois(n,Z + exp(X)/5)
		wmat<-matrix(rgamma(L*n,1),nrow=L,ncol=n)
		iperm<-sample(1:L)

		#Unlinked 
		ipw.est<-rep(L)
		for(l in 1:L){
			w<-wmat[l,]
			pi.hat<-fitted(fastglm(y=Z,x=cbind(1,X),weights=w,family=quasibinomial))
			w0<-(1-Z)/(1-pi.hat)
			w1<-Z/pi.hat
			w<-wmat[iperm[l],]
			ipw.est[l]<-sum(w1*w*Y)/sum(w1*w)-sum(w0*w*Y)/sum(w0*w)
		}
		est.mat[ni,irep,]<-ipw.est	
		e1<-proc.time()[1]
		print(c(n,irep,as.numeric(e1-e0)))
	}
}

hist(apply(est.mat[1,,],1,mean))

ests.ci<-apply(est.mat,1:2,quantile,prob=c(0.025,0.975))

ests.cover<-(ests.ci[1,,]<1 & ests.ci[2,,] > 1)

apply(ests.cover,1,mean)

boxplot(t(apply(est.mat,1:2,mean)))
boxplot(t(apply(est.mat,1:2,var)))

boxplot(t(est.mat[4,1:50,]))


