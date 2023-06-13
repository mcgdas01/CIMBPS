#Compute ATE by BB
set.seed(23498)
nvec <- c(200,500,1000,2000)
nreps<-1000
L<-2000
expit<-function(x){1/(1+exp(-X))}

est.mat<-array(0,c(length(nvec),nreps,L))

for(ni in 1:length(nvec)){
	n<-nvec[ni]
	for(irep in 1:nreps){
		X<- rnorm(n)
    	    	Z <- rbinom(n,1,expit(X))
          	Y <- rpois(n,Z + exp(X)/5)
		wmat<-matrix(rgamma(L*n,1),nrow=L,ncol=n)

		#Linked
		ipw.est<-rep(L)
		for(l in 1:L){
			w<-wmat[l,]
			pi.hat<-fitted(glm(Z~X,weights=w,family=binomial))
			w0<-(1-Z)/(1-pi.hat)
			w1<-Z/pi.hat
			ipw.est[l]<-sum(w1*w*Y)/sum(w1*w)-sum(w0*w*Y)/sum(w0*w)
		}
		est.mat[ni,irep,]<-ipw.est	
		print(c(n,irep))
	}
}

hist(apply(est.mat[1,,],1,mean))

ests.ci<-apply(est.mat,1:2,quantile,prob=c(0.025,0.975))

ests.cover<-(ests.ci[1,,]<1 & ests.ci[2,,] > 1)

apply(ests.cover,1,mean)

boxplot(t(apply(est.mat,1:2,mean)))

boxplot(t(est.mat[4,1:50,]))


