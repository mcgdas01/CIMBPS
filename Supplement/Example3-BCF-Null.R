#install.packages(c('mvnfast','MCMCpack','dbarts','bcf'),dependencies=TRUE)

require(mvnfast)
library('MCMCpack')
library(dbarts)
library(bcf)

expit <- function(x){return(1/(1+exp(-x)))}
Cases <- matrix(NA,nrow = 3,ncol = 5)
Cases[1,] <- c(0,.3,.8,.3,.8)
Cases[2,] <- c(.5,.5,.75,1,1)
Cases[3,] <- seq(0,1.8,len=5)

nvec <- c(200,500,1000,2000)
nrep  <-2000;
nrepBB<-2000;

ests.mat<-array(dim = c(nrep,2,3,4))
ests.mat.bcf<-array(dim = c(nrep,2,3,4,nrepBB))
covered<-array(dim = c(nrep,2,3,4))

ests.mat.var<-array(dim = c(nrep,2,3,4))
ests.mat.025<-array(dim = c(nrep,2,3,4))
ests.mat.975<-array(dim = c(nrep,2,3,4))

ests.matBoot<-array(dim = c(nrep,2,3,4,nrepBB))
ests.matBoot.var<-array(dim = c(nrep,2,3,4,nrepBB))
ests.matBoot.025<-array(dim = c(nrep,2,3,4,nrepBB))
ests.matBoot.975<-array(dim = c(nrep,2,3,4,nrepBB))

betaZ<-0
betaX.y=c(rep(1,4),rep(.75,4))

#Dirichlet weights
dw<-matrix(rgamma(2000*nrepBB,1),nrow=nrepBB)

v0<-10
a0<-1
b0<-1

set.seed(2364)
nsamp<-2000
tot.time<-array(0,c(2,3,4))

for(w in 1:4){
  n <- nvec[w]
  ones<-rep(1,n)
  mones<-matrix(ones)
  for(k in 1:3){
    gam <- Cases[k,]
    for(j in 1:nrep){
      set.seed(100*j + 8768*w + 134599*k)
      X <- cbind(rnorm(n,1,1),rnorm(n,1,1),rnorm(n,-1,1),rnorm(n,-1,1))
      mu.yx <- cbind(X,X[,1]*X[,3],X[,2]*X[,4],X[,1]*X[,4],X[,3]*X[,4])%*%betaX.y
      logit.ps <- cbind(1,X)%*%gam
      mu.z <- 1/(1+exp(-logit.ps))
      Z <- rbinom(n,1,mu.z)
      mu.y <- mu.yx + Z*betaZ
      Y <- rnorm(n,mu.y,1)
	ysq<-sum(Y^2)
      
	pz<-ncol(X)+1
	wv<-n*rep(1/n,n)
	ps.plug<-glm(Z~X,family=quasibinomial,weights=wv)
	#be.mle<-coef(ps.plug)
	psfit<-fitted(ps.plug)

      # two-step
	#Xm<-cbind(1,Z,X,psfit)
	#V0<-diag(v0,ncol(Xm))
	#V0inv<-diag(1/v0,ncol(Xm))

	#Vninv<-t(Xm) %*% Xm + V0inv
	#Vn<-solve(Vninv)
	#vn<-Vn %*% (t(Xm) %*% Y)
	#ests.mat[j,1,k,w]<-vn[2]

      # BCF
	 
	sink("sink-bcf.txt")
	t0<-as.numeric(proc.time()[3])
      fit_default = bcf(Y, Z, x_control = X,x_moderate = mones, nburn = 500, nsim = nrepBB, pihat = matrix(psfit), update_interval = 100000, include_pi = "both")
	tot.time[2,k,w]<-tot.time[2,k,w]+(as.numeric(proc.time()[3])-t0)
	sink()
      # posterior mean
	effects<-rowMeans(fit_default$tau)
      ests.mat.bcf[j,2,k,w,] <- effects
      ests.mat[j,2,k,w] <- mean(effects)
      ests.mat.var[j,2,k,w] <- var(effects)
      # quantiles 	
	qv<-quantile(effects,prob=c(.025,.975))
      ests.mat.025[j,2,k,w] <- qv[1]
      ests.mat.975[j,2,k,w] <- qv[2]

	covered[j,2,k,w]<-as.numeric(betaZ > qv[1] & betaZ < qv[2])
      
      # loop Bayesian bootstrap

	t0<-as.numeric(proc.time()[3])
	X0<-cbind(1,Z,0)
	tX0<-rbind(1,Z,0)

      for(irep in 1:nrepBB){

		wv<-dw[irep,1:n]
		wv<-n*wv/sum(wv)

		ps.plug<-glm(Z~X,family=quasibinomial,weights=wv)
		#be.mle<-coef(ps.plug)
		psfit<-fitted(ps.plug)

		X0[,ncol(X0)]<-psfit
		tX0[nrow(tX0),]<-psfit
		Xm<-X0*sqrt(wv)
		Vninv<-crossprod(Xm)
		#Vn<-solve(Vninv)
		#vn<-Vn %*% (t(X0) %*% (Y*wv))
		vn<-solve(Vninv,tX0 %*% (Y*wv))

		ests.matBoot[j,1,k,w,irep]<-vn[2]

		#if(irep %in% seq(100,nrepBB,len=10)){print(irep);print(timestamp())}
      }
	tot.time[1,k,w]<-tot.time[1,k,w]+(as.numeric(proc.time()[3])-t0)

	ests.mat[j,1,k,w]<-mean(ests.matBoot[j,1,k,w,])


	qv<-quantile(ests.matBoot[j,1,k,w,],prob=c(0.025,0.975))
	covered[j,1,k,w]<-as.numeric(betaZ > qv[1] & betaZ < qv[2])

	cov1<-mean(covered[1:j,1,k,w])
	cov2<-mean(covered[1:j,2,k,w])

	print(sprintf('Scenario %2i n=%4i: replication %4i out of %4i: coverage %5.4f %5.4f Time: %10.2f %10.2f',
			k,n,j,nrep,cov1,cov2,tot.time[1,k,w],tot.time[2,k,w]))

      if(j %in% seq(20,nrep,20)){
         save.image("BCF-analysisBootShort.RData")}
    }
  }
}

bias<-apply(ests.mat,c(2,4,3),mean)
v<-apply(ests.mat,c(2,4,3),var)
rmse<-sqrt(apply((ests.mat)^2,c(2,4,3),mean))
covermat<-round(100*apply(covered,c(2,4,3),mean),1)

####################################################################################

cnames<-c('200','500','1000','2000')

Ex3.tab<-rbind(bias[,,1],rmse[,,1],covermat[,,1])
Ex3.tab<-rbind(Ex3.tab,rbind(bias[,,2],rmse[,,2],covermat[,,2]))
Ex3.tab<-rbind(Ex3.tab,rbind(bias[,,3],rmse[,,3],covermat[,,3]))
colnames(Ex3.tab)<-cnames

Ex3.tab<-data.frame(Method=rep(c('2S','BCF'),9),Ex3.tab)

require(xtable)
print(xtable(Ex3.tab, type = "latex",digits=3), include.rownames=FALSE, file = "Tab-Ex3-BCF-Null-extra.tex")


