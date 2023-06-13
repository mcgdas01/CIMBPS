require(mvnfast)
require(mnormt)

set.seed(29343547)
nreps<-1000
nsamp<-2000
mu<-c(-1,2,0.5)
Sigma<-0.8^abs(outer(1:3,1:3,'-'))

al<-c(1,-1,1,2,-1,2)
be<-c(2,5,1,1,1,5)

sigY<-1
sigZ<-1

expit<-function(x){return(1/(1+exp(-x)))}


nvec<-c(200,500,1000,2000)

ci.matBBnew<-array(0,c(length(nvec),16,nreps,2))
ests.matBBnew<-array(0,c(length(nvec),16,nreps))
coverage.matBBnew<-array(0,c(length(nvec),16,nreps))

nburn<-50
nthin<-5
nits<-nburn+nthin*nsamp

B<-1000

ico<-0
for(nval in 1:length(nvec)){
	n<-nvec[nval]
	t0<-proc.time()[3]
	for(irep in 1:nreps){
		X<-rmvn(n,mu,Sigma)
		Xz<-cbind(1,X,X[,1]*X[,2],X[,2]*X[,3])
		Xt<-cbind(1,X,X[,1]*X[,2],X[,1]*X[,3],X[,2]*X[,3],X[,1]*X[,2]*X[,3])
		etaZ<-Xz%*% al
		Z<-rnorm(n,etaZ,sigZ)
		zsq<-sum(Z^2)

		Xy<-cbind(1,Z,X,X[,2]*X[,3])
		muY<-Xy %*% be
		Y<-rnorm(n)*sigY+muY	

		wmat<-matrix(rgamma(n*B,1),nrow=B,ncol=n)
		wmat<-n*wmat/apply(wmat,1,sum)
		wY<-wmat*matrix(Y,nrow=B,ncol=n,byrow=T)
		wZ<-wmat*matrix(Z,nrow=B,ncol=n,byrow=T)
		al.ests<-matrix(0,nrow=B,ncol=ncol(Xt))
		b.ests<-matrix(0,nrow=B,ncol=n)

		#Method 1: PS True

		swmat<-sqrt(wmat)

		X12y<-cbind(1,Z,etaZ)
		est<-rep(0,B)
		for(b in 1:B){
			X12yw<-X12y*swmat[b,]
			Minv<-crossprod(X12yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X12y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,1,irep,1]<-lui[1]
		ci.matBBnew[nval,1,irep,2]<-lui[2]
		coverage.matBBnew[nval,1,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,1,irep]<-mean(est)

		#Method 2: PS+RA, True
		X13y<-cbind(1,Z,X,etaZ)
		est<-rep(0,B)
		for(b in 1:B){
			X13yw<-X13y*swmat[b,]
			Minv<-crossprod(X13yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X13y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,2,irep,1]<-lui[1]
		ci.matBBnew[nval,2,irep,2]<-lui[2]
		coverage.matBBnew[nval,2,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,2,irep]<-mean(est)


		#Method 3: CF, Parametric

		MZinv<-crossprod(Xt)
		MZ<-solve(MZinv)
		mZ<-MZ %*% (t(Xt) %*% Z)
		cz<-zsq-(t(mZ) %*% (MZinv %*% mZ))[1,1]
		old.sigZ<-1/sqrt(rgamma(B,n/2,cz/2))
		est<-rep(0,B)
		X14y<-cbind(1,Z,0)
		for(b in 1:B){
			old.al<-rmvn(1,mZ,old.sigZ[b]^2*MZ)[1,]
			old.b<-Xt %*% old.al			
			X14y[,3]<-old.b
			X14yw<-X14y*swmat[b,]
			Minv<-crossprod(X14yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X14y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,3,irep,1]<-lui[1]
		ci.matBBnew[nval,3,irep,2]<-lui[2]
		coverage.matBBnew[nval,3,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,3,irep]<-mean(est)

		#Method 4: CF+Ext, Parametric
		est<-rep(0,B)
		X15y<-cbind(1,Z,X,0)
		for(b in 1:B){
			old.al<-rmvn(1,mZ,old.sigZ[b]^2*MZ)[1,]
			old.b<-Xt %*% old.al
			X15y[,ncol(X15y)]<-old.b
			X15yw<-X15y*swmat[b,]
			Minv<-crossprod(X15yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X15y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,4,irep,1]<-lui[1]
		ci.matBBnew[nval,4,irep,2]<-lui[2]
		coverage.matBBnew[nval,4,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,4,irep]<-mean(est)

		#Method 5: 2S, Parametric
		MZinv<-crossprod(Xt)
		MZ<-solve(MZinv)
		mZ<-MZ %*% (t(Xt) %*% Z)
		old.b<-Xt %*% mZ			
		X16y<-cbind(1,Z,old.b)
		est<-rep(0,B)
		for(b in 1:B){
			X16yw<-X16y*swmat[b,]
			Minv<-crossprod(X16yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X16y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,5,irep,1]<-lui[1]
		ci.matBBnew[nval,5,irep,2]<-lui[2]
		coverage.matBBnew[nval,5,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,5,irep]<-mean(est)

		#Method 6: 2S+Ext, Parametric
		est<-rep(0,B)
		X17y<-cbind(1,Z,X,old.b)
		for(b in 1:B){
			X17yw<-X17y*swmat[b,]
			Minv<-crossprod(X17yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X17y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,6,irep,1]<-lui[1]
		ci.matBBnew[nval,6,irep,2]<-lui[2]
		coverage.matBBnew[nval,6,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,6,irep]<-mean(est)

		#Method 7: CF, Unlinked
		est<-rep(0,B)
		X5y<-cbind(1,Z,0)

		lvec<-rep(0,B)

		for(b in 1:B){
			Xtw<-Xt*swmat[b,]
			MZinv<-crossprod(Xtw)
			MZ<-solve(MZinv)
			mZ<-MZ %*% crossprod(Xt,wZ[b,])
			old.al<-mZ
			al.ests[b,]<-old.al
			lvec[b]<-sum((Z- Xt %*% old.al)^2)
			old.b<-Xt %*% old.al
			b.ests[b,]<-old.b
		}

		#Break the link by permutation
		#Take a sample by compiling the estimates and shuffling them

		irow<-sample(1:B)
		for(b in 1:B){
			X5y[,3]<-b.ests[irow[b],]
			X5yw<-X5y*swmat[b,]
			Minv<-crossprod(X5yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X5y,wY[b,])
			est[b]<-m[2]
		}


		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,7,irep,1]<-lui[1]
		ci.matBBnew[nval,7,irep,2]<-lui[2]
		coverage.matBBnew[nval,7,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,7,irep]<-mean(est)

		#Method 8: CF+RA, Unlinked
		est<-rep(0,B)
		X6y<-cbind(1,Z,X,0)
		for(b in 1:B){
			X6y[,ncol(X6y)]<-b.ests[irow[b],]
			X6yw<-X6y*swmat[b,]
			Minv<-crossprod(X6yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X6y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,8,irep,1]<-lui[1]
		ci.matBBnew[nval,8,irep,2]<-lui[2]
		coverage.matBBnew[nval,8,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,8,irep]<-mean(est)

     		#Method 9: 2S, Unlinked
		#old.b<-Xt %*% apply(al.ests,2,mean)
		old.b<-Xt %*% al.ests[which.min(lvec),]

		X7y<-cbind(1,Z,old.b)
		est<-rep(0,B)
		for(b in 1:B){
			X7yw<-X7y*swmat[b,]
			Minv<-crossprod(X7yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X7y,wY[b,])
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,9,irep,1]<-lui[1]
		ci.matBBnew[nval,9,irep,2]<-lui[2]
		coverage.matBBnew[nval,9,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,9,irep]<-mean(est)

		#Method 10: 2S+RA, Unlinked
		X8y<-cbind(1,Z,X,old.b)
		est<-rep(0,B)
		for(b in 1:B){
			X8yw<-X8y*swmat[b,]
			Minv<-crossprod(X8yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X8y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,10,irep,1]<-lui[1]
		ci.matBBnew[nval,10,irep,2]<-lui[2]
		coverage.matBBnew[nval,10,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,10,irep]<-mean(est)


#######################################################################################
		#Method 11: 2S Linked (is it really CF linked)
		#Here 'linked' means linked through the BB, not in a 
		#Change to do this properly in methods 13 and 14?

		#Need to find some way to link the posterior calculations ??

		X10y<-cbind(1,Z,0)
		est<-rep(0,B)
		for(b in 1:B){
			X10y[,3]<-b.ests[b,]
			X10yw<-X10y*swmat[b,]
			Minv<-crossprod(X10yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X10y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,11,irep,1]<-lui[1]
		ci.matBBnew[nval,11,irep,2]<-lui[2]
		coverage.matBBnew[nval,11,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,11,irep]<-mean(est)

		#Method 12: 2S+RA Linked
		X11y<-cbind(1,Z,X,0)
		est<-rep(0,B)
		for(b in 1:B){
			X11y[,ncol(X11y)]<-b.ests[b,]
			X11yw<-X11y*swmat[b,]
			Minv<-crossprod(X11yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X11y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,12,irep,1]<-lui[1]
		ci.matBBnew[nval,12,irep,2]<-lui[2]
		coverage.matBBnew[nval,12,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,12,irep]<-mean(est)

#######################################################################################

		#Method 13: 2S Linked 
		#Here 'linked' means linked through the BB, not in a 
		#Change to do this properly ?

		#Need to find some way to link the posterior calculations ??

		X10y<-cbind(1,Z,0)
		al.ests<-ai.ests<-matrix(0,nrow=B,ncol=ncol(Xt))
		lvec<-rep(0,B)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			wimat<-matrix(rgamma(n*B,1),ncol=n)
			for(bi in 1:B){
				wi<-wimat[bi,]
				wi<-w*(wi/sum(wi))
				Xtw<-Xt*sqrt(wi)
				MZinv<-crossprod(Xtw)
				MZ<-solve(MZinv)
				mZ<-MZ %*% crossprod(Xt,wi*Z)
				ai.ests[bi,]<-mZ
				lvec[bi]<-sum((Z-Xt %*% mZ)^2)
			}	
			#al.ests[b,]<-apply(ai.ests,2,mean)
			al.ests[b,]<-ai.ests[which.min(lvec),]

			old.b<-Xt %*% al.ests[b,]
			b.ests[b,]<-old.b
			X10y[,3]<-old.b
			X10yw<-X10y*swmat[b,]
			Minv<-crossprod(X10yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X10y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,13,irep,1]<-lui[1]
		ci.matBBnew[nval,13,irep,2]<-lui[2]
		coverage.matBBnew[nval,13,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,13,irep]<-mean(est)

		#Method 14: 2S+RA Linked
		X11y<-cbind(1,Z,X,0)
		est<-rep(0,B)
		for(b in 1:B){
			X11y[,ncol(X11y)]<-b.ests[b,]
			X11yw<-X11y*swmat[b,]
			Minv<-crossprod(X11yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X11y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,14,irep,1]<-lui[1]
		ci.matBBnew[nval,14,irep,2]<-lui[2]
		coverage.matBBnew[nval,14,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,14,irep]<-mean(est)

#######################################################################################

		#Method 15: 2S Linked 
		#Here 'linked' means linked through the BB, not in a 
		#Change to do this properly ?

		#Need to find some way to link the posterior calculations ??

		X10y<-cbind(1,Z,0)
		al.ests<-ai.ests<-matrix(0,nrow=B,ncol=ncol(Xt))
		lvec<-rep(0,B)
		est<-rep(0,B)
		wimat<-matrix(rgamma(n*B,1),ncol=n)

		for(b in 1:B){
			w<-wmat[b,]
			wi<-wimat[b,]
			wi<-w*(wi/sum(wi))
			Xtw<-Xt*sqrt(wi)
			MZinv<-crossprod(Xtw)
			MZ<-solve(MZinv)
			mZ<-MZ %*% crossprod(Xt,wi*Z)
			al.ests[b,]<-mZ

			old.b<-Xt %*% al.ests[b,]
			b.ests[b,]<-old.b
			X10y[,3]<-old.b
			X10yw<-X10y*swmat[b,]
			Minv<-crossprod(X10yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X10y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,15,irep,1]<-lui[1]
		ci.matBBnew[nval,15,irep,2]<-lui[2]
		coverage.matBBnew[nval,15,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,15,irep]<-mean(est)

		#Method 16: 2S+RA Linked
		X11y<-cbind(1,Z,X,0)
		est<-rep(0,B)
		for(b in 1:B){
			X11y[,ncol(X11y)]<-b.ests[b,]
			X11yw<-X11y*swmat[b,]
			Minv<-crossprod(X11yw)
			M<-solve(Minv)
			m<-M %*% crossprod(X11y,wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBBnew[nval,16,irep,1]<-lui[1]
		ci.matBBnew[nval,16,irep,2]<-lui[2]
		coverage.matBBnew[nval,16,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBBnew[nval,16,irep]<-mean(est)

		t1<-proc.time()[3]
		if(irep %% 1 == 0) {print(c(nval,n,irep,as.numeric(t1-t0)))}
	}

	cvec<-apply(coverage.matBBnew[nval,,],1,mean)
	print(c(nval,n,cvec))
}

bias<-round(apply(ests.matBBnew,c(1,2),mean)-5,3)
t(bias)

variance<-round(apply(ests.matBBnew,c(1,2),var),3)

rmse<-apply((ests.matBBnew-5)^2,c(1,2),mean)
t(round(sqrt(rmse),3))

cover<-round(100*apply(coverage.matBBnew,c(1,2),mean),3)
t(cover)
round(sqrt(apply(ests.matBBnew,c(1,2),var)*nvec),3)
