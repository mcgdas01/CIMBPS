require(mvnfast)

set.seed(29347)
nreps<-1000
nsamp<-2000
mu<-c(-1,2,0.5)
Sigma<-0.8^abs(outer(1:3,1:3,'-'))

al<-c(1,-1,1,2,-1,2)
be<-c(2,0,1,1,1,5)

sigY<-1
sigZ<-1

expit<-function(x){return(1/(1+exp(-x)))}


nvec<-c(200,500,1000,2000)

ci.matBB<-array(0,c(length(nvec),12,nreps,2))
ests.matBB<-array(0,c(length(nvec),12,nreps))
coverage.matBB<-array(0,c(length(nvec),12,nreps))

nburn<-50
nthin<-5
nits<-nburn+nthin*nsamp

B<-2000

Xall<-rmvn(max(nvec),mu,Sigma)
Xzall<-cbind(1,Xall,Xall[,1]*Xall[,2],Xall[,2]*Xall[,3])
Xtall<-cbind(1,Xall,Xall[,1]*Xall[,2],Xall[,1]*Xall[,3],Xall[,2]*Xall[,3],Xall[,1]*Xall[,2]*Xall[,3])
etaZall<-Xzall%*% al
Zall<-rnorm(max(nvec),etaZall,sigZ)

Xyall<-cbind(1,Zall,Xall,Xall[,2]*Xall[,3])

muYall<-Xyall %*% be

Yall<-rnorm(max(nvec))*sigY+muYall

wmatall<-matrix(rgamma(max(nvec)*B,1),nrow=B,ncol=max(nvec))

ico<-0
for(nval in 1:length(nvec)){
	n<-nvec[nval]
	t0<-proc.time()[3]
	for(irep in 1:nreps){
		ivec<-sample(1:max(nvec),size=n,rep=T)
		X<-Xall[ivec,]
		Xz<-Xzall[ivec,]
		Xt<-Xtall[ivec,]
		etaZ<-etaZall[ivec]

		Z<-Zall[ivec]
		zsq<-sum(Z^2)

		Xy<-Xyall[ivec,]

		muY<-muYall[ivec]

		Y<-Yall[ivec]	

		wmat<-wmatall[,ivec]
		wmat<-n*wmat/apply(wmat,1,sum)
		wY<-wmat*matrix(Y,nrow=B,ncol=n,byrow=T)
		wZ<-wmat*matrix(Z,nrow=B,ncol=n,byrow=T)
		al.ests<-matrix(0,nrow=B,ncol=ncol(Xt))
		b.ests<-matrix(0,nrow=B,ncol=n)




		#Method 1: PS True

		X12y<-cbind(1,Z,etaZ)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X12yw<-X12y*sqrt(w)
			Minv<-crossprod(X12yw)
			M<-solve(Minv)
			m<-M %*% (t(X12y) %*% wY[b,])
			old.be<-m
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,1,irep,1]<-lui[1]
		ci.matBB[nval,1,irep,2]<-lui[2]
		coverage.matBB[nval,1,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,1,irep]<-mean(est)

		#Method 2: PS+RA, True
		X13y<-cbind(1,Z,X,etaZ)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X13yw<-X13y*sqrt(w)
			Minv<-crossprod(X13yw)
			M<-solve(Minv)
			m<-M %*% (t(X13y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,2,irep,1]<-lui[1]
		ci.matBB[nval,2,irep,2]<-lui[2]
		coverage.matBB[nval,2,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,2,irep]<-mean(est)


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
			w<-wmat[b,]
			X14yw<-X14y*sqrt(w)
			Minv<-crossprod(X14yw)
			M<-solve(Minv)
			m<-M %*% (t(X14y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,3,irep,1]<-lui[1]
		ci.matBB[nval,3,irep,2]<-lui[2]
		coverage.matBB[nval,3,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,3,irep]<-mean(est)

		#Method 4: CF+Ext, Parametric
		est<-rep(0,B)
		X15y<-cbind(1,Z,X,0)
		for(b in 1:B){
			old.al<-rmvn(1,mZ,old.sigZ[b]^2*MZ)[1,]
			old.b<-Xt %*% old.al
			X15y[,ncol(X15y)]<-old.b
			w<-wmat[b,]
			X15yw<-X15y*sqrt(w)
			Minv<-crossprod(X15yw)
			M<-solve(Minv)
			m<-M %*% (t(X15y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,4,irep,1]<-lui[1]
		ci.matBB[nval,4,irep,2]<-lui[2]
		coverage.matBB[nval,4,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,4,irep]<-mean(est)

		#Method 5: 2S, Parametric
		MZinv<-crossprod(Xt)
		MZ<-solve(MZinv)
		mZ<-MZ %*% (t(Xt) %*% Z)
		old.b<-Xt %*% mZ			
		X16y<-cbind(1,Z,old.b)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X16yw<-X16y*sqrt(w)
			Minv<-crossprod(X16yw)
			M<-solve(Minv)
			m<-M %*% (t(X16y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,5,irep,1]<-lui[1]
		ci.matBB[nval,5,irep,2]<-lui[2]
		coverage.matBB[nval,5,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,5,irep]<-mean(est)

		#Method 6: 2S+Ext, Parametric
		est<-rep(0,B)
		X17y<-cbind(1,Z,X,old.b)
		for(b in 1:B){
			w<-wmat[b,]
			X17yw<-X17y*sqrt(w)
			Minv<-crossprod(X17yw)
			M<-solve(Minv)
			m<-M %*% (t(X17y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,6,irep,1]<-lui[1]
		ci.matBB[nval,6,irep,2]<-lui[2]
		coverage.matBB[nval,6,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,6,irep]<-mean(est)

		#Method 7: CF, Unlinked
		est<-rep(0,B)
		X5y<-cbind(1,Z,0)

		for(b in 1:B){
			w<-wmat[b,]
			Xtw<-Xt*sqrt(w)
			MZinv<-crossprod(Xtw)
			MZ<-solve(MZinv)
			mZ<-MZ %*% (t(Xt) %*% wZ[b,])
			old.al<-mZ
			al.ests[b,]<-old.al
			old.b<-Xt %*% old.al
			b.ests[b,]<-old.b
		}

		#Break the link by permutation
		irow<-sample(1:B)
		for(b in 1:B){
			X5y[,3]<-b.ests[irow[b],]
			X5yw<-X5y*sqrt(w)
			Minv<-crossprod(X5yw)
			M<-solve(Minv)
			m<-M %*% (t(X5y) %*% (wY[b,]))
			est[b]<-m[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,7,irep,1]<-lui[1]
		ci.matBB[nval,7,irep,2]<-lui[2]
		coverage.matBB[nval,7,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,7,irep]<-mean(est)

		#Method 8: CF+RA, Unlinked
		est<-rep(0,B)
		X6y<-cbind(1,Z,X,0)
		for(b in 1:B){
			w<-wmat[b,]
			X6y[,ncol(X6y)]<-b.ests[irow[b],]
			X6yw<-X6y*sqrt(w)
			Minv<-crossprod(X6yw)
			M<-solve(Minv)
			m<-M %*% (t(X6y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,8,irep,1]<-lui[1]
		ci.matBB[nval,8,irep,2]<-lui[2]
		coverage.matBB[nval,8,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,8,irep]<-mean(est)

     		#Method 9: 2S, Unlinked
		old.b<-Xt %*% apply(al.ests,2,mean)

		X7y<-cbind(1,Z,old.b)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X7yw<-X7y*sqrt(w)
			Minv<-crossprod(X7yw)
			M<-solve(Minv)
			m<-M %*% (t(X7y) %*% (wY[b,]))
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,9,irep,1]<-lui[1]
		ci.matBB[nval,9,irep,2]<-lui[2]
		coverage.matBB[nval,9,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,9,irep]<-mean(est)

		#Method 10: 2S+RA, Unlinked
		X8y<-cbind(1,Z,X,old.b)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X8yw<-X8y*sqrt(w)
			Minv<-crossprod(X8yw)
			M<-solve(Minv)
			m<-M %*% (t(X8y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,10,irep,1]<-lui[1]
		ci.matBB[nval,10,irep,2]<-lui[2]
		coverage.matBB[nval,10,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,10,irep]<-mean(est)

		#Method 11: 2S Linked
		X10y<-cbind(1,Z,0)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X10y[,3]<-b.ests[b,]
			X10yw<-X10y*sqrt(w)
			Minv<-crossprod(X10yw)
			M<-solve(Minv)
			m<-M %*% (t(X10y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}	
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,11,irep,1]<-lui[1]
		ci.matBB[nval,11,irep,2]<-lui[2]
		coverage.matBB[nval,11,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,11,irep]<-mean(est)

		#Method 12: 2S+RA Linked
		X11y<-cbind(1,Z,X,0)
		est<-rep(0,B)
		for(b in 1:B){
			w<-wmat[b,]
			X11y[,ncol(X11y)]<-b.ests[b,]
			X11yw<-X11y*sqrt(w)
			Minv<-crossprod(X11yw)
			M<-solve(Minv)
			m<-M %*% (t(X11y) %*% (wY[b,]))
			old.be<-m
			est[b]<-old.be[2]
		}
		lui<-quantile(est,prob=c(0.025,0.975))
		ci.matBB[nval,12,irep,1]<-lui[1]
		ci.matBB[nval,12,irep,2]<-lui[2]
		coverage.matBB[nval,12,irep]<-as.numeric(lui[1] < be[2] & lui[2] > be[2])
		ests.matBB[nval,12,irep]<-mean(est)

		t1<-proc.time()[3]
		if(irep %% 100 == 0) {print(c(nval,n,irep,as.numeric(t1-t0)))}
	}
	cvec<-apply(coverage.matBB[nval,,],1,mean)
	print(c(nval,n,cvec))
}

bias<-round(apply(ests.matBB,c(1,2),mean)-be[2],3)
t(bias)

variance<-round(apply(ests.matBB,c(1,2),var),3)

rmse<-apply((ests.matBB-be[2])^2,c(1,2),mean)
t(round(sqrt(rmse),3))

cover<-round(100*apply(coverage.matBB,c(1,2),mean),3)
t(cover)
round(sqrt(apply(ests.matBB,c(1,2),var)*nvec),3)
