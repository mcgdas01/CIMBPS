require(mvnfast)
require(mnormt)

nreps<-2000
B<-2000
nsamp<-B

mu<-c(1,1,-1,-1)
Sigma<-diag(rep(1,4))

al.mat<-matrix(c(0,0.3,0.8,0.3,0.8,0.5,0.5,0.75,1,1,0,0.45,0.90,1.35,1.80),nrow=3,byrow=T)

be<-c(0.25,0.25,0.25,0.25,1.5)

true.tau<-0

sigY<-1
sigZ<-1

expit<-function(x){return(1/(1+exp(-x)))}

jt.like<-function(xp,zv,yv,Xa,Xy,sY){

	alv<-xp[1:ncol(Xa)]
	bev<-xp[-c(1:ncol(Xa))]
	eZ<-Xa %*% alv
	bv<-expit(eZ)
	zlike<-sum(zv*eZ-log(1+exp(eZ)))
	eY<-cbind(Xy,bv) %*% bev
	ylike<-sum(dnorm(yv,eY,sY,log=T))
	return(list(like=zlike+ylike,zfit=bv,yfit=eY))
}

jt.like.val<-function(xp,zv,yv,Xa,Xy,sY){
	alv<-xp[1:ncol(Xa)]
	bev<-xp[-c(1:ncol(Xa))]
	eZ<-Xa %*% alv
	bv<-expit(eZ)
	zlike<-sum(zv*eZ-log(1+exp(eZ)))
	eY<-cbind(Xy,bv) %*% bev
	ylike<-sum(dnorm(yv,eY,sY,log=T))
	return(zlike+ylike)
}

log.likeZ<-function(xp,zv,Xa){
	alv<-xp
	eZ<-Xa %*% alv
	bv<-expit(eZ)
	zlike<-sum(zv*eZ-log(1+exp(eZ)))
	return(list(like=zlike,zfit=bv))
}


nvec<-c(200,500,1000,2000)

ci.matBB<-array(0,c(3,length(nvec),12,nreps,2))
ests.matBB<-array(0,c(3,length(nvec),12,nreps))
coverage.matBB<-array(0,c(3,length(nvec),12,nreps))

nburn<-50
nthin<-5
nits<-nburn+nthin*nsamp

set.seed(1327)

ico<-0
for(k in 1:3){
	t0<-proc.time()[3]
	al<-al.mat[k,]
	for(nval in 1:length(nvec)){
		n<-nvec[nval]
		for(irep in 1:nreps){
			X1<-rnorm(n,mu[1],1);X2<-rnorm(n,mu[2],1);X3<-rnorm(n,mu[3],1);X4<-rnorm(n,mu[4],1)
			Xz<-cbind(1,X1,X2,X3,X4)
			Xt<-cbind(Xz,X1*X2,X1*X3,X2*X3,X1*X2*X3)
			etaZ<-Xz%*% al
			pZ<-expit(etaZ)
			Z<-rbinom(n,1,pZ)
		
			Xy<-cbind(Xz[,2:5],X3*X4)
			muY<-Xy %*% be
			Y<-rnorm(n)*sigY+muY	
			ysq<-sum(Y^2)	

			wmat<-matrix(rgamma(n*B,1),nrow=B,ncol=n)
			wmat<-n*wmat/apply(wmat,1,sum)
			swmat<-sqrt(wmat)
			wY<-wmat*matrix(Y,nrow=B,ncol=n,byrow=T)
			al.ests<-matrix(0,nrow=B,ncol=ncol(Xt))
			b.ests<-matrix(0,nrow=B,ncol=n)


			#Method 1: True PS
			imeth<-1
			X1y<-cbind(1,Z,pZ)
			est<-rep(0,B)
			for(b in 1:B){
				X1yw<-X1y*swmat[b,]
				Minv<-crossprod(X1yw)
				M<-solve(Minv)
				m<-M %*% (t(X1y) %*% (wY[b,]))
				old.be<-m
				est[b]<-old.be[2]
			}	
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,1,irep,1]<-lui[1]
			ci.matBB[k,nval,1,irep,2]<-lui[2]
			coverage.matBB[k,nval,1,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,1,irep]<-mean(est)

			next

			#Method 2: True PS Ext
			imeth<-2
			X2y<-cbind(1,Z,X1,X2,X3,X4,pZ)	
			est<-rep(0,B)
			for(b in 1:B){
				X2yw<-X2y*swmat[b,]
				Minv<-crossprod(X2yw)
				M<-solve(Minv)
				m<-M %*% (t(X2y) %*% (wY[b,]))
				old.be<-m
				est[b]<-old.be[2]
			}	
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,2,irep,1]<-lui[1]
			ci.matBB[k,nval,2,irep,2]<-lui[2]
			coverage.matBB[k,nval,2,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,2,irep]<-mean(est)


			#Method 3: CF, Parametric
			imeth<-3
			fZ<-glm(Z~-1+Xt,family=binomial)
			old.al<-coef(fZ)
			old.b<-fitted(fZ)
			X3y<-cbind(1,Z,old.b)
			M<-solve(crossprod(X3y))
			m<-M %*% (t(X3y) %*% Y)
			old.be<-m

			X0<-cbind(1,Z)
			old.sigY<-sqrt(mean((Y - cbind(X0,old.b)%*%old.be)^2))
			imean<-old.al	
			ivar<-vcov(fZ)
			old.val<-logLik(fZ)[1]	
			old.par<-old.al
			old.like<-log.likeZ(old.par,Z,Xt)

			luv<-log(runif(nits))
			par.vals<-rmvn(nits,imean,ivar)
	
			est<-rep(0,nsamp)
			ico<-0
			al.ests<-matrix(0,nrow=B,ncol=ncol(Xt))
			b.ests<-matrix(0,nrow=B,ncol=n)

			for(iter in 1:nits){
				new.par<-par.vals[iter,]
				new.like<-log.likeZ(new.par,Z,Xt)
				new.val<-new.like[[1]]
				new.ival<-dmvn(new.par,imean,ivar,log=T)
				if(luv[iter] < (new.val-new.ival) - (old.val-old.ival)){
					old.par<-new.par
					old.val<-new.val
					old.like<-new.like
					new.ival<-old.ival
				}

				old.b<-old.like$zfit

				if(iter > nburn & iter %% nthin == 0){
					ico<-ico+1;
					al.ests[ico,]<-old.par
					b.ests[ico,]<-old.b
					X3y[,3]<-old.b
					X3yw<-X3y*swmat[ico,]
					M<-solve(crossprod(X3yw))
					m<-M %*% (t(X3y) %*% wY[ico,])
					old.be<-m
					est[ico]<-old.be[2]
				}
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,3,irep,1]<-lui[1]
			ci.matBB[k,nval,3,irep,2]<-lui[2]
			coverage.matBB[k,nval,3,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,3,irep]<-mean(est)


			#Method 4: CF+RA, Parametric
			imeth<-4
			X4y<-cbind(1,Z,X1,X2,X3,X4,old.b)
			est<-rep(0,nsamp)
			for(b in 1:B){
				old.b<-b.ests[b,]
				X4y[,ncol(X4y)]<-old.b
				X4yw<-X4y*swmat[b,]
				M<-solve(crossprod(X4yw))
				m<-M %*% (t(X4y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,4,irep,1]<-lui[1]
			ci.matBB[k,nval,4,irep,2]<-lui[2]
			coverage.matBB[k,nval,4,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,4,irep]<-mean(est)


			#Method 5: 2S,  Parametric
			imeth<-5

			old.b<-expit(Xt %*% apply(al.ests,2,mean))

			X5y<-cbind(1,Z,old.b)
			est<-rep(0,nsamp)
			for(b in 1:B){
				X5yw<-X5y*swmat[b,]
				M<-solve(crossprod(X5yw))
				m<-M %*% (t(X5y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,5,irep,1]<-lui[1]
			ci.matBB[k,nval,5,irep,2]<-lui[2]
			coverage.matBB[k,nval,5,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,5,irep]<-mean(est)


			#Method 6: 2S+RA,  Parametric
			imeth<-6

			X6y<-cbind(1,Z,X1,X2,X3,X4,old.b)
			est<-rep(0,nsamp)
			for(b in 1:B){
				X6yw<-X6y*swmat[b,]
				M<-solve(crossprod(X6yw))
				m<-M %*% (t(X6y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,6,irep,1]<-lui[1]
			ci.matBB[k,nval,6,irep,2]<-lui[2]
			coverage.matBB[k,nval,6,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,6,irep]<-mean(est)


			#Method 7: CF, Unlinked
			imeth<-7
			est<-rep(0,B)
			ivec<-sample(1:B)
			a.samp<-matrix(0,nrow=B,ncol=ncol(Xt))
			b.samp<-matrix(0,nrow=B,ncol=n)
			X7y<-cbind(1,Z,0)
			for(b in 1:B){
				fZ<-glm(Z~-1+Xt,family=binomial,weights=wmat[ivec[b],])
				old.al<-coef(fZ)
				a.samp[b,]<-old.al
				old.b<-fitted(fZ)
				b.samp[ivec[b],]<-old.b
				X7y[,3]<-old.b
				X7yw<-X7y*swmat[b,]
				M<-solve(crossprod(X7yw))
				m<-M %*% (t(X7y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,7,irep,1]<-lui[1]
			ci.matBB[k,nval,7,irep,2]<-lui[2]
			coverage.matBB[k,nval,7,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,7,irep]<-mean(est)
			
			#Method 8: CF+RA, Unlinked
			imeth<-8
			est<-rep(0,B)
			X8y<-cbind(1,Z,X1,X2,X3,X4,0)
			for(b in 1:B){
				X8y[,ncol(X8y)]<-b.samp[ivec[b],]
				X8yw<-X8y*swmat[b,]
				M<-solve(crossprod(X8yw))
				m<-M %*% (t(X8y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,8,irep,1]<-lui[1]
			ci.matBB[k,nval,8,irep,2]<-lui[2]
			coverage.matBB[k,nval,8,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,8,irep]<-mean(est)

			#Method 9: 2S,  Unlinked
			imeth<-9

			old.b<-expit(Xt %*% apply(a.samp,2,mean))

			X9y<-cbind(1,Z,old.b)
			est<-rep(0,nsamp)
			for(b in 1:B){
				X9yw<-X9y*swmat[b,]
				M<-solve(crossprod(X9yw))
				m<-M %*% (t(X9y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,9,irep,1]<-lui[1]
			ci.matBB[k,nval,9,irep,2]<-lui[2]
			coverage.matBB[k,nval,9,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,9,irep]<-mean(est)


			#Method 10: 2S+RA,  Unlinked
			imeth<-10

			X10y<-cbind(1,Z,X1,X2,X3,X4,old.b)
			est<-rep(0,nsamp)
			for(b in 1:B){
				X10yw<-X10y*swmat[b,]
				M<-solve(crossprod(X10yw))
				m<-M %*% (t(X10y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,10,irep,1]<-lui[1]
			ci.matBB[k,nval,10,irep,2]<-lui[2]
			coverage.matBB[k,nval,10,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,10,irep]<-mean(est)

			#Method 11: 2S, Linked
			imeth<-11
			est<-rep(0,B)
			X11y<-cbind(1,Z,0)
			for(b in 1:B){
				X11y[,3]<-b.samp[b,]
				X11yw<-X11y*swmat[b,]
				M<-solve(crossprod(X11yw))
				m<-M %*% (t(X11y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,11,irep,1]<-lui[1]
			ci.matBB[k,nval,11,irep,2]<-lui[2]
			coverage.matBB[k,nval,11,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,11,irep]<-mean(est)
			
			#Method 12: 2S+RA, Linked
			imeth<-12
			est<-rep(0,B)
			X12y<-cbind(1,Z,X1,X2,X3,X4,0)
			for(b in 1:B){
				X12y[,ncol(X12y)]<-b.samp[b,]
				X12yw<-X12y*swmat[b,]
				M<-solve(crossprod(X12yw))
				m<-M %*% (t(X12y) %*% wY[b,])
				old.be<-m
				est[b]<-old.be[2]
			}
			lui<-quantile(est,prob=c(0.025,0.975))
			ci.matBB[k,nval,12,irep,1]<-lui[1]
			ci.matBB[k,nval,12,irep,2]<-lui[2]
			coverage.matBB[k,nval,12,irep]<-as.numeric(lui[1] < true.tau & lui[2] > true.tau)
			ests.matBB[k,nval,12,irep]<-mean(est)

			t1<-proc.time()[3]
			if(irep %% 100 == 0) {print(c(k,nval,n,irep,as.numeric(t1-t0)))}
		}
		cvec<-apply(coverage.matBB[k,nval,,],1,mean)
		print(c(k,nval,n,cvec))
	}
}


round(100*apply(coverage.matBB,c(2:3,1),mean),3)
round(apply(ests.matBB,c(2:3,1),mean),3)
round(apply(ests.matBB,c(2:3,1),var),3)

rmse<-apply((ests.matBB)^2,c(2:3,1),mean)
t(round(sqrt(rmse),3))

round(sqrt(apply(ests.matBB,c(2:3,1),var)*nvec),3)

####################################################################################

Ex3.BB.rmse<-t(round(sqrt(rbind(rmse[,,1],rmse[,,2],rmse[,,3])),3))
rnames<-c('PS           True       ','PS-ext       True','CF           Parametric ',
'CF-ext       Parametric ','2S           Parametric ',
'2S-ext       Parametric ','CF           Unlinked BB','CF-ext       Unlinked BB',
'2S           Unlinked BB','2S-ext       Unlinked BB','2S           Linked B',
'2S-ext       Linked BB  ')
cnames<-rep(c('200','500','1000','2000'),3)
row.names(Ex3.BB.rmse)<-rnames
colnames(Ex3.BB.rmse)<-cnames

print(xtable(as.data.frame(Ex3.BB.rmse), type = "latex",digits=3), file = "Tab-Ex3-BB.tex")

covermat<-100*apply(coverage.matBB,c(2:3,1),mean)

Ex3.BB.cov<-t(rbind(covermat[,,1],covermat[,,2],covermat[,,3]))
row.names(Ex3.BB.cov)<-rnames
colnames(Ex3.BB.cov)<-cnames
print(xtable(as.data.frame(Ex3.BB.cov), type = "latex",digits=1), file = "Tab-Ex3-BB.tex",append=TRUE)



