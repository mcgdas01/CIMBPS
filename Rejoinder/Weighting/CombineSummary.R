#Ests0 - unweighted
#Ests1 - CF weighted
#Ests2 - 2S weighted
#Ests3 - Linked BB
#Ests4 - Unlinked BB
#Estsc - Correct model

apply(apply(ests0,1:2,var),1,mean)
apply(apply(ests1,1:2,var),1,mean)
apply(apply(ests2,1:2,var),1,mean)
apply(apply(ests3,1:2,var),1,mean)
apply(apply(ests4,1:2,var),1,mean)
apply(apply(estsc,1:2,var),1,mean)

#Coverage
all.ests<-array(0,c(6,dim(ests0)))
all.ests[1,,,]<-ests0
all.ests[2,,,]<-ests1
all.ests[3,,,]<-ests2
all.ests[4,,,]<-ests3
all.ests[5,,,]<-ests4
all.ests[6,,,]<-estsc

all.ci<-apply(all.ests,1:3,quantile,c(0.025,0.975))
all.cover<-apply(all.ci,2:4,function(x){x[1]<1 & x[2]>1})

apply(all.cover,1:2,mean)

ate.ests<-apply(all.ests,1:3,mean)

ate.bias<-apply(ate.ests-1,1:2,mean)

round(ate.bias,3)

ate.rmse<-sqrt(apply((ate.ests-1)^2,1:2,mean))

round(ate.rmse,3)