n=100
m=100
u=matrix(sample(c(0,1,2),size=n*m,replace=TRUE),nrow=m)
A=0.9
l=1.7
a=0.5
k=3
library(actuar)
lambda.u=matrix(0,nrow=m,ncol=n)
alpha.u=matrix(0,nrow=m,ncol=n)
for(i in 1:m)
{
lambda.u[i,]=l*(A^u[i,])
alpha.u[i,]=a/(1+u[i,])
}
N=matrix(0,nrow=m,ncol=n)
lambda.u.cap=NULL
for(i in 1:m)
{
	for(j in 1:n)
	N[i,j]=rztpois(n=1,lambda=lambda.u[i,j])
	lambda.u.cap[i]=sum(N[i,])/sum(A^u[i,])
}
lambda.u.cap.samp=matrix(0,nrow=m,ncol=n)
for(i in 1:m)
{
	for(j in 1:n)
	{
	lambda.u.cap.samp[i,j]=lambda.u.cap[i]*(A^u[i,j])
	}
}
N.sim=NULL
expec.N=NULL
shi.c=function(a,l,y)
{
	N.sim=rztpois(n=500,lambda=l)
	for(i in 1:length(y))
		expec.N[i]=mean(((a*y[i])^(N.sim-1))/factorial(N.sim-1))
	return(expec.N)
}
M=2
draw.samples=function(n,a,l,M) 
{
  y=NULL
  count=0
  iter=0
  while(count<n) 
	{
    	iter=iter+1
    	g=rexp(1,rate=a)
    	u=runif(1)
    	ratio=shi.c(a,l,g)/M
    		if (u<=ratio) 
		{
      	count=count+1
      	y[count]=g
    		}
    		if (iter> n * 100) 
		 {
      	   print("Target size not reached. Check if M is too small or PDF is valid.")
      	   break
    		 }
  	}
    return(y[1:count])
}

Y=matrix(0,nrow=m,ncol=n)
for(i in 1:m)
{
	for(j in 1:n)
	Y[i,j]=draw.samples(1,alpha.u[i,j],lambda.u.cap.samp[i,j],M)
}

hist(Y[33,],breaks=20,prob=TRUE,main="Samples via Acceptance-Rejection",xlab="y",col="lightblue")
target.pdf=function(y) 
 {
   a*exp(-a*y)*shi.c(a,l,y)
 }
curve(target.pdf,add=TRUE,col="red",lwd=2)

T=matrix(0,nrow=m,ncol=n)
delta=matrix(0,nrow=m,ncol=n)
c=matrix(rexp(n*m,0.15),nrow=m)
phicap=NULL
for(i in 1:m)
{
	for(j in 1:n)
	{
		T[i,j]=min(Y[i,j],c[i,j])
		if(Y[i,j]<c[i,j]) delta[i,j]=1
		else delta[i,j]=0
	}
	phicap[i]=(n-sum(delta[i,]))/sum(T[i,])
}

un.l=NULL
cen.l=NULL
likelihood=function(a,l,phi,data,delta,cov)
{
   for(i in 1:n)
   {
   if(delta[i]==1)
   un.l[i]=log(a/(1+cov[i]))-((a/(1+cov[i]))*data[i])+log(shi.c(a/(1+cov[i]),l[i],data[i]))-(phi*data[i])
   else
     {
	 integrand=function(x)
	   {
	    (a/(1+cov[i]))*exp(-((a/(1+cov[i]))*x))*shi.c(a/(1+cov[i]),l[i],x)
	   }
	 cen.l[i]=log(integrate(integrand,lower=data[i],upper=Inf)$value)+log(phi)-phi*data[i]
     }
   }
like=sum(na.omit(un.l))+sum(na.omit(cen.l))
return(like)
}

alphacap=sapply(1:m, function(i) {
  current.data=T[i, ]
  current.lambda=lambda.u.cap.samp[i,]
  current.phi=phicap[i]
  current.delta=delta[i,]
  current.cov=u[i,]
  opt.alpha=optimize(f=likelihood,interval=c(0.01,2),maximum=TRUE,l=current.lambda,phi=current.phi,data=current.data,delta=current.delta,cov=current.cov)
  return(opt.alpha$maximum)
})

mse.alpha=(mean(alphacap)-a)^2+var(alphacap)
mse.lambda=(mean(lambda.u.cap)-l)^2+var(lambda.u.cap)
mse.phi=(mean(phicap)-0.15)^2+var(phicap)

inform.matrix.inv=matrix(c(n*var(alphacap),n*cov(alphacap,lambda.u.cap),n*cov(alphacap,phicap),n*cov(alphacap,lambda.u.cap),n*var(lambda.u.cap),n*cov(lambda.u.cap,phicap),n*cov(alphacap,phicap),n*cov(lambda.u.cap,phicap),n*var(phicap)),byrow=T,nrow=3)
inform.matrix=solve(inform.matrix.inv)
Q=NULL
for(i in 1:m)
Q[i]=n*t(c(alphacap[i]-a,lambda.u.cap[i]-l,phicap[i]-0.15))%*%inform.matrix%*%c(alphacap[i]-a,lambda.u.cap[i]-l,phicap[i]-0.15)
theoretical.quantiles=qchisq(ppoints(m),df=3)
observed.quantiles=sort(Q)
plot(theoretical.quantiles,observed.quantiles,main="Chi Square Q-Q plot",xlab="Theoretical Quantiles",ylab="Sample Quantiles",pch=19,col="blue",ylim=c(0,50))
abline(a=0,b=1,col="red",lwd=2)
ks.result=ks.test(Q,"pchisq",df=3)
cor(theoretical.quantiles,observed.quantiles)

alphacap.mean=mean(alphacap)
lambdacap.mean=mean(lambda.u.cap)
phicap.mean=mean(phicap)

surv=NULL
surv.prob=function(a,l,t,u) #--P(YN>t|u)
{
integrand=function(x)
  {
    (a/(1+u))*exp(-(a/(1+u))*x)*shi.c(a/(1+u),l*(A^u),x)
  }
for(i in 1:length(t))
surv[i]=integrate(integrand,lower=t[i],upper=Inf)
return(surv)
}
t=seq(0,60,0.01)
sur1=surv.prob(alphacap.mean,lambdacap.mean,t,u=0)
sur2=surv.prob(alphacap.mean,lambdacap.mean,t,u=1)
sur3=surv.prob(alphacap.mean,lambdacap.mean,t,u=2)
hazard=function(a,l,t,u) #--h(t) of YN 
{
f.y=(a/(1+u))*exp(-(a/(1+u))*t)*shi.c(a/(1+u),l*(A^u),t)
surv=surv.prob(a,l,t,u)
return(f.y/unlist(surv))
}
hazard1=hazard(alphacap.mean,lambdacap.mean,t,u=0)
hazard2=hazard(alphacap.mean,lambdacap.mean,t,u=1)
hazard3=hazard(alphacap.mean,lambdacap.mean,t,u=2)
surv.T=function(a,l,t,u)
{
integrand=function(x)
  {
    (a/(1+u))*exp(-(a/(1+u))*x)*shi.c(a/(1+u),l*(A^u),x)
  }
for(i in 1:length(t))
surv[i]=(integrate(integrand,lower=t[i],upper=Inf)$value)*exp(-(phicap.mean*t[i]))
return(surv)
}
sur1.T=surv.T(alphacap.mean,lambdacap.mean,t,u=0)
sur2.T=surv.T(alphacap.mean,lambdacap.mean,t,u=1)
sur3.T=surv.T(alphacap.mean,lambdacap.mean,t,u=2)

haz1.T=hazard1+phicap.mean
haz2.T=hazard2+phicap.mean
haz3.T=hazard3+phicap.mean

C.star=matrix(rztpois(n*m,k),nrow=m)
N.samp=matrix(0,nrow=m,ncol=n)
for(i in 1:m)
{
	for(j in 1:n)
	{
		N.samp[i,j]=rztpois(1,lambda.u.cap.samp[i,j])
	}
}
nu=matrix(0,nrow=m,ncol=n)
T.star=matrix(0,nrow=m,ncol=n)
for(i in 1:m)
{
	for(j in 1:n)
	{
	T.star[i,j]=min(N.samp[i,j],C.star[i,j])
	if(N.samp[i,j]<C.star[i,j]) nu[i,j]=1
	else nu[i,j]=0
	}
}

like=NULL
likelihood2=function(l,k,data,nu)
{
	for(i in 1:n)
	{
		if(nu[i]==1)
		like[i]=log(dztpois(data[i],l[i]))+log(pztpois(data[i],k,lower.tail=FALSE))
		else 
		like[i]=log(pztpois(data[i],l[i],lower.tail=FALSE)+dztpois(data[i],l[i]))+log(dztpois(data[i],k))
	}
return(sum(like))
}

k.cap=sapply(1:m, function(i) {
  current.data=T.star[i, ]
  current.lambda=lambda.u.cap.samp[i,]
  current.nu=nu[i,]
  opt.k=optimize(f=likelihood2,interval=c(0.01,5),maximum=TRUE,l=current.lambda,data=current.data,nu=current.nu)
  return(opt.k$maximum)
})

mse.lambda=(mean(lambda.u.cap)-l)^2+var(lambda.u.cap)
mse.k=(mean(k.cap)-k)^2+var(k.cap)

lambdacap.mean=mean(lambda.u.cap)
kcap.mean=mean(k.cap)

surv.N=function(t,u)
{
	lambda.u.hat=lambdacap.mean*(A^u)
	sur=pztpois(t,lambda.u.hat, lower.tail = FALSE)
	return(sur)
}
sur1.N=surv.N(t,u=0)
sur2.N=surv.N(t,u=1)
sur3.N=surv.N(t,u=2)

surv.T.star=function(t,u)
{
	lambda.u.hat=lambdacap.mean*(A^u)
	sur=pztpois(t,lambda.u.hat, lower.tail = FALSE)*pztpois(t,kcap.mean, lower.tail = FALSE)
	return(sur)
}

sur1.T.star=surv.T.star(t,u=0)
sur2.T.star=surv.T.star(t,u=1)
sur3.T.star=surv.T.star(t,u=2)
par(mfrow=c(1,3))
plot(t,sur1.T.star,col="darkblue",main="Survival Probability for u=0",xlab="t",ylab="S(t)",lwd=2)
lines(t,sur1.T,col="red",lwd=2)
legend(30,1,c("T*|u","T|u"),col=c("darkblue","red"),lty=1:1,bty="n")
plot(t,sur2.T.star,col="darkblue",main="Survival Probability for u=1",xlab="t",ylab="S(t)",lwd=2)
lines(t,sur2.T,col="red",lwd=2)
legend(30,1,c("T*|u","T|u"),col=c("darkblue","red"),lty=1:1,bty="n")
plot(t,sur3.T.star,col="darkblue",main="Survival Probability for u=2",xlab="t",ylab="S(t)",lwd=2)
lines(t,sur3.T,col="red",lwd=2)
legend(30,1,c("T*|u","T|u"),col=c("darkblue","red"),lty=1:1,bty="n")

hazard.N=function(t,u)
{
	lambda.u.hat=lambdacap.mean*(A^u)
	haz=dztpois(floor(t)+1,lambda.u.hat)/(pztpois(floor(t)+1,lambda.u.hat, lower.tail =FALSE)+dztpois(floor(t)+1,lambda.u.hat))
	return(haz)
}
haz1.N=hazard.N(t,u=0)
haz2.N=hazard.N(t,u=1)
haz3.N=hazard.N(t,u=2)

hazard.C.star=function(t,u)
{
	haz=dztpois(floor(t)+1,kcap.mean)/(pztpois(floor(t)+1,kcap.mean, lower.tail =FALSE)+dztpois(floor(t)+1,kcap.mean))
	return(haz)
}
haz1.C.star=hazard.C.star(t,u=0)
haz2.C.star=hazard.C.star(t,u=1)
haz3.C.star=hazard.C.star(t,u=2)

haz1.T.star=haz1.N+haz1.C.star
haz2.T.star=haz2.N+haz2.C.star
haz3.T.star=haz3.N+haz3.C.star

par(mfrow=c(1,3))
plot(t,haz1.T.star,col="darkblue",main="Hazard rate for u=0",xlab="t",ylab="h(t)",ylim=c(0,max(haz1.T)))
lines(t,haz1.T,col="red",lwd=0.5)
legend("topleft",legend=c("T*|u","T|u"),col=c("darkblue","red"),lty=1:1,bty="n")
plot(t,haz2.T.star,col="darkblue",main="Hazard rate for u=1",xlab="t",ylab="h(t)",ylim=c(0,max(haz2.T.star)+0.3))
lines(t,haz2.T,col="red",lwd=0.5)
legend("topleft",legend=c("T*|u","T|u"),col=c("darkblue","red"),lty=1:1,bty="n")
plot(t,haz3.T.star,col="darkblue",main="Hazard rate for u=2",xlab="t",ylab="h(t)",ylim=c(0,max(haz2.T.star)+0.2))
lines(t,haz3.T,col="red",lwd=0.5)
legend("topleft",legend=c("T*|u","T|u"),col=c("darkblue","red"),lty=1:1,bty="n")

par(mfrow=c(1,2))
plot(t,sur1.T,type="l",col="red",main="Survival Probability of T",xlab="t",ylab="S(t)",lwd=2)
lines(t,sur2.T,type="l",col="blue")
lines(t,sur3.T,type="l",col="green")
legend("topright", legend = c("u=0", "u=1","u=2"),col = c("red","blue","green"),lty = 1,bty = "n" )                           
plot(t,sur1.T.star,col="red",main="Survival Probability of T*",xlab="t",ylab="S(t)",type="s")
lines(t,sur2.T.star,type="s",col="blue")
lines(t,sur3.T.star,type="s",col="green")
legend("topright", legend = c("u=0", "u=1","u=2"),col = c("red","blue","green"),lty = 1,bty = "n" ) 

#PARAMETRIC TEST-------------------------------------------------

S1.hat.l=log(sur1.T/(1-sur1.T))
S2.hat.l=log(sur1.T.star/(1-sur1.T.star))
S1.hat.m=log(sur2.T/(1-sur2.T))
S2.hat.m=log(sur2.T.star/(1-sur2.T.star))
S1.hat.h=log(sur3.T/(1-sur3.T))
S2.hat.h=log(sur3.T.star/(1-sur3.T.star))
S.hat.l=(S1.hat.l+S2.hat.l)/2
S.hat.m=(S1.hat.m+S2.hat.m)/2
S.hat.h=(S1.hat.h+S2.hat.h)/2

num.obs=length(S.hat.l)

z.ho.l=(S1.hat.l-S2.hat.l)/(2^0.5/(num.obs*S.hat.l*(1-S.hat.l))^0.5)
z.ho.m=(S1.hat.m-S2.hat.m)/(2^0.5/(num.obs*S.hat.m*(1-S.hat.m))^0.5)
z.ho.h=(S1.hat.h-S2.hat.h)/(2^0.5/(num.obs*S.hat.h*(1-S.hat.h))^0.5)

na.omit(z.ho.l)>qnorm(1-0.05)
na.omit(z.ho.m)>qnorm(1-0.05)
na.omit(z.ho.h)>qnorm(1-0.05)

#KAPLAN MEIER----------------------------------------------------------
df=data.frame(data=c(T[which(u==0)],T[which(u==1)],T[which(u==2)]),indicator=c(delta[which(u==0)],delta[which(u==1)],delta[which(u==2)]),immunity=c(rep(0,sum(u==0)),rep(1,sum(u==1)),rep(2,sum(u==2))))
km.fit=survfit(Surv(data,indicator)~immunity,data=df)
df.low=subset(df,immunity==0)
df.med=subset(df,immunity==1)
df.hi=subset(df,immunity==2)
km.fit.low=survfit(Surv(data,indicator)~immunity,data=df.low)
km.fit.med=survfit(Surv(data,indicator)~immunity,data=df.med)
km.fit.hi=survfit(Surv(data,indicator)~immunity,data=df.hi)

df.N=data.frame(data=c(T.star[which(u==0)],T.star[which(u==1)],T.star[which(u==2)]),indicator=c(nu[whic
km.fit.low.N=survfit(Surv(data,indicator)~immunity,data=df.N.low)
km.fit.med.N=survfit(Surv(data,indicator)~immunity,data=df.N.med)
km.fit.hi.N=survfit(Surv(data,indicator)~immunity,data=df.N.hi)
par(mfrow=c(1,3))
plot(km.fit.low,conf.int=FALSE,col="red",xlab="t",ylab="S(t)",main="Survival Probability for u=0")
lines(km.fit.low.N,conf.int=FALSE,col="blue")
legend("topright",legend=c("T|u","T*|u"),col=c("red","blue"),lty=1:1,bty="n")
plot(km.fit.med,conf.int=FALSE,col="red",xlab="t",ylab="S(t)",main="Survival Probability for u=1")
lines(km.fit.med.N,conf.int=FALSE,col="blue")
legend("topright",legend=c("T|u","T*|u"),col=c("red","blue"),lty=1:1,bty="n")
plot(km.fit.hi,conf.int=FALSE,col="red",xlab="t",ylab="S(t)",main="Survival Probability for u=2")
lines(km.fit.hi.N,conf.int=FALSE,col="blue")
legend("topright",legend=c("T|u","T*|u"),col=c("red","blue"),lty=1:1,bty="n")
#COX PH MODEL----------------------------------------------------------
cox.fit=coxph(Surv(data,indicator) ~immunity, data = df)
summary(cox.fit)
library(muhaz)
base.haz=muhaz(df.low$data, df.low$indicator)
plot(base.haz$est.grid,base.haz$haz.est, main = "Hazard rate of T",xlab = "t", ylab = "h(t)",type="l",ylim=c(0,max(base.haz$haz.est)),lwd=2)
haz.1=base.haz$haz.est*exp(cox.fit$coef)
lines(base.haz$est.grid,haz.1,col="red",type="l",lwd=2)
haz.2=base.haz$haz.est*exp(cox.fit$coef*2)
lines(base.haz$est.grid,haz.2,col="blue",lwd=2)
legend("topleft",legend=c("u=0","u=1","u=2"),col=c("black","red","blue"),lty=1:1:1,lwd=2:2:2,bty="n")

cox.fit.N=coxph(Surv(data,indicator) ~immunity, data = df.N)
summary(cox.fit.N)
base.haz.N=muhaz(df.N.low$data, df.N.low$indicator)
plot(base.haz.N$est.grid,base.haz.N$haz.est, main = "Hazard rate of T*",xlab = "t", ylab = "h(t)",type="l",ylim=c(0,max(base.haz.N$haz.est)),lwd=2)
haz.1.N=base.haz.N$haz.est*exp(cox.fit.N$coef)
lines(base.haz.N$est.grid,haz.1.N,col="red",type="l",lwd=2)
haz.2.N=base.haz.N$haz.est*exp(cox.fit.N$coef*2)
lines(base.haz.N$est.grid,haz.2.N,col="blue",lwd=2)
legend("topleft",legend=c("u=0","u=1","u=2"),col=c("black","red","blue"),lty=1:1:1,lwd=2:2:2,bty="n")

par(mfrow=c(1,3))
plot(base.haz$est.grid,base.haz$haz.est, main = "Hazard rate for u=0",xlab = "t",col="red", ylab = "h(t)",type="l",ylim=c(0,0.6),lwd=2)
lines(base.haz.N$est.grid,base.haz.N$haz.est,col="blue")
legend("topright",legend=c("T|u","T*|u"),col=c("red","blue"),lty=1:1,lwd=2:2,bty="n")
plot(base.haz$est.grid,haz.1, main = "Hazard rate for u=1",xlab = "t",col="red", ylab = "h(t)",type="l",ylim=c(0,0.6),lwd=2)
lines(base.haz.N$est.grid,haz.1.N,col="blue")
legend("topright",legend=c("T|u","T*|u"),col=c("red","blue"),lty=1:1,lwd=2:2,bty="n")
plot(base.haz$est.grid,haz.2, main = "Hazard rate for u=2",xlab = "t",col="red", ylab = "h(t)",type="l",ylim=c(0,0.6),lwd=2)
lines(base.haz.N$est.grid,haz.2.N,col="blue")
legend("topright",legend=c("T|u","T*|u"),col=c("red","blue"),lty=1:1,lwd=2:2,bty="n")

#NON PARAMETRIC TEST------------------------------------------------

all.data=data.frame(data=c(T,T.star),indicators=c(delta,nu),immunity=rep(u,2),group=c(rep(1,100),rep(2,100)))
all.data.low=subset(all.data,immunity==0)
all.data.med=subset(all.data,immunity==1)
all.data.hi=subset(all.data,immunity==2)
lgr.low=survdiff(Surv(data,indicators)~group,data=all.data.low,rho=0)
lgr.med=survdiff(Surv(data,indicators)~group,data=all.data.med,rho=0)
lgr.hi=survdiff(Surv(data,indicators)~group,data=all.data.hi,rho=0)


