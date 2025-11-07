#MB2306, MB2304, MB2323, MB2330
#Select everything and run at once, requires around 15 minutes to complete
#Numbers are printed on the screen: 1 to 31, 1 to 200 (thrice consecutively), 1 to 100 for tracking progress inside the for loops
#5 pdfs are directly saved in the device, change the directory (see the first line of the program) accrdingly
#Install the "Deriv" package if not already installed

setwd("C:/Users/User/Desktop/SLOE")  #Change this accordingly

library(Deriv)

set.seed(17)

f=function(y0,t0)
{
  if (y0==-1) {log(1+exp(t0))}
  else if (y0==1) {log(1+exp(-t0))}
}
d1=Deriv(f,"t0")
d2=Deriv(d1,"t0")

sloe=function(y,X)  #Input: responses and data matrix, output: eeta
{
  n=nrow(X)
  p=ncol(X)
  ynew=replace(y,(y==0),-1)
  betahat=as.matrix(as.numeric(glm(y~X-1,family="binomial")$coef))  #MLE
  tt=X%*%betahat
  myvec=exp(tt)
  myvec=myvec/((1+myvec)*(1+myvec))
  mat=matrix(0,p,p)
  for (i in 1:p)
  {
    for (j in 1:p)
    {
      mat[i,j]=sum(X[,i]*X[,j]*myvec)  #Hessian
    }
  }
  H0=solve(mat)
  w=qq=S=numeric(n)
  for (i in 1:n)    #The main steps of SLOE
  {
    w[i]=t(as.matrix(X[i,]))%*%H0%*%as.matrix(X[i,])
    qq[i]=w[i]/(1-w[i]*d2(ynew[i],tt[i]))
    S[i]=tt[i]+qq[i]*d1(ynew[i],tt[i])
  }
  eeta=sqrt(var(S)*(n-1)/n)
  return(eeta)
}

dat=function(n,p,m,k,repl,B,type)  #Sample size, No. of Parameters, t distribution df, no. of non-null parameters, bootstrap sample size, pareto/multivariate t
{
  sig=matrix(0,p,p)
  for (i in 1:nrow(sig))
  {
    for (j in 1:ncol(sig))
    {
      sig[i,j]=0.5^(min(abs(i-j),p+1-abs(i-j)))  #Covariance matrix
    }
  }
  eig=eigen(sig)
  pdroot=eig$vec%*%sqrt(diag(eig$val))%*%solve(eig$vec)  #Naive way of simulating from a normal distribution with a specified covariance matrix 
  v=rbinom(1,k,0.5)
  betaa=c(rnorm(v,5,1),rnorm((k-v),-5,1),rep(0,p-k))
  z1=pdroot%*%matrix(rnorm(n*p),p,n)
  z2=sqrt(rchisq(n,m)/(m-2))  #Standardized to make variances 1
  X=t(z1)/(z2*sqrt(p))        #Variance of each column is now 1/p
  if (type=="pareto") {X=matrix(1/(runif(n*p))^(1/5),n,p)}  #Covariates from standard Pareto (M=5,alpha=1) distribution
  Xbeta=X%*%betaa
  y=rbinom(n,1,plogis(Xbeta))
  betahat=glm(y~X-1,family=binomial(link="logit"))$coef
  Xbetahat=X%*%betahat
  myy <<- y
  myX <<- X
  betaa <<- betaa
  betahat <<- betahat
  Xbetahat <<- Xbetahat
  Xbeta <<- Xbeta
  eetatilde=sloe(myy,myX)    #Observed eeta
  
  #1. Eeta-Gamma curve (generates plot after printing 1 to 31)

  myseq=seq(0,1,length=31)    #To find the factor by which the MLE must be scaled down to estimate the actual coefficients
  gamm=numeric(length(myseq))
  J=5
  eetamat=matrix(0,J*length(myseq),2)
  for (i in 1:length(myseq))
  {
print(i)
    mybeta=myseq[i]*betahat
    gamm[i]=sd(Xbetahat)*sqrt((n-1)/n)*myseq[i]
    for (j in 1:J)
    {
      ytemp=rbinom(n,1,plogis(Xbetahat*myseq[i]))  #For each gamma, simulate responses (J times to reduce fluctuation) from the scaled MLE 
      eetamat[J*(i-1)+j,1]=gamm[i]
      eetamat[J*(i-1)+j,2]=sloe(ytemp,myX)
    }
  }

  pdf("gamma_eeta_plot.pdf") 
  plot(eetamat[,1],eetamat[,2],pch=16,xlab="Gamma",ylab="Eta",main="Smoothed LOESS fit of Eta on Gamma")  #A discrete gamma vs eeta plot  
  abline(h=eetatilde,col="red",lwd=2)
  temp0=predict(loess(eetamat[,2]~eetamat[,1]))
  lines(temp0,x=eetamat[,1],col="blue",lwd=2)    #LOESS Smoothed version of the above plot
  dev.off()
  gammahat=eetamat[which.min(abs(temp0-eetatilde)),1]  #Solving for gamma from the monotonically increasing LOESS fit, given the observed eeta

  #2. Scatterplot of the k non null model coefficients and their MLE coungterparts (averaged over M=200 experiments)
  #(generates plot after printing 1 to 200)

  tempsum=rep(0,p)
  M=200  
  arr=numeric(M)
  for (i in 1:M)
  {
print(i)
    y1=rbinom(n,1,plogis(Xbeta))
    tempsum0=glm(y1~myX-1,family=binomial(link="logit"))$coef  #200 MLEs calculated from the original model...benchmark for statistical techniques
    arr[i]=tempsum0[1]
    tempsum=tempsum+tempsum0
  }
  tempsum=tempsum/M

  pdf("normal_qq_plot.pdf")  #Q-Q plot of a coordinate of the averaged out MLE (vs normal quantiles)
  qqnorm(arr,pch=16,cex=0.7)
  qqline(arr,lwd=2,col="red") 
  dev.off()

  pdf("MLE vs actual.pdf")   #Plot of MLE vs actual (non-null) coefficients
  plot(betaa[1:k],tempsum[1:k],main="Average MLE Coefficients vs Actual (non-null) Coefficients",pch=16,xlab="Actual (non-null) coefficients",ylab="Average MLE coefficients")
  abline(lm(tempsum[1:k]~betaa[1:k]),col="red",lwd=2)
  dev.off()

  #3. Failure of Pairs and Parametric; success of Resized
  #(generates plot after printing 1 to 200 twice)

  coefmat1=coefmat2=coefmat3=numeric(repl)
  for (i in 1:repl)
  {
print(i)
    y1=rbinom(n,1,plogis(Xbeta))    #Responses simulated from the actual parameter
    y2=rbinom(n,1,plogis(Xbetahat)) #(Parametric Bootstrap) Responses simulated using MLE
    temp=sample(c(1:n),n,replace=TRUE)  #For Pairs Bootstrap
    result1 <<- glm(y1~X-1,family=binomial(link="logit"))
    result2 <<- glm(y2~X-1,family=binomial(link="logit"))
    result3 <<- glm(y[temp]~X[temp,]-1,family=binomial(link="logit"))
    coefmat1[i]=result1$coef[1]
    coefmat2[i]=result2$coef[1]
    coefmat3[i]=result3$coef[1]
  }
  coefmat1 <<- coefmat1
  coefmat2 <<- coefmat2
  coefmat3 <<- coefmat3

  #4. Bootstrap for estimation of location and scale

  s=myseq[(which.min(abs(temp0-eetatilde))-1)/J]  #The downscaling factor of the MLE we were after
  s <<- s
  betastar=s*betahat  #Betastar should be a good (unbiased) estimate of the actual model coefficients
  betaboot=matrix(0,B,p)
  for (i in 1:B)
  {
print(i)
    ystar=rbinom(n,1,plogis(Xbetahat*s))  #Responses generated from the (unbiasedly) estimated model coefficients
    betaboot[i,]=glm(ystar~X-1,family=binomial(link="logit"))$coef
  }
  sigmahat=numeric(p)
  for (i in 1:p) {sigmahat[i]=sd(betaboot[,i])}  #Estimated bootstrap sd of the MLE distribution
  alphahat=lm(tempsum[1:k]~betastar[1:k],weights=(1/sigmahat^2)[1:k])$coef[2]  #Regression of average MLE coefficients on (unbiasedly) estimated model coefficients
  muhat=betahat[1]/alphahat  #Estimated mean of the MLE distribution

  pdf("failure_and_success.pdf")  #Performance of parametric, paired, resized bootstrap techniques...benchmark being the actual MLE distribution
  plot(density(coefmat1),xlab="Estimated Coordinate of the Parameter Vector",main="Density of a Coordinate of the Estimated Paramater Vector",col="black",lwd=2,xlim=c(3,19),ylim=c(0,0.28))
  lines(density(coefmat2),col="blue",lwd=2)
  lines(density(coefmat3),col="green",lwd=2)
  curve(dnorm(1/sigmahat[1]*(x-muhat))/sigmahat[1],add=TRUE,col="red",lwd=2)
  legend(x="topright",legend=c("Parametric Bootstrap","Pairs Bootstrap","Resized Bootstrap","Exact Distribution of MLE"),lty=1,col=c("blue","green","red","black"),lwd=4,text.font=15,box.lwd=2,title="Procedure",inset=0.02,bg=rgb(1,1,0.5,0.6))
  dev.off()

  pdf("MLE_vs_bootstrapped_MLE.pdf")  #Q-Q plot of standardized MLE vs standardized bootstrapped MLE
  qqplot((arr-mean(arr))/sd(arr),(betaboot[,1]-mean(betaboot[,1]))/sd(betaboot[,1]),pch=16,cex=0.7,main="Q-Q Plot of Standardized MLE vs Standardized Bootstrapped MLE",xlab="Standardized MLE",ylab="Standardized Bootstrap MLE")
  qqline((arr-mean(arr))/sd(arr),lwd=2,col="red")
  dev.off()

  print(gammahat)
  print(sd(Xbeta))
  eetamat <<- eetamat
  alphahat <<- alphahat
  muhat <<- muhat
  sigmahat <<- sigmahat
  betaboot <<- betaboot
}
dat(2000,200,8,20,200,200,"t")
#dat(400,40,8,20,200,200,"pareto")

n=2000
p=200
boundg=boundt=matrix(0,100,2)
countg=countt=0
qa=quantile((betaboot[,1]-s*alphahat*betahat[1])/sigmahat[1],0.05/2,names=FALSE)
qb=quantile((betaboot[,1]-s*alphahat*betahat[1])/sigmahat[1],1-0.05/2,names=FALSE)
for (i in 1:100)  #100 CIs
{
print(i)  
  tempy=rbinom(n,1,plogis(Xbeta))
  tempmle=glm(tempy~myX-1,family=binomial(link="logit"))$coef[1]
  boundg[i,1]=(tempmle-qnorm(0.05/2)*sigmahat[1])/alphahat  #CI for Boot-g (Gaussian approximation)
  boundg[i,2]=(tempmle+qnorm(0.05/2)*sigmahat[1])/alphahat
  boundt[i,1]=(tempmle+qb*sigmahat[1])/alphahat             #CI for Boot-t (Empirical Bootstrap approximation)
  boundt[i,2]=(tempmle+qa*sigmahat[1])/alphahat
  if ((boundg[i,2]<=betaa[1]) && (betaa[1]<=boundg[i,1])) {countg=countg+1}
  if ((boundt[i,2]<=betaa[1]) && (betaa[1]<=boundt[i,1])) {countt=countt+1}
}
print(countg/100)  #Boot-g coverage probabliity
print(countt/100)  #Boot-t coverage probability

print(c(mean(coefmat2),sd(coefmat2)))  #Mean, SD for Parametric Bootstrap
print(c(mean(coefmat3),sd(coefmat3)))  #Mean, SD for Pairs Bootstrap
print(c(muhat,sigmahat[1]))            #Mean, SD for Resized Bootstrap
print(c(mean(coefmat1),sd(coefmat1)))  #Mean, SD for Actual MLE Coefficients

