#adapted from studyPrior package
power_par = function(dat,n,dataPar,priorPars,lik) {
  if(lik=="normal")
  {
    sD = dataPar^2/n
    d = priorPars[2]^2 / (pmax((dat - priorPars[1]) ^ 2, sD + priorPars[2]^2) - sD)
  }
  if(lik=="binomial")
  {
    optfn=function(datC)
    {
      lik.d = function(d) VGAM::dbetabinom.ab(x=datC, size=n, d * (priorPars[1] -1) +1 ,  d *(priorPars[2]-1) +1) #uniform historical prior
      
      opd = BB::spg(par = .005,
                    fn = lik.d,
                    lower=0,
                    upper=1,
                    control=list(maximize=TRUE,
                                 trace=FALSE))
      
      if(opd$convergence!=0) print(opd)
      
      ds = opd$par
      
      return(ds)
    }
    d=sapply(dat,optfn)
  }
  return(d)
}



testDecision= function(dat,prior,alpha,n,n0D,dataPar,priorPars,weights,lik,a0,th0) 
{
  if (prior=="vague" | prior=="info" | prior=="sampling")
  {
    if(lik=="normal")
    {
      post.var = 1/(1/priorPars[2]^2+n/dataPar^2)
      post.mean= post.var*(priorPars[1]/(priorPars[2]^2)+dat/(dataPar^2/n))
      
      dec = (pnorm(th0,post.mean,sqrt(post.var))<=alpha)
    }
    if(lik=="binomial")
    {
      post.a=priorPars[1]+dat
      post.b=priorPars[2]+(n-dat)
      
      dec= (pbeta(th0,post.a,post.b)<=alpha)
    }
  }
  if (prior=="mix") #adapted from RBesT package
  {
    if(lik=="normal")
    {
      post.var = 1/(1/priorPars[2,]^2+n/dataPar^2)
      post.mean = matrix(NA,length(dat),2)
      post.mean[,1]=post.var[1]*(priorPars[1,1]/(priorPars[2,1]^2)+dat/(dataPar^2/n))
      post.mean[,2]=post.var[2]*(priorPars[1,2]/(priorPars[2,2]^2)+dat/(dataPar^2/n))
      dataParPred=matrix(NA,length(dat),2)
      dataParPred[,1] = sqrt(priorPars[2,1]^2 + dataPar^2/n)
      dataParPred[,2] = sqrt(priorPars[2,2]^2 + dataPar^2/n)
      
      margT= exp(log(weights[1]) + dnorm(dat, priorPars[1,1], dataParPred[,1], log=TRUE)) + exp(log(weights[2]) + dnorm(dat, priorPars[1,2], dataParPred[,2], log=TRUE))
      marg1= exp(log(weights[1]) + dnorm(dat, priorPars[1,1], dataParPred[,1], log=TRUE))
      post.weights=cbind(exp(log(marg1)-log(margT)),1-exp(log(marg1)-log(margT)))
      pMix= (1 - post.weights[,2]) * pnorm(th0, post.mean[,1], rep(sqrt(post.var[1]),length(dat))) + post.weights[,2] * pnorm(th0, post.mean[,2], rep(sqrt(post.var[2]),length(dat)))
      dec=(pMix<=alpha)
    }
    if(lik=="binomial")
    {
      post.a=matrix(NA,length(dat),2)
      post.a[,1]=priorPars[1,1]+dat
      post.a[,2]=priorPars[1,2]+dat
      post.b=matrix(NA,length(dat),2)
      post.b[,1]=priorPars[2,1]+(n-dat)
      post.b[,2]=priorPars[2,2]+(n-dat)
      
      marg1=log(weights[1]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,1], shape2=priorPars[2,1], log=TRUE)
      marg2=log(weights[2]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,2], shape2=priorPars[2,2], log=TRUE)
      
      margT= exp(marg1) + exp(marg2)
      post.weights=cbind(exp(marg1-log(margT)),1-exp(marg1-log(margT)))
      pMix= (1 - post.weights[,2]) * pbeta(th0, post.a[,1], post.b[,1]) + post.weights[,2] * pbeta(th0, post.a[,2], post.b[,2])
      dec=(pMix<=alpha)
    }
  }
  if (prior=="power")
  {
    if(lik=="normal")
    {
      prior.sd.EB=dataPar/sqrt(a0*n0D)
      post.var=  1/(1/prior.sd.EB^2+n/dataPar^2)
      post.mean= post.var*(priorPars[1]/(prior.sd.EB^2)+dat/(dataPar^2/n))
      dec=(pnorm(th0,post.mean,sqrt(post.var))<=alpha)
    }
    if(lik=="binomial")
    {
      post.a =  a0 * (priorPars[1] -1) +1 + dat
      post.b =  a0 *(priorPars[2]-1) +1 +(n-dat)
      dec= (pbeta(th0,post.a,post.b)<=alpha)
    }
  }
  return(dec)
}

priorProb= function(prior,n0D,dataPar,priorPars,a0,weights,lik,th0) 
{
  if (prior=="vague" | prior=="info" | prior=="sampling")
  {
    if(lik=="normal")
    {
      pNull= pnorm(th0, priorPars[1], priorPars[2])
    }
    if(lik=="binomial")
    {
      pNull= pbeta(th0, priorPars[1], priorPars[2])
    }
  }
  if (prior=="mix")
  {
    if(lik=="normal")
    {
      pNull= (1 - weights[2]) * pnorm(th0, priorPars[1,1],priorPars[2,1]) + weights[2] * pnorm(th0, priorPars[1,2], priorPars[2,2])
    }
    if(lik=="binomial")
    {
      pNull= (1 - weights[2]) * pbeta(th0, priorPars[1,1], priorPars[2,1]) + weights[2] * pbeta(th0, priorPars[1,2], priorPars[2,2])
    }
  }
  if (prior=="power")
  {
    if(lik=="normal")
    {
      prior.sd.EB=dataPar/sqrt(a0*n0D)
      pNull=pnorm(th0,priorPars[1], prior.sd.EB)
    }
    if(lik=="binomial")
    {
      prior.a.EB =  a0 * (priorPars[1] -1) +1 
      prior.b.EB =  a0 *(priorPars[2]-1) +1 
      pNull=pbeta(th0,prior.a.EB,prior.b.EB)
    }
  }
  return(pNull)
}

postMean= function(dat,prior,n,n0D,dataPar,priorPars,a0,weights,lik) {
  
  if (prior=="vague" | prior=="info" | prior=="sampling")
  {
    if(lik=="normal")
    {
      gamma=dataPar^2/priorPars[2]^2
      post.m=(gamma*priorPars[1]+n*dat)/(gamma+n)
    }
    if(lik=="binomial")
    {
      post.a=priorPars[1]+dat
      post.b=priorPars[2]+(n-dat)
      post.m=(post.a)/(post.a+post.b)
    }
  }
  
  if (prior=="mix") #adapted from RBesT package
  {
    if(lik=="normal")
    {
      post.var = 1/(1/priorPars[2,]^2+n/dataPar^2)
      post.mean = matrix(NA,length(dat),2)
      post.mean[,1]=post.var[1]*(priorPars[1,1]/(priorPars[2,1]^2)+dat/(dataPar^2/n))
      post.mean[,2]=post.var[2]*(priorPars[1,2]/(priorPars[2,2]^2)+dat/(dataPar^2/n))
      dataParPred=matrix(NA,length(dat),2)
      dataParPred[,1] = sqrt(priorPars[2,1]^2 + dataPar^2/n)
      dataParPred[,2] = sqrt(priorPars[2,2]^2 + dataPar^2/n)
      
      margT= exp(log(weights[1]) + dnorm(dat, priorPars[1,1], dataParPred[,1], log=TRUE)) + exp(log(weights[2]) + dnorm(dat, priorPars[1,2], dataParPred[,2], log=TRUE))
      marg1= exp(log(weights[1]) + dnorm(dat, priorPars[1,1], dataParPred[,1], log=TRUE))
      post.weights=cbind(exp(log(marg1)-log(margT)),1-exp(log(marg1)-log(margT)))
      
      post.m=post.weights[,1]*post.mean[,1]+post.weights[,2]*post.mean[,2]
    }
    if(lik=="binomial")
    {
      post.a=matrix(NA,length(dat),2)
      post.a[,1]=priorPars[1,1]+dat
      post.a[,2]=priorPars[1,2]+dat
      post.b=matrix(NA,length(dat),2)
      post.b[,1]=priorPars[2,1]+(n-dat)
      post.b[,2]=priorPars[2,2]+(n-dat)
      
      post.mean=(post.a)/(post.a+post.b)
      
      marg1=log(weights[1]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,1], shape2=priorPars[2,1], log=TRUE)
      marg2=log(weights[2]) + VGAM::dbetabinom.ab(x=dat, size=n, shape1=priorPars[1,2], shape2=priorPars[2,2], log=TRUE)
      
      margT= exp(marg1) + exp(marg2)
      post.weights=cbind(exp(marg1-log(margT)),1-exp(marg1-log(margT)))
      
      post.m=post.weights[,1]*post.mean[,1]+post.weights[,2]*post.mean[,2]
    }
  }
  
  if (prior=="power")
  {
    if(lik=="normal")
    {
      sD=dataPar^2/n
      prior.sd.EB=dataPar/sqrt(a0*n0D)
      post.var=  1/(1/prior.sd.EB^2+n/dataPar^2)
      post.m= post.var*(priorPars[1]/(prior.sd.EB^2)+dat/(dataPar^2/n))
    }
    if(lik=="binomial")
    {
      post.a =  a0 * (priorPars[1] -1) +1 + dat
      post.b =  a0 *(priorPars[2]-1) +1 +(n-dat)
      post.m=(post.a)/(post.a+post.b)
    }
  }
  return(post.m)
}



integrand_Err_MSE_relevance = function(prior=c("vague","info","mix","power","sampling"),n,n0D=NULL,c.alpha1,c.beta1,dataPar=NULL,
                                       priorParsWeights,lik=c('binomial','normal'),truth,decision=c('fixed','normalized'),
                                       alphaIII=0.05,betaIII=0.1,S=5,Fr=400,th0,thR)
{
  if(prior=="mix")
  {
    priorPars=priorParsWeights[2:3,]
    weights=priorParsWeights[1,] 
  }
  else
  {
    priorPars=c(priorParsWeights[2],priorParsWeights[3])
    weights=priorParsWeights[1]
  }
  
  if(lik=="normal")
  {
    #data integration grid
    grid=seq(0-S*dataPar/sqrt(n),0+S*dataPar/sqrt(n),by=dataPar/sqrt(n)/Fr)
    x1t=matrix(rep(grid,length(truth)),length(truth),length(grid), byrow = TRUE)
    x1m=x1t+truth
    x1=as.vector(t(x1m))
    wi=dnorm(grid,0,dataPar/sqrt(n))
    wiM=t(wi/sum(wi))
    
    #effect size
    delta=truth
    sigma=sqrt(2)*dataPar
  }
  
  if(lik=="binomial")
  {
    #data integration grid
    grid=seq(0,n)
    x1=rep(grid,length(truth))
    wi=dbinom(x1,n,rep(truth,each=(n+1)))
    wiM=matrix(wi,length(truth),n+1,byrow=TRUE)
    
    #effect size
    delta=truth-th0
    sigma=sqrt(th0*(1-th0)+truth*(1-truth))
  }
  
  #Phase III power and sample size under truth
  n_t=((qnorm(alphaIII)+qnorm(betaIII))^2)/(delta^2/sigma^2)
  power_t=1-betaIII
  
  #power parameter
  a0=rep(1,length(x1))
  if (prior=="power")
  {
    if(lik=="binomial")
    {
      a0g=power_par(dat=grid,n=n,dataPar=dataPar,priorPars=priorPars,lik=lik)
      a0=rep(a0g,length(truth))
    }
    if(lik=="normal")
    {
      a0=power_par(dat=x1,n=n,dataPar=dataPar,priorPars=priorPars,lik=lik)
    }
  }
  
  #decision threshold
  
  if (decision=="fixed")
  {
    c.alpha=c.alpha1
    c.beta=c.beta1
  }
  if (decision=="normalized")
  {
    p0=priorProb(prior=prior,n0D=n0D,dataPar=dataPar,priorPars=priorPars,a0=a0,weights=weights,lik=lik,th0=th0) 
    c.alpha=c.alpha1/p0
    c.beta=c.beta1/(1-p0)
  }
  
  alpha=c.beta/(c.alpha+c.beta)
  alphaM = matrix(alpha,length(grid),length(truth))
  
  #significance test decision
  if(lik=="binomial")
  {
    significance=testDecision(dat=x1,prior=prior,alpha=alpha,n=n,n0D=n0D,dataPar=dataPar,priorPars=priorPars,weights=weights,lik=lik,a0=a0,th0=th0)
  }
  if(lik=="normal")
  {
    significance= testDecision(dat=x1,prior=prior,alpha=alpha,n=n,n0D=n0D,dataPar=dataPar,priorPars=priorPars,weights=weights,lik=lik,a0=a0,th0=th0) 
  }
  significance_M=matrix(significance,length(grid),length(truth))
  
  #error rate and average alpha
  
  if(lik=="normal")
  {
    error= as.vector(wiM %*% significance_M)
    alphaO= as.vector(wiM %*% alphaM)
  }
  if(lik=="binomial")
  {
    error=diag(wiM %*% significance_M)
    alphaO=diag(wiM %*% alphaM)
  }
  
  errorO=ifelse(truth<=th0,error,(1-error))
  
  #quadratic loss
  post_mean= postMean(dat=x1,prior=prior,n=n,n0D=n0D,dataPar=dataPar,priorPars=priorPars,a0=a0,weights=weights,lik=lik)
  truthV=rep(truth,each=length(grid))
  qV=(post_mean-truthV)^2
  qM=matrix(qV,length(grid),length(truth), byrow = FALSE)
  
  #relevance
  if(lik=="binomial")
  {
    relevance = testDecision(dat=x1,prior=prior,alpha=0.5,n=n,n0D=n0D,dataPar=dataPar,priorPars=priorPars,weights=weights,lik=lik,a0=a0,th0=thR)
  }
  if(lik=="normal")
  {
    relevance=(post_mean>=thR)
  }
  
  #indeterminate
  indeterminateM = matrix(ifelse(significance+relevance==1,1,0),length(grid),length(truth), byrow = FALSE)
  
  
  #actual phase II power and sample size, prob of indeterminate decision and mse
  if(lik=="normal")
  {
    #estimated effect size and Phase III actual sample size
    deltaE=post_mean
    sigmaE=sqrt(2)*dataPar
    n_a=((qnorm(alphaIII)+qnorm(betaIII))^2)/(deltaE^2/sigmaE^2)
    sigmaEn=sqrt(2/n_a)*dataPar
    
    #probability of indeterminate and mse
    indeterminate_p= as.vector(wiM %*% indeterminateM)
    mse= as.vector(wiM %*% qM)
  }
  
  if(lik=="binomial")
  {
    #estimated effect size and Phase III actual sample size
    deltaE=post_mean-th0
    sigmaE=sqrt(th0*(1-th0)+post_mean*(1-post_mean))
    n_a=((qnorm(alphaIII)+qnorm(betaIII))^2)/(deltaE^2/sigmaE^2)
    sigmaEn=sqrt((1/n_a)*(th0*(1-th0)+post_mean*(1-post_mean)))
    
    #probability of indeterminate and mse
    indeterminate_p= diag(wiM %*% indeterminateM)
    mse= diag(wiM %*% qM)
  }
  
  #Phase III actual power
  power_a =pnorm(-qnorm(1-alphaIII)+(rep(delta,each=length(grid))/sigmaEn)) #only for positive delta
  
  #Power loss
  power_lossM= matrix((power_t-power_a*significance*relevance),length(grid),length(truth), byrow = FALSE)
  
  #Sample size gain
  n_gainM= matrix((-rep(n_t,each=length(grid))+n_a*significance*relevance),length(grid),length(truth), byrow = FALSE)
  
  #80% quantile
  ord_pow_loss=apply(power_lossM,2,order)
  ord_n_gain=apply(n_gainM,2,order)
  
  if(lik=="normal")
  {
    power_loss80i= apply(ord_pow_loss,2, function(x) {
    cs=cumsum(t(wiM)[x]) 
    ret=max(which(cs==max(cs[cs<=0.80])))
    x[ret]})
    
    n_gain80i= apply(ord_n_gain,2, function(x) {
      cs=cumsum(t(wiM)[x]) 
      ret=max(which(cs==max(cs[cs<=0.80])))
      x[ret]})
  }
  if(lik=="binomial")
  {
    power_loss80i= sapply(1:length(truth), function(x) {
      col=ord_pow_loss[,x]
      cs=cumsum(t(wiM[x,])[col]) 
      ret=ifelse(min(cs)>0.80,1,max(which(cs==max(cs[cs<=0.80]))))
      #ret=max(which(cs==max(cs[cs<=0.80])))
      col[ret]})
    
    n_gain80i= sapply(1:length(truth), function(x) {
      col=ord_n_gain[,x]
      cs=cumsum(t(wiM[x,])[col]) 
      ret=ifelse(min(cs)>0.80,1,max(which(cs==max(cs[cs<=0.80]))))
      #ret=max(which(cs==max(cs[cs<=0.80])))
      col[ret]})
  }
  power_loss80=power_lossM[cbind(power_loss80i,1:length(truth))]
  n_gain80=n_gainM[cbind(n_gain80i,1:length(truth))]
  
  
  return(cbind('error'=errorO,'alpha'=alphaO,'MSE'=mse,'power_loss'=power_loss80,'pr_indet'=indeterminate_p,'n_gain'=n_gain80))
  
}


require(truncnorm)
priorTr=function(th,prior=c("vague","info","mix","power","sampling"),priorPars,hp=c('null','alt'),weights=NULL,th0) 
{
  if(grepl("mix",prior))
  {
    if (hp=='null')
      out=(weights[1]*dnorm(th,priorPars[1,1],priorPars[2,1])+weights[2]*dnorm(th,priorPars[1,2],priorPars[2,2]))/integrate(function(g) (weights[1]*dnorm(g,priorPars[1,1],priorPars[2,1])+weights[2]*dnorm(g,priorPars[1,2],priorPars[2,2])),-Inf,th0)$value
    if (hp=='alt')
      out=(weights[1]*dnorm(th,priorPars[1,1],priorPars[2,1])+weights[2]*dnorm(th,priorPars[1,2],priorPars[2,2]))/integrate(function(g) (weights[1]*dnorm(g,priorPars[1,1],priorPars[2,1])+weights[2]*dnorm(g,priorPars[1,2],priorPars[2,2])),th0,Inf)$value
  }
  else
  {
    if (hp=='null')
      out= dtruncnorm(th,-Inf,th0,priorPars[1],priorPars[2])
    if (hp=='alt')
      out= dtruncnorm(th,th0,Inf,priorPars[1],priorPars[2])
  }
  out
}



BF=function(dat,n,dataPar,prior=c("vague","info","mix","power","sampling"),priorParsWeights,th0)
{
  if(grepl("mix",prior))
  {
    priorPars=priorParsWeights[2:3,]
    weights=priorParsWeights[1,] 
  }
  else
  {
    a0=ifelse(prior=='power',power_par(dat=dat,n=n,dataPar=dataPar,priorPars=c(priorParsWeights[2],priorParsWeights[3]),lik='normal'),1)
    
    priorPars=c(priorParsWeights[2],priorParsWeights[3]/sqrt(a0))
    weights=priorParsWeights[1]
  }
  
  num=integrate(function(x) (dnorm(dat,x,dataPar/sqrt(n))*priorTr(th=x,prior=prior,priorPars=priorPars,weights=weights,th0=th0,hp='null')),-Inf,th0)$value
  den=integrate(function(x) (dnorm(dat,x,dataPar/sqrt(n))*priorTr(th=x,prior=prior,priorPars=priorPars,weights=weights,th0=th0,hp='alt')),th0,Inf)$value
  
  out=num/den
  out
}
