library(parallel)

rm(list=ls())

# Design function ---------------------------------------------------------

designF = function(index,outcome,dec_cost,max_n,designGrid,th0,thR,prior.listI,n0=NULL,sigma=NULL)
{
  source('Design_Fn.R')
  
  if (outcome=="normal")
  { 
    #sampling prior parameters and pr of null hypothesis
    true.pars=c(designGrid[index,1],sqrt(designGrid[index,2]))
    pSamp=(pnorm(th0, true.pars[1], true.pars[2]))
    sampling = c(1, true.pars[1], true.pars[2])
    
    #parameters for integration grids
    #data
    sl=5
    frl=400
    #sampling prior
    Sh=10
    Frh=400
    
    #integration grid and weights
    yErr=seq(true.pars[1]-Sh*true.pars[2],true.pars[1]+Sh*true.pars[2],by=true.pars[2]/Frh)
    int_gr=list('t1'=yErr[yErr<=th0],'t2'=yErr[yErr>th0],'mse'=yErr,'rel'=yErr[yErr>thR])
    weights.oc=lapply(int_gr, function(x) dnorm(x,true.pars[1],true.pars[2])/sum(dnorm(x,true.pars[1],true.pars[2])))
  }
  
  if (outcome=="binomial")
  {
    #sampling prior parameters and pr of null hypothesis
    true.pars=designGrid[index,]
    pSamp=(pbeta(th0,true.pars[1], true.pars[2]))
    sampling = c(1, true.pars[1], true.pars[2])
    
    #integration grid and weights
    yErr=seq(0,1,0.00005)
    int_gr=list('t1'=yErr[yErr<=th0],'t2'=yErr[yErr>th0],'mse'=yErr,'rel'=yErr[yErr>thR])
    weights.oc=lapply(int_gr, function(x) dbeta(x,true.pars[1],true.pars[2])/sum(dbeta(x,true.pars[1],true.pars[2])))
  }
  
  
  #priors list
  prior.list = c(prior.listI,list("sampling"=sampling))
  
  #test error costs
  c.alpha1 = 0.95
  c.beta1 = 0.05
  
  
  if (dec_cost=="fixed")
  {
    c.alphaS=c.alpha1
    c.betaS=c.beta1
  }
  if (dec_cost=="normalized")
  {
    c.alphaS=c.alpha1/pSamp
    c.betaS=c.beta1/(1-pSamp)
  }
  
  #sample size grid
  n.list=1:max_n
  
  Out_List=lapply(1:length(prior.list),function(y) {sapply(n.list,function(x) 
  {
    conditionalMeas=integrand_Err_MSE_relevance(prior=names(prior.list)[y],
                                                n=x,n0D=n0,c.alpha1=c.alpha1,c.beta1=c.beta1,
                                                dataPar=sigma,priorParsWeights=prior.list[[y]],
                                                lik=outcome,truth=yErr,decision=dec_cost,alphaIII=0.05,betaIII=0.1,
                                                S=sl,Fr=frl,th0=th0,thR=thR)
    
    out=c('type1'=weights.oc$t1 %*% conditionalMeas[yErr<=th0,'error'],
          'type2'=weights.oc$t2 %*% conditionalMeas[yErr>th0,'error'],
          'gamma'=weights.oc$mse %*% conditionalMeas[,'alpha'],
          'ErrSum'=  c.alphaS * (weights.oc$t1 %*% conditionalMeas[yErr<=th0,'error']) *pSamp + c.betaS * (weights.oc$t2 %*% conditionalMeas[yErr>th0,'error']) *(1-pSamp),
          'MSE'=weights.oc$mse %*% conditionalMeas[,'MSE'],
          'PowLoss'=weights.oc$rel %*% conditionalMeas[yErr>thR,'power_loss'],
          'PrIndet'=weights.oc$rel %*% conditionalMeas[yErr>thR,'pr_indet'],
          'NGain'=weights.oc$rel %*% conditionalMeas[yErr>thR,'n_gain'])
    out
  })})
  names(Out_List)= names(prior.list)
  
  cat(index)
  return(c(Out_List,list('Sampling_prior_params'=true.pars)))
  
}


# Run design - Simulations --------------------------------------------------------------

designRunSim = function(outcome=c('binomial','normal'), #type of data outcome
                        dec_cost=c('fixed','normalized'),  #type of test error cost: normalized (costs include prior probabilities of hypotheses) or fixed (costs fixed across prior specifications)
                        max_n, #maximum sample size at which the operating characteristics are computed
                        nCores=1)
{
  outcome=match.arg(outcome)
  dec_cost=match.arg(dec_cost)
  
  if (outcome=='normal')
  {
    #design grid
    true.m.vec=seq(-0.9,1.25,by=0.05)
    true.v.vec=0.02
    designGrid=cbind(rep(true.m.vec,each=length(true.v.vec)),rep(true.v.vec,length(true.m.vec)))
    
    #hypothesis threshold
    th0=0
    
    #relavance threshold
    thR=0.15
    
    #data standard deviation
    sigma = 1 
    
    #historical dat - fitting prior
    n0 = 50
    prior.m = 0.25
    s0 = sqrt(sigma^2/n0)  #historical prior sd
    
    #priors
    info = c(1,prior.m,s0)
    base = c(1,0,10)
    mix = cbind(c(0.5, prior.m, s0), c(0.5,prior.m,10))
    power = c(1, prior.m, s0)
    prior.listI =list("vague"=base, "info"=info, "mix"=mix, "power"=power)
  }
  
  if (outcome=='binomial')
  {
    #design grid
    true.a.vec=seq(1,49,by=2)
    true.b.vec=50-true.a.vec
    
    true.p.vec=true.a.vec/(true.a.vec+true.b.vec)
    designGrid=cbind(true.a.vec,true.b.vec)
    
    #hypothesis threshold
    th0=0.2
    
    #relavance threshold
    thR=0.4
    
    #historical dat - fitting prior
    prior.a = 25
    prior.b = 25
    
    #priors
    info = c(1,prior.a,prior.b)
    base = c(1,1,1)
    mix = cbind(c(0.5, prior.a,prior.b), c(0.5, 1, 1))
    power = c(1,prior.a,prior.b) 
    prior.listI = list("vague"=base, "info"=info, "mix"=mix, "power"=power)
  }
  
  id=1:dim(designGrid)[1]
  desOutF_fixed=mclapply(id, designF, outcome=outcome,dec_cost=dec_cost,max_n=max_n,designGrid=designGrid,th0=th0,thR=thR,prior.listI=prior.listI,n0=n0,sigma=sigma, mc.cores = nCores)
  
  save(desOutF_fixed, file= paste("desOutF_",outcome,"_",dec_cost,".RData",sep=""))
}


# Run design - Data example --------------------------------------------------------------

designRunDataEx = function(dec_cost=c('fixed','normalized'),  #type of test error cost: normalized (costs include prior probabilities of hypotheses) or fixed (costs fixed across prior specifications)
                           max_n, #maximum sample size at which the operating characteristics are computed
                           nCores=1)
{
  outcome='binomial'
  dec_cost=match.arg(dec_cost)
  
  #design grid
  true.a.vec=seq(1,39,by=2)
  true.b.vec=40-true.a.vec
  
  true.p.vec=true.a.vec/(true.a.vec+true.b.vec)
  designGrid=cbind(true.a.vec,true.b.vec)
  
  #hypothesis threshold
  th0=0.075
  
  #relavance threshold
  thR=0.175
  
  #priors
  info = c(1,11,29)
  vague = c(1,0.0811,1)
  mix = cbind(c(0.5, 11,29), c(0.5, 1, 1))
  power = c(1,11,29) 
  prior.listI = list("vague"=info,"mix"=mix,"info"=info,"power"=power)
  
  id=1:dim(designGrid)[1]
  desOutF_fixed=mclapply(id, designF, outcome=outcome,dec_cost=dec_cost,max_n=max_n,designGrid=designGrid,th0=th0,thR=thR,prior.listI=prior.listI,mc.cores = nCores)
  
  save(desOutF_fixed, file= "desOutF_Example.RData")
}


# Run design - Simulations (larger maximum sample size, 800, for Figure 2 and S1) --------------------------------------------------------------

designRunSimLSS = function(outcome='normal',
                        dec_cost='normalized',  
                        max_n= 800)
{
  
    #design grid
    true.m.vec=seq(-0.9,1.25,by=0.05)
    true.v.vec=0.02
    designGrid=cbind(rep(true.m.vec,each=length(true.v.vec)),rep(true.v.vec,length(true.m.vec)))
    
    #hypothesis threshold
    th0=0
    
    #relavance threshold
    thR=0.15
    
    #data standard deviation
    sigma = 1 
    
    #historical dat - fitting prior
    n0 = 50
    prior.m = 0.25
    s0 = sqrt(sigma^2/n0)  #historical prior sd
    
    #priors
    info = c(1,prior.m,s0)
    base = c(1,0,10)
    mix = cbind(c(0.5, prior.m, s0), c(0.5,prior.m,10))
    power = c(1, prior.m, s0)
    prior.listI =list("vague"=base, "info"=info, "mix"=mix, "power"=power)
  
  
  id=which(abs(true.m.vec-prior.m)==min(abs(true.m.vec-prior.m)))
  desOutF_fixed= designF(id,outcome=outcome,dec_cost=dec_cost,max_n=max_n,designGrid=designGrid,th0=th0,thR=thR,prior.listI=prior.listI,n0=n0,sigma=sigma)
  
  save(desOutF_fixed, file= paste("desOutF_",outcome,"_",dec_cost,"LSS.RData",sep=""))
}
