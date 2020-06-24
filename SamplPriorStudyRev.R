library(truncnorm)
library(ggplot2)
library(gridExtra)

rm(list=ls())
source('Design_Fn.R')

SamplPrirorStudy = function ()
{
  #hypothesis threshold
  th0=0
  
  #data variance
  sigma = 1   
  
  #sampling priors considered
  true.m.v=seq(-0.7,0.7,by=0.35)
  true.sd=sqrt(0.02)
  
  #costs
  c.alpha1 = 0.95
  c.beta1 = 0.05
  
  #parameters for integration grids
  sl=5
  frl=400
  
  #sample size
  n.targ=100
  
  
  # compute conditional test error rates ------------------------------------
  
  #evaluation grid ('true' thetas)
  yErr=seq(-3,3,by=true.sd/100)
  
  #conditional test error rates at each evaluation point for each sampling prior
  Err.sampling=sapply(true.m.v, function(m) 
  {
    bf.v=sapply(yErr, function(z) 
    {
      #data grid
      y_grid=seq(z-sl*sigma/sqrt(n.targ),z+sl*sigma/sqrt(n.targ),by=sigma/sqrt(n.targ)/frl)
      norm.dat=sum(dnorm(y_grid,z,sigma/sqrt(n.targ)))
      w.dat=dnorm(y_grid,z,sigma/sqrt(n.targ))/norm.dat
      
      #bayes factor
      bf=sapply(y_grid, function(t) BF(dat=t,n=n.targ,dataPar=sigma,prior='sampling',priorParsWeights=c(1,m,true.sd),th0=th0))
      
      #test error rate
      Err=w.dat%*%(bf<c.beta1/c.alpha1)
      out=ifelse(z<=th0,Err,1-Err)
      return(out)}
    ) })
  
  
  # compute truncated sampling prior densities ------------------------------
  
  yErr0=yErr[yErr<=th0]
  yErr1=yErr[yErr>th0]
  densities=as.vector(sapply(true.m.v,function(z) c(priorTr(th=yErr0,prior='sampling',priorPars=c(z,true.sd),hp='null',th0=th0),
                                                    priorTr(th=yErr1,prior='sampling',priorPars=c(z,true.sd),hp='alt',th0=th0))))
  
  
  
  
  # Figure S2 ---------------------------------------------------------------
  
  
  plotD=data.frame("Theta"=rep(yErr[yErr>-0.9 & yErr<=0.9],length(true.m.v)),
                   "Conditional test errors"=as.vector(Err.sampling[yErr>-0.9 & yErr<=0.9,]),
                   "density"=densities[yErr>-0.9 & yErr<=0.9],
                   "Sampling prior"=as.factor(rep(true.m.v,each=length(yErr[yErr>-0.9 & yErr<=0.9]))))
  
  theme=ggplot()+ theme_light() + 
    scale_color_viridis_d(begin = 0, end = 0.9)+
    xlab(expression(theta))+ geom_vline(xintercept=0,linetype="dashed",colour="gray")
  
  
  g1 <- theme+ geom_line(aes(Theta,Conditional.test.errors,  color=Sampling.prior), 
                         size=0.4, data=plotD) + ylab("Conditional test error rates")+
    ggtitle("Conditional test error rates (BF-based decision)")+ labs(color = "Sampling prior")
  
  g2 <- theme+ geom_line(aes(Theta,density,  color=Sampling.prior), 
                         size=0.4, data=plotD) + ylab("Density")+
    ggtitle("Sampling prior densities (truncated)")+ labs(color = "Sampling prior")
  
  
  
  
  grid.arrange(g1,g2, nrow = 2, ncol=1)
  
  out_dir="~/Figures"
  pdf(file = paste(out_dir,"/","sampling_prior_study.pdf",sep=""), width=7,height=9, pointsize=12)
  grid.arrange(g1,g2, nrow = 2, ncol=1)
  graphics.off()
}