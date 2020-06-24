library(ggplot2)
library(latex2exp)
library(ggpubr)
library(gridExtra)
library(grid)
library(Cairo)

rm(list=ls())


createFigure = function(outcome=c('normal','binomial'),dec_cost=c('normalized','fixed'),opt_sample=c('Nfix', 'Nopt', 'GoalSampling'),type=c("sim","data"))
{
  dec=match.arg(dec_cost)
  outL=match.arg(outcome)
  sampOpt=match.arg(opt_sample)
  type=match.arg(type)
  
  n.targ=100
  c.alpha1 = 0.95
  c.beta1 = 0.05
  
  if (outL=='normal')
  {
    prior.m=0.25
    true.m.vec=seq(-0.9,1.25,by=0.05)
    true.v.vec=0.02
    id=which(abs(true.m.vec-prior.m)==min(abs(true.m.vec-prior.m)))
    n.list=1:300
    sigma = 1   
    n0 = 50
    s0 = sqrt(sigma ^ 2 / n0)
    th0=0
    designGrid=cbind(rep(true.m.vec,each=length(true.v.vec)),rep(true.v.vec,length(true.m.vec)))
    
    load(paste("desOutF_",outL,"_",dec,".RData",sep=""))
    
  }
  
  if (outL=='binomial' & type=="sim")
  {
    prior.m=0.5
    true.a.vec=seq(1,49,by=2)
    true.b.vec=50-true.a.vec
    n.list=1:300
    th0=0.2
    true.m.vec=true.a.vec/(true.a.vec+true.b.vec)
    id=which(true.m.vec==prior.m)
    designGrid=cbind(true.a.vec,true.b.vec)
    
    load(paste("desOutF_",outL,"_",dec,".RData",sep=""))
  }
  
  if (outL=='binomial' & type=="data")
  {
    prior.m=0.275
    true.a.vec=seq(1,39,by=2)
    true.b.vec=40-true.a.vec
    n.list=1:300
    th0=0.075
    true.m.vec=true.a.vec/(true.a.vec+true.b.vec)
    id=which(true.m.vec==prior.m)
    designGrid=cbind(true.a.vec,true.b.vec)
    
    load("desOutF_Example.RData")
  }
  
  #index=dim(designGrid)[1]
  
  design=desOutF_fixed[[id]]
  design[[6]]=NULL
  
  # cost computation vague prior example --------------------------------------------------------
  
  
  if (dec=='normalized' & outL=='normal' & sampOpt=="Nopt")
  {
    load(paste("desOutF_",outL,"_",dec,"_LSS.RData",sep="")) #design computed up to 800 observations for sampling prior of interest
    design=desOutF_fixed
    design[[6]]=NULL
    
    n.listC=1:length(design$vague['type1',])
    x=n.listC
    n_te=2:length(n.listC)
    n_mse=2:length(n.listC)
    
    y=design$vague['ErrSum',]
    lo1 = loess(y~x,span=0.05)
    c_te= -(predict(lo1)[n_te]-predict(lo1)[n_te-1])
    
    y=design$vague['MSE',]
    lo2 = loess(y~x,span=0.005)
    c_mse=-(predict(lo2)[n_mse]-predict(lo2)[n_mse-1])
    
    plotM_oc=data.frame("n"=rep(n.listC,5),
                        "Prior"=rep(c("Vague","Informative","Robust mixture",
                                      "EB power","Sampling"),each=length(n.listC)),
                        "type1"=do.call('c',lapply(design, function(x) x['type1',])),
                        'type2'=do.call('c',lapply(design, function(x) x['type2',])),
                        'ErrSum'=do.call('c',lapply(design, function(x) x['ErrSum',])),
                        'MSE'=do.call('c',lapply(design, function(x) x['MSE',])),
                        'PrIndet'=do.call('c',lapply(design, function(x) x['PrIndet',])),
                        'PowLoss'=do.call('c',lapply(design, function(x) x['PowLoss',])),
                        'NGain'=do.call('c',lapply(design, function(x) x['NGain',])),
                        'gamma'=do.call('c',lapply(design, function(x) x['gamma',])))
    
    plotM_costs=data.frame("n"=c(n_te,n_mse),
                           "OC"=rep(c("Sum of avg. test error rates","Avg. MSE"),each=length(n_te)),
                           "costs"=c(c_te,c_mse))
    
    
    theme1=ggplot()+ xlab("n")+ ylab("")
    
    g1= theme1+ geom_line(aes(n, costs),  color='black', size=0.6, data = subset(plotM_costs, OC=='Sum of avg. test error rates')) + 
      theme(legend.position = "none") +
      ggtitle(TeX('c_n^{SATE}')) 
    
    g2 = theme1+ geom_line(aes(n, costs), color='black', size=0.6, data = subset(plotM_costs, OC=='Avg. MSE')) + 
      theme(legend.position = "none")  +
      ggtitle(TeX('c_n^{AMSE}')) 
    
    g3 = theme1+ geom_line(aes(n, ErrSum),  color='black', size=0.6, data = subset(plotM_oc, Prior=='Vague')) + 
      theme(legend.position = "none")  +
      ggtitle('SATE') +
      xlab("")
    
    g4 = theme1+ geom_line(aes(n, MSE),  color='black', size=0.6, data = subset(plotM_oc, Prior=='Vague')) + 
      theme(legend.position = "none")  +
      ggtitle('AMSE')+
      xlab("")
    
    out_dir="~/Figures"
    cairo_pdf(file = paste(out_dir,"/",outL,"_",dec,"_costExample.pdf", sep = ""), width=6,height=6, pointsize=11)
    grid.arrange(g3,g4,g1,g2, nrow = 2, ncol=2)
    graphics.off()
  }
  
  # cost computation simulation and data --------------------------------------------------------
  
  
  if (dec=='normalized' & sampOpt=='Nopt')
  {
    n.listC=1:length(design$vague['type1',])
    
    if (type=='data')
    {
      n_teV=lapply(design, function(x) max(which(x['type1',]==x['type1',x['type1',]>0.15][length(x['type1',x['type1',]>0.15])])+1,
                                           which(x['type2',]==x['type2',x['type2',]>0.20][length(x['type2',x['type2',]>0.20])])+1))
      n_mseV=lapply(design, function(x) max(which(x['PrIndet',]==x['PrIndet',x['PrIndet',]>0.1][length(x['PrIndet',x['PrIndet',]>0.1])])+1,
                                            which(x['PowLoss',]==x['PowLoss',x['PowLoss',]>0.30][length(x['PowLoss',x['PowLoss',]>0.30])])+1))
    }
    if (type=='sim')
    {
      n_teV=lapply(design, function(x) max(which(x['type1',]==x['type1',x['type1',]>0.10][length(x['type1',x['type1',]>0.10])])+1,
                                           which(x['type2',]==x['type2',x['type2',]>0.20][length(x['type2',x['type2',]>0.20])])+1))
      n_mseV=lapply(design, function(x) max(which(x['PrIndet',]==x['PrIndet',x['PrIndet',]>0.1][length(x['PrIndet',x['PrIndet',]>0.1])])+1,
                                            which(x['PowLoss',]==x['PowLoss',x['PowLoss',]>0.30][length(x['PowLoss',x['PowLoss',]>0.30])])+1))
    }
    
    n_te=ifelse(n_teV$vague==(length(n.listC)+1),length(n.listC),n_teV$vague)
    n_mse=ifelse(n_mseV$vague==(length(n.listC)+1),length(n.listC),n_mseV$vague)
    
    
    w_te=n_te/(n_te+n_mse)
    x=n.listC
    
    y=design$vague['ErrSum',]
    if (outL=='normal')
    {
      lo1 = loess(y~x,span=0.1)
    }
    if (outL=='binomial')
    {
      lo1 = loess(y~x,span=0.7)
    }
    c_te= -(predict(lo1)[n_te]-predict(lo1)[n_te-1])
    #  plot(x,y)
    #  lines(predict(lo1),col=2)
    
    y=design$vague['MSE',]
    lo2 = loess(y~x,span=0.05)
    c_mse=-(predict(lo2)[n_mse]-predict(lo2)[n_mse-1])
    #  plot(x,y)
    #  lines(predict(lo2),col=2)
    
    costs.v=c('c_q'=((1-w_te)*c_te)/(w_te*c_mse),'c_n'=c_te/w_te,'c_te'=c_te,"c_mse"=c_mse,'w_te'=w_te,'n_te'=n_te,'n_mse'=n_mse)
    oc.v=c(design$vague['ErrSum',n_te],design$vague['MSE',n_mse],design$vague['NGain',n_mse])
    out=list('Costs'=costs.v,'OC'=oc.v)
    save(out,file=paste(outL,'_',type,'_costs.RData',sep=''))
  }
  else
  {
    load(paste(outL,'_',type,'_costs.RData',sep=''))
    w_te=out$Costs['w_te']
    c_mse=out$Costs['c_mse']
    c_te=out$Costs['c_te']
    n_te=out$Costs['n_te']
    n_mse=out$Costs['n_mse']
  }
  
  # oc plots for fixed sampling prior ---------------------------------------
  
  if (dec=='normalized' & sampOpt=='Nopt')
  {
    
    plotM_oc=data.frame("n"=rep(n.listC,5),
                        "Prior"=rep(c("Vague","Informative","Robust mixture",
                                      "EB power","Sampling"),each=length(n.listC)),
                        "type1"=do.call('c',lapply(design, function(x) x['type1',])),
                        'type2'=do.call('c',lapply(design, function(x) x['type2',])),
                        'ErrSum'=do.call('c',lapply(design, function(x) x['ErrSum',])),
                        'MSE'=do.call('c',lapply(design, function(x) x['MSE',])),
                        'PrIndet'=do.call('c',lapply(design, function(x) x['PrIndet',])),
                        'PowLoss'=do.call('c',lapply(design, function(x) x['PowLoss',])),
                        'NGain'=do.call('c',lapply(design, function(x) x['NGain',])),
                        'gamma'=do.call('c',lapply(design, function(x) x['gamma',])))
    
    
    # colors=c(
    #   "Vague"= "black",
    #   "Informative"=  "red",
    #   "Robust mixture"="green",
    #   "EB power"="blue",
    #   "Sampling"="magenta"
    # )
    
    ltys=c(
      "Vague"= "solid",
      "Informative"=  "dotted",
      "Robust mixture"="dashed",
      "EB power"="longdash",
      "Sampling"="dotdash"
    )
    
    plotM_oc$Prior=factor(plotM_oc$Prior,levels=c("EB power","Informative","Robust mixture","Vague","Sampling"))
    
    theme=ggplot()+ theme_light() + xlab("n")   +
      ylab("")+scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM_oc$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                        expression(paste(pi[a],":Informative")),
                                                                                        expression(paste(pi[a],":Robust mixture")),
                                                                                        expression(paste(pi[a],":Vague")),
                                                                                        expression(paste(pi[s],":Sampling"))))+
      scale_linetype_manual(values=ltys,labels = c(expression(paste(pi[a],":EB power")),
                                                                                 expression(paste(pi[a],":Informative")),
                                                                                 expression(paste(pi[a],":Robust mixture")),
                                                                                 expression(paste(pi[a],":Vague")),
                                                                                 expression(paste(pi[s],":Sampling"))))
    
    
    if (type=='data')
    {
      theme=ggplot()+ theme_light() + xlab("n")   +
        ylab("")+scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM_oc$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                                     expression(paste(pi[a],":Informative")),
                                                                                                     expression(paste(pi[a],":Robust mixture")),
                                                                                                     expression(paste(pi[a],":RSN")),
                                                                                                     expression(paste(pi[s],":Sampling"))))+
        scale_linetype_manual(values=ltys,labels = c(expression(paste(pi[a],":EB power")),
                                                                                   expression(paste(pi[a],":Informative")),
                                                                                   expression(paste(pi[a],":Robust mixture")),
                                                                                   expression(paste(pi[a],":RSN")),
                                                                                   expression(paste(pi[s],":Sampling"))))
    }
    
    g1 = theme + geom_line(aes(n, type1,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) + theme(legend.title = element_blank())+
      theme(legend.key.width = unit(3, 'lines'),legend.text = element_text(size = 10)) +  guides(color = guide_legend(override.aes = list(size = 0.4),ncol=3))+
      guides(shape = guide_legend(override.aes = list(size = 0.4))) + geom_vline(xintercept=n_te,linetype="dashed",colour="gray") +
      theme(legend.position='top') + ylim(0,NA) +
      ggtitle("Average type I error") +
      ggtitle("Avg. type I error rate") +
      xlab("")
    
    g2 = theme + geom_line(aes(n, type2,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) + 
      theme(legend.position = "none") + geom_vline(xintercept=n_te,linetype="dashed",colour="gray") +
      ggtitle("Avg. type II error rate") + scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1),minor_breaks=c(.1,.3,.5,.7,.9)) +
      xlab("")
    
    g3 = theme + geom_line(aes(n, ErrSum,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) + 
      theme(legend.position = "none") + geom_vline(xintercept=n_te,linetype="dashed",colour="gray")+ ylim(0,NA) +
      ggtitle('SATE') +
      xlab("")
    
    g4 = theme + geom_line(aes(n, MSE,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) + 
      theme(legend.position = "none") + geom_vline(xintercept=n_mse,linetype="dashed",colour="gray")  + ylim(0,NA) +
      ggtitle('AMSE')
    
    
    g5 = theme + geom_line(aes(n, PrIndet,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) + 
      theme(legend.position = "none")+ scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1),minor_breaks=c(.1,.3,.5,.7,.9)) +
      ggtitle(TeX('Avg. P(Indeterm.)')) + geom_vline(xintercept=n_mse,linetype="dashed",colour="gray")
    
    
    g6 = theme + geom_line(aes(n, PowLoss,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) + geom_vline(xintercept=n_mse,linetype="dashed",colour="gray")  + 
      theme(legend.position = "none")+ scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1),minor_breaks=c(.1,.3,.5,.7,.9)) +
      ggtitle(TeX('Avg. power loss (80%)'))
    
    g7 <- theme + geom_line(aes(n, NGain,  color=Prior, linetype=Prior), size=0.6, data = plotM_oc) +  geom_vline(xintercept=n_mse,linetype="dashed",colour="gray")  +
      theme(legend.position = "none")+  ggtitle(TeX('Avg. SS gain (80%)')) +
      ylab("")
    
    
    a= ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_text(size=1)) ,
                         g2,g3,
                         ncol = 3,nrow = 1,common.legend = T,align="v")
    b= grid.arrange(g5,g6,g7,g4, nrow = 1, ncol=4)
    
    out_dir="~/Figures"
    
    if (type=='data')
    {
      cairo_pdf(file = paste(out_dir,"/","data_example_cost.pdf", sep = ""), width=10,height=7.5, pointsize=11)
    }
    else
    {
      cairo_pdf(file = paste(out_dir,"/",outL,"_",dec,"_cost.pdf", sep = ""), width=10,height=7.5, pointsize=11)
    }
    grid.arrange(a,b, nrow = 2, ncol=1)
    graphics.off()
  }
  
  
  # sensitivity analyses plots ----------------------------------------------
  
  if (type=="sim")
  {
    load(paste("desOutF_",outL,"_",dec,".RData",sep=""))
  }
  
  optOC=lapply(desOutF_fixed,function(y) 
  {
    yS=y[[6]]
    y[[6]]=NULL
    
    if (outL=='normal')
      pPs=pnorm(th0,yS[1],sqrt(yS[2]))
    if (outL=='binomial')
      pPs=pbeta(th0,yS[1],yS[2])
    
    lapply(1:5,function(z){
      
      x=y[[z]]
      irisk=w_te*(1/c_te)*x['ErrSum',] + (1-w_te)*(1/c_mse)*x['MSE',] + n.list
      
      if (sampOpt=="Nopt")     
      {
        minirisk=which(irisk==min(irisk))[1]
      }
      if (sampOpt=="Nfix")     
      {
        minirisk=n.targ
      }
      if (sampOpt=="GoalSampling" & type=='sim')     
      {
        n_teV= max(which(x['type1',]==x['type1',x['type1',]>0.1][length(x['type1',x['type1',]>0.1])]),
                   which(x['type2',]==x['type2',x['type2',]>0.20][length(x['type2',x['type2',]>0.20])]))
        n_mseV=max(which(x['PrIndet',]==x['PrIndet',x['PrIndet',]>0.1][length(x['PrIndet',x['PrIndet',]>0.1])]),
                   which(x['PowLoss',]==x['PowLoss',x['PowLoss',]>0.30][length(x['PowLoss',x['PowLoss',]>0.3])]))
        
        n_teV=which(x['ErrSum',]==x['ErrSum',x['ErrSum',]>desOutF_fixed[[id]]$vague['ErrSum',n_te]][length(x['ErrSum',x['ErrSum',]>desOutF_fixed[[id]]$vague['ErrSum',n_te]])])
        n_mseV=which(x['MSE',]==x['MSE',x['MSE',]>desOutF_fixed[[id]]$vague['MSE',n_mse]][length(x['MSE',x['MSE',]>desOutF_fixed[[id]]$vague['MSE',n_mse]])])
        minirisk=ifelse(length(n_teV)==0 & length(n_mseV)==0,1,
                        ifelse(length(n_teV)!=0 & length(n_mseV)==0   & n_teV==length(n.list),length(n.list),
                               ifelse(length(n_teV)!=0 & length(n_mseV)==0   & n_teV!=length(n.list),max(n_teV+1,n_mseV+1),
                                      ifelse(length(n_mseV)!=0  & length(n_teV)==0  & n_mseV==length(n.list),length(n.list),
                                             ifelse(length(n_mseV)!=0  & length(n_teV)==0  & n_mseV!=length(n.list),max(n_teV+1,n_mseV+1),
                                                    ifelse(length(n_teV)!=0 & length(n_mseV)!=0  & (n_teV==length(n.list) | n_mseV==length(n.list)),length(n.list),max(n_teV+1,n_mseV+1)))))))
      }
      
      
      if(dec=="normalized")
      {
        cost=c(c.alpha1/pPs,c.beta1/(1-pPs),((1-w_te)*c_te)/(w_te*c_mse),c_te/w_te)
        names(cost)=NULL
      }
      if(dec=="fixed")
      {
        cost=c(c.alpha1,c.beta1,((1-w_te)*c_te)/(w_te*c_mse),c_te/w_te)
        names(cost)=NULL
      }
      
      outOP=c(x[c('type1','type2','MSE','gamma'),minirisk],'n.opt'=minirisk,'irisk'=irisk[minirisk],'costs'=cost)
    })})
  
  
  
  
  plotM_oc=data.frame("m"=rep(true.m.vec,each=5),
                      "Prior"=rep(c("Vague","Informative","Robust mixture",
                                    "EB power","Sampling"),length(true.m.vec)),
                      "type1"=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){y['type1']}))})),
                      'type2'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){y['type2']}))})),
                      'MSE'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){sqrt(y['MSE'])}))})),
                      'gamma'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){log(y['gamma'])}))})),
                      'opt.n'= do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){y['n.opt']}))})),
                      'irisk'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){y['irisk']}))})),
                      'cost1'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){log(y['costs1'])}))})),
                      'cost2'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){log(y['costs2'])}))})),
                      'costMSE'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){log(y['costs3'])}))})),
                      'costN'=do.call('c',lapply(optOC, function(x) {do.call('c', lapply(x,function(y){log(y['costs4'])}))})))
  
  plotM_costs=data.frame("m"=rep(true.m.vec,4),
                         "OC"=rep(c("type1","type2","MSE","N"),each=length(true.m.vec)),
                         "costs"=c(subset(plotM_oc, Prior=="Sampling")$cost1,subset(plotM_oc, Prior=="Sampling")$cost2,
                                   subset(plotM_oc, Prior=="Sampling")$costMSE,subset(plotM_oc, Prior=="Sampling")$costN))
  
  
  
  ltys=c(
    "Vague"= "solid",
    "Informative"=  "dotted",
    "Robust mixture"="dashed",
    "EB power"="longdash",
    "Sampling"="dotdash"
  )
  
  ltys2=c(
    "type1"= "solid",
    "type2"=  "dotted",
    "MSE"="dashed",
    "N"= "dotdash"
  )
  
  plotM_oc$Prior=factor(plotM_oc$Prior,levels=c("EB power","Informative","Robust mixture","Vague","Sampling"))
  
  theme=ggplot()+ theme_light() + 
    xlab("Sampling prior mean") + geom_vline(xintercept=prior.m,linetype="dashed",colour="gray")  +
    ylab("")+scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM_oc$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                                 expression(paste(pi[a],":Informative")),
                                                                                                 expression(paste(pi[a],":Robust mixture")),
                                                                                                 expression(paste(pi[a],":Vague")),
                                                                                                            expression(paste(pi[s],":Sampling"))))+
    scale_linetype_manual(values=ltys,labels = c(expression(paste(pi[a],":EB power")),
                                                                               expression(paste(pi[a],":Informative")),
                                                                               expression(paste(pi[a],":Robust mixture")),
                                                                               expression(paste(pi[a],":Vague")),
                                                                               expression(paste(pi[s],":Sampling"))))
  
  
  
  theme2=ggplot() + theme_light()+
    xlab("Sampling prior mean") + geom_vline(xintercept=prior.m,linetype="dashed",colour="gray") +
    ylab("") + scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM_costs$OC), labels = c("MSE","Sample size","Type I error","Type II error"))+
    scale_linetype_manual(values=ltys2,labels = c("MSE","Sample size","Type I error","Type II error"))

  
  if (type=='data')
  {
    theme=ggplot()+ theme_light() + 
      xlab("Sampling prior mean") + geom_vline(xintercept=prior.m,linetype="dashed",colour="gray")  +
      ylab("")+scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM_oc$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                                   expression(paste(pi[a],":Informative")),
                                                                                                   expression(paste(pi[a],":Robust mixture")),
                                                                                                   expression(paste(pi[a],":RSN")),
                                                                                                              expression(paste(pi[s],":Sampling"))))+
      scale_linetype_manual(values=ltys,labels = c(expression(paste(pi[a],":EB power")),
                                                                                 expression(paste(pi[a],":Informative")),
                                                                                 expression(paste(pi[a],":Robust mixture")),
                                                                                 expression(paste(pi[a],":RSN")),
                                                                                 expression(paste(pi[s],":Sampling"))))
  }
  
  g1 = theme+ geom_line(aes(m, type1,  color=Prior, linetype=Prior), size=0.5, data = plotM_oc) + 
    theme(legend.title = element_blank())+
    theme(legend.key.width = unit(1.5, 'lines'),legend.text = element_text(size = 10))+ guides(color = guide_legend(override.aes = list(size = 0.4),ncol=3))+
    guides(shape = guide_legend(override.aes = list(size = 0.4)))+
    theme(legend.position=c(0.2,0.8))+
    ggtitle("Avg. type I error rate") +
    xlab("") 
  
  g2 = theme+ geom_line(aes(m, type2,  color=Prior, linetype=Prior), size=0.5, data = plotM_oc) + 
    theme(legend.position = "none")+
    ggtitle("Avg. type II error rate")
  
  
  if (sampOpt=="Nopt")     
  {
    g2=g2 +xlab("")
  }
  
  g3 = theme+ geom_line(aes(m, MSE,  color=Prior, linetype=Prior), size=0.5, data = plotM_oc) + 
    theme(legend.position = "none")+
    ggtitle('Avg. root MSE')  +
    xlab("")
  
  
  g4 = theme+ geom_line(aes(m, opt.n,  color=Prior, linetype=Prior), size=0.5, data = plotM_oc) + 
    theme(legend.position = "none")+
    ggtitle('Optimal sample size') 
  
  if (sampOpt!="GoalSampling")     
  {
    g4=g4 +xlab("")
  }
  
  g5 = theme+ geom_line(aes(m, irisk,  color=Prior, linetype=Prior), size=0.5, data = plotM_oc) + 
    theme(legend.position = "none")+
    ggtitle(TeX('Integrated risk'))
  
  g6=theme2+ 
    geom_line(aes(m, costs, color=OC, linetype=OC), size=0.6, data = plotM_costs) +
    theme(legend.title = element_blank()) +
    theme(legend.position = "top",legend.text = element_text(size = 9)) +
    guides(color = guide_legend(override.aes = list(size = 0.5),ncol=2)) +
    ggtitle(TeX('Log costs')) 
  
  
  
  out_dir="~/Figures"
  
  if (type=='data')
  {
    a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_text(size=1)) ,
                        g3,g2,g4,ncol = 2,nrow = 2,common.legend = T,align="v")
    b=grid.arrange(g6,g5, nrow = 1, ncol=2)
    
    cairo_pdf(file = paste(out_dir,"/","data_example_sens.pdf", sep = ""), width=6,height=8, pointsize=11)
    grid.arrange(a,b, nrow = 2, ncol=1,heights=c(2,1))
    
    graphics.off()
  }
  else
  {
    if (dec=="normalized")     
    {
      if (sampOpt=="Nopt")     
      {
        a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T)) +
                              theme(legend.title=element_text(size=1)),
                            g3,g2,g4,ncol = 2,nrow = 2,common.legend = T,align="v")
        b=grid.arrange(g6,g5, nrow = 1, ncol=2)
        
        cairo_pdf(file = paste(out_dir,"/","/","design_",dec,"_",sampOpt,"_",outL,"_sens.pdf", sep = ""), width=6, height=8, pointsize=11)
        
        grid.arrange(a,b, nrow = 2, ncol=1,heights=c(2,1))
        
        graphics.off()
      }
      
      if (sampOpt=="Nfix")     
      {
        a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_text(size=1)) ,
                            g3,ncol = 2,nrow = 1,common.legend = T,align="v")
        c=grid.arrange(g2,g6, nrow = 1, ncol=2)
        b=grid.arrange(g5, nrow = 1, ncol=1)
        
        cairo_pdf(file = paste(out_dir,"/","/","design_",dec,"_",sampOpt,"_",outL,"_sens.pdf", sep = ""), width=6, height=8, pointsize=11)
        grid.arrange(a,c,b, nrow = 3, ncol=1,heights=c(1,1,1))
        graphics.off()
        
      }
      if (sampOpt=="GoalSampling")     
      {
        a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_text(size=1)) ,
                            g3,g2,g4,ncol = 2,nrow = 2,common.legend = T,align="v")
        
        cairo_pdf(file = paste(out_dir,"/","/","design_",dec,"_",sampOpt,"_",outL,"_sens.pdf", sep = ""), width=6,height=6, pointsize=11)
        grid.arrange(a,nrow = 1, ncol=1)
        graphics.off()
        
      }
    }
    
    if (dec=="fixed")     
    {
      if (sampOpt=="Nopt")     
      {
        a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_text(size=1)) ,
                            g3,g2,g4,ncol = 2,nrow = 2,common.legend = T,align="v")
        b=grid.arrange(g5, nrow = 1, ncol=1)
        
        cairo_pdf(file = paste(out_dir,"/","/","design_",dec,"_",sampOpt,"_",outL,"_sens.pdf", sep = ""), width=6,height=8, pointsize=11)
        
        grid.arrange(a,b, nrow = 2, ncol=1,heights=c(2,1))
        
        graphics.off()
      }
      
      if (sampOpt=="Nfix")     
      {
        a=ggpubr::ggarrange(g1+ guides(col = guide_legend(nrow = 1, byrow = T))+theme(legend.title=element_text(size=1)) ,
                            g3,g2,g5,ncol = 2,nrow = 2,common.legend = T,align="v")
        
        cairo_pdf(file = paste(out_dir,"/","/","design_",dec,"_",sampOpt,"_",outL,"_sens.pdf", sep = ""), width=6,height=6, pointsize=11)
        grid.arrange(a,nrow = 1, ncol=1)
        graphics.off()
        
      }
    }
  }
  
  # bayes factor plot ------------------------------------------------------------
  
  if (outL=="normal" & dec=="normalized" & sampOpt=="Nfix" & type=='sim')     
  {
    
    info = c(1,prior.m,s0)
    base = c(1,0,10)
    mix = cbind(c(0.5, prior.m, s0), c(0.5, prior.m,10))
    power = c(1, prior.m, s0)
    
    #source('Design_Fn.R')
    prior.list = list("vague"=base, "info"=info, "mix"=mix, "power"=power)
    y_grid=seq(-1,1,by=0.05)
    BF.list=lapply(1:length(prior.list), function(z) 
      sapply(y_grid, function(t) BF(dat=t,n=n.targ,dataPar=sigma,prior=names(prior.list)[z],priorParsWeights=prior.list[[z]],th0=th0)))
    names(BF.list)=names(prior.list)
    
    plotM=data.frame("y"=rep(y_grid,4),
                     "lBF"=c(log(unlist(BF.list))),
                     "Prior"=rep(c("Vague","Informative","Robust mixture","EB power"),each=length(y_grid)))
    plotM$Prior=factor(plotM$Prior,levels=c("EB power", "Informative","Robust mixture","Vague"))
    
    
    theme2=ggplot()+ theme_light() + 
      geom_vline(xintercept=prior.m,linetype="dashed",colour="gray")  +
       scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                                   expression(paste(pi[a],":Informative")),
                                                                                                   expression(paste(pi[a],":Robust mixture")),
                                                                                                   expression(paste(pi[a],":Vague"))))+
      scale_linetype_manual(values=ltys,labels = c(expression(paste(pi[a],":EB power")),
                                                                                 expression(paste(pi[a],":Informative")),
                                                                                 expression(paste(pi[a],":Robust mixture")),
                                                                                 expression(paste(pi[a],":Vague"))))
    
    
    g = ggplotGrob(theme2+
                     theme(plot.background = element_rect(colour = "black")) + 
                     geom_line(aes(y, lBF,  color=Prior, linetype=Prior), size=0.3, data =subset(plotM, y>0 & y<0.4)) +
                     theme(plot.background = element_rect(colour = "black")) + theme(legend.position = "none") + 
                     geom_hline(yintercept=log(c.beta1/c.alpha1), size=0.3,linetype="dashed",colour="deepskyblue4",alpha=0.7)+
                     xlab("")+ylab("")+theme(axis.text=element_text(size=5),
                                             axis.title=element_text(size=5))+
                     theme(plot.margin=unit(c(0.1,0.1,0,0),"cm"),
                           panel.spacing = unit(0,"null")))
    
    pbf = theme2 + geom_line(aes(y, lBF, color=Prior, linetype=Prior), size=0.4, data = plotM) +
      annotation_custom(
        grob = g,
        xmin = -0.3,
        xmax = 1,
        ymin = 10,
        ymax = 58
      ) +
      geom_hline(yintercept=log(c.beta1/c.alpha1),linetype="dashed",colour="deepskyblue4", size=0.3,alpha=0.7)+
      xlab(bquote("Data " * bar(y))) +ggtitle("Bayes factor") +
      theme(legend.position="bottom",legend.title = element_blank(),
            legend.text.align = 0,
            legend.text = element_text(size = 7))+ guides(color = guide_legend(override.aes = list(size = 0.5),ncol=2))+
      guides(shape = guide_legend(override.aes = list(size = 0.5)))+
      ylab("log(BF)")
    
    
    pgamma = theme + geom_line(aes(m,gamma,color=Prior, linetype=Prior), size=0.4, data = plotM_oc) +
      ggtitle("Test decision threshold") +
      theme(legend.position="bottom",legend.title = element_blank(),
            legend.text.align = 0,
            legend.text = element_text(size = 7))+ 
      guides(color = guide_legend(override.aes = list(size = 0.5),ncol=2), shape = guide_legend(override.aes = list(size = 0.5)))+
      ylab(expression(log(gamma))) + 
      xlab("Sampling prior mean")
    
    
    cairo_pdf(file = paste(out_dir,"/","design_",dec,"_",sampOpt,"_",outL,"_BF_gamma.pdf", sep = ""), width=7, height=4, pointsize=12)
    grid.arrange(pgamma,pbf,  nrow = 1)
    graphics.off()
    
    
    #sensitivity to historical information location and mixture prior vague component, mixture prior weight
    
    
    info = c(1,prior.m,s0)
    base = c(1,0,10)
    mix = cbind(c(0.5, prior.m, s0), c(0.5, prior.m,10))
    mix_80 = cbind(c(0.8, prior.m, s0), c(0.2, prior.m,10))
    mix_unit = cbind(c(0.5, prior.m, s0), c(0.5, prior.m, 1))
    power = c(1, prior.m, s0)
    
    prior.list = list("vague"=base, "info"=info, "mix"=mix,"mix_80"=mix_80,"mix_unit"=mix_unit, "power"=power)
    y_grid=seq(-1,1,by=0.05)
    BF.list=lapply(1:length(prior.list), function(z) 
      sapply(y_grid, function(t) BF(dat=t,n=n.targ,dataPar=sigma,prior=names(prior.list)[z],priorParsWeights=prior.list[[z]],th0=th0)))
    names(BF.list)=names(prior.list)
    
    plotM=data.frame("y"=rep(y_grid,6),
                     "lBF"=c(log(unlist(BF.list))),
                     "Prior"=rep(c("Vague","Informative","Robust mixture","Robust mixture (0.8)","Robust mixture (unit inform.)","EB power"),each=length(y_grid)))
    plotM$Prior=factor(plotM$Prior,levels=c("EB power", "Informative","Robust mixture","Robust mixture (0.8)","Robust mixture (unit inform.)","Vague"))
    
    ltys2=c(
      "Vague"= "solid",
      "Informative"=  "dotted",
      "Robust mixture"="dashed",
      "Robust mixture (unit inform.)"="twodash",
      "Robust mixture (0.8)"="F1",
      "EB power"="longdash",
      "Sampling"="dotdash"
    )
  
    theme2=ggplot()+ theme_light() + geom_vline(xintercept=prior.m,linetype="dashed",colour="gray") +
      scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                                   expression(paste(pi[a],":Informative")),
                                                                                                   expression(paste(pi[a],":Robust mixture")),
                                                                                                   expression(paste(pi[a],":Robust mixture (unit inform.)")),
                                                                                                   expression(paste(pi[a],":Robust mixture (0.8)")),
                                                                                                   expression(paste(pi[a],":Vague")),
                                                                                                   expression(paste(pi[s],":Sampling"))))+
      scale_linetype_manual(values=ltys2,labels = c(expression(paste(pi[a],":EB power")),
                                                                                 expression(paste(pi[a],":Informative")),
                                                                                 expression(paste(pi[a],":Robust mixture")),
                                                                                 expression(paste(pi[a],":Robust mixture (unit inform.)")),
                                                                                 expression(paste(pi[a],":Robust mixture (0.8)")),
                                                                                 expression(paste(pi[a],":Vague")),
                                                                                 expression(paste(pi[s],":Sampling"))))
    
    
    g = ggplotGrob(theme2+
                     theme(plot.background = element_rect(colour = "black")) + 
                     geom_line(aes(y, lBF,  color=Prior, linetype=Prior), size=0.35, data =subset(plotM, y>0 & y<0.4)) +
                     theme(plot.background = element_rect(colour = "black")) + theme(legend.position = "none") + 
                     geom_hline(yintercept=log(c.beta1/c.alpha1), size=0.35,linetype="dashed",colour="deepskyblue4",alpha=0.7)+
                     xlab("")+ylab("")+theme(axis.text=element_text(size=5),
                                             axis.title=element_text(size=5))+
                     theme(plot.margin=unit(c(0.1,0.1,0,0),"cm"),
                           panel.spacing = unit(0,"null")))
    
    pbf1 = theme2 + geom_line(aes(y, lBF, color=Prior, linetype=Prior), size=0.5, data = plotM) +
      annotation_custom(
        grob = g,
        xmin = -0.3,
        xmax = 1,
        ymin = 9,
        ymax = 58.1
      ) +
      geom_hline(yintercept=log(c.beta1/c.alpha1),linetype="dashed",colour="deepskyblue4", size=0.35,alpha=0.7)+
      xlab(bquote("Data " * bar(y))) +ggtitle(bquote("Bayes factor (" * bar(y)[0] * '=0.25)')) +
      theme(legend.position="bottom",legend.title = element_blank(),
            legend.text.align = 0,
            legend.text = element_text(size = 9))+ guides(color = guide_legend(override.aes = list(size = 0.5),ncol=2))+
      guides(shape = guide_legend(override.aes = list(size = 0.5)))+
      ylab("log(BF)")
    
    ##
    prior.m=0.5
    info = c(1,prior.m,s0)
    base = c(1,0,10)
    mix = cbind(c(0.5, prior.m, s0), c(0.5, prior.m,10))
    mix_80 = cbind(c(0.8, prior.m, s0), c(0.2, prior.m, 10))
    mix_unit = cbind(c(0.5, prior.m, s0), c(0.5, prior.m, 1))
    power = c(1, prior.m, s0)
    
    prior.list = list("vague"=base, "info"=info, "mix"=mix,"mix_unit"=mix_unit,"mix_80"=mix_80, "power"=power)
    y_grid=seq(-1,1,by=0.05)
    BF.list=lapply(1:length(prior.list), function(z) 
      sapply(y_grid, function(t) BF(dat=t,n=n.targ,dataPar=sigma,prior=names(prior.list)[z],priorParsWeights=prior.list[[z]],th0=th0)))
    names(BF.list)=names(prior.list)
    
    plotM=data.frame("y"=rep(y_grid,6),
                     "lBF"=c(log(unlist(BF.list))),
                     "Prior"=rep(c("Vague","Informative","Robust mixture","Robust mixture (unit inform.)","Robust mixture (0.8)","EB power"),each=length(y_grid)))
    plotM$Prior=factor(plotM$Prior,levels=c("EB power", "Informative","Robust mixture","Robust mixture (0.8)","Robust mixture (unit inform.)","Vague"))
    
    theme2=ggplot()+ theme_light() + 
      scale_color_viridis_d(begin = 0, end = 0.9,breaks=levels(plotM$Prior),labels = c(expression(paste(pi[a],":EB power")),
                                                                                       expression(paste(pi[a],":Informative")),
                                                                                       expression(paste(pi[a],":Robust mixture")),
                                                                                       expression(paste(pi[a],":Robust mixture (unit inform.)")),
                                                                                       expression(paste(pi[a],":Robust mixture (0.8)")),
                                                                                       expression(paste(pi[a],":Vague")),
                                                                                       expression(paste(pi[s],":Sampling"))))+
      scale_linetype_manual(values=ltys2,labels = c(expression(paste(pi[a],":EB power")),
                                                    expression(paste(pi[a],":Informative")),
                                                    expression(paste(pi[a],":Robust mixture")),
                                                    expression(paste(pi[a],":Robust mixture (unit inform.)")),
                                                    expression(paste(pi[a],":Robust mixture (0.8)")),
                                                    expression(paste(pi[a],":Vague")),
                                                    expression(paste(pi[s],":Sampling"))))
    
    g = ggplotGrob(theme2+
                     theme(plot.background = element_rect(colour = "black")) + 
                     geom_line(aes(y, lBF,  color=Prior, linetype=Prior), size=0.35, data =subset(plotM, y>0 & y<0.4)) +
                     theme(plot.background = element_rect(colour = "black")) + theme(legend.position = "none") + 
                     geom_hline(yintercept=log(c.beta1/c.alpha1), size=0.35,linetype="dashed",colour="deepskyblue4",alpha=0.7)+
                     xlab("")+ylab("")+theme(axis.text=element_text(size=5),
                                             axis.title=element_text(size=5))+
                     theme(plot.margin=unit(c(0.1,0.1,0,0),"cm"),
                           panel.spacing = unit(0,"null")))
    
    pbf2 = theme2 + 
      geom_vline(xintercept=prior.m,linetype="dashed",colour="gray")+ geom_line(aes(y, lBF, color=Prior, linetype=Prior), size=0.5, data = plotM) +
      annotation_custom(
        grob = g,
        xmin = -0.3,
        xmax = 1,
        ymin = 10,
        ymax = 60
      ) +
      geom_hline(yintercept=log(c.beta1/c.alpha1),linetype="dashed",colour="deepskyblue4", size=0.35,alpha=0.7)+
      xlab(bquote("Data " * bar(y))) +ggtitle(bquote("Bayes factor (" * bar(y)[0] * '=0.5)')) +
      theme(legend.position="bottom",legend.title = element_blank(),
            legend.text.align = 0,
            legend.text = element_text(size = 9))+ guides(color = guide_legend(override.aes = list(size = 0.5),ncol=2))+
      guides(shape = guide_legend(override.aes = list(size = 0.5)))+
      ylab("log(BF)")
    
    
    cairo_pdf(file = paste(out_dir,"/","design_",dec,"_",sampOpt,"_",outL,"_BF_Sensitivity.pdf", sep = ""), width=8.3, height=6, pointsize=12)
    grid.arrange(pbf1,pbf2,  nrow = 1)
    graphics.off()
    
    
  }
}
