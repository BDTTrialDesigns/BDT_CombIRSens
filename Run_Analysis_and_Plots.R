# Design ------------------------------------------------------------------

source('Design.R')

#normal outcome - simulation
designRunSim(outcome='n',dec_cost='normalized',max_n=300,nCores=44)
designRunSim(outcome='n',dec_cost='fixed',max_n=300,nCores=44)
designRunSimLSS() #for Figure 2 and S1

#binomial outcome - simulation
designRunSim(outcome='b',dec_cost='normalized',max_n=300,nCores=25)
designRunSim(outcome='b',dec_cost='fixed',max_n=300,nCores=25)

#binomial outcome - data
designRunDataEx(dec_cost='normalized',max_n=300,nCores=20)



# Figures -----------------------------------------------------------------


source('Figures_Fn_rev.R')
source('Design_Fn.R')

#the following functions produce the figures and save the outputs in '~/Figures':

#normal outcome - simulation
createFigure(outcome='n',dec_cost='n',opt_sample='Nopt',type='sim') 
createFigure(outcome='n',dec_cost='n',opt_sample='Nfix',type='sim')
createFigure(outcome='n',dec_cost='n',opt_sample='G',type='sim') 
createFigure(outcome='n',dec_cost='f',opt_sample='Nopt',type='sim')

#binomial outcome - simulation
createFigure(outcome='b',dec_cost='n',opt_sample='Nopt',type='sim')
createFigure(outcome='b',dec_cost='f',opt_sample='Nopt',type='sim')
createFigure(outcome='b',dec_cost='n',opt_sample='G',type='sim')

#binomial outcome - data
createFigure(outcome='b',dec_cost='n',opt_sample='Nopt',type='data')


# Figure S2 -----------------------------------------------------------------

source('SamplPriorStudy.R')

SamplPrirorStudy()

