source('dep.R')



#   ____________________________________________________________________________
#   Simulate trees                                                          ####


trees<-pbtree(n=4, scale=1, b=.1, d=0, nsim=1000)
#write.tree(trees, file = 'data/trees.tre')
plot(trees)
plot(trees[[3]])

#   ____________________________________________________________________________
#   Sim Histories                                                           ####

trees <- read.tree('data/trees.tre')

#-- Rate Matrix ED_M6
q <- 1
C1 <- initQ(c("r", "b"), q)
A <- initQ(c('A*', 'A'), q)
Tl <- initQ(c('T*', 'T'), q)

step1 <- amaED(A, C1, type = "ql")
Qta <- amaED(Tl, step1, type = "ap")

Qta <- Qta/-diag(Qta)
Qta
Qta <-Qta*.3
Qta
#expm(Qta*100)
# > Qta
# T*  TA*  TAr  TAb
# T*  -0.3  0.3  0.0  0.0
# TA*  0.1 -0.3  0.1  0.1
# TAr  0.1  0.1 -0.3  0.1
# TAb  0.1  0.1  0.1 -0.3


pi <- rep(1/4,4) 
names(pi) <- colnames(Qta)

#----- Simulate hist
Nch=10 # chars to select into matrix
nex.folder <- '../data' # folder for matrix
nex.folder.exp <- '../data_exp' # folder for expanded chars
hist <- list()
i=1
for (i in 1:length(trees)){
  hist.i <- sim.history(trees[[i]], Qta, nsim=100, direction=c("row_to_column"), anc=pi)
  
  fl <- paste0('r', i, '_', Nch, 'ch.nex')
  print(paste0('Sim ', i))
  hist2matrix(hist.i, Nch=Nch, recode.from=c('T*', 'TA*',  'TAr', 'TAb'), recode.to=c(0:3), nex.folder, nex.file=fl)
  
  hist2matrix(hist.i, Nch=Nch, recode.from=c('T*', 'TA*',  'TAr', 'TAb'), recode.to=c('0 1 2 3', '4 5', '6', '7'), nex.folder.exp, nex.file=fl)
  
  hist[[i]] <- hist.i
}


#saveRDS(hist, file = 'data/hist.rds')



#   ____________________________________________________________________________
#   Get RF distance                                                         ####


# Proceed to calculating RF distances
rstudioapi::navigateToFile("calculate_RF_dist.R")






