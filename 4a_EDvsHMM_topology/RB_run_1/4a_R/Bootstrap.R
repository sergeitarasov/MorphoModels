
# To compare the proportions of correctly recovered trees between the original and other models, 
# I used the proportion differential and bootstrapped p-values. 
# The proportion differential $pd$ was calculated as $pd=Prop(ED_6)-Prop(M_i)$, 
# where $Prop$ stands for the proportion of the correctly recovered tress and $M_i$ is one of the models. 
# The p-value for $pd$ was estimated by bootstrapping with 5000 replications.
# 
mj.tb1.final <- readRDS('data/mj_tb1_final.rds')

mj.tb1.final
N.repl <- 5000 
dt.size<- nrow(mj.tb1.final)

# matrix to store pd metrics
pd <- matrix(NA, N.repl, 3)
colnames(pd) <- c('pd.M3', 'pd.Ind', 'pd.SW')
#pd

# Boostraping
for (i in 1:N.repl){
  ss <- sample(1:dt.size, dt.size, replace = TRUE)
  fr <- apply(mj.tb1.final[ss,2:5], 2, table)
  pd[i,1] <- fr[1,2]-fr[1,1]
  pd[i,2] <- fr[1,2]-fr[1,3]
  pd[i,3] <- fr[1,2]-fr[1,4]

}

pd
hist(pd[,1])

#-------Calculate p-values
apply(pd, 2, function(x) length(which(x<=0))/N.repl)

# > apply(pd, 2, function(x) length(which(x<=0))/N.repl)
# pd.M3  pd.Ind pd.SW 
# 0.1902 0.0042 0.0036 

# Proportions of correctly recovered trees
apply(mj.tb1.final[,2:5], 2, table)
# > apply(mj.tb1.final[,2:5], 2, table)
#    M3  M6 Ind  SW
# 0 474 479 458 456
# 1 193 185 191 193
# 2 174 181 187 185
# 3  34  33  35  37
# 4  16  13  20  20

(apply(mj.tb1.final[,2:5], 2, table)/dt.size) %>% round(2)
#    M3   M6  Ind   SW
# 0 0.53 0.54 0.51 0.51
# 1 0.22 0.21 0.21 0.22
# 2 0.20 0.20 0.21 0.21
# 3 0.04 0.04 0.04 0.04
# 4 0.02 0.01 0.02 0.02

