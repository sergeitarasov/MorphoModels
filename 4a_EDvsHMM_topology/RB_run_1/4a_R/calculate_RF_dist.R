source('dep.R')

#   ____________________________________________________________________________
#   MJ cons trees                                                           ####

##  ............................................................................
##  Read trees,    Run 1                                                         ####

rb.out <- '../../Puhti/RB_run_1/output'

mj.M3 <- redTreesDir(rb.path=file.path(rb.out, 'out_ED_M3'), pattern="*mjcon.tre$")
length(mj.M3)
mj.M6 <- redTreesDir(rb.path=file.path(rb.out, 'out_ED_M6'), pattern="*mjcon.tre$")
mj.Ind <- redTreesDir(rb.path=file.path(rb.out, 'out_SMM_ind'), pattern="*mjcon.tre$")
mj.SW <- redTreesDir(rb.path=file.path(rb.out, 'out_SMM_sw_M6'), pattern="*mjcon.tre$")



##  ............................................................................
##  Reoder original tree to match rb ones                                   ####

# original trees
trees <- read.tree('data/trees.tre')
nm0 <- paste0('r', 1:1000)


# reoder trees
nm3 <- str_match(names(mj.M3), "out_r\\s*(.*?)\\s*_")[,2] %>% as.numeric()
length(nm3)
unique(nm3) %>% length()
nm6 <-str_match(names(mj.M6), "out_r\\s*(.*?)\\s*_")[,2] %>% as.numeric()
nm.ind <- str_match(names(mj.Ind), "out_r\\s*(.*?)\\s*_")[,2] %>% as.numeric()
nm.sw <- str_match(names(mj.SW), "out_r\\s*(.*?)\\s*_")[,2] %>% as.numeric()
all(nm3==nm6)
all(nm3==nm.ind)
all(nm3==nm.sw)

# reoder
trees.order <- trees[nm3]

# check roedering
# properly reodered
calcRF(trees.order, mj.M6) %>% table
# not reodered
calcRF(trees, mj.M6) %>% table


##  ............................................................................
##  Run 1. Table                                                            ####

mj.tb1 <- tibble(
  #id=paste0('r', nm3),
  id=nm3,
  M3=calcRF(trees.order, mj.M3),
  M6=calcRF(trees.order, mj.M6),
  Ind=calcRF(trees.order, mj.Ind),
  SW=calcRF(trees.order, mj.SW)
)

mj.tb1
apply(mj.tb1.final[,2:5], 2, table)

##  ............................................................................
##  Read trees, Run 2                                                        ####
rb.out <- '../../Puhti/RB_run_2/output'

mj.M3 <- redTreesDir(rb.path=file.path(rb.out, 'out_ED_M3'), pattern="*mjcon.tre$")
mj.M6 <- redTreesDir(rb.path=file.path(rb.out, 'out_ED_M6'), pattern="*mjcon.tre$")
mj.Ind <- redTreesDir(rb.path=file.path(rb.out, 'out_SMM_ind'), pattern="*mjcon.tre$")
mj.SW <- redTreesDir(rb.path=file.path(rb.out, 'out_SMM_sw_M6'), pattern="*mjcon.tre$")

mj.tb2 <- tibble(
  #id=paste0('r', nm3),
  id=nm3,
  M3=calcRF(trees.order, mj.M3),
  M6=calcRF(trees.order, mj.M6),
  Ind=calcRF(trees.order, mj.Ind),
  SW=calcRF(trees.order, mj.SW)
)

# Both runs should give the same table
apply(mj.tb1[,2:5], 2, table)
# > apply(mj.tb1[,2:5], 2, table)
# M3  M6 Ind  SW
# 0 499 503 481 477
# 1 231 224 226 225
# 2 207 212 223 225
# 3  46  46  45  50
# 4  17  15  25  23
apply(mj.tb2[,2:5], 2, table)
# > apply(mj.tb2[,2:5], 2, table)
# M3  M6 Ind  SW
# 0 503 509 474 476
# 1 223 216 233 227
# 2 209 214 225 225
# 3  48  47  48  50
# 4  17  14  20  22

#saveRDS(mj.tb1, 'data/mj.tb1.rds')
#saveRDS(mj.tb2, 'data/mj.tb2.rds')

# Proportion are different between Run1 and Run2. This happens because some clades may have suppprt very close to .5
# Let's remove different results
rr <- which((mj.tb1$M6-mj.tb2$M6)!=0)
mj.tb1[rr,]

# different rows
diff <- c(
which((mj.tb1$M6-mj.tb2$M6)!=0),
which((mj.tb1$M3-mj.tb2$M3)!=0),
which((mj.tb1$Ind-mj.tb2$Ind)!=0),
which((mj.tb1$SW-mj.tb2$SW)!=0)
)
diff <- unique(diff)

mj.tb1[diff,]

# remove different rows
mj.tb1.final <- mj.tb1[-diff,]
mj.tb2.final <- mj.tb2[-diff,]

apply(mj.tb1.final[,2:5], 2, table)
apply(mj.tb2.final[,2:5], 2, table)

# Nww results are the same
# > apply(mj.tb1.final[,2:5], 2, table)
# M3  M6 Ind  SW
# 0 474 479 458 456
# 1 193 185 191 193
# 2 174 181 187 185
# 3  34  33  35  37
# 4  16  13  20  20


saveRDS(mj.tb1.final, 'data/mj_tb1_final.rds')


##  ............................................................................
##  Bootstrap                                                               ####


# Proceed to Bootstrap
rstudioapi::navigateToFile("Bootstrap.R")



##  ............................................................................
##  MCC TREE                                                               ####

##  ............................................................................
##  Read trees                                                              ####
rb.out <- '../../Puhti/RB_run_1/output'

mcc.M3 <- redTreesDir(rb.path=file.path(rb.out, 'out_ED_M3'), pattern="*mcc.tre$")
length(mcc.M3)
mcc.M6 <- redTreesDir(rb.path=file.path(rb.out, 'out_ED_M6'), pattern="*mcc.tre$")
mcc.Ind <- redTreesDir(rb.path=file.path(rb.out, 'out_SMM_ind'), pattern="*mcc.tre$")
mcc.SW <- redTreesDir(rb.path=file.path(rb.out, 'out_SMM_sw_M6'), pattern="*mcc.tre$")


##  ............................................................................
##  MCC TREE                                                           ####

# calcRF(trees.order, mcc.M3) %>% table

mcc.tb <- tibble(
  #id=paste0('r', nm3),
  id=nm3,
  M3=calcRF(trees.order, mcc.M3),
  M6=calcRF(trees.order, mcc.M6),
  Ind=calcRF(trees.order, mcc.Ind),
  SW=calcRF(trees.order, mcc.SW)
)

apply(mcc.tb[,2:5], 2, table)

# > apply(mcc.tb[,2:5], 2, table)
# M3  M6 Ind  SW
# 0 614 612 606 602
# 2 318 316 310 316
# 4  68  72  84  82

