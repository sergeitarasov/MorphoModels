
hist2matrix <- function(hist, Nch=10, recode.from=c('T*', 'TA*',  'TAr', 'TAb'), recode.to=c(0:3), nex.folder, nex.file='r1_10ch.nex'){
  
  #-- Remove constant chars
  mt <- sapply(hist, function(x) x$states)
  uu <- apply(mt, 2, function(x) length(unique(x)) )
  #table(uu)
  mt <- mt[,which(uu>1)]
  
  # select first Nch chars
  mt <- mt[,1:Nch]
  
  #-- Recode chars
  Qta
  #nex <- mapvalues(mt, c('T*', 'TA*', 'TAr', 'TAb'), c('0 1 2 3', '4 5', '6', '7') )
  nex <- mapvalues(mt, recode.from, recode.to )
  #nex
  
  #--- Write Nex
  taxa <- rownames(nex)
  nex <- as_tibble(nex)
  nex <- add_column(nex, taxa, .before = 1)
  
  path <- file.path(nex.folder, nex.file)
  
  write_nexus_RevBayes(nex, file = path, polym.sep=" and ")
}

redTreesDir <- function(rb.path, pattern="*mcc.tre$"){
  fls <- list.files(path=rb.path, pattern = pattern, recursive = TRUE, full.names = F)
  out.trees <- list()
  
  for (i in fls){
    print(i)
    out.trees[[i]] <- read.nexus(file.path(rb.path, i))
  }
  
  return(out.trees)
}


calcRF <- function(trees1, trees2){
  rf <- c()
  for (i in 1:length(trees1)){
    rf.i <- RF.dist(trees1[[i]], trees2[[i]], normalize = F, check.labels = TRUE, rooted = T)
    rf <- c(rf, rf.i)
  }
  return(rf)
}





