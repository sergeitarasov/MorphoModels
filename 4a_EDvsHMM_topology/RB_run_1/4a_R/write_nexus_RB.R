

library("tibble")
library("dplyr")
library("plyr")
library("magrittr")


#' write nexus for RevBayes
#'
#' @param x tibble data frame with the first column to be taxa, and all columns must be named
#' @param file file name
#' @param polym.sep separator for polymorphic characters
#'
#' @return nexus file
#' @export
#'
#' @examples
#' x<-add_column(mt[,c(1:4)], rds$ph$tip.label, .before = 1)
#' x<-add_column(mt[,c(4)], rds$ph$tip.label, .before = 1)
#' colnames(x)<-NULL
#' x<-x %>% replace(.== c("0 1"), c("0 and 1") )
#' write_nexus_RevBayes(x, file = "test_nex.txt", polym.sep=" and ")

write_nexus_RevBayes<-function (x, file, polym.sep=" ")
{
  # add libraries
  

  # get taxa and covert them to nex format
  taxa<-x[,1]
  taxa<-c(taxa)[[1]]
  taxa<-paste0("'", taxa, "'")
  
  
  # remove first column -- taxa; colums must be named
  x<-select (x,c(-1) )
  
  # recode polymorphisms
  x<-apply(x,2, function(y) gsub(polym.sep, " ", y))
  # insert brackets
  polym.ids<-which(nchar(x)>1)
  x[polym.ids]<-paste0("(", x[polym.ids], ")")
  
  # Make Nexs matrix
  
  paste(x[1,], collapse="")
  MT<-apply(x, 1, function(y) paste(y, collapse="") )
  MT<-cbind(taxa,MT)
  MT<-apply(MT, 1, function(y) paste(y, collapse=" ") )
  MT<-paste(MT, collapse="\n")
  
  # Write nexus
  NCHAR <- paste("NCHAR=", ncol(x), sep = "")
  NTAX <- paste0("NTAX=", length(taxa))
  DATATYPE <- c("DATATYPE = STANDARD")
  GAP<-"GAP = -"
  MISSING<-"MISSING = ?"
  
  # SYMBOLS
  sym<-x[which(nchar(x)==1)] %>% unique()
  sym<-sym[(sym!="?" & sym!="-")]
  sym<-sort(sym)
  SYMBOLS<-paste0("SYMBOLS = ", "\"", paste(sym, collapse = " "), "\"")
  SYMBOLS<- 'SYMBOLS ="0 1 2 3 4 5 6 7 8 9"'
  
  # Write file
  zz <- file(file, "w")
  cat(file=zz, "#NEXUS\n[Data written by write.nexus.RevBayes, ", date(), 
      "]\n")
  cat(file=zz, paste("BEGIN DATA;\n", "DIMENSIONS", NTAX, NCHAR), ";\n", sep="", append = TRUE)
  cat(file=zz, paste("FORMAT DATATYPE = STANDARD", GAP, MISSING, SYMBOLS), ";\n", sep="", append = TRUE)
  cat(file=zz, "MATRIX\n", sep="", append = TRUE)
  cat(file=zz, MT, sep="", append = TRUE)
  cat(file=zz, "\n;\n END;", append = TRUE)
  close(zz)
}

