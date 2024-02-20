#!/usr/bin/env Rscript
#This script was written on 2021_09_21 by Lauren Wasson. It will take in your keyfile and make other keyfiles based on the identifiers
##Usage: Rscript CreateObjectKeyfiles.R <working directory> <keyfile> 

library(dplyr)
args = commandArgs(trailingOnly=TRUE)

setwd(args[1])
keyfile<- read.csv(args[2], fill = TRUE, header = FALSE)
samples<- as.character(unique(keyfile$V1))
objects<- as.character(keyfile$V2)
dictionary<- character(0)

for (i in 1:length(objects)){
  my_string_split <- strsplit(objects, " ")[i]
  dictionary<- append(dictionary,my_string_split)
}
objects<- unique(unlist(dictionary))

for (i in 1:length(objects)) {
  identifier <- objects[i]
  finame<- paste(identifier,"_keyfile.txt", sep="")
  df<- keyfile%>% 
    filter(grepl(identifier, V2))
  write.table(df[1], file = finame, row.names = FALSE, sep = ",", quote = FALSE, col.names  = FALSE )
}
