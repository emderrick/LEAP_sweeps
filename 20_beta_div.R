library(tidyverse)

setwd("/Users/emma/Documents/GitHub/LEAP_sweeps/")
kraken_files <- list.files("data files", pattern = "beta_div_", full.names = T)

for(i in 1: length(kraken_files)){
  kraken_beta <- as.data.frame(read_tsv(kraken_files[i], skip = 18))
  rownames(kraken_beta) <- kraken_beta$x
  kraken_beta <- kraken_beta[, c(2:19)]
  level <- kraken_files[i] %>% str_remove("data files/beta_div_")
  level <- level %>% str_remove(".txt")
  
  for(i in 1:17){
    row = 1 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:16){
    row = 2 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:15){
    row = 3 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:14){
    row = 4 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:13){
    row = 5 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:12){
    row = 6 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:11){
    row = 7 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:10){
    row = 8 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:9){
    row = 9 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:8){
    row = 10 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:7){
    row = 11 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:6){
    row = 12 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:5){
    row = 13 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:4){
    row = 14 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:3){
    row = 15 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:2){
    row = 16 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  for(i in 1:1){
    row = 17 + i
    col = 0 + i
    kraken_beta[row,col] = kraken_beta[col,row]
  }
  
  
  write.csv(kraken_beta, paste("data files/", level, "_beta_div.csv", sep = ""), row.names = F)

}