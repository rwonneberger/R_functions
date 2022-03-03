# Commonly required libraries
library(plyr)
library(data.table)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggfortify)
library(viridis) 
library(Hmisc)
library(janitor)
library(xtable)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(stringr)
library(gplots)
library(reshape2)
library(magrittr)

library(ggpubr)

options("width"=200)


# Functions

# transpose a df and fix the rownames and colnames
transpose_df<-function(df){
  df<-as.data.frame(df)
  df_names <-  df[,1]
  df.T <- as.data.frame(as.matrix(t(df[,-1])))
  colnames(df.T) <- df_names
  
  return(df.T)
}

# transpose a df and fix the rownames and colnames. Add the rownames as the first column. Name first col is given by second function arg
transpose_df_Col1<-function(df, Col1_name){
  df<-as.data.frame(df)
  df_names <-  df[,1]
  df.T <- as.data.frame(as.matrix(t(df[,-1])))
  colnames(df.T) <- df_names
  df.T$Col1<-rownames(df.T)
  df.T%<>%relocate(Col1)
  names(df.T)[1]<-Col1_name
  
  return(df.T)
}


# Round all numeric columns in a dataframe
round_df <- function(df, digits_to_keep) {
  numeric_columns <- sapply(df, mode) == 'numeric'
  df[numeric_columns] <-  round(df[numeric_columns], digits_to_keep)
  
}




# Reformat the per location heritability file from the META-R package to a useable format
format_metar_heri_perloc <- function(df) {
  traits<-df$V1[!is.na(df$V1)] #Make a vector of all traits in the df
  envs<-length(unique(df$Environment[!is.na(df$Environment)] )) +2 #Make a vector with all traits, repeated x times with x being number of environments plus 2 (because of the empty rows inserted between each trait)
  trait_vec<-rep(traits, each=envs) 
  df$Trait<-trait_vec #Add vector to df
  df<-df[!is.na(df$Environment)] #Remove rows with trait names in V1
  df$V1<-NULL #Remove columns
  df$V6<-NULL
  df %<>% relocate(c(Trait, Environment, Heritability, `Genotype Variance`, `Residual Variance`)) #Sort remaining cols
  return(df)
}

# Reformat the across location heritability file from the META-R package to a useable format
format_metar_heri_acrossloc <- function(df) {
  df<-as.data.frame(t(df)) #Transform df and change to data.frame
  df <- data.frame(Trait = row.names(df), df$V1)  #Overwrite df: Add rownames as col and select only h2 column (V1)
  df%<>%filter(Trait != "Environment", Trait!="Statistic") #Remove rows without h2
  df<-df[!is.na(df$df),] #Remove all the BLUE rows. No h2 was calculated for them so we don't want them
  rownames(df)<-NULL #remove row names
  names(df)<-c("Trait", "Heritability") #Rename columns
  return(df)
}







