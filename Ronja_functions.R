# Commonly required libraries
library(plyr)
library(data.table)
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
library(tidyr)
library(readr)
library(purrr)
library(forcats)

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

# Calculate pvals (http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram)
# Function for all-by-all calculation of cor
cor.mtest.all <- function(mat, method, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = method, ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# Function for calculation of cor of only certain columns: The first ones are metabolites, the last ones are the phenotypes
cor.mtest.selected <- function(df, last_metab, method, ...) {
  mat <- as.matrix(df)
  nrows <- last_metab
  ncols <- ncol(df)-last_metab
  p.mat<- matrix(NA, nrows, ncols)
  for (i in 1:last_metab) {
    for (j in 1:ncols) {
      tmp <- cor.test(mat[, i], mat[, last_metab+j], method = method, ...)
      p.mat[i, j] <-  tmp$p.value
    }
  }
  colnames(p.mat) <- colnames(df)[(last_metab+1):ncol(df)] 
  rownames(p.mat) <- colnames(df)[1:last_metab]
  p.mat
}




# Round all numeric columns in a dataframe
round_df <- function(df, digits_to_keep) {
  numeric_columns <- sapply(df, mode) == 'numeric'
  df[numeric_columns] <-  round(df[numeric_columns], digits_to_keep)
  
}


#Look at the first few rows and cols of a df/dt
peek<-function(df, x=10, y=10){
  df[1:x, 1:y]
  }


#Look at the last few rows and cols of a df/dt
peekend<-function(df, x=10, y=10){
  df[(nrow(df)-x):nrow(df), (ncol(df)-x):ncol(df)]
  }



# Custom Manhattan plot function with facet_grid
Man<-function(df, x, y, chr, trait = NULL){
  arguments<-as.list(match.call())
  x = eval(arguments$x, df)
  y = eval(arguments$y, df)
  chr = eval(arguments$chr, df)
  if(missing(trait)) {
  df1<-data.frame(Pos = x/1000000, logp = -log10( y), chrom = chr)
  ggplot(df1, aes(x=Pos, y=logp)) + theme_bw() + geom_point() +  labs(x="Position (Mbp)",y="-log10(pval)") +   theme(legend.position = "none") + theme(axis.text.x =element_text(angle = 90)) + facet_grid(~chrom, scales="free_x", space = "free_x")
  } else {
    Trait = eval(arguments$trait, df)
    df1<-data.frame(Pos = x/1000000, logp = -log10( y), chrom = chr, Trait=Trait)
    ggplot(df1, aes(x=Pos, y=logp)) + theme_bw() + geom_point() +  labs(x="Position (Mbp)",y="-log10(pval)") +   theme(legend.position = "none") + theme(axis.text.x =element_text(angle = 90)) + facet_grid(Trait~chrom, scales="free_x", space = "free_x")
  }
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







