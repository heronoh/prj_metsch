#Creat gephi csv input files for ortholog network from relations of BLASTp,
#DIAMOND or orthoMCL

######Extract orthologs data from fmt6 and creat gephi csv input files for
######  ortholog network from relations of BLASTp, DIAMOND or orthoMCL

###### AUTHOR: Heron Hilario - heronoh@gmail.com


#load libraries
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(sys)
library(purrr)
library(readr)
library(tidyverse)
install.packages("tidygraph")
