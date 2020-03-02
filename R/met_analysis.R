#this script loads gene-gene files from Blast results between
#orthologous genes of 36 yeast species (Metschnikowia and Clavispora) to:
# - manipulate files and filter out interactions by hit %
# - create sp-gene-gene network
# - export propor files to be analyzed in Gephi

#necessary libraries and options:
library(dplyr)
options(stringsAsFactors = F)

#working directory
dir <- "~/Área de Trabalho/GitHub_projects/met_analysis/"
setwd(dir)

#load original edges and nodes files
edges <- read.csv("data/all_edges.csv")
nodes <- read.csv("data/all_nodes.csv")

#make vector with species codes
sps <- unique(substring(c(edges$Source,edges$Target), 1, 3))

#load color codes for species in differente clades
clade_colors <- read.csv("data/clade_colors_2.csv")
clade_colors$color.VARCHAR <- toupper(clade_colors$color)
colnames(clade_colors)[2] <- "color VARCHAR"

#filter edges files to contain only gene-gene hist with 100% score
#make new nodes object with remaining genes
#create column with clade and color information
edges_02 <- edges %>%
  filter(Source!=Target) %>%
  filter(score >= 100) %>%
  select(Source,Target)

nodes_02 <- nodes %>%
  filter(Id %in% unique(c(edges_02$Source,edges_02$Target)))

nodes_02$clade <- substring(nodes_02$Id,1,3)

nodes_02 <- merge(nodes_02,clade_colors,
                  by.x="clade",by.y="sp")

#write edges and nodes to .csv file to be used in Gephi
write.csv(edges_02,file = "~/Área de Trabalho/heron_analysis/all_edges_02.csv",row.names = F)
write.csv(nodes_02,file = "~/Área de Trabalho/heron_analysis/all_nodes_02.csv",row.names = F)

#create sp-gene-gene edges and nodes
genes <- unique(c(edges_02$Source,edges_02$Target))

edges_sp <- data.frame(Source=genes,
                       Target=substring(genes,1,3))

nodes_sp <- data.frame(Id=sps,Label=sps,Class="sps")
nodes_sp$clade <- substring(nodes_sp$Id,1,3)
nodes_sp <- merge(nodes_sp,clade_colors,
                  by.x="clade",by.y="sp")

edges_03 <- as.data.frame(rbind(edges_02,edges_sp))
nodes_03 <- as.data.frame(rbind(nodes_02,nodes_sp))

#write edges and nodes to .csv file to be used in Gephi
write.csv(edges_03,file = "~/Área de Trabalho/heron_analysis/all_edges_03.csv",row.names = F)
write.csv(nodes_03,file = "~/Área de Trabalho/heron_analysis/all_nodes_03.csv",row.names = F)
