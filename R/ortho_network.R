#Creat gephi csv input files for ortholog network from relations of BLASTp,
#DIAMOND or orthoMCL

######Extract orthologs data from fmt6 and creat gephi csv input files for
######  ortholog network from relations of BLASTp, DIAMOND or orthoMCL

###### AUTHOR: Heron Hilario - heronoh@gmail.com


#load libraries
library(dplyr)
library(tibble)
library(stringr)
options(stringsAsFactors = F)

# .libPaths()

###################### USER DEFINED VARIABLES ##################################

## tabular format output from DIAMOND ou BLAST ou UCLUST, raw
fmt6_file <- c("/home/heron-oh/metschnikowias/2020/orthoMCL/dmd_sim/36levs_dmd_sim.fmt6")

#qseqid sseqid |ppos| pident length mismatch gapopen qstart qend sstart send evalue bitscore
# ppos = percent of positive scoring = % similarity
# in this script ppos has to be the 3rd column

## dir where the intermeiate files and the output  will be saved`
trab_dir <- ("/home/heron-oh/metschnikowias/2020/net/dmd_sim/13levs/")

## min alignment lenght to keep in the analysis
min_align <- 60

## min ID percent to keep in the analysis
min_id <- 85

## min qcov & dbcov (both the same value) to keep in the analysis & query/db fasta file
min_qcov <- 0.7
dmd_qry_db_faa <- c("/home/heron-oh/metschnikowias/2020/orthoMCL/dmd/goodProteins.fasta" )

## color code #HEX relation for clades
color_code_file <- c("/home/heron-oh/metschnikowias/2020/net/dmd/20levs/clade_colors.csv")

## define species (clades) to keep or to remove from analysis
sp2remove <- as.array(c("abe","cer","col","con","cub","dek","dra","flo","haw",
                        "kam","mar","pil","shi","bor","bow","cla","cos","ham",
                        "mat","pal","sce","sma","clv"))
sp2keep <- c("aus","bic","gol","tor","mer","loc","ipo","mau","sim","ari","pro","hib","kip")

length(sp2keep)

                      # sp2remove <- as.array(c("abe","cer","col","con","cub","dek","dra","flo","haw",
                      #                         "kam","mar","mau","mer","pil","shi","sim","bor","bow",
                      #                         "cla","ari","cos","ham","hib","ipo","loc","mat","pal",
                      #                         "pro","sce","sma","kip","clv"))
                      # sp2keep <- c("aus","bic","gol","tor")



##creat node "ortlhog" to connect with genes that have link to other genes?
creat_orth_node <- FALSE

##value to assign to speacies-gene interaction strength (recommended = 100)
sp_gene_interation <- c(100)

##remove genes that only match with themselves
remove_auto_hits <- TRUE


################################################################################
## PART I - removing unwanted species and hits

### I.1 - Creating work dir if needed

if (dir.exists(trab_dir)== FALSE){
  dir.create(trab_dir, showWarnings = TRUE, recursive = TRUE)
  print(paste0("the dir ",trab_dir," was created"))
} else{
  print(paste0("the dir ",trab_dir," was not created because it already exists"))
}
setwd(trab_dir)

### I.2 - Assigning color code to genes and species

color_code <- read.csv(color_code_file)
colnames(color_code) <- c("clade","color VARCHAR")
color_code$`color VARCHAR` <- tolower(color_code$`color VARCHAR`)

### I.3 - remove species and their genes from fmt6 file for lighter netwrok

sps <- c(as.character())
for (sp in sp2remove){
  flag <- as.array(paste0(" -e '",as.character(sp),"'"))
  sps <- append(x = sps, values =  flag)
}

grep_sps_out <- c("13levs")
sp_out_command <- paste0("grep -v ", paste(sps, collapse = "")," ",fmt6_file,
                                          " > ", trab_dir,grep_sps_out,".fmt6")

### !!!! UNCOMMENT to run !!!!
#    system(sp_out_command)

### I.4 - remove alignments shorter than min_align

awk_align_out <- paste0(grep_sps_out,"_algn", min_align, ".fmt6")

# $3 is now similarity%, $4 is now id%
align_out_command <- c(paste0("awk '($5 >= ",min_align,"){print}' ",trab_dir,
                              grep_sps_out,".fmt6", " > ", trab_dir, awk_align_out))


# align_out_command <- c(paste0("awk '($4 >= ",min_align,"){print}' ",trab_dir,
#                               grep_sps_out,".fmt6", " > ", trab_dir, awk_align_out))

### !!!! UNCOMMENT to run !!!!
#    system(align_out_command)

#### I.4.a - count and print number of lines removed after min_align filter

lines_bef <- str_split(string = (system(paste0("wc -l ",trab_dir,grep_sps_out,".fmt6"), intern = TRUE)), pattern = "[ ]", n = 2, simplify = TRUE )[,1]
lines_aft <- str_split(string = (system(paste0("wc -l ",trab_dir,awk_align_out), intern = TRUE)), pattern = "[ ]", n = 2, simplify = TRUE )[,1]

lines_removed <- (as.numeric(lines_bef) - as.numeric(lines_aft))

print(paste0("---> ",lines_removed," lines (hits) were removed for having alignment span shorter than ",min_align,"bps"))

### I.5 - remove relations with ID% (sim%) smaller than min_id


awk_id_out <- paste0(grep_sps_out,"_algn", min_align,"_id",min_id,".fmt6")
awk_id_command <- c(paste0("awk '($3 >= ",min_id,"){print}' ",trab_dir,awk_align_out, " > ", trab_dir, awk_id_out))

### !!!! UNCOMMENT to run !!!!
#   system( awk_id_command )


#### I.5.a - count and print number of lines removed after min_id filter

lines_bef <- str_split(string = (system(paste0("wc -l ",trab_dir,awk_align_out), intern = TRUE)), pattern = "[ ]", n = 2, simplify = TRUE )[,1]
lines_aft <- str_split(string = (system(paste0("wc -l ",trab_dir,awk_id_out), intern = TRUE)), pattern = "[ ]", n = 2, simplify = TRUE )[,1]
lines_removed <- (as.numeric(lines_bef) - as.numeric(lines_aft))

print(paste0("---> ",lines_removed," lines (hits) were removed for having ID % smaller than ",min_id))


### I.6 - remove hits with reciprocal coverage smaller than min_align

m8f_out <- paste0(grep_sps_out,"_algn", min_align,"_id",min_id,"_cov",
                  (100*min_qcov),".fmt6")

m8filter_comand <- paste0("/home/heron-oh/scripts/blastm8_filter.pl -b ",
                          trab_dir, awk_id_out, " --qcov ", min_qcov," --query ",
                          dmd_qry_db_faa, " --dcov ",min_qcov," --db ",
                          dmd_qry_db_faa, " > ", trab_dir,m8f_out  )

### !!!! UNCOMMENT to run !!!!
# system(m8filter_comand)

#### I.6.a - count and print number of lines removed after coverage filter

lines_bef <- str_split(string = (system(paste0("wc -l ",trab_dir,awk_id_out), intern = TRUE)), pattern = "[ ]", n = 2, simplify = TRUE )[,1]
lines_aft <- str_split(string = (system(paste0("wc -l ",trab_dir,m8f_out), intern = TRUE)), pattern = "[ ]", n = 2, simplify = TRUE )[,1]
lines_removed <- (as.numeric(lines_bef) - as.numeric(lines_aft))

print(paste0("---> ",lines_removed," lines (hits) were removed for having query or db coverage smaller than ",min_qcov))


### I.7 - remove or not genes that only had matches with themselves

if (remove_auto_hits == TRUE){

  auto_hit_out <- paste0(grep_sps_out,"_algn", min_align,"_id",min_id,"_qcov",
                         (100*min_qcov),"_noAuto",".fmt6")
  auto_hit_cmd <- paste0("awk '($1 != $2){print}' ",trab_dir,m8f_out," > ",
                         trab_dir,auto_hit_out)
  ### !!!! UNCOMMENT to run !!!!
   system(auto_hit_cmd)

  csv_input <- paste0(trab_dir,auto_hit_out)
  cat(paste0("Genes that only match with themselves (exclusive genes) were REMOVED.
             \n\n Working on file: \n",csv_input))

}else{

  csv_input <- paste0(trab_dir,m8f_out)
  cat(paste0("Genes that only match with themselves (exclusive genes) were KEPT.
             \n\n Working on file: \n",csv_input))

}

################################################################################

## PART II - generate nodes and edges .csv files for gehpi

### II.1 - loading fmt6 into table and removing all but the first 3 columns

fmt6_final <- csv_input

#tab_orths <- read.table("/home/heron-oh/metschnikowias/2020/net/dmd/6levs/6levs_70cov_id80_dmd.fmt6")
tab_orths <- read.table(fmt6_final)

#removing columns, except queryID, subjectID, ID%
tab_orths <- tab_orths[,1:3]

colnames(tab_orths) <- c("Source","Target","score")

#removendo os radicais da espécie, só precisa pro meu resultado do dmd
tab_orths[,1] <- str_split(string = tab_orths[,1], pattern = "[|]", n=2, simplify = TRUE)[,2]
tab_orths[,2] <- str_split(string = tab_orths[,2], pattern = "[|]", n=2, simplify = TRUE)[,2]


### II.2 - generating all_nodes.csv ("clade","Id","Label","Class","color VARCHAR")


all_nodes_genes1 <- unique(tab_orths[,1])
all_nodes_genes2 <-  unique(tab_orths[,2])
all_nodes_genes <- unique(c(all_nodes_genes1,all_nodes_genes2))



all_nodes <- tibble::tibble(clade="",
                            Id=as.character(all_nodes_genes),
                            Label=as.character(all_nodes_genes),
                            Class="gene")


all_nodes$clade <- substring(all_nodes$Id,1,3)


all_nodes <-  dplyr::left_join(all_nodes,color_code, by="clade") #deu um warning

color_code$`color VARCHAR` <- as.character(color_code$`color VARCHAR`)
# Warning message:
# Column `clade` joining character vector and factor, coercing into character vector


### !!!! make backups of all_nodes table so you won't have to run everything
###      till here if you make mistakes downwards !!!!

#       all_nodes_bkp <- all_nodes
#       all_nodes <- all_nodes_bkp


### II.2.a add species as nodes

sp_line <- c("")

for (sp in sp2keep){
  sp_line <- c(sp, sp, sp, "species", color_code$`color VARCHAR`[color_code$clade==sp])
  all_nodes <- rbind(all_nodes,sp_line)
  print(sp_line)
}



### II.2.b create ortholog node if requested

if (creat_orth_node == TRUE){

  orth_node  <- tibble::tibble("clade" = "homology","Id" = "ortholog","Label" = "ortholog" , "Class" = "homology", "color VARCHAR" = "#FFFFFF")
  all_nodes <- dplyr::bind_rows(all_nodes, orth_node)
  print("a node 'ORTHOLOG' was created to center genes with orthologs identified")
}else{
  print("no ortholog/exclusive nodes were created")
}



# all_nodes <- select(all_nodes, -"colo VARCHAR")
# all_nodes <- filter(all_nodes, !is.na(`color VARCHAR`))
tail(all_nodes, n = 15)
head(all_nodes)



# add colum for node size , species and genes

 all_nodes <- all_nodes %>%
  mutate(size = if_else(Class == "gene",2,20))



### II.2.c write the first output file, "all_nodes"

all_nodes_out <- c(paste0(grep_sps_out,"_algn", min_align,"_id",min_id,"_cov",
                          (100*min_qcov),"_all_nodes.csv"))
write.csv(all_nodes,file = paste0(trab_dir,all_nodes_out), row.names = FALSE, quote = FALSE )

# Warning message:
#   In write.csv(all_nodes, file = paste0(trab_dir, all_nodes_out),  :
#                  attempt to set 'col.names' ignored

### II.3 - generating all_edges.csv ("Source","Target","score")
##### all_edges.csv already is tab_orths, with some lines missing

### !!!! make backups of all_edges table so you won't have to run everything
###      till here if you make mistakes downwards !!!!

tail(tab_orths, n = 15)
head(tab_orths)

   #     tab_orths_bkp <- tab_orths
#        tab_orths <- tab_orths_bkp


### II.3.a - appending 'gene-species-relation' lines to all_edges


querys <- unique(tab_orths$Source)
subjs <- unique(tab_orths$Target)
queNsub <- sort(unique(c(querys,subjs)))


for (gene in queNsub){

  sp_rad <- substr(gene,0,3)
  gene_sp_line <- c(gene,sp_rad,sp_gene_interation)
  tab_orths <- rbind(tab_orths,gene_sp_line)
  print(gene_sp_line)
}


### II.3.b - create ortholog node interaction (edge) with all genes, if requested

if (creat_orth_node == TRUE){

  tab_orths_uniq <- tibble::tibble("Source" = as.character(),"Target" =as.character(),"score" = as.numeric())
  # for (line in 1:10){
  for (line in 1:nrow(tab_orths)){
    if (tab_orths$Source[line] != tab_orths$Target[line]) {
      orth_line <- tibble::tibble("Source" = tab_orths$Source[line], "Target" = "ortholog", "score" = 100)
      tab_orths_uniq <- dplyr::bind_rows(tab_orths_uniq, orth_line)
    }
  }

  tab_orths_uniq <- dplyr::distinct(tab_orths_uniq)
  tab_orths_uniq$score <- as.numeric(tab_orths_uniq$score)
  #tab_orths <- filter(tab_orths, Source != "ortholog")
  tab_orths <- dplyr::bind_rows(tab_orths, tab_orths_uniq)
}


tail(tab_orths,n = 15)
head(tab_orths)

# readr::read_rds()
# saveRDS(all_nodes, "all_nodes.rds")
# all_nodes <- readRDS("all_nodes.rds")


### II.3.c write the second output file, "all_edges.csv"


all_edges_out <- c(paste0(grep_sps_out,"_algn", min_align,"_id",min_id,"_cov",
                          (100*min_qcov),"_all_edges.csv"))
write.csv(tab_orths,file = paste0(trab_dir,all_edges_out), row.names = FALSE, quote = FALSE )


### III.4 - PARSER END - load .csv files into gephi and generate network (FORCE ATLAS 2)


