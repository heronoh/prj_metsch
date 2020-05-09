#this file has all the procedures performed on Heron's PhD project with Metchnikowia
##Metschnikowia australis Genome assembly
##All available Metschnikowia Genomes prediction


#install.packages("BiocManager")
#a
#BiocManager::install("Rsamtools")
#install.packages("biomartr")
# install.packages("sf")
# install.packages("adegenet")
#install.packages("ape")
library("biomartr")
library("seqinr")
library("readr")
library("stringr")
library("GenomicRanges")
library("Biostrings")
library("Rsamtools")
library("ggplot2")
library("dplyr")
library("ggrepel")

#library("ape")
#read_csv("/home/heron-oh/metschnikowias/2020/metsch_genomes_IDs.list")
#geno_names <- readLines("/home/heron-oh/metschnikowias/2020/metsch_genomes_IDs.list")

working_dir <- c("/home/heron-oh/metschnikowias/prj_metsh")
setwd(working_dir)



#inicializando as tabelas e variáveis de genomas

genomes_dir <- paste0(working_dir,"/data/genomes/ctgs")
list.files(genomes_dir)
#genomes_dir <- c("/home/heron-oh/metschnikowias/2020/genomas/ctgs/")
busco_dir <- paste0(working_dir,"/data/busco/")

busco_dir <- c("/home/heron-oh/metschnikowias/2020/busco/busco_ceta/")
tRNA_dir <- c("/home/heron-oh/metschnikowias/2020/tRNAscan/")
interpro_dir <- c("/home/heron-oh/metschnikowias/2020/interproscan/v39/")
R_dir <- c("/home/heron-oh/metschnikowias/2020/R/")
genomes_file <- list.files(genomes_dir)
genomes_rad <- substr(list.files(genomes_dir), 0,3 )
genomes_path <- paste0(genomes_dir,genomes_file)

#building genomes_tbl whith columns for all genome informations
genomes_tbl <- tibble::tibble(genomes_rad,genomes_file,
                              contigs=numeric(length = length(genomes_rad)),
                          genome_span=numeric(length = length(genomes_rad)),
                          GC_content=numeric(length = length(genomes_rad)),
                          BUSCO_completeness=numeric(length = length(genomes_rad)),
                          genomes_path,
                          num_CDSs=numeric(length = length(genomes_rad)),
                          gene_density=numeric(length = length(genomes_rad)),
                          repeats=numeric(length = length(genomes_rad)),
                          APSs_RAFP=numeric(length = length(genomes_rad)),
                          AFPs_cryo=numeric(length = length(genomes_rad)),
                          AFPs_dbl=numeric(length = length(genomes_rad)),
                          tRNAs=numeric(length = length(genomes_rad)),
                          interpro_ann=numeric(length = length(genomes_rad)),
                          nr_ann=numeric(length = length(genomes_rad)))

#genomes_tbl <- genomes_tbl %>% dplyr::mutate(color_code=character(length = nrow(genomes_tbl)))
#genomes_tbl <- genomes_tbl %>% dplyr::mutate(repeats=character(length = nrow(genomes_tbl)))
#genomes_tbl <- genomes_tbl %>% dplyr::mutate(APSs_RAFP=numeric(length = length(genomes_rad)))
genomes_tbl <- genomes_tbl %>% dplyr::mutate(nr_ann=numeric(length = length(genomes_rad)))
genomes_tbl <- genomes_tbl %>% dplyr::mutate(interpro_ann=numeric(length = length(genomes_rad)))
#color-code
color_code_tab <- readr::read_tsv(file = "/home/heron-oh/metschnikowias/2020/R/genomas.tsv", col_names = FALSE)

colnames(color_code_tab) <- c("clade","color_code","genomes_rad","species name","NCBI genome_ID")

#genomes_tbl <- dplyr::select(genomes_tbl, -c("clade","color_code","species name","NCBI genome_ID"))

genomes_tbl <- dplyr::left_join(genomes_tbl,color_code_tab, by="genomes_rad")


#genomes_tbl <- tibble::add_column(genomes_tbl, num_CDSs="")

#inicializando as tabelas e variáveis de proteomas
proteomes_dir <- c("/home/heron-oh/metschnikowias/2020/genomas/anot/sem_interpro/prot/")
RAFP_results <- c("/home/heron-oh/metschnikowias/2020/AFPs/class/R-AFPpred/results/preds/")
cryoprotect_results <- c("/home/heron-oh/metschnikowias/2020/AFPs/class/cryoprotect/")
proteomes_files <- list.files(proteomes_dir)
proteomes_rads <- substr(list.files(genomes_dir), 0,3 )
proteome_tbl <- tibble::tibble("proteome_rad"="","protein"="")
proteomes_tbl <- tibble::tibble("proteome_rad"=character(),
                                "protein"=character(),"pred_RAFP"=character(),
                                "score_RAFP"=character(),"pred_cryo"=character())


ann_NR_tbl <- read_table(paste0(R_dir,"ann_NR_v38_counts.tab"), col_names = FALSE)
colnames(ann_NR_tbl) <- c("counts","genomes_rad")


#TODO
# continuar com interproscar
#x <- system("ls -l", intern = TRUE)

#y <- system2("ls", c("-l", "-h"), stdout = TRUE, stderr = "output_error.txt")

ncol(genomes_tbl)
nrow(genomes_tbl)


for (geno in 1:nrow(genomes_tbl)) {
   # geno <- 3
  print(geno)
  #   fasta <- Biostrings::readDNAStringSet(genomes_tbl$genomes_path[geno])
  #   contig_num <- length(fasta)
  #   df_temp <- data.frame(x=c(3,3,3), y = c(5,3,2))
  #   as.character(unlist(fasta[10]))
  #   genomes_tbl$GC_content <- df_temp
  #   genomes_tbl$contigs[geno] <- contig_num
  #   names(fasta[10])
  #   names(fasta)
                  # ##GenomicRanges::
                  # geno_fasta <- Biostrings::readDNAStringSet(genomes_tbl$genomes_path[geno])
                  # gen_path <- genomes_tbl$genomes_path[geno]
                  # ##tamanho do genoma
                  # span <- as.numeric(sum(width(geno_fasta)))
                  # genomes_tbl$genome_span[geno] <- span
                  #
                  # ##numero de contigs
                  # ctgs  <- as.numeric(length(width(geno_fasta)))
                  # genomes_tbl$contigs[geno] <- ctgs
                  #
                  # ##GC content
                  #
                  # geno_fasta <- seqinr::read.fasta(file = genomes_tbl$genomes_path[geno],
                  #                     seqtype = "DNA", as.string = FALSE, forceDNAtolower = FALSE)
                  # #seqinr::getName(geno_fasta)
                  # geno_seqs <- seqinr::getSequence(geno_fasta)
                  # geno_seqs <- unlist(geno_seqs, use.names=FALSE)
                  # gc_cont <- seqinr::GC(geno_seqs)
                  # genomes_tbl$GC_content[geno] <-  (100*gc_cont)
                  #
                  #
                  # #BUSCO completeness
                  #
                  # busco_file <- paste0("/home/heron-oh/metschnikowias/2020/busco/short_summary/short_summary_",
                  #                      genomes_tbl$genomes_rad[geno],".txt")
                  #
                  # busco_short <- readLines(busco_file)
                  #
                  # busco_compl <- (as.numeric(stringr::str_extract(stringr::str_extract(busco_short[10],
                  #                 pattern = "\t....\t"), pattern = "\\d\\d\\d\\d"))/1759)*100
                  #
                  # #fazer compl e dupl e comple sing
                  # #busco <- system2( grep 'C:' busco_file)
                  #
                  # genomes_tbl$BUSCO_completeness[geno] <- busco_compl
                  #
                  # #contagem de CDSs preditas
                  #
                  # all_prots <- seqinr::read.fasta(file = paste0(proteomes_dir,genomes_tbl$genomes_rad[geno] ,"_prot.fasta"))
                  # genomes_tbl$num_CDSs[geno] <- length(all_prots)
                  #
                  # #gene density
                  #
                  # genomes_tbl$gene_density[geno] <- ((genomes_tbl$num_CDSs[geno]/genomes_tbl$genome_span[geno])*1000000)
                  #
                  # #repeats
                  #
                  # repeat_tbl <- readLines(paste0("/home/heron-oh/metschnikowias/2020/repeatmasker/",genomes_tbl$genomes_rad[geno],".tbl"))
                  #
                  # repeats <- grep(x = repeat_tbl, pattern = "bases masked",value = TRUE) %>%
                  #  stringr::str_remove(pattern = ".*\\( ") %>%
                  #    stringr::str_remove(pattern = " %\\)")
                  #    genomes_tbl$repeats[geno] <- as.numeric(repeats)
                  #
                  # #tRNAs
                  #
                  # tRNAs_tbl <- readLines(paste0(tRNA_dir,genomes_tbl$genomes_rad[geno],"_tRNAs_stats.txt"))
                  # tRNAs <- as.numeric(grep (x = tRNAs_tbl, pattern = "decoding Standard", value = TRUE) %>%
                  #   stringr::str_remove(pattern = "tRNAs decoding Standard 20 AA:              "))
                  # genomes_tbl$tRNAs[geno] <- tRNAs
                  #  #genomes_tbl$repeats[geno] <-
                  # # stringr::str_match(repeat_tbl,"bases masked" )
                  # # dplyr::filter(repeat_tbl, "bases masked")
                  #
                  #    #proporção do genoma anotado no interproscan
  genomes_tbl$interpro_ann[geno] <-


}

##características do Proteomas e classificação de AFPs

#
#
# proteomes_dir <- c("/home/heron-oh/metschnikowias/2020/genomas/anot/sem_interpro/prot/")
# RAFP_results <- "/home/heron-oh/metschnikowias/2020/AFPs/class/R-AFPpred/results/preds/"
# cryoprotect_results <- "/home/heron-oh/metschnikowias/2020/AFPs/class/cryoprotect/"
# proteomes_files <- list.files(proteomes_dir)
# proteomes_rads <- substr(list.files(genomes_dir), 0,3 )
#
# proteome_tbl <- tibble::tibble("proteome_rad"="","protein"="","size AA"="")
# proteomes_tbl <- tibble::tibble("proteome_rad"=character(),
#                                 "protein"=character(),"pred_RAFP"=character(),
#                                 "score_RAFP"=character(),"pred_cryo"=character())


for(prots in proteomes_rads){
  #prots <- "bic"
  all_prots <- stringr::str_remove(grep("[>]", readLines(paste0(proteomes_dir,prots,
                                          "_prot.fasta")), value = TRUE),"[>]")

  #prots_len <-

  proteome_tbl <- tibble::tibble("proteome_rad"=prots,"protein"=all_prots,)

  RAFP_class_tbl <- read_tsv(paste0(RAFP_results,prots,".txt"),col_names = FALSE)
  colnames(RAFP_class_tbl) <- c("protein","pred_RAFP","score_RAFP")
  RAFP_class_tbl$score_RAFP <- as.character(RAFP_class_tbl$score_RAFP)
  proteome_tbl <- dplyr::left_join(proteome_tbl,RAFP_class_tbl)

  cryoprotect_class_tbl <- read_csv(paste0(cryoprotect_results,prots,".csv"),col_names = TRUE)
  colnames(cryoprotect_class_tbl) <- c("protein","pred_cryo")
  #cryoprotect_class_tbl <- str_remove(cryoprotect_class_tbl$pred_cryo, "on-")
  proteome_tbl <- dplyr::left_join(proteome_tbl,cryoprotect_class_tbl)

  proteomes_tbl <- rbind(proteomes_tbl, proteome_tbl)

}

#proteomes_tbl_bckp <-proteomes_tbl

##genomes_tbl <- genomes_tbl %>% dplyr::mutate(APSs_cryo=numeric(length = length(genomes_rad)))
proteomes_tbl <- proteomes_tbl_bckp %>%
  dplyr::mutate(RAPF_AFP=logical(length = nrow(proteomes_tbl)),
                cryo_AFP=logical(length = nrow(proteomes_tbl)))

#atribuir TRUE or FALSE para o caráter AFP
for (i in 1:nrow(proteomes_tbl)){
  if (proteomes_tbl$pred_RAFP[i] == "AFP"){
    proteomes_tbl$RAPF_AFP[i] <- TRUE
  }else{
    proteomes_tbl$RAPF_AFP[i] <- FALSE
  }
}


for (i in 1:nrow(proteomes_tbl)){
  if (proteomes_tbl$pred_cryo[i] == "AFP"){
    proteomes_tbl$cryo_AFP[i] <- TRUE
  }else{
    proteomes_tbl$cryo_AFP[i] <- FALSE
  }
}
gelo <- proteomes_tbl[,c(6,7)]

gelo %>% filter(cryo_AFP == TRUE) %>% filter(RAPF_AFP == TRUE)

# order(genomes_tbl$genomes_rad)

for (i in  1:length(genomes_tbl$genomes_rad)){

  AFP_tab <- filter(.data = proteomes_tbl, proteomes_tbl$proteome_rad == genomes_tbl$genomes_rad[i])

  RAFP_all <- nrow(AFP_tab)
  RAFP_AFP <- nrow(filter(.data = AFP_tab, AFP_tab$RAPF_AFP == TRUE))
  RAFP_prop <- (RAFP_AFP/RAFP_all)
  genomes_tbl$APSs_RAFP[i] <- RAFP_prop

  cryo_all <- nrow(AFP_tab)
  cryo_AFP <- nrow(filter(.data = AFP_tab, AFP_tab$cryo_AFP == TRUE))
  cryo_prop <- (cryo_AFP/cryo_all)
  genomes_tbl$AFPs_cryo[i] <- cryo_prop

  #double_AFP <- nrow(filter(.data = AFP_tab, AFP_tab$RAPF_AFP == TRUE) %>% filter(.preserve = (AFP_tab$cryo_AFP == TRUE)))
  AFP_temp <- filter(.data = AFP_tab, AFP_tab$RAPF_AFP == TRUE)
  double_AFP <- nrow(filter(.data = AFP_temp, AFP_temp$cryo_AFP == TRUE))
  double_prop <- (double_AFP/cryo_all)
  genomes_tbl$AFPs_dbl[i] <- double_prop
}




######GRÁFICOS##################


save_dir <- "/home/heron-oh/metschnikowias/2020/figs/"
#dir.create(save_dir)
dir.exists(save_dir)
tabela <- genomes_tbl
ordem_esps <- c("sce","cer",
           "con","mer","mat","mar","loc","flo","bor","cub",
           "sma","ipo","mau","ham","haw","pal","kam","pil",
           "dek","bow","sim","cos","col","ari","pro","dra",
           "abe","shi","hib","cla","kip","aus","bic","gol",
           "tor","clv")

# sort(ordem_esps)
# order(ordem_esps)
# rank(ordem_esps)
# tabela[c(5,4), ]


tabela <- tabela[rank(ordem_esps), ]
tabela <- mutate(tabela, name_rad = paste0(tabela$`species name`," (",tabela$genomes_rad,") ") )

tabela$name_rad <- factor(tabela$name_rad, levels = tabela$name_rad)
#tabela$`species name` <- factor(tabela$`species name`, levels = tabela$`species name`)
#tabela$genomes_rad <- factor(tabela$genomes_rad, levels = tabela$genomes_rad)


### Tamanho dos Genomas ######################
graph_name <- "graph_genome_size"
grafico_size <- ggplot(data = tabela,
                       mapping = aes(x = tabela$name_rad,
                       #mapping = aes(x = paste0(tabela$`species name`," (",tabela$genomes_rad,") "),
                                     y = (tabela$genome_span)/1000000)) +
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$genome_span)/1000000),linetype = 2) +
  # geom_hline(yintercept = ((mean(tabela$genome_span)+sd(tabela$genome_span))/1000000)) +
  # geom_hline(yintercept = ((mean(tabela$genome_span)-sd(tabela$genome_span))/1000000)) +
  ylab("Tamanho do genoma (Mb)") +
  xlab("Espécie") +
  ggtitle("Tamanho dos Genomas") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()

grafico_size
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_size, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_size, device = "svg", width =16 , height = 12, dpi = 300)

dev.off()
### Conteúdo GC ######################


graph_name <- "graph_gc_content"
grafico_gc <- ggplot(data = tabela,
                     mapping = aes(x = tabela$name_rad,
                                   y = tabela$GC_content))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$GC_content)),linetype = 2) +
  ylab("Conteúdo GC (%)") +
  xlab("Espécie") +
  ggtitle("Conteudo GC") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_gc
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_gc, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_gc, device = "svg", width =16 , height = 12, dpi = 300)


### CDSs Preditas################


graph_name <- "graph_cds_preds"
grafico_CDSs <- ggplot(data = tabela,
                       mapping = aes(x = tabela$name_rad,
                                   y = tabela$num_CDSs))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$num_CDSs)),linetype = 2) +
  ylab("CDSs preditas") +
  xlab("Espécie") +
  ggtitle("CDSs preditas") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_CDSs
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_CDSs, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_CDSs, device = "svg", width =16 , height = 12, dpi = 300)

### Densidade Gênica################


graph_name <- "graph_gene_dens"
grafico_gene_dens <- ggplot(data = tabela,
                       mapping = aes(x = tabela$name_rad,
                                     y = tabela$gene_density))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$gene_density)),linetype = 2) +
  ylab("CDSs por Mb") +
  xlab("Espécie") +
  ggtitle("Densidade Gênica") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_gene_dens
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_gene_dens, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_gene_dens, device = "svg", width =16 , height = 12, dpi = 300)

### Repeats################


graph_name <- "graph_repeats"
grafico_repeats <- ggplot(data = tabela,
                            mapping = aes(x = tabela$name_rad,
                                          y = tabela$repeats))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$repeats)),linetype = 2) +
  ylab("% de repetições") +
  xlab("Espécie") +
  ggtitle("Conteúdo de repetições por genoma") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_repeats
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_repeats, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_repeats, device = "svg", width =16 , height = 12, dpi = 300)


### Proporção de AFPs por genoma - RAFP-pred################


graph_name <- "graph_RAFP"
grafico_RAFP <- ggplot2::ggplot(data = tabela,
                          mapping = aes(x = tabela$name_rad,
                                        y = tabela$APSs_RAFP))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$APSs_RAFP)),linetype = 2) +

  ylab("% de CDSs classificadas como AFPs") +
  xlab("Espécie") +
  ggtitle("Proporção de AFPs por genoma - RAFP-pred") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_RAFP
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_RAFP, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_RAFP, device = "svg", width =16 , height = 12, dpi = 300)

### Proporção de AFPs por genoma - cryoprotect################


graph_name <- "graph_cryo"
grafico_cryo <- ggplot(data = tabela,
                       mapping = aes(x = tabela$name_rad,
                                     y = tabela$AFPs_cryo))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$AFPs_cryo)),linetype = 2) +
  ylab("% de CDSs classificadas como AFPs") +
  xlab("Espécie") +
  ggtitle("Proporção de AFPs por genoma - Cryoprotect") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_cryo
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_cryo, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_cryo, device = "svg", width =16 , height = 12, dpi = 300)

### Proporção de AFPs por genoma - RAFP-pred & cryoprotect################


graph_name <- "graph_RAFPandCryo"
grafico_RAFPandCRYO <- ggplot(data = tabela,
                       mapping = aes(x = tabela$name_rad,
                                     y = tabela$AFPs_dbl))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$AFPs_dbl)),linetype = 2) +
  ylab("% de CDSs classificadas como AFPs") +
  xlab("Espécie") +
  ggtitle("Proporção de AFPs por genoma \n RAFP-pred & Cryoprotect") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()
grafico_RAFPandCRYO
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_RAFPandCRYO, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_RAFPandCRYO, device = "svg", width =16 , height = 12, dpi = 300)

### Proporção de tRNAs preditos por genoma - tRNAscan-SE################


graph_name <- "graph_tRNAs"
grafico_tRNAs <- ggplot(data = tabela,
                        mapping = aes(x = tabela$name_rad,
                                      y = tabela$tRNAs))+
  geom_bar(stat = "identity", fill = tabela$color_code, width = 0.8) +
  geom_hline(yintercept = (mean(tabela$tRNAs)),linetype = 2) +
  ylab("Número de tRNAs preditos por genoma") +
  xlab("Espécie") +
  ggtitle("Número de tRNAs preditos\npor genoma(tRNAscan-SE)") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 20)) +
  theme(axis.title.y = element_text(size = 20)) +
  coord_flip()

grafico_tRNAs
ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_tRNAs, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_tRNAs, device = "svg", width =16 , height = 12, dpi = 300)

getwd()
###BUSCO#############




### Genome Span X CDSs #########################





graph_name <- "spanXcds"
grafico_spanXcds <-
  ggplot(data = tabela, mapping = aes(x = tabela$genome_span,
                                y = tabela$num_CDSs, group = tabela$`species name`)) +
    xlab("Tamanho do genoma") +
    ylab("Número de CDSs preditas") +
    scale_x_continuous(breaks = c(12500000,15000000,17500000,20000000,22000000),
                       labels =c("12,5 Mb","15 Mb","17,5 Mb","20 Mb","22 Mb")) +
    scale_y_continuous(breaks = c(5100,5500,6000,6500,7000,7500,7400))+
  geom_point(shape = 21, size = 4, color = "white", fill = tabela$color_code) +
  geom_label_repel(data = tabela, label = tabela$`species name`, label.size = 0.4,
                   label.padding = 0.2, segment.alpha = 0.3, nudge_y = 100, force = 0.1) +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(hjust = 1, size = 14, color = "black")) +
    theme(axis.text.y = element_text(hjust = 1, size = 14, color = "black"))

ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_spanXcds, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_spanXcds, device = "svg", width =16 , height = 12, dpi = 300)

####################################################



