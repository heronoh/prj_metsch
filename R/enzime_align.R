library(readr)
library(stringr)
library(tibble)
library(dplyr)
library(reshape2)
# a partir de grupos de ortólogos do orthoMCL, extrair a anotação de cada

#testando o git




sp2remove <- c("abe","cer","col","con","cub","dek","dra","flo","haw",
               "kam","mar","mau","mer","pil","shi","sim","bor","bow",
               "cla","ari","cos","ham","hib","ipo","loc","mat","pal",
               "pro","sce","sma","kip","clv")
sp2keep <- c("aus","bic","gol","tor")

allSp <- sort(c(sp2keep,sp2remove))

#
#         ordem_esps <- c("sce","cer","con","mer","mat","mar","loc","flo","bor","cub",
#                         "sma","ipo","mau","ham","haw","pal","kam","pil","dek","bow",
#                         "sim","cos","col","ari","pro","dra","abe","shi","hib","cla",
#                         "kip","aus","bic","gol","tor","clv")
#
#         # sort(ordem_esps)
#         # order(ordem_esps)
#         # rank(ordem_esps)
#         # tabela[c(5,4), ]
#
#
#         tabela <- tabela[rank(ordem_esps), ]
#         tabela <- mutate(tabela, name_rad = paste0(tabela$`species name`," (",tabela$genomes_rad,") ") )




### I. - contar quantas vezes cada espécie aparece num cluster, para identificar expansões parálogas

clusters <- read_tsv(file = ("/home/heron-oh/metschnikowias/2020/paralogs/36lev_9380_groups_orthoMCL.list"),
                     col_names = FALSE )

clusters <- dplyr::as_tibble(clusters)
clusters$X1 <- str_replace(string = clusters$X1, pattern = "36metsch", replacement = "clust_" )
clusters$X1 <- str_replace(string = clusters$X1, pattern = " genes", replacement = "")
clusters$X1 <- str_replace(string = clusters$X1, pattern = " taxa", replacement = "")
clusters$X1 <- str_replace(string = clusters$X1, pattern = ":", replacement = "")

colnames(clusters) <- c("clust_ID","clust_genes")

#clusters_bckp <- clusters
#clusters <- clusters_bckp
clusters <- mutate(.data = clusters, num_genes = as.numeric(0), num_taxa = as.numeric(0) )




## add columns to count number of genes in each specie

clusters <- mutate(.data = clusters,
aus = 0, bic = 0, gol = 0, tor = 0,
abe = 0, cer = 0, col = 0, con = 0,
cub = 0, dek = 0, dra = 0, flo = 0,
haw = 0, kam = 0, mar = 0, mau = 0,
mer = 0, pil = 0, shi = 0, sim = 0,
bor = 0, bow = 0, cla = 0, ari = 0,
cos = 0, ham = 0, hib = 0, ipo = 0,
loc = 0, mat = 0, pal = 0, pro = 0,
sce = 0, sma = 0, kip = 0, clv = 0)
# nao consegui com o for
# for (sp in allSp) {
# new_col <- paste(sp)
# clusters <- mutate(.data = clusters, new_col = 0 )
# }


## Count cluster num genes and num taxa

for (clus in 1:nrow(clusters)) {

clusters$num_genes[clus] <-
  str_split_fixed(
    string = str_remove(
      string = str_remove(
        string = clusters$clust_ID[clus],
        pattern = "clust_[[:digit:]]++\\("),
      pattern = "\\)"),
      pattern = ",", n = 2)[,1]

clusters$num_taxa[clus] <-
  str_split_fixed(
    string = str_remove(
      string = str_remove(
        string = clusters$clust_ID[clus],
        pattern = "clust_[[:digit:]]++\\("),
      pattern = "\\)"),
    pattern = ",", n = 2)[,2]

clusters$clust_ID[clus] <- str_replace(string = clusters$clust_ID[clus],
                                      pattern = "\\([:alnum:]++,[:alnum:]++\\)",
                                      replacement = "")
}


for (group in 1:nrow(clusters)) {
  clusters$aus[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(aus\\)")
  clusters$bic[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(bic\\)")
  clusters$abe[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(abe\\)")
  clusters$ari[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(ari\\)")
  clusters$bor[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(bor\\)")
  clusters$bow[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(bow\\)")
  clusters$cer[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(cer\\)")
  clusters$cla[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(cla\\)")
  clusters$clv[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(clv\\)")
  clusters$col[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(col\\)")
  clusters$con[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(con\\)")
  clusters$cos[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(cos\\)")
  clusters$cub[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(cub\\)")
  clusters$dek[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(dek\\)")
  clusters$dra[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(dra\\)")
  clusters$flo[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(flo\\)")
  clusters$gol[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(gol\\)")
  clusters$ham[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(ham\\)")
  clusters$haw[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(haw\\)")
  clusters$hib[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(hib\\)")
  clusters$ipo[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(ipo\\)")
  clusters$kam[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(kam\\)")
  clusters$kip[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(kip\\)")
  clusters$loc[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(loc\\)")
  clusters$mar[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(mar\\)")
  clusters$mat[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(mat\\)")
  clusters$mau[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(mau\\)")
  clusters$mer[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(mer\\)")
  clusters$pal[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(pal\\)")
  clusters$pil[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(pil\\)")
  clusters$pro[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(pro\\)")
  clusters$sce[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(sce\\)")
  clusters$shi[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(shi\\)")
  clusters$sim[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(sim\\)")
  clusters$sma[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(sma\\)")
  clusters$tor[group] <- str_count(string = clusters$clust_genes[group], pattern = "\\(tor\\)")
}





## select only single copy orthologs

clusters$num_genes <- as.numeric(clusters$num_genes)
clusters$num_taxa <- as.numeric(clusters$num_taxa)

dplyr::filter(.data = clusters, clusters$num_genes == 36, clusters$num_taxa == 36)
orths_1p1_tbl <- dplyr::filter(.data = clusters, clusters$num_genes == clusters$num_taxa)

paralogs <- dplyr::filter(.data = clusters, clusters$num_genes != clusters$num_taxa)



paralogs %>%
  filter(aus != 0) %>%
  filter(bic != 0) %>%
  filter(gol != 0)
##count num genes in species








#thomaz net
# teste <- ~if_else(. > 0, 1, 0)
# class(teste)

PA_matrix <- clusters %>%
  select(c(1,5:40)) %>%
  mutate_at(.vars = c(2:ncol(.)),~if_else(. > 0, 1, 0)) %>%
  mutate(rowsums = rowSums(.[,c(2:ncol(.))])) %>%
  filter(rowsums < ncol(.)-2) %>%
  select(-rowsums)

sps <- colnames(PA_matrix)[-1]

sp_sp_genes <- as.data.frame(matrix(nrow = 36,ncol = 36))
colnames(sp_sp_genes) <- sps
rownames(sp_sp_genes) <- sps
for (i in 1:length(sps)){
  for (j in 1:length(sps)){
    sp1 <- sps[i]
    sp2 <- sps[j]
    col_1 <- which(colnames(PA_matrix)==sp1)
    col_2 <- which(colnames(PA_matrix)==sp2)
    genes_1 <- PA_matrix %>%
      dplyr::select(1,col_1) %>%
      rename(sp=2) %>%
      filter(sp>0) %>%
      pull(clust_ID) %>%
      as.character()
    genes_2 <- PA_matrix %>%
      dplyr::select(1,col_2) %>%
      rename(sp=2) %>%
      filter(sp>0) %>%
      pull(clust_ID) %>%
      as.character()
    comm <- intersect(genes_1,genes_2)
    universe <- length(PA_matrix$cluster_ID)
    # pval <- round(phyper(q = length(comm)-1,
    #                      m = length(genes_1),
    #                      n = universe - length(genes_1),
    #                      k = length(genes_2),
    #                      lower.tail = F),
    #               digits = 100000)
    sp_sp_genes[i,j] <- length(comm)
  }
}

sp_sp_melt <- sp_sp_genes %>%
  rownames_to_column("sp") %>%
  melt() %>%
  filter(sp!=variable) %>%
  rename(Source=sp,Target=variable,ncomm=value)

nodes <- color_code
nodes$Label <- nodes$clade
colnames(nodes)[1] <- "Id"

write.csv(nodes,"/home/heron-oh/metschnikowias/2020/net/dmd/100/sp_sp_nodes.csv",quote = F,row.names = F)
write.csv(sp_sp_melt,"/home/heron-oh/metschnikowias/2020/net/dmd/100/sp_sp_edges.csv",quote = F,row.names = F)


mutate(.data = clusters, genes = clust_counts[1,], genes = clust_counts[,2])

sps <- c(as.character(allSp))
for (sp in sps){
  print(sp)
  # flag <- as.array(paste0(" -e '",as.character(sp),"'"))
  # sps <- append(x = sps, values =  flag)
}









sps <- c(as.character(allSp))
for (sp in sps){
  flag <- as.array(paste0(" -e '",as.character(sp),"'"))
  sps <- append(x = sps, values =  flag)
}

grep_sps_out <- c("4levs")
sp_out_command <- paste0("grep -v ", paste(sps, collapse = "")," ",fmt6_file,
                         " > ", trab_dir,grep_sps_out,".fmt6")
