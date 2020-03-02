#Extract from orthoMCL groups the counts for genes of each genome and plot one genoma against other

#load libraries
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(plotly)
options(stringsAsFactors = F)





Refazer gráfico de genes compartilhados ordenando pelo numeri de copias no genoma

http://www.transcriptionfactor.org/index.cgi?Home



####################################################################################
##species list
sps <- as.array(c("clv","tor","gol","bic","aus","kip","cla","hib","shi","abe","dra","pro","ari",
                      "col","cos","sim","bow","dek","pil","kam","pal","haw","ham","mau","ipo","sma",
                      "cub","bor","flo","loc","mar","mat","mer","con","cer","sce")
)


#Read orthoMCL out and extract gene counts for genome

orthMCL_file <- c("/home/heron-oh/metschnikowias/2020/orthoMCL/dmd/36metsch_groups.txt")

orthMCL_lines <- readLines(orthMCL_file)

fig_save_dir <- c("/home/heron-oh/metschnikowias/2020/figs/")

#Creat table with orths groups and count species genes in each

orths_tbl <- tibble::tibble(
  Orth_group = as.character(
    str_pad(
      str_remove(
        str_remove(
          str_extract(orthMCL_lines, pattern = "36metsch.+\\:"),
            pattern =c("36metsch")
        ),pattern = "\\:"
      ),width = 4,pad = 0)),
  orths_genes = str_remove(orthMCL_lines, pattern = "36metsch.+\\: "),
  all_genes = as.numeric(0),
  aus = as.numeric(0),
  bic = as.numeric(0),
  abe = as.numeric(0),
  ari = as.numeric(0),
  bor = as.numeric(0),
  bow = as.numeric(0),
  cer = as.numeric(0),
  cla = as.numeric(0),
  clv = as.numeric(0),
  col = as.numeric(0),
  con = as.numeric(0),
  cos = as.numeric(0),
  cub = as.numeric(0),
  dek = as.numeric(0),
  dra = as.numeric(0),
  flo = as.numeric(0),
  gol = as.numeric(0),
  ham = as.numeric(0),
  haw = as.numeric(0),
  hib = as.numeric(0),
  ipo = as.numeric(0),
  kam = as.numeric(0),
  kip = as.numeric(0),
  loc = as.numeric(0),
  mar = as.numeric(0),
  mat = as.numeric(0),
  mau = as.numeric(0),
  mer = as.numeric(0),
  pal = as.numeric(0),
  pil = as.numeric(0),
  pro = as.numeric(0),
  sce = as.numeric(0),
  shi = as.numeric(0),
  sim = as.numeric(0),
  sma = as.numeric(0),
  tor = as.numeric(0),
                            all_other_median = as.numeric(0),
                            all_other_sd = as.numeric(0))



for (group in 1:nrow(orths_tbl)) {
  orths_tbl$all_genes[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "\\|")
  orths_tbl$aus[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "aus\\|")
  orths_tbl$bic[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "bic\\|")
  orths_tbl$all_genes[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "\\|")
  orths_tbl$abe[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "abe\\|")
  orths_tbl$ari[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "ari\\|")
  orths_tbl$bor[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "bor\\|")
  orths_tbl$bow[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "bow\\|")
  orths_tbl$cer[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "cer\\|")
  orths_tbl$cla[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "cla\\|")
  orths_tbl$clv[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "clv\\|")
  orths_tbl$col[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "col\\|")
  orths_tbl$con[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "con\\|")
  orths_tbl$cos[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "cos\\|")
  orths_tbl$cub[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "cub\\|")
  orths_tbl$dek[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "dek\\|")
  orths_tbl$dra[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "dra\\|")
  orths_tbl$flo[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "flo\\|")
  orths_tbl$gol[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "gol\\|")
  orths_tbl$ham[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "ham\\|")
  orths_tbl$haw[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "haw\\|")
  orths_tbl$hib[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "hib\\|")
  orths_tbl$ipo[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "ipo\\|")
  orths_tbl$kam[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "kam\\|")
  orths_tbl$kip[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "kip\\|")
  orths_tbl$loc[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "loc\\|")
  orths_tbl$mar[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "mar\\|")
  orths_tbl$mat[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "mat\\|")
  orths_tbl$mau[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "mau\\|")
  orths_tbl$mer[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "mer\\|")
  orths_tbl$pal[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "pal\\|")
  orths_tbl$pil[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "pil\\|")
  orths_tbl$pro[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "pro\\|")
  orths_tbl$sce[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "sce\\|")
  orths_tbl$shi[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "shi\\|")
  orths_tbl$sim[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "sim\\|")
  orths_tbl$sma[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "sma\\|")
  orths_tbl$tor[group] <- str_count(string = orths_tbl$orths_genes[group], pattern = "tor\\|")
}

#orths_tbl_bckp <- orths_tbl
# orths_tbl <- orths_tbl_bckp

tabela_long <- 0

#Creat table in the LONG format to use median and sd functions (observations as rows (avoid NAs))
tabela_long <- tidyr::gather(orths_tbl, key = "species", value = "gene_num", -c(Orth_group, orths_genes, all_genes))

#tabela_long_bckp <- tabela_long
#tabela_long <- tabela_long_bckp


# fazer a média para as ouras espécies
tabela_long <- tabela_long_bckp %>%
  filter(!(species %in% c("aus","bic","gol"))) %>%
  group_by(Orth_group, all_genes) %>%
  # summarise(all_other_mean = ifelse(mean(gene_num) > 0 & mean(gene_num) < 1,1,mean(gene_num))) %>%
  summarise(all_other_mean = round(mean(gene_num))) %>%
  # summarise(all_other_mean = ifelse(mean(gene_num) > 0,1,mean(gene_num)), all_other_mad = stats::mad(gene_num)) %>%
  ungroup()

unique(tabela_long_bckp$gene_num)
mean(tabela_long_bckp$gene_num)

# add aus & bic observations to final table
orths_mean_tbl <- left_join(tabela_long, dplyr::select(orths_tbl, c(Orth_group, aus, bic, gol)), by = "Orth_group")
# sort final table by number of genes in cluster and MAD
orths_mean_tbl <- arrange(orths_mean_tbl, all_genes, all_other_mean)
# fix final table order at current order
orths_mean_tbl <- mutate(orths_mean_tbl, Orth_group = factor(Orth_group, levels = Orth_group))


#BACKUP
orths_mean_tbl_bckp <- orths_mean_tbl
#orths_mean_tbl <- orths_mean_tbl_bckp

# Filter informative groups: those with different orths numbers in aus and bic
#so para o LOLLIPOP PLOT
      # #n'ao precisa, quero todos
      #         # orths_mean_tbl <- orths_mean_tbl_bckp %>%
      #         #   # filter(all_other_mad != 0) %>%
      #         #   filter(aus != bic) %>%
      #         #   filter(aus != 0) %>%
      #         #   filter(bic != 0)
      #
      # # Calculate median +- MAD for plot median deviation
      # #so para o LOLLIPOP PLOT
      #             # orths_mean_tbl <- orths_mean_tbl_bckp %>%
      #             #   mutate(mad_plus = all_other_median + all_other_mad) %>%
      #               # mutate(mad_minus = all_other_median - all_other_mad)
      #
      # # Assign only >0 values to MAD min values
      # for (line in 1:nrow(orths_mean_tbl)){
      #   if (orths_mean_tbl$mad_minus[line] <0) {
      #     orths_mean_tbl$mad_minus[line] = 0
      #   }
      # }



#tidyr::gather()

#library(viridis)
#library(MASS)


# Plot orths count for each group in aus, bic and the median for the other genomes
            plot_shared_paralogs <-  orths_mean_tbl %>%
              ggplot() +
              geom_hline(yintercept = 0, linetype = "dashed", color ="grey") +
              geom_errorbar(aes(x = Orth_group, ymin = mad_minus, ymax = mad_plus), alpha = 0.6, width = 0.5) +
              geom_segment(aes(x = Orth_group,  xend = Orth_group, y = all_other_median, yend = aus), alpha = 0.4, color = "navyblue") +
              geom_segment(aes(x = Orth_group,  xend = Orth_group, y = all_other_median, yend = bic), alpha = 0.4, color = "dodgerblue3") +
              geom_point(aes(x = Orth_group, y = all_other_median), shape = 23, color = "black", fill = "gray", size = 5) +
              geom_point(aes(x = Orth_group, y = aus), color = "black", fill = "navyblue", size = 3, alpha = 0.8, shape = 21) +
              geom_point(aes(x = Orth_group, y = bic), color = "black", fill = "dodgerblue3", size = 3, alpha = 0.8, shape = 21) +
              scale_y_continuous(trans = 'sqrt', breaks = scales::pretty_breaks(n = max(orths_mean_tbl$bic, orths_mean_tbl$aus))) +
              xlab(label = "ID do grupo de ortólogos") +
              ylab(label = "Número de genes por genoma") +
              ggtitle(label = "Número de cópias parálogas em \nM. australis, M. bicuspidata e demais genomas\n") +
              theme_classic(base_size = 35) +
              theme(plot.title = element_text(hjust = 0.5)) +
              theme(axis.text.x = element_text(hjust = 1, size = 16, color = "black")) +
              theme(axis.text.y = element_text(hjust = 1, size = 16, color = "black", face = "italic")) +
              theme(axis.title.x = element_text(size = 28)) +
              theme(axis.title.y = element_text(size = 28)) +
              coord_flip()

            # mutate(orths_tbl_bckp, density = (orths_tbl_bckp$aus, orths_tbl_bckp$bic))





### plot de paralogos com cópias diferentes em aus e bic


graph_name <- "paralogs_aus_bic"

        # orths_mean_tbl%>%
  # filter(aus != 0 ) %>%
  #   filter(bic != 0) %>%
para_aus_bic <-
  orths_mean_tbl %>%
    filter((bic + aus) >=1 ) %>%
      ggplot( aes( y = aus, x = bic, label = Orth_group)) +
        theme_bw(base_size = 23) +
        scale_fill_gradientn(
          name = "Número de\ncluster\nde ortólogos", guide = "legend",trans = "log1p",
          colors =  c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082"),
          breaks = c(0,1,2,5,10,25,50,100,500,1000,2000,4000),
          labels = c("0","1","2-5","6-10","11-25","26-50","51-100","101-500","501-1000","1001-2000","2001-4000","> 4000") ) +
   # stat_binhex(binwidth = c(1,1)) +
        stat_bin2d( binwidth = c(0.2,0.2)) +
        coord_fixed( xlim = c(0,(max(orths_mean_tbl$bic)),
                     ylim = c(0,(max(orths_mean_tbl$aus)+1)))) +
        geom_abline(slope = 1, alpha = 0.3) +
        xlab( label = "Número de genes por cluster em M. bicuspidata") +
        ylab( label = "Número de genes por cluster em M. australis") +
        scale_x_continuous(breaks = c(unique(orths_mean_tbl$bic))) +
        scale_y_continuous(breaks = c(unique(orths_mean_tbl$aus))) +
        ggtitle( label = "Relação do número de genes por cluster de ortólogos \nentre os genomas de M. australis e M. bicuspidata") +
        theme( plot.title = element_text(size = 35, hjust = 0.5)) +
        theme( axis.text.x = element_text(hjust = 1, size = 23, color = "black")) +
        theme( axis.text.y = element_text(hjust = 1, size = 23, color = "black", face = "italic")) +
        theme( axis.title.x = element_text(size = 28)) +
        theme( axis.title.y = element_text(size = 28))
para_aus_bic

ggsave(filename =  paste0(save_dir,graph_name,".pdf"), plot = para_aus_bic, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,graph_name,".svg"), plot = para_aus_bic, device = "svg", width =16 , height = 12, dpi = 300)



### plot de paralogos com cópias diferentes de 0 em aus e gol
graph_name <- "paralogs_aus_gol"

# orths_mean_tbl%>%
# filter(aus != 0 ) %>%
#   filter(bic != 0) %>%
para_aus_gol <-
  orths_mean_tbl%>%
  filter((gol + aus) >=1 ) %>%
  ggplot( aes(y=aus, x=gol,label =Orth_group)) +
  theme_bw(base_size = 23) +
  scale_fill_gradientn(
    name = "Número de\ncluster\nde ortólogos", guide = "legend",trans = "log1p",
    colors =  c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082"),
    breaks = c(0,1,2,5,10,25,50,100,500,1000,2000,4000),
    labels = c("0","1","2-5","6-10","11-25","26-50","51-100","101-500","501-1000","1001-2000","2001-4000","> 4000") ) +
  # stat_binhex(binwidth = c(1,1)) +
  stat_bin2d(binwidth = c(0.35,0.35)) +
  coord_fixed(xlim = c(0,(max(orths_mean_tbl$gol)),
              ylim = c(0,(max(orths_mean_tbl$aus)+1)))) +
  geom_abline(slope = 1, alpha = 0.3) +
  xlab(label = "Número de genes por cluster em M. golubevii") +
  ylab(label = "Número de genes por cluster\nem M. australis") +
  scale_x_continuous(breaks = c(unique(orths_mean_tbl$gol))) +
  scale_y_continuous(breaks = c(unique(orths_mean_tbl$aus))) +
  ggtitle(label = "Relação do número de genes por cluster de ortólogos \n entre os genomas de M. australis e M. golubevii") +
  theme(plot.title = element_text(size = 35, hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 23, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 23, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 28)) +
  theme(axis.title.y = element_text(size = 28))
para_aus_gol


ggsave(filename =  paste0(save_dir,graph_name,".pdf"), plot = para_aus_gol, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,graph_name,".svg"), plot = para_aus_gol, device = "svg", width =16 , height = 12, dpi = 300)





#############################################################
### paralogos bic e gol
### plot de paralogos com cópias diferentes de 0 em aus e gol
graph_name <- "paralogs_bic_gol"

# orths_mean_tbl%>%
# filter(bic != 0 ) %>%
#   filter(bic != 0) %>%
para_bic_gol <-
  orths_mean_tbl%>%
  filter((gol + bic) >=1 ) %>%
  ggplot( aes(y=bic, x=gol,label =Orth_group)) +
  theme_bw(base_size = 23) +
  scale_fill_gradientn(
    name = "Número de\ncluster\nde ortólogos", guide = "legend",trans = "log1p",
    colors =  c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082"),
    breaks = c(0,1,2,5,10,25,50,100,500,1000,2000,4000),
    labels = c("0","1","2-5","6-10","11-25","26-50","51-100","101-500","501-1000","1001-2000","2001-4000","> 4000") ) +
  # stat_binhex(binwidth = c(1,1)) +
  stat_bin2d(binwidth = c(0.35,0.35)) +
  coord_fixed(xlim = c(0,(max(orths_mean_tbl$gol)),
              ylim = c(0,(max(orths_mean_tbl$bic)+1)))) +
  geom_abline(slope = 1, alpha = 0.3) +
  xlab(label = "Número de genes por cluster em M. golubevii") +
  ylab(label = "Número de genes por cluster em M. bicuspidata") +
  scale_x_continuous(breaks = c(unique(orths_mean_tbl$gol))) +
  scale_y_continuous(breaks = c(unique(orths_mean_tbl$bic))) +
  ggtitle(label = "Relação do número de genes por cluster de ortólogos \n entre os genomas de M. bicuspidata e M. golubevii") +
  theme(plot.title = element_text(size = 35, hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 23, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 23, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 28)) +
  theme(axis.title.y = element_text(size = 28))
para_bic_gol


ggsave(filename =  paste0(save_dir,graph_name,".pdf"), plot = para_bic_gol, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,graph_name,".svg"), plot = para_bic_gol, device = "svg", width =16 , height = 12, dpi = 300)

############################################################
### paralogos aus e media outras
graph_name <- "paralogs_aus_MED"

# orths_mean_tbl%>%
# filter(aus != 0 ) %>%
#   filter(bic != 0) %>%
para_aus_MED <-
  orths_mean_tbl%>%
  filter((all_other_mean + aus) >=1 ) %>%
  ggplot( aes(y=aus, x=all_other_mean,label =Orth_group)) +
  theme_bw(base_size = 23) +
  scale_fill_gradientn(
    name = "Número de\ncluster\nde ortólogos", guide = "legend",trans = "log1p",
    colors =  c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082"),
    breaks = c(0,1,2,5,10,25,50,100,500,1000,2000,4000),
    labels = c("0","1","2-5","6-10","11-25","26-50","51-100","101-500","501-1000","1001-2000","2001-4000","> 4000") ) +
  # stat_binhex(binwidth = c(1,1)) +
  stat_bin2d(binwidth = c(0.15,0.15)) +
  coord_fixed(xlim = c(0,(max(orths_mean_tbl$all_other_mean)),
              ylim = c(0,(max(orths_mean_tbl$aus)+1)))) +
  geom_abline(slope = 1, alpha = 0.3) +
  xlab(label = "Número de genes por cluster em média das Metschnikowia EG") +
  ylab(label = "Número de genes por cluster em M. australis") +
  scale_x_continuous(breaks = c(unique(orths_mean_tbl$all_other_mean))) +
  scale_y_continuous(breaks = c(unique(orths_mean_tbl$aus))) +
  ggtitle(label = "Relação do número de genes\npor cluster de ortólogos entre os genomas de\nM. australis e a média das Metschnikowia E-G") +
  theme(plot.title = element_text(size = 35, hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 23, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 23, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 28)) +
  theme(axis.title.y = element_text(size = 28))
para_aus_MED


ggsave(filename =  paste0(save_dir,graph_name,".pdf"), plot = para_aus_MED, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,graph_name,".svg"), plot = para_aus_MED, device = "svg", width =16 , height = 12, dpi = 300)




############################################################
###paralogos bic e media outras

graph_name <- "paralogs_bic_MED"

# orths_mean_tbl%>%
# filter(bic != 0 ) %>%
#   filter(bic != 0) %>%
para_bic_MED <-
  orths_mean_tbl%>%
  filter((all_other_mean + bic) >=1 ) %>%
  ggplot( aes(y=bic, x=all_other_mean,label =Orth_group)) +
  theme_bw(base_size = 23) +
  scale_fill_gradientn(
    name = "Número de\ncluster\nde ortólogos", guide = "legend",trans = "log1p",
    colors =  c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082"),
    breaks = c(0,1,2,5,10,25,50,100,500,1000,2000,4000),
    labels = c("0","1","2-5","6-10","11-25","26-50","51-100","101-500","501-1000","1001-2000","2001-4000","> 4000") ) +
  # stat_binhex(binwidth = c(1,1)) +
  stat_bin2d(binwidth = c(0.3,0.3)) +
  coord_fixed(xlim = c(0,(max(orths_mean_tbl$all_other_mean)+5)),
              ylim = c(0,(max(orths_mean_tbl$bic)+1))) +
  geom_abline(slope = 1, alpha = 0.3) +
  xlab(label = "Número de genes por cluster em média das Metschnikowia E-G") +
  ylab(label = "Número de genes por cluster em M. bicuspidata") +
  scale_x_continuous(breaks = c(unique(orths_mean_tbl$all_other_mean))) +
  scale_y_continuous(breaks = c(unique(orths_mean_tbl$bic))) +
  ggtitle(label = "Relação do número de genes\npor cluster de ortólogos entre os genomas\nde M. bicuspidata e a média das Metschnikowia E-G") +
  theme(plot.title = element_text(size = 35, hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 23, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 23, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 28)) +
  theme(axis.title.y = element_text(size = 28))
para_bic_MED


ggsave(filename =  paste0(save_dir,graph_name,".pdf"), plot = para_bic_MED, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,graph_name,".svg"), plot = para_bic_MED, device = "svg", width =16 , height = 12, dpi = 300)


###########################################################
### gol e media
graph_name <- "paralogs_gol_MED"

# orths_mean_tbl%>%
# filter(gol != 0 ) %>%
#   filter(gol != 0) %>%
para_gol_MED <-
  orths_mean_tbl%>%
  filter((all_other_mean + gol) >=1 ) %>%
  ggplot( aes(y=gol, x=all_other_mean,label =Orth_group)) +
  theme_bw(base_size = 23) +
  scale_fill_gradientn(
    name = "Número de\ncluster\nde ortólogos", guide = "legend",trans = "log1p",
    colors =  c("#FF0000","#FF7F00","#FFFF00","#00FF00","#0000FF","#4B0082"),
    breaks = c(0,1,2,5,10,25,50,100,500,1000,2000,4000),
    labels = c("0","1","2-5","6-10","11-25","26-50","51-100","101-500","501-1000","1001-2000","2001-4000","> 4000") ) +
  # stat_binhex(binwidth = c(1,1)) +
  stat_bin2d(binwidth = c(0.5,0.5)) +
  coord_fixed(xlim = c(0,(max(orths_mean_tbl$all_other_mean)+5)),
              ylim = c(0,(max(orths_mean_tbl$gol)+1))) +
  geom_abline(slope = 1, alpha = 0.3) +
  xlab(label = "Número de genes por cluster\nem média nas Metschnikowia E-G") +
  ylab(label = "Número de genes por cluster em M. golubevii") +
  scale_x_continuous(breaks = c(unique(orths_mean_tbl$all_other_mean))) +
  scale_y_continuous(breaks = c(unique(orths_mean_tbl$gol))) +
  ggtitle(label = "Relação do número de genes por cluster\nde ortólogos entre os genomas de\nM. golubevii e a média das Metschnikowia E-G") +
  theme(plot.title = element_text(size = 35, hjust = 0.5)) +
  theme(axis.text.x = element_text(hjust = 1, size = 23, color = "black")) +
  theme(axis.text.y = element_text(hjust = 1, size = 23, color = "black", face = "italic")) +
  theme(axis.title.x = element_text(size = 28)) +
  theme(axis.title.y = element_text(size = 28))
para_gol_MED


ggsave(filename =  paste0(save_dir,graph_name,".pdf"), plot = para_gol_MED, device = "pdf", width =16 , height = 12, dpi = 300)
ggsave(filename =  paste0(save_dir,graph_name,".svg"), plot = para_gol_MED, device = "svg", width =16 , height = 12, dpi = 300)










###########################################################

  ggsave(paste0(fig_save_dir,"aus_bic_shared_paralogs.svg"), device = 'svg',
         plot = plot_shared_paralogs, width = 32, height = 32, units = "cm")







################################################################################
  library("UpSetR")
  UpSetR::upset()



upset_tbl <-  select(orths_tbl_bckp, c(Orth_group,orths_genes))


upset_tbl <- mutate(.data = upset_tbl, genomes_gens="")


for (cluster in 1:nrow(upset_tbl)) {

  upset_tbl$genomes_gens[cluster] <- str_remove_all(upset_tbl$orths_genes[cluster], pattern = "\\|[:alnum:]++")

}


upset_tbl_orths <- orths_tbl_bckp[,c(1,3:39)]

#upset_tbl_orths_bckp <- upset_tbl_orths

upset_tbl_orths <- upset_tbl_orths_bckp %>%
  filter(aus != bic) %>%
  filter(aus != 0) %>%
  filter(bic != 0)





for (col in 3:ncol(upset_tbl_orths)) {
  upset_tbl_orths[,col][upset_tbl_orths[,col]>=1] <- 1

}


upset_tbl_orths$Orth_group <- as.factor(upset_tbl_orths$Orth_group)

upset_tbl_orths <- as.data.frame(upset_tbl_orths)

pdf(file="/home/heron-oh/metschnikowias/2020/figs/UPSET_TESTE.pdf", onefile=FALSE) # or other device

 upset(upset_tbl_orths, sets = rev(c("clv","tor","gol","bic","aus","kip","cla","hib",
                                 "shi","abe","dra","pro","ari","col","cos","sim",
                                 "bow","dek","pil","kam","pal","haw","ham","mau",
                                 "ipo","sma","cub","bor","flo","loc","mar","mat",
                                 "mer","con","cer","sce")),
      nsets = 36,
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "black",
      matrix.color = "black", sets.bar.color = "blue", mainbar.y.label = "Número de clusters compartilhados",
      sets.x.label = "Número de genes", shade.color = "gray", point.size = 1,
      text.scale = c(2, 1.8, 2, 1.8, 2, 2), line.size = 0.5,keep.order = TRUE,group.by = "degree")




dev.off()





movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )




upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=scatter_plot, x="ReleaseDate", y="AvgRating"),
                                                         list(plot=scatter_plot, x="ReleaseDate", y="Watches"),list(plot=scatter_plot, x="Watches", y="AvgRating"),
                                                         list(plot=histogram, x="ReleaseDate")), ncols = 2))



