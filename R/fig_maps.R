.libPaths()
install.packages("maps")
library(ggplot2)
library(dplyr)
library(readr)
require(maps)
require(viridis)
library(ggrepel)
# theme_set(
# theme_void()
# )


geno_coord_tbl <- read_tsv("/home/heron-oh/metschnikowias/2020/R/genomas_coord.tsv", col_names = FALSE)

lat_long_tbl <- geno_coord_tbl[,c(2,3,4,6)]

colnames(lat_long_tbl) <- c("color_code","genome_rad","species","coords")


lat_long_tbl$color_code <- toupper(lat_long_tbl$color_code)

lat_long_tbl <- mutate(lat_long_tbl,lat="", long="")


for (line in 1:nrow(lat_long_tbl)) {
  lat_long_tbl$lat[line] <- unlist(strsplit(lat_long_tbl$coords[line],split = "," )  )[1]
  lat_long_tbl$long[line] <- unlist(strsplit(lat_long_tbl$coords[line],split = "," )  )[2]
}


world_map <- ggplot2::map_data("world")
mapa <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill="lightgray", colour = "white") +
  theme_bw()


# library(ggrepel)
#   # teste2 <-
#   teste +
#     geom_jitter(data = lat_long_tbl, mapping = aes(y = as.numeric(lat), x=as.numeric(long), group = color_code),
#                 shape = 21, size  = 4, color = "white" , fill = lat_long_tbl$color_code, inherit.aes = FALSE ) +
#    geom_label_repel(data = lat_long_tbl, mapping = aes(y = as.numeric(lat), x=as.numeric(long), label = genome_rad ),inherit.aes = FALSE)
#   #


#theme()
graph_name <- "geo_levs"
# teste2 <-
grafico_geo_loc <- mapa  +
  geom_label_repel(data = lat_long_tbl,
                   mapping = aes(y = as.numeric(lat), x=as.numeric(long), label = genome_rad ),
                   inherit.aes = FALSE, segment.alpha = 0.3) +
  geom_jitter(data = lat_long_tbl,
              mapping = aes(y = as.numeric(lat), x=as.numeric(long), group = color_code),
              shape = 21, size  = 4, color = "white" , fill = lat_long_tbl$color_code, alpha = 0.7,
              inherit.aes = FALSE, width = 2, height = 2)
  #add golubevii brasileira
  grafico_geo_loc +
    geom_label_repel(mapping = aes(y = -19.0061817, x=-57.660229, label = "gol'"), segment.alpha = 0.3)+
    geom_jitter(mapping = aes(y = -19.0061817, x=-57.660229),
              shape = 21, size  = 4, color = "white" , fill = "#2BCFDF", alpha = 0.7,
              width = 2, height = 2) +
  theme_void()

save_dir <- c("/home/heron-oh/metschnikowias/2020/figs")

ggsave(filename =  paste0(save_dir,"/",graph_name,".pdf"), plot = grafico_geo_loc, device = "pdf", width =16 , height = 9, dpi = 300)
ggsave(filename =  paste0(save_dir,"/",graph_name,".svg"), plot = grafico_geo_loc, device = "svg", width =16 , height = 9, dpi = 300)



