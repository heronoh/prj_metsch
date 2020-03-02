#extract GO anottations from interproscan results

library(stringr)
library(purrr)


interpro_dir <- c("/home/heron-oh/metschnikowias/2020/interproscan/v39/")
genome_rads <- readLines("/home/heron-oh/metschnikowias/2020/36_metsch_IDs.list")
out_dir <- c("/home/heron-oh/metschnikowias/2020/GO_interpro/")

getwd()



for (geno in genome_rads) {
grep_GO_cmd  <- paste0("grep 'GO' ",interpro_dir,geno,"/",geno,".tsv")
lines_grep_GO  <- system(grep_GO_cmd,intern = TRUE)
GOs_xtr <- str_extract_all(lines_grep_GO,pattern = "GO:[:digit:]++")

uniq_GOs <- unique(unlist(GOs_xtr))
assign(x = paste0(geno,"_GOs"),value = uniq_GOs)
# write.csv(uniq_GOs,file = paste0(out_dir,geno,"_GOs.txt"), col.names = FALSE,row.names = FALSE, quote = FALSE )
}

cat(setdiff(aus_GOs,bic_GOs))
cat(setdiff(bic_GOs,aus_GOs))
