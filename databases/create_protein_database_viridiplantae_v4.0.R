#!/usr/bin/env Rscript
library(Biostrings)
domains <- readAAStringSet("/mnt/raid/454_data/databases/protein_domains/Viridiplantae_v4.0/Viridiplantae_v4.0_ALL_protein-domains.fasta")
element_names <- sapply(strsplit(names(domains), split="__"), function(x) paste(x[-1], collapse="__"))
classification <- readLines("/mnt/raid/454_data/databases/protein_domains/Viridiplantae_v4.0/Viridiplantae_v4.0_ALL_classification")
classification <- gsub("/", "_", classification)
names(classification) <- sapply(strsplit(classification, split="\t"), "[[", 1)
classification_formated <- sapply(sapply(strsplit(classification, "\t"), "[", -1), paste, collapse="/")
domain_name <- gsub("__.+", "", names(domains))
full_names <- paste0(names(domains), "#", classification_formated[element_names], ":", domain_name)
names(domains) <- full_names
writeXStringSet(domains, "/mnt/raid/users/petr/workspace/repex_tarean_github/databases/protein_database_viridiplantae_v4.0.fasta")

# later create diamond database and blast database - use singularity container to do that
# diamond v0.9.30.131 and  Package: blast 2.9.0, build Sep 30 2019 01:57:31

singularity_bin <- "/home/petr/data/miniforge3/envs/singularity/bin/singularity"
container_path <- "/mnt/raid/users/petr/workspace/repex_tarean_containers/repex_tarean_0.3.11-579a65d.sif"
database_dir <- "/mnt/raid/users/petr/workspace/repex_tarean_github/databases"
cmd <- paste(singularity_bin, "exec", "-B", database_dir,
             container_path, "diamond makedb --in",
             file.path(database_dir, "protein_database_viridiplantae_v4.0.fasta"),
             "--db", file.path(database_dir, "protein_database_viridiplantae_v4.0.fasta"))
system(cmd)

cmd <- paste(singularity_bin, "exec", "-B", database_dir,
             container_path, "makeblastdb -in",
             file.path(database_dir, "protein_database_viridiplantae_v4.0.fasta"),
             "-dbtype prot")
system(cmd)

library(data.tree)

add_weight <- function(i, name){
    if (is.null(i$parent[[name]])){
        i$parent[[name]] <- i[[name]]
    }else{
        i$parent[[name]] <- i[[name]] + i$parent[[name]]
    }
    if (i$parent$level == 1){
        return()
    }else{
        add_weight(i$parent, name)
    }
}


cls_string <- c(
    "All/contamination",
    "All/organelle/plastid",
    "All/organelle/mitochondria",
    "All/repeat/rDNA/45S_rDNA/18S_rDNA",
    "All/repeat/rDNA/45S_rDNA/25S_rDNA",
    "All/repeat/rDNA/45S_rDNA/5.8S_rDNA",
    "All/repeat/rDNA/5S_rDNA",
    "All/repeat/satellite",
    "All/repeat/mobile_element/Class_I/SINE",
    "All/repeat/mobile_element/Class_II/Subclass_1/TIR/MITE"
)
cls_full_name <- c(
    "contamination",
    "organelle/plastid",
    "organelle/mitochondria",
    "45S_rDNA/18S_rDNA",
    "45S_rDNA/25S_rDNA",
    "45S_rDNA/5.8S_rDNA",
    "5S_rDNA/5S_rDNA",
    "satellite",
    "Class_I/SINE",
    "Class_II/Subclass_1/TIR/MITE"
)

df1 <-  data.frame(pathString = cls_string, full_name=cls_full_name, stringsAsFactors = FALSE, nhits = 0,domains=0, prop=0, mean_weight=0, total_weight=0)
df2 <-  data.frame(pathString = paste("All/repeat/mobile_element",unique(classification_formated),sep="/"), full_name= unique(classification_formated), stringsAsFactors = FALSE,  nhits = 0, domains = 0, prop =0, mean_weight=0, total_weight = 0)
cls_tree <-  as.Node (rbind(df1,df2))
saveRDS(cls_tree, "/mnt/raid/users/petr/workspace/repex_tarean_github/databases/classification_tree_viridiplantae_v4.0.rds")