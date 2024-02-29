#!/usr/bin/env Rscript
# this script is used to convert data.tree object `to updated version
# compativle with data.tree v 1.1.0
library(data.tree)
library(ape)
setwd("databases")
rds_list <- c("classification_tree_viridiplantae_v3.0.rds",
              "classification_viridiplantae_tree.rds",
              "classification_tree_metazoa_v3.rds",
                "classification_tree_metazoa_v2.rds")

for (rds in rds_list) {
    print(rds)
    dt <- readRDS(rds)
    newick <- ToNewick(dt)
    tr <- read.tree(text=newick)
    dt2 <- as.Node(tr)
    out <- gsub("[.]rds$","_2.rds", rds)
    ab <- seq_along(dt2$leaves)
    for (i in seq_along(dt2$leaves)){
        print(i)
        print(dt2$leaves[[i]])
        try({
        dt2$leaves[[i]]$full_name <- dt$leaves[[i]]$full_name
        })
    }

    saveRDS(dt2, out)
}

 class(dt2$leaves[[i]]$full_name)