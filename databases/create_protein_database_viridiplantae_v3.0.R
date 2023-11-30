#!/usr/bin/env Rscript
## prepare protein database in format suitable fro repeatexplorer search
library(Biostrings)
domains = readAAStringSet("/mnt/raid/454_data/databases/protein_domains/new_protein_domains_prelim/coded01/ALL_protein-domains_05.fasta")
                                        # this cannot be used - __ is also in element id!!!
                                        # element_names = gsub("^.+__","",names(domains))  #

# this shou be version 3
domains = readAAStringSet("/mnt/raid/454_data/databases/protein_domains/Viridiplantae_v001/Viridiplantae_v001_ALL_protein-domains.fasta")
                                        # this cannot be used - __ is also in element id!!!
                                        # element_names = gsub("^.+__","",names(domains))  #

element_names = sapply(strsplit(names(domains),split="__"),function(x)paste(x[-1],collapse="__"))

classification = readLines("/mnt/raid/454_data/databases/protein_domains/Viridiplantae_v001/Viridiplantae_v001_ALL_classification")
## classification contain slash in categories - must be replaced with underscore
classification = gsub("/","_",classification)

names(classification) = sapply (strsplit(classification, split="\t"),"[[", 1)
classification_formated = sapply (sapply(strsplit(classification, "\t"), "[",-1), paste, collapse="/")
domain_name = gsub("__.+","",names(domains))
table(domain_name)
full_names = paste0(names(domains),"#", classification_formated[element_names],':', domain_name)
head(full_names)
names(domains) = full_names
writeXStringSet(domains,"/mnt/raid/users/petr/workspace/repex_tarean/databases/protein_database_viridiplantae_v3.0.fasta")


library(data.tree)
library(treemap)
data(GNI2014)
class(GNI2014)
head(classification_formated)

## compile all classification together
dna_dat = readDNAStringSet("/mnt/raid/users/petr/workspace/repex_tarean/databases/dna_database.fasta")
dna_dat = readDNAStringSet("/mnt/raid/users/petr/workspace/repex_tarean/databases/dna_database_masked.fasta")

add_weight = function(i, name){
    if (is.null(i$parent[[name]])){
        i$parent[[name]] = i[[name]]
    }else{
        i$parent[[name]] = i[[name]] + i$parent[[name]]
    }
    if (i$parent$level == 1){
        return()
    }else{
        add_weight(i$parent, name)
    }
}


cls_string = c(
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
cls_full_name =c(
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



df1 = data.frame(pathString = cls_string, full_name=cls_full_name, stringsAsFactors = FALSE, nhits = 0,domains=0, prop=0, mean_weight=0, total_weight=0)
df2 = data.frame(pathString = paste("All/repeat/mobile_element",unique(classification_formated),sep="/"), full_name= unique(classification_formated), stringsAsFactors = FALSE,  nhits = 0, domains = 0, prop =0, mean_weight=0, total_weight = 0)
cls_tree = as.Node (rbind(df1,df2))
saveRDS(object = cls_tree, file = "/mnt/raid/users/petr/workspace/repex_tarean/databases/classification_tree_viridiplantae_v3.0.rds")

print(cls_tree, "nhits", 'domains')

names(cls_tree)
cls_tree$leaves[[3]]$mean_weight

cls_tree$leaves[[1]]$mean_weight



## add to nodes
add_value_to_nodes = function(tr, name="weight"){
  for (i in Traverse(tr)){
      w = sum(sapply(i$leaves,"[[",name))
      i[[name]] = w
  }
  return(tr)
}


tr2 = add_value_to_nodes(cls_tree, name="nhits")
print(tr2,"nhits")

## add_nhits
for (i in cls_tree$leaves){
    cls_string = i$full_name
    ## go to the root:
    nhits = i$nhits
    prop = i$prop
    mean_weight = i$mean_weight
    add_weight(i, "nhits")
}




