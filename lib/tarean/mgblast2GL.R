#!/usr/bin/env Rscript
## get script dir:
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- normalizePath(dirname(script.name))
fout <- commandArgs(T)[[2]]

source(paste(script.dir, "/methods.R", sep=''))
suppressPackageStartupMessages (library(igraph))
fin <- commandArgs(T)[[1]]

colcls = rep("NULL", 12)
colcls[c(1,5,11)] = c("character","character","numeric")
cat("loading mgblast table\n")
df = read.table(pipe(paste("cut -f1,5,11 ",fin)), sep="\t",comment.char="", as.is=TRUE, header= FALSE, colClasses = c("character","character","numeric"))

cat("creating graph\n")
GL = list()
colnames(df) =  c("V1", "V2", "weight")
GL$G = graph.data.frame(df , directed = FALSE)
print(summary(GL$G))
cat("calculating ogdf layouts\n")
try({
    L1 <- OGDFlayout(GL$G, alg=c("fmmm"))
})
cat("calculating fruchterman reingold layouts\n")
set.seed(vcount(GL$G))
L2 = layout.fruchterman.reingold(GL$G,dim=3)
if (class(L1) != "try-error"){
    GL$L <- cbind(L1[[1]],L2)
}else{
    GL$L <- L2
}
cat("saving output\n")
save(GL, file=fout)
