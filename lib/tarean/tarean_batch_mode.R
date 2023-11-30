#!/usr/bin/env Rscript
library(optparse, quiet = TRUE)
library(parallel)
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name,"",
                   initial.options[grep(file.arg.name, initial.options)]
)
script.dir <- normalizePath(dirname(script.name))
oridir=getwd()
options(OGDF = paste0(script.dir,"/OGDF/runOGDFlayout2015.5"))
CPU =  detectCores()
source(paste(script.dir,"/","methods.R", sep=''))
source(paste(script.dir,"/","logo_methods.R", sep=''))
source(paste(script.dir,"/","htmlheader.R", sep=''))

option_list = list(
    make_option(c('-i', '--input_sequences_list'),
                action='store',type='character',
                help='list of fasta sequences file for tarean analysis'
                ),
    make_option(c('-o', '--output_dir'),
                action='store',type='character',
                help='output directory',
                default="./kmer_analysis"),
    make_option(c('-t', '--tRNA_database'),
                action='store',type='character',
                help='path to tRNA database',
                default=NULL),
    make_option(c('-p', '--parallel'),
                action='store_true',
                type='logical',
                help='run in parallel (faster but can exhaust RAM)',
                default=FALSE),
    make_option(c('-N', '--not_paired'),
                action='store_true',
                type='logical',
                help='reads are not paired',
                default=FALSE)

    )

description = paste (strwrap(" put decription here"), collapse ="\n")
epilogue = paste (strwrap(" put epilogue here"), collapse ="\n")
parser=OptionParser(
    option_list=option_list,
    epilogue=epilogue,
    description=description,
    )

opt = parse_args(parser, args=commandArgs(TRUE))
paired = !opt$not_paired
print(opt)
dir.create(opt$output_dir)
fl = readLines(opt$input_sequences_list)
## reorder to avoid running large top graphs at once
ord = sample(seq_along(fl), length(fl))


index=0
info=list()
save.image(paste0(opt$output_dir,"/info.RData")) # for debugin purposes
if (opt$parallel){
    cat("processing in parallel")
    info=mcmapply(
        FUN=tarean,
        input_sequences = fl[ord],
        output_dir = paste0(opt$output_dir,"/",sprintf("%04d",ord)),
        min_kmer_length = 11,
        max_kmer_length = 27,
        CPU = CPU,
        sample_size = 30000,
        reorient_reads = TRUE,
        tRNA_database_path = opt$tRNA_database,
        paired = paired,
        include_layout=FALSE,
        mc.cores=round(1+detectCores()/9),
        mc.set.seed = TRUE,
        mc.preschedule = FALSE,
        SIMPLIFY = FALSE
    )
}else{
    for (i in fl){
        index = index + 1
        dirout=paste0(opt$output_dir,"/",sprintf("%04d",index))
        try({
            info[[i]] = tarean(i, dirout, 11, 27, CPU, 30000, TRUE, opt$tRNA_database, include_layout=FALSE)
            cat("-----------------------------------------------------\n")
            print(info[[i]])
        })
    }
}
save(info, file = paste0(opt$output_dir,"/info.RData"))
save.image("tmp.RData")
## export as csv table
## 'graph_info' is always include:

tr_info = data.frame(do.call(rbind, info[sapply(info,length)>1]))
if (nrow(tr_info)>0){
    ## TR detected
    graph_info = data.frame (do.call(rbind, lapply(info, "[[", "graph_info")))
    graph_info$source=rownames(graph_info)
    tr_info$graph_info=NULL
    tr_info$source = rownames(tr_info)
    graph_tr_info = merge(graph_info, tr_info, all=TRUE, by='source')
    if (any(sapply(graph_tr_info,class)=='list')){
        for (i in colnames(graph_tr_info)){
            graph_tr_info[,i] = unname(unlist(graph_tr_info[,i]))
        }
    }
    write.table(graph_tr_info, file=paste0(opt$output_dir,"/info.csv"), row.names=FALSE,sep="\t", quote= TRUE)
}else{
    ## TR not detected
    graph_info = data.frame (do.call(rbind, lapply(info, function(x) unlist(x[['graph_info']]))))
    graph_info$source=rownames(graph_info)
    write.table(graph_info, file=paste0(opt$output_dir,"/info.csv"), row.names=FALSE,sep="\t", quote = FALSE)
}


