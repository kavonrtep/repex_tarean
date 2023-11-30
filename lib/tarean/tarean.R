#!/usr/bin/env Rscript
library(optparse, quiet = TRUE)
library(parallel)
if (interactive()){
  ## define functions only and exit
  ## assume that working directory was changes with source( chdir=TRUE)!!!
  script.dir=normalizePath('.')
  source('methods.R')
  source('logo_methods.R')
  source('htmlheader.R')
  options(OGDF = paste0(script.dir,"/OGDF/runOGDFlayout"))
  
}else{
  ## get options from command line
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dir <- normalizePath(dirname(script.name))
  oridir=getwd()
  ## parse arguments
  option_list = list(
    make_option(c('-i', '--input_sequences'),action='store',type='character',help='fasta file with input sequences',default=NA),
    make_option(c('-o', '--output_dir'),action='store',type='character',help='output directory',default="./kmer_analysis"),
    make_option(c('-m', '--min_kmer_length'),action='store',type='numeric',help='min kmer length',default=11),
    make_option(c('-x', '--max_kmer_length'),action='store',type='numeric',help='min kmer length',default=27),
    make_option(c('-n', '--cpu'),action='store',type='numeric',help='number of cpu to use',default=NULL),
    make_option(c('-s', '--sample_size'),action='store',type='numeric',help='number of sequences to use for analysis, is set to 0 all sequences are used',default=10000),
    make_option(c('-r', '--reorient_reads'),action='store_true',type='logical',help='number of cpu to use',default=FALSE),
    make_option(c('-l', '--no_layout'),action='store_true',type='logical',help='do not calculate graph layout',default=FALSE),
    make_option(c('-p', '--paired'),action='store_true',type='logical',help='reads are paired',default=FALSE),
    make_option(c('-t', '--tRNA_database='), action='store',type='character',help='path to tRNA database, is set PBS detection is performed',default=NULL)
    
  )

  description = paste (strwrap(" put decription here"), collapse ="\n")
  epilogue = paste (strwrap(" put epilogue here"), collapse ="\n")
  parser=OptionParser(
    option_list=option_list,
    epilogue=epilogue,
    description=description,
    )
  opt = parse_args(parser, args=commandArgs(TRUE))
  ## as Rscript
  options(OGDF = paste0(script.dir,"/OGDF/runOGDFlayout"))
  CPU = ifelse(is.null(opt$cpu), detectCores(), opt$cpu)
  source(paste(script.dir,"/","methods.R", sep=''))
  source(paste(script.dir,"/","logo_methods.R", sep=''))
  source(paste(script.dir,"/","htmlheader.R", sep=''))
  ## set number of CPU to use


  
  ## run tarean:
  tarean(
    opt$input_sequences,
    opt$output_dir,
    opt$min_kmer_length,
    opt$max_kmer_length,
    CPU,
    opt$sample_size,
    opt$reorient_reads,
    opt$tRNA_database,
    !opt$no_layout,
    paired = opt$paired
    )
}
