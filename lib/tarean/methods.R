#!/user/bin/env Rscript

suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(hwriter))
suppressPackageStartupMessages(library(R2HTML))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))

max_ORF_length = function(s) {
  ## check all frames
  L = 0
  for (i in 1:3) {
    L = max(L, nchar(unlist(strsplit(as.character(translate(subseq(s, i))), "*", 
                                     fixed = TRUE))))
    L = max(L, nchar(unlist(strsplit(as.character(translate(subseq(reverseComplement(s), 
                                                                   i))), "*", fixed = TRUE))))
  }
  return(L)
}

kmers2graph = function(kmers, mode = "strong", prop = NULL) {
  kmerLength = nchar(kmers[1, 1])
  if (ncol(kmers) == 2) {
    kmers$size = kmers[, 2]/sum(kmers[, 2])
  }
  colnames(kmers) = c("name", "count", "size")
  if (!is.null(prop)) {  # tohle se nepouziva(prop je null), a je to asi spatne - filtuje se to pred tridenim!!
    p = cumsum(kmers$size)
    kmers = kmers[p < prop, ]
  }
  N = dim(kmers)[1]
  kmers = kmers[order(kmers$size), ]
  ## convert kmers to fasta file
  kms = data.frame(kmer = substring(kmers$name, 1, kmerLength - 1), ids = 1:nrow(kmers),stringsAsFactors = FALSE)
  kme = data.frame(kmer = substring(kmers$name, 2), ide = 1:nrow(kmers), stringsAsFactors = FALSE)

  ## df = merge(kms,kme, by = 'kmer',all=FALSE)[,c(2,3,1)]
  df = inner_join(kme,kms, by = 'kmer')[,c(2,3)]

  ## names(kms) = seq_along(kms)
  ## kme = substring(kmers$name, 2)
  ## names(kme) = seq_along(kme)
  ## ## use new blast!
  ## database = tempfile()
  ## query = tempfile()
  ## output = tempfile()
  ## writeXStringSet(DNAStringSet(kms), filepath = database, format = "fasta")
  ## writeXStringSet(DNAStringSet(kme), filepath = query, format = "fasta")
  ## cmd = paste("makeblastdb -in", database, "-dbtype nucl")
  ## system(cmd, ignore.stdout = TRUE)
  ## cmd = paste("blastn -outfmt '6 qseqid sseqid pident'  -strand plus -dust no -perc_identity 100 -query ", 
  ##     query, "-db", database, "-word_size", kmerLength - 1, "-out", output)
  ## system(cmd)
  ## df = try({
  ##   read.table(output, as.is = TRUE)
  ## })
  ## if (class(df) == "try-error"){
  ##   print("creation of kmer graph failed")
  ##   print(query)
  ##   print(output)
  ##   print(database)
  ##   return(NULL)
  ## }
  ## unlink(query)
  ## unlink(paste(database, "*", sep = ""))
  ## unlist(output)
  gm_mean = function(x, na.rm = TRUE) {
    exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
  }
  
  whg = apply(cbind(kmers[df[, 1], 2], V2 = kmers[df[, 2], 2]), 1, gm_mean)
  G = graph.data.frame(data.frame(V1 = kmers$name[df[, 1]], V2 = kmers$name[df[, 
                                                                               2]], weight = whg), vertices = kmers[, 1:3])
                                        # separate to connected components:
  ccs = clusters(G, mode = mode)$membership
  sel_cls = which(tabulate(ccs) > 1)
  Gs = list()
  for (i in seq_along(sel_cls)) {
    Gs[[i]] = induced.subgraph(G, vids = which(ccs %in% sel_cls[i]))
  }
  ## reorder!!!
  Gs = Gs[order(sapply(Gs, vcount), decreasing = TRUE)]
  return(Gs)
}


OGDFlayout = function(G, ncol = NULL, alg = "fmmm", OGDF = getOption("OGDF")) {
  ## is ogdf binary available?
  if (is.null(OGDF)) {
    OGDF = Sys.getenv("OGDF")
    if ("" == OGDF) {
      options(warn = -1)
      OGDF = system("which runOGDFlayout", intern = TRUE)
      options(warn = 0)
      if (length(OGDF) == 0) {
        cat("path to runOGDFlayout not found\n")
        return(NULL)
      }
      
    }
  }
  if (is.null(ncol)) {
    if (is.null(E(G)$weight)) {
      el = data.frame(get.edgelist(G, names = TRUE), rep(1, ecount(G)))
    } else {
      el = data.frame(get.edgelist(G, names = TRUE), E(G)$weight)
    }
    ncol = paste(tempfile(pattern = as.character(Sys.getpid())), ".layout", sep = "")
    write.table(el, file = ncol, row.names = FALSE, col.names = FALSE, sep = "\t", 
                quote = FALSE)
  } else {
                                        # copy ncol:
    ncol_tmp = paste(tempfile(pattern = as.character(Sys.getpid())), ".layout", 
                     sep = "")
    file.copy(ncol, ncol_tmp)
    ncol = ncol_tmp
  }
  algopt = c("fmmm", "sm", "fme")
  if (!(alg %in% c("fmmm", "sm", "fme") && TRUE)) {
    stop("alg must by :", algopt, "\n")
    
  }
  
                                        # output file:
  Louts = list()
  layout_file = tempfile(pattern = as.character(Sys.getpid()))
  for (i in alg) {
    cmd = paste(OGDF, "-alg", i, "-iff layout -off layout -out", layout_file, 
                ncol)
    system(cmd, intern = TRUE)
    L = read.table(layout_file, skip = ecount(G))
    L = L[match(V(G)$name, L$V2), ]
    Lout = as.matrix(L[, -(1:2)])
    unlink(layout_file)
    Louts[[i]] = Lout
    
  }
                                        # clean up
  unlink(ncol)
  return(Louts)
}

xcolor_code = c(A = "#FF0000", C = "#00FF00", T = "#0000FF", G = "#FFFF00")

kmers2color = function(s, position = NULL) {
  if (is.null(position)) {
    position = round(nchar(s[1])/2, 0)
    ## position = 1
    position = nchar(s[1])
  }
  color_code = c(A = "#FF0000", C = "#00FF00", T = "#0000FF", G = "#FFFF00")
  color_base = substring(s, position, position)
  colors = color_code[color_base]
  names(colors) = color_base
  return(colors)
}
get_sequence = function(g, v, position = NULL) {
  s = V(g)$name
  if (is.null(position)) {
    position = round(nchar(s[1])/2, 0)
    ## position = 1
    position = nchar(s[1])
  }
  nt = paste(substring(s[v], position, position), collapse = "")
  return(nt)
}


get_mimimal_cc = function(km, thr = 20, min_coverage = 0.45, step = 2, start = NULL) {
  if (is.null(start)) {
    i = sum(cumsum(km$freq) < 0.5)
  } else {
    i = sum(cumsum(km$freq) < start)
  }
  continue = TRUE
  while (continue) {
    if (i > nrow(km)) {
      i = nrow(km)
      continue = FALSE
      step = 1
    }
    GG = kmers2graph(km[1:i, ])
    if (length(GG) > 0) {
      if (vcount(GG[[1]]) > thr) {
        if (sum(V(GG[[1]])$size) >= min_coverage) {
          GG[[1]]$input_coverage = sum(km$freq[1:i])
          GG[[1]]$L = OGDFlayout(GG[[1]])[[1]]
          return(GG[[1]])
        }
      }
    }
    i = round(i * step)
    
  }
  if (length(GG) == 0 | is.null(GG)) {
    return(NULL)
  }
  
  GG[[1]]$input_coverage = sum(km$freq[1:i])
  GG[[1]]$L = OGDFlayout(GG[[1]])[[1]]
  return(GG[[1]])
}



paths2string = function(paths) {
  pathstring = sapply(lapply(lapply(paths, as_ids), substring, 1, 1), paste, collapse = "")
  return(pathstring)
}


align_paths = function(paths, G) {
  shift = rep(NA, length(paths))
  thr = 0  # minimal length
  tr_paths = list()
  Seqs = list()
  centre_node = as.numeric(names(sort(table(unlist(paths)), decreasing = TRUE)))[[1]]
  
  for (i in seq_along(paths)) {
    if (centre_node %in% paths[[i]]) {
      S = which(paths[[i]] %in% centre_node)
      shift[i] = S
      if (S == 1) {
        tr_paths[[i]] = paths[[i]]
      } else {
        tr_paths[[i]] = c(paths[[i]][S:length(paths[[i]])], paths[[i]][1:(S - 
                                                                          1)])
      }
      Seqs[[i]] = get_sequence(G, tr_paths[[i]])
    } else {
      shift[i] = NA
    }
  }
  paths_n = lapply(paths, as.numeric)
  tr_paths_n = do.call(cbind, lapply(tr_paths, as.numeric))
  new_shift = shift
  for (i in which(is.na(shift))) {
    score = numeric(length(paths_n))
    for (S in seq_along(paths_n[[i]])) {
      if (S == 1) {
        path_tmp_n = paths_n[[i]]
      } else {
        path_tmp_n = c(paths_n[[i]][S:length(paths_n[[i]])], paths_n[[i]][1:(S - 
                                                                             1)])
      }
      score[S] = sum(tr_paths_n == path_tmp_n)
    }
    if (sum(score) != 0) {
      S = which.max(score)
      new_shift[i] = S
      if (S == 1) {
        tr_paths[[i]] = paths[[i]]
      } else {
        tr_paths[[i]] = c(paths[[i]][S:length(paths[[i]])], paths[[i]][1:(S - 
                                                                          1)])
      }
      Seqs[[i]] = get_sequence(G, tr_paths[[i]])
    }
  }
  shift = new_shift
                                        # try to shift based on the sequences itself
  return(list(Seqs = Seqs[!is.na(shift)], tr_paths = tr_paths[!is.na(shift)], shift = shift[!is.na(shift)]))
}

make_consensus = function(paths_info, G) {
  include = !is.na(paths_info$shift)
  ## make alignments using mafft
  aln = mafft(unlist(paths_info$Seqs[include]))
  CM = calculate_consensus_matrix(aln = aln, tr_paths = paths_info$tr_paths[include], 
                                  G = G)
  gaps = get_gaps_from_alignment(aln)
  CMnorm = CM/rowSums(CM)
  bases = colnames(CM)
  consensus = sapply(apply(CMnorm, 1, function(x) bases[which(x > 0.2)][order(x[x > 
                                                                                0.2], decreasing = TRUE)]), paste, collapse = "")
  consensus2 = gsub("-", "", paste0(substring(consensus, 1, 1), collapse = ""), 
                    fixed = TRUE)
  number_of_SNP = sum(rowSums(CM > 0) > 1)
  SNP_positions = which(rowSums(CM > 0) > 1)
  number_of_position_with_indels = sum(colSums(do.call(rbind, strsplit(as.character(aln), 
                                                                       "")) == "-") > 0)
  indel_positions = which(colSums(do.call(rbind, strsplit(as.character(aln), "")) == 
                                  "-") > 0)
  if (length(SNP_positions) > 0) {
    variable_sites = unlist(c(c(mapply(paste, strsplit((sapply(apply(CMnorm, 
                                                                     1, function(x) bases[which(x > 0.2)]), paste, collapse = ""))[SNP_positions], 
                                                       ""), SNP_positions, sep = "_")), paste("-", indel_positions, sep = "_")))
  } else {
    variable_sites = NULL
  }
  variable_positions = unique(SNP_positions, indel_positions)
  return(list(aln = aln, CM = CM, CMnorm = CMnorm, consensus = consensus, consensus2 = consensus2, 
              number_of_SNP = number_of_SNP, SNP_positions = SNP_positions, number_of_position_with_indels = number_of_position_with_indels, 
              indel_positions = indel_positions, variable_positions = variable_positions, 
              variable_sites = variable_sites, gaps = gaps))
}


estimate_monomer = function(G, weights = NULL, limit = NULL) {
  if (is.null(G)) {
    return(NULL)
  }
  ## estimate monomer from kmer based graph
  V(G)$id = 1:vcount(G)
  GS = induced_subgraph(G, vids = which(degree(G) == 2))  ## keep only vertices without branching
  cls = clusters(GS)$membership
  
  
  ids = mapply(FUN = function(x, y) x[which.max(y)], split(V(GS)$id, cls), split(V(GS)$size, 
                                                                                 cls))  ## from each branch use only one vertex with larges size!
  
  
  ids = ids[order(V(G)$size[ids], decreasing = TRUE)]
  ids_size = V(G)$size[ids]
  N50 = sum(cumsum(ids_size)/sum(ids_size) < 0.5)
  if (length(ids) > 10000) {
    ids = ids[1:N50]
  }
  ## use only large vertices in search!  how many?
  el = get.edgelist(G, names = FALSE)
  node_use = numeric(vcount(G))
  LL = numeric()
  ## W=numeric()
  i = 0
  paths = list()
  if (is.null(weights)) {
    weights = (max(E(G)$weight) - (E(G)$weight) + 1)
    weights = E(G)$weight^(-3)
  }
  included = rep(FALSE, vcount(G))
  W_total = sum(V(G)$size)
  
  coverage = numeric()
  t0 = c()
  i = 0
  j = 0
  for (n in ids) {
    j = j + 1
    t0[j] = Sys.time()
    if (included[n]) {
      next
    }
    m = which(el[, 1] %in% n)
    i = i + 1
    s = get.shortest.paths(G, el[m, 2], el[m, 1], weights = weights, output = "vpath")
    included[as.numeric(s$vpath[[1]])] = TRUE
    paths[[i]] = s$vpath[[1]]
    LL[i] = (length(s$vpath[[1]]))
  }
  
  ## evaluate if paths should be divided to variants - by length and path weight
  paths_clusters = split(paths, LL)
  paths_clusters_tr = mclapply(paths_clusters, FUN = align_paths, G = G, mc.cores = getOption("CPU"))
  
  ## paths_clusters_tr = lapply(paths_clusters, FUN = align_paths, G = G) consensus
  paths_consensus = mclapply(paths_clusters_tr, make_consensus, G = G, mc.cores = getOption("CPU"))
  
  ## evaluate weight for individual paths:
  for (v in seq_along(paths_consensus)) {
    p = paths_clusters_tr[[v]]$tr_paths
    ## clean
    p = p[!sapply(p, function(x) anyNA(x) | is.null(x))]
    L = sapply(p, length)
    p_groups = split(p, L)
    w_groups = sapply(p_groups, function(x) sum(V(G)$size[unique(c(sapply(x, 
                                                                          as.numeric)))]))
    total_score = sum(V(G)$size[unique(c(unlist(sapply(p, as.numeric))))])
    LW = data.frame(`Length estimate` = unique(L), weight = w_groups, stringsAsFactors = FALSE)
    LW = LW[order(LW$weight, decreasing = TRUE), ]
    rownames(LW) = NULL
    paths_consensus[[v]]$total_score = total_score
    paths_consensus[[v]]$length_variant_score = LW
    
  }
  return(list(estimates = paths_consensus, paths = paths_clusters_tr))
}



detect_pbs = function(dimers_file, tRNA_database_path, reads_file, output) {
  ## read_file contain oriented reads!
  min_offset = 10
  max_end_dist = 2
  min_aln_length = 30
  max_pbs_search = 30
  thr_length = 12
  end_position = 23
  insertion_proportion_threshold <- 0.15
  ## read_file contain oriented reads!  for testing
  if (FALSE) {
    library(Biostrings)
    thr_length = 12
    end_position = 23
    dimers_file = "consensus_dimer.fasta"
    reads_file = "reads_oriented.fas"
    tRNA_database_path = "/mnt/raid/users/petr/workspace/repex_tarean/databases/tRNA_database2.fasta"
    
  }
  ## find read which are half aligned do reference dimer
  dimers_seq = readDNAStringSet(dimers_file)
  cmd = paste("makeblastdb -in", dimers_file, "-dbtype nucl")
  system(cmd, ignore.stdout = TRUE)
  output = paste0(reads_file, "blast_out.cvs")
  columns = c("qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "sstrand", 
              "slen", "pident", "length")
  cmd = paste("blastn -outfmt '6", paste(columns, collapse = " "), "' -perc_identity 90 -query ", 
              reads_file, "-db", dimers_file, "-word_size", 7, "-dust no -num_alignments 99999 -strand plus ", 
              "-out", output)
  system(cmd)
  blastdf = read.table(output, as.is = TRUE, col.names = columns, comment.char = "")
  unlink(output)
  blastdf = blastdf[blastdf$length >= min_aln_length, ]
  
  ## expand two whole read to see unligned reads parts
  blastdf_expand = blastdf
  blastdf_expand$qstart = blastdf$qstart - (blastdf$qstart - 1)
  blastdf_expand$sstart = blastdf$sstart - (blastdf$qstart - 1)
  blastdf_expand$qend = blastdf$qend + (blastdf$qlen - blastdf$qend)
  blastdf_expand$send = blastdf$send + (blastdf$qlen - blastdf$qend)
  pS = blastdf$sstart
  pE = blastdf$send
  pSF = blastdf_expand$sstart
  pEF = blastdf_expand$send
  cond1 = pS - pSF >= min_offset & blastdf$qend >= (blastdf$qlen - max_end_dist) & 
    blastdf$sstart >= max_pbs_search
  cond2 = pEF - pE >= min_offset & blastdf$qstart <= max_end_dist & blastdf$send <= 
    blastdf$slen - max_pbs_search
  
  ## coverage of alignments: evaluate coverage at site with breaks, it is neccessary
  ## that there must be at least 20% of interupted alignments at given position to
  ## be considered for additional search
  coverage_profiles = subject_coverage(blastdf)
  
  
  ## extract flanking sequences - cca 50nt, it should contain tRNA sequence search
  ## Left
  fin = tempfile()
  fout = tempfile()
  scoreL = scoreR = 0
  
  
  ## left side unique insertion sites
  if (any(cond1)) {
    insertion_sites = ddply(blastdf[cond1, ], .(sseqid, sstart), nrow)
    ## check coverage
    insertion_sites$coverage = 0
    for (i in 1:nrow(insertion_sites)) {
      insertion_sites$coverage[i] = coverage_profiles[[insertion_sites$sseqid[i]]][insertion_sites$sstart[[i]]]
    }
    insertion_OK_left <- insertion_sites[with(insertion_sites, V1/coverage > 
                                                               insertion_proportion_threshold), ]
    
    if (nrow(insertion_OK_left) > 0) {
      s = ifelse(insertion_OK_left[, "sstart"] - max_pbs_search < 1, 1, insertion_OK_left[, 
                                                                                          "sstart"] - max_pbs_search)
      Lpart = subseq(dimers_seq[match(insertion_OK_left$sseqid, names(dimers_seq))], 
                     s, insertion_OK_left$sstart)
      
      ## names are CONSENSUSID__READID_position
      names(Lpart) = paste0(names(Lpart), "__", insertion_OK_left$V1, "__", 
                            insertion_OK_left$sstart)
      ## check presence TG dinucleotide
      TG = vcountPattern("TG", subseq(dimers_seq[match(insertion_OK_left$sseqid, 
                                                       names(dimers_seq))], insertion_OK_left$sstart - 2, insertion_OK_left$sstart + 
                                                                                                          2))
      if (any(TG > 0)) {
        writeXStringSet(Lpart[TG > 0], filepath = fin)
        cmd = paste("blastn -outfmt '6", paste(columns, collapse = " "), 
                    "' -perc_identity 83 -query ", fin, "-db", tRNA_database_path, 
                    "-word_size", 7, "-strand plus -dust no -max_target_seqs 10000 -evalue 100", 
                    "-out", fout)
        
        system(cmd, ignore.stdout = TRUE)
        df = read.table(fout, as.is = TRUE, col.names = columns, comment.char = "")
        filter1 = df$length >= thr_length
        filter2 = ifelse(df$send > df$sstart, df$send, df$sstart) >= end_position
        df_pass = df[filter1 & filter2, , drop = FALSE]
        df_pass_L = df_pass[!duplicated(df_pass$qseqid), , drop = FALSE]
        scoreL = get_score(df_pass_L)
        write.table(df_pass_L, file = paste0(output, "_L.csv"), row.names = FALSE, 
                    sep = "\t")
      }
    }
  }
  if (any(cond2)) {
    ## search Right
    insertion_sites = ddply(blastdf[cond2, ], .(sseqid, send, slen), nrow)
    ## check coverage
    insertion_sites$coverage = 0
    for (i in 1:nrow(insertion_sites)) {
      insertion_sites$coverage[i] = coverage_profiles[[insertion_sites$sseqid[i]]][insertion_sites$send[[i]]]
    }
    insertion_OK_right <- insertion_sites[with(insertion_sites, V1/coverage > 
                                                                insertion_proportion_threshold), ]
    
    if (nrow(insertion_OK_right) > 0) {
      s = ifelse(insertion_OK_right$send + max_pbs_search > insertion_OK_right$slen, 
                 insertion_OK_right$slen, insertion_OK_right$send + max_pbs_search)
      Rpart = subseq(dimers_seq[match(insertion_OK_right$sseqid, names(dimers_seq))], 
                     insertion_OK_right$send, s)
      names(Rpart) = paste0(names(Rpart), "__", insertion_OK_right$V1, "__", 
                            insertion_OK_right$send)
      
      ## check presence CA dinucleotide
      CA = vcountPattern("CA", subseq(dimers_seq[match(insertion_OK_right$sseqid, 
                                                       names(dimers_seq))], insertion_OK_right$send - 2, insertion_OK_right$send + 
                                                                                                         2, ))
      if (any(CA > 0)) {
        writeXStringSet(Rpart[CA > 0], filepath = fin)
        cmd = paste("blastn -outfmt '6", paste(columns, collapse = " "), 
                    "' -perc_identity 83 -query ", fin, "-db", tRNA_database_path, 
                    "-word_size", 7, "-strand minus -dust no -max_target_seqs 10000 -evalue 100", 
                    "-out", fout)
        
        system(cmd, ignore.stdout = TRUE)
        df = read.table(fout, as.is = TRUE, col.names = columns, comment.char = "")
        filter1 = df$length >= thr_length
        filter2 = ifelse(df$send > df$sstart, df$send, df$sstart) >= end_position
        df_pass = df[filter1 & filter2, , drop = FALSE]
        df_pass_R = df_pass[!duplicated(df_pass$qseqid), , drop = FALSE]
        write.table(df_pass_R, file = paste0(output, "_R.csv"), row.names = FALSE, 
                    sep = "\t")
        scoreR = get_score(df_pass_R)
      }
    }
  }
  unlink(fin)
  unlink(fout)
  return(max(scoreL, scoreR))
}


subject_coverage = function(blastdf) {
  ## calculate coverage for all blast subjects
  coverage_profiles <- by(blastdf, INDICES = blastdf$sseqid, FUN = function(x) {
    as.numeric(coverage(IRanges(start = x$sstart, end = x$send)))
  })
  return(coverage_profiles)
  
}

get_score = function(x) {
  ## keep best tRNA
  if (nrow(x) == 0) {
    return(0)
  }
  xm = data.frame(do.call(rbind, strsplit(x$qseqid, "__")), stringsAsFactors = FALSE)
  xm$score = as.numeric(sapply(strsplit(xm[, 1], "_"), "[[", 4))
  xm$AA = gsub("-.+$", "", gsub("^.+__", "", x$sseqid))
  best_score = max(by(xm$score, INDICES = xm$AA, FUN = sum))
  return(best_score)
}



dotter = function(seq1, seq2 = NULL, params = "") {
  if (is.null(seq2)) {
    seq2 = seq1
  }
  library(Biostrings)
  if (class(seq1) != "DNAStringSet") {
    seq1 = BStringSet(seq1)
  }
  if (class(seq2) != "DNAStringSet") {
    seq2 = BStringSet(seq2)
  }
  sf1 = tempfile("seq1")
  writeXStringSet(seq1, file = sf1)
  sf2 = tempfile("seq2")
  writeXStringSet(seq2, file = sf2)
  system(paste("dotter", params, sf1, sf2), wait = FALSE)
  Sys.sleep(2)
  unlink(c(sf1, sf2))
  return(NULL)
}



dotter2 = function(seq1, seq2 = NULL, params = NULL) {
  if (is.null(seq2)) {
    seq2 = seq1
  }
  if (is.null(params)) {
    params = " -windowsize 30 -threshold 45 "
  }
  library(Biostrings)
  if (class(seq1) != "DNAStringSet") {
    seq1 = DNAStringSet(seq1)
  }
  if (class(seq2) != "DNAStringSet") {
    seq2 = DNAStringSet(seq2)
  }
  L1 = nchar(seq1)
  L2 = nchar(seq2)
  
  tmpdat1 = tempfile()
  
  dir.create(tmpdat1)
  tmpdat2 = tempfile()
  dir.create(tmpdat2)
  oridir = getwd()
  
  
  seq1_merged = DNAStringSet(paste(seq1, collapse = ""))
  seq2_merged = DNAStringSet(paste(seq2, collapse = ""))
  seq2rc_merged = reverseComplement(DNAStringSet(paste(seq2, collapse = "")))
  
  
  sf1 = tempfile("seq1")
  writeXStringSet(seq1_merged, filepath = sf1)
  sf2 = tempfile("seq2")
  sf2rc = tempfile("seq2rc")
  writeXStringSet(seq2_merged, filepath = sf2)
  writeXStringSet(seq2rc_merged, filepath = sf2rc)
  
  cmd1 = paste("dotmatcher -graph data -asequence ", sf1, "-bsequence", sf2, params)
  cmd2 = paste("dotmatcher -graph data -asequence ", sf1, "-bsequence", sf2rc, 
               params)
  setwd(tmpdat1)
  output1 = system(cmd1, intern = TRUE)
  setwd(tmpdat2)
  output2 = system(cmd2, intern = TRUE)
  setwd(oridir)
  
  
  
  fout1 = strsplit(tail(output1, n = 1), split = " ")[[1]][2]
  rawdat1 = readLines(paste(tmpdat1, "/", fout1, sep = ""))
  rawdat1 = rawdat1[1:min(grep("^Rectangle", rawdat1))]
  
  if (length(rawdat1[grep("^Line", rawdat1)]) == 0) {
    coord1 = NULL
  } else {
    coord1 = apply(sapply(strsplit(rawdat1[grep("^Line", rawdat1)], " "), "[", 
                          c(3, 5, 7, 9, 11)), 1, as.numeric)
    coord1 = matrix(coord1, ncol = 5)
  }
  
  fout2 = strsplit(tail(output2, n = 1), split = " ")[[1]][2]
  rawdat2 = readLines(paste(tmpdat2, "/", fout2, sep = ""))
  rawdat2 = rawdat2[1:min(grep("^Rectangle", rawdat2))]
  
  if (length(rawdat2[grep("^Line", rawdat2)]) == 0) {
    coord2 = NULL
  } else {
    coord2 = apply(sapply(strsplit(rawdat2[grep("^Line", rawdat2)], " "), "[", 
                          c(3, 5, 7, 9, 11)), 1, as.numeric)
    coord2 = matrix(coord2, ncol = 5)
  }
  unlink(sf1)
  unlink(sf2)
  unlink(sf2rc)
  unlink(tmpdat1, recursive = TRUE)
  unlink(tmpdat2, recursive = TRUE)
  
  N1 = sum(nchar(seq1))
  N2 = sum(nchar(seq2))
  op = par(xaxs = "i", yaxs = "i", mar = c(5, 2, 6, 10), las = 1)
  on.exit(par(op))
  plot(c(1, N1), c(1, N2), type = "n", xlab = "", ylab = "", axes = FALSE)
  if (!is.null(coord1)) {
    segments(x0 = coord1[, 1], y0 = coord1[, 2], x1 = coord1[, 3], y1 = coord1[, 
                                                                               4])
  }
  if (!is.null(coord2)) {
    segments(x0 = coord2[, 1], y0 = N2 - coord2[, 2], x1 = coord2[, 3], y1 = N2 - 
                                                                          coord2[, 4])
  }
  abline(v = c(0, cumsum(L1)), col = "green")
  abline(h = c(0, cumsum(L2)), col = "green")
  box()
  axis(1, at = c(1, cumsum(L1))[-(length(L1) + 1)], labels = names(seq1), hadj = 0, 
       cex.axis = 0.7)
  axis(4, at = c(1, cumsum(L2))[-(length(L2) + 1)], labels = names(seq2), hadj = 0, 
       cex.axis = 0.7)
  invisible(list(coord1, coord2))
}


CM2sequence = function(x) {
  bases = c(colnames(x), "N")
  x = cbind(x, 0)
  colnames(x) = bases
  seqs = paste(bases[apply(x, 1, which.max)], collapse = "")
  
  return(seqs)
}

## in alingment calculate significance from kmers
calculate_consensus_matrix = function(aln, tr_paths, G) {
  bases = c("A", "C", "G", "T", "-")
  positions = lapply(strsplit(as.character(aln), split = ""), function(x) which(!x %in% 
                                                                                "-"))
  base_matrix = do.call(rbind, strsplit(as.character(aln), split = ""))
  weights = matrix(0, nrow = length(aln), ncol = nchar(aln))
  kmer_rel_proportions = (V(G)$size/table(factor(unlist(tr_paths), levels = 1:vcount(G))))
  for (i in seq_along(positions)) {
    weights[i, positions[[i]]] = kmer_rel_proportions[tr_paths[[i]]]
  }
  names(kmer_rel_proportions) = V(G)$name
  ## get weights for gaps by approximation
  fitgaps = function(y) {
    if (sum(y == 0) == 0) {
      return(y)
    } else {
      y0 = rep(y[y != 0], 3)
      x0 = which(rep(y, 3) != 0)
      fitted = approx(x0, y0, xout = seq_along(rep(y, 3)), rule = 1)
    }
    return(fitted$y[-seq_along(y)][seq_along(y)])
  }
  weights_with_gaps = t(apply(weights, 1, fitgaps))
  ## weights_with_gaps = get_gap_weights(weights, aln,G) ## this step take wey too
  ## long time!!!-not so important TODO handle gaps differently - more effectively
  
  ## get consensus matrix
  CM = sapply(1:nchar(aln[[1]]), function(i) {
    sapply(split(weights_with_gaps[, i], factor(base_matrix[, i], levels = bases)), 
           sum)
  })
  return(t(CM)[, 1:5])
}


get_gaps_from_alignment = function(aln) {
  as.character(aln)
  gaps_positions = unique(do.call(rbind, str_locate_all(as.character(aln), "-+")))
  return(gaps_positions)
}

plot_kmer_graph = function(G, L = NULL, vertex.size = NULL, ord = NULL, upto = NULL, 
                           highlight = NULL) {
  if (!is.null(G$L) & is.null(L)) {
    L = G$L
  }
  if (is.null(L)) {
    ## L=layout.kamada.kawai(G)
    L = OGDFlayout(G)[[1]]
  }
  clr = kmers2color(V(G)$name)
  if (!is.null(highlight)) {
    clr[highlight] = "#00000080"
  }
  if (!is.null(ord)) {
    clr[ord[1:upto]] = paste(clr[ord[1:upto]], "30", sep = "")
  }
  if (is.null(vertex.size)) {
    vertex.size = rescale(V(G)$size, to = c(0.5, 6))
    
  }
  plot(G, layout = L, vertex.label = "", vertex.size = vertex.size, edge.curved = FALSE, 
       vertex.color = clr, vertex.frame.color = "#00000020", edge.width = 2, edge.arrow.mode = 1, 
       edge.arrow.size = 0.2)
  
}

rglplot_kmer_graph = function(G, L = NULL, vertex.size = 4) {
  if (is.null(L)) {
    set.seed(vcount(G))
    L = layout.kamada.kawai(G)
  }
  
  rglplot(G, layout = L, vertex.label = "", vertex.size = vertex.size, edge.curved = FALSE, 
          vertex.color = kmers2color(V(G)$name), edge.width = 2, edge.arrow.mode = 1, 
          edge.arrow.size = 0.5)
}




mafft = function(seqs, params = "--auto --thread 1 ") {
  if (length(seqs) < 2) {
    return(seqs)
  }
  infile = tempfile()
  if (class(seqs) == "character") {
    seqs = DNAStringSet(seqs)
  }
  writeXStringSet(seqs, file = infile)
  outfile = tempfile()
  cmd = paste("mafft --quiet --nuc ", params, infile, "2> /dev/null > ", outfile)
  system(cmd, intern = TRUE, ignore.stderr = FALSE)
  aln = readDNAStringSet(outfile)
  unlink(c(outfile, infile))
  return(aln)
}

mgblast = function(databaseSeq, querySeq) {
  params = " -p 85 -W18 -UT -X40 -KT -JF -F \"m D\" -v100000000 -b100000000 -D4 -C 30 -D 30 "
                                        # no dust filtering:
  paramsDF = " -p 85 -W18 -UT -X40 -KT -JF -F \"m D\" -v100000000 -b100000000 -D4 -F F"
  
  database = tempfile()
  query = tempfile()
  output = tempfile()
  if (class(databaseSeq) == "character") {
    database = databaseSeq
    do_not_delete_database = TRUE
  } else {
    writeXStringSet(databaseSeq, filepath = database, format = "fasta")
    do_not_delete_database = FALSE
  }
  if (class(querySeq) == "character") {
    query = querySeq
    do_not_delete_query = TRUE
  } else {
    writeXStringSet(querySeq, filepath = query, format = "fasta")
    do_not_delete_query = FALSE
  }
  ## create database:
  cmd = paste("formatdb -i", database, "-p F")
  system(cmd)
  cmd = paste("mgblast", "-d", database, "-i", query, "-o", output, params)
  system(cmd)
  if (file.info(output)$size == 0) {
                                        # no hist, try wthou dust masker
    cmd = paste("mgblast", "-d", database, "-i", query, "-o", output, paramsDF)
    system(cmd)
  }
  blastOut = read.table(output, sep = "\t", header = FALSE, as.is = TRUE, comment.char = "")
  unlink(output)
  if (!do_not_delete_query) {
    unlink(query)
  }
  if (!do_not_delete_database) {
    unlink(paste(database, "*", sep = ""))
  }
  colnames(blastOut) = c("query", "q.length", "q.start", "q.end", "subject", "s.length", 
                         "s.start", "s.end", "pid", "weight", "e.value", "strand")
  blastOut
}


estimate_sample_size = function(NV, NE, maxv, maxe) {
  ## density
  d = (2 * NE)/(NV * (NV - 1))
  eEst = (maxv * (maxv - 1) * d)/2
  nEst = (d + sqrt(d^2 + 8 * d * maxe))/(2 * d)
  if (eEst >= maxe) {
    N = round(nEst)
    E = round((N * (N - 1) * d)/2)
    
  }
  if (nEst >= maxv) {
    N = maxv
    E = round((N * (N - 1) * d)/2)
    
  }
  return(N)
}




mgblast2graph = function(blastfile, seqfile, graph_destination, directed_graph_destination, 
                         oriented_sequences, paired = TRUE, repex = FALSE, image_file = NULL, image_file_tmb = NULL, 
                         include_layout = TRUE, pair_completeness = NULL, satellite_model_path = NULL, 
                         maxv = 40000, maxe = 5e+08, seqfile_full = seqfile) {
  cat("loading blast results\n")
  if (repex) {
    cln = c("query", "subject", "weight", "q.length", "q.start", "q.end", "s.length", "s.start", 
            "s.end", "sign")
    colcls = c("character", "character", "numeric", "numeric", "numeric", "numeric", 
               "numeric", "numeric", "numeric", "NULL", "NULL", "character")
    
  } else {
    cln = c("query", "q.length", "q.start", "q.end", "subject", "s.length", "s.start", 
            "s.end","weight", "sign")
    colcls = c("character", "numeric", "numeric", "numeric", "character", "numeric", 
               "numeric", "numeric", "NULL", "numeric", "NULL", "character")
  }
  if (class(blastfile) == "data.frame") {
    blastTable = blastfile
    colnames(blastTable)[12] = "sign"
  } else {
    blastTable = read.table(blastfile, sep = "\t", as.is = TRUE, header = FALSE, 
                            colClasses = colcls, comment.char = "")
    colnames(blastTable) = cln
  }
  ## check for duplicates!
  key = with(blastTable, ifelse(query > subject, paste(query, subject), paste(subject, 
                                                                              query)))
  if (any(duplicated(key))) {
    blastTable = blastTable[!duplicated(key), ]
  }
  seqs = readDNAStringSet(seqfile)
  ## calculate pair completeness for reads before sampling
  if (is.null(pair_completeness)) {
    if (paired) {
      if (seqfile_full != seqfile) {
        seqs_full = readDNAStringSet(seqfile_full)
        pair_counts = tabulate(table(gsub(".$", "", names(seqs_full))))
        rm(seqs_full)
      } else {
        pair_counts = tabulate(table(gsub(".$", "", names(seqs))))
        
      }
      pair_completeness = 1 - pair_counts[1]/sum(pair_counts)
    }else{
      pair_completeness = 0
    }
  }
  NV = length(seqs)  # vertices
  NE = nrow(blastTable)  # nodes
  if (maxv < NV | maxe < NE) {
    ## Sample if graph is large
    V_sample_size = estimate_sample_size(NV, NE, maxv, maxe)
    set.seed(NV)
    seqs = sample(seqs, V_sample_size)
    blastTable = blastTable[blastTable$query %in% names(seqs) & blastTable$subject %in% 
                            names(seqs), ]
  }
  
  blastTable$sign = ifelse(blastTable$sign == "+", 1, -1)
  vnames = unique(c(blastTable$query, blastTable$subject))
  vindex = seq_along(vnames)
  names(vindex) = vnames
  ## correct orientation
  cat("creating graph\n")
  G = graph.data.frame(blastTable[, c("query", "subject", "weight", "sign")], directed = FALSE, 
                       vertices = vnames)
  if (include_layout) {
    ## save temporarily modified blastTable if large
    tmpB = tempfile()
    save(blastTable, file = tmpB)
    rm(blastTable)
    ## 
    cat("graph layout calculation ")
    if (ecount(G) > 2e+06) {
      cat("using fruchterman reingold\n")
      set.seed(vcount(G))
      L = layout.fruchterman.reingold(G, dim = 3)
    } else {
      cat("using OGDF & frucherman reingold\n")
      Ltmp = OGDFlayout(G, alg = c("fmmm"))
      set.seed(vcount(G))
      L = cbind(Ltmp[[1]][, 1:2], layout.fruchterman.reingold(G, dim = 2))
      
    }
    
    GL = list(G = G, L = L)
    save(GL, file = graph_destination)
    if (!is.null(image_file)) {
      cat("exporting graph figure")
      png(image_file, width = 900, height = 900, pointsize = 20)
      plot(GL$G, layout = GL$L[, 1:2], vertex.size = 2, vertex.color = "#000000A0", 
           edge.color = "#00000040", vertex.shape = "circle", edge.curved = FALSE, 
           vertex.label = NA, vertex.frame.color = NA, edge.width = 1)
      dev.off()
      ## create thunmbs:
      system(paste("convert ", image_file, " -resize 100x100 ", image_file_tmb))
      
    }
    rm(GL)
    gc()
    load(tmpB)
    unlink(tmpB)
  }
  
  if (!is.connected(G)) {
    cat("Warning - graph is not connected\n")
    cc = clusters(G, "weak")
    mb = as.numeric(factor(membership(cc), levels = as.numeric(names(sort(table(membership(cc)), 
                                                                          decreasing = TRUE)))))
    selvi = which(mb == 1)  # largest component
    cat("using largest component: \nsize =", round(length(selvi)/vcount(G) * 
                                                   100, 1), "% ;", length(selvi), "reads", "\n")
    G = induced.subgraph(G, vids = selvi)
    blastTable = blastTable[blastTable$query %in% V(G)$name & blastTable$subject %in% 
                                                      V(G)$name, ]
  }
  
  Gmst = minimum.spanning.tree(G)
                                        # create alternative trees - for case that unly suboptima solution is found
  set.seed(123)
  Gmst_alt = list()
  Gmst_alt[[1]] = Gmst
  for (i in 2:6){
    E(G)$weight = runif(ecount(G), 0.1,1)
    Gmst_alt[[i]] = minimum.spanning.tree(G)
  }

  rm(G)
  gc()
  ## six attempts to reorient reads
  flip_names_all = list()
  prop_of_notfit=numeric()
  thr_fit=0.001
  for (ii in 1:6){
    Gmst = Gmst_alt[[ii]]


    blastTable_mst = as.data.frame(get.edgelist(Gmst, names = TRUE))
    colnames(blastTable_mst) = c("query", "subject")
    blastTable_mst$sign = E(Gmst)$sign

    d = dfs(Gmst, root = 1, father = TRUE, order.out = TRUE)
    esign = E(Gmst)$sign
    rc = rep(FALSE, vcount(Gmst))
    j = 0
    p = 0
    for (v in d$order[-1]) {
      j = j + 1
      f = as.numeric(d$father)[v]
      if (is.na(f)) {
        next
      }
      eid = get.edge.ids(Gmst, c(v, f))
      if (esign[eid] == -1) {
        ie = unique(c(which(blastTable_mst$query %in% V(Gmst)$name[v]), which(blastTable_mst$subject %in% 
                                                                              V(Gmst)$name[v])))  # incident edges
        esign[ie] = esign[ie] * -1
        rc[v] = TRUE
        p = p + 1
        if (p > 50) {
          p = 0
          cat("\r")
          cat(j, " : ", sum(esign)/length(esign))
        }
      }
    }
    cat("\t")
    flip_names_all[[ii]] = flip_names = V(Gmst)$name[rc]

    if (!exists("blastTable")) {
      load(tmpB)
      unlink(tmpB)
    }
    ## add nofit later!!
    s.mean = with(blastTable, (s.start + s.end)/2)
    s.mean_corrected = ifelse(blastTable$subject %in% flip_names, blastTable$s.length - 
                                                                  s.mean, s.mean)
    q.mean = with(blastTable, (q.start + q.end)/2)
    q.mean_corrected = ifelse(blastTable$query %in% flip_names, blastTable$q.length - 
                                                                q.mean, q.mean)
    V1 = ifelse(s.mean_corrected > q.mean_corrected, blastTable$subject, blastTable$query)
    V2 = ifelse(s.mean_corrected > q.mean_corrected, blastTable$query, blastTable$subject)
    sign_final = ifelse(V1 %in% flip_names, -1, 1) * ifelse(V2 %in% flip_names, -1, 
                                                            1) * blastTable$sign
    nofit = unique(c(V1[sign_final == "-1"], V2[sign_final == "-1"]))
    prop_of_notfit[ii] = length(nofit)/vcount(Gmst)
    ## check if correctly oriented
    cat("prop notfit", prop_of_notfit[[ii]],"\n")
    if (prop_of_notfit[[ii]]<thr_fit){
      ## OK
      break
    }
  }
  if (!prop_of_notfit[[ii]]<thr_fit){
    ## get best solution
    ii_best = which.min(prop_of_notfit)
    if (ii != ii_best){   # if the last solution was not best, get the best
      flip_names = flip_names_all[[ii_best]]
      s.mean = with(blastTable, (s.start + s.end)/2)
      s.mean_corrected = ifelse(blastTable$subject %in% flip_names, blastTable$s.length - 
                                                                    s.mean, s.mean)
      q.mean = with(blastTable, (q.start + q.end)/2)
      q.mean_corrected = ifelse(blastTable$query %in% flip_names, blastTable$q.length - 
                                                                  q.mean, q.mean)
      V1 = ifelse(s.mean_corrected > q.mean_corrected, blastTable$subject, blastTable$query)
      V2 = ifelse(s.mean_corrected > q.mean_corrected, blastTable$query, blastTable$subject)
      sign_final = ifelse(V1 %in% flip_names, -1, 1) * ifelse(V2 %in% flip_names, -1, 
                                                              1) * blastTable$sign
      nofit = unique(c(V1[sign_final == "-1"], V2[sign_final == "-1"]))
      
    }
  }

  ## exclude all nofit
  
  df2 = data.frame(V1, V2, sign = blastTable$sign, sign_final = sign_final, stringsAsFactors = FALSE)
  rm(blastTable)
  gc()
  vertices = data.frame(name = vnames, reverse_complement = vnames %in% flip_names)
  G = graph.data.frame(df2, vertices = vertices)
  vcount_ori = vcount(G)
  G = induced.subgraph(G, vids = which(!V(G)$name %in% nofit))
  
  G$escore_mst = sum(esign)/length(esign)
  G$escore = sum(sign_final == 1)/length(sign_final)
  cc = clusters(G, "strong")
  mb = as.numeric(factor(membership(cc), levels = as.numeric(names(sort(table(membership(cc)), 
                                                                        decreasing = TRUE)))))
  names(mb) = V(G)$name
  G$loop_index = max(cc$csize)/vcount(G)
  G$prop_of_notfit=min(prop_of_notfit)
  if (is.infinite(G$loop_index) | min(prop_of_notfit)>0.9){ G$loop_index = 0}
  G$coverage = vcount(G)/vcount_ori
  ## check sign in all edges
  V(G)$membership = mb
  save(G, file = directed_graph_destination)
  ## remove nofit
  seqs = seqs[!names(seqs) %in% nofit]
  seqs[names(seqs) %in% flip_names] = reverseComplement(seqs[names(seqs) %in% flip_names])
  names(seqs) = paste(names(seqs), mb[names(seqs)])
  writeXStringSet(seqs, filepath = oriented_sequences)
  ## calculate satellite probability
  if (is.null(satellite_model_path)) {
    pSAT = isSAT = NULL
  } else {
    satellite_model = readRDS(satellite_model_path)
    pSAT = get_prob(G$loop_index, pair_completeness, model = satellite_model)
    isSAT = isSatellite(G$loop_index, pair_completeness, model = satellite_model)
  }
  
  ## get larges cc
  output = list(escore = G$escore, coverage = G$coverage, escore_mts = G$escore_mst, 
                loop_index = G$loop_index, pair_completeness = pair_completeness, graph_file = graph_destination, 
                oriented_sequences = oriented_sequences, vcount = vcount(G), ecount = ecount(G), 
                satellite_probability = pSAT, satellite = isSAT)
  ## clean up
  all_objects = ls()
  do_not_remove = "output"
  rm(list = all_objects[!(all_objects %in% do_not_remove)])
  gc(verbose = FALSE, reset = TRUE)
  return(list2dictionary(output))
}

list2dictionary = function(l){
    dict = "{"
    q='"'
    for (i in 1:length(l)){
        if (class(l[[i]])=="character" | is.null(l[[i]])){
            q2 = "'''"
        }else{
            q2 = ''
        }
        dict = paste0(
            dict,
            q,names(l)[i],q,":",q2, l[[i]], q2,", "
        )
    }
    dict = paste0(dict, "}")
    return(dict)
}

wrap = Vectorize(function(s, width = 80) {
  i1 = seq(1, nchar(s), width)
  i2 = seq(width, by = width, length.out = length(i1))
  return(paste(substring(s, i1, i2), collapse = "\n"))
})


tarean = function(input_sequences, output_dir, min_kmer_length = 11, max_kmer_length = 27, 
                  CPU = 2, sample_size = 10000, reorient_reads = TRUE, tRNA_database_path = NULL, 
                  include_layout = TRUE, paired = TRUE, lock_file=NULL) {
  options(CPU = CPU)
  time0 = Sys.time()
  dir.create(output_dir)
  input_sequences_copy = paste(output_dir, "/", basename(input_sequences), sep = "")
  
  if (!file.copy(input_sequences, input_sequences_copy, overwrite = TRUE)) {
    cat(paste("cannot copy", input_sequences, " to", output_dir), "\n")
    stop()
  }
  
  lock_file = waitForRAM(lock_file = lock_file)
  pair_completeness = NULL
  ## sampling
  if (sample_size != 0) {
    s = readDNAStringSet(input_sequences_copy)
    N = length(s)
    ## pair completness must be calculated before sampling!
    if (N > sample_size) {
      set.seed(N)
      writeXStringSet(sample(s, sample_size), filepath = input_sequences_copy)
      if (paired) {
        pair_counts = tabulate(table(gsub(".$", "", names(s))))
        pair_completeness = 1 - pair_counts[1]/sum(pair_counts)
      }
    }
    rm(s)
  }
  
  if (reorient_reads) {
    input_sequences_oriented = paste0(input_sequences_copy, "oriented.fasta")
    graph_file = paste0(input_sequences_copy, "_graph.RData")
    GLfile = paste0(input_sequences_copy, "_graph.GL")
    cat("reorientig sequences\n")
    
    blastTable = mgblast(input_sequences_copy, input_sequences_copy)
    
    graph_info = mgblast2graph(blastTable, input_sequences_copy, GLfile, graph_file, 
                               input_sequences_oriented, include_layout = include_layout, paired = paired, 
                               pair_completeness = pair_completeness)
    
    ## interupt if it does not look like tandem repeat at all # soft threshold!
    if (is.null(graph_info$pair_completeness)) {
      if (graph_info$loop_index <= 0.4) {
        cat("CLUSTER DID NOT PASS THRESHOLD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return(list2dictionary(list(graph_info = graph_info)))
      }
    } else {
      ## for paired reads:
      if (graph_info$loop_index < 0.7 | graph_info$pair_completeness < 0.4) {
        cat("CLUSTER DID NOT PASS THRESHOLD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return(list2dictionary(list(graph_info = graph_info)))
      }
    }
    
  } else {
    input_sequences_oriented = input_sequences_copy
    graph_info = NULL
  }
  
  
  
  kmer_counts = list()
  kmer_lengths = seq(min_kmer_length, max_kmer_length, by = 4)
  for (i in kmer_lengths) {
    ## use pythonis.null(l[[i]])) function - faster
    cmd = paste(script.dir, "/kmer_counting.py ", input_sequences_oriented, " ", 
                i, sep = "")
    cat("calculation kmers of length ", i, "\n")
    f = system(cmd, intern = TRUE)
    x = read.table(f, as.is = TRUE, sep = "\t")
    kmer_counts[[as.character(i)]] = data.frame(x, freq = x$V2/sum(x$V2))
  }
  
  ## number of kmers:
  N50 = sapply(kmer_counts, function(x) {
    sum(cumsum(x$freq) < 0.5)
  })
  
  N70 = sapply(kmer_counts, function(x) {
    sum(cumsum(x$freq) < 0.7)
  })
  
  
  time1 = Sys.time()
  ggmin = mclapply(kmer_counts, get_mimimal_cc, start = 0.5, mc.cores = CPU)
  time2 = Sys.time()
  cat("kmers graphs created ")
  print(time2 - time1)
  names(ggmin) = names(kmer_counts)
  
  
  ## estimate monomer
  monomers = mclapply(ggmin, estimate_monomer, mc.cores = CPU)
  
  names(monomers) = names(kmer_counts)
  monomers = monomers[!sapply(monomers, is.null)]
  ## error handling:
  error = sapply(monomers, class) == "try-error"
  if (any(error)) {
    cat("\nError in monomers estimation: ")
    cat("calculation failed for monomers length ", names(monomers)[error], "\n")
    print(monomers[error])
    if (any(!error)) {
      monomers = monomers[!error]
    } else {
      stop("monomer estimation failed")
    }
  }
  
  ## summary - make function!!
  total_score = list()
  k = 0
  length_best = numeric()
  score_bn = numeric()
  consensus = character()
  
  for (i in names(monomers)) {
    for (v in seq_along(monomers[[i]]$estimates)) {
      k = k + 1
      total_score[[k]] = c(kmer = as.numeric(i), variant = v, total_score = monomers[[i]]$estimates[[v]]$total_score)
      score_bn[k] = min(rowSums(monomers[[i]]$estimates[[v]]$CM))
      length_best[k] = monomers[[i]]$estimates[[v]]$length_variant_score[1, 
                                                                         1]
      consensus[[k]] = monomers[[i]]$estimates[[v]]$consensus2
    }
  }
  summary_table = as.data.frame(do.call(rbind, total_score))
  summary_table$monomer_length = length_best
  summary_table$consensus_length = nchar(consensus)
  summary_table$score_bn = score_bn
  summary_table$consensus = paste("<pre>", wrap(consensus, width = 80), "<pre>", 
                                  sep = "")
  consensus_DS = DNAStringSet(consensus)
  names(consensus_DS) = with(summary_table, paste0(kmer, "_", variant, "_sc_", 
                                                   signif(total_score), "_l_", monomer_length))
  
  ## filter by score - keep
  
  ## reorder by score
  consensus_DS = consensus_DS[order(summary_table$total_score, decreasing = TRUE)]
  summary_table = summary_table[order(summary_table$total_score, decreasing = TRUE), 
                                ]
  rownames(summary_table) = NULL
  N = nrow(summary_table)
  ## concatenate concensus(ie. make dimer head 2 tail) for pbs detection and other
  ## make something like 'pseudo contig' multimer for mapping - min length 200 bp

  ## searches
  consensus_DS_dimer = DNAStringSet(paste0(consensus_DS, consensus_DS))
  tarean_contigs = DNAStringSet(sapply(consensus_DS,function(x)
    ifelse(nchar(x)<200,
           paste(rep(as.character(x),round(300/nchar(as.character(x))+1)),collapse=''),
           as.character(x)))
    )

  names(consensus_DS_dimer) = names(consensus_DS)
                                        # save sequences:
  consensus_DS_dimer_file = paste0(output_dir, "/consensus_dimer.fasta")
  consensus_DS_file = paste0(output_dir, "/consensus.fasta")
  tarean_contig_file = paste0(output_dir, "/tarean_contigs.fasta")
  writeXStringSet(consensus_DS, consensus_DS_file)
  writeXStringSet(tarean_contigs, tarean_contig_file)
  writeXStringSet(consensus_DS_dimer, consensus_DS_dimer_file)
  if (is.null(tRNA_database_path)) {
    pbs_score = -1
  } else {
    pbs_blast_table = paste0(output_dir, "/pbs_detection")
    pbs_score = detect_pbs(dimers_file = consensus_DS_dimer_file, tRNA_database_path = tRNA_database_path, 
                           reads = input_sequences_oriented, output = pbs_blast_table)
  }
  ## search of open reading frames get the length of the longest
  orf_l = max_ORF_length(consensus_DS_dimer)
  
  
  dir.create(paste0(output_dir, "/img"), recursive = TRUE)
  summary_table$monomer_length_graph = numeric(N)
  summary_table$graph_image = character(N)
  summary_table$monomer_length_logo = numeric(N)
  summary_table$logo_image = character(N)
  ## export graph nd consensus estimate to cluster directory type of output may
  ## change in future
  save(ggmin, file = paste0(output_dir, "/ggmin.RData"))
  save(monomers, file = paste0(output_dir, "/monomers.RData"))
  
  cat("generating HTML output\n")
  for (i in 1:N) {
    kmer = as.character(summary_table$kmer[i])
    variant = summary_table$variant[i]
    ## export graph
    fout_link = paste0("img/graph_", kmer, "mer_", variant, ".png")
    fout = paste0(output_dir, "/", fout_link)
    ## summary_table$monomer_length_graph[i] = summary_table$monomer_length[i]
    ## summary_table$monomer_length_logo[[i]] = nrow(monomers[[kmer]]$estimates[[variant]]$CM)
    summary_table$monomer_length[[i]] = length(monomers[[kmer]]$estimates[[variant]]$consensus)

    if (i <= 10) {
      png(fout, width = 800, height = 800)
      plot_kmer_graph(ggmin[[kmer]], highlight = unlist(monomers[[kmer]]$paths[[variant]]$tr_paths))
      dev.off()
      summary_table$graph_image[i] = hwriteImage(fout_link, link = fout_link, 
                                                 table = FALSE, width = 100, height = 100)
      ## export logo
      png_link = paste0("img/logo_", kmer, "mer_", variant, ".png")
      fout = paste0(output_dir, "/", png_link)
      png(fout, width = 1200, height = round(summary_table$monomer_length[i] * 
                                             1) + 550)
      try(plot_multiline_logo(monomers[[kmer]]$estimates[[variant]]$CM, W = 100))
      dev.off()
      ## export corresponding position probability matrices
      ppm_file = paste0(output_dir, '/ppm_', kmer, "mer_", variant, ".csv")
      ppm_link = paste0('ppm_', kmer, "mer_", variant, ".csv")
      write.table(monomers[[kmer]]$estimates[[variant]]$CM,
                  file = ppm_file,
                  col.names = TRUE, quote = FALSE,
                  row.names = FALSE, sep="\t")
      summary_table$logo_image[i] = hwriteImage(png_link, link = ppm_link, 
                                                table = FALSE, width = 200, height = 100)
    }
    
  }
  
  ## html_report = HTMLReport()
  
  htmlfile = paste0(output_dir, "/report.html")
  cat(htmlheader, file = htmlfile)
  included_columns = c('kmer', 'variant', 'total_score', 'consensus_length','consensus', 'graph_image', 'logo_image')
  summary_table_clean = summary_table[,included_columns]
  colnames(summary_table_clean) = c('k-mer length',
                                    'Variant index',
                                    'k-mer coverage score',
                                    'Consensus length',
                                    'Consensus sequence',
                                    'k-mer based graph',
                                    'Sequence logo')
  HTML(summary_table_clean, file = htmlfile, sortableDF = TRUE)
  HTMLEndFile(file = htmlfile)
  time4 = Sys.time()
  print(time4 - time0)
  if (!is.null(lock_file)){
    print("------removing-lock--------")
    removelock(lock_file)
  }

  print(list(htmlfile = htmlfile, TR_score = summary_table$total_score[1],
             TR_monomer_length = as.numeric(summary_table$consensus_length[1]),
             TR_consensus = summary_table$consensus[1], pbs_score = pbs_score, graph_info = graph_info,
             orf_l = orf_l, tarean_contig_file = tarean_contig_file))
  return(list2dictionary(list(htmlfile = htmlfile, TR_score = summary_table$total_score[1],
              TR_monomer_length = as.numeric(summary_table$consensus_length[1]),
              TR_consensus = summary_table$consensus[1], pbs_score = pbs_score, graph_info = graph_info,
              orf_l = orf_l, tarean_contig_file = tarean_contig_file)))
}


## graph loop index stability
loop_index_instability = function(G) {
  N = 50
  s = seq(vcount(G), vcount(G)/10, length.out = N)
  p = seq(1, 0.1, length.out = N)
  li = numeric()
  for (i in seq_along(s)) {
    print(i)
    gs = induced_subgraph(G, sample(1:vcount(G), s[i]))
    li[i] = max(clusters(gs, "strong")$csize)/vcount(gs)
  }
  instability = lm(li ~ p)$coefficient[2]
  return(instability)
}

isSatellite = function(x, y, model) {
  p = get_prob(x, y, model)
  if (p > model$cutoff) {
    return("Putative Satellite")
  } else {
    return("")
  }
}

get_prob = function(x, y, model) {
  pm = model$prob_matrix
  N = ncol(pm)
  i = round(x * (N - 1)) + 1
  j = round(y * (N - 1)) + 1
  p = pm[i, j]
  return(p)
}


detectMemUsage = function() {
  con = textConnection(gsub(" +", " ", readLines("/proc/meminfo")))
  memInfo = read.table(con, fill = TRUE, row.names = 1)
  close(con)
  memUsage = 1 - (memInfo["MemFree", 1] + memInfo["Cached", 1])/memInfo["MemTotal", 
                                                                        1]
  return(memUsage)
}


makelock<-function(lockfile,lockmsg,CreateDirectories=TRUE){
    lockdir=dirname(lockfile)
    if(!file.exists(lockdir)){
        if(CreateDirectories) dir.create(lockdir,recursive=TRUE)
        else stop("Lock Directory for lockfile ",lockfile," does not exist")
    } 
    if(missing(lockmsg)) lockmsg=paste(system('hostname',intern=TRUE),Sys.getenv("R_SESSION_TMPDIR"))
    if (file.exists(lockfile)) return (FALSE)
                                        # note the use of paste makes the message writing atomic
    cat(paste(lockmsg,"\n",sep=""),file=lockfile,append=TRUE,sep="")
    firstline=readLines(lockfile,n=1)
    if(firstline!=lockmsg){
                                        # somebody else got there first
        return(FALSE)
    } else return(TRUE)
}


removelock<-function(lockfile){
  if(unlink(lockfile)!=0) {
    warning("Unable to remove ",lockfile)
    return (FALSE)
  }
  return (TRUE)
}


waitForRAM = function(p = 0.5,lock_file=NULL) {
  if (detectMemUsage() < p) {
    return(NULL)
    ## check lock file:
  } else {
    cat("waiting for RAM \n")
    free_count = 0
    while (TRUE) {
        if (makelock(lock_file)){
            print("---------locking--------")
            return(lock_file)
        }
      if (detectMemUsage() < p) {
        cat("RAM freed \n")
        return(NULL)
      }
      Sys.sleep(5)
      if (evaluate_user_cpu_usage() == 'free'){
        free_count = free_count + 1
      }else{
        free_count = 0
      }
      if (detectMemUsage() < 0.8 & free_count > 100){
        cat("RAM not free but nothing else is running \n")
        return(NULL)
      }
    }
  }
}

lsmem = function() {
  g = globalenv()
  out_all = envs = list()
  envs = append(envs, g)
  total_size = numeric()
  while (environmentName(g) != "R_EmptyEnv") {
    g <- parent.env(g)
    envs = append(envs, g)
  }
  for (e in envs) {
    
    obj = ls(envir = e)
    if (length(obj) == 0) {
      break
    }
    obj.size = list()
    for (i in obj) {
      obj.size[[i]] = object.size(get(i, envir = e))
    }
    out = data.frame(object = obj, size = unlist(obj.size), stringsAsFactors = FALSE)
    out = out[order(out$size, decreasing = TRUE), ]
    out_all = append(out_all, out)
    total_size = append(total_size, sum(out$size))
  }
  return(list(objects = out_all, total_size = total_size))
} 

evaluate_user_cpu_usage = function(){
  user = Sys.info()["user"]
  a = sum(as.numeric (system(paste ("ps -e -o %cpu -u", user), intern = TRUE)[-1]))
  s = substring (system(paste ("ps -e -o stat -u", user), intern = TRUE)[-1],1,1)
  if (a<5 & sum(s %in% 'D')==0 & sum(s%in% 'R')<2){
    status = 'free'
  }else{
    status = 'full'
  }
  return(status)
}
