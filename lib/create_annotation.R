#!/usr/bin/env Rscript
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")  # this is necessary for handling unicode characters (data.tree package)
suppressPackageStartupMessages(library(data.tree))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(R2HTML))
suppressPackageStartupMessages(library(hwriter))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(plotrix))
suppressPackageStartupMessages(library(png))

source("htmlheader.R")
source("config.R")  # load tandem ranks info
source("utils.R")  
DT_OPTIONS = options = list(pageLength = 1000, lengthMenu = c(10,50,100,1000,5000,10000))
WD = getwd()   # to get script directory when run from Rserve
HTMLHEADER = htmlheader  ## header (character) loaded from htmlheader.R
htmlheader = gsub("Superclusters summary","TAREAN summary", htmlheader)

evaluate_LTR_detection = function(f){
	NO_LTR=NULL
	if (length(readLines(f)) == 11 ){
		return(NO_LTR)
	}
	df=read.table(f,as.is=TRUE,sep="\t", skip=11,fill=TRUE)
	if (ncol(df) != 23){
		#df is smaller if no pbs is detected!
		return(NO_LTR)
	}
  df=df[!df$V13 == "",,drop=FALSE]
	if (nrow(df)==0){
		return(NO_LTR)
	}
	# criteria:
	df_part=df[df$V15 >=12 & df$V20 == 23 & df$V21<df$V20,,drop=FALSE]
	if (nrow(df_part) == 0){
		return(NO_LTR)
	}
  PBS_type = gsub("_","",(str_extract_all(df_part$V13, pattern="_([A-Z][a-z]{2})", simplify=TRUE))) %>%
      paste(collapse=" ")
	return(PBS_type)
}



## annotate superclusters
select_reads_id = function(index, search_type = c("cluster","supercluster")) {
    ## select read if base on the supecluster index need database connection
    ## HITSORTDB!
    search_type = match.arg(search_type)
    x = dbGetQuery(HITSORTDB,
                   paste0("SELECT vertexname FROM vertices WHERE vertexindex IN ",
                          "(SELECT vertexindex  FROM communities ",
                          "WHERE ", search_type,"=\"", index,
                          "\")"))
    return(x$vertexname)
}


get_reads_annotation = function(reads_id) {
    ## select annotation from tables in SEQDB which has name in format *_database
    annot_types = grep("_database", dbListTables(SEQDB), value = TRUE)
    annot = list()
    for (i in annot_types) {
        query = paste0("SELECT * FROM ", i, " WHERE name IN (", paste0("\"", reads_id, 
            "\"", collapse = ", "), ")")
        annot[[i]] = dbGetQuery(SEQDB, query)
    }
    return(annot)
}

supercluster_size = function(supercluster) {
    x = dbGetQuery(HITSORTDB, paste0("SELECT count(*) FROM vertices WHERE vertexindex IN ", 
        "(SELECT vertexindex  FROM communities ", "WHERE supercluster=\"", supercluster, 
        "\")"))
    return(x$"count(*)")
}


cluster_annotation = function(cluster, search_type = c("cluster", "supercluster")){
    ## searcheither for cluster or supercluster annotation
    ## read annotation from sqlite databases database is access though SEQDB (sequence
    ## annotation) and HITSORTDB - clustering information
    search_type = match.arg(search_type)
    reads_id = select_reads_id(cluster, search_type)
    annot = get_reads_annotation(reads_id)
    return(annot)
}

get_tarean_info = function(cluster, search_type = c("cluster", "supercluster")){
    search_type = match.arg(search_type)
    if (search_type == "cluster") {
        search_type = "[index]"
    }
    tarean_info = dbGetQuery(HITSORTDB,
                         paste0(
                             "SELECT [index], supercluster, satellite_probability, tandem_rank, size_real, satellite FROM cluster_info WHERE ",
                             search_type, " = ", cluster))
    nhits = sum(tarean_info$size_real[tarean_info$tandem_rank %in% 1:2])
    proportion = nhits/sum(tarean_info$size_real)
    return( list(
        nhits = nhits,
        proportion = proportion,
        tarean_info = tarean_info)
        )

}

get_ltr_info = function(supercluster){
  ltr_info = dbGetQuery(HITSORTDB,
                        paste0(
                          "SELECT [index], supercluster, ltr_detection, size_real FROM cluster_info WHERE ",
                          "supercluster", " = ", supercluster))
  return(ltr_info)
}

summarize_annotation = function(annot, size = NULL){
    db_id = do.call(rbind, annot)$db_id
    cl_string = gsub("^.+/","",gsub(":.+$", "", gsub("^.+#", "", db_id)))
    weight = do.call(rbind, annot)$bitscore
    domain = str_replace(str_extract(str_replace(db_id, "^.+#", ""), ":.+"), ":", 
                                                          "")
    domain[is.na(domain)] = ""
    domain_table = table(cl_string, domain)
    dt = as.data.frame(domain_table)
    ## remove no count
    dt = dt[dt$Freq != 0, ,drop = FALSE]
    rownames(dt) = paste(dt$cl_string,dt$domain)
    if (!is.null(size)){
        dt$proportion = dt$Freq / size

    }
    return(dt)
}

annot2colors = function(annot, ids, base_color = "#00000050"){
    annot_table = do.call(rbind, annot)
    domains = str_replace(str_extract(str_replace(annot_table$db_id, "^.+#", ""), ":.+"), ":","")
    other_annot = str_replace(annot_table$db_id, "^.+#", "")
    complete_annot = ifelse(is.na(domains), other_annot, domains)
    names(complete_annot) = annot_table$name
    unique_complete_annot = names (table(complete_annot) %>% sort(decreasing = TRUE))
    color_key = unique_colors[seq_along(unique_complete_annot)]
    names(color_key) = unique_complete_annot
    color_table = rep(base_color, length(ids))
    names(color_table) = ids
    color_table[names(complete_annot)] = color_key[complete_annot]
    return(list( color_table = color_table, legend = color_key))
}


 read_annotation_to_tree = function(supercluster, TREE_TEMPLATE = CLS_TREE) {
     ## CLS_TREE is datatree containing 'taxonomy' information of repeats, attribute of
     ## each node/ leave if filled up from annotation annotation is added to leaves
     ## only!! Inner node are NA or NULL
     annot0 = cluster_annotation(supercluster, search_type = "supercluster")
     ## keep only built in annotation
     annot = list(protein_databse = annot0$protein_database, dna_database = annot0$dna_database)
     annot$dna_database$bitscore = annot$dna_database$bitscore / 2 #DNA bitscore corection, blastx - more weight
     annot = do.call(rbind, annot)
     ## remove duplicated hits, keep best
     annot_clean = annot[order(annot$bitscore, decreasing = TRUE), ]
     annot_clean = annot_clean[!duplicated(annot_clean$name), ]
     tarean_info = get_tarean_info(supercluster, search_type = "supercluster")
     ltr_info = get_ltr_info(supercluster)
     db_id = annot_clean$db_id
     cl_string = gsub(":.+$", "", gsub("^.+#", "", db_id))
     weight = annot_clean$bitscore
     domain = str_replace(str_extract(str_replace(db_id, "^.+#", ""), ":.+"), ":", 
        "")
     domain[is.na(domain)] = "NA"
     mean_weight_table = by(weight, INDICES = list(cl_string, domain), FUN = function(x) signif(mean(x)))
     mean_weight_table[is.na(mean_weight_table)] = 0
     total_weight_table = by(weight, INDICES = list(cl_string, domain), FUN = sum)
     total_weight_table[is.na(total_weight_table)] = 0
     cl_table = table(cl_string)
     domain_table = table(cl_string, domain)
     cls_tree = Clone(TREE_TEMPLATE)
     cls_tree$size = supercluster_size(supercluster)
     cls_tree$ltr_info = ltr_info
     for (i in cls_tree$leaves) {
         if (i$full_name %in% names(cl_table)) {
             i$nhits = cl_table[[i$full_name]]
             i$domains = domain_table[i$full_name, ]
             names(i$domains) = colnames(domain_table)
             i$mean_weight = mean_weight_table[i$full_name, ]
             i$total_weight = total_weight_table[i$full_name, ]
             i$proportion = signif(i$nhits/cls_tree$size, 3)
         } else {
             if (i$name == "satellite"){
               i$nhits = tarean_info$nhits
                 i$tandem_rank = tarean_info$tarean_info$tandem_rank
                 i$proportion = tarean_info$proportion
                 i$tarean_info = tarean_info$tarean_info
             }else{
                 i$proportion = 0
             }
         }
     }
     ## create domain string for printing:
     for (i in Traverse(cls_tree)) {
         if (is.null(i$domains)) {
             i$features = ""
         } else if (sum(i$domains) == 0) {
             i$features = ""
         } else {
             dom = i$domains[!names(i$domains) == "NA"]
             i$features = gsub(" *\\(\\)", "", pasteDomains(dom))
         }
     }
     ## add LTR info:
     if  (any(!cls_tree$ltr_info$ltr_detection %in% "None")){
       tr = FindNode(cls_tree, "LTR")
       tr$features = "LTR/PBS"
      }
     return(cls_tree)
 }

pasteDomains = function(dom){
    d = dom[dom != 0]
    if (length(d) == 0){
        return("")
    }else{
        paste0(d, " (", names(d), ")", sep = ", ", collapse="")
    }
    
}

rescale = function(x, newrange) {
    # taken from plotrix package
    if (missing(x) | missing(newrange)) {
        usage.string <- paste("Usage: rescale(x,newrange)\n", "\twhere x is a numeric object and newrange is the new min and max\n", 
            sep = "", collapse = "")
        stop(usage.string)
    }
    if (is.numeric(x) && is.numeric(newrange)) {
        xna <- is.na(x)
        if (all(xna)) 
            return(x)
        if (any(xna)) 
            xrange <- range(x[!xna]) else xrange <- range(x)
        if (xrange[1] == xrange[2]) 
            return(x)
        mfac <- (newrange[2] - newrange[1])/(xrange[2] - xrange[1])
        return(newrange[1] + (x - xrange[1]) * mfac)
    } else {
        warning("Only numeric objects can be rescaled")
        return(x)
    }
}

add_value_to_nodes = function(tr, value_names = c("nhits", "domains", "total_weight", 
    "proportion")) {
    ## propagate values from leaves to nodes, assume that nodes are not defined
    ## (either NA or NULL) leaves coud be be either numeric of list, if list values
    ## are summed up based on the names
    for (vn in value_names) {
        for (i in Traverse(tr)) {
          w = add_leaves_value(i$leaves, vn)
          i[[vn]] = w
        }
    }
}

add_leaves_value = function(x, name) {
    ## check if annotation is multivalue or single value, propagates values from
    ## leaves to root
    n = unique(unlist(lapply(lapply(x, "[[", name), names)))
    if (is.null(n)) {
      ## no names given, we can sum everything to one value
        total = sum(unlist(sapply(x, "[[", name)), na.rm = TRUE)
        return(total)
    } else {
        total = numeric(length(n))
        names(total) = n
        xvals = lapply(x, "[[", name)
        N = 0
        for (i in xvals) {
            for (j in names(i)) {
                total[j] = total[j] + i[j]
                N = N + 1
            }
        }
        return(total)
    }
}


trmap = function(tr, xl = 0, yb = 0, xr = 1, yt = 1, first = TRUE, vertical = TRUE) {
    # treemap for data.tree plotting - experimental
    if (first) {
        plot(xl:xr, yb:yt, type = "n", axes = FALSE, xlab = "", ylab = "")
    }
    if (tr$nhits > 0) {
        width = 2 * (6 - tr$level)
        rect(xl, yb, xr, yt, col = "#00008805", lwd = width, border = "#00000080")
        if (tr$level != 1) {
            text(x = (xl + xr)/2, y = (yb + yt)/2, labels = tr$name, cex = width/4)
        }
    } else {
        return(NULL)
    }
    Nchildren = length(tr$children)
    cat("children:", tr$name, tr$level, "\n")
    if (Nchildren == 0) {
        return(NULL)
    }
    sizes = sapply(tr$children, "[[", "nhits")
    rng = c(sizes, 0) / sum(sizes)
    if (vertical) {
        ## x does not change
        xl2 = rep(xl, Nchildren)
        xr2 = rep(xr, Nchildren)
        yb2 = rescale(rng, c(yb, yt))[1:(Nchildren + 1)]
        yt2 = rescale(rng, c(yb, yt))[-1]
    } else {
        yb2 = rep(yb, Nchildren)
        yt2 = rep(yt, Nchildren)
        xr2 = rescale(rng, c(xl, xr))[1:(Nchildren + 1)]
        xl2 = rescale(rng, c(xl, xr))[-1]
    }
    for (i in 1:Nchildren) {
        trmap(tr$children[[i]], xl2[i], yb2[i], xr2[i], yt2[i], first = FALSE, vertical = !vertical)
    }
    
}

filter_tree = function(tr) {
    ## keep only nodes with positive nhits values
    ## must be used on trees where leaves were already propagated
    ## using add_values_to_nodes
    tr_filtered = Clone(tr)
    Prune(tr_filtered, function(x) x$nhits > 0 | containLTR(x))
    return(tr_filtered)
}

containLTR = function(tr){
  for (i in Traverse(tr)){
    if (i$features == "LTR/PBS"){
      return(TRUE)
    }
  }
  return(FALSE)
}

filter_tree2 = function(tr) {
    tr_filtered = Clone(tr)
    Prune(tr_filtered, function(x) x$parent$nhits > 0)
    return(tr_filtered)
}

## this is for testing purposesm in futurem each node will have its function
## named find_best_hit, it will work like decisin tree.
## functions will be set manually based on training data
formatWidth =  function(x, w){
    if (length (dim(x)) > 1){
        for (i in seq_along(x)){
            delta = nchar(x[,i], type="bytes") - nchar(x[,i], type="char")
            ## delta correction is neccessary for utf-8 correct formating
            x[,i] = sprintf(paste0("%-",w[i] + delta,"s"), x[,i])
        }
        return(x)
    }else{
        delta = nchar(x, type="bytes") - nchar(x, type="char")
        sprintf(paste0("%",w + delta,"s"),
            x) %>% return
    }
}


find_best_hit_repeat = function(cls_tree){
  ## three children:
  ## repeat
  ##  |--rDNA
  ##  |--satellite      this need special handling, rank 1 or 2
  ##  |--mobile_element
  ##
  ## in this case we require that satellite has proportion 0.95
  ## params of good hit
  min_prop_of_positive = 0.7
  min_prop_second_best = 0.9
  min_proportion_x_nhits = 2.5  # based on 0.05 x 50
  min_satellite_proportion = 0.95
  if (isLeaf(cls_tree)) {
    return(cls_tree)
  }
  cond1 = cls_tree$nhits * cls_tree$proportion < min_proportion_x_nhits
  if (cond1){
    return(cls_tree)
  }
  nhits = sapply(cls_tree$children, "[[", "nhits")
  all_hits = cls_tree$root$nhits
  name_of_best_hit = cls_tree$children[[which.max(nhits)]]$name
  best_hit_child = cls_tree$children[[which.max(nhits)]]
  second_max_hits = ifelse(length(nhits) == 1, 0, max(nhits[-which.max(nhits)]))
  cond2 = max(nhits) / sum(nhits) > min_prop_of_positive
  cond3 = max(nhits) / (second_max_hits + max(nhits)) > min_prop_second_best
  cond_satellite = best_hit_child$proportion > min_satellite_proportion & name_of_best_hit == "satellite"
  cond_other = ! name_of_best_hit == "satellite"
  cat(cond2, cond3, cond_satellite, cond_other, "\n")
  if (cond2 & cond3 & (cond_satellite | cond_other)) {
    ## clear case satellite or other
    if ("find_best_hit" %in% names(best_hit_child)){
      ## custom rules for node is defined
      best_hit = best_hit_child$find_best_hit(cls_tree$children[[which.max(nhits)]])
    }else{
                                        # use generic rules
      best_hit = find_best_hit(cls_tree$children[[which.max(nhits)]])
    }
  }else{
    ## do more specific tests for rDNA, or mobile element
    ## rDNA o mobile_elements must be above threshold!
    cond_sat_rank = any(cls_tree$satellite$tandem_rank == 2) & cls_tree$satellite$nhits != 0
    cond_rdna = cls_tree$rDNA$nhits *  cls_tree$rDNA$proportion > min_proportion_x_nhits
    cond_mobile = cls_tree$mobile_element$nhits * cls_tree$mobile_element$proportion  > min_proportion_x_nhits
    cat(cond_sat_rank, "\n")
    cat(cond_rdna,"\n")
    cat(cond_mobile, "\n")
    if (cond_sat_rank  & (cond_rdna | cond_mobile) ){
      cls_tree$satellite$nhits = 0
      best_hit = find_best_hit(cls_tree)
    }else{
      return(cls_tree)
    }
  }
  return(best_hit)
}

find_best_hit = function(cls_tree){
    ## general params of good hit
    min_prop_of_positive = 0.7
    min_prop_second_best = 0.8
    min_proportion_x_nhits = 2.5   # based on 0.05 x 50
    if (isLeaf(cls_tree)) {
        return(cls_tree)
    }


    cond1 = cls_tree$nhits * cls_tree$proportion < min_proportion_x_nhits

    if (cond1){
        ## return if proportions is lower than threshold 
        return(cls_tree)
    }

    nhits = sapply(cls_tree$children, "[[", "nhits")
    best_hit_child = cls_tree$children[[which.max(nhits)]]
  ## more special cases:
    second_max_hits = ifelse(length(nhits) == 1, 0, max(nhits[-which.max(nhits)]))
    cond2 = max(nhits) / sum(nhits) > min_prop_of_positive
    cond3 = max(nhits) / (second_max_hits + max(nhits)) > min_prop_second_best
    if (cond2 & cond3) {
        if ("find_best_hit" %in% names(best_hit_child)){
            ## custom rules for node is defined
            best_hit = best_hit_child$find_best_hit(cls_tree$children[[which.max(nhits)]])
        }else{
            # use generic rules
            best_hit = find_best_hit(cls_tree$children[[which.max(nhits)]])
        }
    } else {
        return(cls_tree)
    }
    return(best_hit)
}


get_annotation_groups = function(tr){
    tr_all = Clone(tr)
    Prune(tr_all, pruneFun=function(...)FALSE) ## prune everithing -> keep root
    tr_all$name = "Unclassified repeat (No evidence)"
    tr_contamination = Clone(tr)$contamination
    tr_organelle = Clone(tr)$organelle
    tr_repeat = Clone(tr)$"repeat"
    tr_repeat$name = "Unclassified_repeat (conflicting evidences)"
    return(
        list(
            tr_repeat = tr_repeat,
            tr_organelle = tr_organelle,
            tr_all = tr_all,
            tr_contamination = tr_contamination
        )
    )
}



format_tree = function(cls_tree, ...) {
    df = ToDataFrameTree(cls_tree, ...)
    ## try to get rid off unicode characters
    df[,1] = gsub("\U00B0", "'", df[,1], fixed = TRUE, useBytes = TRUE)
    df[,1] = gsub("\U00A6", "|", df[,1], fixed = TRUE, useBytes = TRUE)

    df[,1] = gsub("<U+00B0>", "'", df[,1], fixed = TRUE, useBytes = TRUE)
    df[,1] = gsub("<U+00A6>", "|", df[,1], fixed = TRUE, useBytes = TRUE)

    df[,1] = gsub("\xc2\xb0", "'", df[,1], fixed = TRUE, useBytes = TRUE)
    df[,1] = gsub("\xc2\xa6", "|", df[,1], fixed = TRUE, useBytes = TRUE)





    colnames(df)[1] = "                                               "  ## originally levelName
    if ("proportion" %in% colnames(df)){
        df$proportion = signif(df$proportion,2)
    }
    if ("Proportion[%]" %in% colnames(df)){
        df$"Proportion[%]" = round(df$"Proportion[%]", 2)
    }
    ## format header
    w1 = apply(df, 2, function(x) max(nchar(x)))
    w2 = nchar(colnames(df))
    # use whatever with is bigger for formating
    w = ifelse(w1 > w2, w1, w2)
    
    out = character()
    header = mapply(FUN = function(x, y) formatWidth(x, w = y), colnames(df), w) %>%
        paste0(collapse=" | ")
    hr_line = gsub(".","-",header)
    # create output
    # classification lines
    class_lines = formatWidth(df, w) %>%
        apply(1, FUN = function(x) paste0(x,collapse = " | ")) %>%
        paste(collapse = "\n")
    paste(
        header,
        hr_line,
        class_lines,
        sep="\n") %>% return
}


make_final_annotation_template = function(
                                          annot_attributes = list(
                                              Nreads = 0,
                                              Nclusters = 0,
                                              Nsuperclusters = 0,
                                              "Proportion[%]" = 0,
                                              cluster_list = c())
                                          )
{
    ct = Clone(CLS_TREE)
    for (i in Traverse(ct)){
        for (j in names(annot_attributes)){
            i[[j]] = annot_attributes[[j]]
        }
    }
    return(ct)
}


common_ancestor = function(tr1, tr2){
  a1 = tr1$Get('name', traversal = "ancestor")
  a2 = tr2$Get('name', traversal = "ancestor")
  ancestor = intersect(a1,a2)[1]
  return(ancestor)
}

create_all_superclusters_report = function(max_supercluster = 100,
                                           paths,
                                           libdir,
                                           superclusters_dir,
                                           seqdb, hitsortdb,
                                           classification_hierarchy_file,
                                           HTML_LINKS)
{
    ## connect to sqlite databases
    ## seqdb and hitsortdb are path to sqlite files
    HTML_LINKS = nested2named_list(HTML_LINKS)
    paths = nested2named_list(paths)
    report_file = paths[['supercluster_report_html']]
    ## create SEQDB, HTISORTDB and CLS_TREE
    connect_to_databases(seqdb, hitsortdb, classification_hierarchy_file)

    ## append specific classication rules
    CLS_TREE$"repeat"$find_best_hit = find_best_hit_repeat

    ###
  cat("Superclusters summary\n", file = report_file)
  cat("---------------------\n\n", file = report_file, append = TRUE)
    empty = rep(NA,max_supercluster)
    SC_table = data.frame(
        Supercluster = 1 : max_supercluster, "Number of reads" = empty, Automatic_annotation = empty,
        Similarity_hits = empty, TAREAN_annotation = empty,
        Clusters = empty, stringsAsFactors = FALSE, check.names = FALSE
    )
    SC_csv = SC_table  # table as spreadsheet - no html formatings
    final_cluster_annotation = make_final_annotation_template()

  for (sc in 1 : max_supercluster) {
        cat("supercluster: ",sc,"\n")
        cls_tree = read_annotation_to_tree(sc)
        sc_summary = get_supercluster_summary(sc)
        add_value_to_nodes(cls_tree)
        cls_tree_filtered = filter_tree(cls_tree)
        best_hit = find_best_hit(cls_tree)
        ## exception to decision tree put here.
        ## All -> class I ?
        if (best_hit$name == "All"){
          ## check LTR
          if (any(!(unique(cls_tree$ltr_info$ltr_detection) %in% "None"))){
            ## if LTR found - move classification to Class I
            
            LTR = FindNode(best_hit, "LTR")
            nhits_without_class_I = best_hit$nhits - LTR$nhits
            prop_without_class_I = best_hit$proportion - LTR$proportion
            cond3 = nhits_without_class_I * prop_without_class_I < 1
            cond2 = best_hit$nhits * best_hit$proportion < 0.5  # other hits weak
            cond1 = LTR$nhits >= (0.95 * best_hit$nhits) # hits are pre in sub class_I

            if (cond1 | cond2 | cond3){
              best_hit = LTR
            }
          }
        }else{
          ## check is conflict os best hit  with LTR/PBS exists
          if (any(!(unique(cls_tree$ltr_info$ltr_detection) %in% "None"))){
            LTR = FindNode(cls_tree, "LTR")
            ca_name = common_ancestor(LTR, best_hit)
            if (ca_name != "LTR"){
              ## reclassify
              best_hit = FindNode(cls_tree, ca_name)
            }
          }
        }
  
        best_hit_name = best_hit$name
        best_hit_path = paste(best_hit$path, collapse = "/")
        ## add best hit do database:
        sqlcmd = paste0(
            "UPDATE cluster_info SET supercluster_best_hit = \"", best_hit_path, "\" ",
            " WHERE supercluster = ", sc
        )
        dbExecute(HITSORTDB, sqlcmd)
        # add annotation annotation summary stored in final_cluster_annotation - summing up
        best_hit_node = FindNode(final_cluster_annotation, best_hit_name)
        best_hit_node$Nsuperclusters = best_hit_node$Nsuperclusters + 1
        best_hit_node$Nclusters = best_hit_node$Nclusters + sc_summary$Nclusters
        best_hit_node$Nreads = best_hit_node$Nreads + sc_summary$Nreads
        best_hit_node$"Proportion[%]" = best_hit_node$"Proportion[%]" + sc_summary$proportion*100
        best_hit_node$cluster_list = append(best_hit_node$cluster_list, sc_summary$cluster_list)
        best_hit_label = ifelse(best_hit_name == "All", "NA", best_hit_name)
        SC_csv$Automatic_annotation[sc] = best_hit_label
        SC_csv$"Number of reads"[sc] = SC_table$"Number of reads"[sc] = cls_tree$size
        SC_csv$Similarity_hits[sc] = format_tree(cls_tree_filtered,
                                                             "nhits", "proportion", "features")
        SC_table$Similarity_hits[sc] = format_tree(cls_tree_filtered,
                                                               "nhits", "proportion", "features") %>%
            preformatted

        SC_table$Automatic_annotation[sc] = best_hit_label

        G = get_supercluster_graph(
            sc=sc, seqdb = seqdb,hitsortdb = hitsortdb,
            classification_hierarchy_file = classification_hierarchy_file
        )
        ## append tarean annotation
        clist = paste(V(G)$name, collapse =",")
        cl_rank = dbGetQuery(HITSORTDB,
                             paste("SELECT [index], tandem_rank FROM cluster_info WHERE [index] IN (",clist,
                                   ") AND tandem_rank IN (1,2)"))
        ## add information about TR clusters is any
        if (nrow(cl_rank) > 0){
            tarean_annot = paste(
                cl_rank$index, "-",
                RANKS_TANDEM[cl_rank$tandem_rank]
            )
            SC_csv$TAREAN_annotation[sc] = paste(tarean_annot,collapse="\n")
            SC_table$TAREAN_annotation[sc] = mapply(
                hwrite, tarean_annot,
                link = sprintf(HTML_LINKS$ROOT_TO_TAREAN,cl_rank$index)) %>% paste(collapse="<br>")
        }
        ## add links to clusters
        SC_csv$Clusters[sc] = paste((V(G)$name), collapse = ",")
        links = mapply(hwrite, V(G)$name, link = sprintf(HTML_LINKS$ROOT_TO_CLUSTER, as.integer(V(G)$name)))
        SC_table$Clusters[sc] = paste(links, collapse =", ")
        create_single_supercluster_report(G, superclusters_dir)
    }



    ## add html links
    SC_table$Supercluster = mapply(hwrite, SC_table$Supercluster,
                                   link = sprintf(HTML_LINKS$ROOT_TO_SUPERCLUSTER,
                                                  SC_table$Supercluster))



    write.table(SC_csv, file = paths[["superclusters_csv_summary"]],
                sep = "\t", row.names = FALSE)

    
    DT_instance = datatable(SC_table, escape = FALSE, options = DT_OPTIONS)
    saveWidget(DT_instance, file = normalizePath(report_file),
               normalizePath(libdir), selfcontained = FALSE)
    add_preamble(normalizePath(report_file),
                 preamble='<h2>Supercluster annotation</h2> <p><a href="documentation.html#superclust"> For table legend see documentation. <a> </p>')

    annot_groups = get_annotation_groups(final_cluster_annotation)
    saveRDS(object=annot_groups, file = paths[['repeat_annotation_summary_rds']])
    final_cluster_annotation_formated = paste(
        sapply(annot_groups, format_tree, 'Proportion[%]', 'Nsuperclusters', 'Nclusters', "Nreads"),
        collapse = "\n\n\n"
    )
    ## TODO add singleton counts, total counts and extra text - csv output
    ## export annotation summary
    html_summary = start_html(paths[['summarized_annotation_html']], gsub("PAGE_TITLE", "Repeat Annotation Summary", HTMLHEADER))
    html_summary("Repeat annotation summary", HTML.title, HR=2)
    html_summary("This table summarizes the automatic annotations of superclusters that should be verified and manually corrected if necessary. Thus, the table should not be used as the final output of the analysis without critical evaluation.", HTML.title, HR=3)
    html_summary(preformatted(final_cluster_annotation_formated), cat)



    return()
}
## for testing
if (FALSE){
    create_supercluster_report(1:100, report_file, seqdb, hitsortdb, class_file)
}

create_single_supercluster_report = function(G, superclusters_dir){
    sc_dir = paste0(superclusters_dir, "/dir_SC",sprintf("%04d", G$name))
    htmlfile = paste0(sc_dir, "/index.html")
    dir.create(sc_dir)
    title = paste("Supercluster no. ",G$name)
    htmlheader = gsub("PAGE_TITLE", title, HTMLHEADER)
    cat(htmlheader, file = htmlfile)
    file.copy(paste0(WD,"/style1.css"), sc_dir)
    HTML.title(title, file = htmlfile)
    HTML("Simmilarity hits (only hits with proportion above 0.1% in at least one cluster are shown) ", file = htmlfile)
    img_file = paste0(sc_dir,"/SC",sprintf("%0d",G$name), ".png")
    png(filename = img_file, width = 1400, height=1200)
    coords = plot_supercluster(G = G)
    dev.off()
    ## basename - make link relative
    href = paste0("../../clusters/dir_CL",sprintf("%04d", as.numeric(V(G)$name)),"/index.html")
    html_insert_image(basename(img_file), htmlfile,coords=coords, href = href)

    if (is_comparative()){
        HTML("Comparative analysis", file = htmlfile)
        img_file = paste0(sc_dir,"/SC",sprintf("%0d",G$name), "comparative.png")
        png(filename = img_file, width = 1400, height=1200)
        coords = plot_supercluster(G = G, "comparative")
        mtext("comparative analysis")
        dev.off()
        ## basename - make link relative
        href = paste0("../../clusters/dir_CL",sprintf("%04d", as.numeric(V(G)$name)),"/index.html")
        html_insert_image(basename(img_file), htmlfile,coords=coords, href = href)

    }
}

html_insert_image = function(img_file, htmlfile, coords=NULL, href=NULL) {
    img_tag = sprintf('<p> <img src="%s" usemap="#clustermap" > \n</p>\n<br>\n', img_file)
    cat(img_tag, file = htmlfile, append = TRUE)
    if (!is.null(coords)){
        formated_links = character()
        for (i in seq_along(href)){
            formated_links[i] = sprintf('<area shape="circle" coords="%f,%f,%f" href="%s" >\n',
                                        coords[i,1],coords[i,2],coords[,3],href[i])
        }
        cat('<map name="clustermap" >\n',
            formated_links,"</map>\n", file=htmlfile, append = TRUE)

    }
}

html_insert_floating_image = function(img_file,
                                      htmlfile,
                                      width=NULL,
                                      title= "",
                                      footer = ""
                                      ){
    if (is.null(width)){
        width=""
    }
    tag = paste(
        '<div class="floating_img">\n',
        '  <p>', title,'</p>\n',
        '   <a href=',img_file, ">",
        '    <img src="',img_file, '" alt="image" width="',width,'">\n',
        '   </a>',
        '  <p> ',footer,'</p>\n',
        '</div>\n'
       ,sep = "")
    cat(tag, file = htmlfile, append = TRUE)
}

get_supercluster_graph = function(sc, seqdb, hitsortdb, classification_hierarchy_file){
    ## sc - superclusted index
    ## seqdb, hitsortdb - path to sqlite databases
    ## classificationn_hierarchy_file  - path data.tree rds file
    ## SIZE of image
    SIZE = 1000
    connect_to_databases(seqdb, hitsortdb, classification_hierarchy_file)
    clusters = dbGetQuery(HITSORTDB,
                          paste0("SELECT cluster, size FROM superclusters WHERE supercluster=",sc)
                          ) %>% as.data.frame
    ## string for sql query:
    cluster_list = paste0( "(", paste0(clusters$cluster, collapse = ","),")")
    supercluster_ncol = dbGetQuery(HITSORTDB,
                                   paste("SELECT c1,c2,w, k FROM cluster_mate_bond where c1 IN ",
                                         cluster_list, " AND c2 IN ", cluster_list, "AND k > 0.1" 
                                         )
                                   )
    if (nrow(supercluster_ncol)>0){
        ## at least two clusters
        G = graph.data.frame(supercluster_ncol, directed=FALSE)
        L = layout_with_kk(G)
        ## layout is rescaled for easier linking html 
        L[,1] = scales::rescale(L[,1], to = c(1,SIZE) )
        L[,2] = scales::rescale(L[,2], to = c(1,SIZE) )
        G$L = L
    }else{
        ## only one cluster in supercluster
        G = graph.full(n = 1)
        V(G)$name = as.character(clusters$cluster)
        G$L = matrix(c(SIZE/2, SIZE/2), ncol=2)
    }
    G = set_vertex_attr(G, "label", value = paste0("CL",names(V(G))))
    G$name=sc ## names is supercluster
    # create annotation of individual clusters, will be attached as node attribudes
    annot = get_cluster_annotation_summary(clusters)
    ## annotation of nodes(clusters)
    G = set_vertex_attr(G, "annotation", value = annot$clusters[names(V(G))])
    G = set_vertex_attr(G, "size",
                        value = clusters$size[match(names(V(G)), as.character(clusters$cluster))])
    G$annotation = annot$supercluster
    G$clusters = clusters
    # TODO if comparative analysis - add proportion of species
    if (is_comparative()){
        comparative_counts = get_cluster_comparative_counts(cluster_list)
        NACluster = comparative_counts$supercluster
        # some clusters do not have comparative data - they are bellow threshold!
        NACluster[,] = NA
        counts_adjusted = comparative_counts$clusters[names(V(G))]
        for (i in names(V(G))){
            if(is.null(counts_adjusted[[i]])){
                ## get count from database
                seqid = dbGetQuery(HITSORTDB,
                                   paste(
                                       "SELECT vertexname FROM vertex_cluster where cluster = ",
                                       i)
                                   )$vertexname
                codes = get_comparative_codes()$prefix
                x = table(factor(substring(seqid, 1,nchar(codes[1])), levels = codes))
                counts_adjusted[[i]] = data.frame(
                    Freq = as.numeric(x),
                    proportion = as.numeric(x/sum(x)),
                    row.names = names(x))
                
            }
        }

        G = set_vertex_attr(G, "comparative_counts",
                            value = counts_adjusted[V(G)$name]
                            )
        # for whole supercluster
        G$comparative = comparative_counts$supercluster
    }

    return(G)
}

get_cluster_comparative_counts = function(cluster_list){
    counts = dbGetQuery(
        HITSORTDB,
        paste0("SELECT * FROM comparative_counts WHERE clusterindex IN",
               cluster_list
               )
    )
    comparative_counts = apply(counts, 1, function(x)
        y = data.frame(Freq = x[-1], proportion = x[-1]/sum(x[-1]))
    )
    names(comparative_counts) = counts$clusterindex
    total_counts = data.frame(
        Freq = colSums(counts[,-1]),
        proportion = colSums(counts[,-1])/sum(counts[,-1]))
    return(
        list(
            clusters = comparative_counts,
            supercluster = total_counts
            )
    )

}

get_cluster_annotation_summary = function(clusters){
    ## clusters - table of clusters, col: cluster, size
    annot_list = apply(clusters,1,FUN = function(x)summarize_annotation(cluster_annotation(x[1]),x[2]))
    ## if empty - not annotation at all
    if (sum(sapply(annot_list, nrow)) == 0){
        return(list (clusters = NULL,
                 superclusters = NULL
                 ))
    }
    names(annot_list) = as.character(clusters$cluster)
    annot_all = do.call(rbind,annot_list)
    total = by(annot_all$Freq,INDICES=with(annot_all, paste(cl_string,domain)), FUN=sum, simplify=TRUE)
    total_df = data.frame(Freq = c(total), proportion = c(total) /sum(clusters$size))
    ## for ploting it need to contain all annot categories as in s upercluster
    annot_list_full = list()
    for (i in seq_along(annot_list)){
        annot_list_full[[i]] = total_df
        for (j in rownames(total_df)){
            if (!j %in% rownames(annot_list[[i]])){
                annot_list_full[[i]][j,c('Freq','proportion')] = c(0,0)
            }else{
                annot_list_full[[i]][j,c('Freq','proportion')] = annot_list[[i]][j,c('Freq','proportion')]
            }
        }
    }
    names(annot_list_full) = names(annot_list)
    return(list (clusters = annot_list_full,
                 supercluster = total_df
                 )
           )
}

plot_supercluster = function(G,color_scheme=c("annotation","comparative")){
    ## create plot in coords y: 1-1200, x: 1 - 1400
    ## this is fixed for href to work in exact image areas!
    color_scheme = match.arg(color_scheme)
    LIMS0 = matrix(c (1,1000,1,1000), ncol = 2)
    OFFSET_H = 100; OFFSET_W = 200
    lims = LIMS0 + c(0,OFFSET_W * 2, 0, OFFSET_H * 2)
    ## use full range
    par(mar=c(0,0,0,0),xaxs="i", yaxs="i")
    ## proportion of repeats
    if (color_scheme == "annotation"){
        prop = sapply (V(G)$annotation, FUN = function(x)x$proportion)
        brv = c("#BBBBBB",unique_colors)
    }else{
        prop = sapply (V(G)$comparative_counts, FUN = function(x)x$proportion)
        brv = unique_colors
    }
    ## handle special cases - single cluster with annot or single annotation group
    if (is.null(dim(prop))){
        prop=matrix(prop, nrow = 1)
    }
    if (length(prop)==0){
        ## special case - not annotation at all
        prop = matrix(rep(1,vcount(G)),ncol = vcount(G), dimnames=list("NAN", NULL))
    }else{
        if (color_scheme == "annotation"){
            rownames(prop) = rownames(G$annotation)
        }else{
            rownames(prop) = rownames(G$comparative)
        }
        ## reorder by prop
        if (color_scheme =="annotation"){
            prop = prop[order(prop[,1], decreasing = TRUE),,drop=FALSE]
            NAN = 1 - colSums(prop)
            prop = rbind(NAN, prop)
            brv = c("#BBBBBB",unique_colors)
        }else{
        }
    }
    L = G$L  + c(OFFSET_H)
  ## for href - necessary convert coordinater - image is flipped vertically
    include = rowSums(prop>0.001, na.rm = TRUE)>=1
    include[1] = TRUE
    prop = prop[include, , drop=FALSE]
    plot(0,0,type='n', , xlim = lims[,1], ylim = lims[,2])
    plot_edges(G, L, lwd = 8)
    RADIUS = radius_size(G)
    coords = cbind(G$L[,1,drop=FALSE] + 100, 1100 - G$L[,2,drop=FALSE], RADIUS)
    pieScatter(L[,1], L[,2], t(prop), col = brv,
               radius = RADIUS, new.plot = FALSE, border=NA)
    ## add ledend
    legend("topright", legend = rownames(prop),
           col = brv, pch = 15, cex= 1.5, pt.cex= 3)
    ## add cluster labels
    text(x = L[,1], y = L[,2], labels = paste(V(G)$label,"\n",V(G)$size))
    return(coords)
}


radius_size = function(G){
    if (vcount(G) == 1) RADIUS=120
    if (vcount(G) == 2) RADIUS=50
    if (vcount(G) %in% 3:8) RADIUS=40
    if (vcount(G) > 8) RADIUS=30
    return(RADIUS)
}


plot_edges = function(G,L,col="#33000040", lwd = 1){
	e = get.edgelist(G, names = F)
	X0 = L[e[, 1], 1]
	Y0 = L[e[, 1], 2]
	X1 = L[e[, 2], 1]
	Y1 = L[e[, 2], 2]
	segments(X0, Y0, X1, Y1, lwd = lwd,col = col)
}


pieScatter=function(x, y, prop,radius=1, col=NULL, edges=100, new.plot = TRUE, ...){
    if (new.plot){
        plot(x,y,type='n')
    }
    for (i in seq_along(x)){
        p=prop[i,]
        if (length(radius)==1){
            r=radius
        }else{
            r=radius[i]
        }
        include = p != 0
        floating.pie(x[i], y[i],p[include],
                     col=col[include], radius=r,edges=50,...)
    }
}

unique_colors = c("#00FF00", "#0000FF", "#FF0000", 
                  "#FFA6FE", "#FFDB66", "#006401", "#010067", "#95003A",
                  "#007DB5", "#FF00F6", "#FFEEE8", "#774D00", "#90FB92", "#01FFFE",
                  "#0076FF", "#D5FF00", "#FF937E", "#6A826C", "#FF029D",
                  "#FE8900", "#7A4782", "#7E2DD2", "#85A900", "#FF0056",
                  "#A42400", "#00AE7E", "#683D3B", "#BDC6FF", "#263400",
                  "#BDD393", "#00B917", "#9E008E", "#001544", "#C28C9F",
                  "#FF74A3", "#01D0FF", "#004754", "#E56FFE", "#788231",
                  "#0E4CA1", "#91D0CB", "#BE9970", "#968AE8", "#BB8800",
                  "#43002C", "#DEFF74", "#00FFC6", "#FFE502", "#620E00",
                  "#008F9C", "#98FF52", "#7544B1", "#B500FF", "#00FF78",
                  "#FF6E41", "#005F39", "#6B6882", "#5FAD4E", "#A75740",
                  "#A5FFD2", "#FFB167", "#009BFF", "#E85EBE")

if (FALSE){
## For testing
    for ( i in 1:100){
        G = get_supercluster_graph(
            sc=i, seqdb = seqdb,
            hitsortdb = hitsortdb, classification_hierarchy_file = class_file
        )
        png(paste0("/mnt/raid/users/petr/tmp/test.figure_",i,".png"), width = 1400, height=1000)
        plot_supercluster(G)
        dev.off()
        print(i)
    }

}

plotg = function(GG,LL,wlim=NULL,...){
    ## quick ploting function
    e = get.edgelist(GG, names=F)
    w = E(GG)$weight
    if (!is.null(wlim)) {e = e[w > wlim,]; w = w[w > wlim]}
    X0 = LL[e[,1],1]
    Y0 = LL[e[,1],2]
    X1 = LL[e[,2],1]
    Y1 = LL[e[,2],2]
    plot(range(LL[,1]),range(LL[,2]),xlab="",ylab="",axes=FALSE,type="n",...)
    brv = 'grey'
    segments(X0,Y0,X1,Y1,lwd=.5,col=brv)
    points(LL,pch=18,cex=.8,...)
}




get_cluster_info = function(index){
    ## must be connected to hitsort database ( HITSORT)
    x = dbGetQuery(HITSORTDB,
                   paste0("SELECT * FROM cluster_info WHERE [index]=",index))
    return(as.list(x))
}

get_supercluster_info = function(sc){
    ## must be connected to hitsort database ( HITSORT)
    x = dbGetQuery(HITSORTDB,
                   paste0("SELECT * FROM cluster_info WHERE supercluster=",sc))
    return(x)
}

get_supercluster_summary = function(sc){
    info = get_supercluster_info(sc)
    Nclusters = nrow(info)
    Nreads = supercluster_size(sc)
    N = dbGetQuery(SEQDB, "SELECT count(*) from sequences")
    proportion = Nreads/N
    cluster_list = info$index
    return(list(
        Nreads = Nreads,
        Nclusters = Nclusters,
        proportion = proportion,
        cluster_list = cluster_list
    ))
}


get_cluster_connection_info = function(index, search_type=c("pair", "similarity")){
    search_type = match.arg(search_type)
    if (search_type == "pair"){
        x = dbGetQuery(HITSORTDB,
                       sprintf(
                           "SELECT * FROM cluster_mate_connections WHERE c1==%d OR c2=%d",
                           index, index)
                       )
        if (nrow(x) == 0 ){
            ## no connections
            return(NULL)
        }
        other_cl = ifelse(x$c1 == index, x$c2, x$c1)
        out = data.frame(cl = other_cl, N = x$N)
        out$ids = unlist(strsplit(x$ids, split = ", "))
        k = dbGetQuery(HITSORTDB,
                       sprintf(
                           "SELECT * FROM cluster_mate_bond WHERE c1==%d OR c2=%d",
                           index, index)
                       )
        if (k$c1 !=  x$c1 || k$c2 !=x$c2){
            ## tables not in same order
            stop("tables cluster_mate_bond and cluster_mate_connection are not in the same order")
        }
        out$k = k$k
        return(out)
    }else{
        x = dbGetQuery(HITSORTDB,
                       sprintf(
                           "SELECT * FROM cluster_connections WHERE c1==%d OR c2=%d",
                           index, index)
                       )
        if (nrow(x) == 0 ){
            ## no connections
            return(NULL)
        }
        other_cl = ifelse(x$c1 == index, x$c2, x$c1)
        out = data.frame(cl = other_cl, N = x$"count(*)")
        return(out)
    }
}


format_clinfo=function(x){
    include = c(
        "size",
        "size_real",
        "ecount",
        "supercluster",
        "annotations_summary",
        "ltr_detection",
        "pair_completeness",
        "TR_score",
        "TR_monomer_length",
        "loop_index",
        "satellite_probability",
        "TR_consensus",
        "tandem_rank",
        "orientation_score"
    )
    new_names = c(
      "TR_consensus" = "TAREAN consensus",
      "tandem_rank" = "TAREAN_annotation",
      "ltr_detection" = "PBS/LTR",
      'size' = 'number of reads in the graph',
      'size_real' = 'the total number of reads in the cluster',
      'ecount' = "number of edges of the graph:",
      "annotations_summary" = "similarity based annotation"
    )
    values = unlist(x[include]) %>% gsub("\n","<br>", .)
    include = replace(include, match(names(new_names),include), new_names)
    names(values) = include
    ## replace tandem rank with annotation description
    values["TAREAN_annotation"] = RANKS_TANDEM[values["TAREAN_annotation"]]
    ## create HTML table without header
    desc = df2html(data.frame(names(values), values))
    return(desc)
}

create_cluster_report = function(index,
                                 seqdb, hitsortdb,
                                 classification_hierarchy_file,
                                 HTML_LINKS){
    ## index - index of cluster, integer
    ## seqdb, hitsortdb - sqlite databases
    ## classification_hierarchy_file - rds file with classification
    ## other information about cluster is collected from HITSORTDB
    connect_to_databases(seqdb, hitsortdb, classification_hierarchy_file)
    ## basic info about cluster
    clinfo = get_cluster_info(index)
    HTML_LINKS = nested2named_list(HTML_LINKS)
    PNGWIDTH = 1000  # this must be kept so than link in png work correctly
    PNGHEIGHT = 1000
    LEGENDWIDTH = 300
    PS = 30  ## pointsize for plotting
    MAX_N = 10 ## to show mates or similarities connections

    #################################################################################
    ## create HTML title:
    title = paste("Cluster no. ",index)
    htmlheader = gsub("PAGE_TITLE", title, HTMLHEADER)
    cat(htmlheader, file = clinfo$html_report_main)
    file.copy(paste0(WD,"/style1.css"), clinfo$dir)
    HTML.title(title, file = clinfo$html_report_main)
    ## print basic info to HTML header:
    ## make link back to cluster table
    hwrite("Go back to cluster table <br><br>\n", link=HTML_LINKS$CLUSTER_TO_CLUSTER_TABLE) %>%
        cat(file = clinfo$html_report_main, append =TRUE)
    ## create ling to corresponding supercluster
    cat("Cluster is part of ",
        hwrite (paste("supercluster: ", clinfo$supercluster),
                link=sprintf(HTML_LINKS$CLUSTER_TO_SUPERCLUSTER, clinfo$supercluster)),
        "\n<br><br>\n",
        file = clinfo$html_report_main, append =TRUE)
    HTML.title("Cluster characteristics:", HR = 3, file = clinfo$html_report_main)
    cat(format_clinfo(clinfo), file = clinfo$html_report_main, append =TRUE)
    #################################################################################
    ## add link to image with graph layout with domains
    fs = list(
        graph_domains = paste0(clinfo$html_report_files,"/graph_domains.png"),
        graph_mates = paste0(clinfo$html_report_files,"/graph_mates%04d_%04d.png"),
        graph_base = paste0(clinfo$html_report_files,"/graph_base.png"),
        graph_comparative = paste0(clinfo$html_report_files,"/graph_comparative.png")
    )
    fs_relative = list(
        graph_domains = paste0(basename(clinfo$html_report_files),"/graph_domains.png"),
        graph_comparative = paste0(basename(clinfo$html_report_files),"/graph_comparative.png"),
        graph_mates = paste0(basename(clinfo$html_report_files),"/graph_mates%04d_%04d.png")
    )

    dir.create(clinfo$html_report_files)
    load(clinfo$graph_file)
    annot = cluster_annotation(index)
    annot_summary = summarize_annotation(annot, clinfo$size_real)
    vertex_colors = annot2colors(annot,ids = V(GL$G)$name)

    color_proportions = sort(table(vertex_colors$color_table)) / clinfo$size_real
    ## hits with smaller proportion are not shown!
    exclude_colors = names(color_proportions)[color_proportions < 0.001]
    vertex_colors$color_table[vertex_colors$color_table %in% exclude_colors] = "#00000050"
    vertex_colors$legend = vertex_colors$legend[!vertex_colors$legend %in% exclude_colors]

    if (is_comparative()){
        comp_codes = get_comparative_codes()$prefix
        colors = unique_colors[seq_along(comp_codes)]
        names(colors) = comp_codes
        species = substring(V(GL$G)$name,1, nchar(comp_codes[1]))
        ## create image with highlighted domains
        png(fs$graph_comparative, width = PNGWIDTH + LEGENDWIDTH, height = PNGHEIGHT, pointsize = PS)
        layout(matrix(1:2, ncol = 2), width = c(10, 3))
        plotg(GL$G,GL$L, col = colors[species])
        par(mar = c(0, 0, 0, 0))
        plot.new()
        legend("topleft", col=colors,
                   legend = names(colors),
                   pch = 15, cex = 0.7)
        dev.off()
        HTML.title("comparative analysis:", HR=4, file = clinfo$html_report_main)
        html_insert_image(
            img_file = fs_relative$graph_comparative,
            htmlfile = clinfo$html_report_main)
        ## comparative analysis - calculate statistics:
        V1V2 = substring(get.edgelist(GL$G),1, nchar(comp_codes[1]))
        C_observed = table(data.frame(
            V1=factor(V1V2[,1],levels = comp_codes),
            V2=factor(V1V2[,2],levels = comp_codes)
        ))
        C_observed = as.data.frame.matrix(C_observed + t(C_observed))
        P_observed = as.data.frame.matrix(C_observed/sum(C_observed))
        Edge_propotions = table(factor(V1V2, levels = comp_codes))/(length(V1V2))
        P_expected = matrix(Edge_propotions,ncol=1) %*% matrix(Edge_propotions, nrow=1)

        spec_count_df = data.frame(table(factor(species, levels=comp_codes)))
        spec_counts = df2html(
            spec_count_df,
            header =  c("Species", "Read count")
        )
        edge_count_df = data.frame(species = comp_codes, (C_observed/2)[,])
        edge_count = df2html(
            edge_count_df,
            header = c("", comp_codes)
        )

        P_OE_df = data.frame(species=comp_codes,(P_observed/P_expected[,]))
        P_OE = df2html(
            P_OE_df,
            header = c("", comp_codes)
        )
        ## export csv table
        write.table(edge_count_df, file = paste0(clinfo$dir,"/edge_count.csv"), sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE )
        write.table(spec_count_df, file = paste0(clinfo$dir,"/species_count.csv"), sep="\t", col.names = FALSE, row.names = FALSE , quote = FALSE )
        write.table(P_OE_df, file = paste0(clinfo$dir,"/observed_expected_number_of_edges.csv"), sep="\t", col.names = TRUE, row.names=FALSE, quote = FALSE )

        HTML.title("Comparative analysis - species read counts:", HR =5, file = clinfo$html_report_main)
        cat(spec_counts,
            file = clinfo$html_report_main, append =TRUE)

        # show connection only if more than one species in cluster
        if (length(unique(species))>1){
            HTML.title("comparative analysis - number of edges between species:",
                       HR =5, file = clinfo$html_report_main)
            cat(edge_count,
                file = clinfo$html_report_main, append = TRUE)

            HTML.title("comparative analysis - observed/expected number of edges between species",
                       HR =5, file = clinfo$html_report_main)
            cat(P_OE,
                file = clinfo$html_report_main, append = TRUE)
        }
        
    }


    ## create image with highlighted domains
    png(fs$graph_domains, width = PNGWIDTH + LEGENDWIDTH, height = PNGHEIGHT, pointsize = PS)
    layout(matrix(1:2, ncol = 2), width = c(10, 3))
    plotg(GL$G,GL$L, col = vertex_colors$color_table)
    domains_detected = length(vertex_colors$legend) > 0
    par(mar = c(0, 0, 0, 0))
    plot.new()
    if (domains_detected){
        # domains found
        legend("topleft", col=vertex_colors$legend,
               legend = names(vertex_colors$legend),
               pch = 15, cex = 0.7)
    }
    dev.off()

    HTML.title("protein domains:", HR=4, file = clinfo$html_report_main)
    if (!domains_detected){
        HTML("No protein domains detected", file = clinfo$html_report_main)
    }
    HTML("protein domains:", HR=4, file = clinfo$html_report_main)
    html_insert_image(
        img_file = fs_relative$graph_domains,
        htmlfile = clinfo$html_report_main)

    #############################################################################
    if (nrow(annot_summary) == 0){
    HTML.title("Reads annotation summary", HR = 3, file = clinfo$html_report_main)
        HTML("No similarity hits to repeat databases found", file = clinfo$html_report_main)
    }else{
        HTML(annot_summary, file = clinfo$html_report_main, align = "left")
    }


    ## similarity and mate cluster
    mate_clusters = get_cluster_connection_info(index, search_type="pair")
    similar_clusters = get_cluster_connection_info(index, search_type="similarity")
    ## report mate and similarity clusters
    if (!is.null(similar_clusters)){
        HTML.title("clusters with similarity:", file =clinfo$html_report_main, HR = 3)
        cat(df2html(
            similar_clusters,
            header=c("Cluster","Number of similarity hits"),
            sort_col = "N", scroling = TRUE
            ),
        file =clinfo$html_report_main, append=TRUE)
    }
    if (!is.null(mate_clusters)){
        HTML.title("clusters connected through mates:", file =clinfo$html_report_main, HR = 3)
        cat(df2html(
            mate_clusters[,c('cl','N','k')],
            header = c('Cluster','Number of shared<br> read pairs','k'),
            sort_col = "N", scroling = TRUE
            ),
        file = clinfo$html_report_main,append = TRUE
        )

        ## create base graph images - it will serve as background for
        ## mate clusters plots
        png(fs$graph_base,
            width = PNGWIDTH, height = PNGHEIGHT, pointsize = PS)
        par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
        plotg(GL$G,GL$L, col = "#00000050")
        dev.off()
        ## load base as raster image
        base_image = readPNG(fs$graph_base)

        for (i in order(mate_clusters$N, decreasing = TRUE)){
            mate_ids = unlist(strsplit(mate_clusters$ids[[i]],split=","))
            ## print only graph above MAX_N mates
            if (length(mate_ids) < MAX_N){  # TODO  - use constant
                next
            }
            png(sprintf(fs$graph_mates,index, mate_clusters$cl[i]),
                width = PNGWIDTH, height = PNGHEIGHT, pointsize = PS)
            color_mate = gsub(".$","",V(GL$G)$name) %in% mate_ids %>%
                ifelse("#FF0000FF", "#000000AA")
            par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
            plot(range(GL$L[,1]), range(GL$L[,2]), type = "n", xlab = "", ylab = "", axes = FALSE,
                 main = paste0("CL",index," ----> CL",mate_clusters$cl[i]))
            rasterImage(base_image,
                        range(GL$L[,1])[1], range(GL$L[,2])[1],
                        range(GL$L[,1])[2], range(GL$L[,2])[2]
                        )
            points(GL$L[,1:2], col = color_mate,pch=18,cex=.8)
            dev.off()
            title = paste0("CL",index," ----> CL",mate_clusters$cl[i])
            footer = paste0("No. of shared pairs: :", mate_clusters$N[i])
            html_insert_floating_image(
                img_file = sprintf(fs_relative$graph_mates, index, mate_clusters$cl[i]),
                htmlfile = clinfo$html_report_main, width = 200,
                title = title, footer = footer
            )
        }
    }
}
