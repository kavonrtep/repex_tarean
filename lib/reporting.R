#!/usr/bin/env Rscript
library(R2HTML)
library(hwriter)
library(DT)
library(tools)

source("htmlheader.R")
source("config.R")  # load TANDEM_RANKS
source("utils.R")
DT_OPTIONS = list(pageLength = 1000, lengthMenu = c(10, 50, 100, 1000, 5000, 10000))
HTMLHEADER_TAREAN = gsub("PAGE_TITLE","TAREAN summary", htmlheader)
HTMLHEADER_INDEX = gsub("PAGE_TITLE","Clustering summary", htmlheader)

WD = getwd()  # to get script directory when run from Rserve

reformat_header = function(df){
    H = colnames(df)
    H[H=="TR_score"] = "TAREAN k-mer_coverage"
    H[H=="vcount"] = "|V|"
    H[H=="ecount"] = "|E|"
    H[H=="Genome_Proportion[%]"] = "Proportion[%]"
    H[H=="Proportion_Adjusted[%]"] = "Proportion adjusted[%]"
    H[H=="supercluster"] = "Super_cluster"
    H[H=="size_real"] = "Number of reads"
    H[H=="TR_monomer_length"] = "Consensus_length"
    H[H=="TR_consensus"] = "Consensus"
    H[H=="pbs_score"] = "PBS score"
    H[H=="ltr_detection"] = "LTR detection"
    H[H=="kmer_analysis"] = "TAREAN k-mer analysis"
    
   # H[H=="annotations_summary"] = "Similarity_hits"
    H[H=="annotations_summary"] = "Similarity_hits_[above 0.1%]"
    H[H=="annotations_summary_custom"] = "Similarity_hits_to_custom_database"
    H[H=="loop_index"] = "connected_component_index C"
    H[H=="pair_completeness"] = "pair_completeness_index_P"
    H = gsub("_","<br>",H)
    H=gsub("TR_","",H)
    H = capitalize(H)
    colnames(df) = H
    return(df)
}

reformat4html=function(df){
    for (n in colnames(df)){
        if (class(df[,n]) == 'character'){
            df[,n] = gsub("\n","<br>", df[,n])
        }
        if (class(df[,n]) == 'numeric'){
            df[,n] = signif(df[,n],3)
        }
    }
    return(df)
}

capitalize = function(s){
    paste(toupper(substring(s, 1, 1)),
          substring(s, 2),
          sep="")
}


create_main_reports = function(paths, N_clustering, N_input,N_omit, merge_threshold,
                               paired, consensus_files, custom_db, tarean_mode,
                               HTML_LINKS, pipeline_version_info, max_memory,
                               max_number_reads_for_clustering, mincln){
    ## this create main html index and also tarean report ##
    ## index and tarean html reports are created always
    ## extract all paths and directories
    HTML_LINKS = nested2named_list(HTML_LINKS)
    paths = nested2named_list(paths)
    csvfile = paths[['clusters_info']]
    clusters_summary_csv = paths[['clusters_summary_csv']]
    profrep_classification_csv = paths[['profrep_classification_csv']]
    htmlfile = paths[["tarean_report_html"]]
    html_report_dt = paths[["cluster_report_html"]]
    main_report = paths[["main_report_html"]]
    summarized_annnotation_html = paths[["summarized_annotation_html"]]
    libdir = paths[['libdir']]
    clusters_dir = paths[["clusters__relative"]]
    superclusters_dir = paths[['superclusters__relative']]
    seqdb = paths[['sequences_db']]
    hitsortdb = paths[['hitsort_db']]
    connect_to_databases(seqdb, hitsortdb)
    dfraw = read.table(csvfile, as.is=TRUE, header=TRUE, sep="\t", na.strings = c('None','NA'))
    # table must be updated
    dfraw$supercluster_best_hit = dbGetQuery(HITSORTDB, "SELECT supercluster_best_hit FROM cluster_info")[, 1]
    ## columns to use
    selected_cols = c("index", "size_real","size_adjusted", "vcount","ecount",
                     "loop_index", "pair_completeness",
                    'satellite_probability','satellite',
                    'TR_score','pbs_score','ltr_detection', 'TR_monomer_length',
                    'TR_consensus', "annotations_summary", "supercluster", 'tandem_rank',
                    'supercluster_best_hit')

    ## some columns are added (like Graph_layout, clusters,...)
    ## columns for html report
    selected_cols_tarean = c(
        "Cluster",
        "Proportion[%]",
        "Proportion_Adjusted[%]",
        "size_real",
        'satellite_probability',
        'TR_monomer_length',
        'TR_consensus',
        'Graph_layout',
        'kmer_analysis',
        "loop_index",
        "pair_completeness",
        'TR_score',
        "vcount",
        "ecount",
        'pbs_score',
        "annotations_summary"
    )
    selected_cols_main = c(
        "Cluster",
        "supercluster",
        "Proportion[%]",
        "Proportion_Adjusted[%]",
        "size_real",
        'Graph_layout',
        "annotations_summary",
        'ltr_detection',
        'satellite_probability',
        'TAREAN_annotation',
        'TR_monomer_length',
        'TR_consensus',
        'kmer_analysis',
        "loop_index",
        "pair_completeness",
        'TR_score',
        "ecount",
        "vcount"
    )

    if (custom_db){
        selected_cols_main = c(selected_cols_main, "annotations_summary_custom")
        selected_cols_tarean = c(selected_cols_tarean, "annotations_summary_custom")
        selected_cols = c(selected_cols, "annotations_summary_custom")
    }
    if (is_comparative()){
        prefix_codes = dbGetQuery(SEQDB, "SELECT * FROM prefix_codes")
        species_counts = dbGetQuery(HITSORTDB, "SELECT * FROM comparative_counts")
        superclusters = dbGetQuery(HITSORTDB,paste(
                                                 "SELECT supercluster, cluster FROM superclusters WHERE cluster <=",
                                                 nrow(species_counts))
                                   )
        species_counts = merge(superclusters, species_counts, by.x = "cluster", by.y = "clusterindex")
        ## include commented header with total counts:
        cat("# Total counts:\t\t", paste(prefix_codes$N, collapse="\t"),"\n#\n",
            sep="",
            file = paths[['comparative_analysis_counts_csv']])

        write.table(species_counts, file = paths[['comparative_analysis_counts_csv']],
                    sep = "\t", col.names = TRUE, row.names = FALSE, append=TRUE)
        species_counts_formated = apply(
            species_counts[, prefix_codes$prefix, drop = FALSE],
            1, function(x) paste(prefix_codes$prefix, ":", x, "\n",sep='', collapse=""))
        dfraw$species_counts = species_counts_formated[1:nrow(dfraw)]
        selected_cols = c(selected_cols, "species_counts")
        selected_cols_main = c(selected_cols_main, "species_counts")
    }
    

    df_report = dfraw[,selected_cols]
    ## describe tandem ranks:
    df_report$TAREAN_annotation = RANKS_TANDEM[as.character(df_report$tandem_rank)]
    ## remove Cluster_similarity_hits
    df_report_csv = reformat_df_report(df_report) 
    df_report_csv = df_report_csv[,!colnames(df_report_csv) %in% "Cluster_similarity_hits"]
    df_report_csv$Final_annotation=""

    ## make table for profrep classification
    write.table(
        reformat_df_to_profrep_classification(df_report), file = profrep_classification_csv,
        sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

    df_report$"kmer_analysis" = ifelse(dfraw$putative_tandem, hwrite("report", link = dfraw$html_tarean), "N/A")
    df_report$"Graph_layout" = hwriteImage(dfraw$image_file_tmb, link = dfraw$image_file, table = FALSE)
    df_report$Cluster = paste0("CL", df_report$index)
    df_report$"Proportion[%]" = signif (100 * df_report$size_real / N_clustering, 2)
    df_report$"Proportion_Adjusted[%]" = signif (100 * df_report$size_adjusted / N_clustering, 2)
    if (!tarean_mode){

        df_report$Cluster=sapply(
            df_report$index,
            function(x) hwrite(x, link = sprintf("%s/dir_CL%04d/index.html", clusters_dir, x)))

        df_report$supercluster = sapply(
            df_report$supercluster,
            function(x) hwrite(x, link = sprintf("%s/dir_SC%04d/index.html", superclusters_dir, x)))
    }
    ## TAREAN report
    ## copy tarean output data help to place nad make link to it
    file.copy(paste0(WD,"/style1.css"), dirname(htmlfile))
    file.copy(paste0(WD,"/documentation.html"), dirname(htmlfile))

    tarean_html = start_html(htmlfile, HTMLHEADER_TAREAN)
    tarean_html("Tandem Repeat Analyzer", HTML.title, HR=1)
    tarean_html = start_html(htmlfile, HTMLHEADER_TAREAN)
    tarean_html('Run statistics:', HTML.title, HR=2)
    tarean_html(paste("Number of input reads:", N_input ))
    tarean_html(paste("Number of analyzed reads:", N_clustering))
    tarean_html(paste("Cluster merging:",ifelse(merge_threshold == 0,"No", "Yes")))

    ## export links to consensus sequecnes in fasta files
    tarean_html("Consensus files - fasta format:", HTML.title, HR=2)
    for (i in TANDEM_RANKS[TANDEM_RANKS != 0]){ ## no consensus for rank 0
        if (!is.null (consensus_files[[i]])){
            N = sum(dfraw$tandem_rank == TANDEM_RANKS[i])
            name_string = paste(names(TANDEM_RANKS)[i]," - total ", N, "found" )
            tarean_html(paste("<p>",
                  hwrite(name_string, download = name_string,
                         link = basename(consensus_files[[i]][[1]])),
                  "<br>\n"))

        }
    }
    ## print link to documentation ##
    tarean_html("Documentation", HTML.title, HR=2)
    tarean_html(paste('<p> For the explanation of TAREAN output see',
        ' <a href="documentation.html#tra" > the help section </a> <p>'))
    ## HOW TO CITE section)

    ## PRINT TABLES WITH CLUSTERS
    for (n in names(TANDEM_RANKS)){
      tarean_html(n, HTML.title, HR=2)
      inc <- dfraw$tandem_rank == TANDEM_RANKS[n]
      if (sum(inc > 0)){
          tarean_html(reformat4html(
              reformat_header(
                  df_report[inc, selected_cols_tarean ,drop=FALSE]
              )
          ),
          align = "left", digits = 3)
      }else{
        tarean_html("not found")
      }
    }
    
    ## export table with all cluster 
    cat("",file = html_report_dt)

    DT_instance =  df_report[,selected_cols_main, drop = FALSE] %>%
        reformat_header %>% reformat4html %>% datatable(escape = FALSE, options = DT_OPTIONS) %>%
        formatStyle(columns = seq_along(selected_cols), "font-size" = "12px") %>%
        formatStyle(columns = "Similarity<br>hits<br>[above 0.1%]", "min-width" = "500px")

    saveWidget(DT_instance, file = normalizePath(html_report_dt),
               libdir=normalizePath(libdir) , selfcontained = FALSE)

    add_preamble(normalizePath(html_report_dt),
                 preamble='<h2>Cluster annotation</h2> <p><a href="documentation.html#clust"> For table legend see documentation. <a> </p>')

    ## Main page - Clustering info - global information about clustering
    top_clusters_prop = sum(df_report$size_real)/N_clustering
    clustering_info = summary_histogram(fn = paths[["summary_histogram"]], N_clustering, N_omit, df_report$size_adjusted,
                    top_clusters_prop)
    index_html = start_html(main_report, HTMLHEADER_INDEX)
    index_html("Clustering Summary", HTML.title, HR = 1)
    
    index_html(paste0('<a href="',
                      paths[['summary_histogram__relative']],
                      '"> <img src="', paths[['summary_histogram__relative']],
                      '" width="700" border="1" >',
                      ' </a>'), cat)
  index_html('<p> <b> Graphical summary of the clustering results. </b> Bars represent superclusters, with their heights and widths corresponding to the numbers of reads in the superclusters (y-axis) and to their proportions in all analyzed reads (x-axis), respectively. Rectangles inside the supercluster bars represent individual clusters. If the filtering of abundant satellites was performed, the affected clusters are shown in green, and their sizes correspond to the adjusted values. Blue and pink background panels show proportions of reads that were clustered and remained single, respectively. Top clusters are on the left of the dotted line. </p><hr><br><br>',cat)

    index_html('Run information:', HTML.title, HR = 2)
    index_html(paste("Number of input reads:", N_input ))
    index_html(paste("Number of analyzed reads:", N_clustering))
    if (N_omit != 0){
      index_html(paste("Number of reads removed by automatic filtering of abundant putative satellites:", N_omit))
      index_html(paste("Number of remaining reads after filtering of abundant satellites:", N_clustering - N_omit ))
    }
  
    index_html(
        paste(
            "Proportion of reads in top clusters :",
            signif(100 * sum(df_report$size_real)/N_clustering,2),
            "%"
        ))
    index_html(paste("Cluster merging:",ifelse(merge_threshold == 0,"No", "Yes")))
    index_html(paste("Paired-end reads:",ifelse(paired, "Yes", "No")))
    index_html("Available analyses:", HTML.title, HR=2)
    index_html(paste("<p>",hwrite("Tandem repeat analysis", link = HTML_LINKS$INDEX_TO_TAREAN),"</p>"),cat)

    if (!tarean_mode){
        index_html(paste("<p>", hwrite("Cluster annotation", link = HTML_LINKS$INDEX_TO_CLUSTER_REPORT),"</p>"),cat)
        index_html(paste("<p>", hwrite("Supercluster annotation",
                                       link = HTML_LINKS$INDEX_TO_SUPERCLUSTER_REPORT),"</p>"),cat)
        index_html(paste("<p>", hwrite("Repeat annotation summary", link = HTML_LINKS$INDEX_TO_SUMMARIZED_ANNOTATION),"</p>"),cat)
    }

    index_html("Supplementary files:", HTML.title, HR=2)
    index_html(paste("<p>", hwrite("CLUSTER_TABLE.csv",
                                   link = paths[["clusters_summary_csv__relative"]]),"</p>"),cat)
    if (!tarean_mode){

      index_html(paste("<p>", hwrite("SUPERCLUSTER_TABLE.csv",
                                     link = paths[["superclusters_csv_summary__relative"]]),"</p>"),cat)
      index_html(paste("<p>", hwrite("contigs.fasta",
                                     link = paths[["contigs__relative"]]),"</p>"),cat)
    }
  
    if (is_comparative()) {
      index_html(paste("<p>", hwrite("COMPARATIVE_ANALYSIS_COUNTS.csv",
                                     link = paths[["comparative_analysis_counts_csv__relative"]]),"</p>"),cat)
    }
  print("---PATH---")
  print(paths)

  if (is_comparative()) {
    tryCatch({
      imagemap = plot_rect_map(
        read_counts = paths[['comparative_analysis_counts_csv']],
        cluster_annotation = paths[['profrep_classification_csv']],
        output_file = paths[['comparative_summary_map']]
      )},
      error = function(err){
        print(paste("error while plotting ", err))
      }
    )

    HTML.title("Comparative analysis - Total number of reads in clustering analysis", file = main_report)
    index_html(df2html(
      prefix_codes,
      header = c("Code", "Total read count"), rounding_function = round),
      cat
      )
    HTML.title("Comparative analysis - Number of reads in individual clusters", file = main_report)

    index_html(paste0('<img src="', paths[['comparative_summary_map__relative']],
                      '" usemap ="#clustermap" border="2">'), cat)

    index_html(
      "Bar plot on top shows the size of individual clusters. Size of the rectangles in lower panel is proportional to the number of reads in a cluster for each species. Clusters and species were sorted using hierarchical clustering. Bars and rectangles in the plot are hyperlinked to the individual cluster reports.")
    index_html(imagemap)
  }

  how2cite = readLines(paths[["how_to_cite"]])

  index_html(how2cite, cat, sep="\n")
  index_html("<br><hr>", cat)
  index_html('Details:', HTML.title, HR = 3)
  index_html(pipeline_version_info %>% preformatted, cat)
  index_html(paste0("Minimal number of reads in cluster to be considered top cluster : ", mincln))
  index_html(paste0("Reserved Memory : ", round(max_memory/(1024*1024)), "G"))
  index_html(paste0("Maximum number of processable reads with the reserved memory : ", max_number_reads_for_clustering))


  ## export to csv
  clustering_info$Number_of_analyzed_reads = N_clustering
  write.table(t(as.data.frame(clustering_info)),
              file = clusters_summary_csv, sep="\t", col.names = FALSE)
  cat("\n", file = clusters_summary_csv, append = TRUE)
  write.table(
    df_report_csv, file = clusters_summary_csv,
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = TRUE, append=TRUE)
} 

dummy_function = function(){
    print("dummy function")
}

reformat_df_report = function(df_report){
    # for printing to csv  - this should be consise
    df_report$TR_consensus = gsub("(<pre>)|(</pre>)","",df_report$TR_consensus)
    df_report$tandem_rank = NULL
    ## make suitable order and rename
    if ("annotations_summary_custom" %in% colnames(df_report)){
        custom = "annotations_summary_custom"
    }else{
        custom=character()
    }
    df_out = df_report[,c('index',
                          'supercluster',
                          'size_real',
                          'size_adjusted',
                          'supercluster_best_hit',
                          'TAREAN_annotation',
                          'annotations_summary',
                          custom)
                       ]

    colnames(df_out) = c('Cluster',
                         'Supercluster',
                         'Size',
                         'Size_adjusted',
                         'Automatic_annotation',
                         'TAREAN_annotation',
                         'Cluster_similarity_hits',
                         custom)
    return(df_out)
}

reformat_df_to_profrep_classification = function(df_report){
    CL = df_report$index
    best_hit = df_report$supercluster_best_hit
    ## format conversion(in order):
    replacement = list(
        c("/", "|"),
        c("Ty1_copia", "Ty1/copia"),
        c("Ty3_gypsy", "Ty3/gypsy"),
        c("TatIV_Ogre", "TatIV/Ogre"),
        c("Ogre_Tat", "Ogre/Tat"),
        c("EnSpm_CACTA", "EnSpm/CACTA"),
        c("MuDR_Mutator", "MuDR/Mutator"),
        c("PIF_Harbinger", "PIF/Harbinger"),
        c("Tc1/Mariner", "Tc1/Mariner"),
        c("All|", "")
    )
    for (i in replacement){
        best_hit = gsub(i[1], i[2], best_hit, fixed = TRUE)
    }
    best_hit = gsub("^All", "", best_hit, fixed = FALSE)
    best_hit = ifelse(best_hit == "", paste0("unknown_CL", CL), best_hit)
    output = data.frame(Cluster = CL, classification = best_hit, stringsAsFactors = FALSE)
    return(output)
}
