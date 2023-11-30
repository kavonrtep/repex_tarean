#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
 
CONNECTED = FALSE
if (FALSE) {
    ## for testing
    seqdb = "/mnt/raid/spolecny/petr/RE2/comparative_test/sequences.db"
    hitsortdb = "/mnt/raid/spolecny/petr/RE2/comparative_test/hitsort.db"
    class_file = "/mnt/raid/users/petr/workspace/repex_tarean/databases/classification_tree.rds"
    ## connect to sqlite databases
    SEQDB = dbConnect(RSQLite::SQLite(), seqdb)
    HITSORTDB = dbConnect(RSQLite::SQLite(), hitsortdb)
    CLS_TREE = readRDS(class_file)
}

connect_to_databases = function(seqdb, hitsortdb,classification_hierarchy_file = NULL){
    if (!CONNECTED){
        SEQDB <<- dbConnect(RSQLite::SQLite(), seqdb)
        HITSORTDB <<- dbConnect(RSQLite::SQLite(), hitsortdb)
        if (!is.null(classification_hierarchy_file)){
            CLS_TREE <<- readRDS(classification_hierarchy_file)
        }
        CONNECTED <<- TRUE
    }
}

disconnect_database = function(){
    if (CONNECTED){
        dbDisconnect(SEQDB)
        dbDisconnect(HITSORTDB)
        CONNECTED <<- FALSE
    }
}

nested2named_list = function(x){
    y = as.list(unlist(x[[1]]))
    names(y) = unlist(x[[2]])
    return(y)
}

is_comparative = function(){
    prefix_codes = dbGetQuery(SEQDB,"SELECT * FROM prefix_codes")
    if (nrow(prefix_codes) == 0){
        return(FALSE)
    }else{
        return(TRUE)
    }
}

get_comparative_codes = function(){
    prefix_codes = dbGetQuery(SEQDB,"SELECT * FROM prefix_codes")
    return(prefix_codes)
}

add_preamble = function(html_file, preamble){
  html_content=readLines(html_file)
  modified_html_content = gsub("<body>",
       paste("<body>\n", preamble,"\n"),
       html_content)
  cat(modified_html_content, file = html_file, sep="\n")
}


df2html = function(df, header = NULL, sort_col = NULL, digits = 3, rounding_function=signif, decreasing = TRUE, scroling = FALSE, width = 300){
    if (!is.null(sort_col)){
        df = df[order(df[,sort_col], decreasing = decreasing),]
    }
    if (!is.null(digits)){
        for (i in seq_along(df)){
            if(is.numeric(df[,i])){
                df[,i] = rounding_function(df[,i], digits)
            }
        }
    }
    if (is.null(header)){
        h = ""
    }else{
        h = paste0("    <th>",header,"</th>\n", collapse="") %>%
            paste0(" <tr>\n", .,"  </tr>\n")
    }
    x = apply(df,1,function(x)paste0("    <td>",x,"</td>\n", collapse="")) %>%
        paste0("  <tr>\n", .,"  </tr>\n", collapse = "")
    if (scroling){
        cols = paste0('<col width="',rep(round(100/ncol(df)),ncol(df)),'%">\n',collapse ="")
        height = min(200, 22 * nrow(df))
        out = paste0(
            '<table cellspacing="0" cellpadding="0" border="0" width="',width,'">\n',
            '  <tr>\n',
            '    <td>\n',
            '      <table cellspacing="0" cellpadding="1" border="1" width="', width,'" >\n',
            cols,
            h,
            '      </table>\n',
            '   </td>\n',
            ' </tr>\n',
            ' <tr>\n',
            '   <td>\n',
            '     <div style="width:',width,'px; height:',height,'px; overflow:auto;">\n',
            '       <table cellspacing="0" cellpadding="1" border="1" width="',width,'" >\n',
            cols,
            x,
            '       </table>\n',
            '     </div>\n',
            '  </td>\n',
            ' </tr>\n',
            '</table>\n'
        )

    }else{
        out = paste ("<table>\n", h,x, "</table>\n")
    }
    return(out)
}

start_html = function(filename, header){
    cat(header, file = filename)
    html_writer = function(content, fn=HTML, ...){
        fn(content, append = TRUE, file = filename, ...)
    }
}

preformatted = function(x){
    ## make preformatted html text
    return(
        paste(
        "<pre>\n",
        x,
        "</pre>"
        ,sep="")
    )
}


summary_histogram = function(fn, N_clustering, N_omit=0, size_adjusted=NULL, top_clusters_prop){
    ## assume connection do databases
    communities = dbGetQuery(
        HITSORTDB,
        "SELECT DISTINCT cluster, size, supercluster, supercluster_size FROM communities ORDER BY supercluster" 
    )
    if (N_omit != 0){
      ## adjust communities and cluster sizes:
      cluster_to_adjust = which(
        communities$size[order(communities$cluster)][1:length(size_adjusted)] != size_adjusted
      )
      ## keep original value:
      communities$size_original = communities$size
      superclusters_to_adjust = unique(communities$supercluster[communities$cluster %in% cluster_to_adjust])
      for (cl in cluster_to_adjust){
        communities[communities$cluster == cl,'size'] = size_adjusted[cl]
      }
      for (cl in superclusters_to_adjust){
        communities[communities$supercluster == cl,'supercluster_size'] =
          sum(communities[communities$supercluster == cl,'size'])
      }
    }else{
      cluster_to_adjust=NULL
    }
     singlets = N_clustering - sum(communities$size)

     supercluster_size = sort(unique(communities[, c('supercluster', 'supercluster_size')])$supercluster_size, decreasing =  TRUE) 


   clid2size = sort(communities$size, decreasing = TRUE)

    cluster_id = split(communities$cluster, communities$supercluster)
    cluster_id_sort = lapply(cluster_id, function(x)x[order(clid2size[x], decreasing = FALSE)])

  cluster_size_unsorted = split(communities$size, communities$supercluster)
    cluster_size_sort = lapply(cluster_size_unsorted, function(x) (sort(x)))
    ## reorder by size of superclusters
    cluster_size_sort_sort = cluster_size_sort[order(sapply(cluster_size_sort, sum), decreasing = TRUE)]
    cluster_id_sort_sort = cluster_id_sort[order(sapply(cluster_size_sort, sum), decreasing = TRUE)]

  
    Nmax = max(sapply(cluster_size_sort_sort, length))
    M = cbind(
        sapply(cluster_size_sort_sort, function(x)y = c(x, rep(0, 1 + Nmax - length(x)))),
        c(1, rep(0, Nmax)) 
    )

    Mid =  cbind(
        sapply(cluster_id_sort_sort, function(x)y = c(x, rep(0, 1 + Nmax - length(x)))),
        c(1, rep(0, Nmax)) 
    )

    recolor = matrix(ifelse(Mid %in% cluster_to_adjust,TRUE,FALSE), ncol=ncol(Mid))
    indices = which(recolor, arr.ind = TRUE)

  
    png(fn, width = 1200, height = 700, pointsize = 20)

    plot(0,
         xlim = c(0, sum(c(supercluster_size, singlets))),
         ylim = c(0, max(supercluster_size) * 1.2),
         type = "n", yaxs = 'i', axes = FALSE,
         xlab = "Proportion of reads [%]", ylab = "Number of reads",
         main = paste(N_clustering, "reads total"))

    rect(0, 0,
         sum(supercluster_size),
         max(supercluster_size) * 1.2,
         col = "#0000FF10")

    rect(sum(supercluster_size), 0,
         sum(supercluster_size) + singlets,
         max(supercluster_size) * 1.2,
         col = "#FFAAFF10")

    barplot(M,
            width = c(supercluster_size, singlets),
            space = 0, ylim = c(0, max(supercluster_size) * 1.2),
            col = "#AAAAAA", names.arg = rep("", ncol(M)), add = TRUE,
    )


    for (i in seq_along(indices[,1])){
      y1 = sum(M[1:indices[i,'row'],indices[i,'col']])
      x1 = sum(M[,1:indices[i,'col']])
      if(indices[i,'row'] == 1){
        y0=0
      }else{
        y0 = sum(M[1:(indices[i,'row']-1),indices[i,'col']])
      }
      if (indices[i,'col']==1){
        x0=0
      }else{
        x0 = sum(M[,1:(indices[i,'col']-1)])
      }
      rect(x0,y0,x1,y1, col="#88FF88")
    }
    abline(v=top_clusters_prop * sum(c(supercluster_size, singlets)), col="#00000088", lwd=3, lty=3)

    text(sum(supercluster_size) / 2,
         max(supercluster_size) * 1.05,
         labels = paste0(sum(supercluster_size), " reads in\n",
                      length(supercluster_size), " supeclusters (", nrow(communities), " clusters)")
         )

    
    text(sum(supercluster_size) + singlets / 2,
         max(supercluster_size) * 1.05,
         labels = paste(singlets, "singlets"))
    
    axis(1,at=seq(0,N_clustering,length.out=11),label=seq(0,100,by=10))
  dev.off()
  clustering_info = list(
    Number_of_reads_in_clusters = sum(supercluster_size),
    Number_of_clusters =  nrow(communities),
    Number_of_superclusters = length(supercluster_size),
    Number_of_singlets = singlets
  )
  return(clustering_info)
}


rectMap=function(x,scale.by='row',col=1,xlab="",ylab="",grid=TRUE,axis_pos=c(1,4),cexx=NULL,cexy=NULL){
  if (scale.by=='row'){
                                        #x=(x)/rowSums(x)
    x=(x)/apply(x,1,max)
  }
  if (scale.by=='column'){
    x=t(t(x)/apply(x,2,max))
  }
  nc=ncol(x)
  nr=nrow(x)
  coords=expand.grid(1:nr,1:nc)
  plot(coords[,1],coords[,2],type='n',axes=F,xlim=range(coords[,1])+c(-.5,.5),ylim=range(coords[,2])+c(-.5,.5),xlab=xlab,ylab=ylab)
  axis(axis_pos[1],at=1:nr,labels=rownames(x),lty=0,tick=FALSE,line=0,cex.axis=0.5/log10(nr))
  axis(axis_pos[2],at=1:nc,labels=colnames(x),lty=0,tick=FALSE,las=2,line=0 ,hadj=0, cex.axis=0.7)
  axis(2,at=1:nc,labels=colnames(x),lty=0,tick=FALSE,las=2,line=0 ,hadj=1, cex.axis=0.7)

  mtext(side = 1, "Cluster id", las=1, line = 3, cex = 0.5)
  line = 1.5 + log10(nr)
  mtext(side = 2, "Proportions of individual samples", las =0, line = line, cex = 0.5)
  s=c(x)/2  # to get it proportional
  w = c(x)/2
  rect(coords[,1]-0.5,coords[,2]-s,coords[,1]+0.5,coords[,2]+s,col=col,border=NA)
  if (grid){
    abline(v=0:(nr)+.5,h=0:(nc)+.5,lty=2,col="#60606030")
  }
  box(col="#60606030",lty=2)
}

plot_rect_map = function(read_counts,cluster_annotation, output_file,Xcoef=1,Ycoef=1){
  counts = read.table(read_counts,header=TRUE,as.is=TRUE)
  annot = read.table(cluster_annotation, sep="\t",header=FALSE,as.is=TRUE)
  N = nrow(annot)
  colnames(annot) = c("cluster", "Automatic.classification")
  annot$number.of.reads = rowSums(counts[1 : nrow(annot) ,-1])
  unique_repeats = names(sort(table(c(annot$Automatic.classification,rep('nd',N))),decreasing = TRUE))

  M = as.matrix(counts[1:N,-(1:2)])
  rownames(M) = paste0("CL",rownames(M))
  Mn1=(M)/apply(M,1,max)
  Mn2=M/max(M)
  Mn2=M/apply(M,1,sum)

  ord1 = hclust(dist(Mn1),method = "ward.D")$order
  ord2 = hclust(dist(t(Mn2)))$order
  wdth = (400 + N*10 ) * Xcoef
  hgt = (600 + ncol(M)*50) * Ycoef
  ptsize = round((wdth*hgt)^(1/4))
  png(output_file, width=wdth,height=hgt, pointsize = ptsize)  # was 50
  ploting_area_width = 3 + log10(N)*3
  ploting_area_sides = 1
  layout(matrix(c(4,2,3,4,1,3),ncol=3,byrow = TRUE),
         width=c(ploting_area_sides,ploting_area_width,ploting_area_sides),
         height=c(3,ncol(M)*0.5))
  par(xaxs='i', yaxs = 'i')
  par(las=2,mar=c(4,0,0,0),cex.axis=0.5)
  rectMap(Mn2[ord1,ord2],scale.by='none',col=1, grid=TRUE)
  par(las=2,mar=c(1,0,1,0), mgp = c(2,0.5,0))
  barplot(annot$number.of.reads[ord1], col = 1)
  mtext(side = 2, "Cluster size", las = 3, line = 2, cex = 0.5)
  par(mar=c(0,0,10,0))
  plot.new()
  st = dev.off()
  ## calculate coordinated if boxes to create hyperlink
  X0 = wdth/(ploting_area_sides * 2 + ploting_area_width)* ploting_area_sides
  X1 = wdth/(ploting_area_sides * 2 + ploting_area_width)*(ploting_area_sides + ploting_area_width)
  L = round(seq(X0,X1, length.out = N + 1)[1:N])
  R = round(seq(X0,X1, length.out = N + 1)[2:(N + 1)])
  cn = rownames(Mn2[ord1,ord2])
  cluster_links = paste0(
    "seqclust/clustering/clusters/dir_CL",
    sprintf("%04d", as.integer(substring(cn,3 ))),
    "/index.html")
  coords = paste0(L, ",", 1, ",", R, ",", hgt)
  clustermap = paste0(
    '\n<map name="clustermap"> \n',
    paste0(
      '<area shape="rect"\n      coords="',coords, '"\n',
      '      href="', cluster_links, '"\n',
      '      title="', cn, '"/>\n',
      collapse = ""),
    "</map>\n")
  return(clustermap)
}
