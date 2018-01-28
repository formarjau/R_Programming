#Computes one-tailed fisher test from differential expressoin meta-analysis results.

IntersectionPvals <- function(ResultA,ResultB,FDRA,FDRB){
  ResultA <- ResultA[intersect(rownames(ResultA),rownames(ResultB)),]
  ResultB <- ResultB[intersect(rownames(ResultA),rownames(ResultB)),]
  Common <- length(intersect(rownames(ResultA),rownames(ResultB)))
  ResultAfilt <- ResultA[as.numeric(as.character(ResultA$FDR)) < FDRA,]
  ResultAUp <- ResultAfilt[as.numeric(as.character(ResultAfilt$zval)) > 0 ,]
  ResultADown <- ResultAfilt[as.numeric(as.character(ResultAfilt$zval)) < 0 ,]
  ResultBfilt <- ResultB[as.numeric(as.character(ResultB$FDR))< FDRB,]
  ResultBUp <- ResultBfilt[as.numeric(as.character(ResultBfilt$zval)) > 0 ,]
  print(ResultBUp)
  ResultBDown <- ResultBfilt[as.numeric(as.character(ResultBfilt$zval)) < 0 ,]
  UPUP <- length(intersect(as.character(ResultAUp$symb),as.character(ResultBUp$symb)))
  counts = (matrix(data=c(UPUP,nrow(ResultAUp) - UPUP,nrow(ResultBUp)- UPUP,Common + UPUP - nrow(ResultAUp) - nrow(ResultBUp)),nrow=2))
  testUPUP <- fisher.test(counts,alternative="greater")$p.value
  DOWNDOWN <- length(intersect(as.character(ResultADown$symb),as.character(ResultBDown$symb)))
  counts = (matrix(data=c(DOWNDOWN,nrow(ResultADown) - DOWNDOWN,nrow(ResultBDown)-DOWNDOWN,Common + DOWNDOWN - nrow(ResultADown) - nrow(ResultBDown)),nrow=2))
  testDOWNDOWN <- fisher.test(counts,alternative="greater")$p.value
  UPDOWN <- length(intersect(as.character(ResultAUp$symb),as.character(ResultBDown$symb)))
  counts = (matrix(data=c(UPDOWN,nrow(ResultAUp) - UPDOWN,nrow(ResultBDown) - UPDOWN,Common + UPDOWN - nrow(ResultAUp) - nrow(ResultBDown)),nrow=2))
  testUPDOWN <- fisher.test(counts,alternative="greater")$p.value
  DOWNUP <- length(intersect(as.character(ResultADown$symb),as.character(ResultBUp$symb)))
  counts = (matrix(data=c(DOWNUP,nrow(ResultADown) - DOWNUP,nrow(ResultBUp) - DOWNUP,Common  + DOWNUP - nrow(ResultADown) - nrow(ResultBUp)),nrow=2))
  testDOWNUP <- fisher.test(counts,alternative="greater")$p.value
  listOut <- list(testUPUP,testDOWNDOWN,testUPDOWN,testDOWNUP)
  return(listOut)
}

#Funtion that given a directory with a set of meta-analysis results reads them and prepares thom from the next step.

Create_Results_List <- function(directory){
  list_out <- list()
  names <- c()
  files <- dir(directory,full.names = T)
  for(i in 1:length(files)){
    print(i)
    temp_file <- get(load(file=files[i]))
    temp_file <- transform_results(temp_file)
    names <- c(names,gsub("_.*","",strsplit(files[i],"\\/")[[1]][length(strsplit(files[i],"\\/")[[1]])]))
    list_out[[i]] <- temp_file
  }
  names(list_out) <- names
  return(list_out)
}


#Transform results from the 8 item list format to a data.frame with the entrez IDs transformed into symbols.

transform_results <- function(Meta_Result){
  library(org.Hs.eg.db)
  library(annotate)
  df <- data.frame(Meta_Result[[1]],Meta_Result[[2]],Meta_Result[[3]],Meta_Result[[4]],Meta_Result[[5]],Meta_Result[[6]],Meta_Result[[7]],Meta_Result[[8]])
  colnames(df) <- names(Meta_Result)
  rownames(df) <- names(Meta_Result[[1]])
  df$Entrez <- names(Meta_Result[[1]])
  df$symb <- unlist(lookUp(names(Meta_Result[[1]]),'org.Hs.eg', 'SYMBOL'))
  return(df)
}

Compute_pvalues <- function(list,reference,threshold = 0.05){
  list_out_genes <- list()
  list_out_pvals <- list()
  for(i in 1:length(list)){
    if(i != reference){
      name_out <- paste(names(list)[reference],names(list)[i],sep="_")
      list_out_genes[[name_out]] <- IntersectionUpUpDownDown(list[[reference]],list[[i]],threshold,threshold)
      list_out_pvals[[name_out]] <- IntersectionPvals(list[[reference]],list[[i]],threshold,threshold)
    }
  }
  list_out <- list(list_out_genes,list_out_pvals)
  return(list_out)
}


# Counting number of genes in intersection.

IntersectionUpUpDownDown <- function(ResultA,ResultB,FDRA,FDRB){
  ResultAfilt <- ResultA[as.numeric(as.character(ResultA$FDR)) < FDRA,]
  ResultAUp <- ResultAfilt[as.numeric(as.character(ResultAfilt$zval)) > 0 ,]
  #print(ResultAUp)
  ResultADown <- ResultAfilt[as.numeric(as.character(ResultAfilt$zval)) < 0 ,]
  #print(ResultADown)
  ResultBfilt <- ResultB[as.numeric(as.character(ResultB$FDR))< FDRB,]
  ResultBUp <- ResultBfilt[as.numeric(as.character(ResultBfilt$zval)) > 0 ,]
  #print(ResultBUp)
  ResultBDown <- ResultBfilt[as.numeric(as.character(ResultBfilt$zval)) < 0 ,]
  #print(ResultBDown)
  print(as.character(ResultAUp$symb))
  print(as.character(ResultBUp$symb))
  print(as.character(ResultAUp$symb))
  print(as.character(ResultBDown$symb))
  UPUP <- intersect(as.character(ResultAUp$symb),as.character(ResultBUp$symb))
  DOWNDOWN <- intersect(as.character(ResultADown$symb),as.character(ResultBDown$symb))
  UPDOWN <- intersect(as.character(ResultAUp$symb),as.character(ResultBDown$symb))
  DOWNUP <- intersect(as.character(ResultADown$symb),as.character(ResultBUp$symb))
  listOut <- list(UPUP,DOWNDOWN,UPDOWN,DOWNUP)
  return(listOut)
}

# Transforming pvals in a matrix

pvalsMatrix <- function(ListOfIntersecUpUpDownDown){
  salida <- c()
  for (i in 1:length(ListOfIntersecUpUpDownDown)){
    for (j in 1:length(ListOfIntersecUpUpDownDown[[i]])){
      salida <- c(salida,ListOfIntersecUpUpDownDown[[i]][[j]])
    }
  }
  print(salida)
  salidaMat <- matrix(salida,ncol=4,byrow=TRUE)
  colnames(salidaMat) <- c("upup","downdown","updown","downup")
  rownames(salidaMat) <- names(ListOfIntersecUpUpDownDown)
  return(salidaMat)
}

#Cunting number of genes in each intersection

countingDirections <- function(ListOfIntersecUpUpDownDown){
  salida <- c()
  for (i in 1:length(ListOfIntersecUpUpDownDown)){
    for (j in 1:length(ListOfIntersecUpUpDownDown[[i]])){
      salida <- c(salida,length(ListOfIntersecUpUpDownDown[[i]][[j]]))
    }
  }
  #print(salida)
  salidaMat <- matrix(salida,ncol=4,byrow=TRUE)
  colnames(salidaMat) <- c("upup","downdown","updown","downup")
  rownames(salidaMat) <- names(ListOfIntersecUpUpDownDown)
  return(salidaMat)
}

#Correcting pvalues.

CorrecDfPval <- function(pvalsdf){
  sal <- data.frame(matrix(p.adjust(pvalsdf,method = "fdr"),nrow = length(p.adjust(pvalsdf,method = "fdr"))/4,ncol = 4,byrow = FALSE))
  colnames(sal) <- colnames(pvalsdf)
  rownames(sal) <- rownames(pvalsdf)
  return(sal)
}

#Ordering matrices.

order_matrices <- function(corr_pval_mat,out_gen_mat){
  order_1 <- order(corr_pval_mat[,4])
  corr_pval_mat <- corr_pval_mat[order_1,]   
  out_gen_mat <- out_gen_mat[order_1,]
  order_2 <- order(corr_pval_mat[,3])
  corr_pval_mat <- corr_pval_mat[order_2,]   
  out_gen_mat <- out_gen_mat[order_2,]
  order_3 <- order(corr_pval_mat[,2])
  corr_pval_mat <- corr_pval_mat[order_3,]   
  out_gen_mat <- out_gen_mat[order_3,]
  order_4 <- order(corr_pval_mat[,1])
  corr_pval_mat <- corr_pval_mat[order_4,]   
  out_gen_mat <- out_gen_mat[order_4,] 
  list_sal <- list(corr_pval_mat,out_gen_mat)
  return(list_sal)
}

#Plotting tables from previous results.

plot_results <- function(list_of_result_tables,file){
  library(gplots)
  text <- paste(as.matrix(list_of_result_tables[[2]])," (",as.matrix(format(list_of_result_tables[[1]], digits = 3, scientific = 5)),")",sep="")
  text <- matrix(text,ncol = 4)
  colnames(text) <- colnames(list_of_result_tables[[2]])
  rownames(text) <- rownames(list_of_result_tables[[2]])
  my_palette <- colorRampPalette(c("coral2","white","darkolivegreen3"))(n = 59) #Set the color 
  col_breaks = c(seq(max(-log(as.matrix(list_of_result_tables[[1]]))),4.6,length=20),seq(4.59,2.99,length=20),seq(2.98,0.01,length=20))
  col_breaks <- rev(col_breaks)
  pdf(file = file)
  heatmap.2(-log(as.matrix(list_of_result_tables[[1]])),notecol="black",cellnote = text ,trace="none",breaks=col_breaks,col=my_palette, dendrogram = "none", Rowv = FALSE, Colv = FALSE,labRow = rownames(list_of_result_tables[[1]]) ,labCol = colnames(list_of_result_tables[[1]]),key = FALSE,lwid=c(0.05,0.2), lhei=c(0.1,1),cexRow=1,cexCol = 1.3,margins = c(7,11),srtCol=45)
  dev.off()
}


#Generate pathway associated heatmaps.
