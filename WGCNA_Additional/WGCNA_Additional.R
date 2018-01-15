###########################
####PLOTTING FUNCTIONS#####
###########################

#Plot heatmap ordering by PER3 expression levels and by correlation values.

plot_Heat_Ord_By_Gene <- function(corr_results,exprs_results){
  library(circlize)
  library(org.Hs.eg.db)
  library(annotate)
  library(ComplexHeatmap)
  exprs_results <- t(exprs_results)
  exprs_results <- t(scale(t(exprs_results)))
  ordered <- corr_results[order(corr_results[,3],decreasing = T),2]
  exprs_results <- exprs_results[ordered,]
  exprs_results <- exprs_results[,order(exprs_results[1,])]
  #exprs_results <- t(scale(t(exprs_results)))
  rownames(exprs_results) <- paste(unlist(lookUp(rownames(exprs_results), 'org.Hs.eg', 'SYMBOL')),sep=" ")
  Heatmap(exprs_results, col = colorRamp2(c(min(exprs_results),mean(exprs_results),max(exprs_results)), c("blue","white","red")),cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = TRUE, show_column_names = FALSE,row_names_gp = gpar(fontsize = 5))
}

#,"Cor:",round(corr_results[ordered,3],digits = 2)

#Funcions useful when working with WGCNA to generate cleaner code.

plot_hard_tresh <- function(hard_tresh_ob,file){
  cex1=0.7
  thresholds1 = hard_tresh_ob[,1]
  pdf(file = file,width = 15)
  par(mfrow=c(1,2))
  plot(RdichotTable[,1], -sign(RdichotTable[,4])*RdichotTable[,3],xlab="Hard
       Threshold tau",ylab="Scale Free Topology Model Fit,signed R^2", type="n")
  text(RdichotTable[,1], -sign(RdichotTable[,4])*RdichotTable[,3] ,
       labels=thresholds1,cex=cex1)
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.85,col="red")
  plot(RdichotTable[,1], RdichotTable[,6],xlab="Hard Threshold tau",ylab="Mean
       Connectivity", type="n")
  text(RdichotTable[,1], RdichotTable[,6] , labels=thresholds1, cex=cex1) 
  dev.off()
}

plot_soft_thresh <- function(sft,file){
  
  # Plot the results:
  pdf(file = file,width = 10,height = 8)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
}

#Plot TOM clusters 

plot_TOM_Clusters <- function(){
  
}


##############################################
###FIND NEIGHBOURS IN A CORRELATION MATRIX####
##############################################

#Function that given a correlation matrix gets the direct neigbours of the selected genes

get_direct_neighbours <- function(cor_table,gene,threshold,translate = T){
  library(org.Hs.eg.db)
  library(annotate)
  bool <- abs(cor_table[gene,]) > threshold
  #rep(gene,length(colnames(cor_table)[bool]))
  df_sal <- data.frame(rep(gene,length(colnames(cor_table)[bool])),colnames(cor_table)[bool],cor_table[gene,bool])
  colnames(df_sal) <- c("GENE_A","GENE_B","CORRELATION")
  if(translate == T){
    df_sal$GENE_A_SYMB <- unlist(lookUp(df_sal[,1], 'org.Hs.eg', 'SYMBOL'))
    df_sal$GENE_B_SYMB <- unlist(lookUp(df_sal[,2], 'org.Hs.eg', 'SYMBOL'))
  }
  return(df_sal)
}

get_direct_neighbours_and_first <- function(cor_table,gene,threshold){
  library(org.Hs.eg.db)
  library(annotate)
  bool <- abs(cor_table[gene,]) > threshold
  df_sal <- data.frame(rep(gene,length(colnames(cor_table)[bool])),colnames(cor_table)[bool],cor_table[gene,bool])
  colnames(df_sal) <- c("GENE_A","GENE_B","CORRELATION")
  for(i in 1:nrow(df_sal)){
    df_temp <- get_direct_neighbours(cor_table,df_sal[i,2],threshold,translate = F)
    df_sal <- rbind(df_sal,df_temp) 
  }
  df_sal_filt <- data.frame(t(c(NA,NA,NA)))
  colnames(df_sal_filt) <- c("GENE_A","GENE_B","CORRELATION")
  for(i in 1:nrow(df_sal)){
    if(!(paste(df_sal[i,1],df_sal[i,2],sep="_") %in% paste(df_sal_filt[,1],df_sal_filt[,2],sep="_")) & !(paste(df_sal[i,1],df_sal[i,2],sep="_") %in% paste(df_sal_filt[,2],df_sal_filt[,1],sep="_"))){
      df_sal_filt <- rbind(df_sal_filt,df_sal[i,])
      colnames(df_sal_filt) <- c("GENE_A","GENE_B","CORRELATION")
    }
  }
  df_sal_filt <- df_sal_filt[-1,]
  df_sal_filt$GENE_A_SYMB <- unlist(lookUp(df_sal_filt[,1], 'org.Hs.eg', 'SYMBOL'))
  df_sal_filt$GENE_B_SYMB <- unlist(lookUp(df_sal_filt[,2], 'org.Hs.eg', 'SYMBOL'))
  return(df_sal_filt)
}

#################################
####FAST ENRICHMENT FUNCTIONS####
#################################

enrich_with_fgsea <- function(corr_list,pathway){
  library(GSEABase)
  library(fgsea)
  paths <- gmtPathways(pathway)
  fgseaRes <- fgsea(pathways = paths, stats = corr_list,minSize=15,maxSize=10000,nperm=10000)
  return(fgseaRes)
  }
