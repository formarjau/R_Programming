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

#Plot the module conservation data.

plot_Module_Preservation_Statistics <- function(file,mp){
	pdf(file = file)
	# Module labels and module sizes are also contained in the results
	modColors = rownames(mp$preservation$observed[[ref]][[test]])
	moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
	# leave grey and gold modules out
	plotMods = !(modColors %in% c("grey", "gold"));
	# Text labels for points
	text = modColors[plotMods];
	# Auxiliary convenience variable
	plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
	#Main titles for the plot
	mains = c("Preservation Median rank", "Preservation Zsummary");
	# Start the plot
	sizeGrWindow(10, 5);
	#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
	pdf(file = file, width=12, height=12)
	par(mfrow = c(1,2))
	par(mar = c(4.5,4.5,2.5,1))
	for (p in 1:2)
	{
	min = min(plotData[, p], na.rm = TRUE);
	max = max(plotData[, p], na.rm = TRUE);
	# Adjust ploting ranges appropriately
	if (p==2)
		{
		if (min > -max/10) min = -max/10
		ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
		} else
		ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
		plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
		main = mains[p],
		cex = 2.4,
		ylab = mains[p], xlab = "Module size", log = "x",
		ylim = ylim,
		xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
		labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
	# For Zsummary, add threshold lines
		if (p==2)
		{
		abline(h=0)
		abline(h=2, col = "blue", lty = 2)
		abline(h=10, col = "darkgreen", lty = 2)
	}
}
# If plotting into a file, close it
dev.off();
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


# Perform Reactome Enrichment in all the modules create plots and save results.

 perform_enrichment_in_all_modules <- function(exprs,colors,MEs,directory){
 	library(org.Hs.eg.db)
	library(annotate)
	library(ReactomePA)
	geneModuleMembership = as.data.frame(cor(exprs, MEs, use = "p"));
	modNames = substring(names(MEs), 3);
	names(geneModuleMembership) = paste("MM", modNames, sep="");
	gene_module_assign <- data.frame(colors,colnames(exprs),unlist(lookUp(colnames(exprs), 'org.Hs.eg', 'SYMBOL')),row.names = colnames(exprs));
	colnames(gene_module_assign) <- c("Color","Entrez","Symbol");
	list_modules <- list()
	colors <- unique(gene_module_assign$Color)
	for(i in 1:length(colors)){
		genes <- gene_module_assign[gene_module_assign$Color == colors[i],]
		paths <- enrichPathway(gene=genes$Entrez,pvalueCutoff=0.05, readable=T)
		paths_df <- data.frame(paths) 
		try({
		pdf(file=paste(directory,colors[i],"_dotplot.pdf",sep=""),width = 15,height = 15)
		print(dotplot(paths,showCategory=20))
		dev.off()
		pdf(file=paste(directory,colors[i],"_enrichMap.pdf",sep=""),width = 15,height = 15)
		enrichMap(paths, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.8,n = 40)
		list_temp <- list(genes,paths,paths_df)
		dev.off()})
		print(colors[i])
		list_modules[[colors[i]]] <- list_temp
	}
	return(list_modules)
 }

############################################
####CREATING GENE-WISE ANNOTATION TABLES####
############################################

construct_my_gene_table <- function(my_gene_res,go_type = 5){
  ids <- c()
  terms <- c()
  for(i in 1:length(my_gene_res[[go_type]])){
    ids <- c(ids,my_gene_res[[go_type]][[i]]$id)
    terms <- c(terms,my_gene_res[[go_type]][[i]]$term)
  }
  terms <- unique(terms)
  mat <- matrix(0,nrow  = length(unique(my_gene_res[[1]])),ncol = length(terms))
  df <- data.frame(mat)
  rownames(df) <-  unique(my_gene_res[[1]])
  colnames(df) <- terms
  for(i in 1:length(my_gene_res[[go_type]])){
    temp_row <- my_gene_res[[1]][[i]]
    print(temp_row)
    temp_col <- unique(my_gene_res[[go_type]][[i]]$term)
    print(temp_col)
    df[temp_row,temp_col] <- 1
  }
  return(df)
}

##########################################
##OVERREPRESENTATION ENRICHMENT ANALYSIS##
##########################################

enrich_with_enrichr <- function(list_of_genes){
  library("enrichR")
  query <- enrichr(list_of_genes,databases = c("KEGG_2016","Reactome_2016",))
  return(query)
  }


##############################
#PERFORMING GSEA FUNCTIONS####
##############################

#Step one analysis

perform_Step_One <- function(exprs,directory,tagg,threshold_fit = 0.85,minModuleSize = 30){
	print("Selecting Power...")
	powers = c(c(1:20), seq(from = 12, to=20, by=2))
	sft = pickSoftThreshold(exprs, powerVector = powers, verbose = 5)
	save(file=paste(directory,tagg,"_soft_thresh.Rda"),sft)
	softPower <- sft$fitIndices[sft$fitIndices$SFT.R.sq > threshold_fit,][1,1]
	print(softPower)
	print("Computing adjacency matrix...")
	adjacency = adjacency(exprs, power = softPower);
	save(file=paste(directory,tagg,"_adjacency_Healthy.Rda"),adjacency)
	print("Computing TOM...")
	TOM = TOMsimilarity(adjacency);
	save(file=paste(directory,tagg,"_TOM_Healthy.Rda"),TOM)
	dissTOM = 1-TOM
	geneTree = hclust(as.dist(dissTOM), method = "average");
	minModuleSize = minModuleSize;
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
	dynamicColors = labels2colors(dynamicMods)
	number_of_modules <- length(table(dynamicColors))
	number_of_genes_in_modules <- table(dynamicColors)

	pdf(file = paste(directory,tagg,"_geneTree_Mod_Cols.pdf"))
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
	dev.off()
	print("Computing Module eigengenes...")
	MEList = moduleEigengenes(exprs, colors = dynamicColors)
	MEs = MEList$eigengenes
	MEDiss = 1-cor(MEs);
	METree = hclust(as.dist(MEDiss), method = "average");
	pdf(file = paste(directory,tagg,"_METree.pdf"))
	plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
	MEDissThres = 0.25
	abline(h=MEDissThres, col = "red")
	dev.off()
	save(exprs,softPower,adjacency,TOM,dissTOM,geneTree,minModuleSize,dynamicColors,number_of_modules,number_of_genes_in_modules,MEList,MEs,MEDiss,METree,file= paste(directory,tagg,"Step_1_Data.Rdata",sep=""))
}

#Step two ananlysis

perform_Step_Two <- function(exprs,dynamicColors,tagg,MEDissThres = 0.25,directory){
	print("Merging close modules...")
	merge = mergeCloseModules(exprs, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	pdf(file=paste(directory,tagg,"Module_Gene_Content.pdf"))
	barplot(table(merge$colors),col = names(table(merge$colors)),cex.names = .7,las = 2)
	dev.off()
	mergedColors = merge$colors;
	mergedMEs = merge$newMEs;
	pdf(file =paste(directory,tagg,"_Tree_Merged_Colors.pdf"))
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
	dev.off()
	moduleColors = mergedColors
	colorOrder = c("grey", standardColors(50));
	moduleLabels = match(moduleColors, colorOrder)-1;
	MEs = mergedMEs;
	save(merge,mergedColors,mergedMEs,moduleColors,colorOrder,moduleLabels,MEs,file=paste(directory,tagg,"Step_2_Data.Rdata",sep="\t"))
}

perform_Step_Three <- function(exprs,directory,tagg,dat_traits,moduleColors){
	nGenes = ncol(exprs);
	nSamples = nrow(exprs);
	datTraits = dat_traits
	MEs0 = moduleEigengenes(exprs, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)
	moduleTraitCor = cor(MEs, datTraits, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
	save(nGenes,nSamples,datTraits,MEs0,MEs,moduleTraitCor, moduleTraitPvalue, file = paste(directory,tagg,"Step_3_Data.Rdata"))
	pdf(file = 	paste(directory,tagg,"Module_Trait_Associations.pdf"))
	textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
	signif(moduleTraitPvalue, 1), ")", sep = "");
	dim(textMatrix) = dim(moduleTraitCor)
	par(mar = c(6, 8.5, 3, 3));
	# Display the correlation values within a heatmap plot
	labeledHeatmap(Matrix = moduleTraitCor,
	xLabels = names(datTraits),
	yLabels = names(MEs),
	ySymbols = names(MEs),
	colorLabels = FALSE,
	colors = greenWhiteRed(50),
	textMatrix = textMatrix,
	setStdMargins = FALSE,
	cex.text = 0.5,
	zlim = c(-1,1),
	main = paste("Module-trait relationships"))
	dev.off()
	}