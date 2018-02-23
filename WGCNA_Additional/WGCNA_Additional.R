################################
####FUNCTIONS TO IMPORT DATA####
################################

retrieve_expression_matrix <- function(list_dirs){
  for(i in 1:length(list_dirs)){
    print(i)
    exp_file_temp <- dir(list_dirs[i],full.names = T)[grepl("All",dir(list_dirs[i],full.names = T))]
    if(i == 1 ){
      df_exp <- get(load(file = exp_file_temp))[[1]]
      print(dim(df_exp))
      intersected <- rownames(df_exp)
    }else{
      df_exp_temp <- get(load(file = exp_file_temp))[[1]]
      print(dim(df_exp_temp))
      intersected <- intersect(intersected,rownames(df_exp_temp))
      df_exp <- cbind(df_exp[intersected,],df_exp_temp[intersected,])
    }
  }
  return(df_exp)
}

#Create df_exp by scaling each gene within each dataset.

retrieve_expression_matrix_scaling <- function(list_dirs){
  for(i in 1:length(list_dirs)){
    print(i)
    exp_file_temp <- dir(list_dirs[i],full.names = T)[grepl("All",dir(list_dirs[i],full.names = T))]
    if(i == 1 ){
      df_exp <- get(load(file = exp_file_temp))[[1]]
      print(dim(df_exp))
      intersected <- rownames(df_exp)
      df_exp <- t(scale(t(df_exp)))
    }else{
      df_exp_temp <- get(load(file = exp_file_temp))[[1]]
      df_exp_temp <- t(scale(t(df_exp_temp)))
      intersected <- intersect(intersected,rownames(df_exp_temp))
      df_exp <- cbind(df_exp[intersected,],df_exp_temp[intersected,])
    }
  }
  return(df_exp)
}


library(plyr)
library(GEOquery)
retrieve_pheno_matrix <- function(list_dirs){
  for(i in 1:length(list_dirs)){
    print(i)
    pheno_file_temp <- dir(list_dirs[i],full.names = T)[grepl("Eset",dir(list_dirs[i],full.names = T))]
    print(pheno_file_temp)
    if(i == 1 ){
      df_pheno <- pData(get(load(file = pheno_file_temp)))
      df_pheno$pCh_Sample_Names <- rownames(df_pheno)
    }else{
      df_pheno_temp <- pData(get(load(file = pheno_file_temp)))
      df_pheno_temp$pCh_Sample_Names <- rownames(df_pheno_temp)
      df_pheno <- rbind.fill(df_pheno,df_pheno_temp)
    }
  }
  rownames(df_pheno) <- df_pheno$pCh_Sample_Names
  return(df_pheno)
}


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

 perform_enrichment_in_all_modules <- function(exprs,colors,MEs,directory,n_dotplot = 20,n_enrichMap = 20){
 	library(org.Hs.eg.db)
	library(annotate)
	library(ReactomePA)
	geneModuleMembership = as.data.frame(cor(exprs, MEs, use = "p"));
	#print(head(geneModuleMembership))
	modNames = substring(names(MEs), 3);
	names(geneModuleMembership) = paste("MM", modNames, sep="");
	gene_module_assign <- data.frame(colors,colnames(exprs),unlist(lookUp(colnames(exprs), 'org.Hs.eg', 'SYMBOL')),row.names = colnames(exprs));
	colnames(gene_module_assign) <- c("Color","Entrez","Symbol");
	list_modules <- list()
	colors <- unique(gene_module_assign$Color)
	for(i in 1:length(colors)){
		genes <- gene_module_assign[gene_module_assign$Color == colors[i],]
		paths <- enrichPathway(gene=genes$Entrez,pvalueCutoff=0.05, readable=T)
		try({
		paths_df <- as.data.frame(paths)
		if(nrow(paths_df) > 0){
		write.table(file = paste(directory,colors[i],"Enrichment.csv",sep=""),paths_df,quote = FALSE,sep="\t")
		pdf(file=paste(directory,colors[i],"_dotplot.pdf",sep=""),width = 15,height = 15)
		print(dotplot(paths,showCategory=n_dotplot))
		dev.off()
		pdf(file=paste(directory,colors[i],"_enrichMap.pdf",sep=""),width = 15,height = 15)
		enrichMap(paths, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.8,n = n_enrichMap)
		list_temp <- list(genes,paths,paths_df)
		dev.off()
		print(colors[i])
		list_modules[[colors[i]]] <- list_temp
			}
		})
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

#Transform pheno to numeric

transform_Pheno <- function(pheno){
	for(i in 1:ncol(pheno)){
		pheno[,i] <- as.numeric(pheno[,i])
	}
}

#Step one analysis

perform_Step_One <- function(exprs,directory,tagg,threshold_fit = 0.85,minModuleSize = 30){
	print("Selecting Power...")
	powers = c(c(1:20), seq(from = 12, to=20, by=2))
	sft = pickSoftThreshold(exprs, powerVector = powers, verbose = 5)
	save(file=paste(directory,tagg,"_soft_thresh.Rda",sep=""),sft)
	softPower <- sft$fitIndices[sft$fitIndices$SFT.R.sq > threshold_fit,][1,1]
	print(softPower)
	print("Computing adjacency matrix...")
	adjacency = adjacency(exprs, power = softPower);
	save(file=paste(directory,tagg,"_adjacency_Healthy.Rda",sep=""),adjacency)
	print("Computing TOM...")
	TOM = TOMsimilarity(adjacency);
	save(file=paste(directory,tagg,"_TOM_Healthy.Rda",sep=""),TOM)
	dissTOM = 1-TOM
	geneTree = hclust(as.dist(dissTOM), method = "average");
	minModuleSize = minModuleSize;
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
	dynamicColors = labels2colors(dynamicMods)
	number_of_modules <- length(table(dynamicColors))
	number_of_genes_in_modules <- table(dynamicColors)

	pdf(file = paste(directory,tagg,"_geneTree_Mod_Cols.pdf",sep=""))
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
	dev.off()
	print("Computing Module eigengenes...")
	MEList = moduleEigengenes(exprs, colors = dynamicColors)
	MEs = MEList$eigengenes
	MEDiss = 1-cor(MEs);
	METree = hclust(as.dist(MEDiss), method = "average");
	pdf(file = paste(directory,tagg,"_METree.pdf",sep=""))
	plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
	MEDissThres = 0.25
	abline(h=MEDissThres, col = "red")
	dev.off()
	save(exprs,softPower,adjacency,TOM,dissTOM,geneTree,minModuleSize,dynamicColors,number_of_modules,number_of_genes_in_modules,MEList,MEs,MEDiss,METree,file= paste(directory,tagg,"_Step_1_Data.Rdata",sep=""))
}

#Step two ananlysis

perform_Step_Two <- function(exprs,dynamicColors,geneTree,tagg,MEDissThres = 0.25,directory){
	print("Merging close modules...")
	merge = mergeCloseModules(exprs, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	pdf(file=paste(directory,tagg,"_Module_Gene_Content.pdf",sep=""))
	barplot(table(merge$colors),col = names(table(merge$colors)),cex.names = .7,las = 2)
	dev.off()
	mergedColors = merge$colors;
	mergedMEs = merge$newMEs;
	pdf(file =paste(directory,tagg,"_Tree_Merged_Colors.pdf",sep=""))
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
	dev.off()
	moduleColors = mergedColors
	colorOrder = c("grey", standardColors(50));
	moduleLabels = match(moduleColors, colorOrder)-1;
	MEs = mergedMEs;
	save(merge,mergedColors,mergedMEs,moduleColors,colorOrder,moduleLabels,MEs,file=paste(directory,tagg,"_Step_2_Data.Rdata",sep=""))
}

#Step three analysis

perform_Step_Three <- function(exprs,directory,tagg,dat_traits,moduleColors){
	nGenes = ncol(exprs);
	nSamples = nrow(exprs);
	datTraits = dat_traits
	MEs0 = moduleEigengenes(exprs, moduleColors)$eigengenes
	MEs = orderMEs(MEs0)
	moduleTraitCor = cor(MEs, datTraits, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
	save(nGenes,nSamples,datTraits,MEs0,MEs,moduleTraitCor, moduleTraitPvalue, file = paste(directory,tagg,"_Step_3_Data.Rdata",sep=""))
	pdf(file = 	paste(directory,tagg,"Module_Trait_Associations.pdf",sep=""))
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

#Step four analysis

perform_Step_Four <- function(datExpr,moduleColors,tagg,directory){
	MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
	MEs_Ord = orderMEs(MEs)
	library(corrplot)
	cor_to_plot <- cor(MEs_Ord)
	colnames(cor_to_plot) <- gsub("ME","",colnames(cor_to_plot))
	rownames(cor_to_plot) <- gsub("ME","",rownames(cor_to_plot))
	pdf(file=paste(directory,tagg,"Eigengene_Correlations.pdf",sep=""))
	corrplot(cor_to_plot, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 90,tl.cex = 0.7)
	dev.off()
	pdf(file=paste(directory,tagg,"Eigengene_Associations.pdf",sep=""))
	plotEigengeneNetworks(MEs_Ord, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
	plotDendrograms = TRUE, xLabelsAngle = 90,plotHeatmaps = T)
	save(cor_to_plot, file = paste(directory,tagg,"_Step_4_Data.Rdata",sep=""))
	dev.off()
}


create_Gene_Module_Memberships <- function(exprs,moduleColors,genes,tagg,directory,dat_traits){
	library(org.Hs.eg.db)
	library(annotate)
	nGenes = ncol(exprs);
	nSamples = nrow(exprs);
	MEs0 = moduleEigengenes(exprs, moduleColors)$eigengenes
	colnames(MEs0)
	#Cual es la diferencia entre ME y MEs. 
	MEs = orderMEs(MEs0)
	colnames(MEs)
	modNames = substring(names(MEs), 3)
	geneModuleMembership = as.data.frame(cor(exprs, MEs, use = "p"))
	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
	names(geneModuleMembership) = paste("MM", modNames, sep="");
	names(MMPvalue) = paste("p.MM", modNames, sep="");
	#print(head(geneModuleMembership))
	list_of_gene_trait_significance <- list()
	for(i in 1:ncol(dat_traits)){
		print(table(is.na(as.numeric(dat_traits[,i]))))
		geneTraitSignificance_temp = as.data.frame(cor(exprs, as.numeric(dat_traits[,i]), use = "p"));
		GSPvalue_temp = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_temp), nSamples));
		names(geneTraitSignificance_temp) = paste("GS.", colnames(dat_traits)[i], sep="");
		names(GSPvalue_temp) = paste("p.GS.", colnames(dat_traits)[i], sep="");
		temp_list <- list(geneTraitSignificance_temp,GSPvalue_temp)
		list_of_gene_trait_significance[[i]] <- temp_list
	}
	names(list_of_gene_trait_significance) <- colnames(dat_traits)
	gene_module_assign <- data.frame(moduleColors,colnames(exprs),unlist(lookUp(colnames(exprs), 'org.Hs.eg', 'SYMBOL')),row.names = colnames(exprs))
	colnames(gene_module_assign) <- c("Color","Entrez","Symbol")
	
	for(k in 1:length(genes)){
	directory_temp <- paste(directory,tagg,"_",genes[k],"/",sep="")
  if(!dir.exists(directory_temp)){
			  dir.create(directory_temp)
			}
  pdf(file = paste(directory_temp,"gene_Module_Membership.pdf",sep=""))
		barplot(as.numeric(geneModuleMembership[genes[k],]),col = gsub("MM","",names(geneModuleMembership[genes[k],])),cex.names = .7,las = 2)
		dev.off()
	for(i in 1:ncol(geneModuleMembership)){
	  for(j in 1:length(list_of_gene_trait_significance)){
	    module_temp <- colnames(geneModuleMembership)[i]
	    trait_temp <- names(list_of_gene_trait_significance)[j]
	    pdf(file = paste(directory_temp,colnames(geneModuleMembership)[i],"_Vs_",names(list_of_gene_trait_significance)[j],"_",genes[k],"_","gene_MM_Vs_TS.pdf",sep=""))
	    plot(geneModuleMembership[,i],list_of_gene_trait_significance[[j]][[1]][,1],col = ifelse(rownames(geneModuleMembership) == genes[k],"red","black"),cex =  ifelse(rownames(geneModuleMembership) == genes[k],1,1),xlim = c(-1,1),ylim = c(-1,1),ylab = trait_temp,xlab = module_temp)
	    points(geneModuleMembership[genes[k],i],list_of_gene_trait_significance[[j]][[1]][genes[k],1],col="red",cex = 3,pch = 20)
	    text(geneModuleMembership[genes[k],i],list_of_gene_trait_significance[[j]][[1]][genes[k],1]+0.1,labels = genes[k],col="red",cex = 2)
	    dev.off()
	  }
	}
}
save(geneModuleMembership,MMPvalue,list_of_gene_trait_significance,gene_module_assign,file = paste(directory_temp,"/",tagg,"_Single_Gene_Analysis_Data.Rdata",sep=""))
}


######################
###AD HOC FUNCTIONS###
######################

curate_breast_pheno_data <- function(exprs,Pheno){
	intersected <- intersect(colnames(exprs),rownames(Pheno))
	exprs <- exprs[,intersected]
	Pheno <- Pheno[intersected,]
	selected <- c("pCh_Status","pCh_Individual","pCh_Study","pCh_Subtype_PAM50","WGCNA_Basal","WGCNA_Her2","WGCNA_LumA","WGCNA_LumB","WGCNA_Normal","WGCNA_Others","pCh_ER","pCh_PR","pCh_HER2","pCh_Grade","pCh_Age")
	Pheno <- Pheno[,selected]
	pCh_ER <- Pheno$pCh_ER
	pCh_ER[grepl("\\+",pCh_ER)] <- "POS"
	pCh_ER[grepl("\\-",pCh_ER)] <- "NEG"
	pCh_ER[grepl("0",pCh_ER)] <- "NEG"
	pCh_ER[grepl("1",pCh_ER)] <- "POS"
	pCh_ER[grepl("ER_Negative",pCh_ER)] <- "NEG"
	pCh_ER[grepl("ER_Positive",pCh_ER)] <- "POS"
	pCh_ER[grepl("negative",pCh_ER)] <- "NEG"
	pCh_ER[grepl("positive",pCh_ER)] <- "POS"
	pCh_ER[grepl("pos",pCh_ER)] <- "POS"
	pCh_ER[grepl("neg",pCh_ER)] <- "NEG"
	pCh_ER[grepl("Positive",pCh_ER)] <- "POS"
	pCh_ER[grepl("Negative",pCh_ER)] <- "NEG"
	pCh_ER[grepl("3",pCh_ER)] <- NA
	pCh_ER[grepl("4",pCh_ER)] <- NA
	pCh_ER[grepl("5",pCh_ER)] <- NA
	pCh_ER[grepl("6",pCh_ER)] <- NA
	pCh_ER[grepl("7",pCh_ER)] <- NA
	pCh_ER[grepl("8",pCh_ER)] <- NA
	pCh_ER[grepl("N/A",pCh_ER)] <- NA
	pCh_ER[grepl("NA",pCh_ER)] <- NA
	pCh_ER[grepl("UNDET",pCh_ER)] <- NA
	pCh_ER[grepl("not_analyzed",pCh_ER)] <- NA
	Pheno$pCh_ER <- pCh_ER
	pCh_PR <- Pheno$pCh_PR
	pCh_PR[grepl("\\+",pCh_PR)] <- "POS"
	pCh_PR[grepl("\\-",pCh_PR)] <- "NEG"
	pCh_PR[grepl("0",pCh_PR)] <- "NEG"
	pCh_PR[grepl("1",pCh_PR)] <- "POS"
	pCh_PR[grepl("PR_Negative",pCh_PR)] <- "NEG"
	pCh_PR[grepl("PR_Positive",pCh_PR)] <- "POS"
	pCh_PR[grepl("negative",pCh_PR)] <- "NEG"
	pCh_PR[grepl("positive",pCh_PR)] <- "POS"
	pCh_PR[grepl("pos",pCh_PR)] <- "POS"
	pCh_PR[grepl("neg",pCh_PR)] <- "NEG"
	pCh_PR[grepl("Positive",pCh_PR)] <- "POS"
	pCh_PR[grepl("Negative",pCh_PR)] <- "NEG"
	pCh_PR[grepl("2",pCh_PR)] <- NA
	pCh_PR[grepl("3",pCh_PR)] <- NA
	pCh_PR[grepl("4",pCh_PR)] <- NA
	pCh_PR[grepl("5",pCh_PR)] <- NA
	pCh_PR[grepl("6",pCh_PR)] <- NA
	pCh_PR[grepl("7",pCh_PR)] <- NA
	pCh_PR[grepl("8",pCh_PR)] <- NA
	pCh_PR[grepl("N/A",pCh_PR)] <- NA
	pCh_PR[grepl("NA",pCh_PR)] <- NA
	pCh_PR[grepl("UNDET",pCh_PR)] <- NA
	pCh_PR[grepl("not_analyzed",pCh_PR)] <- NA
	Pheno$pCh_PR <- pCh_PR
	pCh_HER2 <- Pheno$pCh_HER2
	pCh_HER2[grepl("\\+",pCh_HER2)] <- "POS"
	pCh_HER2[grepl("\\-",pCh_HER2)] <- "NEG"
	pCh_HER2[grepl("0",pCh_HER2)] <- "NEG"
	pCh_HER2[grepl("1",pCh_HER2)] <- "POS"
	pCh_HER2[grepl("HER2_Negative",pCh_HER2)] <- "NEG"
	pCh_HER2[grepl("HER2_Positive",pCh_HER2)] <- "POS"
	pCh_HER2[grepl("negative",pCh_HER2)] <- "NEG"
	pCh_HER2[grepl("positive",pCh_HER2)] <- "POS"
	pCh_HER2[grepl("pos",pCh_HER2)] <- "POS"
	pCh_HER2[grepl("neg",pCh_HER2)] <- "NEG"
	pCh_HER2[grepl("Positive",pCh_HER2)] <- "POS"
	pCh_HER2[grepl("Negative",pCh_HER2)] <- "NEG"
	pCh_HER2[grepl("2",pCh_HER2)] <- NA
	pCh_HER2[grepl("3",pCh_HER2)] <- NA
	pCh_HER2[grepl("4",pCh_HER2)] <- NA
	pCh_HER2[grepl("5",pCh_HER2)] <- NA
	pCh_HER2[grepl("6",pCh_HER2)] <- NA
	pCh_HER2[grepl("7",pCh_HER2)] <- NA
	pCh_HER2[grepl("8",pCh_HER2)] <- NA
	pCh_HER2[grepl("unk",pCh_HER2)] <- NA
	pCh_HER2[grepl("normal",pCh_HER2)] <- NA
	pCh_HER2[grepl("N/A",pCh_HER2)] <- NA
	pCh_HER2[grepl("NA",pCh_HER2)] <- NA
	pCh_HER2[grepl("UNDET",pCh_HER2)] <- NA
	pCh_HER2[grepl("not_analyzed",pCh_HER2)] <- NA
	Pheno$pCh_HER2 <- pCh_HER2
	WGCNA_Healthy <- ifelse(Pheno$pCh_Status == "NT",1,0)
	Pheno$WGCNA_Healthy <- WGCNA_Healthy
	out <- list(exprs,Pheno)
	return(out)
}