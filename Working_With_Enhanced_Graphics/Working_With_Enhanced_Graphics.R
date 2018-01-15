#Given a list of enrichments compute a network and attribute table.

#Given list of enrichment directoris. Give one as reference. Give a threshold. Read them. Leve only those associated with reference. Modify names.

#Getting dir of a particular report.

get_GSEA_result_in_dir <- function(dir,pos_or_neg = "POS"){
	files <- dir(dir,full.names = T)
	if(pos_or_neg == "POS"){
		files <- files[grepl("report",files) & grepl("pos",files) & grepl("xls",files)]
		#print(files)
	}else{
		files <- files[grepl("report",files) & grepl("neg",files) & grepl("xls",files)]
		#print(files)
	}
	return(files)
}

#Reading GSEA results

read_GSEA_results <- function(file,FDR_threshold){
	res <- read.csv(file,sep="\t")
	res <- res[res$FDR.q.val < FDR_threshold,]
	return(res)
}

#Generating table

generating_Table <- function(directories,reference,threshold){
	temp_pos_ref <- get_GSEA_result_in_dir(directories[reference],"POS")
	temp_pos_ref <- read_GSEA_results(temp_pos_ref,0.05)
	temp_neg_ref <- get_GSEA_result_in_dir(directories[reference],"NEG")
	temp_neg_ref <- read_GSEA_results(temp_neg_ref,0.05)
	temp_all_ref <- rbind(temp_pos_ref,temp_neg_ref)
	refpaths <- temp_all_ref[,1]
	for(i in 1:length(directories)){
		disease <- gsub(".*\\/","",directories[i])
		if(i == 1){
			temp_pos <- get_GSEA_result_in_dir(directories[i],"POS")
			temp_pos <- read_GSEA_results(temp_pos,0.05)
			temp_pos <- temp_pos[,-1]
			temp_pos$disease <- rep(disease,nrow(temp_pos))
			temp_pos$dir <- rep("UP",nrow(temp_pos))
			temp_neg <- get_GSEA_result_in_dir(directories[i],"NEG")
			temp_neg <- read_GSEA_results(temp_neg,0.05)
			temp_neg <- temp_neg[,-1]
			temp_neg$dir <- rep("DOWN",nrow(temp_neg))
			temp_neg$disease <- rep(disease,nrow(temp_neg))
			temp_all <- rbind(temp_pos,temp_neg)
		}	
		else{
			temp_pos <- get_GSEA_result_in_dir(directories[i],"POS")
			temp_pos <- read_GSEA_results(temp_pos,0.05)
			temp_pos <- temp_pos[,-1]
			temp_pos$disease <- rep(disease,nrow(temp_pos))
			temp_pos$dir <- rep("UP",nrow(temp_pos))
			temp_neg <- get_GSEA_result_in_dir(directories[i],"NEG")
			temp_neg <- read_GSEA_results(temp_neg,0.05)
			temp_neg <- temp_neg[,-1]
			temp_neg$disease <- rep(disease,nrow(temp_neg))
			temp_neg$dir <- rep("DOWN",nrow(temp_neg))
			temp_all <- rbind(temp_all,temp_pos,temp_neg)
		}
	}
	temp_all <- temp_all[temp_all[,1] %in% refpaths,]
	return(temp_all)
}

filter_Table <- function(table_data,group_A,group_B,group_A_name,group_B_name,reference,reference_name){
	diseases <- table_data$disease
	group <- c()
	for(i in 1:length(diseases)){
		if(diseases[i] %in% group_A){
			group <- c(group,group_A_name)
		}else if(diseases[i] %in% group_B){
			group <- c(group,group_B_name)
		}else if(diseases[i] %in% reference){
			group <- c(group,reference_name)
		}else{
			group <- c(group,"NO_GROUP")
		}
	}
	table_data$group <- group
	return(table_data)
}

generate_Pie_Chart <- function(list_of_lists,df_conversion,diseases_att){
	path <- c()
	pie_command <- c()
	for(i in 1:length(out_c)){
  	path <- c(path,names(out_c)[i])
  	pie_command <- c(pie_command,paste("piechart: attributelist=",'"',paste(ordered(list_of_lists[[i]]$disease),collapse = ","),'"'," colorlist=",'"',paste(df_conversion[as.character(ordered(list_of_lists[[i]]$disease)),2],collapse = ","),'"',sep=""))
  	print(names(out_c)[i])
  	print(paste("piechart: valueslist=",'"',paste(ordered(list_of_lists[[i]]$disease),collapse = ","),'"'," colorlist=",'"',paste(df_conversion[as.character(ordered(list_of_lists[[i]]$disease)),2],collapse = ","),'"',sep=""))
		}
	out <- data.frame(path,pie_command)
	colnames(out) <- c("path","pie_command")
	df <- data.frame(matrix(0,nrow = nrow(out),ncol = length(diseases_att)))
	colnames(df) <- diseases_att
	for(i in 1:length(diseases_att)){
		for(j in 1:nrow(out)){
			if(grepl(diseases_att[i],out[j,2])){
				df[j,i] <- 1
			}
		}
	}
	out <- cbind(out,df)
	return(out)
}