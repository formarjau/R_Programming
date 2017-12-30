###################################
###FUNCTIONS TO WORK WITH METAQC###
###################################

#Functions used in order to performed quality control of expression based studies.

#Common_Set_Genes: This function returns the common genes found in n studies. Since it works with the output files that I generate in my pipeline.
#This files contain a list with an expression metrix as first element and a factor which deffines classes for the samples in the seccond element.

Common_Set_Genes <- function(vec_of_studies){
	for(i in 1:length(vec_of_studies)){
		if(i == 1){
			intersected <- rownames(get(load(file=vec_of_studies[[i]]))[[1]])
		}else{
			intersected <- intersect(intersected,rownames(get(load(file=vec_of_studies[[i]]))[[1]]))
		}
	}
	closeAllConnections()
	return(intersected)
}

#Prepare_Meta_QC: This functions prepares a list of absolute paths to the individual studies to a list which contains expression data and facotrs for each study
#included in the meta-analysis in the format that metaQC requieres in order to perform the quality control.

Prepare_Meta_QC <- function(vec_of_studies,pos_study_name = 4){
	print("JOJOJJ")
	list_out <- list()
	intersected <- Common_Set_Genes(vec_of_studies)
	names_ds <- c()
	for(i in 1:length(vec_of_studies)){
		print(i)
		print(vec_of_studies[i])
		names_ds <- c(names_ds,strsplit(vec_of_studies[i],"\\/")[[1]][pos_study_name])
		print(vec_of_studies[i])
		temp <- get(load(file=vec_of_studies[i]))
		exp_temp <- temp[[1]]
		colnames(exp_temp) <- (as.numeric(temp[[2]])+ 1)
		exp_temp <- exp_temp[intersected,]
		list_out[[i]] <- exp_temp
	}
	closeAllConnections()
	names(list_out) <- names_ds
	return(list_out)
}