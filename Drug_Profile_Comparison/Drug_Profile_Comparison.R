
compare_drug_profiles <- function(drug_prof_A,drug_prof_B,threshold = 0.05,disease_Names = c("Dis_A","Dis_B")){
	drug_prof_A <- drug_prof_A[drug_prof_A$padj < threshold,]
	mim_or_rev <- ifelse(drug_prof_A$NES > 0,"MIMIC","REVERSE")
	drug_prof_A$mim_or_rev <- mim_or_rev
	colnames(drug_prof_A) <- paste("Dis_A_",colnames(drug_prof_A),sep="")
	drug_prof_B <- drug_prof_B[drug_prof_B$padj < threshold,]
	mim_or_rev <- ifelse(drug_prof_B$NES > 0,"MIMIC","REVERSE")
	drug_prof_B$mim_or_rev <- mim_or_rev
	colnames(drug_prof_B) <- paste("Dis_B_",colnames(drug_prof_B),sep="")
	merged <- merge(drug_prof_A,drug_prof_B,by.x = "Dis_A_pathway", by.y = "Dis_B_pathway")
	merged <- rbind(merged[merged$Dis_A_mim_or_rev == "MIMIC" & merged$Dis_B_mim_or_rev == "MIMIC",],merged[merged$Dis_A_mim_or_rev == "REVERSE" & merged$Dis_B_mim_or_rev == "REVERSE",],merged[merged$Dis_A_mim_or_rev == "MIMIC" & merged$Dis_B_mim_or_rev == "REVERSE",],merged[merged$Dis_A_mim_or_rev == "REVERSE" & merged$Dis_B_mim_or_rev == "MIMIC",])
	merged$Dis_A <- rep(disease_Names[1],nrow(merged))
	merged$Dis_B <- rep(disease_Names[2],nrow(merged))
	return(merged)
}

filter_comparisons <- function(result_data){
	result_data <- result_data[,c("Dis_A_pathway","Dis_A_mim_or_rev","Dis_B_mim_or_rev","Dis_A","Dis_B")]
	return(result_data)
}

contruct_net_given_folder <- function(vec_dirs,vec_names,reference,threshold){
	ref_dis <- get(load(file=vec_dirs[reference]))
	ref_dis_name <- vec_names[reference]
	first <- TRUE
	for(i in 1:length(vec_dirs)){
		if(i != reference){
			temp_dis <-  get(load(file=vec_dirs[i]))
			temp_dis_name <- vec_names[i]
			if(first){
				merged <- compare_drug_profiles(ref_dis,temp_dis,threshold = 0.01,c(ref_dis_name,temp_dis_name))
				merged <- filter_comparisons(merged)
				first <- FALSE
			}else{
				print("Here")
				merged_temp <- compare_drug_profiles(ref_dis,temp_dis,threshold = 0.01,c(ref_dis_name,temp_dis_name))
				merged_temp <- filter_comparisons(merged_temp)
				merged <- rbind(merged,merged_temp)
			}	
		}
	}
	return(merged)
}

Drug_Level_Comparisons <- function(list_of_dirs,reference,threshold,out_file,names){
	temp_ref <- get(load(file = list_of_dirs[reference]))
	print(length(list_of_dirs)/2 + 1)
	pdf(file=out_file,heigh = 30)
	par(mfrow=c(length(list_of_dirs)/2 + 1,2))
	for(i in 1:length(list_of_dirs)){
		if(i != reference){
			temp_1 <- get(load(file = list_of_dirs[i]))
			temp_r1 <- cbind(temp_ref$original_distancies,temp_1$original_distancies,temp_ref$fdr,temp_1$fdr)
			temp_r1 <- temp_r1[(temp_r1[,3] < threshold) & (temp_r1[,4] < threshold),]
			plot(temp_r1[,1],temp_r1[,2],main = paste(names[reference],"Vs",names[i]),xlab = names[reference], ylab = names[i])
		}	
	}
	dev.off()	
}
