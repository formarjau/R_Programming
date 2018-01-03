#Funtion that compares GSEA enrichment data

compare_GSEA_Enrichment <- function(folder_A,folder_B,threshold){
	folder_A_content <- dir(folder_A,full.names = TRUE)
	up_file_A <- folder_A_content[grepl("report_for_na_pos",folder_A_content) & grepl("xls",folder_A_content)]
	up_file_A <- read.csv(up_file_A,sep="\t")
	up_file_A <- up_file_A[as.numeric(up_file_A$FDR.q.val) < threshold,]
	up_file_A$dir <- rep("UP",nrow(up_file_A))
	down_file_A <- folder_A_content[grepl("report_for_na_neg",folder_A_content) & grepl("xls",folder_A_content)]
	down_file_A <- read.csv(down_file_A,sep="\t")
	down_file_A <- down_file_A[as.numeric(down_file_A$FDR.q.val) < threshold,]
	down_file_A$dir <- rep("DOWN",nrow(down_file_A))
	data_A <- rbind(up_file_A,down_file_A)
	colnames(data_A) <- paste("data_A_",colnames(data_A),sep="")
	folder_B_content <- dir(folder_B,full.names = TRUE)
	up_file_B <- folder_B_content[grepl("report_for_na_pos",folder_B_content) & grepl("xls",folder_B_content)]
	up_file_B <- read.csv(up_file_B,sep="\t")
	up_file_B <- up_file_B[as.numeric(up_file_B$FDR.q.val) < threshold,]
	up_file_B$dir <- rep("UP",nrow(up_file_B))
	down_file_B <- folder_B_content[grepl("report_for_na_neg",folder_B_content) & grepl("xls",folder_B_content)]
	down_file_B <- read.csv(down_file_B,sep="\t")
	down_file_B <- down_file_B[as.numeric(down_file_B$FDR.q.val) < threshold,]
	print(head(down_file_B))
	down_file_B$dir <- rep("DOWN",nrow(down_file_B))
	data_B <- rbind(up_file_B,down_file_B)
	colnames(data_B) <- paste("data_B_",colnames(data_B),sep="")
	data_All <- merge(data_A,data_B,by.x="data_A_NAME",by.y="data_B_NAME")
	data_All_filt <- data_All[,c(1,13,25)]
	out_list <- list(data_All,data_All_filt)
	return(out_list)
}

