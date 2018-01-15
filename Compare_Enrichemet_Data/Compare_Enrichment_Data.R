#Funtion that compares GSEA enrichment data

compare_GSEA_Enrichment <- function(folder_A,folder_B,threshold){
	folder_A_content <- dir(folder_A,full.names = TRUE)
	up_file_A <- folder_A_content[grepl("report_for_na_pos",folder_A_content) & grepl("xls",folder_A_content)]
	up_file_A <- read.csv(up_file_A,sep="\t")
	up_file_A <- up_file_A[,1:11]
	up_file_A <- up_file_A[as.numeric(up_file_A$FDR.q.val) < threshold,]
	up_file_A$dir <- rep("UP",nrow(up_file_A))
	down_file_A <- folder_A_content[grepl("report_for_na_neg",folder_A_content) & grepl("xls",folder_A_content)]
	down_file_A <- read.csv(down_file_A,sep="\t")
	down_file_A <- down_file_A[,1:11]
	down_file_A <- down_file_A[as.numeric(down_file_A$FDR.q.val) < threshold,]
	down_file_A$dir <- rep("DOWN",nrow(down_file_A))
	data_A <- rbind(up_file_A,down_file_A)
	colnames(data_A) <- paste("data_A_",colnames(data_A),sep="")
	folder_B_content <- dir(folder_B,full.names = TRUE)
	up_file_B <- folder_B_content[grepl("report_for_na_pos",folder_B_content) & grepl("xls",folder_B_content)]
	up_file_B <- read.csv(up_file_B,sep="\t")
	up_file_B <- up_file_B[,1:11]
	up_file_B <- up_file_B[as.numeric(up_file_B$FDR.q.val) < threshold,]
	up_file_B$dir <- rep("UP",nrow(up_file_B))
	down_file_B <- folder_B_content[grepl("report_for_na_neg",folder_B_content) & grepl("xls",folder_B_content)]
	down_file_B <- read.csv(down_file_B,sep="\t")
	down_file_B <- down_file_B[,1:11]
	down_file_B <- down_file_B[as.numeric(down_file_B$FDR.q.val) < threshold,]
	down_file_B$dir <- rep("DOWN",nrow(down_file_B))
	data_B <- rbind(up_file_B,down_file_B)
	colnames(data_B) <- paste("data_B_",colnames(data_B),sep="")
	data_All <- merge(data_A,data_B,by.x="data_A_NAME",by.y="data_B_NAME")
	print(head(data_All))
	data_All_filt <- data_All[,c(1,12,23)]
	out_list <- list(data_All,data_All_filt)
	return(out_list)
}

#Funtions to generate the plots showing the number of pathways deregulated in the same or different directions.

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

#Contructing disease disease tables

Construct_Dis_Dis_Tables <- function(list_dirs,FDR_threshold){
		D_A_Pos <- read_GSEA_results(get_GSEA_result_in_dir(list_dirs[[1]],"POS"),FDR_threshold)
		print(head(D_A_Pos))
		rownames(D_A_Pos) <- D_A_Pos[,1]
		D_A_Neg <- read_GSEA_results(get_GSEA_result_in_dir(list_dirs[[1]],"NEG"),FDR_threshold)
		print(head(D_A_Neg))
		rownames(D_A_Neg) <- D_A_Neg[,1]
		D_B_Pos <- read_GSEA_results(get_GSEA_result_in_dir(list_dirs[[2]],"POS"),FDR_threshold)
		print(head(D_B_Pos))
		rownames(D_B_Pos) <- D_B_Pos[,1]
		D_B_Neg <- read_GSEA_results(get_GSEA_result_in_dir(list_dirs[[2]],"NEG"),FDR_threshold)
		print(head(D_B_Neg))
		rownames(D_B_Neg) <- D_B_Neg[,1]
		all_paths <- unique(c(rownames(D_A_Pos),rownames(D_A_Neg),rownames(D_B_Pos),rownames(D_B_Neg)))
		up_or_down_A <- rep(NA,length(all_paths))
		up_or_down_B <- rep(NA,length(all_paths))
		up_or_down_A[all_paths %in% rownames(D_A_Pos)] <- "UP"
		up_or_down_A[all_paths %in% rownames(D_A_Neg)] <- "DOWN"
		up_or_down_B[all_paths %in% rownames(D_B_Pos)] <- "UP"
		up_or_down_B[all_paths %in% rownames(D_B_Neg)] <- "DOWN"
		out <- data.frame(all_paths,up_or_down_A,up_or_down_B)
    dis_A <- strsplit(list_dirs[[1]],"\\/")[[1]][length(strsplit(list_dirs[[1]],"\\/")[[1]])]
    dis_B <- strsplit(list_dirs[[2]],"\\/")[[1]][length(strsplit(list_dirs[[2]],"\\/")[[1]])]
    colnames(out) <- c("Paths",dis_A,dis_B)
		return(out)
}

#Countion directions.

Count_directions_And_Perc <- function(pair){
  pair <- na.omit(pair)
  up_up_number <- sum(pair[,2] == "UP" & pair[,3] == "UP")
  down_down_number <- sum(pair[,2] == "DOWN" & pair[,3] == "DOWN")
  up_down_number <- sum(pair[,2] == "UP" & pair[,3] == "DOWN")
  down_up_number <- sum(pair[,2] == "DOWN" & pair[,3] == "UP")
  same_direction_number <- sum(pair[,2] == pair[,3])
  same_direction_proportion <- (sum(pair[,2] == pair[,3])/nrow(pair))*100
  different_direction_number <- sum(pair[,2] != pair[,3])
  different_direction_proportion <- (sum(pair[,2] != pair[,3])/nrow(pair))*100
  df_sal <- data.frame(t(c(colnames(pair)[2],colnames(pair)[3],up_up_number,down_down_number,same_direction_number,same_direction_proportion,up_down_number,down_up_number,different_direction_number,different_direction_proportion)))
  colnames(df_sal) <- c("Dis_A","Dis_B","up_up_number","down_down_number","same_direction_number","same_direction_proportion","up_down_number","down_up_number","different_direction_number","different_direction_proportion")
  print(df_sal)
}

#Countion directions for all

count_for_all <- function(vec_dirs,reference,FDR = 0.05){
  count <- 0
  for(i in 1:length(vec_dirs)){
    if(i != reference & count == 0){
      list_comp <- list(vec_dirs[reference],vec_dirs[i])
      temp <- Construct_Dis_Dis_Tables(list_comp,FDR)
      df <- Count_directions_And_Perc(temp)
      count <- count + 1
    }else if(i != reference){
      list_comp <- list(vec_dirs[reference],vec_dirs[i])
      temp <- Construct_Dis_Dis_Tables(list_comp,FDR)
      df <- rbind(df,Count_directions_And_Perc(temp))
    }
  }
  return(df)
}

#Transform data to generate the barplot.

transform_for_bplot <- function(df){
  names <- paste(df[,1],df[,2],sep="_")
  same <- df$same_direction_proportion
  same_val <- rep("same",length(same))
  diff <- df$different_direction_proportion
  diff_val <- rep("diff",length(diff))
  out <- data.frame(names,c(same,diff),c(same_val,diff_val))
  colnames(out) <- c("names","Per","class")
  return(out)
}

#Create a pallete of the desaired colors.


create_pallete <- function(vector_colors){
  hex_sal <- c()
  for(i in 1:length(vector_colors)){
    temp_rgb <- col2rgb(vector_colors[i], alpha = FALSE)
    red <- as.numeric(temp_rgb[1,1])
    green <- as.numeric(temp_rgb[2,1])
    blue <- as.numeric(temp_rgb[3,1])
    temp_hex = rgb(red = red,green = green,blue = blue, maxColorValue=255)
    hex_sal <- c(hex_sal,temp_hex)
  }
  return(hex_sal)
}