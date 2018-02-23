##############################################
#Functions to work with directory structures.#
##############################################

#Create review structure

create_dir_for_review <- function(directory,disease){
	if(!dir.exists(file.path(directory,"Busquedas"))){
		dir.create(file.path(directory,"Busquedas"))
	}
	if(!(file.exists(file.path(directory,"Busquedas",paste("Screen_cetations_v1_",disease,".xlsx",sep=""))))){
		file.create(file.path(directory,"Busquedas",paste("Screen_cetations_v1_",disease,".xlsx",sep="")))
	}
	if(!dir.exists(file.path(directory,"Deatasets_Incluidos_Data"))){
		dir.create(file.path(directory,"Deatasets_Incluidos_Data"))
	}
	if(!(file.exists(file.path(directory,"Deatasets_Incluidos_Data",paste(disease,"_Dataset.xlsx",sep=""))))){
		file.create(file.path(directory,"Deatasets_Incluidos_Data",paste(disease,"_Dataset.xlsx",sep="")))
	}
	if(!dir.exists(file.path(directory,"Dudas"))){
		dir.create(file.path(directory,"Dudas"))
	}
	if(!(file.exists(file.path(directory,"Dudas",paste(disease,"_Dataset.xlsx",sep=""))))){
		file.create(file.path(directory,"Dudas",paste(disease,"_Dataset.xlsx",sep="")))
	}
	if(!dir.exists(file.path(directory,"Observational_Studies_PDFs"))){
		dir.create(file.path(directory,"Observational_Studies_PDFs"))
	}
	if(!dir.exists(file.path(directory,"Observational_Studies_Searches"))){
		dir.create(file.path(directory,"Observational_Studies_Searches"))
	}
	if(!dir.exists(file.path(directory,"Observational_Studies_Searches","Search_1"))){
		dir.create(file.path(directory,"Observational_Studies_Searches","Search_1"))
		file.create(file.path(directory,"Observational_Studies_Searches","Search_1","Obs_Est_Searhc.txt"))
		file.create(file.path(directory,"Observational_Studies_Searches","Search_1","Search_1.pptx"))
		file.create(file.path(directory,"Observational_Studies_Searches","Search_1","Search_1.txt"))
		file.create(file.path(directory,"Observational_Studies_Searches","Search_1","Search_1.xlsx"))
		file.create(file.path(directory,"Observational_Studies_Searches","Search_1","Selected.xlsx"))
	}
	if(!dir.exists(file.path(directory,"Scripts"))){
		dir.create(file.path(directory,"Scripts"))
	}
	if(!dir.exists(file.path(directory,"Supplementary"))){
		dir.create(file.path(directory,"Supplementary"))
	}
	if(!(file.exists(file.path(directory,"Supplementary",paste(disease,"_Supplementary_Anorexia.xlsx",sep=""))))){
		file.create(file.path(directory,"Supplementary",paste(disease,"_Supplementary_",disease,".docx",sep="")))
	}
	if(!dir.exists(file.path(directory,"Systematic_Review_Searches"))){
		dir.create(file.path(directory,"Systematic_Review_Searches"))
	}
	if(!dir.exists(file.path(directory,"Systematic_Review_Searches","Search_1"))){
		dir.create(file.path(directory,"Systematic_Review_Searches","Search_1"))
		file.create(file.path(directory,"Systematic_Review_Searches","Search_1","Search_1.docx"))
		file.create(file.path(directory,"Systematic_Review_Searches","Search_1","Search_1.pptx"))
		file.create(file.path(directory,"Systematic_Review_Searches","Search_1","Selected.xlsx"))
		file.create(file.path(directory,"Systematic_Review_Searches","Search_1","Sys_Rev_Search.txt"))
	}
}

#Add search folder

add_search_folder <- function(directory,type = c("observational","systematic_review")){
	if(type == "observational"){
		temp <- dir(file.path(directory,"Observational_Studies_Searches"))
		print(temp)
		next_search <- max(as.numeric(unlist(lapply(strsplit(temp,"_"),FUN = function(x) {x[2]})))) + 1
		dir.create(file.path(directory,"Observational_Studies_Searches",paste("Search_",next_search,sep="")))
		file.create(file.path(directory,"Observational_Studies_Searches",paste("Search_",next_search,sep=""),"Obs_Est_Searhc.txt"))
		file.create(file.path(directory,"Observational_Studies_Searches",paste("Search_",next_search,sep=""),paste("Search_",next_search,".pptx",sep="")))
		file.create(file.path(directory,"Observational_Studies_Searches",paste("Search_",next_search,sep=""),paste("Search_",next_search,".txt",sep="")))
		file.create(file.path(directory,"Observational_Studies_Searches",paste("Search_",next_search,sep=""),paste("Search_",next_search,".xlsx",sep="")))
		file.create(file.path(directory,"Observational_Studies_Searches",paste("Search_",next_search,sep=""),"Selected.xlsx"))

	}
	else if(type == "systematic_review"){
		temp <- dir(file.path(directory,"Systematic_Review_Searches"))
		print(temp)
		next_search <- max(as.numeric(unlist(lapply(strsplit(temp,"_"),FUN = function(x) {x[2]})))) + 1
		dir.create(file.path(directory,"Systematic_Review_Searches",paste("Search_",next_search,sep="")))
		file.create(file.path(directory,"Systematic_Review_Searches",paste("Search_",next_search,sep=""),paste("Search_",next_search,".docx",sep="")))
		file.create(file.path(directory,"Systematic_Review_Searches",paste("Search_",next_search,sep=""),paste("Search_",next_search,".pptx")))
		file.create(file.path(directory,"Systematic_Review_Searches",paste("Search_",next_search,sep=""),"Selected.xlsx"))
		file.create(file.path(directory,"Systematic_Review_Searches",paste("Search_",next_search,sep=""),"Sys_Rev_Search.txt"))
	}
}