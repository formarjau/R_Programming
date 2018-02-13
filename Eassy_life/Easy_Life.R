#Function that transform NA of a vector to UNDET string value.

NA_To_Undet <- function(vec){
	vec[is.na(vec) | vec == "NA" | vec == "N/A"] <- "UNDET"
	return(vec)
}

Sp_To_UnSc <- function(vec){
	vec <- gsub(" ","_",vec)
	return(vec)
}


