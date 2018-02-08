#Function that transform NA of a vector to UNDET string value.

NA_To_Undet <- function(vec){
	vec[is.na(vec) | vec == "NA"] <- "UNDET"
	return(vec)
}