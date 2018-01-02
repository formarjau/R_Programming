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
	list_out <- list()
	intersected <- Common_Set_Genes(vec_of_studies)
	names_ds <- c()
	for(i in 1:length(vec_of_studies)){
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

#Plotting MetaQC meta-analysis quality control results. It could stop working in further versions of the package.

Plot_MetaQC <- function(., .scale.coord.var=4, isCAQC=FALSE,xl = 8, yl = 8 ) {
  if(is.null(.$.Scores))
    .$RunQC()
  .dat <- apply(.$.Scores, 2, function(s) {
    scale(s)
  })
  
  .dummy <- apply(.$.Scores, 2, function(s) {
    (-log10(.05/length(s)) - mean(s)) / sd(s)
  })
  
  if(isCAQC) {
    .dat <- cbind(.dat[,-match(c('CQCg','AQCg'),colnames(.dat))], CAQCg=rowMeans(.dat[,match(c('CQCg','AQCg'),colnames(.dat))]))
    .dummy <- c(.dummy[-match(c('CQCg','AQCg'),names(.dummy))], CAQCg=mean(.dummy[match(c('CQCg','AQCg'),names(.dummy))]))
    .dat <- cbind(.dat[,-match(c('CQCp','AQCp'),colnames(.dat))], CAQCp=rowMeans(.dat[,match(c('CQCp','AQCp'),colnames(.dat))]))
    .dummy <- c(.dummy[-match(c('CQCp','AQCp'),names(.dummy))], CAQCp=mean(.dummy[match(c('CQCp','AQCp'),names(.dummy))]))
  }
  
  .res <- prcomp(.dat, center=FALSE)
  
  .coord.dummy <- .dummy %*% .res$rotation[,1:2]
  .coord <- .res$x[,1:2]
  .coord <- sweep(.coord, 2, .coord.dummy)
  .coord.var <- sweep(.res$rotation, 2, .res$sdev, "*")[,1:2] #idea from FactoMineR
  
  #force plots be represented to right-upper side
  .sign <- sign(colSums(sign(.coord.var)))
  .sign <- ifelse(.sign>=0,1,-1)
  .coord <- sweep(.coord, 2, .sign, '*') 
  .coord.var <- sweep(.coord.var, 2, .sign, '*')
  
  .pctEig <- (.res$sdev^2/sum(.res$sdev^2)*100)[1:2]
  
  plot(x=.coord[,1],y=.coord[,2],type="n",xlab=bquote(bold(.(sprintf("1st Principal Component (%2.2f%%)",.pctEig[1])))),
       ylab=bquote(bold(.(sprintf("2nd Principal Component (%2.2f%%)",.pctEig[2])))),
       xlim=c(-xl,xl),
       ylim=c(-yl,yl),
       axes=FALSE)
  axis(side=1,lwd=4,tck=-0.02)
  axis(side=2,lwd=4,tck=-0.02)
  box(bty="L",lwd=4)
  
  abline(v=0, lty=2, lwd=2)
  abline(h=0, lty=2, lwd=2)
  
  for (v in 1:nrow(.coord.var)) {
    arrows(0, 0, .coord.var[v, 1]*.scale.coord.var, .coord.var[v, 2]*.scale.coord.var, 
           lwd=3, length = 0.1, angle = 15, code = 2, col=gray(.4)) 
    text(.coord.var[v, 1]*.scale.coord.var, y = .coord.var[v, 2]*.scale.coord.var, 
         labels = bquote(bold(.(rownames(.coord.var)[v]))), pos = 3, cex=1)
  }
  
  points(x=.coord[,1],y=.coord[,2],pch=1,col="black",cex=2.5,lwd=2)
  text(x=.coord[,1],y=.coord[,2],cex=1) 
}