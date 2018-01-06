#Computes one-tailed fisher test from differential expressoin meta-analysis results.

IntersectionPvals <- function(ResultA,ResultB,FDRA,FDRB){
  ResultA <- ResultA[intersect(rownames(ResultA),rownames(ResultB)),]
  ResultB <- ResultB[intersect(rownames(ResultA),rownames(ResultB)),]
  Common <- length(intersect(rownames(ResultA),rownames(ResultB)))
  ResultAfilt <- ResultA[as.numeric(as.character(ResultA$FDR)) < FDRA,]
  ResultAUp <- ResultAfilt[as.numeric(as.character(ResultAfilt$zval)) > 0 ,]
  ResultADown <- ResultAfilt[as.numeric(as.character(ResultAfilt$zval)) < 0 ,]
  ResultBfilt <- ResultB[as.numeric(as.character(ResultB$FDR))< FDRB,]
  ResultBUp <- ResultBfilt[as.numeric(as.character(ResultBfilt$zval)) > 0 ,]
  print(ResultBUp)
  ResultBDown <- ResultBfilt[as.numeric(as.character(ResultBfilt$zval)) < 0 ,]
  UPUP <- length(intersect(as.character(ResultAUp$symb),as.character(ResultBUp$symb)))
  counts = (matrix(data=c(UPUP,nrow(ResultAUp) - UPUP,nrow(ResultBUp)- UPUP,Common + UPUP - nrow(ResultAUp) - nrow(ResultBUp)),nrow=2))
  testUPUP <- fisher.test(counts,alternative="greater")$p.value
  DOWNDOWN <- length(intersect(as.character(ResultADown$symb),as.character(ResultBDown$symb)))
  counts = (matrix(data=c(DOWNDOWN,nrow(ResultADown) - DOWNDOWN,nrow(ResultBDown)-DOWNDOWN,Common + DOWNDOWN - nrow(ResultADown) - nrow(ResultBDown)),nrow=2))
  testDOWNDOWN <- fisher.test(counts,alternative="greater")$p.value
  UPDOWN <- length(intersect(as.character(ResultAUp$symb),as.character(ResultBDown$symb)))
  counts = (matrix(data=c(UPDOWN,nrow(ResultAUp) - UPDOWN,nrow(ResultBDown) - UPDOWN,Common + UPDOWN - nrow(ResultAUp) - nrow(ResultBDown)),nrow=2))
  testUPDOWN <- fisher.test(counts,alternative="greater")$p.value
  DOWNUP <- length(intersect(as.character(ResultADown$symb),as.character(ResultBUp$symb)))
  counts = (matrix(data=c(DOWNUP,nrow(ResultADown) - DOWNUP,nrow(ResultBUp) - DOWNUP,Common  + DOWNUP - nrow(ResultADown) - nrow(ResultBUp)),nrow=2))
  testDOWNUP <- fisher.test(counts,alternative="greater")$p.value
  listOut <- list(testUPUP,testDOWNDOWN,testUPDOWN,testDOWNUP)
  return(listOut)
}