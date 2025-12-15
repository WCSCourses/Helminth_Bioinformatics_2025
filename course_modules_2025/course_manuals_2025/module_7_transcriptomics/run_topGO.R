run_topGO_R = function(ref, genelist, thres = 0.05) {
  
  ref=read.table(file=ref, stringsAsFactors=FALSE)
  names(ref) = c('id', 'go')
  ref.vec = strsplit(ref$go, split=',', fixed=T)
  names(ref.vec) <- ref$id
  all.ids <- ref$id
  
  vec<-genelist
  
  scores <- rep(0, nrow(ref)) # list of scores
  names(scores) <- ref$id
  scores[ref$id %in% vec] <- 1
  
  geneSelectionFun <- function(score){
    return(score >= 1)
  }
  
  GOdataBP <- new("topGOdata",  ontology = 'BP', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = 5, description = '')
  GOdataMF <- new("topGOdata",  ontology = 'MF', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = 5, description = '')
  GOdataCC <- new("topGOdata",  ontology = 'CC', allGenes = scores, annot = annFUN.gene2GO, gene2GO = ref.vec, geneSelectionFun = geneSelectionFun, nodeSize = 5, description = '')
  
  resultTopgoBP <- runTest(GOdataBP,algorithm="weight01",statistic="Fisher")
  resultTopgoMF <- runTest(GOdataMF,algorithm="weight01",statistic="Fisher")
  resultTopgoCC <- runTest(GOdataCC,algorithm="weight01",statistic="Fisher")
  
  resBP<-GenTable( GOdataBP, topGO = resultTopgoBP, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
  resMF<-GenTable( GOdataMF, topGO = resultTopgoMF, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
  resCC<-GenTable( GOdataCC, topGO = resultTopgoCC, orderBy = "topGO", ranksOf = "fisher", topNodes = 50)
  
  #     writeLines(capture.output(genesInTerm(GOdataBP)), "topgo_out_terms_bp.txt")
  # 		writeLines(capture.output(genesInTerm(GOdataMF)), "topgo_out_terms_mf.txt")
  # 		writeLines(capture.output(genesInTerm(GOdataCC)), "topgo_out_terms_cc.txt")
  
  # 		write.table(resBP, "topgo_out_bp.txt", sep = "	", quote = FALSE)
  # 		write.table(resMF, "topgo_out_mf.txt", sep = "	", quote = FALSE)
  # 		write.table(resCC, "topgo_out_cc.txt", sep = "	", quote = FALSE)
  
  resBP  #MF CC
  genesInTerm(GOdataBP)  #MF CC
  vec
  # thres = 0.05
  
  ##### BP #####
  tempBP <- as.data.frame(matrix(nrow = dim(resBP[resBP$topGO <= thres,])[1], ncol = 8))
  tempBP[,2:7] <- resBP[resBP$topGO <= thres,] 
  if(dim(resBP[resBP$topGO <= thres,])[1] > 0) {
    tempBP$V1 <- c("BP")
    
    for(go in tempBP$V2) {
      gene <- genesInTerm(GOdataBP,whichGO = go)[[1]][which(genesInTerm(GOdataBP,whichGO = go)[[1]] %in% vec)] ## To add at the end of each row in resBP
      gene <- paste(gene, sep = "", collapse = ",")
      tempBP[which(tempBP$V2 == go),8] <- gene
    }
  }
  
  ##### MF #####
  tempMF <- as.data.frame(matrix(nrow = dim(resMF[resMF$topGO <= thres,])[1], ncol = 8))
  tempMF[,2:7] <- resMF[resMF$topGO <= thres,] 
  if(dim(resMF[resMF$topGO <= thres,])[1] > 0) {
    tempMF$V1 <- c("MF")
    
    for(go in tempMF$V2) {
      gene <- genesInTerm(GOdataMF,whichGO = go)[[1]][which(genesInTerm(GOdataMF,whichGO = go)[[1]] %in% vec)] ## To add at the end of each row in resMF
      gene <- paste(gene, sep = "", collapse = ",")
      tempMF[which(tempMF$V2 == go),8] <- gene
    }
  }
  
  ##### CC #####
  tempCC <- as.data.frame(matrix(nrow = dim(resCC[resCC$topGO <= thres,])[1], ncol = 8))
  tempCC[,2:7] <- resCC[resCC$topGO <= thres,] 
  if(dim(resCC[resCC$topGO <= thres,])[1] > 0) {
    tempCC$V1 <- c("CC")
    
    for(go in tempCC$V2) {
      gene <- genesInTerm(GOdataCC,whichGO = go)[[1]][which(genesInTerm(GOdataCC,whichGO = go)[[1]] %in% vec)] ## To add at the end of each row in resCC
      gene <- paste(gene, sep = "", collapse = ",")
      tempCC[which(tempCC$V2 == go),8] <- gene
    }		
  }
  
  ##### Final #####
  tempBP
  tempMF
  tempCC
  resFinal <- rbind(tempBP, tempMF, tempCC)
  resFinal[,1:7]
  resFinal[,8]
  colnames(resFinal) <- c("type", colnames(resBP), "Smp")
  
  return(resFinal)
  
  #		write.table(resFinal, "topgo_out_final.txt", sep = "	", quote = FALSE)
  
}

## Combination of scripts from ar11, ap6, aw17
