#utils
# Gene symbol search function
#library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
#gene = select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c("SYMBOL","GENETYPE","ENSEMBLTRANS"))
#pcg = unique(gene$SYMBOL[which(gene$GENETYPE=="protein-coding")])

#' Convert gene Symbol to ENTREZ ID
#'
#' @param a vector of gene symbols
#'
#' @return a vector of ENTREZ IDs
#' @export
#'
#' @examples
#' entrezid <- SYMBOLtoENTREZ("APOE")
SYMBOLtoENTREZ <- function(x){
  vec=AnnotationDbi::mapIds(org.Hs.eg.db,keys=x,column="ENTREZID",keytype="SYMBOL",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

#' Convert ENSEMBL ID to ENTREZ ID
#'
#' @param a vector of ENSEMBL IDs
#'
#' @return a vector of ENTREZ IDs
#' @export
#'
#' @examples
#' entrezid <- ENSEMBLIDtoENTREZ("ENSG00000130203")
ENSEMBLIDtoENTREZ <- function(x){
  vec=AnnotationDbi::mapIds(org.Hs.eg.db,keys=x,column="ENTREZID",keytype="ENSEMBL",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}


