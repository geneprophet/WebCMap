#' Verify if the input data.frame of the query signature contains both the gene_id and zscore columns.
#' @param query_signature \code{data.frame} object
#'
#' @return A \code{logical} object
#' @export
#'
#' @examples
#' data(query_signature)
#' check_query_signature(query_signature)
#'
check_query_signature = function(query_signature){
  return(all(c(is.element("gene_id",colnames(query_signature)),is.element("zscore",colnames(query_signature)))))
}
#' screening candidate signatures
#'
#' @param query_signature input
#' @param K the top K and bottom K genes
#' @param cores the cores used for computation
#'
#' @return a \code{vector}
#' @export
#'
screen_candidate_signature = function(query_signature,K=50,cores=5){
  if(check_query_signature(query_signature)){
    query_signature = query_signature[order(query_signature$zscore,decreasing=T,na.last=NA),]
    query_up = head(query_signature$gene_id,K)
    query_down = tail(query_signature$gene_id,K)
    clus = parallel::makeCluster(5)
    parallel::clusterExport(clus,varlist=c('K','query_signature'),envir = environment())
    res <- parallel::parLapply(clus, query_up, checkExistenceInCmapDown)
    query_up_res <- do.call(c,res)
    res <- parallel::parLapply(clus, query_down, checkExistenceInCmapUp)
    query_down_res <- do.call(c,res)
    parallel::stopCluster(clus)
    all_candidate_signature_index <- intersect(query_up_res,query_down_res)
    return(all_candidate_signature_index)
  }else{
    print("The query signature data frame must have gene_id and zscore columns")
  }
}
