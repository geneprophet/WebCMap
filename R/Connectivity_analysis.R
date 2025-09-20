#' To identify candidate signatures that are negatively correlated with the query signature by computing six kinds of connectivity scores
#'
#' @param signature_index
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#' data(query_signature)
#' K=50
#' res = negative_query_CMap(506400)
#'

negative_query_CMap = function(signature_index) {
  #order and remove NA zscore, two columns: zscore and gene_id
  query_signature = query_signature[order(query_signature$zscore,
                                          decreasing = TRUE,
                                          na.last = NA), ]
  query_up_gene = head(query_signature$gene_id, K)
  query_down_gene = tail(query_signature$gene_id, K)
  query_signature_gene = c(query_up_gene, query_down_gene)
  tryCatch({
    suppressWarnings({
      options(scipen = 999)
      cmap_signature = retrieveCmapSignature(as.numeric(signature_index))
      #order and remove NA modz
      cmap_signature = cmap_signature[order(cmap_signature$modz,
                                            decreasing = TRUE,
                                            na.last = NA), ]
      cmap_up_gene = head(cmap_signature$gene_id, K)
      cmap_down_gene = tail(cmap_signature$gene_id, K)
      cmap_signature_gene = c(cmap_up_gene, cmap_down_gene)
      #1.Spearman,Cosine and Pearson correlation in all overlaped genes and extreme overlaped genes
      # whole metrics
      whole_overlaped_gene = intersect(query_signature$gene_id, cmap_signature$gene_id)
      extreme_overlaped_gene = intersect(query_signature_gene, cmap_signature_gene)
      # duplicated row names
      merge_cmap_query = merge(cmap_signature,
                               query_signature,
                               by = "gene_id",
                               all.X = TRUE)
      Spearman_whole = cor(merge_cmap_query$zscore,
                           merge_cmap_query$modz,
                           method = "spearman")
      Pearson_whole = cor(merge_cmap_query$zscore,
                          merge_cmap_query$modz,
                          method = "pearson")
      Consine_whole = as.vector(lsa::cosine(merge_cmap_query$zscore, merge_cmap_query$modz))
      # #extreme metrics
      # merge_extreme = merge_cmap_query[match(extreme_overlaped_gene,merge_cmap_query$gene_id),]
      # Spearman_extreme = cor(merge_extreme$zscore,merge_extreme$modz,method="spearman")
      # Pearson_extreme = cor(merge_extreme$zscore,merge_extreme$modz,method="pearson")
      # Consine_extreme = as.vector(lsa::cosine(merge_extreme$zscore,merge_extreme$modz))
      if (Spearman_whole < 0 & Pearson_whole < 0 & Consine_whole < 0) {
        #2.Xsum: disease signature: membership of top and bottom K genes; cmap: modz values of top and bottom K genes
        XUpInDisease = intersect(query_up_gene, cmap_signature_gene)
        XDownInDisease = intersect(query_down_gene, cmap_signature_gene)
        XSum = sum(cmap_signature$modz[which(is.element(cmap_signature$gene_id, XUpInDisease))], na.rm = TRUE) - sum(cmap_signature$modz[which(is.element(cmap_signature$gene_id, XDownInDisease))], na.rm = TRUE)
        ## The permutation results is pre-defined and can be retrived from webserver
        permuteNum = 10000
        permuteResult = retrievePermutationResult(signature_index = as.numeric(signature_index),
                                                  method = "XSum")
        permuteResult$score[is.na(permuteResult$score)] <- 0
        ## Compute the p-value based on permutation results
        XSum_pvalue <-
          sum(abs(permuteResult$score) >= abs(XSum)) / permuteNum
        if (XSum < 0 & XSum_pvalue < 0.05) {
          #3.Connection strength score alias zhangscore
          #3.1 Convert the gene expression matrix to ranked list
          refSort = cmap_signature[order(abs(cmap_signature$modz), decreasing =
                                           TRUE), ]
          refSort$refRank = rank(abs(refSort$modz)) * sign(refSort$modz)
          #3.2 Compute the maximal theoretical score
          queryVector = c(rep(1, length(query_up_gene)), rep(-1, length(query_down_gene)))
          names(queryVector) = c(query_up_gene, query_down_gene)
          maxTheoreticalScore = sum(abs(refSort$refRank)[1:length(queryVector)] * abs(queryVector))
          #3.3 The final score
          Connection_strength_score = sum(queryVector * refSort$refRank[match(names(queryVector), refSort$gene_id)], na.rm =
                                            TRUE) / maxTheoreticalScore
          #3.4 permutation to compute p-value
          permuteNum = 10000
          permuteResult = retrievePermutationResult(signature_index = as.numeric(signature_index),
                                                    method = "CSS")
          permuteResult$score[is.na(permuteResult$score)] <- 0
          ## Compute the two tailed p-value based on bootstrap method
          CSS_pvalue = sum(abs(permuteResult$score) >= abs(Connection_strength_score)) / permuteNum
          if (CSS_pvalue < 0.05 & Connection_strength_score < 0) {
            #4.WTCS
            cmap_ranks = cmap_signature$modz
            names(cmap_ranks) = cmap_signature$gene_id
            query_pathway = list()
            query_pathway[["query_up_gene"]] = query_up_gene
            query_pathway[["query_down_gene"]] = query_down_gene
            fgseaRes = fgsea::fgsea(query_pathway, cmap_ranks, gseaParam = 1)
            ES_up = fgseaRes$ES[which(fgseaRes$pathway == "query_up_gene")]
            ES_up_padj = fgseaRes$padj[which(fgseaRes$pathway == "query_up_gene")]
            ES_down = fgseaRes$ES[which(fgseaRes$pathway == "query_down_gene")]
            ES_down_padj = fgseaRes$padj[which(fgseaRes$pathway == "query_down_gene")]
            if (sum(sign(fgseaRes$ES)) == 0) {
              WTCS = (ES_up - ES_down) / 2
            } else{
              WTCS = 0
            }
            if (WTCS < 0) {
              #5.construct the final result data frame
              final_res = data.frame(
                signature_index = signature_index,
                WTCS = WTCS,
                ES_up = ES_up,
                ES_down = ES_down,
                ES_up_padj = ES_up_padj,
                ES_down_padj = ES_down_padj,
                XSum = XSum,
                XSum_pvalue = XSum_pvalue,
                CSS = Connection_strength_score,
                CSS_pvalue = CSS_pvalue,
                Spearman_correlation = Spearman_whole,
                Pearson_correlation = Pearson_whole,
                Consine_similarity = Consine_whole
              )
              return(final_res)
            } else{
              return(NULL)
            }
          }
          else{
            return(NULL)
          }
        } else{
          return(NULL)
        }
      } else{
        return(NULL)
      }
    })
  }, warning = function(w) {
    return(NULL)
  }, error = function(e) {
    return(
      data.frame(
        signature_index = signature_index,
        WTCS = NA,
        ES_up = NA,
        ES_down = NA,
        ES_up_padj = NA,
        ES_down_padj = NA,
        XSum =
          NA,
        XSum_pvalue = NA,
        CSS = NA,
        CSS_pvalue = NA,
        Spearman_correlation =
          NA,
        Pearson_correlation = NA,
        Consine_similarity = NA
      )
    )
  })
}

#' To identify candidate signatures that are positively correlated with the query signature by computing six kinds of connectivity scores
#'
#' @param signature_index
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#' data(query_signature2)
#' query_signature = query_signature2
#' K=50
#' res = positive_query_CMap(519170)
#'
positive_query_CMap = function(signature_index) {
  #order and remove NA zscore, two columns: zscore and gene_id
  query_signature = query_signature[order(query_signature$zscore,
                                          decreasing = TRUE,
                                          na.last = NA), ]
  query_up_gene = head(query_signature$gene_id, K)
  query_down_gene = tail(query_signature$gene_id, K)
  query_signature_gene = c(query_up_gene, query_down_gene)
  tryCatch({
    suppressWarnings({
      options(scipen = 999)
      cmap_signature = retrieveCmapSignature(as.numeric(signature_index))
      #order and remove NA modz
      cmap_signature = cmap_signature[order(cmap_signature$modz,
                                            decreasing = TRUE,
                                            na.last = NA), ]
      cmap_up_gene = head(cmap_signature$gene_id, K)
      cmap_down_gene = tail(cmap_signature$gene_id, K)
      cmap_signature_gene = c(cmap_up_gene, cmap_down_gene)
      #1.Spearman,Cosine and Pearson correlation in all overlaped genes and extreme overlaped genes
      # whole metrics
      whole_overlaped_gene = intersect(query_signature$gene_id, cmap_signature$gene_id)
      extreme_overlaped_gene = intersect(query_signature_gene, cmap_signature_gene)
      # duplicated row names
      merge_cmap_query = merge(cmap_signature,
                               query_signature,
                               by = "gene_id",
                               all.X = TRUE)
      Spearman_whole = cor(merge_cmap_query$zscore,
                           merge_cmap_query$modz,
                           method = "spearman")
      Pearson_whole = cor(merge_cmap_query$zscore,
                          merge_cmap_query$modz,
                          method = "pearson")
      Consine_whole = as.vector(lsa::cosine(merge_cmap_query$zscore, merge_cmap_query$modz))
      # #extreme metrics
      # merge_extreme = merge_cmap_query[match(extreme_overlaped_gene,merge_cmap_query$gene_id),]
      # Spearman_extreme = cor(merge_extreme$zscore,merge_extreme$modz,method="spearman")
      # Pearson_extreme = cor(merge_extreme$zscore,merge_extreme$modz,method="pearson")
      # Consine_extreme = as.vector(lsa::cosine(merge_extreme$zscore,merge_extreme$modz))
      if (Spearman_whole > 0 & Pearson_whole > 0 & Consine_whole > 0) {
        #2.Xsum: disease signature: membership of top and bottom K genes; cmap: modz values of top and bottom K genes
        XUpInDisease = intersect(query_up_gene, cmap_signature_gene)
        XDownInDisease = intersect(query_down_gene, cmap_signature_gene)
        XSum = sum(cmap_signature$modz[which(is.element(cmap_signature$gene_id, XUpInDisease))], na.rm = TRUE) - sum(cmap_signature$modz[which(is.element(cmap_signature$gene_id, XDownInDisease))], na.rm = TRUE)
        ## The permutation results is pre-defined and can be retrived from webserver
        permuteNum = 10000
        permuteResult = retrievePermutationResult(signature_index = as.numeric(signature_index),
                                                  method = "XSum")
        permuteResult$score[is.na(permuteResult$score)] <- 0
        ## Compute the p-value based on permutation results
        XSum_pvalue <-
          sum(abs(permuteResult$score) >= abs(XSum)) / permuteNum
        if (XSum > 0 & XSum_pvalue < 0.05) {
          #3.Connection strength score alias zhangscore
          #3.1 Convert the gene expression matrix to ranked list
          refSort = cmap_signature[order(abs(cmap_signature$modz), decreasing =
                                           TRUE), ]
          refSort$refRank = rank(abs(refSort$modz)) * sign(refSort$modz)
          #3.2 Compute the maximal theoretical score
          queryVector = c(rep(1, length(query_up_gene)), rep(-1, length(query_down_gene)))
          names(queryVector) = c(query_up_gene, query_down_gene)
          maxTheoreticalScore = sum(abs(refSort$refRank)[1:length(queryVector)] * abs(queryVector))
          #3.3 The final score
          Connection_strength_score = sum(queryVector * refSort$refRank[match(names(queryVector), refSort$gene_id)], na.rm =
                                            TRUE) / maxTheoreticalScore
          #3.4 permutation to compute p-value
          permuteNum = 10000
          permuteResult = retrievePermutationResult(signature_index = as.numeric(signature_index),
                                                    method = "CSS")
          permuteResult$score[is.na(permuteResult$score)] <- 0
          ## Compute the two tailed p-value based on bootstrap method
          CSS_pvalue = sum(abs(permuteResult$score) >= abs(Connection_strength_score)) / permuteNum
          if (CSS_pvalue < 0.05 & Connection_strength_score > 0) {
            #4.WTCS
            cmap_ranks = cmap_signature$modz
            names(cmap_ranks) = cmap_signature$gene_id
            query_pathway = list()
            query_pathway[["query_up_gene"]] = query_up_gene
            query_pathway[["query_down_gene"]] = query_down_gene
            fgseaRes = fgsea::fgsea(query_pathway, cmap_ranks, gseaParam = 1)
            ES_up = fgseaRes$ES[which(fgseaRes$pathway == "query_up_gene")]
            ES_up_padj = fgseaRes$padj[which(fgseaRes$pathway == "query_up_gene")]
            ES_down = fgseaRes$ES[which(fgseaRes$pathway == "query_down_gene")]
            ES_down_padj = fgseaRes$padj[which(fgseaRes$pathway == "query_down_gene")]
            if (sum(sign(fgseaRes$ES)) == 0) {
              WTCS = (ES_up - ES_down) / 2
            } else{
              WTCS = 0
            }
            if (WTCS > 0) {
              #5.construct the final result data frame
              final_res = data.frame(
                signature_index = signature_index,
                WTCS = WTCS,
                ES_up = ES_up,
                ES_down = ES_down,
                ES_up_padj = ES_up_padj,
                ES_down_padj = ES_down_padj,
                XSum = XSum,
                XSum_pvalue = XSum_pvalue,
                CSS = Connection_strength_score,
                CSS_pvalue = CSS_pvalue,
                Spearman_correlation = Spearman_whole,
                Pearson_correlation = Pearson_whole,
                Consine_similarity = Consine_whole
              )
              return(final_res)
            } else{
              return(NULL)
            }
          }
          else{
            return(NULL)
          }
        } else{
          return(NULL)
        }
      } else{
        return(NULL)
      }
    })
  }, warning = function(w) {
    return(NULL)
  }, error = function(e) {
    return(
      data.frame(
        signature_index = signature_index,
        WTCS = NA,
        ES_up = NA,
        ES_down = NA,
        ES_up_padj = NA,
        ES_down_padj = NA,
        XSum =
          NA,
        XSum_pvalue = NA,
        CSS = NA,
        CSS_pvalue = NA,
        Spearman_correlation =
          NA,
        Pearson_correlation = NA,
        Consine_similarity = NA
      )
    )
  })
}



#' @title run negative connectivity analysis between query signature and Extend CMap compound-induced signature
#' @param query_signature a \code{data.frame} must include the gene_id and zscore columns.
#' @param K default 50
#' @param cores default 5
#'
#' @return a \code{data.frame} object
#' @export
#'
#' @examples
#' \dontrun{
#'  data(query_signature1)
#'  res = run_negative_CMap(query_signature1,K=50,cores=10)
#' }

run_negative_CMap = function(query_signature,
                             K = 50,
                             cores = 5) {
  #step1: screen candidate signatures by intersecting the query signature with CMap signatures
  assertthat::assert_that(check_query_signature(query_signature))
  candidate_signatures = screen_negative_candidate_signature(query_signature, K =
                                                               K, cores = cores)
  assertthat::assert_that(length(candidate_signatures) > 0)
  ##step2: run negative connectivity analysis
  clus = parallel::makeCluster(cores)
  parallel::clusterExport(clus,
                          varlist = c('K', 'query_signature'),
                          envir = environment())
  res = parallel::parLapply(clus,
                            as.numeric(candidate_signatures),
                            negative_query_CMap)
  all_res = do.call(rbind, res)
  assertthat::assert_that(nrow(all_res) > 0)
  #meta score, top 5% will be assigned as 1
  if (nrow(all_res) >= 20) {
    Meta_score = as.numeric(rank(all_res$WTCS) <= round(nrow(all_res) * 0.05)) +
      as.numeric(rank(all_res$XSum) <= round(nrow(all_res) * 0.05)) +
      as.numeric(rank(all_res$CSS) <= round(nrow(all_res) * 0.05)) +
      as.numeric(rank(all_res$Spearman_correlation) <= round(nrow(all_res) * 0.05)) +
      as.numeric(rank(all_res$Pearson_correlation) <= round(nrow(all_res) * 0.05)) +
      as.numeric(rank(all_res$Consine_similarity) <= round(nrow(all_res) * 0.05))
  } else{
    Meta_score = as.numeric(rank(all_res$WTCS) <= 1) +
      as.numeric(rank(all_res$XSum) <= 1) +
      as.numeric(rank(all_res$CSS) <= 1) +
      as.numeric(rank(all_res$Spearman_correlation) <= 1) +
      as.numeric(rank(all_res$Pearson_correlation) <= 1) +
      as.numeric(rank(all_res$Consine_similarity) <= 1)
  }
  all_res = cbind(all_res[, 1], Meta_score, all_res[, 2:ncol(all_res)])
  colnames(all_res)[1] = "signature_index"
  ##step3: get signature meta information
  signature_meta = parallel::parLapply(clus, all_res$signature_index, retrieveCmapSignatureMeta)
  all_meta = do.call(rbind, signature_meta)
  all_meta = all_meta[, -c(1, 2)]
  parallel::stopCluster(clus)
  results = cbind(all_res, all_meta)
  results = results[order(results$Meta_score, decreasing = TRUE), ]
  return(results)
}



#' @title run positive connectivity analysis between query signature and Extend CMap compound-induced signature
#' @param query_signature a \code{data.frame} must include the gene_id and zscore columns.
#' @param K default 50
#' @param cores default 5
#'
#' @return a \code{data.frame} object
#' @export
#'
#' @examples
#' \dontrun{
#'  data(query_signature2)
#'  res = run_positive_CMap(query_signature2,K=50,cores=10)
#' }
#'
run_positive_CMap = function(query_signature,
                             K = 50,
                             cores = 5) {
  #step1: screen candidate signatures by intersecting the query signature with CMap signatures
  assertthat::assert_that(check_query_signature(query_signature))
  candidate_signatures = screen_positive_candidate_signature(query_signature, K =
                                                               K, cores = cores)
  assertthat::assert_that(length(candidate_signatures) > 0)
  ##step2: run negative connectivity analysis
  clus = parallel::makeCluster(cores)
  parallel::clusterExport(clus,
                          varlist = c('K', 'query_signature'),
                          envir = environment())
  res = parallel::parLapply(clus,
                            as.numeric(candidate_signatures),
                            positive_query_CMap)
  all_res = do.call(rbind, res)
  assertthat::assert_that(nrow(all_res) > 0)
  #meta score, top 5% will be assigned as 1
  if (nrow(all_res) >= 20) {
    Meta_score = as.numeric(rank(all_res$WTCS,na.last = FALSE) > (nrow(all_res) - round(nrow(all_res) * 0.05))) +
      as.numeric(rank(all_res$XSum,na.last = FALSE) > (nrow(all_res) - round(nrow(all_res) * 0.05))) +
      as.numeric(rank(all_res$CSS,na.last = FALSE) > (nrow(all_res) - round(nrow(all_res) * 0.05))) +
      as.numeric(rank(all_res$Spearman_correlation,na.last = FALSE) > (nrow(all_res) - round(nrow(all_res) *0.05))) +
      as.numeric(rank(all_res$Pearson_correlation,na.last = FALSE) > (nrow(all_res) - round(nrow(all_res) * 0.05))) +
      as.numeric(rank(all_res$Consine_similarity,na.last = FALSE) > (nrow(all_res) - round(nrow(all_res) * 0.05)))
  } else{
    Meta_score = as.numeric(rank(all_res$WTCS,na.last = FALSE) >= nrow(all_res)) +
      as.numeric(rank(all_res$XSum,na.last = FALSE) >= nrow(all_res)) +
      as.numeric(rank(all_res$CSS,na.last = FALSE) >= nrow(all_res)) +
      as.numeric(rank(all_res$Spearman_correlation,na.last = FALSE) >= nrow(all_res)) +
      as.numeric(rank(all_res$Pearson_correlation,na.last = FALSE) >= nrow(all_res)) +
      as.numeric(rank(all_res$Consine_similarity,na.last = FALSE) >= nrow(all_res))
  }

  all_res = cbind(all_res[, 1], Meta_score, all_res[, 2:ncol(all_res)])
  colnames(all_res)[1] = "signature_index"
  ##step3: get signature meta information
  signature_meta = parallel::parLapply(clus, all_res$signature_index, retrieveCmapSignatureMeta)
  all_meta = do.call(rbind, signature_meta)
  all_meta = all_meta[, -c(1, 2)]
  parallel::stopCluster(clus)
  results = cbind(all_res, all_meta)
  results = results[order(results$Meta_score, decreasing = TRUE), ]
  return(results)
}
