#' @keywords internal
#' @title WebCMap: An ultra-fast connectivity analysis tool using web technologies to accelerate computations based on the extended Connectivity Map (CMap).
#'
#' @description An ultra-fast connectivity analysis tool using web technologies to accelerate computations based on the extended Connectivity Map (CMap).
#'
#' @name WebCMap
#'
#' @author Hongen Kang  \email{kanghongen@gmail.com}
#'
#' @import parallel
#' @import httr
#' @import org.Hs.eg.db
#' @importFrom jsonlite fromJSON
#' @importFrom assertthat assert_that
#' @importFrom lsa cosine
#' @importFrom utils head tail
#' @importFrom AnnotationDbi mapIds
#' @importFrom stats cor
#' @import fmsb
#' @import fgsea
#' @import ggvenn
#' @import ggpubr

## usethis namespace: start
## usethis namespace: end
NULL
options(scipen = 999)


# data(query_signature)
# system.time({
#   res = run_negative_CMap(query_signature = query_signature,K=50,cores=10)
# })
