#' checkExistenceInCmapUp
#'
#' @param gene_id the ENTREZ id
#'
#' @return candidate signature index
#' @export
#' @examples
#' data(query_signature)
#' K=50
#' checkExistenceInCmapUp(gene_id="348")
#'


checkExistenceInCmapUp = function(gene_id){
  gene_id = as.numeric(gene_id)
  assertthat::assert_that(assertthat::is.number(gene_id),msg = "the gene id is not included in the extended CMap")
  # 最大重试次数
  max_retries <- 10
  # 初始化重试计数器
  retries <- 0
  # 定义一个标志变量，用于控制循环
  success <- FALSE
  # 循环发送请求，直到成功或达到最大重试次数
  while (!success && retries < max_retries) {
    # 发送 GET 请求
    response_up <- httr::GET(paste0("https://ngdc.cncb.ac.cn/cedr/webcmap/api/cmapzscoreupk",K,"?gene_id=",gene_id))
    # 检查状态码
    if (httr::status_code(response_up) == 200) {
      success <- TRUE  # 请求成功，退出循环
    } else {
      retries <- retries + 1  # 增加重试次数
      #等待再重试
      Sys.sleep(retries)
    }
  }
  # 如果最终请求成功，处理响应数据
  if (success) {
    json_data <- rawToChar(response_up[["content"]])
    parsed_data <- jsonlite::fromJSON(json_data)
    if(!is.null(parsed_data$data)){
      return(parsed_data$data$sig_index)
    }else{
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' checkExistenceInCmapDown
#'
#' @param gene_id the ENTREZ id
#'
#' @return candidate signature index
#' @export
#' @examples
#'
#' data(query_signature)
#' K=50
#' checkExistenceInCmapDown(gene_id="348")
#'

checkExistenceInCmapDown = function(gene_id){
  gene_id = as.numeric(gene_id)
  assertthat::assert_that(assertthat::is.number(gene_id),msg = "the gene id is not included in the extended CMap")
  # 最大重试次数
  max_retries <- 10
  # 初始化重试计数器
  retries <- 0
  # 定义一个标志变量，用于控制循环
  success <- FALSE
  # 循环发送请求，直到成功或达到最大重试次数
  while (!success && retries < max_retries) {
    # 发送 GET 请求
    response_down <- httr::GET(paste0("https://ngdc.cncb.ac.cn/cedr/webcmap/api/cmapzscoredownk",K,"?gene_id=",gene_id))
    # 检查状态码
    if (httr::status_code(response_down) == 200) {
      success <- TRUE  # 请求成功，退出循环
    } else {
      retries <- retries + 1  # 增加重试次数
      #等待再重试
      Sys.sleep(retries)
    }
  }
  # 如果最终请求成功，处理响应数据
  if (success) {
    json_data <- rawToChar(response_down[["content"]])
    parsed_data <- jsonlite::fromJSON(json_data)
    if(!is.null(parsed_data$data)){
      return(parsed_data$data$sig_index)
    }else{
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' Retrieve CMap signature from remote server
#'
#' @param signature_index
#'
#' @return a \code{data.frame} object
#' @export
#' @examples
#'
#' K=50
#' retrieveCmapSignature(signature_index=519170)
#'

retrieveCmapSignature = function(signature_index){
  assertthat::assert_that(assertthat::is.number(signature_index), signature_index >= 1, signature_index <= 720216,msg = "the signature is not included in the extended CMap")
  # 最大重试次数
  max_retries <- 10
  # 初始化重试计数器
  retries <- 0
  # 定义一个标志变量，用于控制循环
  success <- FALSE
  # 循环发送请求，直到成功或达到最大重试次数
  while (!success && retries < max_retries) {
    # 发送 GET 请求
    response <- httr::GET(paste0("https://ngdc.cncb.ac.cn/cedr/webcmap/api/cmapcompound?sig_index=",signature_index))
    # 检查状态码
    if (httr::status_code(response) == 200) {
      success <- TRUE  # 请求成功，退出循环
    } else {
      retries <- retries + 1  # 增加重试次数
      #等待再重试
      Sys.sleep(retries)
    }
  }
  # 如果最终请求成功，处理响应数据
  if (success) {
    json_data <- rawToChar(response[["content"]])
    parsed_data <- jsonlite::fromJSON(json_data)
    cmap_signature <- do.call(c,jsonlite::fromJSON(parsed_data$data$signature))
    gene_id = names(cmap_signature)
    modz = as.numeric(cmap_signature)
    cmap_signature = data.frame(gene_id=gene_id,modz=modz)
    return(cmap_signature)
  } else {
    return(NULL)
  }
}

#' Retrieve CMap signature meta information from remote server
#'
#' @param signature_index
#'
#' @return a \code{data.frame} object
#' @export
#' @examples
#'
#' K=50
#' retrieveCmapSignatureMeta(signature_index=519170)
#'

retrieveCmapSignatureMeta = function(signature_index){
  assertthat::assert_that(assertthat::is.number(signature_index), signature_index >= 1, signature_index <= 720216,msg = "the signature is not included in the extended CMap")
  # 最大重试次数
  max_retries <- 10
  # 初始化重试计数器
  retries <- 0
  # 定义一个标志变量，用于控制循环
  success <- FALSE
  # 循环发送请求，直到成功或达到最大重试次数
  while (!success && retries < max_retries) {
    # 发送 GET 请求
    response <- httr::GET(paste0("https://ngdc.cncb.ac.cn/cedr/webcmap/api/cmapsignature?sig_index=",signature_index))
    # 检查状态码
    if (httr::status_code(response) == 200) {
      success <- TRUE  # 请求成功，退出循环
    } else {
      retries <- retries + 1  # 增加重试次数
      #等待再重试
      Sys.sleep(retries)
    }
  }
  # 如果最终请求成功，处理响应数据
  if (success) {
    json_data <- rawToChar(response[["content"]])
    parsed_data <- jsonlite::fromJSON(json_data)
    signature_meta <- parsed_data$data
    # #remove the first and seconde cols: id and sig_index
    # signature_meta <- signature_meta[,-c(1,2)]
    return(signature_meta)
  } else {
    return(NULL)
  }
}


#' retrieve the Permutation result from remote server
#'
#' @param signature_index the index of CMap signature
#' @param method XSum or CSS
#'
#' @return a \code{data.frame}
#' @export
#' @examples
#'
#' K=50
#' retrievePermutationResult(signature_index=519170,method="CSS")
#'

retrievePermutationResult = function(signature_index,method){
  assertthat::assert_that(assertthat::is.number(signature_index), signature_index >= 1, signature_index <= 720216,msg = "the signature is not included in the extended CMap")
  # 最大重试次数
  max_retries <- 10
  # 初始化重试计数器
  retries <- 0
  # 定义一个标志变量，用于控制循环
  success <- FALSE
  # 循环发送请求，直到成功或达到最大重试次数
  while (!success && retries < max_retries) {
    # 发送 GET 请求
    response <- httr::GET(paste0("https://ngdc.cncb.ac.cn/cedr/webcmap/api/permutation_",tolower(method),"_k",K,"?sig_index=",signature_index))
    # 检查状态码
    if (httr::status_code(response) == 200) {
      success <- TRUE  # 请求成功，退出循环
    } else {
      retries <- retries + 1  # 增加重试次数
      #等待再重试
      Sys.sleep(retries)
    }
  }
  # 如果最终请求成功，处理响应数据
  if (success) {
    json_data <- rawToChar(response[["content"]])
    parsed_data <- jsonlite::fromJSON(json_data)
    premutation_result <- do.call(c,jsonlite::fromJSON(parsed_data$data$permutation))
    permutation = names(premutation_result)
    score = as.numeric(premutation_result)
    premutation_result = data.frame(permutation=permutation,score=score)
    return(premutation_result)
  } else {
    return(NULL)
  }
}
