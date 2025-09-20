#' radar plot
#'
#' @param data the selected connectivity analysis result \code{data.frame} to plot
#'
#' @return The ggplot object
#' @export
#' @examples
#' data(res1)
#' data = res1[1,]
#' radar_plot(data)
radar_plot <- function(data) {
  assertthat::assert_that(nrow(data)==1)
  dat = data.frame(
    WTCS = data$WTCS,
    XSum = data$XSum,
    CSS = data$CSS,
    Spearman_correlation = data$Spearman_correlation,
    Pearson_correlation = data$Pearson_correlation,
    Consine_similarity = data$Consine_similarity
  )

  rownames(dat) <- data$signature_index

  # 添加最大和最小行（需根据实际范围调整）
  max_min <- data.frame(
    WTCS = c(1,-1),
    XSum = c(25,-25),
    CSS = c(1,-1),
    Spearman_correlation = c(1,-1),
    Pearson_correlation = c(1,-1),
    Consine_similarity = c(1,-1)
  )
  plot_data <- rbind(max_min, dat)

  # 绘制雷达图
  p = fmsb::radarchart(
    plot_data,
    title = paste0("CMap Name: ", data$cmap_name),
    pcol = "#2D336B", # 数据线颜色
    plwd = 2, # 数据线宽度
    pfcol = rgb(0, 0, 1, 0.3), # 填充颜色
    cglcol = "grey", # 网格线颜色
    cglty = 1, # 网格线类型
    axislabcol = "grey", # 轴标签颜色
    vlcex = 1.2          # 变量标签字号
  )

  # Add a legend
  #legend(x=0.8, y=1.2, legend = paste0("Signature Index: ",data$signature), bty = "n", pch=20 , col="#2D336B" , cex=1.2, pt.cex=3)

  return(p)
}


#' gsea plot
#'
#' @param data the selected connectivity analysis result \code{data.frame} to plot
#' @param query_signature A data frame containing gene signatures
#' @param K Integer, number of extreme genes
#' @return The ggplot object
#' @export
#' @examples
#' data(res1)
#' data(query_signature1)
#' data = res1[1,]
#' gsea_plot(data, query_signature1, K = 50)
#'
gsea_plot <- function(data, query_signature, K = 50) {
  assertthat::assert_that(nrow(data)==1)
  cmap_signature = retrieveCmapSignature(as.numeric(data$signature_index))
  cmap_ranks = cmap_signature$modz
  names(cmap_ranks) = cmap_signature$gene_id
  query_signature = query_signature[order(query_signature$zscore,
                                          decreasing = TRUE,
                                          na.last = NA), ]
  query_up_gene = head(query_signature$gene_id, K)
  query_down_gene = tail(query_signature$gene_id, K)
  query_pathway = list()
  query_pathway[["query_up_gene"]] = query_up_gene
  query_pathway[["query_down_gene"]] = query_down_gene
  # fgseaRes = fgsea::fgsea(query_pathway,cmap_ranks,gseaParam = 1)
  # ES_up = fgseaRes$ES[which(fgseaRes$pathway=="query_up_gene")]
  # ES_up_padj = fgseaRes$padj[which(fgseaRes$pathway=="query_up_gene")]
  # ES_down = fgseaRes$ES[which(fgseaRes$pathway=="query_down_gene")]
  # ES_down_padj = fgseaRes$padj[which(fgseaRes$pathway=="query_down_gene")]

  p1 = fgsea::plotEnrichment(query_pathway[["query_up_gene"]], cmap_ranks) + ggplot2::labs(title =
                                                                                             "Up Regualted Genes of Query Signature", ) + ggplot2::xlab("Rank of CMap Signature") + ggplot2::ylab("Enrichment Score")
  p2 = fgsea::plotEnrichment(query_pathway[["query_down_gene"]], cmap_ranks) + ggplot2::labs(title =
                                                                                               "Down Regualted Genes of Query Signature") + ggplot2::xlab("Rank of CMap Signature") + ggplot2::ylab("Enrichment Score")
  p = ggpubr::ggarrange(p1, p2)
  return(p)
}

#' venn plot
#'
#' @param data the selected connectivity analysis result \code{data.frame} to plot
#' @param query_signature A data frame containing gene signatures
#' @param K Integer, number of extreme genes
#' @return The ggplot object
#' @export
#' @examples
#' data(res1)
#' data(query_signature1)
#' data = res1[1,]
#' venn_plot(data, query_signature1, K=50)
venn_plot <- function(data, query_signature, K = 50) {
  assertthat::assert_that(nrow(data)==1)
  cmap_signature = retrieveCmapSignature(as.numeric(data$signature_index))
  cmap_signature = cmap_signature[order(cmap_signature$modz,
                                        decreasing = TRUE,
                                        na.last = NA), ]
  cmap_up_gene = head(cmap_signature$gene_id, K)
  cmap_down_gene = tail(cmap_signature$gene_id, K)
  query_signature = query_signature[order(query_signature$zscore,
                                          decreasing = TRUE,
                                          na.last = NA), ]
  query_up_gene = head(query_signature$gene_id, K)
  query_down_gene = tail(query_signature$gene_id, K)
  plot_data <- list(
    `CMap Up Gene` = cmap_up_gene,
    `CMap Down Gene` = cmap_down_gene,
    `Query Up Gene` = query_up_gene,
    `Query Down Gene` = query_down_gene
  )
  if (data$WTCS < 0) {
    p1 = ggvenn::ggvenn(plot_data,
                        c("CMap Up Gene", "Query Down Gene"),
                        fill_color = c("#E195AB", "#7886C7"))
    p2 = ggvenn::ggvenn(plot_data,
                        c("CMap Down Gene", "Query Up Gene"),
                        fill_color = c("#7886C7", "#E195AB"))
    p = ggpubr::ggarrange(p1, p2)
    return(p)
  } else{
    p1 = ggvenn::ggvenn(plot_data,
                        c("CMap Up Gene", "Query Up Gene"),
                        fill_color = c("#FFC785", "#48A6A7"))
    p2 = ggvenn::ggvenn(
      plot_data,
      c("CMap Down Gene", "Query Down Gene"),
      fill_color = c("#FFC785", "#48A6A7")
    )
    p = ggpubr::ggarrange(p1, p2)
    return(p)
  }
}
