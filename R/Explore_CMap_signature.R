# library("cmapR")
# library(data.table)
# gct_file = "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_beta_trt_cp_n720216x12328.gctx"
# expanded_cmap = parse_gctx(gct_file)
# #fwrite(expanded_cmap@mat,file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/expanded_cmap_matrix.txt")
# #fwrite(data.table(rid=expanded_cmap@rid),file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/expanded_cmap_rid.txt")
# #fwrite(data.table(cid=expanded_cmap@cid),file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/expanded_cmap_cid.txt")
#
# # gct_file = "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/CMap/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
# # cmap = parse.gctx(gct_file)
#
#
# library(dplyr)
# library(data.table)
# compoundinfo = fread("/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/compoundinfo_beta.txt")
# geneinfo = fread("/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/geneinfo_beta.txt")
# siginfo = fread("/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/siginfo_beta.txt")
# cellinfo = fread("/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cellinfo_beta.txt")
#
# siginfo[which(siginfo$sig_id==expanded_cmap@cid[1]),]
# geneinfo[which(geneinfo$gene_id==expanded_cmap@rid[1]),]
#
# length(intersect(expanded_cmap@cid,siginfo$sig_id))
# length(intersect(expanded_cmap@rid,geneinfo$gene_id))
#
# #data.tabke::frank, rank 1 means the smallest
# frank(expanded_cmap@mat[,1]) + ranked@mat[,1]
#
# rm("expanded_cmap")
# #get the gene id of the first and bottom 50 genes
# cmp_signature_gene = apply(ranked@mat,MARGIN = 2,function(x){
# 	cmp_gene_order = order(x)
# 	cmp_up_gene = ranked@rid[head(cmp_gene_order,50)]
# 	cmp_down_gene = ranked@rid[tail(cmp_gene_order,50)]
# 	res = c(cmp_up_gene,cmp_down_gene)
# 	return(res)
# })
# cmp_signature_gene = as.data.table(cmp_signature_gene)
# colnames(cmp_signature_gene) = ranked@cid
# fwrite(cmp_signature_gene,file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmp_signature_gene_k50.txt")
#
# #get the gene id of the first and bottom 100 genes
# cmp_signature_gene = apply(ranked@mat,MARGIN = 2,function(x){
# 	cmp_gene_order = order(x)
# 	cmp_up_gene = ranked@rid[head(cmp_gene_order,100)]
# 	cmp_down_gene = ranked@rid[tail(cmp_gene_order,100)]
# 	res = c(cmp_up_gene,cmp_down_gene)
# 	return(res)
# })
# cmp_signature_gene = as.data.table(cmp_signature_gene)
# colnames(cmp_signature_gene) = ranked@cid
# fwrite(cmp_signature_gene,file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmp_signature_gene_k100.txt")
#
# #get the gene id of the first and bottom 200 genes
# cmp_signature_gene = apply(ranked@mat,MARGIN = 2,function(x){
# 	cmp_gene_order = order(x)
# 	cmp_up_gene = ranked@rid[head(cmp_gene_order,200)]
# 	cmp_down_gene = ranked@rid[tail(cmp_gene_order,200)]
# 	res = c(cmp_up_gene,cmp_down_gene)
# 	return(res)
# })
# cmp_signature_gene = as.data.table(cmp_signature_gene)
# colnames(cmp_signature_gene) = ranked@cid
# fwrite(cmp_signature_gene,file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmp_signature_gene_k200.txt")
#
# #get the gene id of the first and bottom 500 genes
# cmp_signature_gene = apply(ranked@mat,MARGIN = 2,function(x){
# 	cmp_gene_order = order(x)
# 	cmp_up_gene = ranked@rid[head(cmp_gene_order,500)]
# 	cmp_down_gene = ranked@rid[tail(cmp_gene_order,500)]
# 	res = c(cmp_up_gene,cmp_down_gene)
# 	return(res)
# })
# cmp_signature_gene = as.data.table(cmp_signature_gene)
# colnames(cmp_signature_gene) = ranked@cid
# fwrite(cmp_signature_gene,file="/xtdisk/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmp_signature_gene_k500.txt")
#
#
#
# ############
# library(DBI)
# library(RSQLite)
# library(data.table)
# all_cmap_up = data.frame()
# all_cmap_down = data.frame()
#
# K=50
#
# con1 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp1.db")
# con2 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp2.db")
# con3 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp3.db")
# con4 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp4.db")
# con5 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp5.db")
# con6 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp6.db")
# con7 = dbConnect(RSQLite::SQLite(), "/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/level5_trt_cp7.db")
#
# for (signature_index in 1:720216) {
#     print(signature_index)
#     db_name = ceiling(signature_index/100000)
#     if(db_name==1){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con1,signature_table_name)
#     }else if(db_name==2){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con2,signature_table_name)
#     }else if(db_name==3){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con3,signature_table_name)
#     }else if(db_name==4){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con4,signature_table_name)
#     }else if(db_name==5){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con5,signature_table_name)
#     }else if(db_name==6){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con6,signature_table_name)
#     }else if(db_name>=7){
#         signature_table_name = paste0("signature_",signature_index)
#         cmap_signature = dbReadTable(con7,signature_table_name)
#     }
#
#     cmap_signature$sig_index = signature_index
#     cmap_signature = cmap_signature[order(cmap_signature$modz,decreasing=T,na.last=NA),]
#     up = head(cmap_signature,n=K)
#     down = tail(cmap_signature,n=K)
#     all_cmap_up = rbind(all_cmap_up,up)
#     all_cmap_down = rbind(all_cmap_down,down)
# }
# fwrite(all_cmap_up,"/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/all_cmap_up_zscore.txt",spe="\t")
# fwrite(all_cmap_down,"/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/all_cmap_down_zscore.txt",spe="\t")
#
#
# #MERGE
# library(data.table)
# all_up = data.frame()
# all_down = data.frame()
# for(i in 1:15){
#         file_name = paste0("/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmap_modz_table/all_cmap_up_zscore_",i,".txt")
#         a = fread(file_name)
#         all_up = rbind(all_up,a)
#
#         file_name2 = paste0("/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmap_modz_table/all_cmap_down_zscore_",i,".txt")
#         b = fread(file_name2)
#         all_down = rbind(all_down,b)
# }
# library(org.Hs.eg.db)
# ENTREZtoSYMBOL <- function(x){
#   vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENTREZID",multiVals="first")
#   repl=which(is.na(vec))
#   for(y in repl){vec[y]=x[y]}
#   return(vec)
# }
#
# all_up$gene_name = ENTREZtoSYMBOL(as.character(all_up$gene_id))
# all_down$gene_name = ENTREZtoSYMBOL(as.character(all_down$gene_id))
#
# fwrite(all_up,file="/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmap_modz_table/merged_cmap_up_zscore.txt",sep="\t")
# fwrite(all_down,file="/p300s/jiapl_group/kanghongen/projects/GWAS_Drug/Expanded_CMap_2020/cmap_modz_table/merged_cmap_down_zscore.txt",sep="\t")
