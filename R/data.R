#' Example query signature
#'
#' An example gene signatures for demonstration.
#'
#' @format A data frame with 11,689 rows and 2 variables:
#' \describe{
#'   \item{gene_id}{gene Entrez ID, character}
#'   \item{zscore}{derived from comparing perturbation and control conditions to quantify both the direction and magnitude of transcriptional changes, numeric}
#' }
#' @usage data(query_signature)
#' @source Generated internally
"query_signature"

#' Example query signature
#'
#' An example gene signatures for demonstration.
#'
#' @format A data frame with 11,689 rows and 2 variables:
#' \describe{
#'   \item{gene_id}{gene Entrez ID, character}
#'   \item{zscore}{derived from comparing perturbation and control conditions to quantify both the direction and magnitude of transcriptional changes, numeric}
#' }
#' @usage data(query_signature1)
#' @source Generated internally
"query_signature1"



#' Example query signature
#'
#' An example gene signatures for demonstration.
#'
#' @format A data frame with 12,328 rows and 2 variables:
#' \describe{
#'   \item{gene_id}{gene Entrez ID, character}
#'   \item{zscore}{derived from comparing perturbation and control conditions to quantify both the direction and magnitude of transcriptional changes, numeric}
#' }
#' @usage data(query_signature2)
#' @source Generated internally
"query_signature2"


#' Example result of run_negative_CMap
#'
#' An example result of run_negative_CMap for demonstration.
#'
#' @format A data frame with 601 rows and 56 variables:
#' \describe{
#'   \item{signature_index }{CMap signature_index, numeric}
#'   \item{Meta_score}{Meta_score, numeric}
#'   \item{WTCS}{Weighted Connectivity Score, numeric}
#'   \item{ES_up}{Enrichment score for the upregulated query genes, numeric}
#'   \item{ES_down}{Enrichment score for the downregulated query genes, numeric}
#'   \item{ES_up_padj}{the adjusted p-value for the enrichment score of the upregulated genes in GSEA, numeric}
#'   \item{ES_down_padj}{the adjusted p-value for the enrichment score of the downregulated genes in GSEA, numeric}
#'   \item{XSum}{eXtreme Sum score, numeric}
#'   \item{XSum_pvalue}{the p value of eXtreme Sum score, numeric}
#'   \item{CSS}{Connection Strength Score, numeric}
#'   \item{CSS_pvalue}{The p value of Connection Strength Score, numeric}
#'   \item{Spearman_correlation}{Spearman correlation coefficients, numeric}
#'   \item{Pearson_correlation}{Pearson correlation coefficients, numeric}
#'   \item{Consine_similarity}{Consine similarity Score, numeric}
#'   \item{pert_id}{A unique identifier for a perturbagen that refers to the perturbagen in general, not to any particular batch or sample, character}
#'   \item{bead_batch}{One instantiation of a complete set of beads, which have been coupled to probes at one time, under the same conditions, character}
#'   \item{nearest_dose}{Rounded version of the dose, mapped to the nearest dose in a predefined series, numeric}
#'   \item{pert_dose}{Precise dose used in the experiment, numeric}
#'   \item{pert_dose_unit}{Units of pert_dose, character}
#'   \item{pert_idose}{The concatenation of pert_dose and pert_dose_unit to create a string containing the dose information. We use a standardized dose for a perturbagen treatment. For example, the less common dose of 10.04 is rounded to 10. This enables grouping of signatures by a common dose., character}
#'   \item{pert_itime}{The length of time, expressed as a number, that a perturbagen was applied to the cells; does not include the unit., character}
#'   \item{pert_time}{The concatenation of pert_time and pert_time_unit to create a string containing the length of time that a perturbagen was applied to the cells. We use a standardized time for a perturbagen treatment. For example, if data is made by treating cells with a perturbagen for 5.5 hours, we round that time to the more common treatment time of 6 hours, numeric}
#'   \item{pert_time_unit}{The unit that applies to the pert_time numerical value, character}
#'   \item{cell_mfc_name}{The supplier's name for the cell line. This will most often be equivalent to the cell_iname, character}
#'   \item{pert_mfc_id}{A manufacturer's id for the perturbagen; by convention, for compounds registered with Broad Compound Management this is the full BRD containing both the compound ID and the batch ID, character}
#'   \item{nsample}{Number of individual replicate profiles (level 4 / z-score) that were used to create the signature (level 5 / aggregate z-score), numeric}
#'   \item{cc_q75}{75th quantile of pairwise spearman correlations in landmark space of replicate level 4 profiles, numeric}
#'   \item{ss_ngene}{The number of significantly differentially expressed transcripts that arise from a particular perturbagen treatment. Defined as the number of genes in the signature with |replicate-adjusted modz| >= 2, numeric}
#'   \item{tas}{Transcriptional activity score - computed as the geometric mean of ss_ngene and cc_q75, normalized to range between 0 and 1, numeric}
#'   \item{pct_self_rank_q25}{Self connectivity of replicates expressed as a percentage of total instances in a replicate set, numeric}
#'   \item{wt}{Comma-delimited list of the replicate weightings used to collapse into the level 5 signature. These will sum to 1, character}
#'   \item{median_recall_rank_spearman}{, character}
#'   \item{median_recall_rank_spearman}{The median pairwise recall rank by spearman correlation between replicate level 4 profiles, numeric}
#'   \item{median_recall_rank_wtcs_50}{The median pairwise recall rank by weighted connectivity score between replicate level 4 profiles, numeric}
#'   \item{median_recall_score_spearman}{The median pairwise spearman correlation between replicate level 4 profiles, numeric}
#'   \item{median_recall_score_wtcs_50}{The median pairwise weighted connectiity score between replicate level 4 profiles, numeric}
#'   \item{batch_effect_tstat}{The result of a one-sample t-test comparing the median correlation between signatures from the same plate to 0. Higher values indicate higher correlations between such signatures, and hence potentially a batch effect, numeric}
#'   \item{batch_effect_tstat_pct}{The percentile of batch_effect_tstat relative to all other signatures, numeric}
#'   \item{is_hiq}{Binary indicating whether the given signature was of high techincal and functional quality. Specific requirements are qc_pass == 1 AND (median_recall_rank_spearman <= 5 OR median_recall_rank_wtcs_50 <= 5), character}
#'   \item{qc_pass}{Binary indicating whether the given signature had at least 50% of its replicates flagged as qc_pass. See definition of qc_pass in instinfo for more details, character}
#'   \item{sig_id}{Unique identifier for the signature, character}
#'   \item{pert_type}{Abbreviated designation for perturbagen type, referring to compound or genetic perturbagens that are used in cell treatments to assess gene expression effects. The various pert_types used by CMap are listed in the pert_type sheet., character}
#'   \item{cell_iname}{Curated name for the cell line, character}
#'   \item{det_wells}{Pipe-delimited list of detection wells, which refers to each well of the detection plate in which an L1000 experiment is conducted, character}
#'   \item{det_plates}{Pipe-delimited list of detection plates, the plate of L1000 experiments that, at the end of the assay pipeline, is put through the Luminex scanners to detect the levels of landmark gene amplicons, character}
#'   \item{distil_ids}{Pipe-delmited list of IDs of individual replicate profiles, referred to as level 4 / z-score data, that is used in creating the signature from replicates assayed together on an L1000 plate. The signature is referred to as level 5 / aggregated z-score data, character}
#'   \item{project_code}{An internal code identifying the project to which a signature belongs, character}
#'   \item{cmap_name}{The internal (CMap-designated) name of a perturbagen. By convention, for genetic perturbations CMap uses the HUGO gene symbol, character}
#'   \item{is_exemplar_sig}{is_exemplar_sig, character}
#'   \item{is_ncs_sig}{is_ncs_sig, character}
#'   \item{is_null_sig}{is_null_sig, character}
#'   \item{target}{The symbol of the gene that the comopund targets, character}
#'   \item{moa}{A curated phrase representing the compound's mechanism of action, character}
#'   \item{canonical_smiles}{Canonical SMILES structure, character}
#'   \item{inchi_key}{InChIKey - hashed version of the InChi identifier (see http://www.iupac.org/inchi/), character}
#'   \item{compound_aliases}{Alternative name for the compound, character}
#'
#' }
#' @usage data(res1)
#' @source Generated internally
"res1"


#' Example result of run_negative_CMap
#'
#' An example result of run_positive_CMap for demonstration.
#'
#' @format A data frame with 923 rows and 56 variables:
#' \describe{
#'   \item{signature_index }{CMap signature_index, numeric}
#'   \item{Meta_score}{Meta_score, numeric}
#'   \item{WTCS}{Weighted Connectivity Score, numeric}
#'   \item{ES_up}{Enrichment score for the upregulated query genes, numeric}
#'   \item{ES_down}{Enrichment score for the downregulated query genes, numeric}
#'   \item{ES_up_padj}{the adjusted p-value for the enrichment score of the upregulated genes in GSEA, numeric}
#'   \item{ES_down_padj}{the adjusted p-value for the enrichment score of the downregulated genes in GSEA, numeric}
#'   \item{XSum}{eXtreme Sum score, numeric}
#'   \item{XSum_pvalue}{the p value of eXtreme Sum score, numeric}
#'   \item{CSS}{Connection Strength Score, numeric}
#'   \item{CSS_pvalue}{The p value of Connection Strength Score, numeric}
#'   \item{Spearman_correlation}{Spearman correlation coefficients, numeric}
#'   \item{Pearson_correlation}{Pearson correlation coefficients, numeric}
#'   \item{Consine_similarity}{Consine similarity Score, numeric}
#'   \item{pert_id}{A unique identifier for a perturbagen that refers to the perturbagen in general, not to any particular batch or sample, character}
#'   \item{bead_batch}{One instantiation of a complete set of beads, which have been coupled to probes at one time, under the same conditions, character}
#'   \item{nearest_dose}{Rounded version of the dose, mapped to the nearest dose in a predefined series, numeric}
#'   \item{pert_dose}{Precise dose used in the experiment, numeric}
#'   \item{pert_dose_unit}{Units of pert_dose, character}
#'   \item{pert_idose}{The concatenation of pert_dose and pert_dose_unit to create a string containing the dose information. We use a standardized dose for a perturbagen treatment. For example, the less common dose of 10.04 is rounded to 10. This enables grouping of signatures by a common dose., character}
#'   \item{pert_itime}{The length of time, expressed as a number, that a perturbagen was applied to the cells; does not include the unit., character}
#'   \item{pert_time}{The concatenation of pert_time and pert_time_unit to create a string containing the length of time that a perturbagen was applied to the cells. We use a standardized time for a perturbagen treatment. For example, if data is made by treating cells with a perturbagen for 5.5 hours, we round that time to the more common treatment time of 6 hours, numeric}
#'   \item{pert_time_unit}{The unit that applies to the pert_time numerical value, character}
#'   \item{cell_mfc_name}{The supplier's name for the cell line. This will most often be equivalent to the cell_iname, character}
#'   \item{pert_mfc_id}{A manufacturer's id for the perturbagen; by convention, for compounds registered with Broad Compound Management this is the full BRD containing both the compound ID and the batch ID, character}
#'   \item{nsample}{Number of individual replicate profiles (level 4 / z-score) that were used to create the signature (level 5 / aggregate z-score), numeric}
#'   \item{cc_q75}{75th quantile of pairwise spearman correlations in landmark space of replicate level 4 profiles, numeric}
#'   \item{ss_ngene}{The number of significantly differentially expressed transcripts that arise from a particular perturbagen treatment. Defined as the number of genes in the signature with |replicate-adjusted modz| >= 2, numeric}
#'   \item{tas}{Transcriptional activity score - computed as the geometric mean of ss_ngene and cc_q75, normalized to range between 0 and 1, numeric}
#'   \item{pct_self_rank_q25}{Self connectivity of replicates expressed as a percentage of total instances in a replicate set, numeric}
#'   \item{wt}{Comma-delimited list of the replicate weightings used to collapse into the level 5 signature. These will sum to 1, character}
#'   \item{median_recall_rank_spearman}{, character}
#'   \item{median_recall_rank_spearman}{The median pairwise recall rank by spearman correlation between replicate level 4 profiles, numeric}
#'   \item{median_recall_rank_wtcs_50}{The median pairwise recall rank by weighted connectivity score between replicate level 4 profiles, numeric}
#'   \item{median_recall_score_spearman}{The median pairwise spearman correlation between replicate level 4 profiles, numeric}
#'   \item{median_recall_score_wtcs_50}{The median pairwise weighted connectiity score between replicate level 4 profiles, numeric}
#'   \item{batch_effect_tstat}{The result of a one-sample t-test comparing the median correlation between signatures from the same plate to 0. Higher values indicate higher correlations between such signatures, and hence potentially a batch effect, numeric}
#'   \item{batch_effect_tstat_pct}{The percentile of batch_effect_tstat relative to all other signatures, numeric}
#'   \item{is_hiq}{Binary indicating whether the given signature was of high techincal and functional quality. Specific requirements are qc_pass == 1 AND (median_recall_rank_spearman <= 5 OR median_recall_rank_wtcs_50 <= 5), character}
#'   \item{qc_pass}{Binary indicating whether the given signature had at least 50% of its replicates flagged as qc_pass. See definition of qc_pass in instinfo for more details, character}
#'   \item{sig_id}{Unique identifier for the signature, character}
#'   \item{pert_type}{Abbreviated designation for perturbagen type, referring to compound or genetic perturbagens that are used in cell treatments to assess gene expression effects. The various pert_types used by CMap are listed in the pert_type sheet., character}
#'   \item{cell_iname}{Curated name for the cell line, character}
#'   \item{det_wells}{Pipe-delimited list of detection wells, which refers to each well of the detection plate in which an L1000 experiment is conducted, character}
#'   \item{det_plates}{Pipe-delimited list of detection plates, the plate of L1000 experiments that, at the end of the assay pipeline, is put through the Luminex scanners to detect the levels of landmark gene amplicons, character}
#'   \item{distil_ids}{Pipe-delmited list of IDs of individual replicate profiles, referred to as level 4 / z-score data, that is used in creating the signature from replicates assayed together on an L1000 plate. The signature is referred to as level 5 / aggregated z-score data, character}
#'   \item{project_code}{An internal code identifying the project to which a signature belongs, character}
#'   \item{cmap_name}{The internal (CMap-designated) name of a perturbagen. By convention, for genetic perturbations CMap uses the HUGO gene symbol, character}
#'   \item{is_exemplar_sig}{is_exemplar_sig, character}
#'   \item{is_ncs_sig}{is_ncs_sig, character}
#'   \item{is_null_sig}{is_null_sig, character}
#'   \item{target}{The symbol of the gene that the comopund targets, character}
#'   \item{moa}{A curated phrase representing the compound's mechanism of action, character}
#'   \item{canonical_smiles}{Canonical SMILES structure, character}
#'   \item{inchi_key}{InChIKey - hashed version of the InChi identifier (see http://www.iupac.org/inchi/), character}
#'   \item{compound_aliases}{Alternative name for the compound, character}
#'
#' }
#' @usage data(res2)
#' @source Generated internally
"res2"
