### Over representation analysis -- gProfiler2
gProfiler2_overrepresentation_analysis <- function (significant_genes_FC_ordered, local_path, comparison, filename = '', p_value_threshold = 0.05) {

    significant_genes_FC_ordered <- list(comparison = significant_genes_FC_ordered)
    # Enrichment pathway analysis
    enrichment_results <- gost(query = significant_genes_FC_ordered, 
                    organism = "mmusculus", ordered_query = TRUE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = p_value_threshold, correction_method = "g_SCS", 
                    domain_scope = "annotated", custom_bg = NULL, 
                    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)

    if (!(is.null(enrichment_results[['result']]) )) {

     # Manhattan plot of results
     plot  <- gostplot(enrichment_results, capped = T, interactive=T)
     htmlwidgets::saveWidget((plot), paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2', '.html'))

     # Saving results tablec
     enrichment_results <- enrichment_results[['result']] |> dplyr::select(-c('parents'))
     write.csv(enrichment_results, paste0(local_path, filename, 'Pathway_enrichment_analysis_gprofiler2', '.csv')) 
    }
    
}

gProfiler2_functional_analysis <- function (results,  cluster, comparison, path='./', FC_threshold, p_value_threshold = 0.05) {

     results <- results[which(duplicated(results$genes) == F),]


####################################### UP ########################################

     local_path <- here(path, paste0('Functional_analysis_UP_gProfiler2_', cluster, '_', comparison, '/'))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange > FC_threshold)) |> arrange(padj) |> arrange(desc(log2FoldChange)) |> pull(genes)
     
     if (length(significant_genes) > 2) {
          gProfiler2_overrepresentation_analysis(significant_genes, local_path =  local_path , comparison = comparison, filename = '', p_value_threshold = p_value_threshold)

     }

######################################## DOWN ########################################

     local_path <- here(path, paste0('Functional_analysis_DOWN_gProfiler2_', cluster, '_', comparison, '/'))
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < p_value_threshold) & (log2FoldChange < -1*FC_threshold)) |> arrange(padj) |> arrange(log2FoldChange) |> pull(genes)

     if (length(significant_genes) > 2) {
          gProfiler2_overrepresentation_analysis(significant_genes, local_path =  local_path , comparison = comparison, filename = '', p_value_threshold = p_value_threshold)
     }

     return()
}

gProfiler2_functional_analysis_cluster_identification <- function (scRNAseq, results, identities = 'seurat_clusters', path='./') {

    color_scale <- viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)
     local_path <- paste0(path, 'Cluster_identification_functional_analysis_gProfiler2/')
     unlink(local_path, recursive = T)
     dir.create(local_path)

     all_genes <- Features(scRNAseq[['RNA']]) |>unique()
     #all_genes_entrezid <- Features(scRNAseq[['RNA']]) |>unique() |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique()

     #mouse_database <- msigdbr(species = 'Mus musculus',category = 'C8') |> dplyr::select(gs_name, entrez_gene)

for (cluster in levels(scRNAseq@meta.data |> pull(!!identities))) {


     name <- paste0('Cluster ', cluster, ' - ')

     significant_results_cluster <- results |> 
          filter(cluster == cluster) |>
          pull(gene) |>
          unique()
    
     gProfiler2_overrepresentation_analysis(significant_results_cluster, local_path, comparison = 'clusters', filename =  name )
 

     }
     return()
}