## Initialization
# rm(list = ls())
# #Library the required

# suppressMessages(library(Seurat))
# suppressMessages(library(dplyr))
# suppressMessages(library(cowplot))
# suppressMessages(library(ggplot2))

# suppressMessages(library(stringr))
# suppressMessages(library(Stat2Data))
# suppressMessages(library(tidyverse))
# suppressMessages(library(patchwork))
# # library(pheatmap))
# suppressMessages(library(viridis)) 
# suppressMessages(library(ggplot2))
# suppressMessages(library(ggrepel))
# suppressMessages(library(scRepertoire))
# suppressMessages(library(circlize))
# suppressMessages(library(scales))
# suppressMessages(library(scCustomize))
# suppressMessages(library(DESeq2)) 
# suppressMessages(library(celldex))
# suppressMessages(library(SingleR))
# suppressMessages(library(gridExtra))

# suppressMessages(library(DOSE))
# suppressMessages(library(pathview))
# suppressMessages(library(clusterProfiler))
# suppressMessages(library(org.Mm.eg.db))
# suppressMessages(library(enrichplot))
# suppressMessages(library(msigdbr))
# suppressMessages(library(pheatmap))
# suppressMessages(library(UCell))
# suppressMessages(library(VennDiagram))
# suppressMessages(library(gprofiler2))

# sessionInfo()


## Functions
### Overrepresentation analysis function GO and MSigDB -- ClusterProfiler
GO_overrepresentation_analysis <- function (significant_genes, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 500, filename = '') {

     color_scale <- viridis(n = 4, direction = -1)
     options(enrichplot.colours = color_scale)
     
     enrichment_results <- enrichGO(gene = significant_genes, 
                    universe = all_genes,
                    keyType = "SYMBOL",
                    OrgDb = org.Mm.eg.db, 
                    ont = ontology, 
                    pAdjustMethod = "BH", 
                    minGSSize    = minGSSize,
                    maxGSSize    = maxGSSize,
                    qvalueCutoff = 0.25)

     enrichment_results_table <- as_tibble(enrichment_results)
     write.csv(enrichment_results_table, paste0(local_path, filename,'GO_OverRepresentation_analysis_results_', ontology, '.csv'))

     if (nrow(enrichment_results_table) > 1) {
          ## Add similarity matrix to the termsim slot of enrichment result
          enrichment_results <- enrichplot::pairwise_termsim(enrichment_results)
          enrichment_results <- simplify(enrichment_results, cutoff=0.7, by="p.adjust", select_fun=min)

          dotplot(enrichment_results,
               showCategory=50,
               title = paste0(filename,'GO Overrepresentation analysis ', ontology),
               label_format = 60)
          ggsave(paste0(filename, 'GO overrepresentation_analysis_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)

          ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
          try(emapplot(enrichment_results, showCategory = (nrow(enrichment_results_table-1))) + ggtitle(paste0(filename, 'Overrepresentation analysis ', ontology)))
          ggsave(paste0(filename, 'GO_overrepresentation_analysis_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)
     }
}

GO_GSEA_analysis <- function (results, local_path, ontology = 'ALL') {
     
     color_scale <- viridis(n = 4, direction = -1)
     options(enrichplot.colours = color_scale)

          #### GSEA ####

     fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
     names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(genes)

     gsea_results <- gseGO(geneList     = fold_changes,
               OrgDb        = org.Mm.eg.db,
               ont          = ontology,
               keyType = "SYMBOL",
               minGSSize    = 5,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)               

     gsea_results_table <- as_tibble(gsea_results) 
     write.csv(gsea_results_table, paste0(local_path, 'Gene_Set_Enrichment_Analysis_results_', ontology,'.csv'))

     if (nrow(gsea_results_table) > 1) {
          dotplot(gsea_results,
               showCategory=50,
               title = paste0('GSEA analysis ', ontology),
               label_format = 60)
          ggsave(paste0('GSEA_dotplot_', ontology,'.pdf'), width = 10, height = 18, path = local_path)

          ## Add similarity matrix to thenes,  termsim slot of enrichment result
          gsea_results <- enrichplot::pairwise_termsim(gsea_results)

          ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
        try(emapplot(gsea_results, showCategory = 60) + ggtitle(paste0('GSEA analysis ', ontology)))
          ggsave(paste0('GSEA_network_', ontology,'.pdf'), width = 14, height = 18, path = local_path)

## Make individual GSEA plots
          local_path2 <- paste0(local_path, 'Individual_GSEA_plot_', ontology, '/')
          unlink(local_path2, recursive = T)
          dir.create(local_path2)

          for (pathway in head(gsea_results$ID, 20)) {
               pathway_name <- gsea_results |> as_tibble() |> filter(ID == pathway) |> pull(Description)
               anno <- gsea_results[pathway, c("NES", "pvalue", "p.adjust")]
               lab <- paste0(names(anno), "=",  round(anno, 4), collapse="\n")

               try(p1 <- enrichplot::gseaplot2(gsea_results, geneSetID = pathway, pvalue_table = FALSE, subplots = 1, base_size = 13,title = paste0(pathway, ' ', pathway_name)))

               try(x_position <- ggplot_build(p1)$layout$panel_params[[1]]$x.range[2]*0.75)
               try(y_position <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-(ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-ggplot_build(p1)$layout$panel_params[[1]]$y.range[1])*0.17)

               try(p1 <- p1 + annotate("text", x_position, y_position, label = lab, hjust=0, vjust=0, size = 5))
               try(p2 <- enrichplot::gseaplot2(gsea_results, geneSetID = pathway, pvalue_table = FALSE, subplots = 2, base_size = 13))
               try(p3 <- enrichplot::gseaplot2(gsea_results, geneSetID = pathway, pvalue_table = FALSE, subplots = 3, base_size = 13))

               try(cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.5, 0.5, 1), align = 'v'))    
               try(ggsave(paste0('GSEA ',pathway_name , '.pdf'), path = local_path2, height = 10, width = 8))
               }
          }
}

GO_functional_analysis <- function (results,  cluster, path='./') {

     results <- results[which(duplicated(results$genes) == F),]
#    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    #results <- results[which(duplicated(results$entrezid) == F),]


####################################### UP ########################################

     local_path <- paste0(path, 'GO_functional_analysis_UP_', cluster, '/')
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < 0.05) & (log2FoldChange > 0)) |> arrange(padj) |> pull(genes)
     all_genes <- results |> arrange(padj) |> pull(genes)

     ######################## ORA ########################

     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'ALL')
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'BP')
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'MF')
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'CC')

     #################### GSEA ####################

     GO_GSEA_analysis(results, local_path, ontology = 'ALL')
     GO_GSEA_analysis(results, local_path, ontology = 'BP')
     GO_GSEA_analysis(results, local_path, ontology = 'MF')
     GO_GSEA_analysis(results, local_path, ontology = 'CC')

######################################## DOWN ########################################

     local_path <- paste0(path, 'GO_functional_analysis_DOWN_', cluster, '/')
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < 0.05) & (log2FoldChange < 0)) |> arrange(padj) |> pull(genes)
     all_genes <- results |> arrange(padj) |> pull(genes)

     ######################## ORA ########################

     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'ALL')
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'BP')
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'MF')
     GO_overrepresentation_analysis(significant_genes, all_genes, local_path, ontology = 'CC')

     #################### GSEA ####################

     GO_GSEA_analysis(results, local_path, ontology = 'ALL')
     GO_GSEA_analysis(results, local_path, ontology = 'BP')
     GO_GSEA_analysis(results, local_path, ontology = 'MF')
     GO_GSEA_analysis(results, local_path, ontology = 'CC')

     return()
}




GO_functional_analysis_cluster_identification <- function (scRNAseq, results, identities = 'seurat_clusters', path='./') {

    color_scale <- viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)
     local_path <- paste0(path, 'Cluster_identification_functional_analysis/')
     unlink(local_path, recursive = T)
     dir.create(local_path)

     all_genes <- Features(scRNAseq[['RNA']]) |>unique()
     all_genes_entrezid <- Features(scRNAseq[['RNA']]) |>unique() |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique()

     mouse_database <- msigdbr(species = 'Mus musculus',category = 'C8') |> dplyr::select(gs_name, entrez_gene)

for (cluster in unique(scRNAseq@meta.data |> pull(!!identities))) {

####################################### GO ########################################

     name <- paste0('Cluster ', cluster, ' - ')

     significant_results_cluster <- results |> 
          filter(cluster == cluster) |>
          pull(gene) |>
          unique()

     GO_overrepresentation_analysis(significant_results_cluster, all_genes, local_path, ontology = 'ALL', minGSSize = 5, maxGSSize = 2000, filename = name)  
 
#################### msigdbr ####################

     significant_results_cluster  <-  significant_results_cluster |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) |> unique() 
     enrichment_results <- enricher(gene = significant_results_cluster, 
                    universe = all_genes_entrezid,
                    #keyType = "ENTREZID",
                    #OrgDb = org.Mm.eg.db, 
                    pAdjustMethod = "BH", 
                    minGSSize    = 5,
                    maxGSSize    = 2000,
                    qvalueCutoff = 0.05,
                    TERM2GENE = mouse_database)
                        
        enrichment_results_table <- as_tibble(enrichment_results)
        write.csv(enrichment_results_table, paste0(local_path, name,'MSigDbr_OverRepresentation_analysis_results_', 'C8','.csv'))

        if (nrow(enrichment_results_table) > 0) {
            ## Add similarity matrix to the termsim slot of enrichment result
            enrichment_results <- enrichplot::pairwise_termsim(enrichment_results)
            
            enrichment_results_table <- as_tibble(enrichment_results) 
            write.csv(enrichment_results_table, paste0(local_path, 'MSigDbr_OverRepresentation_analysis_results_', 'C8','.csv'))

            dotplot(enrichment_results,
                showCategory=50,
                title = paste0('MSigDbr Overrepresentation analysis ', 'C8'),
                label_format = 60)
            ggsave(paste0(name,'MSigDbr_overrepresentation_analysis_dotplot_', 'C8','.pdf'), width = 10, height = 18, path = local_path)

            ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
            emapplot(enrichment_results, showCategory = 60) + ggtitle(paste0('Overrepresentation analysis ', 'C8'))
            ggsave(paste0(name,'MSigDbr_overrepresentation_analysis_network_', 'C8','.pdf'), width = 14, height = 18, path = local_path)
        }
     
}
     return()
}
msigdbr_functional_analysis <- function (results,cluster,  path='./') { 
    
    results <- results[which(duplicated(results$genes) == F),]
    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    results <- results[which(duplicated(results$entrezid) == F),] |>
        drop_na(entrezid)

        
    local_path <- paste0(path, 'MSigDb_functional_analysis_', cluster,'/')
    unlink(local_path, recursive = T)
    dir.create(local_path)

    color_scale <- viridis(n = 4, direction = -1)
    options(enrichplot.colours = color_scale)

    for (database in c('H', 'C2', 'C3', 'C8')) {

        mouse_database <- msigdbr(species = 'Mus musculus',category = database) |> dplyr::select(gs_name, entrez_gene)

        #### ORA ####

        significant_genes <- filter(results, padj < 0.05) |> arrange(padj) |> pull(entrezid)
        all_genes <- results |> arrange(padj) |> pull(entrezid)

        enrichment_results <- enricher(gene = significant_genes, 
                        universe = all_genes,
                        #keyType = "ENTREZID",
                        #OrgDb = org.Mm.eg.db, 
                        pAdjustMethod = "BH", 
                        minGSSize    = 5,
                        maxGSSize    = 500,
                        qvalueCutoff = 0.05,
                        TERM2GENE = mouse_database)
                        
        enrichment_results_table <- as_tibble(enrichment_results)
        write.csv(enrichment_results_table, paste0(local_path, 'OverRepresentation_analysis_results_', database,'.csv'))

        if (nrow(enrichment_results_table) > 0) {
            ## Add similarity matrix to the termsim slot of enrichment result
            enrichment_results <- enrichplot::pairwise_termsim(enrichment_results)
            
            enrichment_results_table <- as_tibble(enrichment_results) 
            write.csv(enrichment_results_table, paste0(local_path, 'OverRepresentation_analysis_results_', database,'.csv'))

            dotplot(enrichment_results,
                showCategory=50,
                title = paste0('Overrepresentation analysis ', database),
                label_format = 60)
            ggsave(paste0('ORA_dotplot_', database,'.pdf'), width = 10, height = 18, path = local_path)

            ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
            emapplot(enrichment_results, showCategory = 60) + ggtitle(paste0('Overrepresentation analysis ', database))
            ggsave(paste0('ORA_network_', database,'.pdf'), width = 14, height = 18, path = local_path)
        }

        #### GSEA ####

        fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
        names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(entrezid)

        gsea_results <- GSEA(geneList     = fold_changes,
                    # OrgDb        = org.Mm.eg.db,
                    # ont          = "ALL",
                    # keyType = "SYMBOL",
                    minGSSize    = 5,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE,
                    TERM2GENE = mouse_database)
        
        gsea_results_table <- as_tibble(gsea_results)
        write.csv(gsea_results_table, paste0(local_path, 'GSEA_analysis_results_', database,'.csv'))

        if (nrow(gsea_results_table) > 1) {
            gsea_results_table <- as_tibble(gsea_results) 
            write.csv(gsea_results_table, paste0(local_path, 'GSEA_analysis_results_', database,'.csv'))

            dotplot(gsea_results,
                showCategory=50,
                title = paste0('GSEA analysis ', database),
                label_format = 60)
            ggsave(paste0('GSEA_dotplot_', database,'.pdf'), width = 10, height = 18, path = local_path)

            ## Add similarity matrix to thenes,  termsim slot of enrichment result
            gsea_results <- enrichplot::pairwise_termsim(gsea_results)

            ## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
            emapplot(gsea_results, showCategory = 60) + ggtitle(paste0('GSEA analysis ', database))
            ggsave(paste0('GSEA_network_', database,'.pdf'), width = 14, height = 18, path = local_path)
        }
        
    }

    
    return()
}
pathways_of_interest_analysis <- function (results,  pathways_of_interest, cluster,path='./', FC_threshold = 0.3, p_value_threshold = 0.1, max_overlaps = 1000, label_size = 5, group1, group2, comparison) {

    results <- results[which(duplicated(results$genes) == F),]
    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    results <- results[which(duplicated(results$entrezid) == F),] |>
    drop_na(entrezid)
        
    local_path <- paste0(path, 'pathways_of_interest_', cluster, '/')
    unlink(local_path, recursive = T)
    dir.create(local_path)


    pathways_of_interest <- pathways_of_interest |>
        pivot_longer(cols = everything(),values_to = 'genes', names_to = 'pathway', values_drop_na = TRUE) |>
        arrange(desc(pathway)) 

    pathways_of_interest$genes <- pathways_of_interest |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) 

    pathways_of_interest <- pathways_of_interest |> filter(!is.na(genes) )

    #head(pathways_of_interest)

    fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
    names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(entrezid)

    for (term_of_interest in unique(pathways_of_interest$pathway)) {

        term <- pathways_of_interest |> filter(pathway == term_of_interest)

        gsea_results <- GSEA(geneList     = fold_changes,
                    # OrgDb        = org.Mm.eg.db,
                    # ont          = "ALL",
                    # keyType = "SYMBOL",
                    minGSSize    = 4,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE,
                    TERM2GENE = term)
                    
        anno <- gsea_results[term_of_interest, c("NES", "pvalue", "p.adjust")]
        lab <- paste0(names(anno), "=",  round(anno, 4), collapse="\n")

        p1 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 1, base_size = 13,title = term_of_interest)

        x_position <- ggplot_build(p1)$layout$panel_params[[1]]$x.range[2]*0.75
        y_position <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-(ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-ggplot_build(p1)$layout$panel_params[[1]]$y.range[1])*0.17

        p1 <- p1 + annotate("text", x_position, y_position, label = lab, hjust=0, vjust=0, size = 5)
        p2 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 2, base_size = 13)
        p3 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 3, base_size = 13)     
                
        cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.5, 0.5, 1), align = 'v')    
        ggsave(paste0('GSEA_', term_of_interest, '.pdf'), path = local_path, height = 10, width = 8)

        ########## scatterplot ##########

        ## prepare for visualization
        results <- results %>% 
                        
                        mutate(                        
                            genes_to_label2 = ifelse(entrezid %in% term$genes,  genes ,''),
                            add_label = ifelse(entrezid %in% term$genes,  'YES' ,'NO')
                            ) |>
                        mutate(add_label <- factor(add_label, levels = c('NO', 'YES'))) |>     
                        mutate(fold_change_direction = case_when(
                                                log2FoldChange>= 0 ~ 'UP',
                                                log2FoldChange<= 0 ~ "DOWN",
                                                TRUE ~ 'NO')) |>
                        arrange(add_label)
                         
                            
        my_colors <- c( "#32228b", "gray")
        names(my_colors) <- c("YES", "NO")
        
        # Scatterplot
        limx <- results |> pull(paste0('Avg_', group1)) |> max()
        limy <- results |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*5
        
        results |> 
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), label = genes_to_label2, col = add_label))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.2,
                    show.legend = FALSE,
                    max.overlaps = max_overlaps,
                    max.time = 30,
                    max.iter = 10000000,
                    nudge_x = ifelse(results$fold_change_direction == 'UP', -0.75, 0.75),
                    nudge_y = ifelse(results$fold_change_direction == 'UP', 0.75, -0.75),
                    aes(segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            
            
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=term_of_interest,
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+       
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(local_path, 'Pseudobulk scatter ', term_of_interest, '.pdf'))
    }
    return()
}


pathways_of_interest_analysis2 <- function (results,  pathways_of_interest_table, cluster,path='./', FC_threshold = 0.3, p_value_threshold = 0.1, max_overlaps = 1000, label_size = 5, group1, group2, comparison) {

    results <- results[which(duplicated(results$genes) == F),]
    results$entrezid <-  results |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID)
    results <- results[which(duplicated(results$entrezid) == F),] |>
    drop_na(entrezid)
        
    local_path <- paste0(path, 'pathways_of_interest_', cluster, '/')
    unlink(local_path, recursive = T)
    dir.create(local_path)


    # pathways_of_interest_table <- pathways_of_interest_table |>
    #     pivot_longer(cols = everything(),values_to = 'genes', names_to = 'pathway', values_drop_na = TRUE) |>
    #     arrange(desc(pathway)) 

    # pathways_of_interest_table$genes <- pathways_of_interest_table |> pull(genes) |> bitr(fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db', drop = FALSE) |> pull(ENTREZID) 

    pathways_of_interest_table <- pathways_of_interest_table |> filter(!is.na(genes) )

    #head(pathways_of_interest_table)

    fold_changes <- results |> arrange(desc(log2FoldChange)) |> pull(log2FoldChange)
    names(fold_changes) <- results |> arrange(desc(log2FoldChange)) |> pull(entrezid)

    for (term_of_interest in unique(pathways_of_interest_table$gs_name)) {

        term <- pathways_of_interest_table |> filter(gs_name == term_of_interest)

        gsea_results <- GSEA(geneList     = fold_changes,
                    # OrgDb        = org.Mm.eg.db,
                    # ont          = "ALL",
                    # keyType = "SYMBOL",
                    minGSSize    = 4,
                    maxGSSize    = 500,
                    pvalueCutoff = 1,
                    verbose      = FALSE,
                    TERM2GENE = term)
                    
        anno <- gsea_results[term_of_interest, c("NES", "pvalue", "p.adjust")]
        lab <- paste0(names(anno), "=",  round(anno, 4), collapse="\n")

        p1 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 1, base_size = 13,title = term_of_interest)

        x_position <- ggplot_build(p1)$layout$panel_params[[1]]$x.range[2]*0.75
        y_position <- ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-(ggplot_build(p1)$layout$panel_params[[1]]$y.range[2]-ggplot_build(p1)$layout$panel_params[[1]]$y.range[1])*0.17

        p1 <- p1 + annotate("text", x_position, y_position, label = lab, hjust=0, vjust=0, size = 5)
        p2 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 2, base_size = 13)
        p3 <- enrichplot::gseaplot2(gsea_results, geneSetID = term_of_interest, pvalue_table = FALSE, subplots = 3, base_size = 13)     
                
        cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1.5, 0.5, 1), align = 'v')    
        ggsave(paste0('GSEA_', term_of_interest, '.pdf'), path = local_path, height = 10, width = 8)

        ########## scatterplot ##########

        ## prepare for visualization
        results <- results %>% 
                        
                        mutate(                        
                            genes_to_label2 = ifelse(entrezid %in% term$genes,  genes ,''),
                            add_label = ifelse(entrezid %in% term$genes,  'YES' ,'NO')
                            ) |>
                        mutate(add_label <- factor(add_label, levels = c('NO', 'YES'))) |>     
                        mutate(fold_change_direction = case_when(
                                                log2FoldChange>= 0 ~ 'UP',
                                                log2FoldChange<= 0 ~ "DOWN",
                                                TRUE ~ 'NO')) |>
                        arrange(add_label)
                         
                            
        my_colors <- c( "#32228b", "gray")
        names(my_colors) <- c("YES", "NO")
        
        # Scatterplot
        limx <- results |> pull(paste0('Avg_', group1)) |> max()
        limy <- results |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*5
        
        results |> 
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), label = genes_to_label2, col = add_label))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.2,
                    show.legend = FALSE,
                    max.overlaps = max_overlaps,
                    max.time = 30,
                    max.iter = 10000000,
                    nudge_x = ifelse(results$fold_change_direction == 'UP', -0.75, 0.75),
                    nudge_y = ifelse(results$fold_change_direction == 'UP', 0.75, -0.75),
                    aes(segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            
            
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=term_of_interest,
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+       
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(local_path, 'Pseudobulk scatter ', term_of_interest, '.pdf'))
    }
    return()
}


### Over representation analysis -- gProfiler2
gProfiler2_overrepresentation_analysis <- function (significant_genes_FC_ordered, local_path, comparison, filename = '') {

    significant_genes_FC_ordered <- list(comparison = significant_genes_FC_ordered)
    # Enrichment pathway analysis
    enrichment_results <- gost(query = significant_genes_FC_ordered, 
                    organism = "mmusculus", ordered_query = TRUE, 
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                    user_threshold = 0.05, correction_method = "g_SCS", 
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

gProfiler2_functional_analysis <- function (results,  cluster, comparison, path='./', FC_threshold) {

     results <- results[which(duplicated(results$genes) == F),]


####################################### UP ########################################

     local_path <- paste0(path, 'Functional_analysis_UP_gProfiler2_', cluster, '_', comparison, '/')
     unlink(local_path, recursive = T)
     dir.create(local_path)
     print(colnames(results))

     significant_genes <- results |> filter((padj < 0.05) & (log2FoldChange > FC_threshold)) |> arrange(padj) |> arrange(desc(log2FoldChange)) |> pull(genes)
     
     if (length(significant_genes) > 2) {
          gProfiler2_overrepresentation_analysis(significant_genes, local_path =  local_path , comparison = comparison, filename = '')

     }

######################################## DOWN ########################################

     local_path <- paste0(path, 'Functional_analysis_DOWN_gProfiler2_', cluster, '_', comparison, '/')
     unlink(local_path, recursive = T)
     dir.create(local_path)

     significant_genes <- results |> filter((padj < 0.05) & (log2FoldChange < -1*FC_threshold)) |> arrange(padj) |> arrange(log2FoldChange) |> pull(genes)

     if (length(significant_genes) > 2) {
          gProfiler2_overrepresentation_analysis(significant_genes, local_path =  local_path , comparison = comparison, filename = '')
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
### Pseudobulk function
########## Pseudobulk ##########

pseudobulk <- function (scRNAseq, comparison, group1, group2, cluster='all_clusters', path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold = 0.4, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4')) {

    my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")

    gene_lists_path <- paste0(path, 'gene_lists/')
    figures_path <- paste0(path, 'figures/')

    # unlink(gene_lists_path, recursive = T)
    # unlink(figures_path, recursive = T)
    dir.create(gene_lists_path)
    dir.create(figures_path)
    print(paste('Cluster',cluster))

    Idents(scRNAseq) <- comparison


    ########## Continue here ##########

    print('number of cells in group 1')
    #print(scRNAseq@meta.data |> filter( !!sym(comparison) 'Cluster Annotations'in'Cluster Annotations' c(group1) ) |> nrow())
    print(scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group1 )) |> nrow())
    print('number of cells in group 2')
    print(scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group2 )) |> nrow())
    


    #Checked there are enough cells
    if ((scRNAseq@meta.data |> 
            filter(str_detect( !!as.name(comparison) , group1 )) |> 
            nrow()  < 10) | 
            (scRNAseq@meta.data |> filter(str_detect( !!as.name(comparison) , group2 )) |> nrow()  < 10)) {
        DEG_count <- 'Not enough cells'
        DEG_UP_count <- 'Not enough cells'
        DEG_DOWN_count <- 'Not enough cells'
        return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
    }
    
    counts <- AggregateExpression(scRNAseq, group.by=c(comparison),
                            assays='RNA',
                            slot='counts',
                            return.seurat=FALSE)

    counts <- counts$RNA |> as.data.frame()

    # print(counts)
    # print(colnames( counts ))
    
    # Run DE Analysis
    #Generate sample level metadata
    colData <- data.frame(samples=colnames(counts)) |>
                mutate(condition = ifelse(grepl(group1, samples), group1, group2))
    
    ## Filter
    counts <- counts |> mutate(row_sums=rowSums(counts)) |> filter(row_sums >= 10) |> dplyr::select(-row_sums)
    
    print('Group 1 Length')
    print(nrow(colData |> filter(condition == group1)))
    print('Group 2 Length')
    print(nrow(colData |> filter(condition == group2)))

    #Perform DESeq2
    if ((length(unique(colData$condition)) != 2 ) | (nrow(colData |> filter(condition == group1)) < 2) | (nrow(colData |> filter(condition == group2)) < 2)) {
        DEG_count <- 'Not enough biological replicates per group'
        DEG_UP_count <- 'Not enough biological replicates per group'
        DEG_DOWN_count <- 'Not enough biological replicates per group'
        return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
    }
    
    #Create DESeq2 object
    
    dds <- DESeqDataSetFromMatrix(countData = counts,
                            colData = colData,
                            design = ~condition)

    dds$condition <- factor(dds$condition, levels = c(group1, group2))

    

    ## DESeq2 QC
    rld <- rlog(dds, blind=TRUE) #rlog normalization

    DESeq2::plotPCA(rld, ntop=500, intgroup='condition') #PCA
    ggsave(filename=paste0('Pseudobulk_PCA_', cluster, '.pdf'), path=figures_path) 

    PCA_table <- DESeq2::plotPCA(rld, ntop=500, intgroup='condition', returnData = T) #PCA table
    write.csv(PCA_table, file=paste(path, 'PCA_pseudobulk', cluster, group2, 'vs', group1, '.csv', sep='_'))

    #################### Run DESeq2
    dds <- DESeq(dds)

    #Check the coefficients for the scRNAseq
    resultsNames(dds)

    #Generate results object
    results <- results(dds) |> as.data.frame()
    
    #Get Normalized Counts
    normalized_counts <- counts(dds, normalized = T)
    normalized_counts <- normalized_counts |>
            as.data.frame() |>
            rownames_to_column('genes') |>
            as_tibble() |>
            rowwise() |>
            mutate(
                !!paste0('Avg_', group2) := mean(c_across(contains(group2))),
                !!paste0('Avg_', group1) := mean(c_across(contains(group1))),
            ) |>
            ungroup()

    #Add gene annotations:
    annotations <- read.csv("M:/LPD-MIS members' data/Eduard Ansaldo/Bioinformatics/analysis_templates/annotations.csv")
    results <- results |>
                    rownames_to_column('genes') |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('genes' = 'gene_name')) |>
                    left_join(y = normalized_counts, by = c('genes' = 'genes'))

    head(results)

    results_filtered <- filter(results, padj < p_value_threshold & ((!!sym(paste0('Avg_', group2)) > expression_threshold_for_gene_list) | !!sym(paste0('Avg_', group1)) > expression_threshold_for_gene_list) & (log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)) %>% arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange >  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <  FC_threshold)



    write.csv(results |> arrange(desc(padj)), file= paste(paste(gene_lists_path, 'ALL_GENES_DEG_Analysis', sep=''), 'pseudobulk', cluster, group2, 'vs', group1, '.csv', sep='_'))
    write.csv(results_filtered_UP, file=paste(paste(gene_lists_path, 'DEG_UP', sep=''), 'pseudobulk', cluster, group2, 'vs', group1, '.csv', sep='_'))
    write.csv(results_filtered_DOWN, file=paste(paste(gene_lists_path, 'DEG_DOWN', sep=''), 'pseudobulk', cluster, group2, 'vs', group1, '.csv', sep='_'))

        

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)


    ## prepare for visualization
    results_scatter <- results %>% 
                    drop_na(pvalue) |>
                    mutate(
                        log10_pval = log10(padj+10^-90)*-1,
                        distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2)) ,
                        #bottom_limit = (!!sym(paste0('Avg_', group1))-x_intercept^2)/!!sym(paste0('Avg_', group1)),
                        #top_limit = bottom_limitp
                         )
                        
    results_scatter <- results_scatter |> mutate(
                        genes_to_label_first = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold) & (distance_from_diagonal > distance_from_diagonal_threshold) & ((!!sym(paste0('Avg_', group2)) > 100) | !!sym(paste0('Avg_', group1)) > 100), genes,NA),
                        genes_to_label_second  = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold) & ((!!sym(paste0('Avg_', group2)) > label_threshold) | !!sym(paste0('Avg_', group1)) > label_threshold) & is.na(genes_to_label_first), genes,NA),
                        genes_to_label = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold) & is.na(genes_to_label_first) & is.na(genes_to_label_second), genes,NA),
                        genes_to_label_volcano = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold), genes,NA),
                        diffexpressed = case_when(
                            log2FoldChange>=FC_threshold & padj < p_value_threshold  ~ 'UP', 
                            log2FoldChange<=-1*FC_threshold & padj < p_value_threshold  ~ "DOWN",
                            TRUE ~ 'NO'))  |>
                    mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))) |>
                    arrange(diffexpressed)          
    
    print(head(results_scatter |> filter(diffexpressed != 'NO') |> arrange(desc(distance_from_diagonal)) |> pull(distance_from_diagonal), n = 11))

    #ggplot(results, aes(x= !!paste0('Avg_', group1), y = !!paste0('Avg_', group2), label = genes_to_label, col = diffexpressed))+

    
    # Scatterplot

    limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
    limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
    mylims <- max(limx, limy)*6
       
    results_scatter |> 
        ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
            geom_point(size=1.3
            ) +
            geom_abline(slope = 1, intercept = 0)+
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 30,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label_first,segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label_second, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste('Pseudobulk DEGs in', str_replace(cluster,pattern = '_',replace = ' ') ),
                    x=paste0('Average Normalized Counts in ', group1),
                    y=paste0('Average Normalized Counts in ',  group2))+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=11),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13))+
       scale_x_log10(limits =  c(0.5, mylims))+
       scale_y_log10(limits =  c(0.5, mylims))

    ggsave(paste0(figures_path, 'Pseudobulk scatter DEG in ', cluster, '.pdf'))


    #Filter values which are not significant but with high  FC that would bias the plot visualization
    initial_number_of_genes <- nrow(results_scatter)
    max_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
    min_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
    if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
        min_FC_up_significant <- -3
    }  
    if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
        max_FC_up_significant <- 3    
    }
    results_volcano <- results_scatter %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant)))
    final_number_of_genes <- nrow(results_volcano)
    print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

    # # Volcano Plot
    results_volcano |> 
        arrange(desc(padj)) |>
        ggplot(aes(x=log2FoldChange, y=log10_pval, label=genes_to_label_volcano, col=diffexpressed)) +
        geom_point(size=1.5) +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            aes(segment.size=0.5, segment.alpha=0.8, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        geom_vline(xintercept=FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste('Pseudobulk DEGs in', str_replace(cluster,pattern = '_',replace = ' ') ),
                    x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                    y= 'Log10 Adj. p-value')+
        # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=15),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13),
            )         +
            scale_y_continuous(n.breaks = 8) +
            scale_x_continuous(n.breaks = 8)
    ggsave(paste0(figures_path, 'Pseudobulk volcano DEG in ', cluster, '.pdf'))

########## Overrepresentation analysis ##########

gProfiler2_functional_analysis(results,  cluster = cluster, comparison = comparison, path= path , FC_threshold = FC_threshold)


if (!is.null(pathways_of_interest)) {
    pathways_of_interest_analysis(results = results, pathways_of_interest = pathways_of_interest,  cluster = cluster, path = path, group1 = group1, group2 = group2, comparison = comparison)
    }

########## Plotting individual genes of interest ##########
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]
            
            print(genes_to_plot)
            
        
        results_scatter <- results_scatter |> mutate(
                            genes_to_label_first = ifelse(genes %in% genes_to_plot, genes, NA),
                            diffexpressed = case_when(
                                log2FoldChange>=FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot) ~ 'UP', 
                                log2FoldChange<=-1*FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot)  ~ "DOWN",
                                genes %in% genes_to_plot & log2FoldChange >= 0 & padj < p_value_threshold ~ 'INTEREST UP',
                                genes %in% genes_to_plot & log2FoldChange < 0 & padj < p_value_threshold ~ 'INTEREST DOWN',
                                TRUE ~ 'NO'))  |>
                        mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP', 'INTEREST DOWN', 'INTEREST UP'))) |>
                        arrange(diffexpressed)          
        
        my_colors <- c("green4", "darkorchid4", "gray", 'green3', 'darkorchid3')
        names(my_colors) <- c("DOWN", "UP", "NO", 'INTEREST DOWN', 'INTEREST UP')

        print(head(results_scatter |> filter(diffexpressed != 'NO') |> arrange(desc(distance_from_diagonal)) |> pull(distance_from_diagonal), n = 11))

        #ggplot(results, aes(x= !!paste0('Avg_', group1), y = !!paste0('Avg_', group2), label = genes_to_label, col = diffexpressed))+

        
        # Scatterplot

        limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
        limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*6
        
        results_scatter |> 
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.35,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    max.time = 10,
                    max.iter = 10000000,
                    nudge_x = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', -1, 1),
                    nudge_y = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', 0.75, -0.75),
                    aes(label = genes_to_label_first,segment.size=0.5, segment.alpha=1, segment.curvature=0)) +
            
            scale_colour_manual(values=my_colors)+
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=paste('Pseudobulk DEGs in', str_replace(cluster,pattern = '_',replace = ' ') ),
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(figures_path, 'Pseudobulk scatter ',gene_list,' in ', cluster, ' ', group2, ' vs ', group1, '.pdf'))


        #Filter values which are not significant but with high  FC that would bias the plot visualization
        initial_number_of_genes <- nrow(results_scatter)
        max_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
        min_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
        if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
            min_FC_up_significant <- -3
        }  
        if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
            max_FC_up_significant <- 3    
        }
        results_volcano <- results_scatter %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant)))
        final_number_of_genes <- nrow(results_volcano)
        print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

        # # Volcano Plot
        results_volcano |> 
            #arrange(desc(padj)) |>
            ggplot(aes(x=log2FoldChange, y=log10_pval, label=genes_to_label_first, col=diffexpressed)) +
            geom_point(size=1.5) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = max_overlaps,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', 3 ,-3),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 3, 3),
                #nudge_x = (5),
                aes(segment.size=0.3, segment.alpha=1, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            geom_vline(xintercept= FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
            geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
            geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=paste('Pseudobulk DEGs in', str_replace(cluster,pattern = '_',replace = ' ') ),
                        x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                        y= 'Log10 Adj. p-value')+
            # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=15),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13),
                )         +
                scale_y_continuous(n.breaks = 8) +
                scale_x_continuous(n.breaks = 8)
        ggsave(paste0(figures_path, 'Pseudobulk volcano ',gene_list,' in ', cluster, ' ', group2, ' vs ', group1, '.pdf'))
 
    }}

     return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
   
}                            
## Wilcox DE analysis
########## DEG_FindMarkers ##########

DEG_FindMarkers <- function (scRNAseq, comparison, group1, group2, cluster='all_clusters', path='./', FC_threshold = 0.3, gene_lists_to_plot = NULL, p_value_threshold = 0.05) {

    Idents(scRNAseq) <- comparison 
    print(cluster)
    print(group1)
    print(group2)
    
    DEG_scRNAseq <- FindMarkers(object = scRNAseq, ident.1 = group1, ident.2 = group2)
                                        #min.pct=0.005,
                                        #logfc.threshold=FC_threshold)

    scRNAseq_CPM <- scRNAseq |> AggregateExpression(group.by=c(comparison),
                                        assays = 'RNA',
                                        slot = 'data',
                                        return.seurat=TRUE,
                                        normalization.method='RC',
                                        scale.factor = 1e6)

    counts_CPM <- scRNAseq_CPM |> GetAssayData(assay = 'RNA', slot = 'data') |>
                            as.data.frame() |>
                            rownames_to_column(var = 'gene') |>
                            mutate(
                                !!paste0('Avg_', group2) := !!as.name(group2),
                                !!paste0('Avg_', group1) := !!as.name(group1)
                            )     
    
    DEG_scRNAseq <- DEG_scRNAseq %>% mutate(avg_log2FC = avg_log2FC*-1) 

    #Add gene annotations:
    annotations <- read.csv("M:/LPD-MIS members' data/Eduard Ansaldo/Bioinformatics/analysis_templates/annotations.csv")    
    DEG_scRNAseq <- DEG_scRNAseq |>
                    rownames_to_column('genes') |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('genes' = 'gene_name')) |>
                    left_join(y = counts_CPM, by = c('genes' = 'gene'))

    DEG_scRNAseq_filtered <- filter(DEG_scRNAseq, p_val_adj < 0.05) %>% arrange(p_val_adj)
    DEG_scRNAseq_filtered_UP <- filter(DEG_scRNAseq_filtered, avg_log2FC >  0) 
    DEG_scRNAseq_filtered_DOWN <- filter(DEG_scRNAseq_filtered, avg_log2FC <  0)

    write.csv(DEG_scRNAseq_filtered, file= paste(paste(path, 'DEG', sep=''), 'FindMarkers', cluster, group1, 'vs', group2, '.csv', sep='_'))
    write.csv(DEG_scRNAseq_filtered_UP, file=paste(paste(path, 'DEG_UP', sep=''), 'FindMarkers', cluster, group1, 'vs', group2, '.csv', sep='_'))
    write.csv(DEG_scRNAseq_filtered_DOWN, file=paste(paste(path, 'DEG_DOWN', sep=''), 'FindMarkers', cluster, group1, 'vs', group2, '.csv', sep='_'))

    # Return number of DEGs:
    DEG_count <- nrow(DEG_scRNAseq_filtered)
    DEG_UP_count <- nrow(DEG_scRNAseq_filtered_UP)
    DEG_DOWN_count <- nrow(DEG_scRNAseq_filtered_DOWN)

    ## prepare for visualization
    DEG_scRNAseq <- DEG_scRNAseq %>% mutate(
                        log10_pval = log10(p_val_adj+10^-90)*-1,
                        genes_to_label = ifelse((avg_log2FC >= 0.2 | avg_log2FC  <= -0.2)  & (p_val_adj < 0.05), genes,NA),
                        diffexpressed = case_when(
                                                avg_log2FC>=0.2 & p_val_adj < 0.05  ~ 'UP',
                                                avg_log2FC<=-0.2 & p_val_adj < 0.05  ~ "DOWN",
                                                TRUE ~ 'NO')) |>
                        mutate(diffexpressed  =  factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))) |>
                        arrange(diffexpressed)

    my_colors <- c("green4", "darkorchid4", "gray")    
    names(my_colors) <- c("DOWN", "UP", "NO")

    # Scatterplot
    DEG_scRNAseq |> 
        # arrange(desc(p_val_adj)) |>
        ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), label = genes_to_label, col = diffexpressed))+
            geom_point(size=1.5
            ) +
            geom_abline(slope = 1, intercept = 0)+
            geom_text_repel(
            size=5.5,
            box.padding = 0.7,
            show.legend = FALSE,
            max.overlaps = 6,
            aes(segment.size=0.5, segment.alpha=0.5, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        #  geom_vline(xintercept=0.2, col="lavenderblush2", linetype=2, size=0.5) +
        #  geom_vline(xintercept=-0.2, col="lavenderblush2", linetype=2, size=0.5) +
        #  geom_hline(yintercept=-1*log10(0.05), col="lavenderblush2", linetype=2, size=0.5)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste('FindMarkers DEGs in', cluster ),
                    x=paste0('Average CPM in ', group1),
                    y=paste0('Average CPM in ',  group2),
        # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=20),
            axis.text= element_text(size=13),
            axis.title= element_text(size=15)))+
        scale_y_continuous(trans = 'log10')+
        scale_x_continuous(trans = 'log10')

    ggsave(paste0(path, 'FindMarkers scatter DEG in ', cluster, '.pdf'))


    #Filter values which are not significant but with high  FC that would bias the plot visualization
    initial_number_of_genes <- nrow(DEG_scRNAseq)
    max_FC_up_significant <- DEG_scRNAseq %>% filter(diffexpressed != 'NO') %>% dplyr::select(avg_log2FC) %>% max()
    min_FC_up_significant <- DEG_scRNAseq %>% filter(diffexpressed != 'NO') %>% dplyr::select(avg_log2FC) %>% min()
    DEG_scRNAseq <- DEG_scRNAseq %>% filter(!(diffexpressed == 'NO' & (avg_log2FC <= min_FC_up_significant | avg_log2FC >= max_FC_up_significant)))
    final_number_of_genes <- nrow(DEG_scRNAseq)
    
    print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))
    
    #Volcano Plot
    ggplot(DEG_scRNAseq, aes(x=avg_log2FC, y=log10_pval, label=genes_to_label, col=diffexpressed)) +
    geom_point(size=1.5) +
    geom_text_repel(
        size=5.3,
        box.padding = 1,
        show.legend = FALSE,
        max.overlaps = 10,
        aes(segment.size=0.4, segment.alpha=0.3, segment.curvature=0)) +
    scale_colour_manual(values=my_colors)+
    geom_vline(xintercept=0.2, col="lavenderblush2", linetype=2, size=0.5) +
    geom_vline(xintercept=-0.2, col="lavenderblush2", linetype=2, size=0.5) +
    geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
    theme(text=element_text(size=20), legend.position="none")+
    labs(title=paste('FindMarkers DEGs in', cluster ),
                x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                y= 'Log10 Adj. p-value')+
    # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
    theme_classic(base_size = 28, base_line_size=1) +
    theme(legend.position="none", 
          title = element_text(size=20),
          axis.text= element_text(size=13),
          axis.title= element_text(size=15),
        )
    ggsave(paste(path, 'FindMarkers Differential gene expression in', cluster, '.pdf'))


    
########## Plotting individual genes of interest ##########
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]
            
            print(genes_to_plot)
            
        
        DEG_scRNAseq <- DEG_scRNAseq |> mutate(
                            genes_to_label = ifelse(genes %in% genes_to_plot, genes, NA),
                            diffexpressed = case_when(
                                avg_log2FC>=FC_threshold & p_val_adj < p_value_threshold & !(genes %in% genes_to_plot) ~ 'UP', 
                                avg_log2FC<=-1*FC_threshold & p_val_adj < p_value_threshold & !(genes %in% genes_to_plot)  ~ "DOWN",
                                genes %in% genes_to_plot & avg_log2FC >= 0 & p_val_adj < p_value_threshold ~ 'INTEREST UP',
                                genes %in% genes_to_plot & avg_log2FC < 0 & p_val_adj < p_value_threshold ~ 'INTEREST DOWN',
                                TRUE ~ 'NO'))  |>
                        mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP', 'INTEREST DOWN', 'INTEREST UP'))) |>
                        arrange(diffexpressed)          
        
        my_colors <- c("green4", "darkorchid4", "gray", 'green3', 'darkorchid3')
        names(my_colors) <- c("DOWN", "UP", "NO", 'INTEREST DOWN', 'INTEREST UP')

        #print(head(DEG_scRNAseq |> filter(diffexpressed != 'NO') |> arrange(desc(distance_from_diagonal)) |> pull(distance_from_diagonal), n = 11))

        #ggplot(results, aes(x= !!paste0('Avg_', group1), y = !!paste0('Avg_', group2), label = genes_to_label, col = diffexpressed))+

        
        # Scatterplot

        limx <- DEG_scRNAseq |> pull(paste0('Avg_', group1)) |> max()
        limy <- DEG_scRNAseq |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*6
        
        DEG_scRNAseq |> 
             ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), label = genes_to_label, col = diffexpressed))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    #size=label_size,
                    box.padding = 0.35,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    max.time = 10,
                    max.iter = 10000000,
                    nudge_x = ifelse(DEG_scRNAseq$diffexpressed == 'UP' | DEG_scRNAseq$diffexpressed == 'INTEREST UP', -1, 1),
                    nudge_y = ifelse(DEG_scRNAseq$diffexpressed == 'UP' | DEG_scRNAseq$diffexpressed == 'INTEREST UP', 0.75, -0.75),
                    aes(label = genes_to_label,segment.size=0.5, segment.alpha=1, segment.curvature=0)) +
            
            scale_colour_manual(values=my_colors)+
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=paste('Pseudobulk DEGs in', str_replace(cluster,pattern = '_',replace = ' ') ),
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(path, 'Scatter ',gene_list,' in ', cluster, ' ', group2, ' vs ', group1, '.pdf'))


        #Filter values which are not significant but with high  FC that would bias the plot visualization
        initial_number_of_genes <- nrow(DEG_scRNAseq)
        max_FC_up_significant <- DEG_scRNAseq %>% filter(diffexpressed != 'NO') %>% dplyr::select(avg_log2FC) %>% max(na.rm = T)
        min_FC_up_significant <- DEG_scRNAseq %>% filter(diffexpressed != 'NO') %>% dplyr::select(avg_log2FC) %>% min(na.rm = T)
        if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
            min_FC_up_significant <- -3
        }  
        if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
            max_FC_up_significant <- 3    
        }
        results_volcano <- DEG_scRNAseq %>% filter(!(diffexpressed == 'NO' & (avg_log2FC < min_FC_up_significant | avg_log2FC > max_FC_up_significant)))
        final_number_of_genes <- nrow(results_volcano)
        print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

        # # Volcano Plot
        results_volcano |> 
            #arrange(desc(p_val_adj)) |>
            ggplot(aes(x=avg_log2FC, y=log10_pval, label=genes_to_label, col=diffexpressed)) +
            geom_point(size=1.5) +
            geom_text_repel(
                #size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = Inf,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(DEG_scRNAseq$diffexpressed == 'UP' | DEG_scRNAseq$diffexpressed == 'INTEREST UP', 1 ,-1),
                nudge_y = ifelse(DEG_scRNAseq$diffexpressed == 'UP', 1, 1),
                #nudge_x = (5),
                aes(segment.size=0.3, segment.alpha=1, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            geom_vline(xintercept= FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
            geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
            geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=paste('Pseudobulk DEGs in', str_replace(cluster,pattern = '_',replace = ' ') ),
                        x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                        y= 'Log10 Adj. p-value')+
            # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=15),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13),
                )         +
                scale_y_continuous(n.breaks = 8) +
                scale_x_continuous(n.breaks = 8)
        ggsave(paste0(path, 'Volcano ',gene_list,' in ', cluster, ' ', group2, ' vs ', group1, '.pdf'))
 
    }}
                           

    return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
}
# Bulk functions
bulk_analysis <- function (counts_table, comparison = 'Groups', group1, group2, cell_type, path='./', FC_threshold = 0.3, p_value_threshold = 0.05, max_overlaps = 15, label_size = 5, pathways_of_interest = NULL, label_threshold = 100000, distance_from_diagonal_threshold =0.4, gene_lists_to_plot = NULL, expression_threshold_for_gene_list = 20, colors = c('green4', 'darkorchid4')) {

   my_colors <- c(colors, "gray")
    names(my_colors) <- c("DOWN", "UP", "NO")

    gene_lists_path <- paste0(path, 'gene_lists/')
    figures_path <- paste0(path, 'figures/')

    # unlink(gene_lists_path, recursive = T)
    # unlink(figures_path, recursive = T)
    dir.create(gene_lists_path)
    dir.create(figures_path)
    print(paste(cell_type))

    counts <- tibble(counts_table) |> column_to_rownames('genes')

    # Run DE Analysis
    #Generate sample level metadata
    colData <- data.frame(samples=colnames(counts)) |>
                mutate(condition = ifelse(grepl(group1, samples), group1, group2))
    
    ## Filter
   counts <- counts |> mutate(row_sums=rowSums(counts)) |> filter(row_sums >= 10) |> dplyr::select(-row_sums)
    
    print('Group 1 Length')
    print(nrow(colData |> filter(condition == group1)))
    print('Group 2 Length')
    print(nrow(colData |> filter(condition == group2)))

    #Perform DESeq2
    if ((length(unique(colData$condition)) != 2 ) | (nrow(colData |> filter(condition == group1)) < 2) | (nrow(colData |> filter(condition == group2)) < 2)) {
        DEG_count <- 'Not enough biological replicates per group'
        DEG_UP_count <- 'Not enough biological replicates per group'
        DEG_DOWN_count <- 'Not enough biological replicates per group'
        return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
    }
    
    #Create DESeq2 object
    
    dds <- DESeqDataSetFromMatrix(countData = counts,
                            colData = colData,
                            design = ~condition)

    dds$condition <- factor(dds$condition, levels = c(group1, group2))

    

    ## DESeq2 QC
    rld <- rlog(dds, blind=TRUE) #rlog normalization

    DESeq2::plotPCA(rld, ntop=500, intgroup='condition') #PCA
    ggsave(filename=paste0('Pseudobulk_PCA_', cell_type, '.pdf'), path=figures_path) 

    #################### Run DESeq2
    dds <- DESeq(dds)

    #Check the coefficients for the scRNAseq
    resultsNames(dds)

    #Generate results object
    results <- results(dds) |> as.data.frame()
    
    #Get Normalized Counts
    normalized_counts <- counts(dds, normalized = T)
    normalized_counts <- normalized_counts |>
            as.data.frame() |>
            rownames_to_column('genes') |>
            as_tibble() |>
            rowwise() |>
            mutate(
                !!paste0('Avg_', group2) := mean(c_across(contains(group2))),
                !!paste0('Avg_', group1) := mean(c_across(contains(group1))),
            ) |>
            ungroup()

    #Add gene annotations:
    annotations <- read.csv("M:/LPD-MIS members' data/Eduard Ansaldo/Bioinformatics/analysis_templates/annotations.csv")
    results <- results |>
                    rownames_to_column('genes') |>
                    left_join(y= unique(annotations[,c('gene_name', 'description')]),
                        by = c('genes' = 'gene_name')) |>
                    left_join(y = normalized_counts, by = c('genes' = 'genes'))

    head(results)

    results_filtered <- filter(results, padj < p_value_threshold & ((!!sym(paste0('Avg_', group2)) > expression_threshold_for_gene_list) | !!sym(paste0('Avg_', group1)) > expression_threshold_for_gene_list) & (log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)) %>% arrange(padj)
    results_filtered_UP <- filter(results_filtered, log2FoldChange >  FC_threshold) 
    results_filtered_DOWN <- filter(results_filtered, log2FoldChange <  FC_threshold)



    write.csv(results |> arrange(desc(padj)), file= paste(paste(gene_lists_path, 'ALL_GENES_DEG_Analysis', sep=''), 'pseudobulk', cell_type, group2, 'vs', group1, '.csv', sep='_'))
    write.csv(results_filtered_UP, file=paste(paste(gene_lists_path, 'DEG_UP', sep=''), 'pseudobulk', cell_type, group2, 'vs', group1, '.csv', sep='_'))
    write.csv(results_filtered_DOWN, file=paste(paste(gene_lists_path, 'DEG_DOWN', sep=''), 'pseudobulk', cell_type, group2, 'vs', group1, '.csv', sep='_'))

        

    # Return number of DEGs:
    DEG_count <- nrow(results_filtered)
    DEG_UP_count <- nrow(results_filtered_UP)
    DEG_DOWN_count <- nrow(results_filtered_DOWN)


    ## prepare for visualization
    results_scatter <- results %>% 
                    drop_na(pvalue) |>
                    mutate(
                        log10_pval = log10(padj)*-1,
                        distance_from_diagonal =  (abs((log10(!!sym(paste0('Avg_', group2))+1)) - (log10(!!sym(paste0('Avg_', group1))+1)))/sqrt(2)) ,
                        #bottom_limit = (!!sym(paste0('Avg_', group1))-x_intercept^2)/!!sym(paste0('Avg_', group1)),
                        #top_limit = bottom_limitp
                         )
                        
    results_scatter <- results_scatter |> mutate(
                        genes_to_label_first = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold) & (distance_from_diagonal > distance_from_diagonal_threshold) & ((!!sym(paste0('Avg_', group2)) > 100) | !!sym(paste0('Avg_', group1)) > 100), genes,NA),
                        genes_to_label_second  = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold) & ((!!sym(paste0('Avg_', group2)) > label_threshold) | !!sym(paste0('Avg_', group1)) > label_threshold) & is.na(genes_to_label_first), genes,NA),
                        genes_to_label = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold) & is.na(genes_to_label_first) & is.na(genes_to_label_second), genes,NA),
                        genes_to_label_volcano = ifelse((log2FoldChange >= FC_threshold | log2FoldChange  <= -1*FC_threshold)  & (padj < p_value_threshold), genes,NA),
                        diffexpressed = case_when(
                            log2FoldChange>=FC_threshold & padj < p_value_threshold  ~ 'UP', 
                            log2FoldChange<=-1*FC_threshold & padj < p_value_threshold  ~ "DOWN",
                            TRUE ~ 'NO'))  |>
                    mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP'))) |>
                    arrange(diffexpressed)          
    

    print(head(results_scatter |> filter(diffexpressed != 'NO') |> arrange(desc(distance_from_diagonal)) |> pull(distance_from_diagonal), n = 11))

    #ggplot(results, aes(x= !!paste0('Avg_', group1), y = !!paste0('Avg_', group2), label = genes_to_label, col = diffexpressed))+

    
    # Scatterplot

    limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
    limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
    mylims <- max(limx, limy)*6
       
    results_scatter |> 
        ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
            geom_point(size=1.3
            ) +
            geom_abline(slope = 1, intercept = 0)+
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 30,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label_first,segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label_second, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = 10,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP', -1, 1),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 0.75, -0.75),
                aes(label = genes_to_label, segment.size=0.3, segment.alpha=0.4, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste('Pseudobulk DEGs in', str_replace(cell_type,pattern = '_',replace = ' ') ),
                    x=paste0('Average Normalized Counts in ', group1),
                    y=paste0('Average Normalized Counts in ',  group2))+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=11),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13))+
       scale_x_log10(limits =  c(0.5, mylims))+
       scale_y_log10(limits =  c(0.5, mylims))

    ggsave(paste0(figures_path, 'Pseudobulk scatter DEG in ', cell_type, '.pdf'))


    #Filter values which are not significant but with high  FC that would bias the plot visualization
    initial_number_of_genes <- nrow(results_scatter)
    max_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
    min_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
    if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
        min_FC_up_significant <- -3
    }  
    if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
        max_FC_up_significant <- 3    
    }
    results_volcano <- results_scatter %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant)))
    final_number_of_genes <- nrow(results_volcano)
    print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

    # # Volcano Plot
    results_volcano |> 
        arrange(desc(padj)) |>
        ggplot(aes(x=log2FoldChange, y=log10_pval, label=genes_to_label_volcano, col=diffexpressed)) +
        geom_point(size=1.5) +
        geom_text_repel(
            size=label_size,
            box.padding = 0.35,
            show.legend = FALSE,
            max.overlaps = max_overlaps,
            max.time = 10,
            max.iter = 10000000,
            aes(segment.size=0.5, segment.alpha=0.8, segment.curvature=0)) +
        scale_colour_manual(values=my_colors)+
        geom_vline(xintercept=FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
        geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
        theme(text=element_text(size=20), legend.position="none")+
        labs(title=paste('Pseudobulk DEGs in', str_replace(cell_type,pattern = '_',replace = ' ') ),
                    x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                    y= 'Log10 Adj. p-value')+
        # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
        theme_classic(base_size = 28, base_line_size=1) +
        theme(legend.position="none", 
            title = element_text(size=15),
            axis.text= element_text(size=10),
            axis.title= element_text(size=13),
            )         +
            scale_y_continuous(n.breaks = 8) +
            scale_x_continuous(n.breaks = 8)
    ggsave(paste0(figures_path, 'Pseudobulk volcano DEG in ', cell_type, '.pdf'))

########## Overrepresentation analysis ##########

gProfiler2_functional_analysis(results, cluster =cell_type, comparison = comparison, path= path , FC_threshold = FC_threshold)


if (!is.null(pathways_of_interest)) {
    pathways_of_interest_analysis2(results = results, pathways_of_interest = pathways_of_interest,  cluster = cell_type, path = path, group1 = group1, group2 = group2, comparison = comparison)
    }

########## Plotting individual genes of interest ##########
    if (!is.null(gene_lists_to_plot)) {
        for (gene_list in names(gene_lists_to_plot)) {
            genes_to_plot <- gene_lists_to_plot[[gene_list]]
            
            print(genes_to_plot)
            
        
        results_scatter <- results_scatter |> mutate(
                            genes_to_label_first = ifelse(genes %in% genes_to_plot, genes, NA),
                            diffexpressed = case_when(
                                log2FoldChange>=FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot) ~ 'UP', 
                                log2FoldChange<=-1*FC_threshold & padj < p_value_threshold & !(genes %in% genes_to_plot)  ~ "DOWN",
                                genes %in% genes_to_plot & log2FoldChange >= 0 & padj < p_value_threshold ~ 'INTEREST UP',
                                genes %in% genes_to_plot & log2FoldChange < 0 & padj < p_value_threshold ~ 'INTEREST DOWN',
                                TRUE ~ 'NO'))  |>
                        mutate(diffexpressed = factor(diffexpressed, levels = c('NO', 'DOWN', 'UP', 'INTEREST DOWN', 'INTEREST UP'))) |>
                        arrange(diffexpressed)          
        
        my_colors <- c("green4", "darkorchid4", "gray", 'green3', 'darkorchid3')
        names(my_colors) <- c("DOWN", "UP", "NO", 'INTEREST DOWN', 'INTEREST UP')

        print(head(results_scatter |> filter(diffexpressed != 'NO') |> arrange(desc(distance_from_diagonal)) |> pull(distance_from_diagonal), n = 11))

        #ggplot(results, aes(x= !!paste0('Avg_', group1), y = !!paste0('Avg_', group2), label = genes_to_label, col = diffexpressed))+

        
        # Scatterplot

        limx <- results_scatter |> pull(paste0('Avg_', group1)) |> max()
        limy <- results_scatter |> pull(paste0('Avg_', group2)) |> max()
        mylims <- max(limx, limy)*6
        
        results_scatter |> 
            ggplot(aes(x = !!sym(paste0('Avg_', group1)), y = !!sym(paste0('Avg_', group2)), col = diffexpressed))+
                geom_point(size=1.3
                ) +
                geom_abline(slope = 1, intercept = 0)+
                geom_text_repel(
                    size=label_size,
                    box.padding = 0.35,
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    max.time = 10,
                    max.iter = 10000000,
                    nudge_x = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', -1, 1),
                    nudge_y = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', 0.75, -0.75),
                    aes(label = genes_to_label_first,segment.size=0.5, segment.alpha=1, segment.curvature=0)) +
            
            scale_colour_manual(values=my_colors)+
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=paste('Pseudobulk DEGs in', str_replace(cell_type,pattern = '_',replace = ' ') ),
                        x=paste0('Average Normalized Counts in ', group1),
                        y=paste0('Average Normalized Counts in ',  group2))+
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=11),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13))+
        scale_x_log10(limits =  c(0.5, mylims))+
        scale_y_log10(limits =  c(0.5, mylims))

        ggsave(paste0(figures_path, 'Pseudobulk scatter ',gene_list,' in ', cell_type, ' ', group2, ' vs ', group1, '.pdf'))


        #Filter values which are not significant but with high  FC that would bias the plot visualization
        initial_number_of_genes <- nrow(results_scatter)
        max_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% max(na.rm = T)
        min_FC_up_significant <- results_scatter %>% filter(diffexpressed != 'NO') %>% dplyr::select(log2FoldChange) %>% min(na.rm = T)
        if (min_FC_up_significant > -3 | is.na(min_FC_up_significant)) {
            min_FC_up_significant <- -3
        }  
        if (max_FC_up_significant  < 3 | is.na(max_FC_up_significant)) {
            max_FC_up_significant <- 3    
        }
        results_volcano <- results_scatter %>% filter(!(diffexpressed == 'NO' & (log2FoldChange < min_FC_up_significant | log2FoldChange > max_FC_up_significant)))
        final_number_of_genes <- nrow(results_volcano)
        print(paste('Removed', initial_number_of_genes-final_number_of_genes, 'non-significant genes that would bias the plot visualization'))

        # # Volcano Plot
        results_volcano |> 
            #arrange(desc(padj)) |>
            ggplot(aes(x=log2FoldChange, y=log10_pval, label=genes_to_label_first, col=diffexpressed)) +
            geom_point(size=1.5) +
            geom_text_repel(
                size=label_size,
                box.padding = 0.35,
                show.legend = FALSE,
                max.overlaps = max_overlaps,
                max.time = 10,
                max.iter = 10000000,
                nudge_x = ifelse(results_scatter$diffexpressed == 'UP' | results_scatter$diffexpressed == 'INTEREST UP', 3 ,-3),
                nudge_y = ifelse(results_scatter$diffexpressed == 'UP', 3, 3),
                #nudge_x = (5),
                aes(segment.size=0.3, segment.alpha=1, segment.curvature=0)) +
            scale_colour_manual(values=my_colors)+
            geom_vline(xintercept= FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
            geom_vline(xintercept=-FC_threshold, col="lavenderblush2", linetype=2, size=0.5) +
            geom_hline(yintercept=-1*log10(p_value_threshold), col="lavenderblush2", linetype=2, size=0.5)+
            theme(text=element_text(size=20), legend.position="none")+
            labs(title=paste('Pseudobulk DEGs in', str_replace(cell_type,pattern = '_',replace = ' ') ),
                        x=paste('Average log2 FC (', group2, '/', group1, ')', sep=''),
                        y= 'Log10 Adj. p-value')+
            # coord_cartesian(xlim=c(-1.15, 1.15), ylim = c(0, 95), expand = FALSE)+
            theme_classic(base_size = 28, base_line_size=1) +
            theme(legend.position="none", 
                title = element_text(size=15),
                axis.text= element_text(size=10),
                axis.title= element_text(size=13),
                )         +
                scale_y_continuous(n.breaks = 8) +
                scale_x_continuous(n.breaks = 8)
        ggsave(paste0(figures_path, 'Pseudobulk volcano ',gene_list,' in ', cell_type, ' ', group2, ' vs ', group1, '.pdf'))
 
    }}

     return(c(all_count=DEG_count, UP_count=DEG_UP_count, DOWN_count=DEG_DOWN_count))
   
}                            