#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(kableExtra)
library(stargazer)
library(readr)
library(deconica)
library(tibble)
library(FactoMineR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(fgsea)
library(knitr)
library(tidyr)
library(tidyverse)
library(data.table)

library(factoextra)
library(cluster)
library(NbClust)

library(shiny)
library(plotly)
library(stringr)

library(kableExtra)
library(GGally)

# Define server logic required to draw a histogram
function(input, output, session) {
########################################################FONCTION /OBJET############################################################
  
  
########################################################ACCUEIL####################################################################
  CancerSEA_sig <- read_table("CancerSEA.txt")
  
  # Render the table using renderUI
  output$cancerSEA_sig_table <- renderText({
    
    CancerSEA_sig %>%
      mutate(across(everything(), ~replace_na(.x, "")))%>%
      kable(format = "html", 
            escape = FALSE) %>%
      kable_styling(bootstrap_options = c("condensed"),
                    full_width = TRUE,
                    position = "center") %>%
      scroll_box(width = "100%", height = "500px") %>%
      row_spec(0,color = "orchid", bold = TRUE) %>%
      as.character()
 
  })
  
  Hallmark_sig <- read_table("Cancer_Gene_Census_Hallmarks_Of_Cancer.txt")
  
  output$Hallmarks_sig_table <- renderText({
    Hallmark_sig %>%
      mutate(across(everything(), ~replace_na(.x, "")))%>%
      kable(format = "html", 
            escape = FALSE) %>%
      kable_styling(bootstrap_options = c( "condensed"),
                    full_width = TRUE,
                    position = "center") %>%
      scroll_box(width = "100%", height = "500px") %>%
      row_spec(0,color = "orchid", bold = TRUE) %>%
      as.character()
  })
  
  XBP1s_sig <- c( "ASS1", "C3", "CCL20", "COL4A6", "CXCL2", "CXCL5", "CXCL8", "IFI44L", "IL1B", "IL6", "KCNN2", "MMP1", "MMP12", "MMP3", "PLA2G4A", "PPP4R4", "SERPINB2", "TFPI2", "ZNF804A")
  RIDD_sig <- c("ANGPT1", "CFH", "CFI", "CLEC3B", "COL3A1", "COL8A1", "DACH1", "DCN", "FHL1", "GAS1", "LUM", "OXTR", "PLAC8", "RGS4", "TAGLN", "TGFB2", "THBS1", "TIMP3", "TMEM255A")
  df_38 <- data.frame("XBP1s" = XBP1s_sig, "RIDD" = RIDD_sig)
  
  output$XBP1s_RIDD_19_sig_table <- renderText({
    df_38 %>%
      mutate(across(everything(), ~replace_na(.x, "")))%>%
      kable(format = "html", 
            escape = FALSE) %>%
      kable_styling(bootstrap_options = c("condensed"),
                    full_width = TRUE,
                    position = "center") %>%
      scroll_box(width = "100%", height = "500px") %>%
      row_spec(0,color = "orchid", bold = TRUE) %>%
      as.character()
  })
  
  IRE1_PUR <-read_table("annotation_allprobes.txt")
  
  XBP1_PUR <- filter(IRE1_PUR,Component == "XBP1s")
  XBP1_PUR<- XBP1_PUR$Affy_geneIDs
  
  RIDD_PUR <- filter(IRE1_PUR,Component == "RIDD")
  RIDD_PUR<- RIDD_PUR$Affy_geneIDs
  RIDD_PUR <- c(RIDD_PUR, rep(" ", 3))
  
  df_pur <- data.frame("XBP1s" = XBP1_PUR, "RIDD" = RIDD_PUR)
  
  output$XBP1s_RIDD_pur_sig_table <- renderText({
    df_pur %>%
      mutate(across(everything(), ~replace_na(.x, "")))%>%
      kable(format = "html", 
            escape = FALSE) %>%
      kable_styling(bootstrap_options = c("condensed"),
                    full_width = TRUE,
                    position = "center") %>%
      scroll_box(width = "100%", height = "500px") %>%
      row_spec(0,color = "orchid", bold = TRUE) %>%
      as.character()
  })
########################################################DECONVOLUTION####################################################################    
  
  TCGA_counts <- read_tsv("TCGA-GBM.htseq_counts.tsv")
  TCGA_counts <- column_to_rownames(TCGA_counts, var = colnames(TCGA_counts)[1])
  
  
  output$resume_TCGA  <- renderDT({
    datatable(TCGA_counts, options = list(pageLength = 5))
      
  })  
  
  code_to_show <- '
  TCGA_ica<- run_fastica (
  TCGA_final,
  overdecompose = FALSE,
  with.names = FALSE,
  gene.names = row.names(TCGA_final),
  samples = colnames(TCGA_final),
  n.comp = 25,
  R = TRUE
)
  '
  
  # Use renderPrint to display the code
  output$CodeCleanChunk <- renderPrint({
    cat(code_to_show)
  })
  
  # Nettoyer les Ensembl IDs dans TCGA_counts
  TCGA_counts <- TCGA_counts %>%
    rownames_to_column("ensembl_id")
  
  TCGA_counts$ensembl_id <- sub("\\..*", "", TCGA_counts$ensembl_id)
  
  # Lecture et nettoyage
  correspondance <- read.table("ENSG_ID.txt", header = TRUE, sep = ",")
  
  # Nettoyer les IDs s'ils contiennent une version
  correspondance$Gene.stable.ID <- sub("\\..*", "", correspondance$Gene.stable.ID)
  
  # Garder seulement les colonnes utiles
  correspondance_clean <- correspondance %>%
    dplyr::select(Gene.stable.ID, Gene.name) %>%
    dplyr::filter(Gene.name != "")  # On enlève les lignes sans nom de gène
  
  # Jointure
  TCGA_joined <- TCGA_counts %>%
    left_join(correspondance_clean, by = c("ensembl_id" = "Gene.stable.ID"))
  
  # Supprimer les lignes sans nom de gène
  TCGA_filtered <- TCGA_joined %>%
    filter(!is.na(Gene.name))
  
  # Supprimer les doublons si même nom de gène apparaît plusieurs fois
  TCGA_unique <- TCGA_filtered %>%
    group_by(Gene.name) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  # Finalisation
  TCGA_final <- TCGA_unique %>%
    column_to_rownames("Gene.name")
  
  TCGA_ica<- run_fastica (
    TCGA_final,
    overdecompose = FALSE,
    with.names = FALSE,
    gene.names = row.names(TCGA_final),
    samples = colnames(TCGA_final),
    n.comp = 25,
    R = TRUE
  )
########################################################TAB XBP1 ET RIDD####################################################################  
  
  basis.list_XBP1_RIDD_pur <- c(
    list(XBP1_pur = data.frame(gene = XBP1_PUR),
         RIDD_pur = data.frame(gene = RIDD_PUR))
  )
  
  basis.list_XBP1_RIDD <- c(
    list(XBP1 = data.frame(gene = XBP1s_sig),
         RIDD = data.frame(gene = RIDD_sig))
  )
  
  ex_gsea <- data.frame(
    Description = c("Up regulated genes with a low pvalue",
                    "Up regulated genes with a high pvalue",
                    "Down regulated genes with a high pvalue",
                    "Down regulated genes with a low pvalue")
  )
  
  output$exemple_GSEA_XR <- renderUI({
    HTML(
      ex_gsea %>%
        kable(format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "orchid", bold = TRUE) %>%
        as.character()
    )
  })
  
  # On crée une liste pour mettre tout nos résultats
  all_fgsea_XBP1_RIDD <- list()
  
  # S'assurer que les noms des stats sont corrects
  TCGA_ica$names <- rownames(TCGA_final)
  
  for (i in 1:ncol(TCGA_ica$S)) {
    component_name <- paste0("ICA_", i)
    
    stats <- TCGA_ica$S[, i]
    names(stats) <- TCGA_ica$names
    stats <- stats[!is.na(stats)]  # filtrer les NA
    
    for (pathway_name in names(basis.list_XBP1_RIDD)) {
      # Extraire les gènes de la signature
      genes <- basis.list_XBP1_RIDD[[pathway_name]]$gene
      
      # Créer un sous-ensemble de pathways comme liste
      pathway_list <- list()
      pathway_list[[pathway_name]] <- genes
      
      # Vérifier qu'il y a assez de gènes
      if (length(intersect(names(stats), genes)) >= 5) {
        fgsea_res_XBP1_RIDD <- fgsea(
          pathways = pathway_list,
          stats = stats,
          minSize = 5,
          maxSize = 20000
        )
        
        # Ajouter le nom de la composante et du pathway
        fgsea_res_XBP1_RIDD$Component <- component_name
        fgsea_res_XBP1_RIDD$Pathway <- pathway_name
        
        all_fgsea_XBP1_RIDD[[length(all_fgsea_XBP1_RIDD) + 1]] <- fgsea_res_XBP1_RIDD
      }
    }
  }
  
  
  # Combiner en un seul data.frame
  fgsea_all_results_XBP1_RIDD <- bind_rows(all_fgsea_XBP1_RIDD)
  
  
  ## Pour XBP1s pur et RIDD pur
  all_fgsea_XBP1_RIDD_pur <- list()
  # S'assurer que les noms des stats sont corrects
  TCGA_ica$names <- rownames(TCGA_final)
  for (i in 1:ncol(TCGA_ica$S)) {
    component_name <- paste0("ICA_", i)
    
    stats <- TCGA_ica$S[, i]
    names(stats) <- TCGA_ica$names
    stats <- stats[!is.na(stats)]  # filtrer les NA
    
    for (pathway_name in names(basis.list_XBP1_RIDD_pur)) {
      # Extraire les gènes de la signature
      genes <- basis.list_XBP1_RIDD_pur[[pathway_name]]$gene
      
      # Créer un sous-ensemble de pathways comme liste
      pathway_list <- list()
      pathway_list[[pathway_name]] <- genes
      
      # Vérifier qu'il y a assez de gènes
      if (length(intersect(names(stats), genes)) >= 5) {
        fgsea_res_XBP1_RIDD_pur <- fgsea(
          pathways = pathway_list,
          stats = stats,
          minSize = 5,
          maxSize = 20000
        )
        
        # Ajouter le nom de la composante et du pathway
        fgsea_res_XBP1_RIDD_pur$Component <- component_name
        fgsea_res_XBP1_RIDD_pur$Pathway <- pathway_name
        
        all_fgsea_XBP1_RIDD_pur[[length(all_fgsea_XBP1_RIDD_pur) + 1]] <- fgsea_res_XBP1_RIDD_pur
      }
    }
  }
  # Combiner en un seul data.frame
  fgsea_all_results_XBP1_RIDD_pur <- bind_rows(all_fgsea_XBP1_RIDD_pur)
  
  
  ## On combine le pur et ceux d'IRE1 38
  fgsea_combined <- bind_rows(fgsea_all_results_XBP1_RIDD_pur, fgsea_all_results_XBP1_RIDD)
  
  
  
  p_combined <- ggplot(fgsea_combined, aes(x = Component, y = Pathway)) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size(range = c(1, 6)) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(color = "NES", size = "-log10(padj)")
  
  output$result_GSEA_IRE1_38 <- renderUI({
    HTML(fgsea_all_results_XBP1_RIDD %>%
      kable(format = "html", 
            escape = FALSE) %>%
      kable_styling(bootstrap_options = c("condensed"),
                    full_width = TRUE,
                    position = "center") %>%
      scroll_box(width = "100%", height = "500px") %>%
      row_spec(0,color = "orchid", bold = TRUE) %>%
      as.character())
  })
  
  output$result_GSEA_IRE1_pur <- renderUI({
    HTML(fgsea_all_results_XBP1_RIDD_pur %>%
           kable(format = "html", 
                 escape = FALSE) %>%
           kable_styling(bootstrap_options = c("condensed"),
                         full_width = TRUE,
                         position = "center") %>%
           scroll_box(width = "100%", height = "500px") %>%
           row_spec(0,color = "orchid", bold = TRUE) %>%
           as.character())
  })

  output$result_GSEA_IRE1_graph <- renderPlotly({
    
    
    ggplotly(p_combined)
  })

  
  #Ici on fait le score d'abondance des 10 gènes les plus exprimés dans chaque composante
  TCGA_final_marker <- generate_markers(TCGA_ica,10)
  score <- get_scores(TCGA_ica$log.counts,TCGA_final_marker)
  
  # On prend les résultats les plus significatif (pval<1%)
  fgsea_all_results_XBP1_RIDD_significatif <- filter(fgsea_all_results_XBP1_RIDD,padj<0.05)
  
  #on associe à chaque composante la pathway adapté (on la associé grace au test de GSEA )
  # Alors comme des signatures peuvent être associé à plusieurs composantes on va nommer les signatures de la facon suivante (sig_UP_1,sig_UP_2,...)
  correspondance_Comp_XBP1_RIDD <- fgsea_all_results_XBP1_RIDD_significatif %>%
    mutate(Direction = ifelse(NES > 0, "UP", "DOWN"),
           BaseName = paste0(Pathway, "_", Direction)) %>%
    group_by(BaseName) %>%
    mutate(Nom = paste0(BaseName, "_", row_number())) %>%
    ungroup() %>%
    select(Component, Nom)
  
  score <- as.data.frame(score)
  
  # On ne garde que les composantes qu'on peut associer aux signatures
  score_long_XBP1_RIDD <- score %>%
    tibble::rownames_to_column(var = "Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Component",
      values_to = "Abundance"
    ) %>%
    # Adaptation du nom des composantes pour matcher
    mutate(Component = paste0("ICA_", gsub("IC", "", Component))) 
  
  score_long_joined_XBP1_RIDD <- score_long_XBP1_RIDD %>%
    left_join(correspondance_Comp_XBP1_RIDD, by = "Component", relationship = "many-to-many") %>%
    filter(!is.na(Nom))
  
  # IDEE : Additionner le score des meme pathways ensemble 
  
  score_grouped_XBP1_RIDD <- score_long_joined_XBP1_RIDD %>%
    mutate(BaseName = sub("_[0-9]+$", "", Nom)) %>%  # supprime le suffixe _1, _2, etc.
    group_by(Sample, BaseName) %>%
    summarise(Score = sum(Abundance), .groups = "drop")
  
  
  
  # Filtrage des données
  filtered_data_38 <- reactive({
    req(input$selected_basenames)
    score_grouped_XBP1_RIDD[score_grouped_XBP1_RIDD$BaseName %in% input$selected_basenames, ]
  })
  
  output$result_ab_IRE138_graph <- renderPlotly({
    p_ab_38<-ggplot(filtered_data_38(), aes(x = Sample, y = Score, fill = BaseName)) +
      geom_bar(stat = "identity") +
      labs(
        title = "Scores d'abondance par pathway d'IRE1 38 au sein des échantillons TCGA",
        x = "Échantillon",
        y = "Score d'abondance"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    ggplotly(p_ab_38)
  })
  
  fgsea_all_results_XBP1_RIDD_pur_significatif <- filter(fgsea_all_results_XBP1_RIDD_pur,padj<0.05)
  
  correspondance_Comp_XBP1_RIDD_pur <- fgsea_all_results_XBP1_RIDD_pur_significatif %>%
    mutate(Direction = ifelse(NES > 0, "UP", "DOWN"),
           BaseName = paste0(Pathway, "_", Direction)) %>%
    group_by(BaseName) %>%
    mutate(Nom = paste0(BaseName, "_", row_number())) %>%
    ungroup() %>%
    select(Component, Nom)
  
  score <- as.data.frame(score)
  
  # On ne garde que les composantes qu'on peut associer aux signatures
  score_long_XBP1_RIDD_pur <- score %>%
    tibble::rownames_to_column(var = "Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Component",
      values_to = "Abundance"
    ) %>%
    # Adaptation du nom des composantes pour matcher
    mutate(Component = paste0("ICA_", gsub("IC", "", Component))) 
  
  score_long_joined_XBP1_RIDD_pur <- score_long_XBP1_RIDD_pur %>%
    left_join(correspondance_Comp_XBP1_RIDD_pur, by = "Component",relationship = "many-to-many") %>%
    filter(!is.na(Nom))
  
  # IDEE : ADDITIONNER LE SCORE Des meme pathwayensemble 
  
  score_grouped_XBP1_RIDD_pur<- score_long_joined_XBP1_RIDD_pur %>%
    mutate(BaseName = sub("_[0-9]+$", "", Nom)) %>%  # supprime le suffixe _1, _2, etc.
    group_by(Sample, BaseName) %>%
    summarise(Score = sum(Abundance), .groups = "drop")
  
  ggplot(score_grouped_XBP1_RIDD_pur, aes(x = Sample, y = Score, fill = BaseName)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Scores d'abondance par pathway",
      x = "Échantillon",
      y = "Score d'abondance"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
  
  # Filtrage des données
  filtered_data_pur <- reactive({
    req(input$selected_basenames)
    score_grouped_XBP1_RIDD_pur[score_grouped_XBP1_RIDD_pur$BaseName %in% input$selected_basenames, ]
  })
  
  output$result_ab_IRE1pur_graph <- renderPlotly({
    p_ab_pur<-ggplot(filtered_data_pur(), aes(x = Sample, y = Score, fill = BaseName)) +
      geom_bar(stat = "identity") +
      labs(
        title = "Scores d'abondance pour XBP1s pur et RIDD pur au sein des échantillons TCGA",
        x = "Échantillon",
        y = "Score d'abondance"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    ggplotly(p_ab_pur)
  })
  
  score_grouped_XBP1_RIDD_wide <- score_grouped_XBP1_RIDD %>%
    pivot_wider(names_from = BaseName, values_from = Score)
  
  score_grouped_XBP1_RIDD_pur_wide <- score_grouped_XBP1_RIDD_pur %>%
    pivot_wider(names_from = BaseName, values_from = Score)
  
  score_XBP1_RIDD_pur_impur <-merge(score_grouped_XBP1_RIDD_pur_wide,score_grouped_XBP1_RIDD_wide,by = "Sample")
  
  output$acp_ire1 <- renderPlot({
    
    acp <- PCA(score_XBP1_RIDD_pur_impur[, -1], scale.unit = TRUE, graph = TRUE)
    fviz_pca_ind(acp)  # individus
    p <- fviz_pca_var(acp, 
                      repel = TRUE,    # évite que les labels se chevauchent
                      col.var = "black")  # couleur unie pour plus de lisibilité
    
    # Rendre interactif avec plotly
    p
  })
  
  numeric_data_H <- score_XBP1_RIDD_pur_impur[, -1]
  data_scaled_H <- scale(numeric_data_H)
  
  # Déterminer le nombre optimal de clusters en utilisant l'indice de silhouette
  silhouette_info <- fviz_nbclust(data_scaled_H, kmeans, method = "silhouette")
  nb_groupe_silhouette <- which.max(silhouette_info$data$y)
  I.intra = sapply(1:20,FUN=function(k) kmeans(numeric_data,centers=k,nstart=50)$tot.withinss)
  
  diff_I.intra <- diff(I.intra)
  diff_ratio <- diff_I.intra / I.intra[1:(length(I.intra) - 1)]
  nb_groupe_elbow <- which.max(diff_ratio) + 1
  plot(I.intra,type="b",xlab="nb groupes",ylab="inertie intra") 
  res_nbclust <- NbClust(data_scaled_H, distance = "euclidean", 
                         min.nc = 2, max.nc = 15, method = "kmeans", index = "all")
  
  # Le nombre optimal de clusters selon la majorité des indices
  nb_groupe <- res_nbclust$Best.nc[1]
  
  
  output$k_XR <- renderText({nb_groupe})
  
  output$kmeans_ire1 <-renderPlot({
    numeric_data_pur <- score_XBP1_RIDD_pur_impur[, -1]  # Exclude the first column (Sample)
    
    set.seed(123)  # Pour la reproductibilité
    km <- kmeans(numeric_data_pur, centers = nb_groupe, nstart = 50)
    
    clusplot(numeric_data_pur,km$cluster,color = TRUE, shade = TRUE, 
             labels = nb_groupe, lines = 0, main = "Clusplot XBP1s RIDD kmeans")
  })
  nb_groupe_XR_CAH <- NbClust(data_scaled_H, 
                              distance = "euclidean", 
                              min.nc = 2, 
                              max.nc = 10, 
                              method = "ward.D2",   # ou "complete", "average", etc.
                              index = "all") 
  
  output$k_XR_CAH <- renderText({
    nb_groupe_XR_CAH$Best.nc[1]
  })
  output$cah_SEA <- renderPlot({
    d <- dist(score_grouped_XBP1_RIDD_wide[, -1], method = "euclidean")
    h <- hclust(d,method = "complete")
    plot(h, main = "Dendrogramme de la CAH XBP1 et RIDD")
    rect.hclust(h, k = nb_groupe_XR_CAH$Best.nc[1], border = 'darkorchid')
    groupes_cah <- cutree(h, k = nb_groupe_XR_CAH$Best.nc[1])
    clusplot(score_grouped_XBP1_RIDD_wide[, -1], groupes_cah, color = TRUE, shade = TRUE, 
             labels = nb_groupe_XR_CAH$Best.nc[1], lines = 0, main = "Clusplot XBP1s RIDD CAH")
  })
########################################################TAB CANCER SEA####################################################################
  output$exemple_GSEA_CS <- renderUI({
    HTML(
      ex_gsea %>%
        kable(format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "orchid", bold = TRUE) %>%
        as.character()
    )
  })
  
  
  # Convertir CancerSEA_sig (data.frame ou matrice) en liste de data.frames
  CancerSEA_sig_l <- lapply(CancerSEA_sig, function(genes) {
    genes <- na.omit(genes)
    data.frame(gene = genes)
  })
  
  # Attribuer les noms
  names(CancerSEA_sig_l) <- colnames(CancerSEA_sig)
  
  basis.list_CancerSEA <- CancerSEA_sig_l
  
  all_fgsea_CancerSEA <- list()
  
  # S'assurer que les noms des stats sont corrects
  TCGA_ica$names <- rownames(TCGA_final)
  
  for (i in 1:ncol(TCGA_ica$S)) {
    component_name <- paste0("ICA_", i)
    
    stats <- TCGA_ica$S[, i]
    names(stats) <- TCGA_ica$names
    stats <- stats[!is.na(stats)]  # filtrer les NA
    
    for (pathway_name in names(basis.list_CancerSEA)) {
      # Extraire les gènes de la signature
      genes <- basis.list_CancerSEA[[pathway_name]]$gene
      
      # Créer un sous-ensemble de pathways comme liste
      pathway_list <- list()
      pathway_list[[pathway_name]] <- genes
      
      # Vérifier qu'il y a assez de gènes
      # On fait le gsea
      if (length(intersect(names(stats), genes)) >= 5) {
        fgsea_res_CancerSEA <- fgsea(
          pathways = pathway_list,
          stats = stats,
          minSize = 5,
          maxSize = 20000
        )
        
        # Ajouter le nom de la composante et du pathway
        fgsea_res_CancerSEA$Component <- component_name
        fgsea_res_CancerSEA$Pathway <- pathway_name
        
        all_fgsea_CancerSEA[[length(all_fgsea_CancerSEA) + 1]] <- fgsea_res_CancerSEA
      }
    }
  }
  
  # Combiner en un seul data.frame
  # On met les résultats du gsea dans un tableau Kable pour une éventuelle meilleur lecture des résultats
  fgsea_all_results_CancerSEA <- bind_rows(all_fgsea_CancerSEA)
  
  output$result_GSEA_SEA <- renderUI({
    HTML(fgsea_all_results_CancerSEA %>%
           kable(format = "html", 
                 escape = FALSE) %>%
           kable_styling(bootstrap_options = c("condensed"),
                         full_width = TRUE,
                         position = "center") %>%
           scroll_box(width = "100%", height = "500px") %>%
           row_spec(0,color = "orchid", bold = TRUE) %>%
           as.character())
  })
  
  output$result_GSEA_SEA_graph <- renderPlotly({
    p_CancerSEA <- ggplot(fgsea_all_results_CancerSEA, aes(x = Component, y = Pathway)) +
      geom_point(aes(size = -log10(padj), color = NES)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size(range = c(1, 6)) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(color = "NES", size = "-log10(padj)")
    ggplotly(p_CancerSEA)
  })
  
  fgsea_all_results_CancerSEA_significatif <- filter(fgsea_all_results_CancerSEA,padj<0.05)
  
  correspondance_Comp_CancerSEA <- fgsea_all_results_CancerSEA_significatif %>%
    mutate(Direction = ifelse(NES > 0, "UP", "DOWN"),
           BaseName = paste0(Pathway, "_", Direction)) %>%
    group_by(BaseName) %>%
    mutate(Nom = paste0(BaseName, "_", row_number())) %>%
    ungroup() %>%
    select(Component, Nom)
  
  score <- as.data.frame(score)
  
  # On ne garde que les composantes qu'on peut associer aux signatures
  score_long_CancerSEA <- score %>%
    tibble::rownames_to_column(var = "Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Component",
      values_to = "Abundance"
    ) %>%
    # Adaptation du nom des composantes pour matcher
    mutate(Component = paste0("ICA_", gsub("IC", "", Component))) 
  
  score_long_joined_CancerSEA <- score_long_CancerSEA %>%
    left_join(correspondance_Comp_CancerSEA, by = "Component",relationship = "many-to-many") %>%
    filter(!is.na(Nom))
  
  
  # IDEE : ADDITIONNER LE SCORE Des meme pathwayensemble 
  
  score_grouped_CancerSEA <- score_long_joined_CancerSEA %>%
    mutate(BaseName = sub("_[0-9]+$", "", Nom)) %>%  # supprime le suffixe _1, _2, etc.
    group_by(Sample, BaseName) %>%
    summarise(Score = sum(Abundance), .groups = "drop")
  
  all_basenames_SEA <- unique(score_grouped_CancerSEA$BaseName)
  # Observers pour les boutons
  observeEvent(input$select_all, {
    updatePickerInput(session, "selected_basenames",
                             selected = all_basenames_SEA)
  })
  
  observeEvent(input$clear_all, {
    updatePickerInput(session, "selected_basenames",
                             selected = character(0))
  })
  
  observeEvent(input$select_up, {
    up_labels <- all_basenames_SEA[str_detect(all_basenames_SEA, "UP$")]
    updatePickerInput(session, "selected_basenames",
                             selected = up_labels)
  })
  
  observeEvent(input$select_down, {
    down_labels <- all_basenames_SEA[str_detect(all_basenames_SEA, "DOWN$")]
    updatePickerInput(session, "selected_basenames",
                             selected = down_labels)
  })
  
  # Filtrage des données
  filtered_data_SEA <- reactive({
    sel <- input$selected_basenames
    if (is.null(sel) || length(sel) == 0) {
      return(NULL)
    }
    score_grouped_CancerSEA[score_grouped_CancerSEA$BaseName %in% sel, ]
  })
  
  
  # Plot interactif
  output$result_ab_SEA_graph <- renderPlotly({
    data <- filtered_data_SEA()
    
    if (is.null(data) || nrow(data) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = "Aucune donnée sélectionnée"))
    }
    
    p <- ggplot(data, aes(x = Sample, y = Score, fill = BaseName)) +
      geom_bar(stat = "identity") +
      labs(
        title = "Scores d'abondance par pathway",
        x = "Échantillon",
        y = "Score d'abondance"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    ggplotly(p)
  })
  
  score_grouped_CancerSEA_wide <- score_grouped_CancerSEA %>%
    pivot_wider(names_from = BaseName, values_from = Score)
  acp <- PCA(score_grouped_CancerSEA_wide[, -1], scale.unit = TRUE, graph = FALSE)
  
  output$ind_SEA <- renderPlot({fviz_pca_ind(acp)}) # individus
    
  output$var_SEA <- renderPlot({fviz_pca_var(acp)})  # variables
  
  numeric_data_SEA <- score_grouped_CancerSEA_wide[, -1]
  data_scaled_SEA <- scale(numeric_data)
  
  # Déterminer le nombre optimal de clusters en utilisant l'indice de silhouette
  res_nbclust_SEA <- NbClust(numeric_data_SEA, distance = "euclidean", 
                         min.nc = 2, max.nc = 15, method = "kmeans", index = "silhouette")
  
  # Le nombre optimal de clusters selon la majorité des indices
  nb_groupe_SEA <- res_nbclust_SEA$Best.nc[1]
  
  output$k_SEA <- renderText({
    nb_groupe_SEA
  })
  
  output$kmeans_SEA <- renderPlot({
    numeric_data <- score_grouped_CancerSEA_wide[, -1]  # Exclude the first column (Sample)
    
    set.seed(123)  # Pour la reproductibilité
    km <- kmeans(numeric_data_SEA, centers = nb_groupe_SEA, nstart = 50)
    print(km)
    clusplot(numeric_data_SEA,km$cluster,color = TRUE, shade = TRUE, 
             labels = nb_groupe_SEA, lines = 0, main = "Clusplot pathways de CancerSEA kmeans")
  })
  
  k_CAH_SEA <- NbClust(data_scaled_SEA, 
                       distance = "euclidean", 
                       min.nc = 2, 
                       max.nc = 10, 
                       method = "ward.D2",   # ou "complete", "average", etc.
                       index = "silhouette") 
  
  nb_groupe_SEA_CAH <- k_CAH_SEA$Best.nc[1] 
  
  output$k_SEA_CAH <- renderText({
    nb_groupe_SEA_CAH
  })
    
  output$cah_ire1 <- renderPlot({
    d <- dist(score_grouped_CancerSEA_wide[, -1], method = "euclidean")
    h <- hclust(d,method = "complete")
    plot(h, main = "Dendrogramme de la CAH de CancerSEA")
    rect.hclust(h, k = nb_groupe_SEA_CAH, border = 'darkorchid')
    groupes_cah <- cutree(h, k = nb_groupe_SEA_CAH)
    clusplot(score_grouped_CancerSEA_wide[, -1], groupes_cah, color = TRUE, shade = TRUE, 
             labels = nb_groupe_SEA_CAH, lines = 0, main = "Clusplot de CancerSEA CAH")
  })
  
  
########################################################TAB HALLMARKS####################################################################
  output$exemple_GSEA_HM <- renderUI({
    HTML(
      ex_gsea %>%
        kable(format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "orchid", bold = TRUE) %>%
        as.character()
    )
  })
  
  Hallmark_sig_l <- lapply(Hallmark_sig, function(genes) {
    genes <- na.omit(genes)
    data.frame(gene = genes)
  })
  
  # Attribuer les noms
  names(Hallmark_sig_l) <- colnames(Hallmark_sig)
  
  basis.list_Hallmark <- Hallmark_sig_l
  
  all_fgsea_Hallmark <- list()
  # S'assurer que les noms des stats sont corrects
  TCGA_ica$names <- rownames(TCGA_final)
  for (i in 1:ncol(TCGA_ica$S)) {
    component_name <- paste0("ICA_", i)
    
    stats <- TCGA_ica$S[, i]
    names(stats) <- TCGA_ica$names
    stats <- stats[!is.na(stats)]  # filtrer les NA
    
    for (pathway_name in names(basis.list_Hallmark)) {
      # Extraire les gènes de la signature
      genes <- basis.list_Hallmark[[pathway_name]]$gene
      
      # Créer un sous-ensemble de pathways comme liste
      pathway_list <- list()
      pathway_list[[pathway_name]] <- genes
      
      # Vérifier qu'il y a assez de gènes
      if (length(intersect(names(stats), genes)) >= 5) {
        fgsea_res_Hallmark <- fgsea(
          pathways = pathway_list,
          stats = stats,
          minSize = 5,
          maxSize = 20000
        )
        
        # Ajouter le nom de la composante et du pathway
        fgsea_res_Hallmark$Component <- component_name
        fgsea_res_Hallmark$Pathway <- pathway_name
        
        all_fgsea_Hallmark[[length(all_fgsea_Hallmark) + 1]] <- fgsea_res_Hallmark
      }
    }
  }
  # Combiner en un seul data.frame
  fgsea_all_results_Hallmark <- bind_rows(all_fgsea_Hallmark)
  
  output$result_GSEA_Hallmark <- renderUI({
    HTML(fgsea_all_results_Hallmark %>%
           kable(format = "html", 
                 escape = FALSE) %>%
           kable_styling(bootstrap_options = c("condensed"),
                         full_width = TRUE,
                         position = "center") %>%
           scroll_box(width = "100%", height = "500px") %>%
           row_spec(0,color = "orchid", bold = TRUE) %>%
           as.character())
  })
  
  output$result_GSEA_Hallmark_graph <- renderPlotly({
    p_Hallmark <- ggplot(fgsea_all_results_Hallmark, aes(x = Component, y = Pathway)) +
      geom_point(aes(size = -log10(padj), color = NES)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size(range = c(1, 6)) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(color = "NES", size = "-log10(padj)")
    ggplotly(p_Hallmark)
  })
  
  fgsea_all_results_Hallmark_significatif <- filter(fgsea_all_results_Hallmark,padj<0.05)
  
  correspondance_Comp_Hallmark <- fgsea_all_results_Hallmark_significatif %>%
    mutate(Direction = ifelse(NES > 0, "UP", "DOWN"),
           BaseName = paste0(Pathway, "_", Direction)) %>%
    group_by(BaseName) %>%
    mutate(Nom = paste0(BaseName, "_", row_number())) %>%
    ungroup() %>%
    select(Component, Nom)
  
  score <- as.data.frame(score)
  
  # On ne garde que les composantes qu'on peut associer aux signatures
  score_long_Hallmark <- score %>%
    tibble::rownames_to_column(var = "Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Component",
      values_to = "Abundance"
    ) %>%
    # Adaptation du nom des composantes pour matcher
    mutate(Component = paste0("ICA_", gsub("IC", "", Component))) 
  
  score_long_joined_Hallmark <- score_long_Hallmark %>%
    left_join(correspondance_Comp_Hallmark, by = "Component",relationship = "many-to-many") %>%
    filter(!is.na(Nom))
  
  # IDEE : ADDITIONNER LE SCORE Des meme pathwayensemble 
  
  score_grouped_Hallmark<- score_long_joined_Hallmark %>%
    mutate(BaseName = sub("_[0-9]+$", "", Nom)) %>%  # supprime le suffixe _1, _2, etc.
    group_by(Sample, BaseName) %>%
    summarise(Score = sum(Abundance), .groups = "drop")
  
  # Exemple : unique basenames utilisés pour les boutons
  all_basenames_Hallmarks <- unique(score_grouped_Hallmark$BaseName)

  observeEvent(input$select_all, {
    updatePickerInput(session, "selected_basenames",
                             selected = all_basenames_Hallmarks)
  })
  
  observeEvent(input$clear_all, {
    updatePickerInput(session, "selected_basenames",
                             selected = character(0))
  })
  
  observeEvent(input$select_up, {
    up_labels <- all_basenames_Hallmarks[str_detect(all_basenames_Hallmarks, "UP$")]
    updatePickerInput(session, "selected_basenames",
                             selected = up_labels)
  })
  
  observeEvent(input$select_down, {
    down_labels <- all_basenames_Hallmarks[str_detect(all_basenames_Hallmarks, "DOWN$")]
    updatePickerInput(session, "selected_basenames",
                             selected = down_labels)
  })
  
  # Filtrage des données
  filtered_data_Hallmarks <- reactive({
    req(input$selected_basenames)
    score_grouped_Hallmark[score_grouped_Hallmark$BaseName %in% input$selected_basenames, ]
  })
  
  # Plot interactif
  output$result_ab_Hallmarks_graph <- renderPlotly({
    p <- ggplot(filtered_data_Hallmarks(), aes(x = Sample, y = Score, fill = BaseName)) +
      geom_bar(stat = "identity") +
      labs(
        title = "Scores d'abondance par pathway",
        x = "Échantillon",
        y = "Score d'abondance"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    ggplotly(p)
  })
  
  score_grouped_Hallmark_wide <- score_grouped_Hallmark %>%
    pivot_wider(names_from = BaseName, values_from = Score)
  acp <- PCA(score_grouped_Hallmark_wide[, -1], scale.unit = TRUE, graph = FALSE)
  
  output$ind_Hallmarks <- renderPlot({fviz_pca_ind(acp)})
  
  output$var_Hallmarks <- renderPlot({fviz_pca_var(acp)})
  
  numeric_data_H <- score_grouped_Hallmark_wide[, -1]
  data_scaled_H <- scale(numeric_data_H)
  
  # Déterminer le nombre optimal de clusters en utilisant l'indice de silhouette
  fviz_nbclust(data_scaled_H, kmeans, method = "silhouette")
  I.intra = sapply(1:20,FUN=function(k) kmeans(numeric_data,centers=k,nstart=50)$tot.withinss)
  plot(I.intra,type="b",xlab="nb groupes",ylab="inertie intra")
  
  numeric_data_H <- score_grouped_Hallmark_wide[, -1]
  data_scaled_H <- scale(numeric_data_H)
  
  res_nbclust_H <- NbClust(data_scaled_H, distance = "euclidean", 
                           min.nc = 2, max.nc = 20, method = "kmeans", index = "silhouette")
  
  # Le nombre optimal de clusters selon la majorité des indices
  nb_groupe_H <- res_nbclust_H$Best.nc[1]
  output$k_H <- renderText({nb_groupe_H})
  
  output$kmeans_Hallmarks <- renderPlot({
    numeric_data_H <- score_grouped_Hallmark_wide[, -1]  # Exclude the first column (Sample)
    
    set.seed(123)  # Pour la reproductibilité
    km <- kmeans(numeric_data_H, centers = nb_groupe_H, nstart = 50)
    print(km)
    clusplot(numeric_data_H,km$cluster,color = TRUE, shade = TRUE, 
             labels = nb_groupe_H, lines = 0, main = "Clusplot pathways kmeans")
  })
  
  nb_groupe_H_CAH <- NbClust(data_scaled_H, 
                              distance = "euclidean", 
                              min.nc = 2, 
                              max.nc = 20, 
                              method = "ward.D2",   
                              index = "silhouette") 
  
  output$k_H_CAH <- renderText({
    nb_groupe_H_CAH$Best.nc[1]
  })
  
  output$cah_Hallmarks <- renderPlot({
    d <- dist(score_grouped_Hallmark_wide[, -1], method = "euclidean")
    h <- hclust(d,method = "complete")
    plot(h, main = "Dendrogramme de la CAH des Hallmark")
    rect.hclust(h, k = nb_groupe_H_CAH$Best.nc[1], border = 'darkorchid')
    groupes_cah <- cutree(h, k = nb_groupe_H_CAH$Best.nc[1])
    clusplot(score_grouped_Hallmark_wide[, -1], groupes_cah, color = TRUE, shade = TRUE, 
             labels = nb_groupe_H_CAH$Best.nc[1], lines = 0, main = "Clusplot de Hallmark CAH")
  })
  
########################################################TAB QUALITE SIGNATURE####################################################################
  
  df <- score_XBP1_RIDD_pur_impur
  
  score_XBP1_RIDD_pur_impur_long <- score_XBP1_RIDD_pur_impur %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Pathway_Direction",
      values_to = "Score"
    )
  
  # Séparer en deux colonnes : Pathway et Direction
  score_XBP1_RIDD_pur_impur_long <- score_XBP1_RIDD_pur_impur_long %>%
    separate(Pathway_Direction, into = c("Pathway", "Direction"), sep = "_(?=UP|DOWN)", remove = FALSE)
  
  # Supprimer la colonne Pathway_Direction
  score_XBP1_RIDD_pur_impur_long <- score_XBP1_RIDD_pur_impur_long %>% select(-Pathway_Direction)
  
  # Créer le tableau de contingence
  table_contingence <- table(score_XBP1_RIDD_pur_impur_long$Direction, score_XBP1_RIDD_pur_impur_long$Pathway)
  
  # Affichage
  output$contingence_occurence <- renderUI({
    HTML(
      table_contingence %>%
        kable(format = "html",
              escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "white", bold = TRUE) %>%
        as.character()
    )
  })
  
  
  # Reshape en long
  score_XBP1_RIDD_pur_impur_long <- score_XBP1_RIDD_pur_impur %>%
    pivot_longer(cols = -Sample, names_to = "Pathway_Direction", values_to = "Score") %>%
    separate(Pathway_Direction, into = c("Pathway", "Direction"), sep = "_(?=UP|DOWN)")
  
  # Créer la table de contingence avec la **somme des scores**
  table_sum_scores <- score_XBP1_RIDD_pur_impur_long %>%
    group_by(Direction, Pathway) %>%
    summarise(Sum_Score = sum(Score), .groups = "drop") %>%
    pivot_wider(names_from = Pathway, values_from = Sum_Score)
  
  # Affichage
  output$contingency_scores<- renderUI({
    HTML(
      table_sum_scores %>%
        kable(format = "html", escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "white", bold = TRUE) %>%
        as.character()
    )
  })
  
  df_active_38 <- score_grouped_XBP1_RIDD %>%
    filter(BaseName %in% c("XBP1_UP", "RIDD_DOWN")) %>%
    group_by(Sample) %>%
    summarise(IRE1_actif = sum(Score, na.rm = TRUE))
  
  # on range tout en fonction du ratio croissant pour une question de lisibilité
  df_active_38 <-df_active_38  %>%
    arrange(IRE1_actif) %>%
    mutate(Sample = factor(Sample, levels = Sample))
  
  p_act_38<- ggplot(df_active_38, aes(x = Sample, y = IRE1_actif, group = 1)) +
    geom_point( size = 2, color ="darkorchid") +
    geom_line() +
    theme_minimal() +
    labs(title = "score de l'activité d'IRE1 par échantillon basé sur les signatures XBP1s et RIDD issues de la signature gloable d'IRE1 38",
         x = "Échantillon",
         y = "score de l'activité d'IRE1") +
    theme(axis.text.x = element_blank())
  
  output$activite_38 <- renderPlotly({
    ggplotly(p_act_38)
  })
  
  df_active_pur <- score_grouped_XBP1_RIDD_pur %>%
    filter(BaseName %in% c("XBP1_pur_UP", "RIDD_pur_DOWN")) %>%
    group_by(Sample) %>%
    summarise(IRE1_actif = sum(Score, na.rm = TRUE))
  
  # on range tout en fonction du ratio croissant pour une question de lisibilité
  df_active_pur <-df_active_pur  %>%
    arrange(IRE1_actif) %>%
    mutate(Sample = factor(Sample, levels = Sample))
  
  p_act_pur <- ggplot(df_active_pur, aes(x = Sample, y = IRE1_actif, group = 1)) +
    geom_point( size = 2, color ="darkorchid") +
    geom_line() +
    theme_minimal() +
    labs(title = "score de l'activité d'IRE1 par échantillon basés sur les signatures XBP1s pur et RIDD pur",
         x = "Échantillon",
         y = "score de l'activitéd'IRE1") +
    theme(axis.text.x = element_blank())
  
  output$activite_pur <- renderPlotly({
    ggplotly(p_act_pur)
  })
  
  # Créer un dataframe combiné pour les deux ensembles de données
  df_active_38$Group <- "Groupe IRE1 38"
  df_active_pur$Group <- "Groupe Pur"
  combined_df <- rbind(df_active_38, df_active_pur)
  combined_df$Group <- factor(combined_df$Group, 
                              levels = c("Groupe IRE1 38", "Groupe Pur"),
                              labels = c("XBP1s 19 et RIDD 19", "XBP1s pur et RIDD pur"))
  
  # Créer le boxplot
  output$comparaison_activite <- renderPlot({
    ggplot(combined_df, aes(x = Group, y = IRE1_actif, fill = Group)) +
      geom_boxplot() +
      labs(
        title = "Comparaison des scores IRE1_actif entre les deux types de signatures",
        x = "Type de signature",
        y = "Score IRE1_actif",
        fill = "Type de signature"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right"
      )
  })
#################TRUC MATHEO########################
  score_XR_pur <- read.csv(file = "score_grouped_XBP1_RIDD_pur.csv", sep = ",")
  score_XR_pur <- score_XR_pur[,-1]
  colnames(score_XR_pur) <- c("Sample_pur", "BaseName_pur", "Score_pur")
  score_XR <- read.csv(file = "score_grouped_XBP1_RIDD.csv", sep = ",")
  score_XR <- score_XR[,-1]
  signature_genes <- read.table(file = "annotation_allprobes.txt", header = T)
  
  #load("TCGA.RData") # tcga_cohort_stupp_xiang_226 tcga_peset_554 TCGA_counts tcga_geset_554 xref_tcga_peset_554
  score_XR <- score_XR[score_XR$BaseName != "RIDD_UP",]
  
  scores <- cbind(score_XR, score_XR_pur)
  
  
  
  scores_transformed_impur <- scores %>%
    pivot_wider(
      id_cols = Sample,
      names_from = BaseName,
      values_from = c(Score)
    )
  scores_transformed_impur$status = "impur"
  scores_transformed_pur <- scores %>%
    pivot_wider(
      id_cols = Sample,
      names_from = BaseName,
      values_from = c(Score_pur)
    )
  scores_transformed_pur$status = "pur"
  scores_transformed <- rbind(scores_transformed_pur,scores_transformed_impur)
  
  XBP1_DOWN_l <- score_XR$Score[score_XR$BaseName == "XBP1_DOWN"]
  XBP1_UP_l <- score_XR$Score[score_XR$BaseName == "XBP1_UP"]
  
  
  
  XBP1_DOWN_pur_l <- score_XR_pur$Score[score_XR_pur$BaseName == "XBP1_pur_DOWN"]
  XBP1_UP_pur_l <- score_XR_pur$Score[score_XR_pur$BaseName == "XBP1_pur_UP"]
  
  
  
  
  ## Quand les deux vecteur ont les rang dans le mm ordre c'est -1
  ## Quand les deux vecteur ont les rang dans l'ordre inverse c'est 1
  # Si la valeur de l'indice est proche de 0 cela signifie que les rangs sont sans lien clair 
  
  
  
  
  indice_inversion_vecteurs <- function(v1, v2) {
    if (length(v1) != length(v2)) {
      stop("Les deux vecteurs doivent avoir la même longueur.")
    }
    
    r1 <- rank(-v1, ties.method = "min")
    r2 <- rank(-v2, ties.method = "min")
    
    indice <- -cor(r1, r2, method = "spearman")
    return(indice)
  }
  ind_X_19 <- indice_inversion_vecteurs(XBP1_DOWN_l, XBP1_UP_l)
 
  output$my_box1 <- renderText({
    (ind_X_19)
  })
  
  ind_X_pur <- indice_inversion_vecteurs(XBP1_DOWN_pur_l, XBP1_UP_pur_l)
  output$my_box2 <- renderText({
    ind_X_pur
  })
  
  XBP1_DOWN_l_dc <- quantile(XBP1_DOWN_l, probs = seq(0,1,0.1))
  XBP1_UP_l_dc <- quantile(XBP1_UP_l, probs = seq(0,1,0.1))
  
  ind_X_19_dc <- indice_inversion_vecteurs(XBP1_DOWN_l_dc, XBP1_UP_l_dc)
  
  #output$my_box1_dc <- renderText({
  #  (ind_X_19_dc)
  #})
  
  XBP1_DOWN_pur_l_dc <- quantile(XBP1_DOWN_pur_l, probs = seq(0,1,0.1))
  XBP1_UP_pur_l_dc <- quantile(XBP1_UP_pur_l, probs = seq(0,1,0.1))
  
  ind_X_pur_dc <- indice_inversion_vecteurs(XBP1_DOWN_pur_l_dc, XBP1_UP_pur_l_dc)
  
  #output$my_box2_dc <- renderText({
  #  ind_X_pur_dc
  #})
  
  df_inv_rang_38 <- data.frame("DOWN" = XBP1_DOWN_l,
                      "UP"= XBP1_UP_l)
  output$inversion_rang_38 <- renderPlot(
    ggparcoord(df_inv_rang_38,
                  columns = 1:2, 
               showPoints = TRUE,
               mapping = aes(color = NULL)) + 
      geom_line(color = "orchid2") +
      geom_point(color = "orchid2") +
      theme_minimal()
  )
  
  df_inv_rang_pur <- data.frame("DOWN" = XBP1_DOWN_pur_l,
                                "UP" = XBP1_UP_pur_l)
  
  output$inversion_rang_pur <- renderPlot(
    ggparcoord(df_inv_rang_pur,
                  columns = 1:2, 
                  showPoints = TRUE,
                  mapping = aes(color = NULL)) + 
      geom_line(color = "orchid2") +
      geom_point(color = "orchid2") +
      theme_minimal()
  )
  
  #df_inv_rang_38_dc <- data.frame("DOWN" = XBP1_DOWN_l_dc,
  #                             "UP"= XBP1_UP_l_dc)
  #output$inversion_rang_38_dc <- renderPlot(
  #  ggparcoord(df_inv_rang_38_dc,
  #                columns = 1:2, 
  #                showPoints = TRUE,
  #                mapping = aes(color = NULL)) + 
  #    geom_line(color = "orchid2") +
  #    geom_point(color = "orchid2") +
  #    theme_minimal()
  #)
  
  df_inv_rang_pur_dc <- data.frame("DOWN" = XBP1_DOWN_pur_l_dc,
                                  "UP"= XBP1_UP_pur_l_dc)
  
  #output$inversion_rang_pur_dc <- renderPlot(
  #  ggparcoord(df_inv_rang_pur_dc,
  #             columns = 1:2, 
  #             showPoints = TRUE,
  #             mapping = aes(color = NULL)) + 
  #    geom_line(color = "orchid2") +
  #    geom_point(color = "orchid2") +
  #    theme_minimal()
  #)
  
  # Filter the data
  filtered_dt <- fgsea_combined[padj < 0.05 & abs(NES) > 2,]
  
  
  
  # Étape 1 : filtrer les pathways avec NES fort et padj significatif
  filtered_df <- fgsea_combined %>%
    filter(abs(NES) > 1, padj < 0.05)
  
  # Étape 2 : passer au format large pour comparaison
  wide_df <- filtered_df %>%
    pivot_wider(
      id_cols = Component,
      names_from = Pathway,
      values_from = NES
    )
  
  # Étape 3 : appliquer la logique conditionnelle
  valid_components <- wide_df %>%
    filter(
      (!is.na(XBP1_pur) & !is.na(RIDD_pur) &
         ((XBP1_pur > 0 & RIDD_pur < 0) | (XBP1_pur < 0 & RIDD_pur > 0))) |
        (!is.na(XBP1) & !is.na(RIDD) &
           ((XBP1 > 0 & RIDD < 0) | (XBP1 < 0 & RIDD > 0)))
    ) %>%
    pull(Component)
  
  # Étape 4 : extraire les lignes complètes du jeu filtré
  result <- filtered_df %>%
    filter(Component %in% valid_components)
  
  
  output$gsea_ire1_actif <- renderPlotly({
    p<-ggplot(result, aes(x = Component, y = Pathway)) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size(range = c(1, 6)) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(color = "NES", size = "-log10(padj)")
    ggplotly(p)
  })
  #on associe à chaque composante la pathway adapté (on la associé grace au test de GSEA )
  # Alors comme des signatures peuvent être associé à plusieurs composantes on va nommer les signatures de la facon suivante (sig_UP_1,sig_UP_2,...)
  correspondance_Comp_XBP1_RIDD_pur_19 <- result %>%
    mutate(Direction = ifelse(NES > 0, "UP", "DOWN"),
           BaseName = paste0(Pathway, "_", Direction)) %>%
    group_by(BaseName) %>%
    mutate(Nom = paste0(BaseName, "_", row_number())) %>%
    ungroup() %>%
    select(Component, Nom)
  
  score <- as.data.frame(score)
  
  # On ne garde que les composantes qu'on peut associer aux signatures
  score_long_XBP1_RIDD_pur_19  <- score %>%
    tibble::rownames_to_column(var = "Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "Component",
      values_to = "Abundance"
    ) %>%
    # Adaptation du nom des composantes pour matcher
    mutate(Component = paste0("ICA_", gsub("IC", "", Component))) 
  
  score_long_joined_XBP1_RIDD_pur_19  <- score_long_XBP1_RIDD_pur_19  %>%
    left_join(correspondance_Comp_XBP1_RIDD_pur_19 , by = "Component") %>%
    filter(!is.na(Nom))
  
  # IDEE : Additionner le score des meme pathways ensemble 
  
  score_grouped_XBP1_RIDD_pur_19  <- score_long_joined_XBP1_RIDD_pur_19  %>%
    mutate(BaseName = sub("_[0-9]+$", "", Nom)) %>%  # supprime le suffixe _1, _2, etc.
    group_by(Sample, BaseName) %>%
    summarise(Score = sum(Abundance), .groups = "drop")
  
  p_pur_19 <-ggplot(score_grouped_XBP1_RIDD_pur_19 , aes(x = Sample, y = Score, fill = BaseName)) +
    geom_bar(stat = "identity") +
    labs(
      title = "Scores d'abondance par pathway",
      x = "Échantillon",
      y = "Score d'abondance"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
  
  output$score_ire1actif <- renderPlotly({
    ggplotly(p_pur_19 )
  })
  
  score_grouped_XBP1_RIDD_pur_19
  ##OCCURENCE 
  score_grouped_XBP1_RIDD_pur_19$Regulation <- sub(".*_(UP|DOWN)", "\\1", score_grouped_XBP1_RIDD_pur_19$BaseName)
  score_grouped_XBP1_RIDD_pur_19$Gene <- sub("(XBP1|RIDD|XBP1_pur|RIDD_pur).*", "\\1", score_grouped_XBP1_RIDD_pur_19$BaseName)
  
  
  # Créer un tableau de contingence avec la somme des scores
  contingency_table <- table(score_grouped_XBP1_RIDD_pur_19$Regulation,score_grouped_XBP1_RIDD_pur_19$Gene)
  
  
  # Afficher le tableau de contingence
  output$contingence_occurence_actif <- renderUI({
    HTML(
      contingency_table %>%
        kable(format = "html",
              escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "white", bold = TRUE) %>%
        as.character()
    )
  })# Now column 1 is the Direction column
  ##SCORE
  
  
  # Créer un tableau de contingence avec la somme des scores
  formatted_table <- dcast(as.data.table(score_grouped_XBP1_RIDD_pur_19), Regulation ~ Gene, value.var = "Score", fun.aggregate = sum)
  
  # Afficher le tableau formaté
  output$contingence_score_actif <- renderUI({
    HTML(
      formatted_table  %>%
        kable(format = "html",
              escape = FALSE) %>%
        kable_styling(bootstrap_options = c("condensed"),
                      full_width = TRUE,
                      position = "center") %>%
        row_spec(0, color = "white", bold = TRUE) %>%
        as.character()
    )
  })
  
  
  
  
}