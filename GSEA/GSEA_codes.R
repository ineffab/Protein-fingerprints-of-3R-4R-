#All required packages 
library(GOSemSim)
library(AnnotationHub)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggarchery)
library(enrichplot)
library(stats)
library(DOSE)
library(ReactomePA)
library(tidyr)
library(dplyr)
library(tibble)
library(pheatmap)
library(ComplexHeatmap)
library(ggtree)
library(stringr)
library(ggplot2)
library(ggrepel)
library(scales)
library(cowplot)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(corrplot)
library(RColorBrewer)

#setup the working directory that has all your input files and would have your outputs stored
setwd("E:/AAT/OneDrive/1 - ONGOING PUBLICATION/EN PREPARATION/4- Espourteille et al-2023/WGCNA/PSP_PICK_WGCNA")

########### PREPROCESSING ###################################################

#load the datafile contaning the proteomic data with sample details
data_mito_PSP = read.csv("PiD_fold_change_GSEA.csv")

#load background list (usually all the proteins names in your data)
#background = read.csv("Background_list.csv")


#Converting the gene symbol ids to entrez ids
entrez_mapping <- bitr(data_mito_PSP$Proteins,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

# Merge the Entrez IDs back into your dataframe
df_merged <- merge(data_mito_PSP, entrez_mapping, by.x = "Proteins", by.y = "SYMBOL", all.x = TRUE)

## assume 1st column is ID
## 2nd column is FC (Fold change)

## feature 1: numeric vector: FC
geneList = df_merged$Fold_change

## feature 2: named vector: entrez ID
names(geneList) = as.charcater(df_merged$ENTREZID)

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

#remove NA's 
geneList <- geneList[!is.na(names(geneList))]


########### GO ENRICHMENT ###################################################
#GO classification first; initiate with BP  
#level 3 is better than level 2 because we are not getting the specific terms

#Biological Pathways
ego_bp <- gseGO(geneList   = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                minGSSize    = 2,
                maxGSSize    = 100,
                pvalueCutoff = 0.05,
                verbose      = FALSE)

#write.csv(ego_bp@result, "GSEA_mito_BP_all.csv")

#Molecular Function
ego_mf <- gseGO(geneList   = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "MF",
                minGSSize    = 2,
                maxGSSize    = 100,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
#write.csv(ego_mf@result, "GSEA_mito_MF_all.csv")

#Cellular Components
ego_cc <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "CC",
                minGSSize    = 2,
                maxGSSize    = 100,
                pvalueCutoff = 0.05,
                verbose      = FALSE)
#write.csv(ego_cc@result, "GSEA_mito_CC_all.csv")

#KEGG
ego_kk <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  minGSSize    = 2,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)
#write.csv(ego_kk@result, "GSEA_mito_KEGG_all.csv")

#Wiki Pathways
ego_wp = gseWP(geneList, organism = "Homo sapiens")
#write.csv(ego_wp@result, "GSEA_mito_Wiki_all.csv")

#Reactome
ego_react = gsePathway(geneList, pvalueCutoff = 0.2, pAdjustMethod = "BH", verbose=FALSE)

########### GROUP ANNOTATION ###################################################

# Pick top N from each
extract_top_gsea <- function(gsea_obj, db_name, top_n) {
  gsea_obj@result %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::select(ID, Description, NES, p.adjust, core_enrichment) %>%
    dplyr::mutate(Database = db_name)
}


df_go_bp     <- extract_top_gsea(ego_bp, "GOBP", 50)
df_go_mf     <- extract_top_gsea(ego_mf, "GOMF", 50)
df_go_cc     <- extract_top_gsea(ego_cc, "GOCC", 50)
df_kegg   <- extract_top_gsea(ego_kk, "KEGG",  50)
df_react  <- extract_top_gsea(ego_react, "Reactome", 50)
df_wp     <- extract_top_gsea(ego_wp, "WP", 50)

combined_df <- bind_rows(df_go_bp,df_go_mf,df_go_cc, df_kegg, df_react, df_wp)
write.csv(combined_df,"Enrich_50_PiD_ombine_all.csv")

# Create a lowercase version of the term description
df <- combined_df %>%
  mutate(Description_lc = tolower(Description))

#I tried to assign the name of module manually here 
# Define a function to assign modules
assign_module <- function(desc) {
  if (str_detect(desc, "synapse|synaptic|neuro|glutamate|gaba|presynaptic|exocytosis|recycling")) {
    return("Synaptic regulation")
  } else if (str_detect(desc, "immune|b cell|t cell|immunoglobulin|inflammatory|cytokine|neutrophil|complement")) {
    return("Immune Response and Regulation")
  } else if (str_detect(desc, "peptidase|hydrolase|endopeptidase|enzyme inhibitor|protease")) {
    return("Peptidase activity")
  } else if (str_detect(desc, "extracellular matrix|ecm|integrin|collagen|adhesion|basement membrane")) {
    return("Extracellular matrix and adhesion")
  } else if (str_detect(desc, "mitochondria|oxidoreductase|tryptophan|metabolism|vitamin d")) {
    return("Mitochondria and metabolism")
  } else if (str_detect(desc, "endoplasmic reticulum|er lumen|golgi|secretory|retrograde")) {
    return("Endoplasmic reticulum and secretory pathway")
  } else {
    return("Unclassified")
  }
}

# Apply the function to assign modules
df$Functional_Module <- sapply(df$Description_lc, assign_module)

write.csv(df, "Assign_PiD_50.csv")

########### MANUAL ANNOTATING ###################################################
#As not all the terms are described well by the keywords mentioned above, 
#it is highly recommended to go through the file manually
#and annotate the suitable and reasonable term with the categories used before. 
#And the manually annotated file would be used further


########### VISUALIZATION ###################################################

#Manually_assigned_unclassified terms
Assigned_data <- read.csv("Assign_PiD_50.csv")

# Optional: filter to remove "Unclassified"
df_1 <- Assigned_data %>% filter(Functional_Module != "Unclassified")

# STEP 1: Process data
agg_df <- df_1 %>%
  mutate(
    GeneSetSize = stringr::str_count(core_enrichment, "/") + 1
  ) %>%
  group_by(Functional_Module, Database) %>%
  summarise(
    Avg_NES = mean(NES, na.rm = TRUE),
    Avg_FDR = mean(p.adjust, na.rm = TRUE),
    Avg_GeneSetSize = mean(GeneSetSize),
    .groups = "drop"
  ) %>%
  mutate(
    NegLog10_FDR = -log10(Avg_FDR),
    ScaledGeneSetSize = log1p(Avg_GeneSetSize)  # visual scaling
  )

# STEP 2: Fix consistent order for factors
# Replace with your desired fixed order
database_order <- c("GOBP","GOCC","GOMF","KEGG","WP", "Reactome")  # customize this
module_order <- rev(sort(unique(agg_df$Functional_Module)))  # or define manually

agg_df$Database <- factor(agg_df$Database, levels = database_order)
agg_df$Functional_Module <- factor(agg_df$Functional_Module, levels = module_order)

# Define common scales manually outside the plotting function

# Define gene set size legend: use full global range you want
gene_legend_values <- c(5, 10, 25, 50, 100)
gene_legend_scaled <- log1p(gene_legend_values)

# STEP 4: Final bubble plot
ggplot(agg_df, aes(x = Database, y = Functional_Module)) +
  geom_point(
    aes(size = ScaledGeneSetSize, fill = Avg_NES, alpha = NegLog10_FDR),
    shape = 21, color = "black", stroke = 1.1
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Average NES"
  ) +
  scale_size_continuous(
    name = "Gene Set Size",
    breaks = gene_legend_scaled,
    labels = gene_legend_values,
    range = c(5, 17),
  ) +
  scale_alpha(
    name = "-log10(FDR)",
    range = c(0.4, 1),
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_line(),
    axis.text.x = element_text(colour = "black", face = "bold", angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(colour = "black", face = "bold", size = 15),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    legend.position = "right",
    legend.box = "vertical",
    legend.justification = "center",
    legend.direction = "vertical"
  ) +
  labs(
    title = "Gene Set Enrichment Analysis of BD-EV Proteins in PiD"
  )

#Put dot plot to show mitochondrial and apototic terms. 

specific <- filter(df_1, df_1$Functional_Module == "Endoplasmic reticulum and secretory pathway")

# Prepare plot data
specific <- specific %>%
  mutate(
    NegLog10FDR = -log10(p.adjust),
    Term = factor(Description, levels = rev(Description))  # To control order in plot
  )%>%
  mutate(
    GeneSetSize = stringr::str_count(core_enrichment, "/") + 1
  )

# Split into two dataframes
#df_mito <- specific %>% filter(Functional_Module == "Mitochondria and metabolism")
#df_apop <- specific %>% filter(Functional_Module == "Apoptosis")


#Dot Plot
ggplot(specific, aes(x = NES, y = Term)) +
  geom_point(aes(size = GeneSetSize, color = NegLog10FDR)) +
  scale_color_gradient(low = "blue" ,high = "red", name = "-log10(FDR)") +
  scale_size(range = c(4, 10), name = "Gene Set Size") +
  #facet_wrap(~ Functional_Module, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", size = 14, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    title = "Enriched Endoplasmic reticulum Terms"
  )


# Bar plot
ggplot(df_mito, aes(x = NES, y = Term, fill = NegLog10FDR)) +
  geom_col(width = 0.6) +
  scale_fill_gradient(low = "midnightblue", high = "red", name = "-log10(FDR)") +
  #facet_wrap(~ Functional_Module, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  ) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    title = "Top Mitochondrial Terms in PSP"
  )


###facets by database

ggplot(agg_df, aes(x = Functional_Module, y = Avg_NES)) +
  geom_point(
    aes(size = ScaledGeneSetSize, fill = Avg_NES, alpha = NegLog10_FDR),
    shape = 21, color = "black", stroke = 1.1
  ) +
  facet_wrap(~ Database, scales = "free_y") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Average NES"
  ) +
  scale_size_continuous(
    name = "Gene Set Size",
    breaks = gene_legend_scaled,
    labels = gene_legend_values,
    range = c(5, 17)
  ) +
  scale_alpha(range = c(0.4, 1), name = "-log10(FDR)") +
  coord_flip() +  # Optional: makes modules on y-axis, easier to read
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_line(),
    axis.text.x = element_text(colour = "black", face = "bold", size = 12),
    axis.text.y = element_text(colour = "black", face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5)
  ) +
  labs(
    title = "GSEA of BD-EV Proteins in PSP by Database"
  )


########### EXTRACTION ###################################################

# Function to extract top geneSetID
get_top_term_info <- function(gsea_obj, by = "p.adjust") {
  df <- gsea_obj@result
  if (nrow(df) == 0) return(NULL)
  top <- df %>% arrange(.data[[by]]) %>% slice(1)
  return(list(ID = top$ID[1], Description = top$Description[1]))
}


# # Function to generate GSEA plot with title
# get_gsea_plot <- function(gsea_obj, term_info, db_label) {
#   plot_out <- gseaplot2(
#     gsea_obj,
#     geneSetID = term_info$ID,
#     title = paste0(db_label, ": ", term_info$Description)
#   )
#   if (inherits(plot_out, "gg")) {
#     return(plot_out)
#    } else if (is.list(plot_out)) {
#     return(plot_out[[1]])
#   } else {
#     stop("Unknown gseaplot2 output format")
#   }
# }

get_gsea_plot <- function(gsea_obj, term_info, db_label) {
  plots <- gseaplot2(
    gsea_obj,
    geneSetID = term_info$ID,
    subplots = 1:2,
    title = paste0(db_label, ": ", term_info$Description)
  )
  
  if (is.list(plots)) {
    # Combine enrichment score and ranked metric into one full GSEA plot
    return(plots[[1]] / plots[[2]] + patchwork::plot_layout(heights = c(2, 1)))
  } else {
    return(plots)
  }
}
 

# Extract top terms
top_go_bp    <- get_top_term_info(ego_bp, "p.adjust")
top_go_mf<- get_top_term_info(ego_mf, "p.adjust")
top_go_cc    <- get_top_term_info(ego_cc, "p.adjust")
top_kegg  <- get_top_term_info(ego_kk, "p.adjust")
top_react <- get_top_term_info(ego_react, "p.adjust")
top_wp    <- get_top_term_info(ego_wp, "p.adjust")

# Create plots safely
p_go_bp     <- get_gsea_plot(ego_bp, top_go_bp, "GOBP")
p_go_mf     <- get_gsea_plot(ego_mf, top_go_mf, "GOMF")
p_go_cc     <- get_gsea_plot(ego_cc, top_go_cc, "GOCC")
p_kegg   <- get_gsea_plot(ego_kk, top_kegg, "KEGG")
p_react  <- get_gsea_plot(ego_react, top_react, "Reactome")
p_wp     <- get_gsea_plot(ego_wp, top_wp, "WP")

(p_go_bp | p_go_mf) / (p_go_cc | p_kegg) / (p_react | p_wp) + plot_layout(guides = "collect")


# # Combine with cowplot
# cowplot::plot_grid(
#   cowplot::ggdraw(p_go_bp),
#   cowplot::ggdraw(p_go_mf),
#   cowplot::ggdraw(p_go_cc),
#   cowplot::ggdraw(p_kegg),
#   cowplot::ggdraw(p_react),
#   cowplot::ggdraw(p_wp),
#   ncol = 3,
#   labels = LETTERS[1:6]
# )