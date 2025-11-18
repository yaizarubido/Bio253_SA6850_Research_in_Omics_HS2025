####HEATMAP + Dendogramm

####Spalten orden so das alle PSN zusammen und dann alle TBS zusammen!
#####Spalten nach Klone ordnen
##darstellen als heatmap und Dendogram uf de rows???? J.Grossman
##STX --> den log 2 nehmen
##STX nicht im log2
####delta vom delta nehmen comparing within clones!


######################################################
### HEATMAP ohne Dendrogram --> sortiert nach Clone###
######################################################
library(tidyverse)
library(pheatmap)

### 1) Clone-Reihenfolge bestimmen
clone_order <- df %>%
    distinct(Clones) %>%      # alle Clone einmal
    pull(Clones)

### 2) Samples in Clone-Reihenfolge sortieren
sample_order <- df %>%
    mutate(Clones = factor(Clones, levels = clone_order)) %>%
    arrange(Clones) %>%       # nach Clone sortieren
    pull(Sample)

### 3) Matrix nach Clone-Reihenfolge sortieren
matrix_prot <- matrix_ready[sample_order, ]

### 4) HEATMAP ohne Dendrogram zeichnen
pheatmap(
    matrix_prot,
    scale = "none",
    cluster_rows = FALSE,      # ❌ kein Row-Dendrogram
    cluster_cols = FALSE,      # ❌ kein Col-Dendrogram
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = "Proteomics STX operon (crtM – crtO), sorted by Clone",
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    fontsize = 12,
    border_color = "grey90"
)

##################################################################
###Heatmap ohen STX, ohne dendogram --> SORTIERT NACH TSB/PASN####
##################################################################
library(tidyverse)
library(pheatmap)

### Matrix_ready verwenden (crtM–crtO)
# rows = Samples
# cols = genes crtM–crtO
# values = log2FC

### 2) Samples sortieren (PASN/PSAN zuerst, dann TSB)
sample_info <- df %>%
    distinct(Sample, Condition)

PSN_samples <- sample_info %>% filter(Condition %in% c("PASN","PSAN")) %>% pull(Sample)
TSB_samples <- sample_info %>% filter(Condition == "TSB") %>% pull(Sample)

order_samples <- c(PSN_samples, TSB_samples)

matrix_prot <- matrix_ready[order_samples, ]

### 3) HEATMAP ohne Dendrogram zeichnen
pheatmap(
    matrix_prot,
    scale = "none",                # log2FC unverändert anzeigen
    cluster_rows = FALSE,          # ❌ kein Row-Dendrogram
    cluster_cols = FALSE,          # ❌ kein Col-Dendrogram
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = "Proteomics STX operon (crtM – crtO) log2FC",
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    fontsize = 12,
    border_color = "grey90"
)

#####################################################################
### Heatmap + Dendrogram  FEEDBACK JONAS GROSSMAN   ###
#####################################################################
library(tidyverse)
library(pheatmap)

genes <- c("crtM", "crtN", "crtP", "crtQ", "crtO")

### 1) Matrix vorbereiten
matrix_op <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    select(Sample, all_of(genes)) %>%
    column_to_rownames("Sample")

### 2) Gruppen trennen
samples_6850 <- rownames(matrix_op)[str_detect(rownames(matrix_op), "^6850")]
samples_JE   <- rownames(matrix_op)[str_detect(rownames(matrix_op), "^JE")]

### 3) Optionales Clustering pro Gruppe
cluster_order <- function(mat) {
    hc <- hclust(dist(mat))
    rownames(mat)[hc$order]
}

order_6850 <- cluster_order(matrix_op[samples_6850, ])
order_JE   <- cluster_order(matrix_op[samples_JE, ])

### 4) Leerzeile erzeugen (mit richtigen Namen!)
empty_row <- as.data.frame(matrix(rep(NA, length(genes)), nrow = 1))
colnames(empty_row) <- genes
rownames(empty_row) <- " "

### 5) Endmatrix
matrix_plot <- rbind(
    matrix_op[order_6850, ],
    empty_row,
    matrix_op[order_JE, ]
)

### 6) Heatmap
pheatmap(
    matrix_plot,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = "none",
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    fontsize = 12,
    main = "Proteomics STX operon (crtM–crtO)\n6850 vs JE strains",
    border_color = NA,
    na_col = "white",     # Leerzeile sichtbar
    treeheight_row = 0,
    treeheight_col = 0
)

###Noch den Graph machen Vergelich pairwise TSB vs. PASN Vergleich within clones!
