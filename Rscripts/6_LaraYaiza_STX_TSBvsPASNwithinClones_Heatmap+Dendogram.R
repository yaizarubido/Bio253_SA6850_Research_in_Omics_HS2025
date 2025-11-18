### OMICS-2025 — Δ Heatmaps (TSB–PASN)

library(tidyverse)
library(readxl)
library(pheatmap)
library(grid)
library(gridExtra)

#### Daten einlesen
pfad <- "C:/Users/rubid/Documents/1_Prio_Biomedizin/AAA/00_AAA_BIO253/Bio253_SA6850_Research_in_Omics_HS2025-main"
datei <- "tabel_STX.xlsx"
full_path <- file.path(pfad, datei)

df <- read_excel(full_path, sheet = "Tabelle2", skip = 3)

genes <- c("crtM","crtN","crtP","crtQ","crtO")

### 1) Missing → NA → numeric
df <- df %>%
    mutate(
        across(all_of(genes), as.character),
        across(all_of(genes), ~na_if(.,"Missing")),
        across(all_of(genes), as.numeric)
    )

### 2) Delta-Funktion
compute_delta <- function(df, strain_name) {

    df_strain <- df %>%
        filter(Strain == strain_name) %>%
        select(Strain, Clones, Condition, all_of(genes)) %>%
        pivot_longer(all_of(genes),
                     names_to="Gene",
                     values_to="Value") %>%
        pivot_wider(names_from=Condition, values_from=Value) %>%
        mutate(Delta = TSB - PASN)

    mat <- df_strain %>%
        select(Clones, Gene, Delta) %>%
        pivot_wider(names_from = Gene, values_from = Delta) %>%
        column_to_rownames("Clones") %>%
        as.matrix()

    # fürs Clustering → NA ersetzen
    mat_clust <- mat
    mat_clust[is.na(mat_clust)] <- 0

    list(plot = mat, clust = mat_clust)
}

### 3) Delta-Matrizen
m6850 <- compute_delta(df, "6850")
mJE2  <- compute_delta(df, "JE2")

### 4) pheatmap → gtable extrahieren
extract_grob <- function(mat_list, title) {
    ph <- pheatmap(
        mat_list$plot,
        clustering_distance_rows = dist(mat_list$clust),
        clustering_distance_cols = dist(t(mat_list$clust)),
        color = colorRampPalette(c("#2166ac","white","#b2182b"))(100),
        na_col = "black",
        main = title,
        silent = TRUE
    )

    # pheatmap-Objekt enthält $gtable
    return(ph$gtable)
}

g1 <- extract_grob(m6850, "Δ (TSB – PASN) — Strain 6850")
g2 <- extract_grob(mJE2,  "Δ (TSB – PASN) — Strain JE2")

### 5) Heatmaps untereinander
gridExtra::grid.arrange(g1, g2, ncol = 1)
