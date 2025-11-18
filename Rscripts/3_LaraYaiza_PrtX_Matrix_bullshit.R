#### Omics BIO253
#### Matrix Proteomix STX_protein expression
### Arrow Matrix with STX-Average

library(tidyverse)
library(readxl)
library(writexl)

#### Daten einlesen
pfad <- "C:/Users/rubid/Documents/1_Prio_Biomedizin/AAA/00_AAA_BIO253/Bio253_SA6850_Research_in_Omics_HS2025-main"
datei <- "LaraYaiza_tabel_STX.xlsx"

full_path <- file.path(pfad, datei)

#### richtige Tabelle einlesen
df <- read_excel(full_path,sheet = "Tabelle2", skip=3)

#### Gen-Spalten definieren
genes <- c("crtM", "crtN", "crtP", "crtQ", "crtO")

#### Sample-Namen generieren (Strain + Clone + Condition)
df <- df %>%
    mutate(Sample = paste(Strain, Clones, Condition, sep = "-"))

#### Missing → NA umwandeln & numerisch machen
df <- df %>%
    mutate(across(all_of(genes), as.character)) %>%   # alles erstmal als text
    mutate(across(all_of(genes), ~ na_if(., "Missing"))) %>%  # Missing -> NA
    mutate(across(all_of(genes), as.numeric))          # zurück zu numeric


#### Matrix mit ZAHLEN direkt erstellen (Samples als Zeilen, Gene als Spalten)
matrix_ready <- df %>%
    select(Sample, all_of(genes)) %>%
    column_to_rownames("Sample")

#### Ausgabe ansehen
matrix_ready

### Matrix mit MEHREREN PFEILEN
arrowify <- function(x, factor = 2) {
    sapply(x, function(v) {
        if (is.na(v)) return("")
        n <- round(abs(v) * factor)
        if (n == 0) return("-")
        if (v > 0) return(paste(rep("↑", n), collapse = ""))
        if (v < 0) return(paste(rep("↓", n), collapse = ""))
    })
}

arrow_matrix <- matrix_ready %>%
    mutate(across(everything(), ~ arrowify(.)))

arrow_matrix

###==============================================================

#############################################################
### Erweiterung: Vertikale Pfeil-Matrix für log2FC-Werte  ###
#############################################################

library(tidyverse)
library(ggplot2)

### 1) Matrix ins lange Format bringen
df_long <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    pivot_longer(
        cols = -Sample,
        names_to = "Gene",
        values_to = "log2FC"
    )

### 2) Pfeilmatrix-Plot (genau EIN vertikaler Pfeil pro Zelle)
pfeil_plot <- ggplot(df_long, aes(x = Gene, y = 0)) +
    geom_segment(
        aes(x = Gene,
            xend = Gene,
            y = 0,
            yend = log2FC,
            size = abs(log2FC),
            color = log2FC),
        arrow = arrow(length = unit(0.25, "cm")),
        lineend = "round",
        na.rm = TRUE
    ) +
    facet_grid(
        rows = vars(Sample),
        switch = "y"
    ) +
    scale_color_gradient2(
        low = "blue",      # negative FC (down)
        mid = "grey80",
        high = "red",      # positive FC (up)
        midpoint = 0
    ) +
    scale_size(
        range = c(0.5, 2.5),
        guide = "none"
    ) +
    geom_hline(
        yintercept = 0,
        color = "black",
        linewidth = 0.3
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 14) +
    theme(
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        panel.spacing.y = unit(1, "lines")
    )

### Plot anzeigen
pfeil_plot

#########==========================================
###NEUER PFEILPLOT - BESSER????

library(ggplot2)
library(dplyr)
library(tidyr)

### Matrix → long format
df_long <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "log2FC") %>%
    mutate(
        Sample_f = factor(Sample, levels = rev(unique(Sample))),  # ordnen
        Gene_f   = factor(Gene, levels = genes)                   # Gene sortieren
    )

### Arrow-Heatmap
ggplot(df_long, aes(x = Gene_f, y = Sample_f)) +

    # Hintergrund-Kacheln
    geom_tile(fill = "grey95", color = "white") +

    # Pfeile
    geom_segment(
        aes(x = Gene_f,
            xend = Gene_f,
            y = as.numeric(Sample_f),
            yend = as.numeric(Sample_f) + (log2FC * 0.4),
            color = log2FC,
            size = abs(log2FC)),
        arrow = arrow(length = unit(0.25, "cm")),
        lineend = "round",
        na.rm = TRUE
    ) +

    scale_color_gradient2(
        low = "#2166ac",
        mid = "grey80",
        high = "#b2182b",
        midpoint = 0,
        name = "log2FC"
    ) +

    scale_size(range = c(0.5, 1.5), guide = "none") +

    # Theme
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.margin = margin(20, 20, 20, 20)
    )

####====================================================
library(ggplot2)
library(dplyr)
library(tidyr)

df_long <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "log2FC") %>%
    mutate(
        Sample_f = factor(Sample, levels = rev(unique(Sample))),
        Gene_f = factor(Gene, levels = c("crtM", "crtN", "crtP", "crtQ", "crtO")),
        log2FC_clamped = ifelse(abs(log2FC) < 0.1, 0, log2FC)  # Avoid tiny arrows
    )

max_len <- max(abs(df_long$log2FC), na.rm = TRUE)

ggplot(df_long, aes(Gene_f, Sample_f)) +

    # Soft background tiles
    geom_tile(fill = "grey97", color = "white") +

    # Centered arrows with normalized lengths
    geom_segment(
        aes(
            x = as.numeric(Gene_f),
            xend = as.numeric(Gene_f),
            y = as.numeric(Sample_f),
            yend = as.numeric(Sample_f) +
                (log2FC_clamped / max_len) * 0.45,  # normalized length
            color = log2FC,
            size = abs(log2FC_clamped)
        ),
        arrow = arrow(length = unit(0.25, "cm"),
                      ends = "last",
                      type = "closed"),
        lineend = "round",
        na.rm = TRUE
    ) +

    # Colors: symmetrical, clean diverging palette
    scale_color_gradient2(
        low = "#4575b4",   # blue
        mid = "grey80",
        high = "#d73027",  # red
        midpoint = 0,
        name = "log2FC"
    ) +

    # Control arrow thickness
    scale_size(range = c(0.4, 1.8), guide = "none") +

    # Cleaner theme
    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
    )
######===================================================
####PERFECT MATRIX

df_long <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "log2FC") %>%
    mutate(
        Sample_f = factor(Sample, levels = rev(unique(Sample))),
        Gene_f = factor(Gene, levels = c("crtM","crtN","crtP","crtQ","crtO")),
        sign = case_when(
            is.na(log2FC) ~ "",     # no arrow for NA
            log2FC > 0 ~ "↑",
            log2FC < 0 ~ "↓",
            TRUE ~ "•"
        ),
        mag = ifelse(
            is.na(log2FC),
            0,  # NA → invisible
            scales::rescale(abs(log2FC), to = c(4, 10))
        )
    )

ggplot(df_long, aes(Gene_f, Sample_f)) +

    geom_tile(fill = "grey97", color = "white") +

    # Draw only non-NA arrows
    geom_text(
        data = df_long %>% filter(!is.na(log2FC)),
        aes(label = sign, color = log2FC, size = mag),
        fontface = "bold"
    ) +

    scale_color_gradient2(
        low = "#2166ac",
        high = "#b2182b",
        mid = "grey80",
        midpoint = 0,
        name = "log2FC"
    ) +

    scale_size_identity() +

    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.margin = margin(10, 20, 10, 10)
    )
######=====================================================
###add STX-Average als letzte Spalte

### 1) Long format for arrows
df_long <- matrix_ready %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Gene", values_to = "log2FC") %>%
  mutate(
    Sample_f = factor(Sample, levels = rev(unique(Sample))),
    Gene_f = factor(Gene, levels = c("crtM","crtN","crtP","crtQ","crtO")),
    sign = case_when(
      is.na(log2FC) ~ "",
      log2FC > 0 ~ "↑",
      log2FC < 0 ~ "↓",
      TRUE ~ "•"
    ),
    mag = ifelse(
      is.na(log2FC),
      0,
      scales::rescale(abs(log2FC), to = c(4, 10))
    )
  )

### 2) Extract STX values and create a new “STX” column as a pseudo-gene
df_stx <- df %>%
  select(Sample, STX = `Phenotype STX Average`) %>%
  mutate(
    Sample_f = factor(Sample, levels = rev(unique(df_long$Sample))),
    Gene_f = factor("STX", levels = c("crtM","crtN","crtP","crtQ","crtO","STX")),
    sign = "",        # no arrows
    log2FC = NA,
    mag = 0
  )

### 3) Combine arrow data + STX column
df_plot <- bind_rows(
  df_long %>% mutate(Gene_f = factor(Gene_f, levels = c("crtM","crtN","crtP","crtQ","crtO","STX")),
                     STX = NA),
  df_stx
)

### 4) Final plot
ggplot(df_plot, aes(Gene_f, Sample_f)) +

  geom_tile(fill = "grey97", color = "white") +

  # draw arrows
  geom_text(
    data = df_plot %>% filter(!is.na(log2FC)),
    aes(label = sign, color = log2FC, size = mag),
    fontface = "bold"
  ) +

  # draw STX numbers in new column
  geom_text(
    data = df_plot %>% filter(Gene_f == "STX"),
    aes(label = round(STX, 2)),
    color = "black",
    fontface = "bold",
    size = 4
  ) +

  scale_color_gradient2(
    low = "#2166ac",
    high = "#b2182b",
    mid = "grey80",
    midpoint = 0,
    name = "log2FC"
  ) +

  scale_size_identity() +

  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

####Spalten orden so das alle PSN zusammen und dann alle TBS zusammen!
##darstellen als heatmap und dendogram uf de rows???? J.Grossman
##STX --> den log 2 nehmen
####===================================================================
#############################################
### Heatmap + Dendrogram + log2(STX)   FEEDBACK JONAS GROSS   ###
#############################################

library(tidyverse)
library(pheatmap)

###-----------------------------------------
### 1) STX log2-transformieren
###-----------------------------------------

df2 <- df %>%
    select(Sample, STX = `Phenotype STX Average`) %>%
    mutate(STX_log2 = log2(STX))

###-----------------------------------------
### 2) Matrix_ready erweitern mit STX_log2
###-----------------------------------------

matrix2 <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    left_join(df2 %>% select(Sample, STX_log2), by = "Sample") %>%
    column_to_rownames("Sample")

### Jetzt Matrix enthält:
### crtM, crtN, crtP, crtQ, crtO, STX_log2

###-----------------------------------------
### 3) Reihenfolge der Samples: zuerst PSN, dann TSB
###-----------------------------------------

# Condition aus Sample herauslesen
sample_info <- df %>%
    distinct(Sample, Condition)

PSN_samples <- sample_info %>% filter(Condition %in% c("PASN","PSAN")) %>% pull(Sample)
TSB_samples <- sample_info %>% filter(Condition == "TSB") %>% pull(Sample)

order_samples <- c(PSN_samples, TSB_samples)

# Matrix in gewünschte Reihenfolge
matrix2 <- matrix2[order_samples, ]

###-----------------------------------------
### 4) Heatmap mit Row-Dendrogram
###-----------------------------------------

pheatmap(
    matrix2,
    scale = "none",              # KEINE weitere Skalierung
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    clustering_distance_cols = "euclidean",
    main = "Proteomix STX operon + STX phenotype (log2)",
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    fontsize = 12,
    border_color = NA,
    treeheight_row = 60,
    treeheight_col = 40,
    cluster_cols = FALSE  # CRT-Gene in fixer Reihenfolge
)

