#### Omics BIO253 - Proteomic STX
#### Matrix MIT EINZELNEN PFEILEN (diverser Länge) + STX-AVERAGE
#############################################
###TODO - POSSIBLE IMPROVEMENTS:
######add a black vertical line separating the STX-Average values (or even more spacing)
######add/improve axes + add evtl. a titel
library(tidyverse)
library(readxl)
library(writexl)

#### Daten einlesen
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
file <- "resources/LaraYaiza_tabel_STX.xlsx"          # relative folder path!

#### richtige Tabelle einlesen
df <- read_excel(file, sheet = "Tabelle2", skip=3)

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

#### Ausgabe ansehen --> MATRIX MIT ZAHLEN
matrix_ready

############################################################################
### Erweiterung: Vertikale Pfeil-Matrix für log2FC-Werte  -->ALPHABETICAL###
############################################################################

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
#########################################################
### ADDITION STX-Average als letzte Spalte hinzufügen ###
##########################################################

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
    select(Sample, STX = `STX Average`) %>%
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

##################################################################
#########################################################
### ADDITION — Clean STX column (clone-level) + spacing
##########################################################

### 1) Long format for arrows (same as before)
df_long <- matrix_ready %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Gene", values_to = "log2FC") %>%
  left_join(df %>% select(Sample, Clones), by = "Sample") %>%
  mutate(
    Sample_f = factor(Sample, levels = rev(unique(Sample))),
    Gene_f = factor(Gene, levels = c("crtM","crtN","crtP","crtQ","crtO")),
    sign = ifelse(is.na(log2FC), "", ifelse(log2FC > 0, "↑", "↓")),
    mag  = ifelse(is.na(log2FC), 0, scales::rescale(abs(log2FC), to = c(4,10)))
  )

### 2) Build clone-level STX table
df_stx <- df %>%
  distinct(Clones, STX = `STX Average`) %>%    # only one row per clone
  left_join(df %>% select(Sample, Clones), by = "Clones") %>%  # match Sample rows
  mutate(
    Gene_f = "STX",
    Sample_f = factor(Sample, levels = rev(unique(df_long$Sample))),
    sign = "",         # no arrows here
    log2FC = NA,
    mag = 0
  )

### 3) Combine arrow data + STX column
df_plot <- bind_rows(
  df_long %>% mutate(Gene_f = factor(Gene_f, levels = c("crtM","crtN","crtP","crtQ","crtO","STX"))),
  df_stx  %>% mutate(Gene_f = factor(Gene_f, levels = c("crtM","crtN","crtP","crtQ","crtO","STX")))
)

### 4) Plot with extra spacing to STX
ggplot(df_plot, aes(Gene_f, Sample_f)) +
  geom_tile(fill = "grey97", color = "white") +

  # arrows for genes
  geom_text(
    data = df_plot %>% filter(Gene_f != "STX"),
    aes(label = sign, color = log2FC, size = mag),
    fontface = "bold"
  ) +

  # clone-level STX numbers
  geom_text(
    data = df_plot %>% filter(Gene_f == "STX"),
    aes(label = round(STX, 2)),
    color = "black",
    fontface = "bold",
    size = 4
  ) +

  scale_color_gradient2(
    low = "#2166ac", high = "#b2182b", mid = "grey80",
    midpoint = 0, name = "log2FC"
  ) +
  scale_size_identity() +

  # add spacing before STX
  scale_x_discrete(
    expand = expansion(add = c(0.1, 1.4))   # <-- makes STX far right separated
  ) +

  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(10, 40, 10, 10)   # extra right margin
  )

#########################################################
### FINAL: STX nur einmal pro Klon + grosser Abstand
##########################################################

# 1) Long format for arrow matrix
df_long <- matrix_ready %>%
    rownames_to_column("Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "log2FC") %>%
    left_join(df %>% select(Sample, Clones, Condition), by = "Sample") %>%
    mutate(
        Sample_f = factor(Sample, levels = rev(unique(Sample))),
        Gene_f = factor(Gene, levels = c("crtM","crtN","crtP","crtQ","crtO")),
        sign = ifelse(is.na(log2FC), "", ifelse(log2FC > 0, "↑", "↓")),
        mag  = ifelse(is.na(log2FC), 0, scales::rescale(abs(log2FC), to = c(4,10)))
    )

# 2) STX-One-Per-Clone: keep only STX on TSB row
df_stx <- df %>%
    select(Sample, Clones, Condition, STX = `STX Average`) %>%
    mutate(
        keep_STX = ifelse(Condition == "TSB", TRUE, FALSE)
    ) %>%
    mutate(
        STX_show = ifelse(keep_STX, STX, NA),
        Gene_f = "STX",
        Sample_f = factor(Sample, levels = rev(unique(df_long$Sample))),
        sign = "",
        log2FC = NA,
        mag = 0
    ) %>%
    select(Sample, Sample_f, Gene_f, STX_show, sign, log2FC, mag)

# 3) Combine gene arrows + STX
df_plot <- bind_rows(
    df_long %>% mutate(Gene_f = factor(Gene_f,
                                       levels = c("crtM","crtN","crtP","crtQ","crtO","STX"))),
    df_stx  %>% mutate(Gene_f = factor(Gene_f,
                                       levels = c("crtM","crtN","crtP","crtQ","crtO","STX")))
)

# 4) Plot with strong spacing before STX
ggplot(df_plot, aes(Gene_f, Sample_f)) +
    geom_tile(fill = "grey97", color = "white") +

    # gene arrows
    geom_text(
        data = df_plot %>% filter(Gene_f != "STX"),
        aes(label = sign, color = log2FC, size = mag),
        fontface = "bold"
    ) +

    # STX only once per clone (TSB row)
    geom_text(
        data = df_plot %>% filter(Gene_f == "STX"),
        aes(label = round(STX_show, 2)),
        color = "black",
        fontface = "bold",
        size = 5
    ) +

    scale_color_gradient2(
        low = "#2166ac", high = "#b2182b", mid = "grey80",
        midpoint = 0, name = "log2FC"
    ) +
    scale_size_identity() +

    # BIG spacing before STX
    scale_x_discrete(
        expand = expansion(add = c(0.1, 3.0))   # <<--- much larger gap after crtO
    ) +

    theme_minimal(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(10, 60, 10, 10)   # extra margin right
    )
