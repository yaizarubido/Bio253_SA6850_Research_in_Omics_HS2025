#### Omics BIO253 - Proteomic STX
#### Matrix with numbers and MULTIPLE ARROW Matrix
###

library(tidyverse)
library(readxl)
library(writexl)

#### Daten einlesen
setwd("../Bio253_SA6850_Research_in_Omics_HS2025")
#path <- "..Bio253_SA6850_Reseach_in_Omics_HS2025/resources"        # relative folder path!
file <- "resources/LaraYaiza_tabel_STX.xlsx"                        # Dateiname Resource
# full_path <- file.path(path, file)
# print(full_path)
# print(file.exists(full_path))

#### richtige Tabelle einlesen
df <- read_excel(file,sheet = "Tabelle2",skip=3)

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

#### Ausgabe ansehen mit ZAHLEN
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

###################################################################
### Erweiterung: Pfeil-Matrix mit sortierter Sample-Reihenfolge ###
### (oben TSB, unten PASN)                                      ###
###################################################################

# 1) Condition-Information extrahieren
sample_info <- df %>%
    select(Sample, Condition)

# 2) Reihenfolge definieren: TSB zuerst, dann PASN / PSAN
sample_order <- sample_info %>%
    mutate(order_group = case_when(
        Condition == "TSB" ~ 1,
        Condition %in% c("PASN", "PSAN") ~ 2,
        TRUE ~ 3
    )) %>%
    arrange(order_group, Sample) %>%
    pull(Sample)

# 3) arrow_matrix in dieser Reihenfolge neu sortieren
arrow_matrix_sorted <- arrow_matrix %>%
    slice(match(sample_order, rownames(arrow_matrix)))

# 4) Ausgabe zeigen
arrow_matrix_sorted
