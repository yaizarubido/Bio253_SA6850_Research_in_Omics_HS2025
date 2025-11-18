###BIO253-OMICS - Proteomix-STX-Production in Staphylococus aureus
###
###


library(tidyverse)
library(readxl)
library(pheatmap)
library(grid)
library(gridExtra)

###Daten einlesen
pfad <- "C:/Users/rubid/Documents/1_Prio_Biomedizin/AAA/00_AAA_BIO253/Bio253_SA6850_Research_in_Omics_HS2025-main"
datei <- "tabel_STX.xlsx"
full_path <- file.path(pfad, datei)

df <- read_excel(full_path, sheet = "Tabelle2", skip = 3)

###Spaltennamen überprüfen:
#names(df)

#Missing values in numeric values umwandlen! --> Missin gmit NA ersetzen!
df[df == "Missing"] <- NA
numeric_cols <- c("crtM", "crtN", "crtP", "crtQ", "crtO",
                  "STX Average", "Growth rate", "H2O2 survival",
                  "sigma-B", "rsbU", "rsbV", "rsbW")

df <- df %>%
    mutate(across(all_of(numeric_cols), ~as.numeric(.)))

##data strucutre check
#str(df)
#head(df)

#### delta = PASN-TSB berechenen für STX-Production & STX-Regulation PrtX
delta_df <- df %>%
    pivot_wider(
        names_from = Condition,
        values_from = c(crtM, crtN, crtP, crtQ, crtO,
                        `sigma-B`, rsbU, rsbV, rsbW,
                        `Growth rate`)
    ) %>%
    mutate(
        delta_crtM = crtM_PASN - crtM_TSB,
        delta_crtN = crtN_PASN - crtN_TSB,
        delta_crtP = crtP_PASN - crtP_TSB,
        delta_crtQ = crtQ_PASN - crtQ_TSB,
        delta_crtO = crtO_PASN - crtO_TSB,

        delta_sigmaB = `sigma-B_PASN` - `sigma-B_TSB`,
        delta_rsbU   = rsbU_PASN - rsbU_TSB,
        delta_rsbV   = rsbV_PASN - rsbV_TSB,
        delta_rsbW   = rsbW_PASN - rsbW_TSB,

        delta_growth = `Growth rate_PASN` - `Growth rate_TSB`
    )

### EINDEUTIGE Zeilen-ID erzeugen: Strain_Clone
delta_df <- delta_df %>%
    mutate(Clone_ID = paste(Strain, Clones, sep = "_"))

###Heatmap-Matrix erstellen
heat_data <- delta_df %>%
    select(Clone_ID, starts_with("delta_")) %>%
    column_to_rownames("Clone_ID")

###Spalten mit NA entfernen --> sonst HeatMap verfälscht!????
heat_data_clean <- heat_data %>%
    select(where(~!any(is.na(.))))
##CHECK WHICH GOT "DELETED"
#print("Entfernte Spalten ween NA:")
#print(setdiff(colnames(heat_data),colnames(heat_data_clean)))
###crtQ, crtO, rsbU --> keine delta berechenbar!

###Heatmap erzeugen
pheatmap(
    heat_data_clean,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 11,
    main = "Δ Proteomics + Growth rate (PASN – TSB)"
)


################################################################
###WITHOUT GROWTHRATE - ONLY GENES###
################################################################
heat_data_clean_noGrowth <- heat_data_clean %>%
    select(-delta_growth)

pheatmap(
    heat_data_clean_noGrowth,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 11,
    main = "Δ Proteomics (PASN – TSB, ohne Growth rate)"
)
###MACHT DAS SINN REGULATORY UND PRODUCTION IN EIN MAP???
###ODER SEPARATE PRODUCTION & REGULATORY MAP?

##############################################################
###ONLY STX_PRODUCTION GENES###
##############################################################
heat_stx <- heat_data_clean %>%
    select(delta_crtM, delta_crtN, delta_crtP)

pheatmap(
    heat_stx,
    scale = "row",
    main = "Δ STX-Biosynthesis Proteins"
)

############################################################
###ONLY STX-REGULATORY GENES###
###########################################################
###sigmaB, rsbU, rsbV, rsbW

heat_reg <- heat_data %>%
    select(delta_sigmaB, delta_rsbU, delta_rsbV, delta_rsbW)

# NA-freie Spalten sicherstellen
heat_reg <- heat_reg %>%
    select(where(~ !any(is.na(.))))

pheatmap(
    heat_reg,
    scale = "row",
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    fontsize = 11,
    main = "Δ STX-Regulation (σB & rsb genes)"
)
