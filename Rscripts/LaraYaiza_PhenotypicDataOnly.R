###bio253-omics
###PHENOTYPIC DATA

#What can we look at?
#Correlations betwenn STX and H2O2survival and effect on growth rate difference!
#Does higher staphyloxanthin correlate with ROS resistance?
#Does more staphyloxanthin (antioxidant pigment) lead to better peroxide resistance? --> allready well-establishes in the literature --> can confirm this!
#Do pigment-rich clones handle P. aeruginosa supernatant better (or worse)?

### BIO253-OMICS -- Phenotypic Data Analysis (FINAL WORKING VERSION)

### BIO253-OMICS -- Phenotypic Data Analysis (FINAL WORKING VERSION)

library(tidyverse)
library(readxl)
library(GGally)

### 1) Load data
pfad <- "C:/Users/rubid/Documents/1_Prio_Biomedizin/AAA/00_AAA_BIO253/Bio253_SA6850_Research_in_Omics_HS2025-main"
datei <- "tabel_STX.xlsx"
full_path <- file.path(pfad, datei)

df <- read_excel(full_path, sheet = "Tabelle2", skip = 3)

### 2) Clean numeric columns
df[df == "Missing"] <- NA

numeric_cols <- c("crtM", "crtN", "crtP", "crtQ", "crtO",
                  "STX Average", "Growth rate", "H2O2 survival",
                  "sigma-B", "rsbU", "rsbV", "rsbW")

df <- df %>%
    mutate(across(all_of(numeric_cols), ~as.numeric(.)))

### 3) Create delta_df (now **correct**)
delta_df <- df %>%
    pivot_wider(
        names_from = Condition,
        values_from = `Growth rate`
    ) %>%
    mutate(delta_growth = PASN - TSB) %>%    # <-- CORRECT COLUMN NAMES
    select(Clones, Strain, delta_growth)

### 4) Build phenotype dataframe
pheno <- df %>%
    distinct(Clones, Strain, `STX Average`, `H2O2 survival`) %>%
    left_join(delta_df, by = c("Clones", "Strain"))

print("Phenotype dataframe:")
print(pheno)

### 5) Correlation STX vs H2O2
cor.test(pheno$`STX Average`, pheno$`H2O2 survival`, method="pearson")
cor.test(pheno$`STX Average`, pheno$`H2O2 survival`, method="spearman")

ggplot(pheno, aes(x = `STX Average`, y = `H2O2 survival`, color = Strain)) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE, color="black") +
    theme_minimal()

### 6) Correlation STX vs delta growth
cor.test(pheno$`STX Average`, pheno$delta_growth, method="pearson")
cor.test(pheno$`STX Average`, pheno$delta_growth, method="spearman")

ggplot(pheno, aes(x = `STX Average`, y = delta_growth, color = Strain)) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE, color="black") +
    theme_minimal()

### 7) Scatter matrix
ggpairs(pheno,
        columns = c("STX Average", "H2O2 survival", "delta_growth"),
        aes(color = Strain))

#JE2 clones have higher STX and higher σᵇ activation
# → They may show better oxidative resistance.
#6850 clones have lower STX
# → Possibly lower H₂O₂ survival.
#Δ Growth (PASN – TSB) differs strongly between strains
# → Could correlate weakly or not at all with STX.
# You might find:
#STX ↔ H₂O₂ survival: positive correlation
#STX ↔ Δ Growth: weak or no correlation
