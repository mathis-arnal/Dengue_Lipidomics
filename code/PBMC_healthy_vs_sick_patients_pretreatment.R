# Healthy vs sick patients data preparation

library(readr)
library(dplyr)

## Merge the two datasets

file1 <- "data/pretreatment/healthy_vs_sick_patients/Healthy_patients.tsv"
file2 <- "data/pretreatment/healthy_vs_sick_patients/MISSE_Global_Lipidomics_results_selected_patients_D0.tsv"

df1 <- read_tsv(file1, col_types = cols(.default = "c"))
df2 <- read_tsv(file2, col_types = cols(.default = "c"))

# Verification du nom de la premiere colonne (feature)
col_name1 <- colnames(df1)[1]
col_name2 <- colnames(df2)[1]
cat(sprintf("\nNom de la colonne feature dans fichier 1: '%s'\n", col_name1))
cat(sprintf("Nom de la colonne feature dans fichier 2: '%s'\n", col_name2))

cat(sprintf("Fichier 1: %d lignes, %d colonnes\n", nrow(df1), ncol(df1)))
cat(sprintf("Fichier 2: %d lignes, %d colonnes\n", nrow(df2), ncol(df2)))

# Nettoyage des donnees AVANT la fusion
cat("\n", strrep("=", 60), "\n")
cat("NETTOYAGE DES DONNEES\n")
cat(strrep("=", 60), "\n")

# 1. Remplacer les virgules par des points-virgules dans la colonne feature
virgules_in_features_df1 <- sum(grepl(",", df1[[col_name1]]), na.rm = TRUE)
virgules_in_features_df2 <- sum(grepl(",", df2[[col_name2]]), na.rm = TRUE)

cat("\nVirgules detectees dans les noms de features:\n")
cat(sprintf("  - Fichier 1: %d lignes concernees\n", virgules_in_features_df1))
cat(sprintf("  - Fichier 2: %d lignes concernees\n", virgules_in_features_df2))

if (virgules_in_features_df1 > 0) {
  df1[[col_name1]] <- gsub(",", ";", df1[[col_name1]])
  cat("  -> Virgules remplacees par ';' dans fichier 1\n")
}

if (virgules_in_features_df2 > 0) {
  df2[[col_name2]] <- gsub(",", ";", df2[[col_name2]])
  cat("  -> Virgules remplacees par ';' dans fichier 2\n")
}

# 2. Remplacer les virgules par des points dans les colonnes numeriques
numeric_cols_df1 <- setdiff(colnames(df1), col_name1)
numeric_cols_df2 <- setdiff(colnames(df2), col_name2)

cat("\nConversion des virgules decimales en points:\n")
cat(sprintf("  - Fichier 1: %d colonnes numeriques a traiter\n", length(numeric_cols_df1)))
cat(sprintf("  - Fichier 2: %d colonnes numeriques a traiter\n", length(numeric_cols_df2)))

df1 <- df1 %>%
  mutate(across(all_of(numeric_cols_df1), ~ as.numeric(gsub(",", ".", .))))
df2 <- df2 %>%
  mutate(across(all_of(numeric_cols_df2), ~ as.numeric(gsub(",", ".", .))))

cat("  -> Conversion terminee\n")

cat("\nApercu apres nettoyage:\n")
cat("Fichier 1:\n")
print(head(df1, 3))
cat("\nFichier 2:\n")
print(head(df2, 3))

# Fusion des deux dataframes sur la colonne feature (inner join)
# Cela ne garde que les features presentes dans les deux fichiers
df_merged <- inner_join(df1, df2, by = setNames(col_name2, col_name1),
                         suffix = c("_healthy", "_D0"))

cat(sprintf("\nApres fusion: %d lignes, %d colonnes\n", nrow(df_merged), ncol(df_merged)))
cat(sprintf("Nombre de features communes: %d\n", nrow(df_merged)))

# Sauvegarde du resultat
output_file <- "data/pretreatment/healthy_vs_sick_patients/healthy_sick_lipidomics.tsv"
write_tsv(df_merged, output_file)

cat(sprintf("\nFichier fusionne sauvegarde: %s\n", output_file))

# Bilan des lignes non fusionnees
cat("\n", strrep("=", 60), "\n")
cat("BILAN DES LIGNES NON FUSIONNEES\n")
cat(strrep("=", 60), "\n")

features_df1 <- unique(df1[[col_name1]])
features_df2 <- unique(df2[[col_name2]])
features_merged <- unique(df_merged[[col_name1]])

only_in_df1 <- setdiff(features_df1, features_df2)
only_in_df2 <- setdiff(features_df2, features_df1)

cat(sprintf("\nFeatures uniquement dans Healthy_patients.tsv: %d\n", length(only_in_df1)))
if (length(only_in_df1) > 0) {
  cat("Exemples:\n")
  for (feature in only_in_df1) cat(sprintf("  - %s\n", feature))
}

cat(sprintf("\nFeatures uniquement dans MISSE_Global_Lipidomics_results_selected_patients_D0.tsv: %d\n",
            length(only_in_df2)))
if (length(only_in_df2) > 0) {
  cat("Exemples:\n")
  for (feature in only_in_df2) cat(sprintf("  - %s\n", feature))
}

cat("\nResume:\n")
cat(sprintf("  - Total features dans df1: %d\n", length(features_df1)))
cat(sprintf("  - Total features dans df2: %d\n", length(features_df2)))
cat(sprintf("  - Features communes (fusionnees): %d\n", length(features_merged)))
cat(sprintf("  - Features perdues (total): %d\n", length(only_in_df1) + length(only_in_df2)))
cat(sprintf("  - Taux de fusion: %.1f%%\n",
            length(features_merged) / max(length(features_df1), length(features_df2)) * 100))

## Group information table file preparation

# Liste des patients sick avec leur groupe (mild ou severe)
sick_patients <- c(
  "BS-082D0" = "mild",
  "BS-364D0" = "mild",
  "KT-193D0" = "severe",
  "KT-247D0" = "severe",
  "KT-312D0" = "severe",
  "KT-313D0" = "mild",
  "KT-347D0" = "severe",
  "KT-412D0" = "mild",
  "KT-417D0" = "severe",
  "KT-445D0" = "mild",
  "KT-522D0" = "severe",
  "KT-525D0" = "severe",
  "KT-537D0" = "mild",
  "KT-538D0" = "mild",
  "KT-539D0" = "mild",
  "KT-565D0" = "severe",
  "KT-663D0" = "mild",
  "KT-695D0" = "mild",
  "KT-705D0" = "severe",
  "KT-723D0" = "mild"
)

# Recuperation des noms de colonnes (patients) de df1, en excluant la colonne feature
healthy_patients <- setdiff(colnames(df1), col_name1)

# Creation du dataframe de groupes : patients sick puis patients healthy
df_groups <- bind_rows(
  tibble(
    sample_name = names(sick_patients),
    label_name = names(sick_patients),
    group = unname(sick_patients),
    pair = "NA"
  ),
  tibble(
    sample_name = healthy_patients,
    label_name = healthy_patients,
    group = "healthy",
    pair = "NA"
  )
)

# Sauvegarde du fichier
group_output_file <- "data/pretreatment/healthy_vs_sick_patients/group_information_table_healthy_vs_sick_patients_D0.tsv"
write_tsv(df_groups, group_output_file)

cat(sprintf("\nFichier de groupes cree: %s\n", group_output_file))
cat(sprintf("Nombre total d'echantillons: %d\n", nrow(df_groups)))
cat(sprintf("  - Patients healthy: %d\n", length(healthy_patients)))
cat(sprintf("  - Patients mild: %d\n", sum(sick_patients == "mild")))
cat(sprintf("  - Patients severe: %d\n", sum(sick_patients == "severe")))

cat("\nApercu du fichier de groupes:\n")
print(head(df_groups, 10))
