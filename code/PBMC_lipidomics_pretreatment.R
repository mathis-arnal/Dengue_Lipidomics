# Lipidomics

library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)


# Prepare dataset for LipidSig

# Read the "Normalization" sheet (values already divided by number of cells, i.e.
# ng per M of cells), not the default "Approached_Quantification" sheet (ISTD-corrected
# ng only, not yet normalized to cell count) - mirrors the plasma pipeline reading
# "Normalisation_µL_plasma" rather than "Approached_Quantification"/"Raw datas".
# This sheet has a different layout than Approached_Quantification (no ISTD column, an
# extra "M of cells" row), so it's parsed positionally rather than via skip=/rename().
raw_pbmc <- read_excel(
  "data/raw_data/260204_PBMC_DMI_Global_Results.xlsx",
  sheet = "Normalization", col_names = FALSE
)

pbmc_sample_names <- as.character(raw_pbmc[4, 3:ncol(raw_pbmc)])

pbmc_family <- raw_pbmc[[1]]
for (i in seq_along(pbmc_family)) if (i > 1 && is.na(pbmc_family[i])) pbmc_family[i] <- pbmc_family[i - 1]

pbmc_data_rows <- 8:nrow(raw_pbmc)
abundance <- raw_pbmc[pbmc_data_rows, ]
abundance[[1]] <- pbmc_family[pbmc_data_rows]
colnames(abundance) <- c("family", "feature", pbmc_sample_names)
abundance <- abundance %>% mutate(across(all_of(pbmc_sample_names), as.numeric))

# column_order <- c(
#   "BS-064", "JV-035", "JV-048", "KT-515", "KT-926", "JV-148", "JV-157",
#   "KT-193", "KT-247", "KT-312", "KT-347", "KT-352", "KT-417", "KT-434",
#   "KT-522", "KT-525", "KT-565", "KT-612", "KT-705", "KT-718", "JV-071",
#   "BS-082", "BS-336", "BS-346", "BS-351", "BS-364", "BS-377", "JV-138",
#   "KT-313", "KT-412", "KT-445", "KT-537", "KT-538", "KT-539", "KT-663",
#   "KT-695", "KT-716", "KT-723", "KT-771", "KT-805", "KT-880", "KT-974",
#   "BS-671"
# )

# "KT-880 D0" -> "KT-880D0" (drop the space before the timepoint suffix)
colnames(abundance) <- sub(" (D0|D03|D10|D60)$", "\\1", colnames(abundance))

# Drop D0 for healthy patients 
# "BS-126-RAT-05-M2D0" -> "BS-126-RAT-05-M2"
colnames(abundance) <- sub("(-M[0-9]+)(D0|D03|D10|D60)$", "\\1", colnames(abundance))

colnames(abundance)

# change (1,2) DG 32:0 to DG 32:0
abundance$feature <- gsub("\\(\\d+,\\d+\\) DG", "DG", abundance$feature)
# change Cer 34:1,O2 to Cer 34:1;O2
abundance$feature <- gsub(",", ";", abundance$feature)

# Get rid of the total rows, which are not needed 
abundance <- abundance %>% dplyr::filter(!grepl("Total", feature))

write_tsv(abundance, "data/PBMC/PBMC_lipid_abundance_data_D0.tsv")


# ------------------------------------------------------------------------------------------------
# Comparison of the lipid names between plasma samples and PBMC samples
PBMC_lipid_names <- abundance$feature
plasma_lipid_names  <- read_tsv("data/lipidsig_datasets/healthy_vs_sick_patients/healthy_sick_lipidomics.tsv")  %>% 
pull(feature)

print(sprintf("Number of lipids in PBMC dataset: %d", length(PBMC_lipid_names)))
print(sprintf("Number of lipids in plasma dataset: %d", length(plasma_lipid_names)))
matched_lipids <- intersect(PBMC_lipid_names, plasma_lipid_names)
print(sprintf("Number of matched lipids between PBMC and plasma datasets: %d", length(matched_lipids)))
print("Difference in PBMC dataset:")
setdiff(PBMC_lipid_names, plasma_lipid_names)
print("Difference in plasma dataset:")
setdiff(plasma_lipid_names, PBMC_lipid_names)

# ------------------------------------------------------------------------------------------------
# Create the group info 

## Group information table file preparation
sick_patients_df <- read.xlsx("data/raw_data/Final_Sample Shipment_ Montpellier_patient list_Dorothee_08082024 (2).xlsx",startRow = 2 )
# Liste des patients sick avec leur groupe (mild ou severe)

sick_patients_status <- sick_patients_df %>% 
  dplyr::select(`Patient.code`, `1997.Classification`) %>% 
  tidyr::drop_na(`1997.Classification` ) %>% 
  mutate(status = ifelse(`1997.Classification`  %in% c("DSS", "DHF"), "severe", "mild"),
         patient = paste0(`Patient.code`,"D0")) %>% 
  dplyr::select(patient, status)

# Fix sample and patient code mismatch

sick_patients_status$patient <- sub("KT-312D0", "KT-310D0", sick_patients_status$patient)
sick_patients_status$patient <- sub("JV-035D0", "BS-035D0", sick_patients_status$patient)

healthy_patients_status <- read.xlsx("data/raw_data/Final_Sample Shipment_ Montpellier_patient list_Dorothee_08082024 (2).xlsx",
                              sheet = "HD controls") %>% 
  dplyr::select(X1) %>% 
  mutate(patient = X1,
         status = "healthy"
  ) %>% 
  dplyr::select(patient, status)


all_patients <-rbind(sick_patients_status, healthy_patients_status)

# Creation du dataframe de groupes
df_groups <- 
  tibble(
    sample_name = all_patients$patient,
    label_name = all_patients$patient,
    group = all_patients$status,
    pair = "NA"
  )

# Sauvegarde du fichier
group_output_file <- "data/PBMC/group_information_table_healthy_vs_sick_patients_D0_all.tsv"
write_tsv(df_groups, group_output_file)

cat(sprintf("\nFichier de groupes cree: %s\n", group_output_file))
cat(sprintf("Nombre total d'echantillons: %d\n", nrow(df_groups)))
cat(sprintf("  - Patients healthy: %d\n", sum(df_groups$group == "healthy")))
cat(sprintf("  - Patients mild: %d\n", sum(df_groups$group == "mild")))
cat(sprintf("  - Patients severe: %d\n", sum(df_groups$group == "severe")))

cat("\nApercu du fichier de groupes:\n")
print(head(df_groups, 10))


# ------------------------------------------------------------------------------------------------
# Create the day of fever csv 
sick_patients_fever <- sick_patients_df %>%
  dplyr::select(`Patient.code`, `Day.of.Fever`) %>%
  rename(patient = `Patient.code`, day_of_fever = `Day.of.Fever`) %>%
  tidyr::drop_na(day_of_fever) %>%
  dplyr::select(patient, day_of_fever)

write_tsv(sick_patients_fever, "data/sick_patients_day_of_fever.tsv")  
