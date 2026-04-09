library(readxl)
library(dplyr)
library(purrr)
library(stringr)

input_file <- "data/raw/raw_Data_morpho.xlsx"
output_rds <- "data/derived/dat_clean.rds"
output_csv <- "data/derived/dat_clean.csv"

sheet_map <- tibble(sheet = excel_sheets(input_file)) %>%
  mutate(
    dataset = if_else(str_detect(str_to_lower(sheet), "^simulated"), "simulated", "known"),
    region = if_else(str_detect(str_to_lower(sheet), "cranial"), "cranial", "pelvic")
  )

dat_clean <- pmap_dfr(sheet_map, function(sheet, dataset, region) {
  read_excel(input_file, sheet = sheet) %>%
    mutate(dataset = dataset, region = region, sheet = sheet, .before = 1)
})

dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
saveRDS(dat_clean, output_rds)
write.csv(dat_clean, output_csv, row.names = FALSE)

print(dplyr::count(dat_clean, dataset, region))
