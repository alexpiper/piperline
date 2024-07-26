# Parse input arguments 
options <- commandArgs(trailingOnly = TRUE)
#options[options =="NA"] <- NA_character_
options <- as.list(options)

names(options) <- c(
  "sample_sheet",
  "run_parameters",
  "pcr_primers",
  "for_primer_seq",
  "rev_primer_seq",
  "target_gene",
  "max_primer_mismatch",
  "read_min_length",
  "read_max_length",
  "read_max_ee",
  "read_trunc_length",
  "read_trim_left",
  "read_trim_right",
  "asv_min_length",
  "asv_max_length",
  "high_sensitivity",
  "concat_unmerged",
  "genetic_code",
  "coding",
  "phmm",
  "idtaxa_db",
  "ref_fasta",
  "idtaxa_confidence",
  "run_blast",
  "blast_min_identity",
  "blast_min_coverage",
  "target_kingdom",
  "target_phylum",
  "target_class",
  "target_order",
  "target_family",
  "target_genus",
  "target_species",
  "min_sample_reads",
  "min_taxa_reads",
  "min_taxa_ra",
  "threads",
  "start_fresh"
)

print(options)

#library(renv)
library(pak)
library(targets)
library(tarchetypes)

suppressWarnings(suppressMessages(source("_targets_packages.R")))
suppressWarnings(suppressMessages(source("R/functions.R")))
suppressWarnings(suppressMessages(source("R/themes.R")))

# get number of primers
primer_length <- length(options$pcr_primers %>% str_split_1(",|;"))

# Check that number of provided parameters are either 1, or match the number of provided primers
param_lengths <- options[!names(options) %in% c("sample_sheet", "run_parameters", "pcr_primers")] %>%
  purrr::map_dbl(~length(.x %>% str_split_1(",|;")))

if(any(!param_lengths %in% c(1, primer_length))){
  stop(names(param_lengths[param_lengths %in% c(1, primer_length)]), " parameters must have either 1 argument, or the same as pcr_primers (", primer_length,")")
}

# Get sample sheet & runparameters
sample_sheet <- options$sample_sheet%>%
  str_split_1(",")
run_parameters <- options$run_parameters%>%
  str_split_1(",")

message("Creating new params file from input arguments")

# Create parameters sheet
params <- options[!names(options) %in% c("sample_sheet", "run_parameters")] %>%
  purrr::map(~{
    #if(is.na(.x)){return(NA)}
    .x %>%
      str_split_1(",|;") %>%
      str_remove("\\[.*\\]")}) %>%
  bind_rows()

write_csv(params, "sample_data/loci_params.csv")

# Create samplesheets
message("Creating new sample data file from input arguments")

# Create samplesheet containing samples and run parameters for all runs
samdf <- create_samplesheet(SampleSheet = sample_sheet, runParameters = run_parameters, template = "V4") %>%
  distinct()

# Check that sample_ids contain fcid, if not; attatch
samdf <- samdf %>%
  mutate(sample_id = case_when(
    !str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
    TRUE ~ sample_id
  ))

# Validate samples match the sample sheet
fastqFs <- purrr::map(list.dirs("data", recursive=FALSE),
                      list.files, pattern="_R1_", full.names = TRUE) %>%
  unlist() %>%
  str_remove(pattern = "^(.*)\\/") %>%
  str_remove(pattern = "(?:.(?!_S))+$")

fastqFs <- fastqFs[!str_detect(fastqFs, "Undetermined")]

missing_fastqFs <- setdiff(fastqFs, samdf$sample_id)
if (length(missing_fastqFs) > 0) warning("The fastq file/s: ", missing_fastqFs, " are not in the sample sheet")

missing_sample_ids <- setdiff(samdf$sample_id, fastqFs)
if (length(missing_sample_ids) > 0) {
  warning("The fastq file: ", missing_sample_ids, " is missing, dropping from samplesheet \n")
  samdf <- samdf %>% filter(!sample_id %in% missing_sample_ids)
}

# Add primer details to the sample sheet based on string matching

# Function to dynamically create case_when logic for primers
create_case_when <- function(input) {
  case_when_expr <- map2_chr(names(input), input, ~ sprintf("str_detect(sample_id, '%s') ~ '%s'", .x, .y))
  case_when_expr <- paste(case_when_expr, collapse = ", ")
  case_when_expr <- paste0("case_when(", case_when_expr, ", TRUE ~ NA_character_)")
  rlang::parse_expr(case_when_expr)
}

# Helper function to create named vectors
create_named_vector <- function(input) {
  input_values <- str_split_1(input, ",") %>% str_remove("\\[.*\\]")
  input_names <- str_split_1(input, ",") %>% str_remove("].*$") %>% str_remove("^.*\\[")
  setNames(input_values, input_names)
}

# Define named vectors for patterns and their corresponding primer sequences
pattern_to_primers <- create_named_vector(options$pcr_primers)
pattern_to_for <- create_named_vector(options$for_primer_seq)
pattern_to_rev <- create_named_vector(options$rev_primer_seq)

# Add primers to sample sheet based on sample_name matching
samdf <- samdf %>%
  mutate(pcr_primers = !!create_case_when(pattern_to_primers))%>%
  mutate(for_primer_seq = !!create_case_when(pattern_to_for))%>%
  mutate(rev_primer_seq = !!create_case_when(pattern_to_rev))

# Write out sample tracking sheet
write_csv(samdf, "sample_data/Sample_info.csv")

# Check if start_fresh is set
if(any(isTRUE(as.logical(options$start_fresh)))){
  message("Start fresh parameter is set - removing all prior targets before running")
  tar_destroy(destroy="all", ask=FALSE)
} else {
  message("Continuing pipeline from last run")
}

# run pipeline
tar_make(script = "_targets.R")
