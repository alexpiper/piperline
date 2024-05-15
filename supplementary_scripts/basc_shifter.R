
# commands modified from: https://github.com/jackscanlan/piperline/blob/main/jack_notes/basc_run.md

# Parse input arguments 
options <- commandArgs(trailingOnly = TRUE)

#library(renv)
library(pak)
library(targets)
library(tarchetypes)

suppressWarnings(suppressMessages(source("_targets_packages.R")))
suppressWarnings(suppressMessages(source("R/functions.R")))
suppressWarnings(suppressMessages(source("R/themes.R")))

# Params to add in step_add_parameters
# Positions of arguments are handled by the SLURM script argument passing and should be consistent here
params <- tibble(
  # Primer parameters
  pcr_primers = options[1],
  for_primer_seq = options[2],
  rev_primer_seq = options[3],

  target_gene= options[4],
  max_primer_mismatch= options[5],
  
  # Read filtering
  read_min_length = options[6],
  read_max_length = options[7],
  read_max_ee = options[8],
  read_trunc_length = options[9],
  read_trim_left = options[10],
  read_trim_right = options[11],
  
  # ASV filtering
  asv_min_length = options[12],
  asv_max_length = options[13],
  high_sensitivity = options[14],
  concat_unmerged = options[15],
  genetic_code = options[16],
  coding = options[17],
  phmm = options[18],
  
  # Taxonomic assignment
  idtaxa_db = options[19],
  ref_fasta = options[20],
  idtaxa_confidence = options[21],
  run_blast= options[22],
  blast_min_identity = options[23],
  blast_min_coverage = options[24],
  target_kingdom = options[25],
  target_phylum = options[26],
  target_class = options[27],
  target_order = options[28],
  target_family = options[29],
  target_genus = options[30],
  target_species= options[31],
  
  # Sample & Taxon filtering
  min_sample_reads = options[32],
  min_taxa_reads= options[33],
  min_taxa_ra = options[34],
  
  # General pipeline parameters
  threads = options[35],
)

write_csv(params, "sample_data/loci_params.csv")

runs <- dir("data/") #Find all directories within data
SampleSheet <- list.files(paste0("data/", runs), pattern= "SampleSheet", full.names = TRUE)
runParameters <- list.files(paste0("data/", runs), pattern= "[Rr]unParameters.xml", full.names = TRUE)

# Create samplesheet containing samples and run parameters for all runs
samdf <- create_samplesheet(SampleSheet = SampleSheet, runParameters = runParameters, template = "V4") %>%
  distinct()

# Check that sample_ids contain fcid, if not; attatch
samdf <- samdf %>%
  mutate(sample_id = case_when(
    !str_detect(sample_id, fcid) ~ paste0(fcid,"_",sample_id),
    TRUE ~ sample_id
  ))

# Check that samples match samplesheet
fastqFs <- purrr::map(list.dirs("data", recursive=FALSE),
                      list.files, pattern="_R1_", full.names = TRUE) %>%
  unlist() %>%
  str_remove(pattern = "^(.*)\\/") %>%
  str_remove(pattern = "(?:.(?!_S))+$")

# Filter undetermined reads from sample sheet
fastqFs <- fastqFs[!str_detect(fastqFs, "Undetermined")]

# Check for fastq files that are missing from samplesheet
if (length(setdiff(fastqFs, samdf$sample_id)) > 0) {warning("The fastq file/s: ", setdiff(fastqFs, samdf$sample_id), " are not in the sample sheet") }

# Check for sample_ids that dont have a corresponding fastq file
if (length(setdiff(samdf$sample_id, fastqFs)) > 0) {
  warning(paste0("The fastq file: ",
                 setdiff(samdf$sample_id, fastqFs),
                 " is missing, dropping from samplesheet \n")) 
  samdf <- samdf %>%
    filter(!sample_id %in% setdiff(samdf$sample_id, fastqFs))
}


# Add primers to sample sheet
samdf <- samdf %>%
  mutate(pcr_primers = options[1],
         for_primer_seq = options[2],
         rev_primer_seq = options[3]
  )

# Write out sample tracking sheet
write_csv(samdf, "sample_data/Sample_info.csv")

# run pipeline
tar_make(script = "_targets.R")
