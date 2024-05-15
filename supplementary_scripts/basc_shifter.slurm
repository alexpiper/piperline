#!/bin/bash

# The name of the job:
#SBATCH --job-name="piperline"

# The project ID which this job should run under:
#SBATCH --account="pathogens"

# The partition in which to submit the job:
#SBATCH --partition="batch"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# The total amount of memory in megabytes in the job:
#SBATCH --mem=40GB

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=7-0:0:00

# Output errors and logs into same file
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.out

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi


# Argument parsing
while [ $# -gt 0 ] ; do
  case $1 in
    --pcr_primers) pcr_primers="$2" ;;
    --for_primer_seq) for_primer_seq="$2" ;;
    --rev_primer_seq) rev_primer_seq="$2" ;;
    --target_gene) target_gene="$2" ;;
    --max_primer_mismatch) max_primer_mismatch="$2" ;;
    --read_min_length) read_min_length="$2" ;;
    --read_max_length) read_max_length="$2" ;;
    --read_max_ee) read_max_ee="$2" ;;
    --read_trunc_length) read_trunc_length="$2" ;;
    --read_trim_left) read_trim_left="$2" ;;
    --read_trim_right) read_trim_right="$2" ;;
    --asv_min_length) asv_min_length="$2" ;;
    --asv_max_length) asv_max_length="$2" ;;
    --high_sensitivity) high_sensitivity="$2" ;;
    --concat_unmerged) concat_unmerged="$2" ;;
    --genetic_code) genetic_code="$2" ;;
    --coding) coding="$2" ;;
    --phmm) phmm="$2" ;;
    --idtaxa_db) idtaxa_db="$2" ;;
    --ref_fasta) ref_fasta="$2" ;;
    --idtaxa_confidence) idtaxa_confidence="$2" ;;
    --run_blast) run_blast="$2" ;;
    --blast_min_identity) blast_min_identity="$2" ;;
    --blast_min_coverage) blast_min_coverage="$2" ;;
    --target_kingdom) target_kingdom="$2" ;;
    --target_phylum) target_phylum="$2" ;;
    --target_class) target_class="$2" ;;
    --target_order) target_order="$2" ;;
    --target_family) target_family="$2" ;;
    --target_genus) target_genus="$2" ;;
    --target_species) target_species="$2" ;;
    --min_sample_reads) min_sample_reads="$2" ;;
    --min_taxa_reads) min_taxa_reads="$2" ;;
    --min_taxa_ra) min_taxa_ra="$2" ;;
    --threads) threads="$2" ;;
    
  esac
  shift
done

# Print parameters
echo Parameters used for analysis:
echo pcr_primers=$pcr_primers
echo for_primer_seq=$for_primer_seq
echo rev_primer_seq=$rev_primer_seq
echo target_gene=$target_gene
echo max_primer_mismatch=$max_primer_mismatch
echo read_trunc_length=$read_trunc_length
echo read_min_length=$read_min_length
echo read_max_length=$read_max_length
echo read_max_ee=$read_max_ee
echo read_trim_left=$read_trim_left
echo read_trim_right=$read_trim_right
echo asv_min_length=$asv_min_length
echo asv_max_length=$asv_max_length
echo high_sensitivity=$high_sensitivity
echo concat_unmerged=$concat_unmerged
echo genetic_code=$genetic_code
echo coding=$coding
echo phmm=$phmm
echo idtaxa_db=$idtaxa_db
echo ref_fasta=$ref_fasta
echo idtaxa_confidence=$idtaxa_confidence
echo run_blast=$run_blast
echo blast_min_identity=$blast_min_identity
echo blast_min_coverage=$blast_min_coverage
echo target_kingdom=$target_kingdom
echo target_phylum=$target_phylum
echo target_class=$target_class
echo target_order=$target_order
echo target_family=$target_family
echo target_genus=$target_genus
echo target_species=$target_species
echo min_sample_reads=$min_sample_reads
echo min_taxa_reads=$min_taxa_reads
echo min_taxa_ra=$min_taxa_ra
echo threads=$threads

# Create parameters file

module load shifter

# Update shifter image
# shifterimg pull docker:jackscanlan/piperline-multi:0.0.1

# Run pipeline in shifter
echo Running pipeline
shifter --image jackscanlan/piperline-multi:0.0.1 --env=R_LIBS_USER=/usr/local/lib/R/site-library \
	Rscript --vanilla ./supplementary_scripts/basc_shifter.R \
  $pcr_primers \
  $for_primer_seq \
  $rev_primer_seq \
  $target_gene \
  $max_primer_mismatch \
  $read_min_length \
  $read_max_length \
  $read_max_ee \
  $read_trunc_length \
  $read_trim_left \
  $read_trim_right \
  $asv_min_length \
  $asv_max_length \
  $high_sensitivity \
  $concat_unmerged \
  $genetic_code \
  $coding \
  $phmm \
  $idtaxa_db \
  $ref_fasta \
  $idtaxa_confidence \
  $run_blast \
  $blast_min_identity \
  $blast_min_coverage \
  $target_kingdom \
  $target_phylum \
  $target_class \
  $target_order \
  $target_family \
  $target_genus \
  $target_species \
  $min_sample_reads \
  $min_taxa_reads \
  $min_taxa_ra \
  $threads 

exit

# Output useful job stats
/usr/local/bin/showJobStats.scr 
