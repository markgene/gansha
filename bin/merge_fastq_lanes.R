library(optparse)
library(dplyr)
library(stringr)

message("# Given an input directory of FASTQs, create the command to merge the same samples loaded on different lanes.")

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="input directory"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="output directory"),
  make_option(c("-r", "--recursive"), action="store", default=FALSE, type='character',
              help="output directory")
)
opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$input)) {
  stop("Require input directory")
}

if (is.na(opt$output)) {
  stop("Require output directory")
}

if (!dir.exists(opt$output)) dir.create(opt$output, recursive = TRUE)

# Test
# opt <-
#   list(input = "/rgs01/scratch/users/jchen4/projects/193662", 
#        output = "/rgs01/scratch/users/jchen4/pipeline/cr/dev/merged_fastq_lanes",
#        recursive = FALSE)
# fastq_files <-
#   c(
#     "/rgs01/scratch/users/jchen4/projects/193662/1974039_HA_S13_L001_R1_001.fastq.gz",
#     "/rgs01/scratch/users/jchen4/projects/193662/1974039_HA_S13_L001_R2_001.fastq.gz",
#     "/rgs01/scratch/users/jchen4/projects/193662/1974039_HA_S13_L002_R1_001.fastq.gz",
#     "/rgs01/scratch/users/jchen4/projects/193662/1974039_HA_S13_L002_R2_001.fastq.gz",
#     "/rgs01/scratch/users/jchen4/projects/193662/1974040_HC_S14_L001_R1_001.fastq.gz",
#     "/rgs01/scratch/users/jchen4/projects/193662/1974040_HC_S14_L001_R2_001.fastq.gz"
#   )

fastq_files <- list.files(path = opt$input, recursive = opt$recursive, pattern = "_001.fastq.gz$", full.names = TRUE)

# Create sample data frame
data.frame(Fullpath = fastq_files) %>%
  dplyr::mutate(Filename = basename(Fullpath)) %>%
  dplyr::mutate(Sample = stringr::str_remove(Filename, "_L00[0-9]{1}_R[12]{1}_001.fastq.gz$")) %>%
  dplyr::mutate(R12 = stringr::str_replace(Filename, "^.*_L00[0-9]{1}_(R[12]{1})_001.fastq.gz$", "\\1")) %>%
  dplyr::mutate(Lane = stringr::str_replace(Filename, "^.*_(L00[0-9]{1})_(R[12]{1})_001.fastq.gz$", "\\1")) -> sample_df
# Check if each sample was loaded on the same lanes.
sample_df %>%
  dplyr::arrange(Sample, R12, Lane) %>%
  dplyr::group_by(Sample, R12) %>%
  dplyr::summarise(
    Lanes = paste(Lane, sep = ", ", collapse = ", ")
  ) %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(is_good = ifelse(length(unique(Lanes)) == 1, TRUE, FALSE)) -> check_df

if (!all(check_df$is_good)) {
  bad_samples <- check_df$Sample[!check_df$is_good]
  stop(paste("Some samples have different sets of lanes, e.g. a sample has R1 on L001 and R2 on L002.\n", bad_samples, collapse =  " "))
}

sample_df %>%
  dplyr::arrange(Sample, R12, Lane) %>%
  dplyr::group_by(Sample, R12) %>%
  dplyr::summarise(
    Infile = paste(Fullpath, collapse = " "),
    Outfile = unique(paste0(Sample, "_", R12, ".fq.gz")),
    Outpath = file.path(opt$output, Outfile),
    Cmd = paste("cat", Infile, ">", Outpath, collapse = " ")
  ) %>%
  dplyr::ungroup() -> res_df

cat(res_df$Cmd, sep = "\n")
