options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("seqinr", quietly = TRUE)) install.packages("seqinr"); library(seqinr)
if (!requireNamespace("coil", quietly = TRUE)) install.packages("coil"); library(coil)
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel"); library(parallel)

args <- commandArgs(trailingOnly = TRUE)

input_fasta <- "input.fasta"
input_csv <- "metadata.csv"
output_fasta <- "output.fasta"
output_csv <- "output_metadata.csv"
cores <- 1

for (i in seq(1, length(args), by = 2)) {
  key <- args[i]
  value <- args[i + 1]
  if (key == "-i") input_fasta <- value
  if (key == "-m") input_csv <- value
  if (key == "-o") output_fasta <- value
  if (key == "-c") output_csv <- value
  if (key == "-n") cores <- as.numeric(value)
}

options(mc.cores = min(cores, parallel::detectCores(logical = TRUE)))

barcodes <- read.fasta(input_fasta, as.string = TRUE)

parsed_names_data <- lapply(1:length(barcodes), function(i) {
  unlist(strsplit(names(barcodes)[[i]], "\\|"))
})

barcodes <- data.frame(
  name = sapply(parsed_names_data, function(x) x[[1]]),
  sequence = unname(unlist(barcodes))
)

process_sequence <- function(i) {
  tryCatch({
    result <- coi5p_pipe(barcodes$sequence[i], 
                         name = barcodes$name[i], 
                         trans_table = 5,
                         triple_translate = TRUE)
    cat("Sequence", i, "of", length(barcodes$name), "processed.\n")
    return(result)
  }, error = function(e) {
    cat("Error in sequence", i, ":", e$message, "\n")
    return(NULL)
  })
}

barcodes$coi_output <- mclapply(1:length(barcodes$name), process_sequence)

full_coi <- flatten_coi5p(barcodes$coi_output)
write.csv(full_coi, "coil_output.csv")
retained_coi <- full_coi[full_coi$indel_likely == FALSE & full_coi$stop_codons == FALSE, ]
write.fasta(sequences = as.list(retained_coi$raw), names = retained_coi$name, file.out = output_fasta)

metadata <- read.csv(input_csv)
filtered_metadata <- metadata[metadata$processid %in% retained_coi$name, ]
write.csv(filtered_metadata, file = output_csv)
