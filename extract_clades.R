suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: extract_clades.R <tree_file> <subclade_file> <output_dir>\n", 
       "Example: ./extract_clades.R tree.tre subclades.txt subtrees")
}

tree_file <- args[1]
subclade_file <- args[2]
output_dir <- args[3]

extract_clades <- function(tree_file, subclade_file, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  tree <- read.tree(tree_file)
  subclade_data <- readLines(subclade_file)
  clades <- list()
  current_clade <- NULL
  subclade_name <- NULL
  for (line in subclade_data) {    
    if (grepl("^Subclade [0-9]+:", line)) {
      if (!is.null(current_clade) && length(current_clade) > 1) {
        clades[[subclade_name]] <- current_clade
      }
      subclade_name <- sub("Subclade ([0-9]+):", "\\1", line)
      if (nchar(subclade_name) == 0) {
        stop("Subclade name extraction failed for line: ", line)
      }
      current_clade <- c()
    } else {
      current_clade <- c(current_clade, trimws(line))
    }
  }

  if (!is.null(current_clade) && length(current_clade) > 1) {
    clades[[subclade_name]] <- current_clade
  }

  for (subclade_name in names(clades)) {
    clade_tips <- clades[[subclade_name]]
    valid_tips <- intersect(tree$tip.label, clade_tips)

    if (length(valid_tips) > 1) {
      mrca_node <- getMRCA(tree, valid_tips)
      root_length <- tree$edge.length[which(tree$edge[, 2] == mrca_node)]
      subtree <- drop.tip(tree, setdiff(tree$tip.label, valid_tips))

      if (length(root_length) == 1 && !is.na(root_length)) {
        subtree$root.edge <- root_length
      } else {
        cat("Warning: No root length found for subclade_", subclade_name, ". Root edge not set.\n")
      }

      output_file <- file.path(output_dir, paste0("subclade_", subclade_name, ".tre"))
      write.tree(subtree, file = output_file)
      cat("Saved subtree to:", output_file, "\n")
    } else {
      cat("Skipping clade: subclade_", subclade_name, "- only one sequence, skipping.\n")
    }
  }
}

# Run the function
extract_clades(tree_file, subclade_file, output_dir)
