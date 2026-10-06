#!/usr/bin/env Rscript

# Build the ratePlacer reference files for one tree from a reference alignment and a tree.
#
# Usage: Rscript make_reference.R ALIGNMENT.fasta TREE.newick TREE_NUM OUT_DIR [phy]
#
#   ALIGNMENT.fasta  Reference multiple sequence alignment. Sequence names are the FASTA
#                    header up to the first whitespace and must match the tree's tip labels.
#   TREE.newick      Tree of the reference sequences (e.g. RAxML GTRGAMMA). It is rooted
#                    at the midpoint, so it may be unrooted.
#   TREE_NUM         Tree number (0, 1, ...), used as the prefix of every output file.
#   OUT_DIR          Output directory (created if missing). Use the same directory for all
#                    trees; it is ratePlacer's -p directory.
#   phy              Optional: also write the alignment in PHYLIP format (needs seqmagick).
#
# Outputs, in OUT_DIR:
#   <TREE_NUM>_reference.txt    alignment + ultrametric tree with tips numbered 1..n (gfl/ratePlacer input)
#   <TREE_NUM>_parameter.txt    GTR+Gamma parameters (gfl/ratePlacer input)
#   <TREE_NUM>_pared.fa         alignment in tip-number order: the k-th sequence is leaf k
#   <TREE_NUM>_ultrametric.tree fitted ultrametric tree with the original tip names
#   <TREE_NUM>_pared.phy        (only with 'phy') alignment in PHYLIP format

suppressPackageStartupMessages({
  library(phangorn)
  library(Biostrings)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript make_reference.R ALIGNMENT.fasta TREE.newick TREE_NUM OUT_DIR [phy]", call. = FALSE)
}

fastaFile <- args[1]
treeFile <- args[2]
treeNum <- args[3]
outDir <- args[4]
writePhy <- length(args) >= 5 && args[5] == "phy"

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
outFile <- function(suffix) file.path(outDir, paste0(treeNum, "_", suffix))

# 1) Load MSA and tree

nt_dna <- readDNAStringSet(fastaFile)
# FASTA names are the header up to the first whitespace, as in most sequence tools
names(nt_dna) <- sub("[[:space:]].*$", "", names(nt_dna))
if (anyDuplicated(names(nt_dna))) {
  stop("Duplicate sequence names in ", fastaFile, ": ",
       paste(head(unique(names(nt_dna)[duplicated(names(nt_dna))]), 5), collapse = ", "), call. = FALSE)
}
nt_seqs <- as.phyDat(nt_dna)
rm(nt_dna); invisible(gc())

tree <- midpoint(read.tree(treeFile))

# 2) Pare the alignment to the sequences in the tree

missing_tips <- setdiff(tree$tip.label, names(nt_seqs))
if (length(missing_tips) > 0) {
  stop(length(missing_tips), " tree tip(s) have no sequence in ", fastaFile, ", e.g. ",
       paste(head(missing_tips, 5), collapse = ", "), call. = FALSE)
}
n_dropped <- sum(!(names(nt_seqs) %in% tree$tip.label))
if (n_dropped > 0) {
  message("Dropping ", n_dropped, " sequence(s) that are not in the tree")
}
nt_seqs_pared <- nt_seqs[which(names(nt_seqs) %in% tree$tip.label)]
rm(nt_seqs); invisible(gc())

# 3) Coerce the tree to be ultrametric

tree_ultra <- phangorn:::minEdge(tree, tau = 1e-5, enforce_ultrametric = TRUE)
rm(tree); invisible(gc())

# 4) Fit GTR+Gamma on the ultrametric tree
# The branch lengths are not meaningful yet, so they are re-optimized along with the model.
# k must stay 4: gfl and ratePlacer use 4 discrete gamma categories.

bf_pared <- baseFreq(nt_seqs_pared)
fit_ultra <- pml(tree_ultra, data = nt_seqs_pared, k = 4, bf = bf_pared, rate = 1.0)
rm(tree_ultra); invisible(gc())

fitGTR_ultra <- optim.pml(fit_ultra, model = "GTR", optRooted = T, optQ = T, optGamma = TRUE, optBf = TRUE,  optRate = FALSE, optNni = FALSE, rearrangement = "none", control = pml.control(tau=1e-5, maxit = 10))
rm(fit_ultra, bf_pared); invisible(gc())

# Write the fitted tree with the original tip names (without node labels)
newTree <- fitGTR_ultra$tree
newTree$node.label <- NULL
write.tree(newTree, file = outFile("ultrametric.tree"))
rm(newTree)

# 5) Write the parameter file: gamma shape, base frequencies (A C G T), GTR rates (AC AG AT CG CT GT)

param_file <- outFile("parameter.txt")
if (file.exists(param_file)) {
  file.remove(param_file)
}

write(fitGTR_ultra$shape, param_file, append = T)
write(sprintf("%.8f", fitGTR_ultra$bf), param_file, ncolumns=length(fitGTR_ultra$bf), append = T)	#like this due to rounding issues
write(fitGTR_ultra$Q, param_file, ncolumns = length(fitGTR_ultra$Q), append = T)

# 6) Write the reference file

## Tips in the order they appear in the Newick string; leaf k of the reference file is the k-th of these.
## write.tree writes tips in cladewise edge order, so read the order from the edge matrix.
get_tips_in_newick_order <- function(phy) {
  phy <- reorder(phy, "cladewise")
  tip_nodes <- phy$edge[, 2][phy$edge[, 2] <= Ntip(phy)]
  phy$tip.label[tip_nodes]
}

tips_in_newick_order <- get_tips_in_newick_order(fitGTR_ultra$tree)
tips_in_newick_order_indices <- match(tips_in_newick_order, names(nt_seqs_pared))
nt_seqs_reordered <- nt_seqs_pared[tips_in_newick_order_indices]
rm(nt_seqs_pared); invisible(gc())

##write reordered MSA as FASTA file
write.phyDat(nt_seqs_reordered, file = outFile("pared.fa"), format = "fasta")

##optionally write PHY file (enable with 5th argument "phy")
if (writePhy) {
  library(seqmagick)  # Only load if needed
  fas2phy(outFile("pared.fa"), outfile = outFile("pared.phy"), type = "sequential")
}

##convert phyDat directly to character matrix
nt_seqs_matrix <- as.character(nt_seqs_reordered)
rm(nt_seqs_reordered); invisible(gc())

ref_file <- outFile("reference.txt")
if (file.exists(ref_file)) {
  file.remove(ref_file)
}

##write sequences
seq_strings <- apply(nt_seqs_matrix, 1, paste, collapse = "")
writeLines(c(paste(nrow(nt_seqs_matrix), ncol(nt_seqs_matrix)), seq_strings), ref_file)
rm(nt_seqs_matrix, seq_strings); invisible(gc())

##number tree tips as they show in the newick tree
infile_tree <- fitGTR_ultra$tree
rm(fitGTR_ultra); invisible(gc())

label_map <- setNames(seq_along(tips_in_newick_order), tips_in_newick_order)
infile_tree$tip.label <- label_map[infile_tree$tip.label]
rm(tips_in_newick_order, label_map); invisible(gc())

##add tree to the reference file
infile_tree$node.label <- NULL
write.tree(infile_tree, file = ref_file, append = TRUE)
rm(infile_tree); invisible(gc())
