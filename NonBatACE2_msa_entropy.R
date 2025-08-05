library(Biostrings)
library(msa)
library(seqinr)

nonbat_fasta_file <- "non-bat_ACE2.fasta"
nonbat_aa_seqences <- readAAStringSet(nonbat_fasta_file)
aligned2 <- msa(nonbat_aa_seqences, method = "Muscle")
aligned2_seqinr <- msaConvert(aligned2, type = "seqinr::alignment")
num_seqs <- length(aligned2_seqinr$nam)
seq_lengths <- unique(nchar(aligned2_seqinr$seq))
cat("Number of sequences aligned:", num_seqs, "\n")
cat("Length of aligned sequences:", paste(seq_lengths, collapse = ", "), "\n")
write.fasta(sequences = as.list(aligned2_seqinr$seq),
            names = aligned2_seqinr$nam,
            file.out = "non-bat_ACE2_alignment_for_ggmsa.fasta")
library(ggmsa)
library(ggplot2)
plot_faceted <- ggmsa("non-bat_ACE2_alignment_for_ggmsa.fasta",
                      color = "Zappo_AA") + 
  ggtitle("Non-Bat ACE2 Protein Alignment (Faceted View)") +
  facet_msa(field = 200) + 
  theme_minimal() 

print(plot_faceted)

ggsave("faceted_ACE2_alignment.png",
       plot = plot_faceted, 
       width = 15,        
       height = 8,         
       units = "in",        
       dpi = 300) 

#Per-Site Entropy (Sequence Variability)
alignment <- read.alignment("non-bat_ACE2_alignment_for_ggmsa.fasta", format ="fasta")
alignment_matrix2 <- as.matrix.alignment(alignment)
alignment_matrix2 <- apply(alignment_matrix2, 2, as.character)
entropy_without_gaps2 <- function(column) {
  column <- column[column != "-" & column != "." & column != "*"]
  if(length(column) == 0) return(0)
  freqs <- table(column)/length(column)
  -sum(freqs*log2(freqs))
} 

ent_values2 <- apply(alignment_matrix2, 2, entropy_without_gaps2)
plot(ent_values, type ="h",
     main = "Per-Site Entropy (Gaps Excluded)",
     xlab = "Alignment Position",
     ylab = "Entropy")
#Save entropy values for later analysis to compare with non-bat sequences
write.csv(ent_values2, "non-bat_entropy_values.csv", row.names = FALSE)

