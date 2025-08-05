R Project: Bat vs Non-Bat ACE2 Analysis
This project compares ACE2 protein sequences between bat and non-bat species using entropy and conservation metrics.
Objective: Identify key amino acid residues that differ significantly between groups, possibly influencing viral receptor binding (e.g., SARS-CoV-2).
Methods used: Shannon Entropy (diversity at each amino acid site)
- Conservation score analysis (0–11 scale)
- Custom Divergence Score based on:
  - Consensus differences
  - Entropy differences
  - Conservation differences
  - Tools used: R packages: `readxl`, `msa`, `Biostrings`, `ggplot2`
- Jalview for sequence visualization
- 
