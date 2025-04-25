Directory contains scripts used to benchmark scGND.

Cell imbalances were introduced in the datasets to study how well the tools handle integrating data which contains imbalances. These were introduced by removing a randomly selected cell type in varied percentages (25%, 50%, 75% and 100%) and the UMAP and PCA embeddings were extracted. These embeddings were then scored using KBet and comparison across the tools was carried out using this method.

Information on the datasets used:
1. PBMC10K:
   *  scRNASeq data obtained from PBMC and contains ~10K cells.
   *  Additionally there was PBMC10K_rm25pct, PBMC10K_rm50pct, PBMC10K_rm75pct and PBMC10K_rm100pct which had varied cell imbalances introduced.
2. PBMC23K:
   *  scRNASeq data obtained from PBMC and contains ~23K cells.
   *  Additionally there was PBMC23K_rm25pct, PBMC23K_rm50pct, PBMC23K_rm75pct and PBMC23K_rm100pct which had varied cell imbalances introduced.
3. BMMC:
   *  scRNASeq data obtained from BMMC and contains ~50K cells.
   *  Additionally there was BMMC_rm25pct, BMMC_rm50pct, BMMC_rm75pct and BMMC_rm100pct which had varied cell imbalances introduced.

The file name contains the tool used to perform standard scRNASeq steps until data integration to extract the required embeddings.
