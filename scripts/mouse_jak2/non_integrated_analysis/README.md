This folder contains scripts for scRNASeq analysis on non-integrated data since collaborator wanted to check if integration was affecting comparative studies between samples.

1. R_240617_check_non-integrated_5samples.Rmd:
   *  Checking quality of non-integrated samples.
     
2. R_240617_MS_2-6_non_integrated_Clustering.Rmd:
   *  Performed clustering on non-integrated data.
   *  Retaining only Cluster 0 since that has the information of interest with 4 different sub-clusters based on sample conditions (JAK2V617F_HMb1_1, JAK2V617F_Hamster_IgG C57BL6J_Hamster_IgG, C57BL6J_HMb1_1 and, C57BL6J_HMb1_1).
   *  Further isolated JAK2 since that was of interest. Identified DEGs for both conditions.
    
3. R_240618_DAVID_Gene_Enrichment.Rmd:
   *  Performed gene enrichment analysis using DAVID by identifying top 500 DEGs for both JAK2 conditions.
   *  Plotted results using treemap and scatterplot.
   *  Additionally plotted pathway enrichment plot for comparison with ClusterProfiler results.

4. R_240620_ClusterProfiler_Gene_Enrichment.Rmd:
   *  Performed gene enrichment using ClusterProfiler on all positive markers identified during clustering for JAK2 conditions.
   *  Plotted pathway enrichment plot for comparision with DAVID results.
