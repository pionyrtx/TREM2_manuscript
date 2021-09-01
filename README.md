# TREM2_manuscript
Novel code and cell IDs used for the manuscript titled "Targeting TREM2 on tumor associated macrophages enhances immunotherapy"

### monomac_final_object_cells.rds
Character vector of cell IDs for human DTC monocyte-macrophage object in the paper. 
Sample IDs are appended to cell ID to maintain uniqueness.

### cd8cd4_final_object_cells.rds
Character vector of the cell IDs used for the cd8cd4 T cell object in the paper.
Sample IDs are appended to cell ID to maintain uniqueness.

### ct26_myeloid_object_cells.rds
Character vector of the cell IDs used for the myeloid subset of the CT26 4-arm experiment in the paper. 
Cell IDs are appended to the sample ID (iso, t2, pd1, or combo) to maintain uniqueness.  

### ct26_lymphoid_object_cells.rds
Character vector of the cell IDs used for the lymphoid subset of the CT26 4-arm experiment in the paper. 
Cell IDs are appended to the sample ID (iso, t2, pd1, or combo) to maintain uniqueness.  

### ovarian_obj_trem2_paper_version_cells.rds
Character vector of the cell IDs used for the untreated ovarian DTC object in the paper.   

### two_object_sample_correlation.R
Function used for creating the correlation heatmaps in the paper, this also performs the Monte-Carlo permutation test for significance and output p-values.   
