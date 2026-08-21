# WTAgingXMetabolism           

**Date created:** 8/21/26    
**Last Updated:** 8/21/26    
**Author:** Chloe Lucido     
**Project:** Wildtype Aging Cerebral Metabolism         
**Purpose:** Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data        

**Notes:**              
- MCC denotes scripts run using the supercomputer

## Scripts:      

### **RCTD_SPLIT/**      

**- _adataobj 02.py_** creating adata object to run RCTD-py (same format, cell ID prefixes, etc. as merged seurat obj)     
**- _RCTDPY 03.py_** running RCTD-py using adata obj (spatial_obj) from script 02 and reference obj (ABC 10xv3 scRNA seq.) and saving outputs as parquest files for reconstruct_rctd_from_rctdpy() in SPLIT        
**- _reconstruct_rctdpy 04.R_** converting rectd-py results to rctd spacexr S4 obj and RCTD post-processing            
**- _SPLIT 05.R_** running spatially aware SPLIT and SPLIT shift             
**- _SPLIT_SCT 06.R_** clean, re-attach images to split shift obj and run SCT per slide          

#### Reference RPCA       

**- _MCC_SPLIT_refRPCA 07._** Running reference-based RPCA integration and processing                   
**- _MCC_SPLIT_refRPCA 08._** Cluster at different creating feature plots, UMAPs and imagedimplots for annotation of RCTD SPLIT refRPCA obj         
         

**Important Notes:**    
