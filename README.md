# WTAgingXMetabolism           

**Date created:** 4/22/26    
**Last Updated:** 6/12/26    
**Author:** Chloe Lucido     
**Project:** Wildtype Aging Cerebral Metabolism         
**Purpose:** Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data        

**Notes:**              
- MCC denotes scripts run using the supercomputer 

## Scripts:      

**- _01._** Creating xenium objs for each slide (7 objs total)    

### **LogNormalize/**         
          
**- _02._** Merging all objs and processing (NormalizeData->FindVariableFeatures->ScaleData)    
**- _03._** Annotating all cells on log normalize and scaled obj                   
**- _04._** Subsetting, subclustering, and annotating neuronal subtypes      
**- _05._** Subsetting, subclustering, and annotating vascular mural cell types      
**- _06._** Metabolic Scoring        

### **SCT/**    

**- _MCC_SCT_02._** Merging all slide objs and processing (SCTransform) - NO INTEGRATION      
**- _SCT_03._** Finding markers and annotating SCT (not integrated) obj

#### Reference RPCA  

**- _MCC_int 02_1_** Creating object list split by slide and SCT on each for RPCA integration                
**- _MCC_int 02_2_v2_** Integrating using reference-based RPCA              

#### Harmony integration    

**- _MCC_harmonyint_03_v2_** Integrating SCTransformed object using harmony         


### **SpaNorm/**      

**- _SpaNorm 02._** Running SpaNorm per sample, and creating sample objlist with corresponding SpaNorm assay attached and merged obj with merged SpaNorm assay attached     
**- _MCC_SpaNormprocessing_noint 03._** Processing merged obj using SpaNorm assay (scaleData -> dim reductions) - NO INTEGRATION        

#### Reference RPCA       

**- _MCC_SpaNorm_RPCA 03._** Running reference-based RPCA integration on SpaNorm obj split by sample, processing and findmarkers           
**- _MCC_SpaNorm_RPCA 04._** Creating feature plots, UMAPs and imagedimplots for annotation of SpaNorm RPCA integrated obj    
**- _MCC_SpaNorm_RPCA 05._** Subclustering clusters 9 and 18 and annotating         



**Important Notes:**    

**R session and packages Info**   
