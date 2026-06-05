# WTAgingXMetabolism           

**Date created:** 4/22/26    
**Last Updated:** 5/5/26    
**Author:** Chloe Lucido     
**Project:** Wildtype Aging Cerebral Metabolism         
**Purpose:** Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data        

Scripts:        
**- 01.** Creating xenium objs for each slide (7 objs total)    
**NORMALIZATION METHODS**      
> Normalize and Scaledata 
     
**- 02.** Merging all objs and processing (NormalizeData->FindVariableFeatures->ScaleData)     
**- MCC 02.** Merging all objs and processing (SCTransform) *ONLY TO BE RUN ON SUPERCOMPUTER*    
**INTEGRATION SCRIPTS**       
**- MCC_int 02_1** Creating object list split by slide and SCT on each for RPCA integration *ONLY TO BE RUN ON SUPERCOMPUTER*        
**- MCC_int 02_2** Integrating using RPCA *ONLY TO BE RUN ON SUPERCOMPUTER*      
**- MCC_harmonyint 03.** Integrating using harmony method *ONLY TO BE RUN ON SUPERCOMPUTER*         
**ANNOTATION SCRIPTS**        
**- SCT 03.** Finding markers and annotating cells on SCTransformed obj         
**- 03.** Annotating all cells on log normalize and scaled obj              
**- 04.** Subsetting, subclustering, and annotating neuronal subtypes      
**- 05.** Subsetting, subclustering, and annotating vascular mural cell types      
 

**Important Notes:**    

**R session and packages Info**   
