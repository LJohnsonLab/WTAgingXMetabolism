# WTAgingXMetabolism           

**Date created:** 4/22/26    
**Last Updated:** 5/5/26    
**Author:** Chloe Lucido     
**Project:** Wildtype Aging Cerebral Metabolism         
**Purpose:** Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data        

## Scripts:      

**- _01._** Creating xenium objs for each slide (7 objs total)    

### **NORMALIZATION AND SCALING SCRIPTS**         

#### Normalize and Scaledata           
**- _02._** Merging all objs and processing (NormalizeData->FindVariableFeatures->ScaleData)    

#### SCTransform       
**- _MCC 02._** Merging all objs and processing (SCTransform) *ONLY TO BE RUN ON SUPERCOMPUTER*      

#### SpaNorm      
**- _SpaNorm 02._** Running SpaNorm per slide, and creating slide objlist with corresponding SpaNorm assay attached and merged obj with merged SpaNorm assay attached         
**- _MCC SpaNormprocessing_noint 03._** Processing merged obj using SpaNorm assay (scaleData -> dim reductions) NO INTEGRATION        

### **INTEGRATION SCRIPTS**   

#### reference RPCA       
**- _MCC_int 02_1_** Creating object list split by slide and SCT on each for RPCA integration *ONLY TO BE RUN ON SUPERCOMPUTER*        
**- _MCC_int 02_2_** Integrating using RPCA *ONLY TO BE RUN ON SUPERCOMPUTER*      
**- _MCC_SpaNorm_int 03._**         

#### Harmony integration         
**- _MCC_harmonyint 03._** Integrating using harmony method *ONLY TO BE RUN ON SUPERCOMPUTER*     



### **ANNOTATION SCRIPTS**        
**- _SCT 03._** Finding markers and annotating cells on SCTransformed obj         
**- _03._** Annotating all cells on log normalize and scaled obj              
**- _04._** Subsetting, subclustering, and annotating neuronal subtypes      
**- _05._** Subsetting, subclustering, and annotating vascular mural cell types      
**- _06._** Metabolic Scoring 

**Important Notes:**    

**R session and packages Info**   
