# WTAgingXMetabolism           

**Date created:** 7/7/26    
**Last Updated:** 7/7/26    
**Author:** Chloe Lucido     
**Project:** Wildtype Aging Cerebral Metabolism         
**Purpose:** Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data        

**Notes:**              
- MCC denotes scripts run using the supercomputer

## Scripts:      

### **SpaNorm/**      

**- _SpaNorm 02._** Running SpaNorm per sample, and creating sample objlist with corresponding SpaNorm assay attached and merged obj with merged SpaNorm assay attached     
**- _MCC_SpaNormprocessing_noint 03._** Processing merged obj using SpaNorm assay (scaleData -> dim reductions) - NO INTEGRATION        

#### Reference RPCA       

**- _MCC_SpaNorm_RPCA 03._** Running reference-based RPCA integration on SpaNorm obj split by sample, processing and findmarkers           
**- _MCC_SpaNorm_RPCA 04._** Creating feature plots, UMAPs and imagedimplots for annotation of SpaNorm RPCA integrated obj    
**- _MCC_SpaNorm_RPCA 05._** Subclustering clusters 9, 18, 1, and 3       
**- _MCC_SpaNorm_RPCA 06._** renaming clusters and adding old (from 6/11 and before) XE regional annotations        
**- _MCC_SpaNorm_RPCA 07._** making pretty figures for ICBEM and calculating cell proportions          
**- _SpaNorm_RPCA 08._** stats on cell proportions and making cell proportions line graphs         
**- _MCC_SpaNorm_RPCA 09._** metabolic pathway scoring (global and cell-specific z-scores and mean expression)
**- _SpaNorm_RPCA 10._** metabolic pathway score figures       
**- _MCC_SpaNorm_RPCA 11._** subsetting, subclustering and annotating astrocytes          
**- _MCC_SpaNorm_RPCA 12._** astro obj: metabolic pathway scoring and making pretty figs (dim plots)
**- _SpaNorm_RPCA 13._** astro obj: metabolic pathway figures

**Important Notes:**    
