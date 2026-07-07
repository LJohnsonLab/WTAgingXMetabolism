# WTAgingXMetabolism           

**Date created:** 7/7/26    
**Last Updated:** 7/7/26    
**Author:** Chloe Lucido     
**Project:** Wildtype Aging Cerebral Metabolism         
**Purpose:** Analysis of 4 different age groups(16-, 36-, 59-, and 92-wks; n=3M / 3F per group) of WT mice Xenium ST data        

**Notes:**              
- MCC denotes scripts run using the supercomputer

## Scripts:      

### **SCT/**    

**- _MCC_SCT_02._** Merging all slide objs and processing (SCTransform) - NO INTEGRATION      
**- _SCT_03._** Finding markers and annotating SCT (not integrated) obj

#### Reference RPCA  

**- _MCC_int 02_1_** Creating object list split by slide and SCT on each for RPCA integration                
**- _MCC_int 02_2_v2_** Integrating using reference-based RPCA              

#### Harmony integration    

**- _MCC_harmonyint_03_v2_** Integrating SCTransformed object using harmony         


**Important Notes:**    
