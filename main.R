# main.R

## Process data files ##
source("scripts/1_process_gpr.R")
source("scripts/2_processed_data.R")
gc()  
rm(list = ls())  
## VAR-MERTON MODEL ##
source("scripts/3_zfactor.R", local = TRUE)
gc()  
rm(list = ls())  
# OIRF or PS GIRF is equivalent.
OIRF = TRUE
if (OIRF){
  source("scripts/4_oirf.R")  
}else{
  source("scripts/4_PSGIRF.R")
  }

