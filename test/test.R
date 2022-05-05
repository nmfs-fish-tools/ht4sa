

# load packages
library(ht4sa)
library(Rmpi)
library(r4ss)

#source("../R/ht4sa.R")


#_____________________________________________________________________________________________________________________________
# set directory paths
proj_dir = "/Users/mattadmin/ht4sa/test/"
dir_utility = "/Users/mattadmin/SAP/Code/Utilities/"
dir_model = paste0(proj_dir,"test_case/")
dir_model_source = "input_files"
model_name_stem = "input_files/"
dir_input_files = paste0(dir_model,dir_model_source,"/")
dir_ss = paste0(proj_dir,"input_files/")
dir.create(paste0(proj_dir,"Summary/",model_name_stem),recursive=TRUE,showWarnings=FALSE)

print(dir_input_files)
#_____________________________________________________________________________________________________________________________
# set-up testing options
    testing_options_df = expand.grid(steepness=c(0.45,0.6,0.75),
                                     sigmaR=c(0.3,0.4,0.5),
                                     stringsAsFactors=FALSE)


#run cases in parallel
b<- run_ht4sa_ss_MPI(testing_options_df,
                             proj_dir,
                             dir_utility,
                             dir_model,
                             dir_model_source,
                             model_name_stem,
                             dir_input_files,
                             dir_ss)
