

# load packages
library(ht4sa)
library(Rmpi)
library(r4ss)
library(pointr)
library(tidyverse)

newPointer=function(inputValue){ 
  object=new.env(parent=globalenv()) 
  object$value=inputValue 
  class(object)='pointer'
  
  return(object) 
} 

updatePointerValue=function (object, ...) { # create S3 generic 
  UseMethod("updatePointerValue")  
}  

install.packages("pointr")

#source("../R/ht4sa.R")
setClass("EnsembleUnit", 
         representation(
           element = "character",
           indices = "vector", #column-row
           values="vector"))

set_ensemble_element<-function(input = list(), 
                               element = EnsembleUnit,
                               index){
  
  if(length(element@indices)== 1){
    input[[element@element]][[element@indices[1]]]<-element@values[index]
  }else if(length(element@indices)== 2){
    input[[element@element]][[element@indices[1]]][[element@indices[2]]]<-element@values[index]
    }
  
  return(input)
}

set_ensemble_element<-function(input = list(), 
                               element = EnsembleUnit,
                               value = numeric){
  
  if(length(element@indices)== 1){
    input[[element@element]][[element@indices[1]]]<-value
  }else if(length(element@indices)== 2){
    input[[element@element]][[element@indices[1]]][[element@indices[2]]]<-value
  }
  
  return(input)
}

build_parameter_sets<-function(count = as.integer(), 
                               current = as.integer(), 
                               working = vector,
                               source,
                               combos){


  if (current == 1) {
    for (i in 1:length(source[[1]])) {
      v=c(0)
      v[1] = source[[1]][i]
      combos<-build_parameter_sets(count, current + 1, v, source, combos)
    }
  } else if (current < length(source)) {
    v = working;
    append(v,0)
    for (i in 1:length(source[[current]])) {
      v[length(v)] = source[[current]][i]
      combos <-build_parameter_sets(count, current + 1, v, source, combos)
    }
  } else {
    v = working;
    append(v,0)
    for (i in 1:length(source[[current]])) {
      v[length(v)] = source[[current]][i]
      print(source[current][i])
      combos[[length(combos)+1]]<-v
      count<-count+1
    }
    
  }
  return(combos)
}

source<-list(c(1,2,3), c(4,5,6))
source[[length(source)+1]]<-c(0.1,.3)
source
combos<-list()
working<-rep(0,2)
working
cc<-build_parameter_sets(1,1,working,source,combos)
print(cc)
class(combos)

c<-expand.grid(combos)#choose(combos,2)
print(c)
ensemble_units<-list()

e<-new("EnsembleUnit", 
       element = "SR_parms", 
       indices = c(2,4),
       values=c(1.1,2.0,3.4))



ensemble_units<-append(ensemble_units, e)
print(ensemble_units)

e2<-new("EnsembleUnit", 
       element = "SR_parms", 
       indices = c(1,4),
       values=c(1.1,2.0,3.4))
ensemble_units<-append(ensemble_units, e2)
print(ensemble_units)
#_____________________________________________________________________________________________________________________________
# set directory paths
proj_dir = "/Users/mattadmin/ht4sa/test/"
dir_utility = "/Users/mattadmin/SAP/Code/Utilities/"
dir_model = paste0(proj_dir,"test_case/")
dir_model_source = "input_files"
model_name_stem = "input_files/"
dir_input_files = paste0(dir_model,dir_model_source,"/")
dir_ss = paste0(dir_model,"input_files/")
dir.create(paste0(proj_dir,"Summary/",model_name_stem),recursive=TRUE,showWarnings=FALSE)

print(dir_input_files)

tmp_ctl = SS_readctl(
  file = paste0(dir_ss, "control.ss"),
  use_datlist = TRUE,
  datlist = paste0(dir_ss, "data.ss")
)
print(tmp_ctl)

print(tmp_ctl["SR_parms"])

x <- tmp_ctl["SR_parms"] %>% slice(4)
tmp_ctl[["SR_parms"]][[2]][[4]]
tmp_ctl[[e@element]][[e@indices[1]]][[e@indices[2]]]

tmp_ctl<-set_ensemble_element(tmp_ctl, e, 1)
tmp_ctl[[e@element]][[e@indices[1]]][[e@indices[2]]]

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
