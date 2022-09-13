library(Rcpp)
library(ht4sa)

#load the ht4sa module
ht4sa <- Rcpp::Module("ht4sa", PACKAGE = "ht4sa")

#create the ensemble object
ensemble <- new(ht4sa$ht4sa_ensemble)

#create an interface object
beverton_holt_1 <- new(ht4sa$beverton_holt)
#add it to ht4sa recruitment list
ht4sa_recruitment_models[[beverton_holt_1$get_id()]] <-
  beverton_holt_1


beverton_holt_2 <- new(ht4sa$beverton_holt)
ht4sa_recruitment_models[[beverton_holt_2$get_id()]] <-
  beverton_holt_2

beverton_holt_3 <- new(ht4sa$beverton_holt)
ht4sa_recruitment_models[[beverton_holt_3$get_id()]] <-
  beverton_holt_3

ricker_1 <- new(ht4sa$ricker)
ht4sa_recruitment_models[[ricker_1$get_id()]] <- ricker_1

ricker_2 <- new(ht4sa$ricker)
ht4sa_recruitment_models[[ricker_2$get_id()]] <- ricker_2

logistic_selectivity_1 <- new(ht4sa$logistic_selectivity)
ht4sa_selectivity_models[[logistic_selectivity_1$get_id()]] <-
  logistic_selectivity_1

logistic_selectivity_2 <- new(ht4sa$logistic_selectivity)
ht4sa_selectivity_models[[logistic_selectivity_2$get_id()]] <-
  logistic_selectivity_2

double_logistic_selectivity_1 <-
  new(ht4sa$double_logistic_selectivity)
ht4sa_selectivity_models[[double_logistic_selectivity_1$get_id()]] <-
  double_logistic_selectivity_1

double_normal_selectivity_1 <- new(ht4sa$double_normal_selectivity)
ht4sa_selectivity_models[[double_normal_selectivity_1$get_id()]] <-
  double_normal_selectivity_1

von_bertalanffy<-new(ht4sa$von_bertalanffy)
ht4sa_growth_models[[von_bertalanffy$get_id()]] <-von_bertalanffy

# schnute<-new(ht4sa$schnute)
# ht4sa_growth_models[[schnute$get_id()]] <-schnute

print(length(ht4sa_recruitment_models))
#extract ids and create combinations
models <- ht4sa_create_models()




cat("recruitment models\n")
for (i in 1:length(ht4sa_recruitment_models)) {
  cat("id = ")
  cat(ht4sa_recruitment_models[[i]]$get_id())
  cat(", ss_type_id = ")
  cat(ht4sa_recruitment_models[[i]]$ss_id)
  cat(", name = ")
  cat(ht4sa_recruitment_models[[i]]$name)
  cat(", category = ")
  cat(ht4sa_recruitment_models[[i]]$category)
  cat("\n")
}
cat("\nselectivity models\n")
for (i in 1:length(ht4sa_selectivity_models)) {
  cat("id = ")
  cat(ht4sa_selectivity_models[[i]]$get_id())
  cat(", ss_type_id = ")
  cat(ht4sa_selectivity_models[[i]]$ss_id)
  cat(", name = ")
  cat(ht4sa_selectivity_models[[i]]$name)
  cat(", category = ")
  cat(ht4sa_selectivity_models[[i]]$category)
  cat("\n")
}

cat("Combinations: \n")
print(models)

for (i in 1:nrow(models)) {
  names<-colnames(models) 

  cat("\nModel ")
  cat(i)
  cat(" of ")
  cat(nrow(models))
  cat(":\n")
  path<-paste(as.character(as.vector(models[i,])), collapse = "_")
  cat("Path to input: /User/matthew/ht4sa/ss/input/")
  cat(path)
  cat("\n")
 
  if ("recruitment" %in% colnames(models))
  {
    cat(ht4sa_recruitment_models[[models$recruitment[[i]]]]$name)
    cat("(id = ")
    cat(ht4sa_recruitment_models[[models$recruitment[[i]]]]$get_id())
    cat(")\n")
  }
  
  if ("selectivity" %in% colnames(models))
  {
    cat(ht4sa_selectivity_models[[models$selectivity[[i]]]]$name)
    cat("(id = ")
    cat(ht4sa_selectivity_models[[models$selectivity[[i]]]]$get_id())
    cat(")\n")
  }
  
  if ("growth" %in% colnames(models))
  {
    cat(ht4sa_growth_models[[models$growth[[i]]]]$name)
    cat("(id = ")
    cat(ht4sa_growth_models[[models$growth[[i]]]]$get_id())
    cat(")\n")
  }
}
