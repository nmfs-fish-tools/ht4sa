library(Rcpp)
library(ht4sa)

ht4sa <- Rcpp::Module("ht4sa", PACKAGE= "ht4sa")

ensemble<-new(ht4sa$ht4sa_ensemble)

recruitment_models<-list()
selectivity_models<-list()

beverton_holt_1<-new(ht4sa$beverton_holt)
ensemble$add_recruitment_unit(beverton_holt_1$get_id())
recruitment_models[[beverton_holt_1$get_id()]]<-beverton_holt_1


beverton_holt_2<-new(ht4sa$beverton_holt)
ensemble$add_recruitment_unit(beverton_holt_2$get_id())
recruitment_models[[beverton_holt_2$get_id()]]<-beverton_holt_2

beverton_holt_3<-new(ht4sa$beverton_holt)
ensemble$add_recruitment_unit(beverton_holt_3$get_id())
recruitment_models[[beverton_holt_3$get_id()]]<-beverton_holt_3

ricker_1<-new(ht4sa$ricker)
ensemble$add_recruitment_unit(ricker_1$get_id())
recruitment_models[[ricker_1$get_id()]]<-ricker_1

ricker_2<-new(ht4sa$ricker)
ensemble$add_recruitment_unit(ricker_2$get_id())
recruitment_models[[ricker_2$get_id()]]<-ricker_2

logistic_selectivity_1<-new(ht4sa$logistic_selectivity)
ensemble$add_selectivity_unit(logistic_selectivity_1$get_id())
selectivity_models[[logistic_selectivity_1$get_id()]]<-logistic_selectivity_1

logistic_selectivity_2<-new(ht4sa$logistic_selectivity)
ensemble$add_selectivity_unit(logistic_selectivity_2$get_id())
selectivity_models[[logistic_selectivity_2$get_id()]]<-logistic_selectivity_2

double_logistic_selectivity_1<-new(ht4sa$double_logistic_selectivity)
ensemble$add_selectivity_unit(double_logistic_selectivity_1$get_id())
selectivity_models[[double_logistic_selectivity_1$get_id()]]<-double_logistic_selectivity_1

double_normal_selectivity_1<-new(ht4sa$double_normal_selectivity)
ensemble$add_selectivity_unit(double_normal_selectivity_1$get_id())
selectivity_models[[double_normal_selectivity_1$get_id()]]<-double_normal_selectivity_1


models<-expand.grid(recruitment = 
                      ensemble$get_recruitment_units(),
                    selectivity = ensemble$get_selectivity_units())


cat("recruitment models\n")
for(i in 1:length(recruitment_models)){
  cat("id = ")
  cat(recruitment_models[[i]]$get_id())
  cat(", ss_type_id = ")
  cat(recruitment_models[[i]]$ss_id)
  cat(", name = ")
  cat(recruitment_models[[i]]$name)
  cat(", category = ")
  cat(recruitment_models[[i]]$category)
  cat("\n")
}
cat("\nselectivity models\n")
for(i in 1:length(selectivity_models)){
  cat("id = ")
  cat(selectivity_models[[i]]$get_id())
  cat(", ss_type_id = ")
  cat(selectivity_models[[i]]$ss_id)
  cat(", name = ")
  cat(selectivity_models[[i]]$name)
  cat(", category = ")
  cat(selectivity_models[[i]]$category)
  cat("\n")
}

cat("Combinations: \n")
print(models)
