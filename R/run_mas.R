#' Function to run Metapopulation Assessment System (https://github.com/nmfs-fish-tools/r4MAS)
#' @name run_mas
#' @description Function to read simulated true values from the operating model and to run MAS
#' @param maindir Main working directory
#' @param subdir Estimation model working directory
#' @param om_sim_num Number of iterations from the operating model
#' @param casedir Case working directory
#' @param em_bias_cor Use bias correction in the estimation model or not?
#' @export
run_mas <- function(
  maindir = maindir, 
  subdir = "MAS", 
  om_sim_num = NULL, 
  casedir = casedir, 
  em_bias_cor = em_bias_cor, 
  input_filename = NULL) {
  
  if (!("r4MAS" %in% installed.packages()[, "Package"])) stop("Please install r4MAS!")

  setwd(file.path(casedir, "output", subdir))
  unlink(list.files(file.path(casedir, "output", "MAS"), full.names = TRUE), recursive = TRUE)
  sapply(1:om_sim_num, function(x) dir.create(file.path(casedir, "output", subdir, paste("s", x, sep = ""))))

  for (om_sim in 1:om_sim_num) {
    load(file=file.path(casedir, "output", "OM", paste("OM", om_sim, ".RData", sep="")))

    library(r4MAS)
    r4mas <- Rcpp::Module("rmas", PACKAGE = "r4MAS")

    nyears <- om_input$nyr
    nseasons <- 1
    nages <- om_input$nages
    ages <- om_input$ages
    area1 <- new(r4mas$Area)
    area1$name <- "area1"
    
    recruitment <- new(r4mas$BevertonHoltRecruitment)
    recruitment$R0$value <- om_input$R0 / 1000
    recruitment$R0$estimated <- TRUE
    recruitment$R0$phase <- 1
    recruitment$h$value <- om_input$h
    recruitment$h$estimated <- FALSE
    recruitment$h$phase <- 3
    recruitment$h$min <- 0.2001
    recruitment$h$max <- 1.0
    recruitment$sigma_r$value <- om_input$logR_sd
    recruitment$sigma_r$estimated <- FALSE
    recruitment$sigma_r$min <- 0
    recruitment$sigma_r$max <- 1.0
    recruitment$sigma_r$phase <- 2
    recruitment$estimate_deviations <- TRUE
    recruitment$constrained_deviations <- TRUE
    recruitment$deviations_min <- -15.0
    recruitment$deviations_max <- 15.0
    recruitment$deviation_phase <- 2
    recruitment$SetDeviations(rep(0.0, times = nyears))
    recruitment$use_bias_correction <- em_bias_cor
    
    growth <- new(r4mas$VonBertalanffyModified)
    empirical_weight <- rep(om_input$W.kg, times = om_input$nyr)
    survey_empirical_weight <- replicate(nages * nyears, 1.0)
    growth$SetUndifferentiatedCatchWeight(empirical_weight)
    growth$SetUndifferentiatedWeightAtSeasonStart(empirical_weight)
    growth$SetUndifferentiatedWeightAtSpawning(empirical_weight)
    growth$SetUndifferentiatedSurveyWeight(survey_empirical_weight)
    
    maturity <- new(r4mas$Maturity)
    maturity$values <- om_input$mat.age 
    
    natural_mortality <- new(r4mas$NaturalMortality)
    natural_mortality$SetValues(om_input$M.age)
    
    # Only 1 area in this model
    movement <- new(r4mas$Movement)
    movement$connectivity_females <- c(0.0)
    movement$connectivity_males <- c(0.0)
    movement$connectivity_recruits <- c(0.0)
    
    initial_deviations <- new(r4mas$InitialDeviations)
    initial_deviations$values <- rep(0.0, times = om_input$nages)
    initial_deviations$estimate <- TRUE
    initial_deviations$phase <- 2
    
    population <- new(r4mas$Population)
    for (y in 1:(nyears))
    {
      population$AddMovement(movement$id, y)
    }
    population$AddNaturalMortality(natural_mortality$id, area1$id, "undifferentiated")
    population$AddMaturity(maturity$id, area1$id, "undifferentiated")
    population$AddRecruitment(recruitment$id, 1, area1$id)
    population$SetInitialDeviations(initial_deviations$id, area1$id, "undifferentiated")
    population$SetGrowth(growth$id)
    population$sex_ratio <- 0.5
    
    # Catch index values and observation errors
    catch_index <- new(r4mas$IndexData)
    catch_index$values <- em_input$L.obs$fleet1
    catch_index$error <- rep(em_input$cv.L$fleet1, times = om_input$nyr)
    # Catch composition data
    catch_comp <- new(r4mas$AgeCompData)
    catch_comp$values <- as.vector(t(em_input$L.age.obs$fleet1))
    catch_comp$sample_size <- rep(em_input$n.L$fleet1, nyears * nseasons)
    # Likelihood component settings
    fleet_index_comp_nll <- new(r4mas$Lognormal)
    fleet_index_comp_nll$use_bias_correction <- FALSE
    fleet_age_comp_nll <- new(r4mas$Multinomial)
    # Fleet selectivity settings
    if (om_input$sel_fleet$fleet1$pattern == 1) {
      fleet_selectivity <- new(r4mas$LogisticSelectivity)
      fleet_selectivity$a50$value <- om_input$sel_fleet$fleet1$A50.sel
      fleet_selectivity$a50$estimated <- TRUE
      fleet_selectivity$a50$phase <- 2
      fleet_selectivity$a50$min <- 0.0
      fleet_selectivity$a50$max <- max(om_input$ages)
      fleet_selectivity$slope$value <- 1 / om_input$sel_fleet$fleet1$slope.sel
      fleet_selectivity$slope$estimated <- TRUE
      fleet_selectivity$slope$phase <- 2
      fleet_selectivity$slope$min <- 0.0001
      fleet_selectivity$slope$max <- 5
    }
    
    if (om_input$sel_fleet$fleet1$pattern == 2) {
      fleet_selectivity <- new(r4mas$DoubleLogisticSelectivity)
      
      fleet_selectivity$alpha_asc$value <- om_input$sel_fleet$fleet1$A50.sel1
      fleet_selectivity$alpha_asc$estimated <- TRUE
      fleet_selectivity$alpha_asc$phase <- 2
      fleet_selectivity$alpha_asc$min <- 0.0
      fleet_selectivity$alpha_asc$max <- max(om_input$ages)
      
      fleet_selectivity$beta_asc$value <- om_input$sel_fleet$fleet1$slope.sel1
      fleet_selectivity$beta_asc$estimated <- TRUE
      fleet_selectivity$beta_asc$phase <- 2
      fleet_selectivity$beta_asc$min <- 0.0001
      fleet_selectivity$beta_asc$max <- max(om_input$ages)
      
      fleet_selectivity$alpha_desc$value <- om_input$sel_fleet$fleet1$A50.sel2
      fleet_selectivity$alpha_desc$estimated <- TRUE
      fleet_selectivity$alpha_desc$phase <- 2
      fleet_selectivity$alpha_desc$min <- 0.0
      fleet_selectivity$alpha_desc$max <- max(om_input$ages)
      
      fleet_selectivity$beta_desc$value <- om_input$sel_fleet$fleet1$slope.sel2
      fleet_selectivity$beta_desc$estimated <- TRUE
      fleet_selectivity$beta_desc$phase <- 2
      fleet_selectivity$beta_desc$min <- 0.0001
      fleet_selectivity$beta_desc$max <- max(om_input$ages)
      
    }
    
    
    # Fishing mortality settings
    fishing_mortality <- new(r4mas$FishingMortality)
    fishing_mortality$estimate <- TRUE
    fishing_mortality$phase <- 1
    fishing_mortality$min <- 0.0
    fishing_mortality$max <- 4
    fishing_mortality$SetValues(om_output$f)
    # Create the fleet
    fleet <- new(r4mas$Fleet)
    fleet$AddIndexData(catch_index$id, "undifferentiated")
    fleet$AddAgeCompData(catch_comp$id, "undifferentiated")
    fleet$SetIndexNllComponent(fleet_index_comp_nll$id)
    fleet$SetAgeCompNllComponent(fleet_age_comp_nll$id)
    fleet$AddSelectivity(fleet_selectivity$id, 1, area1$id)
    fleet$AddFishingMortality(fishing_mortality$id, 1, area1$id)
    
    # Survey index values and observation errors
    survey_index <- vector(mode = "list", length = om_input$survey_num)
    
    for (i in 1:om_input$survey_num) {
      survey_index[[i]] <- new(r4mas$IndexData)
      # survey_index[[i]]$id <- i
      survey_index[[i]]$values <- em_input$survey.obs[[i]]
      survey_index[[i]]$error <- rep(em_input$cv.survey[[i]], times = om_input$nyr)
    }
    
    # Survey composition
    survey_comp <- vector(mode = "list", length = om_input$survey_num)
    for (i in 1:om_input$survey_num) {
      survey_comp[[i]] <- new(r4mas$AgeCompData)
      # survey_comp[[i]]$id <- i
      survey_comp[[i]]$values <- as.vector(t(em_input$survey.age.obs[[i]]))
      survey_comp[[i]]$sample_size <- rep(em_input$n.survey[[i]], times = om_input$nyr)
    }
    
    # Likelihood component settings
    
    survey_index_comp_nll <- vector(mode = "list", length = om_input$survey_num)
    survey_age_comp_nll <- vector(mode = "list", length = om_input$survey_num)
    for (i in 1:om_input$survey_num) {
      survey_index_comp_nll[[i]] <- new(r4mas$Lognormal)
      survey_index_comp_nll[[i]]$use_bias_correction <- FALSE
      # survey_index_comp_nll[[i]]$id <- i
      
      survey_age_comp_nll[[i]] <- new(r4mas$Multinomial)
      # survey_age_comp_nll[[i]]$id <- i
    }
    
    # Survey selectivity settings
    survey_selectivity <- vector(mode = "list", length = om_input$survey_num)
    
    for (i in 1:om_input$survey_num) {
      selectivity_option <- om_input$sel_survey[[i]]$pattern
      
      if (selectivity_option == 1) {
        survey_selectivity[[i]] <- new(r4mas$LogisticSelectivity)
        # survey_selectivity[[i]]$id <- i
        survey_selectivity[[i]]$a50$value <- om_input$sel_survey[[i]]$A50.sel
        survey_selectivity[[i]]$a50$estimated <- TRUE
        survey_selectivity[[i]]$a50$phase <- 2
        survey_selectivity[[i]]$a50$min <- 0
        survey_selectivity[[i]]$a50$max <- max(om_input$ages)
        
        survey_selectivity[[i]]$slope$value <- 1 / om_input$sel_survey[[i]]$slope.sel
        survey_selectivity[[i]]$slope$estimated <- TRUE
        survey_selectivity[[i]]$slope$phase <- 2
        survey_selectivity[[i]]$slope$min <- 0.0001
        survey_selectivity[[i]]$slope$max <- max(om_input$ages)
      }
      
      if (selectivity_option == 2) {
        survey_selectivity[[i]] <- new(r4mas$DoubleLogisticSelectivity)
        # survey_selectivity[[i]]$id <- i
        
        survey_selectivity[[i]]$alpha_asc$value <- om_input$sel_survey[[i]]$A50.sel1
        survey_selectivity[[i]]$alpha_asc$estimated <- TRUE
        survey_selectivity[[i]]$alpha_asc$phase <- 2
        survey_selectivity[[i]]$alpha_asc$min <- 0.0
        survey_selectivity[[i]]$alpha_asc$max <- max(om_input$ages)
        
        survey_selectivity[[i]]$beta_asc$value <- om_input$sel_survey[[i]]$slope.sel1
        survey_selectivity[[i]]$beta_asc$estimated <- TRUE
        survey_selectivity[[i]]$beta_asc$phase <- 2
        survey_selectivity[[i]]$beta_asc$min <- 0.0001
        survey_selectivity[[i]]$beta_asc$max <- max(om_input$ages)
        
        survey_selectivity[[i]]$alpha_desc$value <- om_input$sel_survey[[i]]$A50.sel2
        survey_selectivity[[i]]$alpha_desc$estimated <- TRUE
        survey_selectivity[[i]]$alpha_desc$phase <- 2
        survey_selectivity[[i]]$alpha_desc$min <- 0.0
        survey_selectivity[[i]]$alpha_desc$max <- max(om_input$ages)
        
        survey_selectivity[[i]]$beta_desc$value <- om_input$sel_survey[[i]]$slope.sel2
        survey_selectivity[[i]]$beta_desc$estimated <- TRUE
        survey_selectivity[[i]]$beta_desc$phase <- 2
        survey_selectivity[[i]]$beta_desc$min <- 0.0001
        survey_selectivity[[i]]$beta_desc$max <- max(om_input$ages)
        
      }
      
    }
    
    # Create the survey
    survey <- vector(mode = "list", length = om_input$survey_num)
    for (i in 1:om_input$survey_num) {
      survey[[i]] <- new(r4mas$Survey)
      
      survey[[i]]$AddIndexData(survey_index[[i]]$id, "undifferentiated")
      survey[[i]]$AddAgeCompData(survey_comp[[i]]$id, "undifferentiated")
      survey[[i]]$SetIndexNllComponent(survey_index_comp_nll[[i]]$id)
      survey[[i]]$SetAgeCompNllComponent(survey_age_comp_nll[[i]]$id)
      survey[[i]]$AddSelectivity(survey_selectivity[[i]]$id, 1, area1$id)
      
      survey[[i]]$q$value <- em_input$survey_q[[i]]
      survey[[i]]$q$min <- 0
      survey[[i]]$q$max <- 10
      survey[[i]]$q$estimated <- TRUE
      survey[[i]]$q$phase <- 1
    }
    
    mas_model <- new(r4mas$MASModel)
    mas_model$compute_variance_for_derived_quantities<-FALSE
    mas_model$nyears <- nyears
    mas_model$nseasons <- nseasons
    mas_model$nages <- nages
    mas_model$extended_plus_group <- max(om_input$ages)
    mas_model$ages <- ages
    mas_model$catch_season_offset <- 0.0
    mas_model$spawning_season_offset <- 0.0
    mas_model$survey_season_offset <- 0.0
    mas_model$AddPopulation(population$id)
    mas_model$AddFleet(fleet$id)
    for (i in 1:om_input$survey_num) {
      mas_model$AddSurvey(survey[[i]]$id)
    }
    # mas_model$tolerance <- 0.0001
    mas_model$max_iterations <- 10000

    # Run MAS
    mas_model$Run()
    output_file <- file.path(casedir, "output", subdir, paste("s", om_sim, sep = ""), paste("s", om_sim, ".json", sep = ""))
    write(mas_model$GetOutput(), file = toString(output_file))
    mas_model$Reset()
  }
}
