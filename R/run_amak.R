#' @export
run_amak <- function(maindir=NULL, subdir="AMAK", om_sim_num=NULL, casedir=casedir, input_filename=NULL){
  setwd(file.path(casedir, "output", subdir))
  unlink(list.files(file.path(casedir, "output", "AMAK"), full.names = TRUE), recursive = TRUE)
  sapply(1:om_sim_num, function(x) dir.create(file.path(casedir, "output", subdir, paste("s", x, sep=""))))

  file.copy(file.path(maindir, "em_input", paste("amak_", input_filename, ".dat", sep="")), file.path(maindir, "em_input", "amak.dat"), overwrite = T)
  file.copy(file.path(maindir, "em_input", paste("amak_data_", input_filename, ".dat", sep="")), file.path(maindir, "em_input", "amak_data.dat"), overwrite = T)

  ctlf <- file.path(maindir, "em_input", "amak.dat")
  datf <- file.path(maindir, "em_input", "amak_data.dat")

  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    load(file=file.path(casedir, "output", "OM", paste("OM", om_sim, ".RData", sep="")))

    if(modify_input == "all") {

    }

    if(modify_input == "partial") {
      char.lines <- readLines(ctlf)
      char.lines[grep("#SigmaR", char.lines)+1] <- gsub(gsub(" .*$", "", char.lines[grep("#SigmaR", char.lines)+1]), om_input$logR_sd, char.lines[grep("#SigmaR", char.lines)+1]) #Replace the value before 1st space
      #char.lines[grep("#catchability", char.lines)+1] <- gsub(gsub(" .*$", "", char.lines[grep("#catchability", char.lines)+1]), em_input$survey_q$survey1, char.lines[grep("#catchability", char.lines)+1])
      writeLines(char.lines, con=file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.dat"))

      char.lines <- readLines(datf)
      char.lines[grep("#styr\t", char.lines)+1] <- as.character(paste(om_output$year[1], collapse="\t"))
      char.lines[grep("#endyr\t", char.lines)+1] <- as.character(paste(om_output$year[length(om_output$year)], collapse="\t"))

      char.lines[grep("#catch\t", char.lines)+1] <- as.character(paste(em_input$L.obs$fleet1, collapse="\t"))
      char.lines[grep("#catch_cv", char.lines)+1] <- as.character(paste(rep(em_input$cv.L$fleet1, length(em_input$L.obs$fleet1)), collapse="\t"))
      char.lines[grep("#sample_ages_fsh", char.lines)+1] <- as.character(paste(rep(em_input$n.L$fleet1, length(em_input$L.obs$fleet1)), collapse="\t"))
      for (i in 1:nrow(em_input$L.age.obs$fleet1)){
        char.lines[grep("#page_fsh", char.lines)+i]<-as.character(paste(em_input$L.age.obs$fleet1[i,],collapse="\t"))
      }

      for (i in 1:nrow(em_input$L.age.obs$fleet1)){
      char.lines[grep("#wt_age_fsh", char.lines)+i] <- as.character(paste(om_input$W.mt, collapse="\t"))
      }

      for (survey_id in 1:om_input$survey_num){
        temp <- em_input$survey.obs[[survey_id]]
        char.lines[grep("#biom_ind", char.lines)+survey_id]<-as.character(paste(temp, collapse="\t"))

        temp <- em_input$cv.survey[[survey_id]]*em_input$survey.obs[[survey_id]]
        char.lines[grep("#biom_cv", char.lines)+survey_id]<-as.character(paste(temp, collapse="\t"))

        char.lines[grep("#sample_ages_ind", char.lines)+survey_id] <- as.character(paste(rep(em_input$n.survey[[survey_id]], length(em_input$survey.obs[[survey_id]])), collapse="\t"))
      }

      temp <- do.call(rbind, em_input$survey.age.obs)
      for (i in 1:nrow(temp)){
        char.lines[grep("#page_ind", char.lines)+i]<-as.character(paste(temp[i,],collapse="\t"))
      }

      for (i in 1:nrow(temp)){
      char.lines[grep("#wt_age_ind", char.lines)+i] <- as.character(paste(rep(1, om_input$nages), collapse="\t"))
      }

      writeLines(char.lines, con=file.path(casedir, "output", subdir, paste("s", om_sim, sep=""),  "amak_data.dat"))
    }
  }

  cl <- ifelse(detectCores()==1, detectCores(), detectCores()-2)
  registerDoParallel(cl)

  foreach (om_sim = 1:om_sim_num) %dopar% {
    setwd(file.path(casedir, "output", subdir, paste("s", om_sim, sep="")))
    file.copy(file.path(maindir, "em_input", "amak.exe"), file.path(casedir,"output", subdir, paste("s", om_sim, sep=""), "amak.exe"), overwrite = T)
    system(paste(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.exe"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.dat"), sep = " "), show.output.on.console = FALSE)
    file.remove(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.exe"))
    file_list <- list.files(path = getwd())
    file.remove(c(file_list[!(file_list %in% c(list.files(path = getwd(), pattern = c(".rep")),
                                               list.files(path = getwd(), pattern = c(".par")),
                                               list.files(path = getwd(), pattern = c(".std")),
                                               list.files(path = getwd(), pattern = c(".rdat")),
                                               list.files(path = getwd(), pattern = c(".cov")),
                                               list.files(path = getwd(), pattern = c(".dat"))))]))
  }
  #stopCluster(cl)

  # for (om_sim in 1:om_sim_num){
  #   setwd(file.path(casedir, "output", subdir, paste("s", om_sim, sep="")))
  #   file.copy(file.path(maindir, "em_input", "amak.exe"), file.path(casedir,"output", subdir, paste("s", om_sim, sep=""), "amak.exe"), overwrite = T)
  #   system(paste(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.exe"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.dat"), sep = " "), show.output.on.console = FALSE)
  #   file.remove(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "amak.exe"))
  # }
}


