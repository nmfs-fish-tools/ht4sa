#' @export
run_ss <- function(maindir=maindir, subdir="SS", om_sim_num=NULL, casedir=casedir, em_bias_cor=em_bias_cor, input_filename=NULL, initial_equilibrium_F=TRUE){
  if(!("r4ss" %in% installed.packages()[,"Package"])) install.packages("r4ss")
  library(r4ss)

  setwd(file.path(casedir, "output", subdir))
  file.copy(file.path(maindir, "em_input", paste("control_", input_filename, ".ss", sep="")), file.path(maindir, "em_input", "control.ss"), overwrite = T)
  file.copy(file.path(maindir, "em_input", paste("data_", input_filename, ".ss", sep="")), file.path(maindir, "em_input", "data.ss"), overwrite = T)
  file.copy(file.path(maindir, "em_input", paste("wtatage_", input_filename, ".ss", sep="")), file.path(maindir, "em_input", "wtatage.ss"), overwrite = T)

  unlink(list.files(file.path(casedir, "output", "SS"), full.names = TRUE), recursive = TRUE)
  sapply(1:om_sim_num, function(x) dir.create(file.path(casedir, "output", subdir, paste("s", x, sep=""))))

  ss_data <- SS_readdat_3.30(file.path(maindir, "em_input", "data.ss"), verbose=FALSE)

  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    load(file=file.path(casedir,"output", "OM", paste("OM", om_sim, ".RData", sep="")))
    #ss_ctl <- SS_readctl_3.30(file.path(maindir, "em_input", "control.ss"), verbose=FALSE, nseas=1, N_areas=1, Nages=om_input$nages, Ngenders=1, Npopbins=NA, Nfleet=om_input$fleet_num, Nsurveys=om_input$survey_num, Do_AgeKey=FALSE, N_tag_groups=NA, N_CPUE_obs=c(0,0,9,12))

    if(modify_input == "all") {

    } #incomplement
    if(modify_input == "partial") {
      if (initial_equilibrium_F==FALSE){
        ss_data$catch <- as.data.frame(cbind(om_input$year,
                                             rep(1, times=om_input$nyr),
                                             rep(1, times=om_input$nyr),
                                             em_input$L.obs$fleet1,
                                             rep(sqrt(log(1+em_input$cv.L$fleet1^2)), times=length(em_input$L.obs$fleet1))))
      } else {
        ss_data$catch <- as.data.frame(rbind(c(-999, 1, 1, em_input$L.obs$fleet1[1]/om_input$nages, sqrt(log(1+em_input$cv.L$fleet1^2))),
                                             cbind(om_input$year,
                                                   rep(1, times=om_input$nyr),
                                                   rep(1, times=om_input$nyr),
                                                   em_input$L.obs$fleet1,
                                                   rep(sqrt(log(1+em_input$cv.L$fleet1^2)), times=length(em_input$L.obs$fleet1))))
        )
      }


      names(ss_data$catch) <- c("year", "seas", "fleet", "catch", "catch_se")

      if (om_input$survey_num == 1) {
       ss_data$CPUE <- as.data.frame(cbind(om_input$year, rep(1, times=om_input$nyr), rep(2, times=om_input$nyr), em_input$survey.obs$survey1, rep(sqrt(log(1+em_input$cv.survey$survey1^2)), times=length(em_input$survey.obs$survey1))))
       names(ss_data$CPUE) <- c("year", "seas", "index", "obs", "se_log")
       ss_data$agecomp[, 9:ncol(ss_data$agecomp)] <- rbind(cbind(rep(em_input$n.L$fleet1, length(em_input$L.obs$fleet1)), em_input$L.age.obs$fleet1), cbind(rep(em_input$n.survey$survey1, length(em_input$survey.obs$survey1)), em_input$survey.age.obs$survey1))
      }
      if (om_input$survey_num == 2) {
        ss_data$CPUE <- as.data.frame(rbind(
        cbind(om_input$year, rep(1, times=om_input$nyr), rep(2, times=om_input$nyr), em_input$survey.obs$survey1, rep(sqrt(log(1+em_input$cv.survey$survey1^2)), times=length(em_input$survey.obs$survey1))),
        cbind(om_input$year, rep(1, times=om_input$nyr), rep(3, times=om_input$nyr), em_input$survey.obs$survey2, rep(sqrt(log(1+em_input$cv.survey$survey2^2)), times=length(em_input$survey.obs$survey2)))))
       names(ss_data$CPUE) <- c("year", "seas", "index", "obs", "se_log")
       ss_data$agecomp[, 9:ncol(ss_data$agecomp)] <- rbind(
        cbind(rep(em_input$n.L$fleet1, length(em_input$L.obs$fleet1)), em_input$L.age.obs$fleet1),
        cbind(rep(em_input$n.survey$survey1, length(em_input$survey.obs$survey1)), em_input$survey.age.obs$survey1),
        cbind(rep(em_input$n.survey$survey2, length(em_input$survey.obs$survey2)), em_input$survey.age.obs$survey2))
      }


      # ss_ctl$Q_parms$INIT <- log(em_input$survey_q$survey1*1000)
      # ss_ctl$Q_parms$PRIOR <- log(em_input$survey_q$survey1*1000)
      # ss_ctl$SRparm$INIT[which(rownames(ss_ctl$SRparm)=="SR_sigmaR")] <- om_input$logR_sd
      # ss_ctl$SRparm$PRIOR[which(rownames(ss_ctl$SRparm)=="SR_sigmaR")] <- om_input$logR_sd
      # ss_ctl$F_setup[1] <- om_output$f[1]
      # if(em_bias_cor==FALSE) {
      #   ss_ctl$max_bias_adj <- 0
      # }else{
      #   ss_ctl$max_bias_adj <- 1
      # }

    }
    SS_writedat_3.30(datlist = ss_data, outfile = file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "data.ss"), overwrite = T, verbose = FALSE)
    #SS_writectl_3.30(ctllist = ss_ctl, outfile = file.path(maindir, "output", subdir, paste("s", om_sim, sep=""), "control.ss"), overwrite = T, verbose = F)

    file.copy(file.path(maindir, "em_input", "control.ss"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "control.ss"), overwrite = T)
    file.copy(file.path(maindir, "em_input", "starter.ss"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "starter.ss"), overwrite = T)
    file.copy(file.path(maindir, "em_input", "forecast.ss"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "forecast.ss"), overwrite = T)
    file.copy(file.path(maindir, "em_input", "wtatage.ss"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "wtatage.ss"), overwrite = T)
  }

  cl <- ifelse(detectCores()==1, detectCores(), detectCores()-2)
  registerDoParallel(cl)

  foreach (om_sim = 1:om_sim_num) %dopar% {
    setwd(file.path(casedir, "output", subdir, paste("s", om_sim, sep="")))
    file.copy(file.path(maindir, "em_input", "ss.exe"), file.path(casedir,"output", subdir, paste("s", om_sim, sep=""), "ss.exe"), overwrite = T)
    system(paste(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "ss.exe"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "data.ss"), sep = " "), show.output.on.console = FALSE)
    file.remove(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "ss.exe"))
    # file_list <- list.files(path = getwd())
    # file.remove(c(file_list[!(file_list %in% c(list.files(path = getwd(), pattern = c(".rep")),
    #                                            list.files(path = getwd(), pattern = c(".par")),
    #                                            list.files(path = getwd(), pattern = c(".std")),
    #                                            list.files(path = getwd(), pattern = c(".rdat")),
    #                                            list.files(path = getwd(), pattern = c(".cov")),
    #                                            list.files(path = getwd(), pattern = c("Report.sso")),
    #                                            list.files(path = getwd(), pattern = c("ss.cor")),
    #                                            list.files(path = getwd(), pattern = c(".ss"))))]))
  }
  #stopCluster(cl)

  # for (om_sim in 1:om_sim_num){
  #   # setwd(file.path(maindir, "output", subdir, paste("s", om_sim, sep="")))
  #   # system(paste(file.path(maindir, "em_input", "ss.exe"), " data.ss", sep = ""), show.output.on.console = FALSE)
  #
  #   setwd(file.path(casedir, "output", subdir, paste("s", om_sim, sep="")))
  #   file.copy(file.path(maindir, "em_input", "ss.exe"), file.path(casedir,"output", subdir, paste("s", om_sim, sep=""), "ss.exe"), overwrite = T)
  #   system(paste(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "ss.exe"), file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "data.ss"), sep = " "), show.output.on.console = FALSE)
  #   file.remove(file.path(casedir, "output", subdir, paste("s", om_sim, sep=""), "ss.exe"))
  # }
}
