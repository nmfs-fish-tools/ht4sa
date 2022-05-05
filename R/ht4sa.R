




library(r4ss)
library(Rmpi)



# Load the R MPI package if it is not already loaded.
if (!is.loaded("mpi_initialize")) {
  library("Rmpi")
}

.Last <- function() {
  if (is.loaded("mpi_initialize")) {
    if (mpi.comm.size(1) > 0) {
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}



#run on mpi child
#' @export
run_ss_child <- function(testing_options_df,
                         begin,
                         end,
                         proj_dir,
                         dir_utility,
                         dir_model,
                         dir_model_source,
                         model_name_stem,
                         dir_input_files,
                         dir_ss,
                         TESTING = FALSE) {
  
  #hard code for now
  exe_path<-"/Users/mattadmin/ht4sa/bin"
  print(Sys.info()['sysname'])
  print(begin[mpi.comm.rank()])
  print(end[mpi.comm.rank()])
  
  #_____________________________________________________________________________________________________________________________
  # loop over options
  #TESTING = FALSE
  for (i in begin[mpi.comm.rank()]:end[mpi.comm.rank()])
  {
    print(i)
    #_____________________________________________________________________________________________________________________________
    # define directory structure
    model_name = paste0(testing_options_df[i,], collapse = "_")
    dir_run = paste0(dir_model, model_name_stem, model_name, "/")
    dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)
    dir_plot = paste0(proj_dir,
                      "Plot/Model_Runs/",
                      model_name_stem,
                      model_name,
                      "/")
    dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)
    print(dir_run)
    #_____________________________________________________________________________________________________________________________
    # transfer files
    # intial files
    FileList = c("starter.ss", "forecast.ss", "data.ss", "control.ss")
    file.copy(paste0(dir_input_files, FileList), dir_run, overwrite =
                TRUE)
    # ss files
    FileList = list.files(dir_ss)
    file.copy(paste0(dir_ss, FileList), dir_run, overwrite =
                TRUE)
    
    
    #_____________________________________________________________________________________________________________________________
    # modify input files
    # update the US HI longline selex patterns
    # make same for both 17 & 18
    tmp_ctl = SS_readctl(
      file = paste0(dir_run, "control.ss"),
      use_datlist = TRUE,
      datlist = paste0(dir_run, "data.ss")
    )
    
    # Update steepness
    tmp_ctl = SS_readctl(
      file = paste0(dir_run, "control.ss"),
      use_datlist = TRUE,
      datlist = paste0(dir_run, "data.ss")
    )
    # change SRR configuration from survival function to BH SRR
    tmp_ctl$SR_function = 3
    tmp_ctl$SR_parms$PR_type = 0
    tmp_ctl$SR_parms = tmp_ctl$SR_parms[-3, ]
    rownames(tmp_ctl$SR_parms)[2] = "SR_BH_steep"
    tmp_ctl$SR_parms$LO[2] = 0.2
    tmp_ctl$SR_parms$HI[2] = 1
    tmp_ctl$SR_parms$INIT[2] = testing_options_df$steepness[i]
    tmp_ctl$SR_parms$PRIOR[2] = testing_options_df$steepness[i]
    tmp_ctl$SR_parms$PR_SD[2] = 0.05
    #_____________________________________________________________________________________________________________________________
    # Update sigmaR
    tmp_ctl$SR_parms$INIT[3] = testing_options_df$sigmaR[i]
    #_____________________________________________________________________________________________________________________________
    # write-out files
    SS_writectl(
      tmp_ctl,
      paste0(dir_run, "control.ss"),
      version = "3.30",
      overwrite = TRUE
    )
    rm(list = c("tmp_ctl"))
    #_____________________________________________________________________________________________________________________________
    # (FOR TESTING ONLY) modify starter.ss
    if (TESTING)
    {
      tmp_starter = SS_readstarter(file = paste0(dir_run, "starter.ss"),
                                   verbose = TRUE)
      tmp_starter$last_estimation_phase = 0
      SS_writestarter(
        tmp_starter,
        dir = dir_run,
        file = "starter.ss",
        overwrite = TRUE,
        verbose = TRUE,
        warn = TRUE
      )
      rm(list = "tmp_starter")
    }
    
    #_____________________________________________________________________________________________________________________________
    # run SS
    if (TESTING)
    {
      system(paste0(
        "powershell cd ",
        dir_run,
        " ; ./ss_win.exe -nohess 3>&1 2>&1 > out.log"
      ))
      tmp = readLines(paste0(dir_run, "out.log"), skipNul = TRUE)
      if (length(grep("!!  Run has completed  !!", tmp, fixed = TRUE)) ==
          0) {
        rm(list = c("tmp"))
        stop("Error: Bad model settings.")
      }
    } else{

      # system(paste0("powershell cd ",dir_run," ; ./ss_win.exe"))
      if (print(Sys.info()['sysname']) == "Darwin") {
        system(paste0("cd ", dir_run))
        setwd(dir_run)
        system(paste0(exe_path,"/ss_osx"))
      } else if (print(Sys.info()['sysname']) == "Linux") {
        system(paste0("cd ", dir_run))
        setwd(dir_run)
        system(paste0(exe_path,"/ss_linux"))
      } else if (print(Sys.info()['sysname']) == "Windows") {

      }
    }
  }
  
  return (c(begin[mpi.comm.rank()], end[mpi.comm.rank()]))
}

#Run model in parallel
#' @export
run_ht4sa_ss_MPI <- function(testing_options_df,
                             proj_dir,
                             dir_utility,
                             dir_model,
                             dir_model_source,
                             model_name_stem,
                             dir_input_files,
                             dir_ss,
                             TESTING = FALSE) {
  id <- mpi.comm.rank(comm = 0)
  ns <-  mpi.universe.size() - 1
  
  
  nsims <- nrow(testing_options_df)
  if (ns > nsims) {
    ns <- nsims
  }
  begin <- rep(0, ns)
  end <- rep(0, ns)
  
  #create scenario segments
  print(paste("id = ", id))
  if (id == 0) {
    if (nsims < ns) {
      segments <- 1
      for (i in 1:nsims) {
        begin[i] <- as.integer((i - 1)  + 1)
        end[i] <- as.integer(i + 1)
      }
    } else{
      segments <- nsims / ns
      for (i in 1:ns) {
        if (i < ns) {
          begin[i] <- as.integer((i - 1) * segments + 1)
          end[i] <- as.integer(i * segments)
        } else{
          begin[i] <- as.integer((i - 1) * segments + 1)
          end[i] <- nsims
        }
      }
    }
    print("launching children....")
    mpi.spawn.Rslaves(nslaves = ns)
    
    mpi.bcast.Rfun2slave()
    #pass the function to all slaves
    
    # show(ht4sa::run_ss_child)
    mpi.bcast.Robj2slave(obj = r4ss::SS_readctl)
    mpi.bcast.Robj2slave(obj = SS_writectl)
    mpi.bcast.Robj2slave(run_ss_child)
    
    
    print("broadcasting run....")
    #execute all children
    x <- mpi.remote.exec(
      run_ss_child,
      testing_options_df,
      begin,
      end,
      proj_dir,
      dir_utility,
      dir_model,
      dir_model_source,
      model_name_stem,
      dir_input_files,
      dir_ss,
      TESTING
    )
    
    print("done broadcasting run....")
    print(x)
    # Tell all slaves to close down, and exit the program
    mpi.close.Rslaves(dellog = FALSE)
  }
  # save summary
  # save(
  #  selex_results_df,
  # file = paste0(
  #  proj_dir,
  # "Summary/",
  #model_name_stem,
  # "selex_results_df.RData"
  #)
  #)
  mpi.quit()
}


#run models sequentially
#' @export
run_ht4sa_ss_local <- function(testing_options_df,
                               proj_dir,
                               dir_utility,
                               dir_model,
                               dir_model_source,
                               model_name_stem,
                               dir_input_files,
                               dir_ss,
                               TESTING = FALSE) {
  # loop over options
  for (i in 1:nrow(testing_options_df))
  {
    #_____________________________________________________________________________________________________________________________
    # define directory structure
    model_name = paste0(testing_options_df[i,], collapse = "_")
    dir_run = paste0(dir_model, model_name_stem, model_name, "/")
    dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)
    dir_plot = paste0(proj_dir,
                      "Plot/Model_Runs/",
                      model_name_stem,
                      model_name,
                      "/")
    dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)
    print(dir_run)
    #_____________________________________________________________________________________________________________________________
    # transfer files
    # intial files
    FileList = c("starter.ss", "forecast.ss", "data.ss", "control.ss")
    file.copy(paste0(dir_input_files, FileList), dir_run, overwrite =
                TRUE)
    # ss files
    FileList = list.files(dir_ss)
    file.copy(paste0(dir_ss, FileList), dir_run, overwrite =
                TRUE)
    
    
    #_____________________________________________________________________________________________________________________________
    # modify input files
    # update the US HI longline selex patterns
    # make same for both 17 & 18
    tmp_ctl = SS_readctl(
      file = paste0(dir_run, "control.ss"),
      use_datlist = TRUE,
      datlist = paste0(dir_run, "data.ss")
    )
    
    # Update steepness
    tmp_ctl = SS_readctl(
      file = paste0(dir_run, "control.ss"),
      use_datlist = TRUE,
      datlist = paste0(dir_run, "data.ss")
    )
    # change SRR configuration from survival function to BH SRR
    tmp_ctl$SR_function = 3
    tmp_ctl$SR_parms$PR_type = 0
    tmp_ctl$SR_parms = tmp_ctl$SR_parms[-3, ]
    rownames(tmp_ctl$SR_parms)[2] = "SR_BH_steep"
    tmp_ctl$SR_parms$LO[2] = 0.2
    tmp_ctl$SR_parms$HI[2] = 1
    tmp_ctl$SR_parms$INIT[2] = testing_options_df$steepness[i]
    tmp_ctl$SR_parms$PRIOR[2] = testing_options_df$steepness[i]
    tmp_ctl$SR_parms$PR_SD[2] = 0.05
    #_____________________________________________________________________________________________________________________________
    # Update sigmaR
    tmp_ctl$SR_parms$INIT[3] = testing_options_df$sigmaR[i]
    #_____________________________________________________________________________________________________________________________
    # write-out files
    SS_writectl(
      tmp_ctl,
      paste0(dir_run, "control.ss"),
      version = "3.30",
      overwrite = TRUE
    )
    rm(list = c("tmp_ctl"))
    #_____________________________________________________________________________________________________________________________
    # (FOR TESTING ONLY) modify starter.ss
    if (TESTING)
    {
      tmp_starter = SS_readstarter(file = paste0(dir_run, "starter.ss"),
                                   verbose = TRUE)
      tmp_starter$last_estimation_phase = 0
      SS_writestarter(
        tmp_starter,
        dir = dir_run,
        file = "starter.ss",
        overwrite = TRUE,
        verbose = TRUE,
        warn = TRUE
      )
      rm(list = "tmp_starter")
    }
    
    #_____________________________________________________________________________________________________________________________
    # run SS
    if (TESTING)
    {
      system(paste0(
        "powershell cd ",
        dir_run,
        " ; ./ss_win.exe -nohess 3>&1 2>&1 > out.log"
      ))
      tmp = readLines(paste0(dir_run, "out.log"), skipNul = TRUE)
      if (length(grep("!!  Run has completed  !!", tmp, fixed = TRUE)) ==
          0) {
        rm(list = c("tmp"))
        stop("Error: Bad model settings.")
      }
    } else{
      print(paste0(paste0(
        "powershell cd ", dir_run, " ; ./ss_win.exe"
      )))
      # system(paste0("powershell cd ",dir_run," ; ./ss_win.exe"))
      if (print(Sys.info()['sysname']) == "Darwin") {
        # system(paste0("cd ",dir_run))
        # system("./ss_osx")
      } else if (print(Sys.info()['sysname']) == "Linux") {
        # system(paste0("cd ",dir_run))
        # system("./ss_osx")
      } else if (print(Sys.info()['sysname']) == "Windows") {
        # system(paste0("cd ",dir_run))
        # system("./ss_osx")
      }
    }
  }
  # save summary
  #save(
  # selex_results_df,
  #file = paste0(
  # proj_dir,
  #"Summary/",
  #model_name_stem,
  #"selex_results_df.RData"
  #)
  #)
}
