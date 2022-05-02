


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
run_ss_child <- function(selex_options_df,
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
  #_____________________________________________________________________________________________________________________________
  # loop over options
  #TESTING = FALSE
  for (i in begin[mpi.comm.rank()]:end[mpi.comm.rank()])
  {
    #_____________________________________________________________________________________________________________________________
    # define directory structure
    model_name = paste0(selex_options_df[i,], collapse = "_")
    dir_run = paste0(dir_model, model_name_stem, model_name, "/")
    dir.create(dir_run, recursive = TRUE, showWarnings = FALSE)
    dir_plot = paste0(proj_dir,
                      "Plot/Model_Runs/",
                      model_name_stem,
                      model_name,
                      "/")
    dir.create(dir_plot, recursive = TRUE, showWarnings = FALSE)
    
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
    
    # change selex type
    tmp_ctl$size_selex_types$Pattern[17] = 27
    tmp_ctl$size_selex_types$Male[17] = selex_options_df$male_option[i]
    tmp_ctl$size_selex_types$Special[17] = selex_options_df$nodes[i]
    tmp_ctl$size_selex_types$Pattern[18] = 27
    tmp_ctl$size_selex_types$Male[18] = selex_options_df$male_option[i]
    tmp_ctl$size_selex_types$Special[18] = selex_options_df$nodes[i]
    
    # update selex params
    idx_f17 = grep("F17", rownames(tmp_ctl$size_selex_parms))
    tmp_ctl$size_selex_parms = tmp_ctl$size_selex_parms[-idx_f17,]
    idx_f18 = grep("F18", rownames(tmp_ctl$size_selex_parms))
    tmp_ctl$size_selex_parms = tmp_ctl$size_selex_parms[-idx_f18,]
    idx_f16_tail = tail(grep("F16", rownames(tmp_ctl$size_selex_parms)), n =
                          1)
    # if sex selex offset
    if (selex_options_df$male_option[i] == 2)
    {
      new_f17 = matrix(0,
                       nrow = (3 + 2 * selex_options_df$nodes[i]) + 4,
                       ncol = 14)
    } else {
      new_f17 = matrix(0,
                       nrow = 3 + 2 * selex_options_df$nodes[i],
                       ncol = 14)
    }
    colnames(new_f17) = colnames(tmp_ctl$size_selex_parms)
    new_f17[1, 1:7] = c(-1, 13, selex_options_df$spline_option[i], 0, 0, 0,-1)
    if (selex_options_df$cubic_spline[i])
    {
      new_f17[2, 1:7] = c(0, 2e30, 1e30, 0, 0, 0,-1)
      new_f17[3, 1:7] = c(0, 2e30, 1e30, 0, 0, 0,-1)
    } else{
      new_f17[2, 1:7] = c(-0.001, 1, 0.13, 0, 0, 0,-1)
      new_f17[3, 1:7] = c(-1, 0.001,-0.03, 0, 0, 0,-1)
    }
    # node placement
    tmp_nodes = selex_options_df$nodes[i]
    tmp_nodes_loc = floor(seq(
      from = 50,
      to = 250,
      length.out = tmp_nodes
    ) / 5) * 5
    for (j in 1:tmp_nodes)
    {
      new_f17[3 + j, 1:7] = c(0, 300, tmp_nodes_loc[j], 0, 0, 0,-1)
    }
    # selex values at nodes
    for (j in 1:tmp_nodes)
    {
      if (j == 1)
      {
        new_f17[(3 + tmp_nodes) + j, 1:7] = c(-9, 7, 0, 0, 0, 0,-1)
      } else {
        new_f17[(3 + tmp_nodes) + j, 1:7] = c(-9, 7, 0, 0, 0, 0, 2)
      }
    }
    # if sex selex offset
    if (selex_options_df$male_option[i] == 2)
    {
      new_f17[nrow(new_f17) - 3, 1:7] = c(195, 205, 200, 0, 0, 0,-1)
      new_f17[nrow(new_f17) - 2, 1:7] = c(-1, 1, 0, 0, 0, 0,-1)
      new_f17[nrow(new_f17) - 1, 1:7] = c(-9, 7, 0, 0, 0, 0, 3)
      new_f17[nrow(new_f17), 1:7] = c(-9, 7, 0, 0, 0, 0, 3)
    }
    new_f17 = as.data.frame(new_f17)
    
    tmp_ctl$size_selex_parms = rbind(tmp_ctl$size_selex_parms[1:idx_f16_tail,],
                                     new_f17,
                                     new_f17,
                                     tmp_ctl$size_selex_parms[(idx_f16_tail + 1):nrow(tmp_ctl$size_selex_parms),])
    
    
    SS_writectl(
      tmp_ctl,
      paste0(dir_run, "control.ss"),
      version = "3.30",
      overwrite = TRUE
    )
    rm(list = "tmp_ctl",
       "idx_f17",
       "idx_f18",
       "idx_f16_tail",
       "tmp_nodes",
       "tmp_nodes_loc")
    
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
    setwd(dir_run)
    if (TESTING)
    {
      system("powershell ./ss_win.exe -nohess 3>&1 2>&1 > out.log")
      tmp = readLines("out.log", skipNul = TRUE)
      if (length(grep("!!  Run has completed  !!", tmp, fixed =
                      TRUE)) == 0) {
        rm(list = c("tmp"))
        stop("Error: Bad selex settings.")
      }
    } else{
      system("powershell ./ss_win.exe")
    }
    
    if (!TESTING)
    {
      #_____________________________________________________________________________________________________________________________
      # make plots
      r4ss_model_run = SS_output(dir_run)
      # SS_plots(replist = r4ss_model_run, pdf = FALSE, png = TRUE, html = TRUE, printfolder = "", dir = dir_plot)
      # SS_plots(replist = r4ss_model_run, pdf = TRUE, png = FALSE, html = FALSE, printfolder = "", dir = dir_plot)
      
      #_____________________________________________________________________________________________________________________________
      # extract metrics: LF fit
      selex_results_df$selex_par[i] = sum(new_f17$PHASE > 0) +
        1
      selex_results_df$f17_LL[i] = r4ss_model_run$likelihoods_by_fleet$F17_USA_Lonline_DP[10]
      selex_results_df$f18_LL[i] = r4ss_model_run$likelihoods_by_fleet$F18_USA_Lonline_DP[10]
      
      #_____________________________________________________________________________________________________________________________
      # compress SS output
      tmp_files = list.files()
      system(paste0(
        "powershell tar -czf ss_files.tar.gz ",
        paste0(tmp_files, collapse = " ")
      ))
      
      #_____________________________________________________________________________________________________________________________
      # clean-up
      file.remove(tmp_files)
      rm(
        list = c(
          "tmp_out",
          "pointer",
          "new_f17",
          "r4ss_model_run",
          "tmp_files",
          "model_name"
        )
      )
    }
  }
  
  return (c(begin[mpi.comm.rank()], end[mpi.comm.rank()]))
}

#Run model in parallel
run_ss_ht4sa_MPI <- function(selex_options_df,
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
  id <- mpi.comm.rank(comm = 0)
  ns <-  mpi.universe.size() - 1
  
  
  nsims <- nrow(selex_options_df)
  
  begin <- rep(0, ns)
  end <- rep(0, ns)
  
  #create scenario segments
  print(paste("id = ", id))
  if (id == 0) {
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
    print("launching children....")
    mpi.spawn.Rslaves(nslaves = ns)
    
    #pass the function to all slaves
    mpi.bcast.Robj2slave(obj = run_ss_child)
    
    
    print("broadcasting run....")
    #execute all children
    x <- mpi.remote.exec(
      run_ss_child,
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
    
    # Tell all slaves to close down, and exit the program
    mpi.close.Rslaves(dellog = FALSE)
  }
  # save summary
  save(
    selex_results_df,
    file = paste0(
      proj_dir,
      "Summary/",
      model_name_stem,
      "selex_results_df.RData"
    )
  )
}


#run models sequentially
run_ht4sa_ss_local <- function(selex_options_df,
                               proj_dir,
                               dir_utility,
                               dir_model,
                               dir_model_source,
                               model_name_stem,
                               dir_input_files,
                               dir_ss,
                               TESTING = FALSE) {
  # loop over options
  #TESTING = FALSE
  for (i in 1:nrow(selex_options_df))
  {
    #_____________________________________________________________________________________________________________________________
    # define directory structure
    model_name = paste0(selex_options_df[i,], collapse =
                          "_")
    dir_run = paste0(dir_model, model_name_stem, model_name, "/")
    dir.create(dir_run, recursive =
                 TRUE, showWarnings = FALSE)
    dir_plot = paste0(proj_dir,
                      "Plot/Model_Runs/",
                      model_name_stem,
                      model_name,
                      "/")
    dir.create(dir_plot, recursive =
                 TRUE, showWarnings = FALSE)
    
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
    
    # change selex type
    tmp_ctl$size_selex_types$Pattern[17] = 27
    tmp_ctl$size_selex_types$Male[17] = selex_options_df$male_option[i]
    tmp_ctl$size_selex_types$Special[17] = selex_options_df$nodes[i]
    tmp_ctl$size_selex_types$Pattern[18] = 27
    tmp_ctl$size_selex_types$Male[18] = selex_options_df$male_option[i]
    tmp_ctl$size_selex_types$Special[18] = selex_options_df$nodes[i]
    
    # update selex params
    idx_f17 = grep("F17", rownames(tmp_ctl$size_selex_parms))
    tmp_ctl$size_selex_parms = tmp_ctl$size_selex_parms[-idx_f17,]
    idx_f18 = grep("F18", rownames(tmp_ctl$size_selex_parms))
    tmp_ctl$size_selex_parms = tmp_ctl$size_selex_parms[-idx_f18,]
    idx_f16_tail = tail(grep("F16", rownames(tmp_ctl$size_selex_parms)), n =
                          1)
    # if sex selex offset
    if (selex_options_df$male_option[i] ==
        2)
    {
      new_f17 = matrix(0,
                       nrow = (3 + 2 * selex_options_df$nodes[i]) + 4,
                       ncol = 14)
    } else {
      new_f17 = matrix(0,
                       nrow = 3 + 2 * selex_options_df$nodes[i],
                       ncol = 14)
    }
    colnames(new_f17) = colnames(tmp_ctl$size_selex_parms)
    new_f17[1, 1:7] = c(-1, 13, selex_options_df$spline_option[i], 0, 0, 0,-1)
    if (selex_options_df$cubic_spline[i])
    {
      new_f17[2, 1:7] = c(0, 2e30, 1e30, 0, 0, 0,-1)
      new_f17[3, 1:7] = c(0, 2e30, 1e30, 0, 0, 0,-1)
    } else{
      new_f17[2, 1:7] = c(-0.001, 1, 0.13, 0, 0, 0,-1)
      new_f17[3, 1:7] = c(-1, 0.001,-0.03, 0, 0, 0,-1)
    }
    # node placement
    tmp_nodes = selex_options_df$nodes[i]
    tmp_nodes_loc = floor(seq(
      from = 50,
      to = 250,
      length.out = tmp_nodes
    ) / 5) * 5
    for (j in 1:tmp_nodes)
    {
      new_f17[3 + j, 1:7] = c(0, 300, tmp_nodes_loc[j], 0, 0, 0,-1)
    }
    # selex values at nodes
    for (j in 1:tmp_nodes)
    {
      if (j == 1)
      {
        new_f17[(3 + tmp_nodes) + j, 1:7] = c(-9, 7, 0, 0, 0, 0,-1)
      } else {
        new_f17[(3 + tmp_nodes) + j, 1:7] = c(-9, 7, 0, 0, 0, 0, 2)
      }
    }
    # if sex selex offset
    if (selex_options_df$male_option[i] ==
        2)
    {
      new_f17[nrow(new_f17) - 3, 1:7] = c(195, 205, 200, 0, 0, 0,-1)
      new_f17[nrow(new_f17) -
                2, 1:7] = c(-1, 1, 0, 0, 0, 0,-1)
      new_f17[nrow(new_f17) -
                1, 1:7] = c(-9, 7, 0, 0, 0, 0, 3)
      new_f17[nrow(new_f17), 1:7] = c(-9, 7, 0, 0, 0, 0, 3)
    }
    new_f17 = as.data.frame(new_f17)
    
    tmp_ctl$size_selex_parms = rbind(tmp_ctl$size_selex_parms[1:idx_f16_tail,],
                                     new_f17,
                                     new_f17,
                                     tmp_ctl$size_selex_parms[(idx_f16_tail + 1):nrow(tmp_ctl$size_selex_parms),])
    
    
    SS_writectl(
      tmp_ctl,
      paste0(dir_run, "control.ss"),
      version = "3.30",
      overwrite = TRUE
    )
    rm(list = "tmp_ctl",
       "idx_f17",
       "idx_f18",
       "idx_f16_tail",
       "tmp_nodes",
       "tmp_nodes_loc")
    
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
    setwd(dir_run)
    if (TESTING)
    {
      system("powershell ./ss_win.exe -nohess 3>&1 2>&1 > out.log")
      tmp = readLines("out.log", skipNul = TRUE)
      if (length(grep("!!  Run has completed  !!", tmp, fixed =
                      TRUE)) == 0) {
        rm(list = c("tmp"))
        stop("Error: Bad selex settings.")
      }
    } else{
      system("powershell ./ss_win.exe")
    }
    
    if (!TESTING)
    {
      #_____________________________________________________________________________________________________________________________
      # make plots
      r4ss_model_run = SS_output(dir_run)
      # SS_plots(replist = r4ss_model_run, pdf = FALSE, png = TRUE, html = TRUE, printfolder = "", dir = dir_plot)
      # SS_plots(replist = r4ss_model_run, pdf = TRUE, png = FALSE, html = FALSE, printfolder = "", dir = dir_plot)
      
      #_____________________________________________________________________________________________________________________________
      # extract metrics: LF fit
      selex_results_df$selex_par[i] = sum(new_f17$PHASE >
                                            0) + 1
      selex_results_df$f17_LL[i] = r4ss_model_run$likelihoods_by_fleet$F17_USA_Lonline_DP[10]
      selex_results_df$f18_LL[i] = r4ss_model_run$likelihoods_by_fleet$F18_USA_Lonline_DP[10]
      
      #_____________________________________________________________________________________________________________________________
      # compress SS output
      tmp_files = list.files()
      system(paste0(
        "powershell tar -czf ss_files.tar.gz ",
        paste0(tmp_files, collapse = " ")
      ))
      
      #_____________________________________________________________________________________________________________________________
      # clean-up
      file.remove(tmp_files)
      rm(
        list = c(
          "tmp_out",
          "pointer",
          "new_f17",
          "r4ss_model_run",
          "tmp_files",
          "model_name"
        )
      )
    }
  }
  # save summary
  save(
    selex_results_df,
    file = paste0(
      proj_dir,
      "Summary/",
      model_name_stem,
      "selex_results_df.RData"
    )
  )
}
