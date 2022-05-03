

# Nicholas Ducharme-Barth
# 05/02/2022
# htc4sa
# Simple example based on 2017 ISC north pacific blue shark stock assessment
# Set-up small job array
# 1) change steepness
# 2) change sigmaR

# Caveat: model results could be nonsense since this model was tuned to run with the LFSR relationship

#_____________________________________________________________________________________________________________________________
# load packages
	library(r4ss)

#_____________________________________________________________________________________________________________________________
# set directory paths
	proj_dir = "D:/HOME/SAP/htc4sa_pilot/test_case/"
	dir_input_files = paste0(proj_dir,"input_files/")
	dir_model = paste0(proj_dir,"output_files/")
	dir_ss = paste0(proj_dir,"executable/SS3/3.30.19.01.safe/")

#_____________________________________________________________________________________________________________________________
# set-up testing options
	testing_options_df = expand.grid(steepness=c(0.45,0.6,0.75),
								     sigmaR=c(0.3,0.4,0.5),
								     stringsAsFactors=FALSE)

#_____________________________________________________________________________________________________________________________
# extract indices and run DFA to create composite index
	TESTING = FALSE
	for(i in 1:nrow(testing_options_df))
	{	
		#_____________________________________________________________________________________________________________________________
		# define directory structure
			model_name = paste0(testing_options_df[i,],collapse="_")
			dir_run = paste0(dir_model,model_name,"/") 
			dir.create(dir_run,recursive=TRUE,showWarnings=FALSE)

		#_____________________________________________________________________________________________________________________________
		# transfer files
			# intial files
				FileList=c("starter.ss","forecast.ss","data.ss","control.ss")
				file.copy(paste0(dir_input_files,FileList),dir_run,overwrite=TRUE)

			# ss files
				FileList="ss_win.exe"
				file.copy(paste0(dir_ss,FileList),dir_run,overwrite=TRUE)
		# Update starter.ss to read correct control and data files
			tmp_starter = SS_readstarter(file = paste0(dir_run,"starter.ss"), verbose = TRUE)
			tmp_starter$datfile = "data.ss"
			tmp_starter$ctlfile = "control.ss"
			SS_writestarter(tmp_starter, dir = dir_run, file = "starter.ss",overwrite = TRUE, verbose = TRUE, warn = TRUE)
			rm(list="tmp_starter")

		#_____________________________________________________________________________________________________________________________
		# Update steepness
			tmp_ctl = SS_readctl(file=paste0(dir_run,"control.ss"),use_datlist = TRUE,datlist = paste0(dir_run,"data.ss"))
			# change SRR configuration from survival function to BH SRR
				tmp_ctl$SR_function = 3
				tmp_ctl$SR_parms$PR_type = 0
				tmp_ctl$SR_parms = tmp_ctl$SR_parms[-3,]
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
			SS_writectl(tmp_ctl,paste0(dir_run,"control.ss"), version = "3.30", overwrite = TRUE)
			rm(list=c("tmp_ctl"))
		#_____________________________________________________________________________________________________________________________
		# (FOR TESTING ONLY) modify starter.ss
			if(TESTING)
			{
				tmp_starter = SS_readstarter(file = paste0(dir_run,"starter.ss"), verbose = TRUE)
				tmp_starter$last_estimation_phase = 0
				SS_writestarter(tmp_starter, dir = dir_run, file = "starter.ss",overwrite = TRUE, verbose = TRUE, warn = TRUE)
				rm(list="tmp_starter")
			} 

		#_____________________________________________________________________________________________________________________________
		# run SS
			if(TESTING)
			{
				system(paste0("powershell cd ",dir_run," ; ./ss_win.exe -nohess 3>&1 2>&1 > out.log"))
				tmp = readLines(paste0(dir_run,"out.log"), skipNul = TRUE)
				if(length(grep("!!  Run has completed  !!",tmp,fixed=TRUE))==0){rm(list=c("tmp"));stop("Error: Bad model settings.")}
			} else{
				system(paste0("powershell cd ",dir_run," ; ./ss_win.exe"))
			}		
	}

