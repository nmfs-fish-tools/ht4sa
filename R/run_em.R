#' @export
run_em <- function(em_names=NULL, input_list=NULL, em_input_filenames=NULL){

  if (!file.exists(file.path(maindir, "em_input"))) stop ("Missing estimation model input file!")
  if (is.null(em_names)) stop ("Missing EM information!")

  maindir <- input_list$maindir
  om_sim_num <- input_list$om_sim_num
  case_name <- input_list$case_name
  casedir <- file.path(maindir, case_name)
  em_bias_cor <- input_list$em_bias_cor
  initial_equilibrium_F <- input_list$initial_equilibrium_F

  invisible(sapply(em_names, function(x) {
    if (!file.exists(file.path(casedir, "output", x))) dir.create(file.path(casedir, "output", x))
  }))

  if("AMAK" %in% em_names) run_amak(maindir=maindir, om_sim_num=om_sim_num, casedir=casedir, input_filename=em_input_filenames$AMAK)
  if("ASAP" %in% em_names) run_asap(maindir=maindir, om_sim_num=om_sim_num, casedir=casedir, input_filename=em_input_filenames$ASAP)
  if("BAM" %in% em_names) run_bam(maindir=maindir, om_sim_num=om_sim_num, casedir=casedir, em_bias_cor=em_bias_cor, input_filename=em_input_filenames$BAM)
  if("SS" %in% em_names) run_ss(maindir=maindir, om_sim_num=om_sim_num, casedir=casedir, em_bias_cor=em_bias_cor, input_filename=em_input_filenames$SS, initial_equilibrium_F=initial_equilibrium_F)
  if("MAS" %in% em_names) run_mas(maindir=maindir, om_sim_num=om_sim_num, casedir=casedir, em_bias_cor=em_bias_cor)
}
