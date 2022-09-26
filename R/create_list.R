CreateList<-function(nyr, nages, fleet_num, L, survey_num, survey, L.age, survey.age, cv.L, cv.survey, n.L, n.survey){
  
  ## Purpose: Add lognormal observation error to landings and survey, and sampling error to age comps
  ##
  ## INPUT DATA:
  ##   L, survey = error free time series of landings and survey
  ##   L.age, survey.age = error free composition data from which to draw samples
  ## INPUT PARAMETERS:
  ##   cv.L, cv.survey = CVs
  ##   n.L, n.survey = annual sample sizes of age comps
  ##
  ## OUTPUT:
  ##   time series with lognormal error
  ##   comp data with sampling error
  
  #####################################################################################
  #### True landings and survey index ####
  nobs.L <- vector(mode="list", length=fleet_num)
  sd.L <- vector(mode="list", length=fleet_num)
  ln.L <- vector(mode="list", length=fleet_num)
  L.obs <- vector(mode="list", length=fleet_num)
  names(L.obs) <- paste("fleet", 1:fleet_num, sep="")
  
  invisible(sapply(1:fleet_num, function(x){
    nobs.L[[x]] <- length(L[[x]])
    sd.L[[x]] <- sqrt(log(1+cv.L[[x]]^2)) #SD in log space, given CV in arithmetic space
    ln.L[[x]] <- rnorm(nobs.L[[x]], mean=0, sd=sd.L[[x]]) #log error
    L.obs[[x]] <<- L[[x]]*exp(ln.L[[x]]) #multiplicative lognormal error
  }))
  
  nobs.survey <- vector(mode="list", length=survey_num)
  sd.survey <- vector(mode="list", length=survey_num)
  ln.survey <- vector(mode="list", length=survey_num)
  survey.obs <- vector(mode="list", length=survey_num)
  names(survey.obs) <- paste("survey", 1:survey_num, sep="")
  
  invisible(sapply(1:survey_num, function(x){
    nobs.survey[[x]] <- length(survey[[x]])
    sd.survey[[x]] <- sqrt(log(1+cv.survey[[x]]^2))
    ln.survey[[x]] <<- rnorm(nobs.survey[[x]], mean=0, sd=sd.survey[[x]])
    survey.obs[[x]] <<- survey[[x]]*exp(ln.survey[[x]])
  }))
  
  #### Sampling of age comps (rows=years, columns=ages) ####
  
  L.age.obs <- vector(mode="list", length=fleet_num)
  names(L.age.obs) <- paste("fleet", 1:fleet_num, sep="")
  invisible(sapply(1:length(L.age.obs), function(x) L.age.obs[[x]] <<- matrix(0, nrow=nyr, ncol=nages)))
  
  invisible(sapply(1:fleet_num, function(x){
    for (i in 1:length(L.obs[[x]])){
      if(sum(L.age[[x]][i,])==0) {
        probs <- rep(0, length(L.age[[x]][i,]))
        L.age.obs[[x]][i,] <<- rep(0, length(L.age[[x]][i,]))
      } else {
        probs <- L.age[[x]][i,]/sum(L.age[[x]][i,])
        L.age.obs[[x]][i,] <<- rmultinom(n=1, size=n.L[[x]], prob=probs)/n.L[[x]]
      }
    }
  }))
  
  survey.age.obs <- vector(mode="list", length=survey_num)
  names(survey.age.obs) <- paste("survey", 1:survey_num, sep="")
  invisible(sapply(1:length(survey.age.obs), function(x) survey.age.obs[[x]] <<- matrix(0, nrow=nyr, ncol=nages)))
  
  invisible(sapply(1:survey_num, function(x){
    for (i in 1:length(survey.obs[[x]])){
      if (sum(survey.age[[x]][i,])==0) {
        probs <- rep(0, length(survey.age[[x]][i,]))
        survey.age[[x]][i,] <<- rep(0, length(survey.age[[x]][i,]))
      } else {
        probs <- survey.age[[x]][i,]/sum(survey.age[[x]][i,])
        survey.age.obs[[x]][i,] <<- rmultinom(n=1, size=n.survey[[x]], prob=probs)/n.survey[[x]]
      }
    }
  }))
  
  return(list(L.obs=L.obs, survey.obs=survey.obs, L.age.obs=L.age.obs, survey.age.obs=survey.age.obs, n.L=n.L, n.survey=n.survey, survey_q=om_output$survey_q))
  
} # end ObsModel