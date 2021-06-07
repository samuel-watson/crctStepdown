#' Randomisation test confidence interval limit search
#'
#' A multi-variate Robbins-Monroe search process to estimate the upper of lower
#' limits for a confidence set for parameters from a list of fitted model objects.
#'
#' A version of the search process proposed by Garthwaite (1996) adapted for multiple
#' limits. Given a set of estimates of the upper or lower, the process calculates the test
#' statistics for the two-sided null hypotheses that the treatment effects equal these values,
#' and then conducts a single iteration randomisation test of the same null hypotheses. The
#' estimates are then updated based on whether the actual test statistic is higher or lower
#' than the randomisation test statistic it would be compared to under the resampling stepdown
#' approach of Romano & Wolf (2005). At the limits of the confidence set all values should
#' be rejected in a two-sided hypothesis test with a family-wise error rate of alpha, which
#' provides a probabilistic basis for the search process. See Watson (2021) for more details.
#'
#' @param fitlist A list of p fitted mer class model objects
#' @param data A data frame used to fit the models in fitlist
#' @param actual_tr A vector of length p with the original point estimates of the
#' treatment effects
#' @param start A vector of length p with a set of starting values for the search. The uncorrected
#' confidence set limits (tr_eff +/- 2*SE) usually provide a good starting point.
#' @param nsteps Number of steps for the search process.
#' @param alpha Numeric value. The process searches for the 100(1-2*alpha) confidence intervals
#' @param plots Logical value indicating whether to plot the search process. Default to TRUE
#' @param cl_var String indicating the name of the column in data with the IDs of the clusters
#' @param rand_func The name of a function that re-randomises the clusters. The function should
#' produce a data frame that identifies the clusters in the treatment group under the new
#' randomisation scheme. The data frame can either have a single column with name cl_var or
#' two columns of cl_var and t identifying the cluster ID and time period a cluster joins
#' the treatment group.  If NULL then clusters are randomised in a 1:1 ratio to treatment and control
#' @param verbose Logical indicating whether to provide verbose output showing progress and estimates
#' @return A vector of length p with the estimates of the limits
#' @importFrom methods is
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#' @export
conf_int_search <- function(fitlist,
                            data,
                            actual_tr,
                            start,
                            nsteps=1000,
                            alpha=0.025,
                            plots=TRUE,
                            cl_var = "cl",
                            rand_func = NULL,
                            verbose=TRUE){
  if(!is(fitlist,"list"))stop("fitlist should be a list.")
  if(!all(unlist(lapply(fitlist,function(x)I(is(x,"glmerMod")|is(x,"lmerMod"))))))stop("All elements of fitlist should be lme4 model objects")
  p <- length(fitlist)
  if(length(actual_tr)!=p)stop("length(actual_tr)!=length(fitlist)")
  if(length(start)!=p)stop("length(actual_tr)!=length(fitlist)")

  bound <- start
  k <- 2/(1.96*((2*pi)^-0.5)*exp((-1.96^2)/2))
  dfv <- matrix(NA,nrow=nsteps,ncol=length(actual_tr))
  actual_t <- rep(NA,length(fitlist))

  for(i in 1:nsteps){
    if(verbose)cat("\rStep: ",i," of ",nsteps,": ",bound)
    nullfitlist <- list()
    for(j in 1:length(fitlist)){
      nullfitlist[[j]] <- est_null_model(fitlist[[j]],
                                         data,
                                         tr_var="treat",
                                         null_par = bound[j])
    }
    dfv[i,] <- bound
    for(j in 1:length(fitlist)){
      actual_t[j] <- qscore_stat(nullfitlist[[j]],
                                 data,
                                 null_par=bound[j],
                                 tr_assign = "treat")
    }


    #order of statistics
    val <- lme_permute2(nullfitlist,
                        data,
                        null_par = bound,
                        cl_var=cl_var,
                        rand_func=rand_func)
    actual_t <- abs(actual_t)
    val <- abs(val)
    pos_t <- order(actual_t)
    t.st <- rep(NA,length(actual_t))
    p.st <- rep(NA,length(actual_t))
    step <- rep(NA,length(actual_t))
    for(j in 1:length(actual_t)){
      t.st[pos_t[(length(pos_t) - (j-1))]] <- max(actual_t[pos_t[1:(length(pos_t) - (j-1))]])
      p.st[pos_t[(length(pos_t) - (j-1))]] <- max(val[pos_t[1:(length(pos_t) - (j-1))]])
      step[pos_t[(length(pos_t) - (j-1))]] <- k*(actual_tr[pos_t[(length(pos_t) - (j-1))]] -
                                                   bound[pos_t[(length(pos_t) - (j-1))]])
    }
    rjct <- I(t.st > p.st)

    for(j in 1:length(actual_t)){
      if(rjct[j]){
        bound[j] <- bound[j] + step[j]*(alpha)/i
      } else {
        bound[j] <- bound[j] - step[j]*(1-alpha)/i
      }
    }

  }

  if(plots){

    dfv <- data.frame(step=rep(1:nsteps,length(actual_tr)),
                      par=rep(1:length(actual_tr),each=nsteps),
                      value=c(dfv))
    p1 <- ggplot2::ggplot(data=dfv,aes(x=.data$step,y=.data$value))+
      ggplot2::geom_line()+
      ggplot2::facet_wrap(~.data$par, scales = "free_y")
    print(p1)
  }

  return(bound)
}
