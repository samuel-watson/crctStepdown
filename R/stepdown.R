#' Extracts the test statistics
#'
#' Extracts the test statistics from the output of the \code{permute} function.
#' Returns the largest value from a specified subset of rows, each row is the test
#' statistic from a different null hypothesis.
#' @param out Array output by the \code{permute} function
#' @param positions Vector indicating which rows of out to use
#' @return Vector of numeric values of length \code{ncol(out)}
perm_dist <- function(out,
                      positions){
  out <- out[positions,]
  out <- abs(out)
  if(length(positions)>1){
    vals <- apply(out,2,max)
  } else {
    vals <- out
  }
  return(vals)
}

#' Conduct the randomisation-based stepdown procedure
#'
#' For a set of models fit with lme4, the function will conduct the randomisation tests
#' and generate p-values for the null hypotheses of no treatment effect that controls the
#'  family-wise error rate, and generates a 100(1-alpha)% confidence set for the
#'  treatment effect model parameters.
#'
#' @param fitlist A list of models fitted with lme4. All models should be fit using the
#' same data frame.
#' @param tr_var String indicating the name of the column in data that is a binary indicator
#' for whether the observation was under the treatment (1=treatment, 0=control)
#' @param cl_var String specifying the name of the column identifying the clusters/cluster-time
#' @param data A data frame containing the data used to fit the models in fitlist
#' @param alpha Numeric. 100(1-alpha)% confidence intervals are calculated. Default it 0.05
#' @param plots Logical indicating whether to plot permutational distributions and confidence
#' interval search during running of function. Default is TRUE
#' @param n_permute Number of permutations of the randomisation test to run
#' @param nsteps Number of steps of the confidence interval search process
#' @param type Method of correction: options are "rw" = Romano-Wolf randomisation test based stepdown, "h"
#' = Holm standard stepdown, "h" = Holm stepdown using randomisation test, "b" = standard Bonferroni, "br" =
#' Bonerroni using randomisation test, or "none" = randomisation test with no correction.
#' @param rand_func String of the name of a function that re-randomises the clusters. The function should
#' produce a data frame that identifies the clusters in the treatment group under the new
#' randomisation scheme. The data frame can either have a single column with name cl_var or
#' two columns of cl_var and t identifying the cluster ID and time period a cluster joins
#' the treatment group. If NULL then clusters are randomised in a 1:1 ratio to treatment and control
#' @param confint Logical indicating whether to run the confidence interval search process
#' @param sigma optional, estimated covariance matrix of the observations. If provided then the weighted q-score statistic is used.
#' @param verbose Logical indicating whether to provide detailed output
#' @return A data frame with the point estimates, p-values, and confidence intervals
#' @importFrom methods is
#' @examples
#' out <- twoarm_sim()
#' data <- out[[1]]
#'   fit1 <- lme4::glmer(y1 ~ treat + (1|cl) ,
#' data=data,
#' family="poisson")
#'
#' fit2 <- lme4::glmer(y2 ~ treat + (1|cl),
#'                     data=data,
#'                     family="poisson")
#'  stepdown(fitlist=list(fit1,fit2),
#'           data=data,
#'           n_permute = 100,
#'           nsteps=100,
#'           plots=FALSE,
#'           verbose=TRUE)
#' @export
stepdown <- function(fitlist,
                     tr_var = "treat",
                     cl_var = "cl",
                     data,
                     alpha=0.05,
                     plots=TRUE,
                     n_permute=1000,
                     nsteps=1000,
                     type="rw",
                     rand_func=NULL,
                     confint=TRUE,
                     sigma = NULL,
                     verbose=TRUE){

  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")
  if(!cl_var%in%colnames(data))stop("cl_var not in colnames(data)")
  if(!all(unlist(lapply(fitlist,function(x)I(is(x,"glmerMod")|is(x,"lmerMod")|is(x,"glm")|is(x,"lm"))))))stop("All elements of fitlist should be glm, lm, glmer, or lmer model objects")
  if(alpha<=0|alpha>=1)stop("alpha should be between 0 and 1")
  if(!type%in%c("rw","b","h","br","hr","none"))stop("type should be one of rw, b, br, h, hr, or none.")

  if(verbose)cat("Extracting treatment effect estimates.\n")
  tr_eff <- rep(NA,length(fitlist))
  tr_sd <- rep(NA,length(fitlist))
  tr_p <- rep(NA,length(fitlist))
  for(i in 1:length(fitlist)){
    res <- summary(fitlist[[i]])
    tr_eff[i] <- res$coefficients[tr_var,'Estimate']
    tr_sd[i] <- res$coefficients[tr_var,'Std. Error']
    if(ncol(res$coefficients)<4){
      tr_p[i] <- res$coefficients[tr_var,3]
      tr_p[i] <- 2*(1-pnorm(abs(tr_p[i])))
    } else {
      tr_p[i] <- res$coefficients[tr_var,4]
    }

  }

  tr_st <- rep(NA,length(fitlist))
  tr_icc <- rep(NA,length(fitlist))

  if(verbose)cat("Point estimates: ",round(tr_eff,2),"\n")
  if(verbose)cat("Uncorrected SE: ",round(tr_sd,2),"\n")

  if(!is.null(sigma)){
    inv_sigma <- solve(sigma)
  } else {
    inv_sigma <- NULL
  }

  if(type%in%c("rw","br","hr","none")){
    nullfitlist <- list()
    for(i in 1:length(fitlist)){
      nullfitlist[[i]] <- est_null_model(fitlist[[i]],
                                         data=data,
                                         tr_var = tr_var,
                                         null_par = 0)
    }

    for(i in 1:length(fitlist)){
      tr_st[i] <- qscore_stat(nullfitlist[[i]],
                              data=data,
                              tr_var = tr_var,
                              cl_var = cl_var,
                              tr_assign = tr_var,
                              inv_sigma = inv_sigma)
    }


    #cat("ICC: ",round(tr_icc,3),"\n")
    if(verbose)cat("Starting permutations...\n")
    out <- permute(nullfitlist,
                   data,
                   n_permute = n_permute,
                   cl_var = cl_var,
                   tr_var = tr_var,
                   rand_func=rand_func,
                   inv_sigma = inv_sigma)
    anyna <- apply(out,2,function(x)any(is.na(x)))
    out <- out[,!anyna]
    out <- abs(out)
    tr_st <- abs(tr_st)
    if(verbose)cat("Test statistics: ",round(tr_st,2),"\n")
    # print(apply(out,1,sd))

    #first determine corrected p-values
    #test statistics
    if(verbose)cat("Calculating p-values\n")
    tr_p <- rep(NA,length(fitlist))
    #order of statistics
    ord_t <- order(tr_st)
    if(plots)plist <- list()

    if(type=="rw"){
      for(i in 1:length(ord_t)){
        test_stat <- max(tr_st[ord_t[1:(length(ord_t) - (i-1))]])
        vals <- perm_dist(out,ord_t[1:(length(ord_t) - (i-1))])
        tr_p[ord_t[(length(ord_t) - (i-1))]] <- (1+length(vals[vals>test_stat]))/(length(vals)+1)

        if(plots){
          plist[[i]] <- ggplot2::qplot(vals,bins=30) +
            ggplot2::geom_vline(xintercept = test_stat,color="red")+
            ggplot2::labs(x=paste0(i,"th largest statistic"))+
            ggplot2::annotate("text",x = Inf, y= Inf,hjust=1,vjust=1,color="red",
                              label=paste0("p = ",round(tr_p[which(ord_t==(length(ord_t) - (i-1)))],2)))
        }

      }
    }

    if(type=="br"){
      for(i in 1:length(ord_t)){
        vals <- out[i,]
        tr_p[i] <- (1+length(vals[vals>tr_st[i]]))/(length(vals)+1)
        tr_p[i] <- min(tr_p[i]*length(fitlist),1)

        if(plots){
          plist[[i]] <- ggplot2::qplot(vals,bins=30) +
            ggplot2::geom_vline(xintercept = tr_st[i],color="red")+
            ggplot2::annotate("text",x = Inf, y= Inf,hjust=1,vjust=1,color="red",
                              label=paste0("p = ",round(tr_p[i],2)))
        }
      }
    }

    if(type=="hr"){
      for(i in 1:length(ord_t)){
        vals <- out[ord_t[(length(ord_t) - (i-1))],]
        tr_p[ord_t[(length(ord_t) - (i-1))]] <- (1+length(vals[vals>tr_st[ord_t[(length(ord_t) - (i-1))]]]))/(length(vals)+1)
        tr_p[ord_t[(length(ord_t) - (i-1))]] <- min(tr_p[ord_t[(length(ord_t) - (i-1))]]*(length(fitlist)+1-i),1)

        if(plots){
          plist[[i]] <- ggplot2::qplot(vals,bins=30) +
            ggplot2::geom_vline(xintercept = tr_st[ord_t[i]],color="red")+
            ggplot2::labs(x=paste0(i,"th largest statistic"))+
            ggplot2::annotate("text",x = Inf, y= Inf,hjust=1,vjust=1,color="red",
                              label=paste0("p = ",round(tr_p[ord_t[i]],2)))
        }
      }
    }

    if(type=="none"){
      for(i in 1:length(ord_t)){
        vals <- out[i,]
        tr_p[i] <- (1+length(vals[vals>tr_st[i]]))/(length(vals)+1)

        if(plots){
          plist[[i]] <- ggplot2::qplot(vals,bins=30) +
            ggplot2::geom_vline(xintercept = tr_st[i],color="red")+
            ggplot2::annotate("text",x = Inf, y= Inf,hjust=1,vjust=1,color="red",
                              label=paste0("p = ",round(tr_p[i],2)))
        }
      }
    }

    if(plots) print(ggpubr::ggarrange(plotlist = plist))
  }

  if(type=="b"){
    tr_p <- ifelse(tr_p * length(fitlist) > 1, 1, tr_p * length(fitlist) )
  }

  if(type=="h"){
    ord_p <- order(tr_p)
    for(i in 1:length(fitlist)){
      tr_p[ord_p[(length(ord_p) - (i-1))]] <- tr_p[ord_p[(length(ord_p) - (i-1))]]*(length(fitlist)-i+1)
    }
    tr_p <- ifelse(tr_p>1,1,tr_p)
  }

  if(verbose)cat("P-values: ",round(tr_p,2),"\n")

  # next use search algorithm to determine 95% confidence interval
  if(confint){
    if(verbose)cat("Searching for confidence intervals...\n")
    if(verbose)cat("Upper intervals:\n")
    ci_upper <- conf_int_search(fitlist,
                                data,
                                actual_tr=tr_eff,
                                start=tr_eff+2.5*tr_sd,
                                nsteps=nsteps,
                                alpha=alpha,
                                plots = plots,
                                cl_var = cl_var,
                                tr_var = tr_var,
                                rand_func=rand_func,
                                verbose = verbose,
                                type=type,
                                sigma = sigma)

    if(verbose)cat("\nLower intervals:\n")
    ci_lower <- conf_int_search(fitlist,
                                data,
                                actual_tr=tr_eff,
                                start=tr_eff-2.5*tr_sd,
                                nsteps=nsteps,
                                alpha=alpha,
                                plots=plots,
                                cl_var = cl_var,
                                tr_var = tr_var,
                                rand_func=rand_func,
                                verbose = verbose,
                                type=type,
                                sigma = sigma)

  } else {
    ci_upper <- rep(NA,length(tr_p))
    ci_lower <- rep(NA,length(tr_p))
  }


  if(verbose)cat("\nCompleted!\n")
  results <- data.frame(model=1:length(tr_eff),
                        mean = round(tr_eff,3),
                        lower_ci = round(ci_lower,3),
                        upper_ci = round(ci_upper,3),
                        p_value = round(tr_p,3),
                        icc = round(tr_icc,3))

  return(results)
}
