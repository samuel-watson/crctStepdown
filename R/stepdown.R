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
#' For a set of models fit with lme4, base R, or glmmrBase, the function will conduct the randomisation tests
#' and generate p-values for the null hypotheses of no treatment effect that controls the
#'  family-wise error rate, and generates a 100(1-alpha)% confidence set for the
#'  treatment effect model parameters.
#'
#' @param fitlist A list of models fitted with lme4, base R (lm or glm), or glmmrBase. All models should be fit using the
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
#' @param rand_func String of the name of a function that re-randomises the clusters. The function
#' must take the arguments `nJ` for the number of clusters and `nT` for the number of time periods.
#' The function should produce a data frame that identifies the clusters in the treatment group under the new
#' randomisation scheme. The data frame can either have a single column with name cl_var or
#' two columns of cl_var and t identifying the cluster ID and time period a cluster joins
#' the treatment group. If NULL then clusters are randomised in a 1:1 ratio to treatment and control
#' @param confint Logical indicating whether to run the confidence interval search process
#' @param sigma optional, list of estimated covariance matrices of the observations from the models in fitlist.
#' If provided then the weighted q-score statistic is used.
#' @param ci_start_values Optional list. The list should contain named vectors "upper" and/or "lower" that provide
#' a set of starting values for the upper and/or lower confidence interval searches, respectively. Alternatively,
#' a named scalar `scale` can be provided such that the starting values of the confidence interval search procedure
#' are est +/- `scale`*SE.
#' @param verbose Logical indicating whether to provide detailed output
#' @return A data frame with the point estimates, p-values, and confidence intervals
#' @importFrom methods is
#' @examples
#' \dontshow{
#' setParallelCRT(FALSE) # for the CRAN check
#' }
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
                     ci_start_values = NULL,
                     verbose=TRUE){

  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")
  if(!cl_var%in%colnames(data))stop("cl_var not in colnames(data)")
  if(!all(unlist(lapply(fitlist,function(x)I(is(x,"glmerMod")|is(x,"lmerMod")|is(x,"glm")|is(x,"lm")|is(x,"mcml"))))))stop("All elements of fitlist should be glm, lm, glmer (lme4), lmer (lme4), or mcml (glmmrBase) model objects")
  if(alpha<=0|alpha>=1)stop("alpha should be between 0 and 1")
  if(!type%in%c("rw","b","h","br","hr","none"))stop("type should be one of rw, b, br, h, hr, or none.")

  if(verbose&!is.null(sigma))cat("Covariance matrix provided, weighted quasi-score statistic will be used.\n")
  if(verbose)cat("Extracting treatment effect estimates.\n")
  tr_eff <- rep(NA,length(fitlist))
  tr_sd <- rep(NA,length(fitlist))
  tr_p <- rep(NA,length(fitlist))
  for(i in 1:length(fitlist)){
    if(!is(fitlist[[i]],"mcml")){
      res <- summary(fitlist[[i]])
      tr_eff[i] <- res$coefficients[tr_var,'Estimate']
      tr_sd[i] <- res$coefficients[tr_var,'Std. Error']
      if(ncol(res$coefficients)<4){
        tr_p[i] <- res$coefficients[tr_var,3]
        tr_p[i] <- 2*(1-stats::pnorm(abs(tr_p[i])))
      } else {
        tr_p[i] <- res$coefficients[tr_var,4]
      }
    } else {
      res <- fitlist[[i]]
      tr_eff[i] <- res$coefficients[res$coefficients$par == paste0("b_",tr_var),'est']
      tr_sd[i] <- res$coefficients[res$coefficients$par == paste0("b_",tr_var),'SE']
      tr_p[i] <- res$coefficients[res$coefficients$par == paste0("b_",tr_var),'p']
    }
  }

  tr_st <- rep(NA,length(fitlist))
  tr_icc <- rep(NA,length(fitlist))

  if(verbose)cat("Point estimates: ",round(tr_eff,2),"\n")
  if(verbose)cat("Uncorrected SE: ",round(tr_sd,2),"\n")

  form.cl <- paste0("~factor(",cl_var,")-1")
  Z <- model.matrix(as.formula(form.cl),data=data)
  inv_sigma <- list()

  if(!is.null(sigma)){
    w.opt <- TRUE
    for(i in 1:length(sigma)){
      if(is(sigma[[i]],"list")){
        inv_sigma[[i]] <- list()
        for(j in 1:length(sigma[[i]])){
          inv_sigma[[i]][[j]] <- solve(sigma[[i]][[j]])
        }
      } else {
        inv_sigma[[i]] <- solve(sigma[[i]])
      }
    }
  } else {
    w.opt <- FALSE
    for(i in 1:length(fitlist)){
      inv_sigma[[i]] <- as.matrix(1)
    }
  }

  new_rand <- function(){
    if(is.null(rand_func)){
      df_tr <- sample(unique(data[,cl_var]),ceiling(length(unique(data[,cl_var]))/2))
      trnew <- rep(0,nrow(data))
      trnew[data[,cl_var] %in% df_tr] <- 1
    } else {
      df_tr <- do.call(rand_func,list(nJ=length(unique(data[,cl_var])),
                                      nT=length(unique(data$t))))
      trnew <- rep(0,nrow(data))
      if(!is(df_tr,"data.frame")){
        stop("rand_func must produce a data frame")
      } else {
        if("t" %in% colnames(df_tr)){
          for(i in 1:nrow(df_tr)){
            trnew[data[,cl_var]==df_tr$cl[i]&data$t>=df_tr$t[i]] <- 1
          }
        } else {
          trnew[data[,cl_var] %in% df_tr[,cl_var]] <- 1
        }
      }
    }
    return(trnew)
  }

  if(verbose)cat("Starting permutations...\n")
  trmat <- sapply(1:n_permute,function(i)new_rand())
  dtr <- data[,tr_var]
  dtr[dtr==0] <- -1

  if(type%in%c("rw","br","hr","none")){
    nullfitlist <- list()
    Xnull <- list()
    familylist <- list()
    xb <- list()
    ypred <- list()
    y <- list()
    resids <- list()
    family2list <- list()
    trlist <- list()
    qtest <- list()
    for(i in 1:length(fitlist)){
      out <- est_null_model(fitlist[[i]],
                            data=data,
                            tr_var = tr_var,
                            null_par = 0)

      nullfitlist[[i]] <- out$fit
      Xnull[[i]] <- out$X
      familylist[[i]] <- out$family
      if(out$family[[1]] == "gaussian"){
        xb[[i]] <- out$fit$fitted.values
      } else {
        xb[[i]] <- out$fit$linear.predictors
      }
      y[[i]] <- out$y
      ypred[[i]] <- out$family$linkinv(xb[[i]])
      resids[[i]] <- out$y - ypred[[i]]
      family2list[[i]] <- out$family[[2]]

      tr_st[i] <- qscore_impl(as.vector(resids[[i]]),
                              as.vector(dtr),
                              as.vector(xb[[i]]),
                              as.matrix(inv_sigma[[i]]),
                              as.character(family2list[[i]]),
                              Z,
                              w.opt)

      qtest[[i]] <- as.vector(permutation_test_impl(as.vector(resids[[i]]),
                                                    as.matrix(trmat),
                                                    as.vector(xb[[i]]),
                                                    as.matrix(inv_sigma[[i]]),
                                                    as.character(family2list[[i]]),
                                                    Z,
                                                    w.opt,
                                                    as.numeric(n_permute),
                                                    verbose))
    }

    out <- matrix(Reduce(rbind,qtest), nrow=length(fitlist))


    out <- abs(out)
    tr_st <- abs(tr_st)
    if(verbose)cat("Test statistics: ",round(tr_st,2),"\n")


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
                              label=paste0("p = ",round(tr_p[which(ord_t==(length(ord_t) - (i-1)))],2)))+
            ggplot2::theme_bw()+
            ggplot2::theme(panel.grid = ggplot2::element_blank())
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
                              label=paste0("p = ",round(tr_p[i],2)))+
            ggplot2::theme_bw()+
            ggplot2::theme(panel.grid = ggplot2::element_blank())
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
                              label=paste0("p = ",round(tr_p[ord_t[i]],2)))+
            ggplot2::theme_bw()+
            ggplot2::theme(panel.grid = ggplot2::element_blank())
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
                              label=paste0("p = ",round(tr_p[i],2)))+
            ggplot2::theme_bw()+
            ggplot2::theme(panel.grid = ggplot2::element_blank())
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
    trmat <- sapply(1:nsteps,function(i)new_rand())
    if(verbose)cat("Lower\n")

    if(!is.null(ci_start_values)){
      if(c("lower")%in%tolower(names(ci_start_values))){
        ci_lower_start <- ci_start_values[[which(tolower(names(ci_start_values))=="lower")]]
        if(length(ci_lower_start)!=length(tr_eff))stop("Wrong length start values for lower interval bound")
      } else if(c("scale")%in%tolower(names(ci_start_values))) {
        se_scale <- ci_start_values[[which(tolower(names(ci_start_values))=="scale")]]
        ci_lower_start <- tr_eff-se_scale*tr_sd
      }
    } else {
      ci_lower_start <- tr_eff-3*tr_sd
    }

    ci_lower <- confint_search(start = ci_lower_start,
                            b =  tr_eff,
                            n = nrow(data),
                            nmodel = length(fitlist),
                            Xnull_= Xnull,
                            y = y,
                            tr_ = as.vector(data[,tr_var]),
                            new_tr_mat = as.matrix(trmat),
                            invS = inv_sigma,
                            family = familylist,
                            family2 = family2list,
                            Z = Z,
                            type = type,
                            nsteps=nsteps,
                            weight = w.opt,
                            alpha = alpha,
                            verbose = verbose)


    #out <<- ci_lower
    if(verbose)cat("\nUpper\n")

    if(!is.null(ci_start_values)){
      if(c("upper")%in%tolower(names(ci_start_values))){
        ci_upper_start <- ci_start_values[[which(tolower(names(ci_start_values))=="upper")]]
        if(length(ci_upper_start)!=length(tr_eff))stop("Wrong length start values for upper interval bound")
      } else if(c("scale")%in%tolower(names(ci_start_values))) {
        ci_upper_start <- tr_eff+se_scale*tr_sd
      }
    } else {
      ci_upper_start <- tr_eff+3*tr_sd
    }

    ci_upper <- confint_search(start = ci_upper_start,
                            b =  tr_eff,
                            n = nrow(data),
                            nmodel = length(fitlist),
                            Xnull_= Xnull,
                            y = y,
                            tr_ = as.vector(data[,tr_var]),
                            new_tr_mat = as.matrix(trmat),
                            invS = inv_sigma,
                            family = familylist,
                            family2 = family2list,
                            Z = Z,
                            type = type,
                            nsteps=nsteps,
                            weight = w.opt,
                            alpha = alpha,
                            verbose = verbose)

    if(plots){
      dfl <- as.data.frame(ci_lower$values)
      dfu <- as.data.frame(ci_upper$values)
      dfl <-  reshape2::melt(dfl, id.vars=NULL)
      dfu <- reshape2::melt(dfu, id.vars=NULL)
      dfl$iter <- dfu$iter <- 1:nsteps
      dfl$interval <- "lower"
      dfu$interval <- "upper"
      iter <- value <- NULL
      print(ggplot2::ggplot(data=rbind(dfu,dfl),ggplot2::aes(x=iter,y=value))+
              ggplot2::geom_line()+
              ggplot2::facet_grid(interval~variable))+
              ggplot2::theme_bw()+
              ggplot2::theme(panel.grid = ggplot2::element_blank())
    }

  } else {
    ci_upper <- rep(NA,length(tr_p))
    ci_lower <- rep(NA,length(tr_p))
  }


  if(verbose)cat("\nCompleted!\n")
  results <- data.frame(model=1:length(tr_eff),
                        mean = round(tr_eff,3),
                        lower_ci = round(ci_lower$bound,3),
                        upper_ci = round(ci_upper$bound,3),
                        p_value = round(tr_p,3),
                        icc = round(tr_icc,3))

  return(results)
}

#' Extracts the dependent variable name from glm, lm, or mer model
#'
#' @param fit A fitted model object of class glm, lm, or *merMod
#' @return A string with the name of the dependent variable from the model
#' @importFrom methods is
#' @examples
#' out <- twoarm_sim()
#' data <- out[[1]]
#' fit1 <- lme4::glmer(y1 ~ treat + (1|cl) ,
#'                     data=data,
#'                     family="poisson")
#' outname_fit(fit1)
#' @export
outname_fit <- function(fit){
  if(!(is(fit,"glm")|is(fit,"lm")|is(fit,"glmerMod")|is(fit,"lmerMod")))
    stop("Model class should be glm, lm, or merMod")

  if(is(fit,"glm")){
    outv <- strsplit(as.character(fit$formula), " ")[[2]][1]
  } else if(is(fit,"lm")){
    outv <- strsplit(as.character(fit$call[[2]])," ")[[2]][1]
  } else if(grepl("mer",class(fit))){
    outv <- strsplit(as.character(fit@call[[2]])," ")[[2]]
  } else if(is(fit,"fastglm")|is(fit,"fastLm")){
    outv <- fit$outv
  }
  return(outv)
}
