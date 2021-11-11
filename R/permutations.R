#' Generate a new permutation
#'
#' Returns the test statistic from a specified null hypothesis and model under a single new
#' permutation
#'
#' @param fitlist A list of glm model objects fitted under the null hypotheses
#' @param data A data frame containing the data used to fit the models in fitlist
#' @param null_par A vector of the same length as fitlist specifying the value(s) of the
#' treatment effect parameter(s) under the null hypotheses
#' @param cl_var String specifying the name of the column identifying the clusters/cluster-time
#' @param rand_func String of the name of a function that re-randomises the clusters. The function should
#' produce a data frame that identifies the clusters in the treatment group under the new
#' randomisation scheme. The data frame can either have a single column with name cl_var or
#' two columns of cl_var and t identifying the cluster ID and time period a cluster joins
#' the treatment group. If NULL then clusters are randomised in a 1:1 ratio to treatment and control
#' @return A vector of the length of fitlist with the test statistics for each model and null
#' hypothesis
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
#' fitlist <- list(fit1,fit2)
#' nullfitlist <- list()
#' for(i in 1:length(fitlist)){
#'   nullfitlist[[i]] <- est_null_model(fitlist[[i]],
#'                                      data,
#'                                      tr_var = "treat",
#'                                      null_par = 0)
#' }
#' out <- lme_permute2(nullfitlist,
#'                data=data,
#'                cl_var = "cl")
#' @export
lme_permute2 <- function(fitlist,
                         data,
                         null_par=rep(0,length(fitlist)),
                         cl_var = "cl",
                         rand_func= NULL){
  if(!is(fitlist,"list"))stop("fitlist should be a list.")
  if(!is(data,"data.frame"))stop("data should be a data frame")
  if(length(null_par)!=length(fitlist))stop("length(null_par) should be the same as length(fitlist)")
  if(!cl_var%in%colnames(data))stop("cl_var not in colnames(data)")

  #generate a new randomisation scheme
  if(is.null(rand_func)){
    df_tr <- sample(unique(data[,cl_var]),ceiling(length(unique(data[,cl_var]))/2))#gen_rand_order(length(unique(data$cl)),length(unique(data$t)))
    data$treat_perm <- 0
    data[data[,cl_var] %in% df_tr, 'treat_perm'] <- 1
  } else {
    df_tr <- do.call(rand_func,list(nJ=length(unique(data[,cl_var])),
                                    nT=length(unique(data$t))))
    data$treat_perm <- 0
    if(!is(df_tr,"data.frame")){
      stop("rand_func must produce a data frame")
    } else {
      if("t" %in% colnames(df_tr)){
        for(i in 1:nrow(df_tr)){
          data[data[,cl_var]==df_tr$cl[i]&data$t>=df_tr$t[i],'treat_perm'] <- 1
        }
      } else {
        data[data[,cl_var] %in% df_tr[,cl_var], 'treat_perm'] <- 1
      }
    }

  }

  res <- rep(NA,length(fitlist))

  for(i in 1:length(fitlist)){
    res[i] <- qscore_stat(fitlist[[i]],
                          data,
                          null_par = null_par[i],
                          tr_assign = "treat_perm")
  }

  return(res)
}

#' Wrapper function to replicate permutations
#'
#' Calls lme_permute2 to replicate the permutations
#'
#' @param fitlist A list of glm model objects fitted under the null hypotheses
#' @param data A data frame containing the data used to fit the models in fitlist
#' @param n_permute Number of permutations to conduct
#' @param null_pars A vector of the same length as fitlist specifying the value(s) of the
#' treatment effect parameter(s) under the null hypotheses
#' @param cl_var String specifying the name of the column identifying the clusters/cluster-time
#' @param rand_func The name of a function that re-randomises the clusters. The function should
#' produce a data frame that identifies the clusters in the treatment group under the new
#' randomisation scheme. The data frame can either have a single column with name cl_var or
#' two columns of cl_var and t identifying the cluster ID and time period a cluster joins
#' the treatment group.  If NULL then clusters are randomised in a 1:1 ratio to treatment and control
#' @return An array of dimension length(fitlist)*n_permute containing the test statistics
#' for each model and each iteration
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
#' fitlist <- list(fit1,fit2)
#' nullfitlist <- list()
#' for(i in 1:length(fitlist)){
#'   nullfitlist[[i]] <- est_null_model(fitlist[[i]],
#'                                      data,
#'                                      tr_var = "treat",
#'                                      null_par = 0)
#' }
#' out <- permute(nullfitlist,
#'                data=data,
#'                n_permute = 10,
#'                cl_var = "cl")
#' @export
permute <- function(fitlist,
                    data,
                    n_permute=100,
                    null_pars=rep(0,length(fitlist)),
                    cl_var = "cl",
                    rand_func=NULL){
  if(!cl_var%in%colnames(data))stop("cl_var not in colnames(data)")
  res <- replicate(n_permute,lme_permute2(fitlist,
                                          data,
                                          null_par = null_pars,
                                          cl_var=cl_var,
                                          rand_func = rand_func))

  return(res)
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
  if(!(is(fit,"glm")|!is(fit,"lm")|is(fit,"glmerMod")|is(fit,"lmerMod")))
    stop("Model class should be glm, lm, or mer")

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

