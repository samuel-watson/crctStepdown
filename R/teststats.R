#' Estimates null model
#'
#' Given an lme4 model object and the value of the treatment effect parameter under the
#' null hypothesis, the function returns a glm or lm object fitted under the null model
#' with no cluster effects. For linear models (lmer) the value of the null is subtracted
#' from the value of the outcome for those in receipt of the treatment and an lm model
#' is fitted with no treatment effect. For generalised linear models (glmer) the model is
#' refitted as a glm model with the treatment effect specified as an offset.
#'
#' @param fit An lme4 model object
#' @param data The data frame used to fit model fit
#' @param tr_var A string indicating the name of the column in data that is a binary indicator
#' for whether the observation was under the treatment (1=treatment, 0=control)
#' @param null_par Numeric the value of tr_var parameter under the null hypothesis
#' @return An lm or glm model fit under the null model
#' @importFrom methods is
#' @importFrom stats model.matrix as.formula coef
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
#' @export
est_null_model <- function(fit,
                           data,
                           tr_var="treat",
                           null_par){
  if(!(is(fit,"glmerMod")|is(fit,"lmerMod")|is(fit,"glm")|is(fit,"lm")))stop("fit should be glm, lm, glmer, or lmer model object")
  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")

  type <- ifelse(is(fit,"glmerMod")|is(fit,"lmerMod"),"mer","glm")

  if(type=="mer"){
    fixeff <- names(lme4::fixef(fit))
    family <- fit@resp$family
  } else if(type=="glm"){
    fixeff <- names(coef(fit))
    family <- stats::family(fit)
  }

  fixeff <- fixeff[!fixeff%in%c(tr_var,"(Intercept)")]
  outv <- outname_fit(fit)
  #call1 <- fit@call[[1]]

  if(length(fixeff)==0){
    form <- paste0(outv," ~ 1")
  }
  if(length(fixeff)==1){
    form <- paste0(outv," ~ ",fixeff)
  }
  if(length(fixeff)>1){
    form1 <- paste0(fixeff, sep=" + ", collapse = "")
    form1 <- stringr::str_sub(form1,1,nchar(form1)-3)
    form <- paste0(outv," ~ ",form1)
  }

  off1 <- c(data[,tr_var]*null_par)

  X <- model.matrix(object=as.formula(form),data)
  Y <- data[!is.na(data[,outv]),outv]

  if(is(fit,"glmerMod")|is(fit,"glm")){
    #family <- fit@call[[4]]
    if(all(off1==0)){
      f1 <- fastglm::fastglm(X,Y,family=family,method=3)
    } else {
      f1 <- fastglm::fastglm(X,Y,family=family,offset=off1,method=3)
    }
  } else if(is(fit,"lmerMod")|is(fit,"lm")){
    data[,outv] <- data[,outv]-off1
    f1 <- RcppArmadillo::fastLm(X=X,y=Y)
  }

  f1$outv <- outv

  return(f1)
}

#' Calculates the randomisation test statistic
#'
#' Calculates a randomisation test statistic based on the sum of generalised, studentized
#' residuals for a given model fit and null hypothesis of the value of the treatment effect
#' parameter
#' @param fit A glm or lm model fit
#' @param data The data frame used to fit model fit
#' @param null_par Numeric the value of tr_var parameter under the null hypothesis
#' @param tr_var String indicating the name of the column in data that is a binary indicator
#' for whether the observation was under the treatment (1=treatment, 0=control)
#' @param cl_var String specifying the name of the column identifying the clusters/cluster-time
#' @param tr_assign String specifying the treatment assignment under a new permutation. If calculating
#' the test statistics under the original treatment assignment then this should be the same
#' as tr_var
#' @param inv_sigma optional, inverse of the covariance matrices of the observations. If provided then the weighted q-score statistic
#' is calculated
#' @return The value of the test statistic
#' @importFrom methods is
#' @importFrom stats predict.lm predict.glm family
#' @examples
#' out <- twoarm_sim()
#' data <- out[[1]]
#'   fit1 <- glm(y1 ~ treat ,
#'               data=data,
#'               family="poisson")
#' qscore_stat(fit=fit1,
#'             data=data)
#' @export
qscore_stat <- function(fit,
                        data,
                        null_par=0,
                        tr_var = "treat",
                        cl_var="cl",
                        tr_assign="treat",
                        inv_sigma = NULL){
  #if(!(is(fit,"glm")|is(fit,"lm")))stop("fit should be a glm or lm model")
  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")
  if(!cl_var%in%colnames(data))stop("cl_var not in colnames(data)")

  df <- data
  df[,tr_var] <- 0
  outv <- outname_fit(fit)
  if(class(fit)[1] == "glm"){
    pr1 <- predict.glm(fit,data,type="response")
    family <- family(fit)
    xb <- predict.glm(fit,data,type="link")
  } else if(is(fit,"fastglm")|is(fit,"fastLm")){
    pr1 <- fit$fitted.values
    family <- fit$family
    xb <- fit$linear.predictors
  } else {
    pr1 <- predict.lm(fit,data)
    family <- family(fit)
    xb <- pr1
  }

  if(is(fit,"fastLm")|is(fit,"lm")){
    pr1 <- pr1 + null_par*data[!is.na(data[,outv]),tr_var]
  }

  cltmp <- as.numeric(as.factor(data[,cl_var]))-1

  if(is.null(inv_sigma)){
    sc <- qscore(y=data[!is.na(data[,outv]),outv],
                 x=pr1,
                 Tr=data[!is.na(data[,outv]),tr_assign],
                 cl=cltmp,
                 ncl=length(unique(cltmp)))
  } else {
    tr <- data[!is.na(data[,outv]),tr_assign]
    tr[tr==0] <- -1
    g <- get_G(xb,family)
    sc <- qscorew(y=data[!is.na(data[,outv]),outv],
                 x=pr1,
                 Tr=diag(tr),
                 g = g,
                 sigma = inv_sigma,
                 cl=cltmp,
                 ncl=length(unique(cltmp)))
  }

  return(sc)
}

#' Get inverse partial derivative of link function
#'
#' Returns the inverse partial derivative of the link function with respect to the linear predictor
#'
#' @param x vector of predicted model values on the scale of the linear predictor
#' @param family a family object, see \link[stats]{family}. Only identity, log, logit, and probit link functions currently supported.
#' @return A vector of values of the inverse of the partial derivative of the link function with respect to the linear predictor
#' @examples
#' get_G(seq(-1,1,length.out=10),family=binomial())
#' @importFrom stats family dnorm
#' @export
get_G <- function(x,
                  family){

  if(!family[[2]]%in%c("identity","log","logit","probit"))stop("only identity, log, logit, and probit link functions currently supported for
                                                               analyses using covariance matrix")

  if(family[[2]] == "identity"){
    dx <- rep(1,length(x))
  } else if(family[[2]] == "log"){
    dx <- exp(x)
  } else if(family[[2]] == "logit"){
    dx <- (exp(x)/(1+exp(x)))*(1-exp(x)/(1+exp(x)))
  } else if(family[[2]] == "probit"){
    dx <- -1/dnorm(x)
  }
  return(dx)
}
