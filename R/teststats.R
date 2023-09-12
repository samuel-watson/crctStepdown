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
  if(!(is(fit,"glmerMod")|is(fit,"lmerMod")|is(fit,"glm")|is(fit,"lm")|is(fit,"mcml")))stop("fit should be glm, lm, glmer, lmer, or mcml model object")
  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")

  if(!is(fit,"mcml")){
    if(is(fit,"glmerMod")){
      fixeff <- names(lme4::fixef(fit))
      family <- fit@resp$family
    } else if(is(fit,"lmerMod")){
      fixeff <- names(lme4::fixef(fit))
      family <- stats::gaussian()
    } else if(is(fit,"glm")){
      fixeff <- names(coef(fit))
      family <- stats::family(fit)
    } else if(is(fit,"lm")){
      fixeff <- names(coef(fit))
      family <- stats::gaussian()
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
    Y <- data[!is.na(data[,outv]),outv]
  } else {
    if(fit$family == "bernoulli"){
      family <- stats::binomial()
    } else {
      family <- do.call(fit$family,list())
    }
    form1 <- fit$mean_form
    form1 <-  gsub(tr_var,"",form1)
    if(substr(form1,1,1) == "+")form1 <- substr(form1,2,nchar(form1))
    form1 <-  gsub("\\+\\+","+",form1)
    form1 <-  gsub("\\+\\-","-",form1)
    form <- paste0("~",form1)
    Y <- fit$y
    outv <- "y"
  }

  off1 <- c(data[,tr_var]*null_par)
  X <- model.matrix(object=as.formula(form),data)

  if(is(fit,"glmerMod")|is(fit,"glm")|(is(fit,"mcml")&&family[[1]]!="gaussian")){
    #family <- fit@call[[4]]
    if(all(off1==0)){
      f1 <- fastglm::fastglm(X,Y,family=family,method=3)
    } else {
      f1 <- fastglm::fastglm(X,Y,family=family,offset=off1,method=3)
    }
  } else {
    #data[,outv] <- data[,outv]-off1
    Y <- Y - off1
    f1 <- RcppArmadillo::fastLm(X=X,y=Y)
  }

  f1$outv <- outv

  return(list(fit=f1,X=X,y=Y,family=family))
}
