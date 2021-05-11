#' Inverse logit link function
#'
#' Inverse logit link function exp(x)/(1+exp(x))
#' @param x Numerical value
#' @export
h_logit <- function(x){
  return(exp(x)/(1+exp(x)))
}

#' Poisson exponential link function
#'
#' Exponential function
#' @param x Numerical value
#' @export
h_pois <- function(x){
  return(exp(x))
}

#' Identity link function
#'
#' Identity function
#' @param x Numerical value
#' @export
h_norm <- function(x){
  return(x)
}

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
#' @export
est_null_model <- function(fit,
                           data,
                           tr_var="treat",
                           null_par){
  if(!(is(fit,"glmer")|is(fit,"lmer")))stop("fit should be glmer or lmer model object")
  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")

  s1 <- summary(fit)
  fixeff <- rownames(s1$coefficients)
  fixeff <- fixeff[!fixeff%in%c(tr_var,"(Intercept)")]
  outv <- outname_fit(fit)
  call1 <- fit@call[[1]]

  if(length(fixeff)==0){
    form <- paste0(outv," ~ 1")
    #environment(form) <- environment()
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

  if(call1 == "glmer"){
    family <- fit@call[[4]]
    if(all(off1==0)){
      f1 <- glm(form,data,family=family)
    } else {
      f1 <- do.call("glm",list(formula=form,data=data,family=family,offset=off1))
    }
  } else if(call1=="lmer"){
    data[,outv] <- data[,outv]-off1
    f1 <- do.call("lm",list(formula=form,data=data))
  }

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
#' @return The value of the test statistic
#' @export
qscore_stat <- function(fit,
                        data,
                        null_par=0,
                        tr_var = "treat",
                        cl_var="cl",
                        tr_assign="treat"){
  if(!(is(fit,"glm")|is(fit,"lm")))stop("fit should be a glm or lm model")
  if(!is(data,"data.frame"))stop("Data should be a data frame")
  if(!tr_var%in%colnames(data))stop("tr_var not in colnames(data)")
  if(!cl_var%in%colnames(data))stop("cl_var not in colnames(data)")

  df <- data
  df[,tr_var] <- 0
  outv <- outname_fit(fit)
  if(class(fit)[1] == "glm"){
    pr1 <- predict(fit,data,type="link")
  } else {
    pr1 <- predict(fit,data)
    pr1 <- pr1 + null_par*data[,tr_var]
  }
  s1 <- summary(fit)
  call1 <- fit$call[[1]]#fit@call[[1]]
  #family <- nullfitlist[[1]]$family[[1]]
  if(call1=="glm"){
    family <- fit$family[[1]]#fit@call[[4]]
    if(family=="poisson"){
      f2 <- "h_pois"
    }
    if(family=="binomial"){
      f2 <- "h_logit"
    }
  } else {
    f2 <- "h_norm"
  }

  pr2 <- do.call(f2,list(pr1))
  resid <- matrix(data[,outv]-pr2,ncol=1)
  Tval <- ifelse(data[,tr_assign]==1,1,-1)
  score <- sum(Tval*resid, na.rm = TRUE)
  aggscore <- aggregate(Tval*resid,list(data[,cl_var]),sum, na.rm=TRUE)
  scoresq <- sum(aggscore$V1^2, na.rm = TRUE)
  return(score/sqrt(scoresq))
}
