#' Function to generate a stepped-wedge cRCT randomisation allocation
#'
#' Function to generate a stepped-wedge cRCT randomisation allocation. Assumes a baseline
#' and endline period in which no clusters and all clusters have the intervention,
#' respectively.
#'
#' @param nJ Number of clusters
#' @param nT Number of time points
#' @return A data frame with columns cl and t indicating the time pp
#' @examples
#   gen_rand_order(10,7)
#' @export
gen_rand_order <- function(nJ,nT){
  df_tr <- data.frame(cl=sample(1:nJ,nJ),t=NA)
  rem <- nJ %% (nT-2)
  nper <- floor(nJ/(nT-2))
  reps <- rep(c(nper,nper+1),
              c(nT-2-rem, rem))
  reps <- sample(reps,length(reps),replace = FALSE)
  df_tr$t <- rep(2:(nT-1),reps)
  return(df_tr)
}

#' Simulates data from a two-arm parallel cluster randomised trial
#'
#' Simple simulation of two Poisson distributed outcomes for a two-arm
#' parallel cluster randomised trial with no baseline measures. A log-linear model
#' is specified y~Poisson(lambda) with lambda = exp(mu + beta*D + theta) where D is the
#' treatment effect indicator equal to one in clusters with the treatment and zero
#' otherwise, and theta~N(0,sigma^2) is the cluster random effect.
#'
#' @param nJ Vector of two integers with the number of clusters in treatment and
#' control arms
#' @param N Number of individuals per cluster
#' @param mu Vector of two numeric values with the intercept terms for the two models on the log
#' scale
#' @param beta Vector of two numeric values that are the treatment effect parameters in the two models
#' @param sig_cl Vector of two values equal to the variance of the random effect in each model
#' @return A list consisting of: (1) data frame with the cluster IDs (cl), treatment effect indicators (treat),
#' and two outcomes (y1, y2), and (2) the values of the treatment effect parameters used in the simulation.
#' @importFrom stats rnorm rpois
#' @export
twoarm_sim <- function(nJ=c(7,7),
                       N = 20,
                       mu = rep(1,2),
                       beta = c(0,0),
                       sig_cl = rep(0.05,2)){

  tr <- sample(1:sum(nJ),ceiling(sum(nJ)/2))
  data <- data.frame(cl=rep(1:sum(nJ),each=N),
                     id = rep(1:N,sum(nJ)))
  data$treat <- 0
  data[data$cl %in% tr, 'treat'] <- 1
  alpha_1 <- rnorm(sum(nJ),0,sd=sqrt(sig_cl[1]))
  alpha_2 <- rnorm(sum(nJ),0,sd=sqrt(sig_cl[2]))
  data$y1 <- rpois(nrow(data),exp(mu[1] + beta[1]*data$treat + alpha_1[data$cl]))
  data$y2 <- rpois(nrow(data),exp(mu[2] + beta[2]*data$treat + alpha_2[data$cl]))

  return(list(data,beta))
}
