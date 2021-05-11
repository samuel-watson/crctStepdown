test_that("whole package", {
  #generate sim data
  data <- twoarm_sim()
  delta <- data[[2]]
  data <- data[[1]]

  fit1 <- lme4::glmer(y1 ~ treat + (1|cl) ,
                data=data,
                family="poisson")

  fit2 <- lme4::glmer(y2 ~ treat + (1|cl),
                data=data,
                family="poisson")

  res2 <- stepdown(fitlist=list(fit1,fit2),
                            data=data,
                            n_permute = 100,
                            nsteps=100,
                            plots=FALSE,
                            verbose=TRUE)
  pv <- res2$p_value
  lc <- res2$lower_ci

  expect_s3_class(res2,"data.frame")
  expect_type(pv,"double")
  expect_type(lc,"double")
})

