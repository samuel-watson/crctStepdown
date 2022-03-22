# crctStepdown

Permutation test inference for generalised linear mixed models including corrections for multiple testing. The package was originally named after the permutation
based stepdown procedure of Romano and Wolf (2005), but has expanded to include Bonferroni, Holm's stepdown method, and no correction at all. The package provides
functions to estimate p-values and confidence intervals using models fits from either `lme4`, `glm`, or `lm` in R. The permutational test statistic is either unweigted
or weighted by the covariance matrix of observations as described by Braun and Feng (2001). Confidence intervals are estimated using an iterative search procedure. All 
the methods are described in Watson et al (2021).

## Usage
As an example of usage the user can fit two models using `lme4`. The models can be of different families.
```
fit1 <- lme4::lmer(y1 ~ treat + x1 + x2, data=data)
fit2 <- lme4::glmer(y2 ~ treat + x1 + x2, data= data)
```

Then these can be passed to the stepdown function:
```
stepdown(fitlist=list(fit1,fit2),
          data=data,
          n_permute = 1000,
          nsteps=1000,
          plots=TRUE,
          verbose=TRUE,
          type = "rw")
```
where type can be `rw` for Romano-Wolf, `h` or `hr` for standard Holm and Holm using permutation tests respectively, `b` or `br` for standard Bonferroni or Bonferroni
using permutation tests, and `none` for no correction for multiple testing. To use the weighted statistic the user can provide a covariance matrix to the 
argument `sigma`. 

## References
Braun and Feng 2001. [Optimal Permutation Tests for the Analysis of Group Randomized Trials.](https://www.tandfonline.com/doi/abs/10.1198/016214501753382336) _Journal of the American Statistical Association_ 96(456):1424-1432

Romano and Wolf 2005. [Exact and approximate stepdown methods for multiple hypothesis testing.](https://www.tandfonline.com/doi/abs/10.1198/016214504000000539) _Journal of the American Statistical Association_ 100: 94–108

Watson, Akinyemi, and Hemming. 2021 [Multiple testing corrections for p-values and confidence intervals from generalised linear mixed models.](https://arxiv.org/abs/2107.10017) arXiv 2107.10017v2

