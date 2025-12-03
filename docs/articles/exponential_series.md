# Masked Data Likelihood Model: Components with Exponentially Distributed Lifetimes Arranged In Series Configuration

The R package `likelihood.md.series.system` is a framework for
estimating the parameters of latent component lifetimes from *masked
data* in a series system.

## Exponentially Distributed Component Lifetimes

Consider a series system in which the components have exponentially
distributed lifetimes. The $j$ component of the $i$ has a lifetime
distribution given by
$$T_{ij} \sim \operatorname{EXP}\left( \lambda_{j} \right)$$ for
$j = 1,\ldots,m$. Thus, $\lambda = (\lambda_{1},\ldots,\lambda_{m})$.
The random variable $T_{ij}$ has a reliability, pdf, and hazard function
given respectively by $$\begin{aligned}
{R_{j}\left( t|\lambda_{j} \right)} & {= \exp\left( - \lambda_{j}t \right),} \\
{f_{j}\left( t|\lambda_{j} \right)} & {= \lambda_{j}\exp\left( - \lambda_{j}t \right),} \\
{h_{j}\left( \cdot |\lambda_{j} \right)} & {= \lambda_{j}}
\end{aligned}$$ where $t > 0$ is the lifetime and $\lambda_{j} > 0$ is
the failure rate of the $j$-th component.

The lifetime of the series system composed of $m$ components with
exponentially distributed lifetimes has a reliability function given by
$$R_{T_{i}}\left( t_{i}|{\mathbf{λ}} \right) = \exp( - \sum\limits_{j = 1}^{m}\lambda_{j}t_{i})$$
where $t > 0$. A series system with exponentially distributed lifetimes
is also exponentially distributed.

The series system’s failure rate function is given by
$$h\left( \cdot |{\mathbf{λ}} \right) = \sum\limits_{j = 1}^{m}\lambda_{j}$$
whose proof follows from Theorem .

We see that the failure rate ${\mathbf{λ}} = \sum_{j = 1}^{n}\lambda$ is
*constant*, consistent with the the exponential distribution being the
only continuous distribution that has a constant failure rate.

The pdf of the series system is given by
$$f_{T_{i}}\left( t_{i}|{\mathbf{λ}} \right) = (\sum\limits_{j = 1}^{m}\lambda_{j})\exp( - \sum\limits_{j = 1}^{m}\lambda_{j}t_{i})$$
where $t_{i} > 0$ is the lifetime of the system. The conditional
probability that component $k$ is the cause of a system failure at time
$t$ is given by
$$f_{K_{i}|T_{i}}\left( k|t,{\mathbf{λ}} \right) = f_{K_{i}}\left( k|{\mathbf{λ}} \right) = \frac{\lambda_{k}}{\sum\limits_{p = 1}^{m}\lambda_{p}}$$
where $k \in \{ 1,\ldots,m\}$ and $t > 0$. Due to the constant failure
rates of the components, $K_{i}$ and $T_{i}$ are mutually independent.
The joint pdf of $K_{i}$ and $T_{i}$ is given by
$$f_{K_{i},T_{i}}\left( k,t|{\mathbf{λ}} \right) = \lambda_{k}\exp( - \sum\limits_{j = 1}^{m}\lambda_{j}t)$$
where $k \in \{ 1,\ldots,m\}$ and $t > 0$.

## Likelihood Model

In this study, the system is a series system with $m$ components. The
true DGP for the system lifetime is in the exponential series system
family, i.e., the component lifetimes are exponentially and
independently distributed and we denote the true parameter value by
$\theta$.

The principle object of study is \$, which in the case of the
exponential series system family consists of $m$ rate (scale) parameters
for each component lifetime,
${\mathbf{λ}} = \left( \lambda_{1},\ldots,\lambda_{m} \right)\prime$.

We are interested in estimating the $\theta$ from masked data. The
masking comes in two independent forms:

- Censored system failure times, e.g., right-censoring. The system
  failure time is the minimum of the component lifetimes, and it is
  right-censored if the system failure time does not occur during the
  observation period, $$T_{i} = \min\{\tau_{i},T_{i1},\ldots,T_{im}\},$$
  where $\tau_{i}$ is the right-censoring time for the $i$ observation
  and $T_{i1},\ldots,T_{im}$ are the component lifetimes for the $i$th
  system.

- The cause of failure, the failed component, is masked. This masking
  comes in the form of a candidate set $\mathcal{C}_{i}$ that, on
  average, conveys information about the component cause of failure.

The candidate set $\mathcal{C}_{i}$ is a random variable that is a
subset of $\{ 1,\ldots,m\}$. The true DGP for the candidate set model
has a general form that may be denoted by
$$\Pr\{\mathcal{C}_{i} = c_{i}|T_{1} = j,\ldots,T_{m},\theta,\text{other factors}\}.$$

This is a pretty complicated looking model, and we are not even
interested in the DGP for candidate sets, except to the extent that it
affects the sampling distribution of the MLE for $\theta$.

In theory, given some candidate set model, we could construct a joint
likelihood function for the full model and jointly estimate the
parameters of both the candidate set model and $\theta$. In practice,
however, this could be a very challenging task unless we make some
simplifying assumptions about the DGP for candidate sets.

### Candidate set models

In every model we consider, we assume that the candidate set
$\mathcal{C}_{i}$ is only a function of the component lifetimes
$T_{i1},\ldots,T_{im}$, $\theta$, and the right-censoring time
$\tau_{i}$. That is, the candidate set $\mathcal{C}_{i}$ is independent
of any other factors (or held constant for the duration of the
experiment), like ambient temperature, and these factors also have a
neglible effect on the series system lifetime and thus we can ignore
them.

#### Reduced likelihood model

In the Bernoulli candidate set model, we make the following assumptions
about how candidate sets are generated:

- $C_{1}$: The index of the failed component is in the candidate set,
  i.e., $\Pr\{ K_{i} \in \mathcal{C}_{i}\} = 1$, where
  $K_{i} = \arg\min_{j}\{ T_{ij}:j = 1,\ldots,m\}$.

- $C_{2}$: The probability of $C_{i}$ given $K_{i}$ and $T_{i}$ is
  equally probable when the failed component varies over the components
  in the candidate set, i.e.,
  $\Pr\{\mathcal{C}_{i} = c_{i}\left| K_{i} = j,T_{i} = t_{i},\theta\} = \Pr\{ C_{i} = c_{i} \right|K_{i} = j\prime,T_{i} = t_{i}\}$
  for any $j,j\prime \in c_{i}$.

- $C_{3}$: The masking probabilities are conditionally independent of
  $\theta$ given $K_{i}$ and $T_{i}$, i.e.,
  $\Pr\{\mathcal{C}_{i} = c_{i}|K_{i} = j,T_{i} = t_{i}\}$ is not a
  function of $\theta$.

Using these simplifying assumptions, we can arrive at a reduced
likelihood function that only depends on $\theta$ and the observed data
and as long as our candidate set satisfies conditions $C_{1}$, $C_{2}$,
and $C_{3}$, our reduced likelihood function obtains the same MLEs as
the full likelihood function.

We see that
$$\Pr\{\mathcal{C}_{i} = c_{i},|K_{i} = j,T_{i} = t_{i}\} = g\left( c_{i},t_{i} \right),$$
since the probability cannot depend on $j$ by condition $C_{2}$ and
cannot depend on $\theta$ by condition $C_{3}$. Thus, we can write the
likelihood function as
$$L(\theta) = \prod\limits_{i = 1}^{n}f_{T_{i}}\left( t_{i}|\theta \right)g\left( c_{i},t_{i} \right).$$

We show that $g\left( c_{i},t_{i} \right)$ is proportional to
$$g\left( c_{i},t_{i} \right) \propto \sum\limits_{j \in c_{i}}f_{j}\left( t_{i}|\theta_{j} \right)\prod\limits_{l = j,l \neq j}^{m}R_{l}\left( t_{i}|\theta_{l} \right),$$
and thus

Note, however, that different ways in which the conditions are met will
yield MLEs with different sampling distributions, e.g., more or less
efficient estimators.

#### Bernoulli candidate set model \#1

This is a special case of the reduced likelihood model. In this model,
we satisfy conditions $C_{1}$, $C_{2}$, and $C_{3}$, but we include each
of the non-failed components with a fixed probability $p$, $0 < p < 1$.

In the simplest case, $p = 0.5$, and candidate set $c_{i}$ has a
probability given by
$$\Pr\{\mathcal{C}_{i} = c_{i}|K_{i} = j,T_{i} = t_{i}\} = \begin{cases}
(1/2)^{m - 1} & {{\text{if}\mspace{6mu}}{j \in c_{i}}{\mspace{6mu}\text{and}\mspace{6mu}}{c_{i} \subseteq \{ 1,\ldots,m\}}} \\
0 & {{{\text{if}\mspace{6mu}}{j \notin c_{i}}}.}
\end{cases}$$

#### Bernoulli candidate set model \#2

Now, we remove condition $C_{2}$. We still assume conditions $C_{1}$ and
$C_{3}$, but we allow $C_{i}$ to depend on the failed component $K_{i}$,
i.e.,
$$\Pr\{\mathcal{C}_{i} = c_{i}\left| K_{i} = j,T_{i} = t_{i},\theta\} \neq \Pr\{ C_{i} = c_{i} \right|K_{i} = j\prime,T_{i} = t_{i}\}$$
for $j,j\prime \in c_{i}$.

In this case, we can write the likelihood function as
$$L(\theta) = \prod\limits_{i = 1}^{n}f_{T_{i}}(t_{i}\left| \theta)\prod\limits_{j = 1}^{m}\Pr\{ K_{i} = j \right|T_{i} = t_{i}\}\prod\limits_{c_{i} \in \mathcal{C}_{i}}g\left( c_{i},t_{i},j \right).$$

## Simulation

The most straightforward series system to estimate is the series system
with exponentially distributed component lifetimes.

Suppose an exponential series system with $m$ components is
parameterized by the following R code:

``` r
theta <- c(1,     # component 1 failure rate
           1.1,   # 3
           0.95,  # 5
           1.15,  # 6
           1.1)   # 7

m <- length(theta)
```

So, in our study, $\theta = (1,1.1,0.95,1.15,1.1)\prime$. The component
assigned to index $j$ has an exponentially distributed lifetime with a
failure rate $\theta_{j}$, e.g., $\theta_{2} = 1.1$ is the failure rate
of the component indexed by $2$.

Let’s simulate generating the lifetimes of the $m = 5$ components for
this series system:

``` r
set.seed(7231) # set seed for reproducibility
n <- 7500
comp_times <- matrix(nrow=n,ncol=m)
for (j in 1:m)
    comp_times[,j] <- rexp(n,theta[j])
comp_times <- md_encode_matrix(comp_times,"t")
print(comp_times,n=4)
#> # A tibble: 7,500 × 5
#>      t1     t2    t3     t4    t5
#>   <dbl>  <dbl> <dbl>  <dbl> <dbl>
#> 1  2.95 1.67   0.777 0.0476 0.380
#> 2  1.10 2.20   2.59  3.27   0.626
#> 3  1.20 0.0956 3.74  1.80   2.13 
#> 4  1.44 0.0460 0.123 0.522  0.102
#> # ℹ 7,496 more rows
```

Next, we use the function `md_series_lifetime_right_censoring` to
decorate the masked data with the right-censor-censoring time chosen by
the probability $\Pr\{ T_{i} > \tau\} = 0.75$:

``` r
q <- 0.25
tau <- rep(-(1/sum(theta))*log(q),n)
data <- comp_times %>% md_series_lifetime_right_censoring(tau)
print(data,n=4,drop_latent=TRUE)
#> Latent variables:  t1 t2 t3 t4 t5 
#> # A tibble: 7,500 × 2
#>        t delta
#>    <dbl> <lgl>
#> 1 0.0476 TRUE 
#> 2 0.262  FALSE
#> 3 0.0956 TRUE 
#> 4 0.0460 TRUE 
#> # ℹ 7,496 more rows
```

### Masked component cause of failure

We simulate candidate sets using the Bernoulli candidate model with an
appropriate set of parameters to satisfy conditions $C_{1}$, $C_{2}$,
and $C_{3}$:

``` r
p <- .3
data <- data %>% md_bernoulli_cand_c1_c2_c3(p)
print(data[,paste0("q",1:m)],n=4)
#> Latent variables:  q1 q2 q3 q4 q5 t1 t2 t3 t4 t5 
#> # A tibble: 7,500 × 5
#>      q1    q2    q3    q4    q5
#>   <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1   0.3   0.3   0.3   1     0.3
#> 2   0     0     0     0     0  
#> 3   0.3   1     0.3   0.3   0.3
#> 4   0.3   1     0.3   0.3   0.3
#> # ℹ 7,496 more rows
```

Now, to generate candidate sets, we sample from these probabilities:

``` r
data <- data %>% md_cand_sampler()
print(md_boolean_matrix_to_charsets(data,drop_set=TRUE),drop_latent=TRUE,n=6)
#> Latent variables:  q1 q2 q3 q4 q5 t1 t2 t3 t4 t5 
#> # A tibble: 7,500 × 3
#>        t delta x        
#>    <dbl> <lgl> <chr>    
#> 1 0.0476 TRUE  {1, 4}   
#> 2 0.262  FALSE {}       
#> 3 0.0956 TRUE  {1, 2, 5}
#> 4 0.0460 TRUE  {1, 2, 4}
#> 5 0.0290 TRUE  {4}      
#> 6 0.160  TRUE  {3, 5}   
#> # ℹ 7,494 more rows
```

We see that after dropping latent (unobserved) columns, we only have the
right censoring time, right censoring indicator, and the candidate sets.
(Note that this time we showed the candidate sets in a more friendly way
using `md_boolean_matrix_to_charsets`.)

## Likelihood Model

The likelihood model is a statistical model that describes the
distribution of the observed data as a function of the parameters of
interest.

We construct a likelihood model for the masked data model with
exponentially distributed component lifetimes with the following code:

``` r
model <- exp_series_md_c1_c2_c3()
```

## Maximum likelihood estimation

The log-likelihood for our masked data model when we assume Conditions ,
, and is given by
$$\ell(\lambda) = \sum\limits_{i = 1}^{n}\left( 1 - \delta_{i} \right)\log(\sum\limits_{j \in c_{i}}\lambda_{j}) - (\sum\limits_{i = 1}^{n}s_{i})(\sum\limits_{j = 1}^{m}\lambda_{j}).$$

The set of solutions to the MLE equations must be stationary points,
i.e., a point at which the score function of type
$\left. {\mathbb{R}}^{m}\mapsto{\mathbb{R}}^{m} \right.$ is zero. The
$j$-th component of the output of the score function is given by
$$\frac{\partial\ell}{\partial\lambda_{p}} = \sum\limits_{i = 1}^{n}(\sum\limits_{j \in c_{i}}\lambda_{j})^{- 1}1_{\{ p \in c_{i}{\mspace{6mu}\text{and}\mspace{6mu}}\delta_{i} = 0\}} - \sum\limits_{i = 1}^{n}s_{i}.$$

We may find an MLE by solving the maximum likelihood equation , i.e., a
set of (stationary) points satisfying
$$\frac{\partial\ell}{\partial\lambda_{j}}|_{{\widehat{\lambda}}_{j}} = 0$$
for $j = 1,\ldots,m$. We approximate a solution to this problem by using
the iterative Newton-Raphson method as described in Section .

The Newton-Raphson method needs the observed information matrix, which
is a function of $\lambda$ of type
$\left. {\mathbb{R}}^{m}\mapsto{\mathbb{R}}^{m \times m} \right.$. The
$(j,k)$-th element of $J(\lambda)$ is given by
$$\frac{\partial^{2}\ell}{\partial\lambda_{j}\partial\lambda_{k}} = \sum\limits_{i = 1}^{n}(\sum\limits_{j \in c_{i}}\lambda_{j})^{- 2}1_{\{ j,k \in c_{i}{\mspace{6mu}\text{and}\mspace{6mu}}\delta_{i} = 0\}}.$$

### Log-likelihood of $\theta$ given masked data

The reduced log-likelihood function (the log of the kernel of the
likelihood function) is given by
$$\ell\left( \theta|\text{data} \right) = - \left( \sum\limits_{i = 1}^{n}t_{i} \right)\left( \sum\limits_{j = 1}^{m}\theta_{j} \right) + \sum\limits_{i = 1}^{n}\left( 1 - \delta_{i} \right)\log\left( \sum\limits_{j \in c_{i}}\theta_{j} \right).$$

We compute the log-likelihood function using generic dispatch:

``` r
ll <- loglik(model)
```

The returned function `ll(df, par)` evaluates the log-likelihood at
parameter `par` given data `df`. For example, at the true parameter
value:

``` r
ll(data, theta)
#> [1] -1485
```

Note that the implementation uses minimally sufficient statistics, which
improves computational efficiency.

The log-likelihood function contains the maximum amount of information
about parameter $\theta$ given the sample of masked data `data`
satisfying conditions $C_{1}$, $C_{2}$, and $C_{3}$.

Suppose we do not know that $\theta = (1,1.1,0.95,1.15,1.1)\prime$. With
the log-likelihood, we may estimate $\theta$ with $\widehat{\theta}$ by
solving
$$\widehat{\theta} = \operatorname{argmax}_{\theta \in \Omega}\ell(\theta),$$
i.e., finding the point that *maximizes* the log-likelihood on the
observed sample `data`. This is known as *maximum likelihood estimation*
(MLE). We typically solve for the MLE by solving
$$\nabla\ell|_{\theta = \widehat{\theta}} = 0.$$

A popular choice is gradient ascent, which is an iterative method based
on the update rule
$$\theta^{(n + 1)} = \theta^{n} + \eta\ell\left( \theta^{n} \right),$$
where $\eta$ is the learning rate.

We can also obtain the score (gradient) function via generic dispatch:

``` r
grad <- score(model)
```

The score at the true parameter should be close to zero (at the MLE, it
is exactly zero):

``` r
grad(data, theta)
#> [1] -28.488  12.506  -9.824 -24.504 -12.322
```

The `likelihood.model` framework provides analytical score and Hessian
implementations when available, falling back to numerical
differentiation otherwise.

In what follows, we use
[algebraic.mle](https://github.com/queelius/algebraic.mle) to help solve
the MLE equations and display various properties of the solution.

To solve the MLE equation, we use the generic
[`fit()`](https://generics.r-lib.org/reference/fit.html) function, which
dispatches to `fit.likelihood_model` for any object with
`"likelihood_model"` in its class. The
[`fit()`](https://generics.r-lib.org/reference/fit.html) function
returns a solver that uses `optim` internally with the BFGS method by
default.

``` r
# Get the solver from the model using generic dispatch
solver <- fit(model)

# Solve for MLE with initial guess
theta0 <- rep(1, m)
estimate <- solver(data, par = theta0, method = "Nelder-Mead")
```

The result is an `mle` object from the `algebraic.mle` package with rich
accessor methods:

``` r
# Print summary with confidence intervals
print(estimate)
#> Maximum likelihood estimator of type mle_likelihood_model is normally distributed.
#> The estimates of the parameters are given by:
#> [1] 0.9545 1.1430 0.9440 1.1085 1.0881
#> The standard error is  0.04314 0.04496 0.04247 0.04481 0.04465 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>          2.5% 97.5%
#> param1 0.8700 1.039
#> param2 1.0549 1.231
#> param3 0.8607 1.027
#> param4 1.0207 1.196
#> param5 1.0005 1.176
#> The MSE of the individual components in a multivariate estimator is:
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  0.0018608 -0.0002202 -0.0002258 -0.0002585 -0.0002654
#> [2,] -0.0002202  0.0020211 -0.0002331 -0.0002648 -0.0002360
#> [3,] -0.0002258 -0.0002331  0.0018034 -0.0002185 -0.0002453
#> [4,] -0.0002585 -0.0002648 -0.0002185  0.0020078 -0.0002316
#> [5,] -0.0002654 -0.0002360 -0.0002453 -0.0002316  0.0019936
#> The log-likelihood is  -1483 .
#> The AIC is  2976 .
```

We can access specific components of the MLE:

``` r
# Point estimate
theta.hat <- estimate$theta.hat
cat("MLE:", round(theta.hat, 4), "\n")
#> MLE: 0.9545 1.143 0.944 1.109 1.088

# Standard errors
cat("SE:", round(estimate$sigma, 4), "\n")
#> SE: 0.0019 -2e-04 -2e-04 -3e-04 -3e-04 -2e-04 0.002 -2e-04 -3e-04 -2e-04 -2e-04 -2e-04 0.0018 -2e-04 -2e-04 -3e-04 -3e-04 -2e-04 0.002 -2e-04 -3e-04 -2e-04 -2e-04 -2e-04 0.002

# Log-likelihood at MLE
cat("Log-likelihood:", round(estimate$loglike, 4), "\n")
#> Log-likelihood: -1483
```

Recall that the true parameter is
$\theta = (1,1.1,0.95,1.15,1.1)\prime$.

Due to sampling variability, different runs of the experiment will
result in different outcomes, i.e., $\widehat{\theta}$ has a sampling
distribution. We see that $\widehat{\theta} \neq \theta$, but it is
reasonably close. We may measure this sampling variability using the
variance-covariance matrix, bias, mean squared error (MSE), and
confidence intervals.

### Monte Carlo simulation study

To understand the sampling properties of the MLE, we conduct a Monte
Carlo simulation study. We repeatedly:

1.  Generate masked data from the true model
2.  Fit the likelihood model to obtain $\widehat{\theta}$
3.  Compute asymptotic confidence intervals

From the collection of estimates, we can compute:

- **Bias**:
  $\text{E}\left\lbrack \widehat{\theta} \right\rbrack - \theta$
- **Variance**: $\text{Var}\left\lbrack \widehat{\theta} \right\rbrack$
- **MSE**:
  $\text{E}\left\lbrack \left( \widehat{\theta} - \theta \right)^{2} \right\rbrack = \text{Bias}^{2} + \text{Var}$
- **Coverage probability**: Proportion of CIs containing the true value

``` r
set.seed(7231)

# Simulation parameters
B <- 200       # Number of Monte Carlo replications
alpha <- 0.05  # Significance level for CIs

# Storage for results
estimates <- matrix(NA, nrow = B, ncol = m)
se_estimates <- matrix(NA, nrow = B, ncol = m)
ci_lower <- matrix(NA, nrow = B, ncol = m)
ci_upper <- matrix(NA, nrow = B, ncol = m)
converged <- logical(B)
```

``` r
for (b in 1:B) {
    # Generate component lifetimes
    comp_times_b <- matrix(nrow = n, ncol = m)
    for (j in 1:m) {
        comp_times_b[, j] <- rexp(n, theta[j])
    }
    comp_times_b <- md_encode_matrix(comp_times_b, "t")

    # Apply masking pipeline
    data_b <- comp_times_b %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_bernoulli_cand_c1_c2_c3(p) %>%
        md_cand_sampler()

    # Fit model
    tryCatch({
        result_b <- solver(data_b, par = theta0, method = "Nelder-Mead")
        estimates[b, ] <- result_b$theta.hat
        se_estimates[b, ] <- sqrt(diag(result_b$sigma))

        # Asymptotic Wald CIs
        z <- qnorm(1 - alpha/2)
        ci_lower[b, ] <- result_b$theta.hat - z * sqrt(diag(result_b$sigma))
        ci_upper[b, ] <- result_b$theta.hat + z * sqrt(diag(result_b$sigma))
        converged[b] <- result_b$converged
    }, error = function(e) {
        converged[b] <<- FALSE
    })
}

cat("Convergence rate:", mean(converged, na.rm = TRUE), "\n")
#> Convergence rate: 0.99
```

#### Bias, Variance, and MSE

``` r
# Compute summary statistics (only for converged runs)
valid <- converged & !is.na(estimates[, 1])
est_valid <- estimates[valid, , drop = FALSE]

bias <- colMeans(est_valid) - theta
variance <- apply(est_valid, 2, var)
mse <- bias^2 + variance
rmse <- sqrt(mse)

# Create results table
results_df <- data.frame(
    Component = 1:m,
    True = theta,
    Mean_Est = colMeans(est_valid),
    Bias = bias,
    Variance = variance,
    MSE = mse,
    RMSE = rmse,
    Rel_Bias_Pct = 100 * bias / theta
)

knitr::kable(results_df, digits = 4,
             caption = "Monte Carlo Results: Bias, Variance, and MSE",
             col.names = c("Component", "True θ", "Mean θ̂", "Bias",
                          "Variance", "MSE", "RMSE", "Rel. Bias %"))
```

| Component | True θ | Mean θ̂ |    Bias | Variance |    MSE |   RMSE | Rel. Bias % |
|----------:|-------:|-------:|--------:|---------:|-------:|-------:|------------:|
|         1 |   1.00 | 1.0026 |  0.0026 |   0.0020 | 0.0020 | 0.0453 |      0.2627 |
|         2 |   1.10 | 1.0972 | -0.0028 |   0.0023 | 0.0023 | 0.0484 |     -0.2516 |
|         3 |   0.95 | 0.9461 | -0.0039 |   0.0017 | 0.0017 | 0.0414 |     -0.4146 |
|         4 |   1.15 | 1.1516 |  0.0016 |   0.0019 | 0.0019 | 0.0437 |      0.1421 |
|         5 |   1.10 | 1.0967 | -0.0033 |   0.0023 | 0.0023 | 0.0483 |     -0.3024 |

Monte Carlo Results: Bias, Variance, and MSE

#### Confidence Interval Coverage

The asymptotic $(1 - \alpha)$% Wald confidence interval is:
$${\widehat{\theta}}_{j} \pm z_{1 - \alpha/2} \cdot \text{SE}\left( {\widehat{\theta}}_{j} \right)$$

We assess coverage probability: the proportion of intervals that contain
the true parameter value.

``` r
# Compute coverage for each component
coverage <- numeric(m)
for (j in 1:m) {
    valid_j <- valid & !is.na(ci_lower[, j]) & !is.na(ci_upper[, j])
    covered <- (ci_lower[valid_j, j] <= theta[j]) & (theta[j] <= ci_upper[valid_j, j])
    coverage[j] <- mean(covered)
}

# Mean CI width
mean_width <- colMeans(ci_upper[valid, ] - ci_lower[valid, ], na.rm = TRUE)

coverage_df <- data.frame(
    Component = 1:m,
    True = theta,
    Coverage = coverage,
    Nominal = 1 - alpha,
    Mean_Width = mean_width
)

knitr::kable(coverage_df, digits = 4,
             caption = paste0("Coverage Probability of ", 100*(1-alpha), "% Confidence Intervals"),
             col.names = c("Component", "True θ", "Coverage", "Nominal", "Mean Width"))
```

| Component | True θ | Coverage | Nominal | Mean Width |
|----------:|-------:|---------:|--------:|-----------:|
|         1 |   1.00 |   0.9495 |    0.95 |     0.1721 |
|         2 |   1.10 |   0.9343 |    0.95 |     0.1770 |
|         3 |   0.95 |   0.9596 |    0.95 |     0.1694 |
|         4 |   1.15 |   0.9545 |    0.95 |     0.1794 |
|         5 |   1.10 |   0.9444 |    0.95 |     0.1769 |

Coverage Probability of 95% Confidence Intervals

#### Sampling Distribution Visualization

``` r
# Plot sampling distributions
par(mfrow = c(1, min(m, 5)), mar = c(4, 4, 2, 1))
for (j in 1:min(m, 5)) {
    hist(est_valid[, j], breaks = 20, probability = TRUE,
         main = paste0("Component ", j),
         xlab = expression(hat(theta)[j]),
         col = "lightblue", border = "white")
    abline(v = theta[j], col = "red", lwd = 2, lty = 2)
    abline(v = mean(est_valid[, j]), col = "blue", lwd = 2)
    legend("topright", legend = c("True", "Mean Est."),
           col = c("red", "blue"), lty = c(2, 1), lwd = 2, cex = 0.7)
}
```

![](exponential_series_files/figure-html/sampling-dist-plot-1.png)

#### Summary

The Monte Carlo simulation demonstrates that the MLE for the exponential
series system with masked data:

1.  **Is approximately unbiased** - the bias is small relative to the
    true parameter values
2.  **Has well-characterized variance** - consistent with asymptotic
    theory
3.  **Achieves nominal coverage** - the $(1 - \alpha)$% confidence
    intervals contain the true value approximately $(1 - \alpha)$% of
    the time

These properties validate the use of this likelihood model for inference
about component failure rates from masked series system data.

## Sensitivity Analysis

The MLE properties depend on several factors: sample size, masking
probability, and right-censoring rate. In this section, we study how
these factors affect estimation accuracy.

### Effect of Masking Probability

The masking probability $p$ controls how much information about the
failed component is obscured. When $p = 0$, only the failed component is
in the candidate set (perfect information). As $p$ increases, more
non-failed components are included, making estimation more difficult.

``` r
set.seed(7231)

# Use smaller sample for sensitivity study
n_sens <- 500
B_sens <- 100

# Masking probabilities to test
p_values <- seq(0, 0.5, by = 0.1)

# Fixed right-censoring (25% censored)
q_sens <- 0.25
tau_sens <- rep(-(1/sum(theta))*log(q_sens), n_sens)

# Storage
mask_results <- list()
```

``` r
for (p_idx in seq_along(p_values)) {
    p_curr <- p_values[p_idx]
    est_p <- matrix(NA, nrow = B_sens, ncol = m)

    for (b in 1:B_sens) {
        # Generate data
        comp_b <- matrix(nrow = n_sens, ncol = m)
        for (j in 1:m) comp_b[, j] <- rexp(n_sens, theta[j])
        comp_b <- md_encode_matrix(comp_b, "t")

        data_b <- comp_b %>%
            md_series_lifetime_right_censoring(tau_sens) %>%
            md_bernoulli_cand_c1_c2_c3(p_curr) %>%
            md_cand_sampler()

        tryCatch({
            fit_b <- solver(data_b, par = theta0, method = "Nelder-Mead")
            if (fit_b$converged) est_p[b, ] <- fit_b$theta.hat
        }, error = function(e) NULL)
    }

    # Compute statistics
    valid_p <- !is.na(est_p[, 1])
    mask_results[[p_idx]] <- list(
        p = p_curr,
        bias = colMeans(est_p[valid_p, , drop = FALSE]) - theta,
        variance = apply(est_p[valid_p, , drop = FALSE], 2, var),
        mse = (colMeans(est_p[valid_p, , drop = FALSE]) - theta)^2 +
              apply(est_p[valid_p, , drop = FALSE], 2, var)
    )
}
```

``` r
# Extract bias and MSE for plotting
bias_mat <- sapply(mask_results, function(x) mean(abs(x$bias)))
mse_mat <- sapply(mask_results, function(x) mean(x$mse))

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Mean absolute bias vs masking probability
plot(p_values, bias_mat, type = "b", pch = 19, col = "blue",
     xlab = "Masking Probability (p)",
     ylab = "Mean Absolute Bias",
     main = "Bias vs Masking Probability")
grid()

# Mean MSE vs masking probability
plot(p_values, mse_mat, type = "b", pch = 19, col = "red",
     xlab = "Masking Probability (p)",
     ylab = "Mean MSE",
     main = "MSE vs Masking Probability")
grid()
```

![](exponential_series_files/figure-html/masking-sensitivity-plot-1.png)

``` r
mask_df <- data.frame(
    p = p_values,
    Mean_Abs_Bias = bias_mat,
    Mean_MSE = mse_mat,
    Mean_RMSE = sqrt(mse_mat)
)

knitr::kable(mask_df, digits = 4,
             caption = "Effect of Masking Probability on Estimation Accuracy",
             col.names = c("Masking Prob.", "Mean |Bias|", "Mean MSE", "Mean RMSE"))
```

| Masking Prob. | Mean \|Bias\| | Mean MSE | Mean RMSE |
|--------------:|--------------:|---------:|----------:|
|           0.0 |        0.0133 |   0.0143 |    0.1197 |
|           0.1 |        0.0136 |   0.0186 |    0.1365 |
|           0.2 |        0.0065 |   0.0212 |    0.1456 |
|           0.3 |        0.0194 |   0.0303 |    0.1741 |
|           0.4 |        0.0285 |   0.0389 |    0.1973 |
|           0.5 |        0.0240 |   0.0590 |    0.2430 |

Effect of Masking Probability on Estimation Accuracy

As expected, increasing the masking probability generally increases both
bias and MSE. With $p = 0$ (no masking of non-failed components), we
have the most information and achieve the best estimation accuracy. As
$p$ increases toward 0.5, the candidate sets become less informative,
leading to less precise estimates.

### Effect of Right-Censoring Rate

Right-censoring reduces the number of exact failure times observed. When
a system is right-censored, we know the system survived beyond the
censoring time, but we don’t observe the actual failure time or failed
component.

``` r
set.seed(7231)

# Censoring quantiles (proportion surviving past tau)
# q = 0.1 means 10% survive (heavy censoring)
# q = 0.9 means 90% survive (light censoring)
q_values <- c(0.9, 0.7, 0.5, 0.3, 0.1)  # Survival probabilities

# Fixed masking probability
p_cens <- 0.2

# Storage
cens_results <- list()
```

``` r
for (q_idx in seq_along(q_values)) {
    q_curr <- q_values[q_idx]
    tau_curr <- rep(-(1/sum(theta))*log(q_curr), n_sens)
    est_q <- matrix(NA, nrow = B_sens, ncol = m)

    for (b in 1:B_sens) {
        # Generate data
        comp_b <- matrix(nrow = n_sens, ncol = m)
        for (j in 1:m) comp_b[, j] <- rexp(n_sens, theta[j])
        comp_b <- md_encode_matrix(comp_b, "t")

        data_b <- comp_b %>%
            md_series_lifetime_right_censoring(tau_curr) %>%
            md_bernoulli_cand_c1_c2_c3(p_cens) %>%
            md_cand_sampler()

        tryCatch({
            fit_b <- solver(data_b, par = theta0, method = "Nelder-Mead")
            if (fit_b$converged) est_q[b, ] <- fit_b$theta.hat
        }, error = function(e) NULL)
    }

    # Compute statistics
    valid_q <- !is.na(est_q[, 1])
    cens_rate <- 1 - mean(data_b$delta)  # Actual censoring rate
    cens_results[[q_idx]] <- list(
        q = q_curr,
        cens_rate = cens_rate,
        bias = colMeans(est_q[valid_q, , drop = FALSE]) - theta,
        variance = apply(est_q[valid_q, , drop = FALSE], 2, var),
        mse = (colMeans(est_q[valid_q, , drop = FALSE]) - theta)^2 +
              apply(est_q[valid_q, , drop = FALSE], 2, var)
    )
}
```

``` r
# Extract statistics
cens_rates <- sapply(cens_results, function(x) x$cens_rate)
bias_cens <- sapply(cens_results, function(x) mean(abs(x$bias)))
mse_cens <- sapply(cens_results, function(x) mean(x$mse))

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Mean absolute bias vs censoring rate
plot(cens_rates * 100, bias_cens, type = "b", pch = 19, col = "blue",
     xlab = "Censoring Rate (%)",
     ylab = "Mean Absolute Bias",
     main = "Bias vs Censoring Rate")
grid()

# Mean MSE vs censoring rate
plot(cens_rates * 100, mse_cens, type = "b", pch = 19, col = "red",
     xlab = "Censoring Rate (%)",
     ylab = "Mean MSE",
     main = "MSE vs Censoring Rate")
grid()
```

![](exponential_series_files/figure-html/censoring-sensitivity-plot-1.png)

``` r
cens_df <- data.frame(
    Cens_Rate_Pct = round(cens_rates * 100, 1),
    Mean_Abs_Bias = bias_cens,
    Mean_MSE = mse_cens,
    Mean_RMSE = sqrt(mse_cens)
)

knitr::kable(cens_df, digits = 4,
             caption = "Effect of Right-Censoring Rate on Estimation Accuracy",
             col.names = c("Censoring %", "Mean |Bias|", "Mean MSE", "Mean RMSE"))
```

| Censoring % | Mean \|Bias\| | Mean MSE | Mean RMSE |
|------------:|--------------:|---------:|----------:|
|        91.0 |        0.0556 |   0.1745 |    0.4177 |
|        68.2 |        0.0173 |   0.0568 |    0.2383 |
|        52.2 |        0.0118 |   0.0340 |    0.1843 |
|        28.8 |        0.0183 |   0.0245 |    0.1566 |
|         7.8 |        0.0085 |   0.0171 |    0.1308 |

Effect of Right-Censoring Rate on Estimation Accuracy

Higher censoring rates (more systems surviving past the observation
period) lead to increased bias and MSE. This is expected because:

1.  Censored observations contribute less information to the likelihood
2.  For censored systems, we only know the system survived beyond
    $\tau$, not which component would have failed
3.  No candidate set information is available for censored observations

### Practical Recommendations

Based on these sensitivity analyses:

1.  **Sample size**: Larger samples (n \> 500) generally provide stable,
    well-behaved estimates
2.  **Masking probability**: Keep $p$ as low as practically possible.
    Even moderate masking ($p \leq 0.3$) produces acceptable results
    with sufficient sample size
3.  **Censoring**: Heavy censoring (\> 50%) significantly degrades
    estimation quality. Design experiments with adequate follow-up time
    when possible
4.  **Combined effects**: When both masking and censoring are high,
    consider increasing sample size to compensate for information loss
