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
