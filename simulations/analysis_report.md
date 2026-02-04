# Simulation Analysis: MLE for Exponential Series Systems with Masked Data

## Model and Notation

We study maximum likelihood estimation of component failure rates
$\theta = (\lambda_1, \ldots, \lambda_m)$ for a series system with
exponentially distributed component lifetimes, under two forms of masking:

1. **Right-censoring**: system lifetime $T_i = \min(\tau, T_{i1}, \ldots, T_{im})$
2. **Candidate set masking**: the failed component is identified only up to
   a candidate set $\mathcal{C}_i \subseteq \{1, \ldots, m\}$

The candidate set model satisfies conditions C1 (failed component in set),
C2 (uniform probability given component cause), and C3 (independence from
$\theta$), yielding a reduced likelihood that depends only on observed data.

### True Parameters

| Component | $\lambda_j$ |
|-----------|-------------|
| 1         | 1.00        |
| 2         | 1.10        |
| 3         | 0.95        |
| 4         | 1.15        |
| 5         | 1.10        |

System failure rate: $\sum \lambda_j = 5.30$.

---

## Experiment 1: Monte Carlo Study

**Design:** $n = 7500$, $B = 200$ replications, masking probability $p = 0.3$,
survival probability $q = 0.25$ (approximately 25% censoring).
Optimization: Nelder-Mead. 95% Wald confidence intervals.

**Convergence rate:** 99% (198/200).

### Bias, Variance, and MSE

| Component | True $\lambda$ | Mean $\hat\lambda$ | Bias | Variance | RMSE | Rel. Bias % |
|-----------|----------------|---------------------|----------|----------|--------|-------------|
| 1         | 1.00           | 1.000               | 0.0005   | 0.0019   | 0.0439 | 0.05        |
| 2         | 1.10           | 1.100               | 0.0004   | 0.0018   | 0.0426 | 0.04        |
| 3         | 0.95           | 0.949               | -0.0010  | 0.0021   | 0.0463 | -0.11       |
| 4         | 1.15           | 1.151               | 0.0011   | 0.0021   | 0.0453 | 0.10        |
| 5         | 1.10           | 1.107               | 0.0073   | 0.0020   | 0.0452 | 0.67        |

### Confidence Interval Coverage

| Component | True $\lambda$ | Coverage | Nominal | Mean CI Width |
|-----------|----------------|----------|---------|---------------|
| 1         | 1.00           | 0.944    | 0.95    | 0.1725        |
| 2         | 1.10           | 0.965    | 0.95    | 0.1775        |
| 3         | 0.95           | 0.955    | 0.95    | 0.1697        |
| 4         | 1.15           | 0.934    | 0.95    | 0.1797        |
| 5         | 1.10           | 0.950    | 0.95    | 0.1776        |

### Key Findings

- **Essentially unbiased.** All relative biases are below 0.7%. The
  squared-bias contribution to MSE is negligible (e.g., component 1:
  bias^2 = 2e-7 vs variance = 0.0019).

- **Uniform relative precision.** RMSE ranges 0.043--0.046, or roughly
  4--5% of the true rates. Components with larger failure rates have
  proportionally wider CIs, consistent with asymptotic efficiency.

- **Coverage near nominal.** Coverage ranges 93.4%--96.5%. Component 4
  at 93.4% slightly undercovers, a known finite-sample property of Wald
  intervals for rate parameters.

- **CI width as design metric.** Mean widths 0.170--0.180 correspond
  to approximately +/-8.5% of the true rate. Width scales as $1/\sqrt{n}$.

---

## Experiment 2: Effect of Masking Probability

**Design:** $n = 500$, $B = 100$ per $p$ value, $p \in \{0, 0.1, \ldots, 0.5\}$,
fixed censoring at ~25%.

| Masking $p$ | Mean |Bias| | Mean MSE | Mean RMSE |
|-------------|-------------|----------|-----------|
| 0.0         | 0.0092      | 0.0148   | 0.1218    |
| 0.1         | 0.0120      | 0.0187   | 0.1368    |
| 0.2         | 0.0082      | 0.0211   | 0.1451    |
| 0.3         | 0.0152      | 0.0325   | 0.1802    |
| 0.4         | 0.0142      | 0.0396   | 0.1991    |
| 0.5         | 0.0260      | 0.0623   | 0.2495    |

### Key Findings

- MSE increases 4.2x from $p=0$ (0.015) to $p=0.5$ (0.062). RMSE roughly
  doubles from 0.12 to 0.25.

- The degradation is driven entirely by variance, not bias. Mean |bias|
  fluctuates between 0.008 and 0.026 with no systematic trend.

- The marginal cost accelerates: MSE roughly doubles from $p=0$ to $p=0.3$,
  then doubles again from $p=0.3$ to $p=0.5$.

- Residual MSE of 0.015 at $p=0$ reflects the censoring-only contribution.

---

## Experiment 3: Effect of Right-Censoring Rate

**Design:** $n = 500$, $B = 100$ per $q$ value, $q \in \{0.1, 0.2, \ldots, 0.9\}$,
fixed masking at $p = 0.2$.

| Surv. $q$ | Cens. % | Mean |Bias| | Mean MSE | Mean RMSE |
|-----------|---------|-------------|----------|-----------|
| 0.1       | 10.0    | 0.0049      | 0.0193   | 0.1391    |
| 0.2       | 20.0    | 0.0126      | 0.0244   | 0.1561    |
| 0.3       | 30.1    | 0.0135      | 0.0279   | 0.1672    |
| 0.4       | 40.2    | 0.0113      | 0.0333   | 0.1824    |
| 0.5       | 50.2    | 0.0143      | 0.0402   | 0.2005    |
| 0.6       | 60.2    | 0.0185      | 0.0497   | 0.2230    |
| 0.7       | 69.9    | 0.0310      | 0.0590   | 0.2429    |
| 0.8       | 80.2    | 0.0138      | 0.0876   | 0.2959    |
| 0.9       | 90.1    | 0.0341      | 0.1829   | 0.4276    |

### Key Findings

- MSE grows 9.5x from 10% censoring (0.019) to 90% censoring (0.183), far
  exceeding the 4.2x range from masking. Censoring eliminates *both* the
  failure time and the candidate set.

- There is a clear inflection around 50% censoring. MSE is manageable at
  50% (0.040) but jumps to 0.059 at 70% and 0.183 at 90%. The relationship
  is convex.

- Bias remains modest even under extreme censoring (~3% of true rates at
  worst). The MLE is consistent but increasingly inefficient.

- At 90% censoring (~50 of 500 exact failures), the MLE still converges but
  with RMSE ~ 0.43 (~40% of true rates).

---

## Experiment 4: Effect of Sample Size

**Design:** $n \in \{100, 250, 500, 1000, 2500, 5000\}$, $B = 100$,
$p = 0.3$, ~25% censoring.

| $n$ | Conv. Rate | Mean RMSE | Mean Coverage | Mean CI Width |
|------|-----------|-----------|---------------|---------------|
| 100  | 0.97      | 0.3989    | 0.934         | 1.5267        |
| 250  | 0.98      | 0.2576    | 0.916         | 0.9540        |
| 500  | 0.98      | 0.1670    | 0.949         | 0.6788        |
| 1000 | 0.99      | 0.1169    | 0.970         | 0.4800        |
| 2500 | 1.00      | 0.0781    | 0.948         | 0.3042        |
| 5000 | 1.00      | 0.0533    | 0.956         | 0.2150        |

### CI Width Scaling Verification

The $1/\sqrt{n}$ scaling of CI width is confirmed empirically (reference: $n = 500$):

| $n$ | Width Ratio | Theory $\sqrt{500/n}$ | Ratio of Ratios |
|------|------------|------------------------|-----------------|
| 100  | 2.249      | 2.236                  | 1.006           |
| 250  | 1.405      | 1.414                  | 0.994           |
| 500  | 1.000      | 1.000                  | 1.000           |
| 1000 | 0.707      | 0.707                  | 1.000           |
| 2500 | 0.448      | 0.447                  | 1.002           |
| 5000 | 0.317      | 0.316                  | 1.002           |

All ratios are within 0.6% of the theoretical prediction.

### Key Findings

- CI width follows $1/\sqrt{n}$ scaling almost exactly, confirming the
  asymptotic normality of the MLE.

- Coverage is erratic at small $n$: 93.4% at $n=100$, drops to 91.6% at
  $n=250$, then stabilizes near 95% for $n \geq 500$. The Wald interval
  is unreliable below $n \approx 500$.

- Convergence rate improves from 97% at $n=100$ to 100% at $n \geq 2500$.

- At $n = 500$, CI width is ~0.68 (about 65% of the true rate), which is
  too wide for most practical purposes. At $n = 5000$, width shrinks to
  ~0.22 (about 20% of the true rate).

---

## Experiment 5: Joint Masking x Censoring Interaction

**Design:** $n = 500$, $B = 100$, $p \in \{0, 0.2, 0.4\}$,
$q \in \{0.1, 0.3, 0.5, 0.7, 0.9\}$.

### Mean MSE (across components)

|          | q=0.1  | q=0.3  | q=0.5  | q=0.7  | q=0.9  |
|----------|--------|--------|--------|--------|--------|
| p=0.0    | 0.0128 | 0.0170 | 0.0239 | 0.0425 | 0.1235 |
| p=0.2    | 0.0183 | 0.0223 | 0.0377 | 0.0536 | 0.1610 |
| p=0.4    | 0.0308 | 0.0431 | 0.0532 | 0.1041 | 0.2972 |

### Mean RMSE (across components)

|          | q=0.1  | q=0.3  | q=0.5  | q=0.7  | q=0.9  |
|----------|--------|--------|--------|--------|--------|
| p=0.0    | 0.113  | 0.130  | 0.155  | 0.206  | 0.351  |
| p=0.2    | 0.135  | 0.149  | 0.194  | 0.232  | 0.401  |
| p=0.4    | 0.175  | 0.207  | 0.231  | 0.323  | 0.545  |

### Key Findings

- The worst case is $(p = 0.4, q = 0.9)$ with MSE = 0.297, a $23\times$
  degradation from the best case $(p = 0, q = 0.1)$ with MSE = 0.013.

- Censoring dominates at all masking levels. Moving from $q = 0.1$ to
  $q = 0.9$ at $p = 0$ increases MSE by $9.6\times$. Moving from $p = 0$
  to $p = 0.4$ at $q = 0.1$ increases MSE by only $2.4\times$.

- The effects are roughly multiplicative (log-additive). The product of
  marginal effects ($9.6 \times 2.4 = 23.0$) closely matches the joint
  effect ($23.2\times$). This means the interaction between masking and
  censoring is weak --- their effects compound but do not amplify each other.

- At moderate settings $(p = 0.2, q = 0.3)$, RMSE = 0.149, or about
  14% of the true rates --- likely acceptable for most applications.

---

## Practical Recommendations

1. **Prioritize reducing censoring over masking.** Censoring degrades MSE by
   ~10x across its range versus ~4x for masking. Extending the observation
   window is more valuable than improving diagnostic resolution.

2. **Moderate masking is tolerable.** MSE only doubles from $p = 0$ to
   $p = 0.3$. If reducing masking below 0.3 requires expensive diagnostics,
   the incremental benefit may not justify the cost.

3. **The 50% censoring threshold.** Below ~50% censoring, MSE degrades
   gradually. Above 50%, MSE grows convexly. Experimental designs should
   aim for less than 50% censoring.

4. **Minimum sample size: $n \geq 500$.** Below this, Wald coverage is
   unreliable (91.6% at $n=250$ vs nominal 95%). At $n = 500$ with
   moderate masking and censoring, RMSE is ~17% of the true rates.

5. **CI widths for sample size planning.** Width scales as $1/\sqrt{n}$.
   At $n = 7500$ with $p = 0.3$ and 25% censoring, CI widths are ~0.17.
   Halving width to ~0.085 requires $n \approx 30000$.

6. **Joint effects are multiplicative.** When both masking and censoring
   are present, their effects on MSE roughly multiply. Use the marginal
   sensitivity curves to estimate the joint effect for a given $(p, q)$
   combination.

---

## Reproducibility

All results generated by the simulation framework in `simulations/`.
Base seed: 7231. Each replication uses a deterministic seed
derived from the base seed and replication index.

To regenerate:
```bash
cd simulations
Rscript run_all.R --reset
```
