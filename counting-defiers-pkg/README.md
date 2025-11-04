counting-defiers
================

Design-based maximum likelihood estimation (MLE) for identifying the joint distribution of always-takers, compliers, defiers, and never-takers from randomization assignment (Z) and takeup (D) data.

This package provides a design-based framework for estimating latent response type distributions in binary treatment takeup experiments. It supports both aggregated 2×2 count data and individual-level observations.

----------------------------------------------------------------------
Installation
----------------------------------------------------------------------

Install from source:

    pip install .

Once the package is published to PyPI, users will be able to install it directly with:

    pip install counting-defiers

Requirements:
- Python ≥ 3.9  
- tqdm ≥ 4.0 (used for progress bars and time estimates during exhaustive grid search)

----------------------------------------------------------------------
Quick Start
----------------------------------------------------------------------

From aggregated counts:

```python
from counting_defiers import counting_defiers_command

# Inputs correspond to the observed 2×2 table:
# (xI1, xI0, xC1, xC0)
res = counting_defiers_command(50, 11, 23, 31, method="approx", auxiliary=False)
print(res.report())
```

From individual-level data:
```python
from counting_defiers import counting_defiers_from_ZD

# Z: randomized assignment indicator (1 = intervention, 0 = control)
# D: observed take-up indicator (1 = treatment takeup, 0 = did not)
Z = [1, 1, 0, 1, 0, 1, 0, 0]
D = [1, 0, 1, 1, 0, 1, 0, 0]

res = counting_defiers_from_ZD(Z, D, method="exhaustive", auxiliary=True)
print(res.report())
```

Both commands return a CountingDefiersResult object whose .report() method prints a formatted summary of standard statistics, and MLEs. More statistics are reported if `auxiliary = True`.

----------------------------------------------------------------------
Options
----------------------------------------------------------------------

Each command supports the following parameters:

- **method**  
  `"approx"` (fast approximation) or `"exhaustive"` (exhaustive grid search).  
  *Default:* `"approx"`

- **auxiliary**  
  Whether to include auxiliary statistics aside from just the MLE (and creidble sets if `method = "exhaustive"`)  
  *Default:* `"proposed"`

- **level**  
  Credible-set level (e.g. `0.95` for 95%).  
  *Default:* `0.95`

- **show_progress**  
  Whether to display a tqdm progress bar (mainly relevant for exhaustive mode).  
  *Default:* `True`

----------------------------------------------------------------------
Note on Approximation
----------------------------------------------------------------------

When `method="approx"`, the package estimates the maximum likelihood estimates (MLEs) using a fast local search rather than testing every possible joint distribution of always-takers, compliers, defiers, and never-takers.

### Initialization

The approximation algorithm begins by finding the distribution on the "corner" of the estimated Fréchet set with the highest likelihood.  
It also considers two simple joint distributions:

- Only always-takers and never-takers  
- Only compliers and defiers  

These provide plausible starting points for the likelihood search. If either the joint dsitribution with only always and never takers or the joint distribution with only compliers and defiers has the highest likelihood, it is chosen as the MLE.
Otherwise, the algorithm moves to the local cube search.

### Local Cube Search

Around the "corner" of the estimated Fréchet set with the highest likelihood, the algorithm searches a small four-dimensional integer cube centered on it, checking all nearby joint distributions that:

- sum exactly to *n*, and  
- have nonnegative counts.

Empirically, we find that the maximum $\ell_\inf$ distance between the the true MLE and the highest likelihood Fréchet corner increases at a rate just under $O(\sqrt{n})$ (barring MLEs that consist only of always takers and never takers or compliers and defiers). Accordingly, the cube’s half-width (`delta`) scales at the rate $\sqrt{n}$, ensuring the search is both fast and likeliy covers the true MLE.  
Each candidate’s log-likelihood is evaluated, and the algorithm keeps all points within a small numerical tolerance of the best value.

### Edge Expansion

If the best solution lies on the cube’s boundary, the algorithm performs one larger pass with an expanded cube.  
This step captures nearby high-likelihood points the first pass might miss without having to exhaustively enumerate every possible combination.

### Practical Implications

- The approximate method usually matches the true global MLE for moderate sample sizes while running much faster.  
- It returns all tied MLEs found within tolerance.  
- If you need guaranteed global maxima and full credible sets, use `method="exhaustive"` instead.

----------------------------------------------------------------------
Example Output
----------------------------------------------------------------------

Below is a full example of the output printed by `.report()`:

```
Enumerating Joint Distributions: 100%|██████████| 266916/266916 [00:01<00:00, 152926.72Joint Distribution/s]

Standard Statistics
----------------------------------------------
Average Effect              50/61 - 23/54 = 39.37%
95% Confidence Interval     [23.03%, 55.72%]
Fisher's Exact Test p-value 1.552e-05
Intervention Takeup Rate    50/61 = 81.97%
Control Takeup Rate         23/54 = 42.59%
Sample Size                  115


Christy and Kowalski Design-Based Maximum Likelihood Estimates and Auxiliary Statistics
------------------------------------------------------------------------------------
Always takers
  MLE: 28/115 = 24.35%
  95% Smallest Credible Set: [0,63]/115 = [0.00%, 54.78%]
  Largest Possible Support: [0,73]/115 = [0.00%, 63.48%]
  Estimated Frechet Bounds: [28,49]/115 = [24.35%, 42.61%]
  95% SCS within Est. Frechet: [28,39]/115 U [41,49]/115 = [24.35%, 33.91%] U [35.65%, 42.61%]

Compliers
  MLE: 66/115 = 57.39%
  95% Smallest Credible Set: [23,81]/115 = [20.00%, 70.43%]
  Largest Possible Support: [0,81]/115 = [0.00%, 70.43%]
  Estimated Frechet Bounds: [45,66]/115 = [39.13%, 57.39%]
  95% SCS within Est. Frechet: [45,53]/115 U [55,66]/115 = [39.13%, 46.09%] U [47.83%, 57.39%]

Defiers
  MLE: 21/115 = 18.26%
  95% Smallest Credible Set: [0,34]/115 = [0.00%, 29.57%]
  Largest Possible Support: [0,34]/115 = [0.00%, 29.57%]
  Estimated Frechet Bounds: [0,21]/115 = [0.00%, 18.26%]
  95% SCS within Est. Frechet: [0,8]/115 U [10,21]/115 = [0.00%, 6.96%] U [8.70%, 18.26%]

Never takers
  MLE: 0/115 = 0.00%
  95% Smallest Credible Set: [0,32]/115 = [0.00%, 27.83%]
  Largest Possible Support: [0,42]/115 = [0.00%, 36.52%]
  Estimated Frechet Bounds: [0,21]/115 = [0.00%, 18.26%]
  95% SCS within Est. Frechet: [0,11]/115 U [13,21]/115 = [0.00%, 9.57%] U [11.30%, 18.26%]

```

----------------------------------------------------------------------
License
----------------------------------------------------------------------

License information here.

----------------------------------------------------------------------
Project Links
----------------------------------------------------------------------

Homepage: https://github.com/YOURNAME/counting-defiers  
Issues: https://github.com/YOURNAME/counting-defiers/issues

----------------------------------------------------------------------
Citation
----------------------------------------------------------------------

If you use this package in academic work, please cite it as:

citation.



