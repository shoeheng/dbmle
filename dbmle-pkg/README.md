dbmle
================

Design-based maximum likelihood estimation (MLE) for identifying the joint distribution of always-takers, compliers, defiers, and never-takers from randomization assignment (Z) and takeup (D) data.

This package provides a design-based framework for estimating latent response type distributions in experiments with binary assignment and takeup. It supports both aggregated 2×2 count data and individual-level observations. Note that this design-based likelihood works only for experiments using a Bernoulli randomized design or a completely randomized design. 

----------------------------------------------------------------------
Installation
----------------------------------------------------------------------

### Option 1 - Install directly from Github (recommended)

Open a terminal (or Command Prompt on Windows), navigate to your proejct directory, and run:

    python3 -m pip install "git+https://github.com/shoeheng/dbmle@main#subdirectory=dbmle-pkg"

Or on some systems: 

    -m pip install "git+https://github.com/shoeheng/dbmle.git@main#egg=dbmle&subdirectory=dbmle-pkg"

This installs the package and the command-line tool `dbmle`.

### Option 2 - Install from source

Open a terminal (or Command Prompt on Windows) and run:

    git clone https://github.com/shoeheng/dbmle.git
    cd dbmle/dbmle-pkg
    python3 -m pip install .

Requirements:
- Python ≥ 3.9  
- `pip` installed (usually included with Python)

----------------------------------------------------------------------
Usage
----------------------------------------------------------------------

The MLE can be calculated from aggregated counts $(x_{I1},x_{I0},x_{C1},x_{C0})$, which are the counts of subjects who take up in intervention $x_{I1}$,
who do not take up in intervention $x_{I0}$,
who take up in control $x_{C1}$,
and who do not take up in control $x_{C0}$.

### 1. Python 

A sample script in Python using `dbmle` with aggregate counts looks as follows:

```python
from dbmle import dbmle

# Inputs correspond to the observed 2×2 table:
# (xI1, xI0, xC1, xC0) = (50, 11, 23, 31)
res = dbmle(50, 11, 23, 31, method="approx", auxiliary=True)
print(res.report())
```

It can also be calculated from the individual-level data using the `dbmle_from_ZD` command:
```python
from dbmle import dbmle_from_ZD

# Z: randomized assignment indicator (1 = intervention, 0 = control)
# D: observed take-up indicator (1 = treatment takeup, 0 = did not)
Z = [1, 1, 0, 1, 0, 1, 0, 0]
D = [1, 0, 1, 1, 0, 1, 0, 0]

res = dbmle_from_ZD(Z, D, method="exhaustive", auxiliary=True)
print(res.report())
```

Both commands return a DBMLEResult object whose .report() method prints a formatted summary of standard statistics, and MLEs. More statistics are reported if `auxiliary = True`.

### 2. Command Line Usage

Once installed, you can also use dbmle directly in the command line. The first example Python code above would equivalently be

    dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31 --method exhuastive --auxiliary 

To set `auxiliary=False`, simply remove `--auxiliary` from the command line.

### 3. Using `dbmle` in Stata

Stata now supports calling Python directly, meaning one need not leave Stata. After installing `dbmle`, one can run 

    python set exec "C:PATH\TO\PYTHON\python.exe"

to set your Python path in Stata. Once your path is set, run
```stata
python:
from dbmle import dbmle
res = dbmle(2, 1, 1, 2, method="exhaustive", auxiliary=True)
print(res.report())
end
```

Note that if Stata says "restart required", run:

    exit, clear

and reopen Stata.

Alternatively, you can run shell command inside Stata using `!`. After `dbmle` is installed, you can run

    ! dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31 --method exhuastive --auxiliary 

to get the same result.

If this doesn't work, find the directory where the dbmle.exe is kept and create a global. It is typically found under the "Scripts" folder in Python. For example,

    global dbmle "C:\PATH\TO\PYTHON\Scripts\dbmle.exe"

You can then run the above code.

----------------------------------------------------------------------
Options
----------------------------------------------------------------------

Each command supports the following parameters:

- **method**  
  `"approx"` (fast approximation) or `"exhaustive"` (exhaustive grid search).  
  *Default:* `"approx"`

- **auxiliary**  
  Whether to include auxiliary statistics aside from just the MLE (and creidble sets if `method = "exhaustive"`)  
  *Default:* `False`
  **Important:** When `auxiliary=True`, the package automatically performs an **exhaustive** grid search to compute the exact MLE and credible sets, even if `method="approx"` was requested since the likelihood of every joint distribution needs to be calculated for the credible set.  
  
- **level**  
  Credible-set level (e.g. `0.95` for 95%).  
  *Default:* `0.95`

- **show_progress**  
  Whether to display a tqdm progress bar (mainly relevant for exhaustive mode).  
  *Default:* `True`
  
### Additional options for `dbmle_from_ZD(...)`

When using the individual-data level command (`Z` assignment, `D` take-up), you also have controls for invalid data handling:

- **`invalid_policy`**  
  How to handle entries in `Z` or `D` that are not clean binary values (accepted binaries are: `0/1`, `True/False`, `0.0/1.0`, and strings `"0"`/`"1"`).  
  Choices:  
  - `"drop"` – ignore any record with an invalid `Z` or `D` (default).  
  - `"coerce-0"` – coerce invalid values to `0`.  
  - `"coerce-1"` – coerce invalid values to `1`.  
  - `"raise"` – strict mode; raise a `ValueError` on the first invalid entry.  
  *Default:* `"drop"`

- **`warn_on_invalid`**  
  Emit a single summary warning if any entries were dropped or coerced. The warning reports counts and a few example indices.  
  *Default:* `True`

**Notes:**  
- After cleaning (according to `invalid_policy`), both arms must be non-empty (at least one `Z=1` and one `Z=0`), otherwise a `ValueError` is raised.  
- When `auxiliary=True` in `dbmle_from_ZD(...)`, the same override applies: an exhaustive grid search is run to obtain the exact MLE and credible sets, regardless of `method`.
----------------------------------------------------------------------
Note on Approximation
----------------------------------------------------------------------

When `method="approx"`, the package estimates the maximum likelihood estimates (MLEs) using a fast local search rather than testing every possible joint distribution of always-takers, compliers, defiers, and never-takers. For all experiment results resulting from experiments with an equal number of individuals in intervention and control up to a sample size of 200, the appoximation is correct. We also randomly sampled 100 different experiment results for experiments with sample sizes between 500 and 1000, intervention and control not necessarily equally sizes, and find the approximation is correct for all sampled experiment results.

### Initialization

The approximation algorithm begins by considering three candidate joint distributions. The first is finding the distribution on the "corner" of the estimated Fréchet set with the highest likelihood.  
It also considers two simple joint distributions:

- Only always-takers and never-takers: $(\theta_{11},\theta_{10},\theta_{01},\theta_{00})=(x_{I1}+x_{C1},0,0,x_{I0}+x_{C0})$
- Only compliers and defiers: $(\theta_{11},\theta_{10},\theta_{01},\theta_{00})=(0, x_{I1}+x_{C0},x_{I0}+x_{C1},0)$

These provide plausible starting points for the likelihood search. If either the joint dsitribution with only always and never takers or the joint distribution with only compliers and defiers has the highest likelihood, it is chosen as the MLE.
Otherwise, the algorithm moves to the local cube search.

### Local Cube Search

Around the "corner" of the estimated Fréchet set with the highest likelihood, the algorithm searches a small four-dimensional integer cube centered on it, checking all nearby joint distributions that:

- sum exactly to *n*, and  
- have nonnegative counts.

Empirically, we find that the maximum $\ell_\inf$ distance between the the true MLE and the highest likelihood Fréchet corner increases at a rate just around $O(\sqrt{n})$ (barring MLEs that consist only of always takers and never takers or compliers and defiers). Accordingly, the cube’s half-width (`delta`) scales at the rate $0.3n^{0.58}$, a scale emprically decided by fitting a power law relating the maximum $\ell_\inf$ distance to the sample size, ensuring the search is both fast and likeliy covers the true MLE.  
Each candidate’s log-likelihood is evaluated, and the algorithm keeps all points within a small numerical tolerance of the best value.

### Edge Expansion

If the best solution lies on the cube’s boundary, the algorithm performs one larger pass with an expanded cube.  
This step captures nearby high-likelihood points the first pass might miss without having to exhaustively enumerate every possible combination.

----------------------------------------------------------------------
Example Output
----------------------------------------------------------------------

Below is a full example of the output printed by `.report()`, resulting from the first example Python script:

```
UserWarning: Note: When 'auxiliary=True' is specified with method='approx', the function automatically performs an exhaustive grid search to compute the true MLE and 95% credible set. The MLE is exact as well.
  res = dbmle(
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

Homepage: https://github.com/shoeheng/dbmle  
Issues: https://github.com/shoeheng/dbmle/issues

----------------------------------------------------------------------
Citation
----------------------------------------------------------------------

If you use this package in academic work, please cite it as:

citation.



















