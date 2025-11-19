`dbmle`
================

This package provides a design-based maximum likelihood estimate of the numbers of always takers, compliers, defiers, and never takers in the sample of people in an experiment using the method from Christy and Kowalski (2025). It supports both individual-level data and aggregated counts. Note that this design-based likelihood works only for experiments using a Bernoulli randomized design or a completely randomized design. We also provide directions to use this package in Stata. 

----------------------------------------------------------------------
Installation
----------------------------------------------------------------------

Open a terminal or command prompt, navigate to your project directory, and run:

    python3 -m pip install "git+https://github.com/shoeheng/dbmle@main#subdirectory=dbmle-pkg"

Once the package is published on PyPI, the installation command to install the package will be

    pip install dbmle

This installs the package and the command-line tool `dbmle`. 

Requirements:
- Python ≥ 3.9  
- `pip` installed (usually included with Python)

----------------------------------------------------------------------
Usage
----------------------------------------------------------------------

The MLE can be calculated from aggregated counts $(x_{I1},x_{I0},x_{C1},x_{C0})$, which are the counts of subjects who take up in intervention $x_{I1}$,
do not take up in intervention $x_{I0}$,
take up in control $x_{C1}$,
and do not take up in control $x_{C0}$.

### 1. Within Python 

A sample script in Python using `dbmle` with aggregate counts looks as follows:

```python
from dbmle import dbmle

# Inputs correspond to aggregated count data:
# (xI1, xI0, xC1, xC0) = (50, 11, 23, 31), from Johnson and Goldstein (2003)
res = dbmle(50, 11, 23, 31)
print(res.report())
```

It can also be calculated from the individual-level data using the `dbmle_from_ZD` command:
```python
from dbmle import dbmle_from_ZD

# Z: randomized assignment indicator (1 = intervention, 0 = control)
# D: observed takeup indicator (1 = treatment takeup, 0 = did not)
Z = [1, 1, 1, 0, 0, 0]
D = [1, 1, 0, 1, 0, 0]

res = dbmle_from_ZD(Z, D)
print(res.report())
```

Both commands return a DBMLEResult object whose .report() method prints a formatted summary of standard statistics, and MLEs.

### 2. Command Line Usage

Once installed, you can also use dbmle directly in the command line. The first example Python code above would equivalently be

    dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31 

----------------------------------------------------------------------
Using `dbmle` in Stata
----------------------------------------------------------------------

Below is a guide to using Python within Stata. The first step is installing Python. If you already have Python 3.9 or higher installed, skip this step.

### 1. Install Python

Go to https://www.python.org/downloads/ and download the latest Python installer for your operating system. run the installer and check the option "Add Python to PATH" after running the installer.

### 2. Open Stata and tell Stata how to use Python

In the Stata command prompt, run

    python query

You should get an output similar to

    Python Settings
      set python_exec      /path/to/python
      set python_userpath  

    Python system information
      initialized          no
      version              3.12.8
      architecture         64-bit
      library path         /.../lib64/libpython3.12.so.1.0

that is, `set python_exec` should have a valid path to python and `version` should be 3.9 or greater. Once this is done, `dbmle` can be used directly in the command line of stata by typing a command like 

    ! dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31

----------------------------------------------------------------------
Parameters
----------------------------------------------------------------------

Aisde from the data input, each command supports the following parameters:

- **output:** `"basic"`, `"auxiliary"`, _or_ `"approx"`  
  *Default:* `"approx"`  
What statistics are to be calculated and displayed. `"basic"` performs an exhaustive grid search and returns the MLE(s) along with the smallest credible set. `"auxiliary`" returns the statistics that `"basic"` returns along with the largest possible support, estimated Frechet bounds, and the smallest credible set conditional on being within the estimated Frechet set. `"approx"` uses a significantly faster approximation algorithm to calculate the MLE(s) and only returns the MLE(s). All three return a standard statistics table as well.
  
- **level:** _float_  
  *Default:* `0.95`  
  Smallest credible-set level (e.g. `0.95` for 95%).

- **show_progress:** _bool_   
  *Default:* `True`  
  Whether to display a `tqdm` progress bar for the exhaustive grid search (not relevant `output="approx"`). 

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
































