`dbmle`
================

This package provides a design-based maximum likelihood estimate of the numbers of always takers, compliers, defiers, and never takers in the sample of people in an experiment, following [Christy and Kowalski (2025)](https://arxiv.org/abs/2412.16352) "Counting Defiers: A Design-Based Model of an Experiment Can Reveal Evidence Beyond the Average Effect". The package is compatible with experiments that use either a Bernoulli randomized design or a completely randomized design, and results are exact in both cases. This package also reports other statistics from Christy and Kowalski (2025), including credible sets. See below for instructions on how to use this package in Stata.   

----------------------------------------------------------------------
Installation
----------------------------------------------------------------------

Open a terminal or command prompt, navigate to your project directory, and run:

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
# D: observed takeup indicator (1 = takeup, 0 = no takeup)
Z = [1, 1, 1, 0, 0, 0]
D = [1, 1, 0, 1, 0, 0]
# Data is equivalent to stylized example in Figure 1 of Christy and Kowalski (2025)

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

Go to https://www.python.org/downloads/ and download the latest Python installer for your operating system. Run the installer and check the option "Add Python to PATH" after running the installer.

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

that is, `set python_exec` should have a valid path to Python and `version` should be 3.9 or greater. Once this is done, `dbmle` can be used directly in the command line of Stata by typing a command like 

    ! dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31

You can also use Python in a Stata script by starting the script with `python:` and ending with `end`:

```
python:
    from dbmle import dbmle
    res = dbmle(50, 11, 23, 31)
    print(res.report())
end
```

You can find more information about using Python in Stata here: https://www.stata.com/python/

----------------------------------------------------------------------
Parameters
----------------------------------------------------------------------

Aside from the data input, each command supports the following parameters:

- **output:** `"basic"`, `"auxiliary"`, _or_ `"approx"`  
  *Default:* `"basic"`  
What statistics are to be calculated and displayed. `"basic"` performs an exhaustive grid search and returns the MLE(s) along with the smallest credible set. `"auxiliary`" returns the statistics that `"basic"` returns along with the largest possible support, estimated Fréchet bounds, and the smallest credible set conditional on being within the estimated Fréchet set. `"approx"` uses a significantly faster approximation algorithm to calculate the MLE(s) and only returns the MLE(s). All three return a standard statistics table as well.
  
- **level:** _float_  
  *Default:* `0.95`  
  Smallest credible-set level (e.g. `0.95` for 95%).

- **show_progress:** _bool_   
  *Default:* `True`  
  Whether to display a `tqdm` progress bar for the exhaustive grid search (not relevant when `output="approx"`). For the Tappin et al. (2015) data, a sample of 612 where $(x_{I1},x_{I0},x_{C1},x_{C0})=(69, 237, 26, 280)$, the `tqdm` bar looks as follows:

        Enumerating Joint Distributions:   7%|▋         | 2814961/38579155 [00:36<09:53, 60237.47Joint Distribution/s]
  
  The `7%|▋         |` is a visual indicator of how many distributions in the grid have been evaluated so far. The fraction `2814961/38579155` tells you exactly how many distributions have been calculated over how many distributions total there are to calculate. The times `00:36<09:53` tell you the time elapsed so far and the predicted time left to complete the grid search, so in this case, the total run time is expected to be about 10.5 minutes. Finally, `60237.47Joint Distribution/s` tells you how many distributions can be calculated in a second (the predicted time left is based on the history of this value).

Note that there is no parameter for the design of the randomization (Bernoulli randomized or completely randomized) since the results are the same, as their respective design-based likelihoods are proportional to each other.

----------------------------------------------------------------------
Note on Approximation
----------------------------------------------------------------------

When `method="approx"`, the package estimates MLEs using a fast local search rather than testing every possible distribution of always takers, compliers, defiers, and never takers. For all outcomes resulting from experiments with an equal number of individuals in intervention and control up to a sample size of 200, we have verified the approximation is correct.

The approximation begins by considering three candidate joint distributions:

- The endpoint of the estimated Fréchet set with the highest likelihood.  
- A joint distribution with only always takers and never takers, so $x_{I1}+x_{C1}$ always takers and $x_{I0}+x_{C0}$ never takers.
- A joint distribution with only compliers and defiers, so $x_{I1}+x_{C0}$ compliers and $x_{I0}+x_{C1}$ defiers.

If either of the two-type joint distributions attains the highest likelihood, it is immediately returned as the MLE. Otherwise, the algorithm performs a local cube search. In this case, we search a small four-dimensional integer cube around the endpoint of the estimated Fréchet set with the highest likelihood. The width of the search increases with sample size. If the distribution with the highest likelihood lies on the boundary of the cube, we search a slightly larger cube. 

----------------------------------------------------------------------
Example Usages
----------------------------------------------------------------------

To get the basic output, in Python, you can run 

```
from dbmle import dbmle

# (xI1, xI0, xC1, xC0) = (50, 11, 23, 31)
res = dbmle(50, 11, 23, 31)
print(res.report())
```
or equivalently in Stata,

```
python:
    from dbmle import dbmle
    res = dbmle(50, 11, 23, 31)
    print(res.report())
end
```

and in the command line interface,

```
dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31
```

All will give the same output:

```
Standard Statistics
----------------------------------------------
Average Effect              50/61 - 23/54 = 39.37%
95% Confidence Interval     [23.03%, 55.72%]
Fisher's Exact Test p-value 1.552e-05
Intervention Takeup Rate    50/61 = 81.97%
Control Takeup Rate         23/54 = 42.59%
Sample Size                  115


Christy and Kowalski Design-Based Maximum Likelihood Estimates
------------------------------------------------------------------------------------
Always takers
  MLE: 28/115 = 24.35%
  95% Smallest Credible Set: [0,63]/115 = [0.00%, 54.78%]

Compliers
  MLE: 66/115 = 57.39%
  95% Smallest Credible Set: [23,81]/115 = [20.00%, 70.43%]

Defiers
  MLE: 21/115 = 18.26%
  95% Smallest Credible Set: [0,34]/115 = [0.00%, 29.57%]

Never takers
  MLE: 0/115 = 0.00%
  95% Smallest Credible Set: [0,32]/115 = [0.00%, 27.83%]
```

To get the auxiliary table output, you can run

```
from dbmle import dbmle

res = dbmle(50, 11, 23, 31, output="auxiliary")
print(res.report())
```

or in Stata,

```
python:
    from dbmle import dbmle
    res = dbmle(50, 11, 23, 31, output="auxiliary")
    print(res.report())
end
```

or in the command line interface

```
dbmle --xI1 50 --xI0 11 --xC1 23 --xC0 31 --output auxiliary
```

which will all yield the same output:

```
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

Replace `"auxiliary"` with `"approx"` to use the fast MLE approximation, which will return the output

```
Standard Statistics
----------------------------------------------
Average Effect              50/61 - 23/54 = 39.37%
95% Confidence Interval     [23.03%, 55.72%]
Fisher's Exact Test p-value 1.552e-05
Intervention Takeup Rate    50/61 = 81.97%
Control Takeup Rate         23/54 = 42.59%
Sample Size                  115


Christy and Kowalski Design-Based Maximum Likelihood Estimates*
------------------------------------------------------------------------------------
Always takers
  MLE: 28/115 = 24.35%

Compliers
  MLE: 66/115 = 57.39%

Defiers
  MLE: 21/115 = 18.26%

Never takers
  MLE: 0/115 = 0.00%

* MLE estimates obtained from an approximation algorithm implemented 
by the dbmle package (Christy, Kowalski, and Zhang 2025)
```

----------------------------------------------------------------------
Project Links
----------------------------------------------------------------------

- PyPI: https://pypi.org/project/dbmle/
- Source code: https://github.com/shoeheng/dbmle
- Issue tracker: https://github.com/shoeheng/dbmle/issues

----------------------------------------------------------------------
Citation
----------------------------------------------------------------------

If you use `dbmle` in your academic work, please cite Christy and Kowalski (2025) along with this package:

```
@misc{christy2025dbmle,
  author       = {Christy, Neil and Kowalski, Amanda and Zhang, Shuheng},
  title        = {dbmle: Design-Based Maximum Likelihood Estimation for Always Takers, Compliers, Defiers, and Never Takers},
  year         = {2025},
  howpublished = {\url{https://pypi.org/project/dbmle/}},
  note         = {Python package version 0.0.2}
}
```

















