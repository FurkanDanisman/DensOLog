# DensOLog

**DensOLog** is our density-estimation algorithm designed for **binned / rounded / interval-grouped** observations.  
This repository provides code to run DensOLog and compare it against **three baseline approaches** commonly used for binned-data density estimation.

---

## What’s in this repo

1. DensOLog and three baseline (KernSmooth, binnednp, NKDEBD) algorithms.
2. Simulation codes and results to compare methods.
3. Code to compare methods visually.
4. Application of DensOLog on real data and synthetic data. 

### Repository structure

- `Algorithms/`  
  Algorithms of **DensOLog** and the **3 comparison methods**, plus shared utilities (e.g., binning helpers, normalizations, evaluation metrics).

- `Simulations/`  
  Scripts to reproduce simulation studies: data generation, repeated runs, saving results, and summary tables.

- `Plot_Comparison/`  
  Plotting code used to generate comparison figures (overlayed estimated densities, error curves, etc.).

- `Application_Example/`  
  End-to-end example on a real (or illustrative) dataset: data prep → estimation → plots.
  
---

## Methods compared

This repo compares **DensOLog** with the following methods:

### 1) `KernSmooth` (binned kernel density estimation)
We use functions from the **KernSmooth** R package, which implements kernel smoothing and density estimation routines supporting the framework in *Wand & Jones (1995)*.

**Typical usage in this repo:** binned KDE via `bkde()` (and related helpers).

- Package: `KernSmooth` (Ripley; supporting Wand & Jones)
- Reference (book): Wand & Jones (1995)

---

### 2) `binnednp` (kernel estimation for interval-grouped data)
We use routines from the **binnednp** package, which provides kernel density / distribution estimation and bandwidth selection methods for **interval-grouped data**.

- Package: `binnednp`
- Reference: Barreiro-Ures et al. (2019)

> Note: `binnednp` has been archived/removed from CRAN in recent years; if installation is an issue, you can install from an archive snapshot or use the code paths in this repo that do not require it.

---

### 3) Nonlinear kernel density estimation for binned data (Bernoulli)
We include/replicate the method from:

- **Blower & Kelsall (2002)**, *Bernoulli*: “Nonlinear kernel density estimation for binned data: convergence in entropy.”

Bandwidth selection for this method uses selectors available in base R (`stats`), via the `bandwidth` documentation (e.g., `bw.nrd`, `bw.SJ`, etc.), which references classic selectors such as Silverman’s rule-of-thumb and Sheather–Jones.

---

## References

### KernSmooth

	- Wand, M. P., & Jones, M. C. (1995). Kernel Smoothing. Chapman & Hall.
	- Ripley, B. D. (maintainer). KernSmooth: Functions for Kernel Smoothing Supporting Wand & Jones (1995). R package.

###  binnednp
	- Barreiro-Ures, D., Francisco-Fernández, M., Cao, R., Fraguela, B. B., Doallo, R., González-Andújar, J. L., & Reyes, M. (2019). Analysis of interval-grouped data in weed science: The binnednp Rcpp package. (Preprint / and published version in Ecology & Evolution).

###  Nonlinear KDE for binned data
	- Blower, G., & Kelsall, J. E. (2002). Nonlinear kernel density estimation for binned data: convergence in entropy. Bernoulli, 8(4), 423–449.

###  Bandwidth selection (R base stats)
	- R stats documentation: Bandwidth Selectors for Kernel Density Estimation (e.g., bw.nrd, bw.SJ, etc.).
(Includes references such as Silverman (1986), Scott (1992), and Sheather & Jones (1991).)
