# Networks-Decide-for-Us

This repository contains the refactored code for the "Networks Decide for Us" project.

## Refactoring Plan

The goal is to refactor the project to use a modified version of the `netdiffuseR` package (which supports stochastic transmission).

### Requirements

- `netdiffuseR` (stochastic-transmission branch): `devtools::install_github("USCCANA/netdiffuseR", ref = "stochastic-transmission")`

### Tasks

- [ ] Update simulation scripts (`04_*.R`) to use `rdiffnet(..., exposure.mode = "stochastic")`.
- [ ] Update "Rational" agents logic: set threshold to $\epsilon$ (approx 0).
