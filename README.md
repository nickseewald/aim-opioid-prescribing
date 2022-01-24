# Supplementary Code for "Effects of state opioid prescribing laws on use of opioid and other pain treatments among commercially insured U.S. adults" (McGinty, et al. 2022)

This repository contains two files:
- [`analysis_example.R`](analysis_example.R), which contains sample code for performing the augmented synthetic controls method for a single state and a single outcome, and also produces plots and summary statistics.
- [`xaugsynth.R`](xaugsynth.R) contains functions which extend the [augsynth](https://github.com/ebenmichael/augsynth) package to produce plots of observed and predicted potential outcomes for treatment states and their synthetic controls. Also contains code to compute average outcomes in the pre- and post-treatment periods as well as confidence intervals for those averages.

## Usage
Download the two .R files in this repository and store them in the same location before running `analysis_example.R`. `xaugsynth.R` can be used by itself alongside the augsynth package.
