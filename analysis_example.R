## Install "remotes" package to download augsynth from GitHub
if (!("remotes" %in% installed.packages()))
  install.packages('remotes')

## Install augsynth from Github (available at https://github.com/ebenmichael/augsynth).
## A list of package dependencies is available at 
## https://raw.githubusercontent.com/ebenmichael/augsynth/master/DESCRIPTION.
## The below code will install the package as it existed when analyses were performed.
if (!("augsynth" %in% installed.packages()))
  remotes::install_github("ebenmichael/augsynth", dependencies = TRUE,
                          ref = "ab71fc726a32fe084b200be7b21469be84f43500")
## Load augsynth
library(augsynth)

## Load ggplot2 for plotting
library(ggplot2)

## Load data
load("ohio_cohort1.Rdata")

## Source helper functions to extend augsynth; see GitHub for code
source("xaugsynth.R")

## Run augmented synthetic control method using ridge regression to augment, as 
## recommended by Ben-Michael, et al. (2021). 
covsyn <-
  augsynth(
    opioid_any ~ policy |
      age + female + mi + sud + elix + fpl_totalpop + employ_age20_4 + 
      edu_hs_grad + raceeth_blacknohs + raceeth_hispanic,
    unit = egeoloc,
    time = month,
    data = ohio_cohort1,
    progfunc = "ridge",  #use ridge regression to impute control outcomes
    scm = TRUE,          #use SCM weighting function
    fixedeff = TRUE      #include unit fixed effects
  )

## Compute average outcomes in pre- and post-law periods, using code from 
## xaugsynth.R (see GitHub)
print(averageOutcomes(covsyn), digits = 4)

## Plot for Online Supplement C, using code from xaugsynth.R (see GitHub)
potentialplot(covsyn, 
              xlab = "Month (July 1, 2009-June 30, 2013)",
              ylab = "% of patients with any opioid Rx, per month",
              main = "Ohio: Proportion patients with any opioid Rx, per month",
              ylim = c(0, 0.1),
              legendParams = list(legend = c("Ohio", "Synthetic Control"),
                                  horiz = T, cex = 1))

## Plot for Online Supplement D
plot(covsyn) +
  labs(x = "Month (July 1, 2009-June 30, 2013",
       y = "Percentage point change") + 
  expand_limits(x = c(0, 50), y = c(-5, 5))


