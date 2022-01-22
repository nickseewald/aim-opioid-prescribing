# Install gplots if not already installed; used for color matching confidence
# band to synthetic control line in potentialplot()
if (!("gplots" %in% installed.packages()))
  install.packages('gplots')

# Function for plotting potential outcomes
## augsynth: object of class "augsynth", returned from augsynth() function
## synthParams: optional *named* list of base R graphical parameters to control
##              the line for the synthetic control outcome
## ciParams: optional *named* list of base R graphical parameters to control the
##           confidence interval
## legendParams: optional *named* list of base R graphical parameters to control
##               the plot legend
## ...: optional arguments used to plot the treated unit line passed to plot()
##
## Invisibly returns a data frame containing all information needed to plot 
## using other software (e.g., Excel) if exported
potentialplot <- function(augsynth, 
                          synthParams = list(col = "blue", lty = 2),
                          ciParams = list(border = NA), 
                          legendParams = list(x = "bottom", horiz = T, cex = .8),
                          ...) {
  args      <- list(...)
  summ      <- summary(augsynth)
  att       <- summ$att[summ$att$Time >= augsynth$t_int, ]
  synth     <- predict(augsynth)
  combData  <- with(augsynth$data, cbind(X, y))
  txOutcome <- combData[augsynth$data$trt == 1, ]
  
  ci   <- data.frame(
    "lb" = txOutcome[augsynth$data$time >= augsynth$t_int] - att$upper_bound,
    "ub" = txOutcome[augsynth$data$time >= augsynth$t_int] - att$lower_bound
  )
  
  if (is.null(args$xlab))
    args$xlab <- "Time"
  if (is.null(args$ylab))
    args$ylab <- "Potential Outcome"
  if (is.null(args$ylim))
    args$ylim <- range(c(synth, txOutcome, ci$lb, ci$ub))
  
  # Set some default graphical parameters in case they're not provided
  if (is.null(synthParams$lty))
    synthParams$lty <- 2
  if (is.null(synthParams$lwd))
    synthParams$lwd <- 1
  if (is.null(ciParams$col))
    ciParams$col <- paste0(gplots::col2hex(synthParams$col), "33")
  if (is.null(ciParams$border))
    ciParams$border <- NA
  if (is.null(legendParams$legend))
    legendParams$legend <- c("Treated Unit", "Synthetic Control")
  if (is.null(args$col))
    args$col <- "black"
  if (is.null(args$lty))
    args$lty <- 1
  if (is.null(args$lwd))
    args$lwd <- 1
  if (is.null(legendParams$x))
    legendParams$x <- "bottom"
  if (is.null(legendParams$cex))
    legendParams$cex <- .8

  # Start with a blank plot
  plot(txOutcome ~ augsynth$data$time,
       type = "n",...)
  
  # Add confidence band around post-period synthetic control
  do.call(polygon, c(ciParams, "x" = list(c(att$Time, rev(att$Time))),
          "y" = list(c(ci$lb, rev(ci$ub)))))
  
  # Add vertical line at treatment time
  abline(v = augsynth$t_int, lty = "dotted")
  
  # Add treated line
  do.call(lines, c(txOutcome ~ augsynth$data$time, args))
  
  # Add synthetic control line
  do.call(lines, c(synth ~ augsynth$data$time, synthParams))
  
  # Add legend
  do.call(legend, c(legendParams,
                    lty = list(c(args$lty, synthParams$lty)),
                    lwd = list(c(args$lwd, synthParams$lwd)),
                    col = list(c(args$col, synthParams$col))))
  
  # Create data frame to be returned invisibly
  dat <- data.frame("time" = augsynth$data$time,
                    "treated_unit" = txOutcome,
                    "synth_ctrl" = synth)
  dat <- merge(dat, 
               cbind("time" = augsynth$data$time[augsynth$data$time >=
                                                   augsynth$t_int], 
                     ci), by = "time", all = T)
  
  invisible(dat)
}

# Function to compute averages by group (and average ATT)
averageOutcomes <- function(augsynth, conf.level = .95) {
  
  # Extract some information from the augsynth object
  synth <- predict(augsynth)
  att   <- summary(augsynth, alpha = 1 - conf.level)$att
  X     <- augsynth$data$X[augsynth$data$trt == 1, ]
  y     <- augsynth$data$y[augsynth$data$trt == 1, ]
  time  <- augsynth$data$time
  t_int <- augsynth$t_int
  T0    <- sum(time < t_int)
  T1    <- sum(time >= t_int)
  
  # Compute appropriate t critical values for confidence bounds
  tstar0 <- qt((1 + conf.level) / 2, df = T0 - 1)
  tstar1 <- qt((1 + conf.level) / 2, df = T1 - 1)
  
  # Compute an "average" SE for the overall/average ATT by assuming individual
  # (i.e., time-specific) ATTs have approximately t-distributed CIs and solving
  # for each SE, then pooling
  att_se <- with(subset(att, Time >= t_int),
                 sqrt(mean((c((Estimate - lower_bound) / qnorm((1 + conf.level) / 2),
                        (upper_bound - Estimate) / qnorm((1 + conf.level) / 2))^2))))
  
  # Construct a data frame to return
  d <- data.frame("Estimate" = c(mean(X), mean(y), 
                                 mean(synth[time < t_int]),
                                 mean(synth[time >= t_int]),
                                 mean(att$Estimate[time >= t_int])),
                  "lower_bound" = c(
                    mean(X) - tstar0 * sd(X) / sqrt(T0),
                    mean(y) - tstar1 * sd(y) / sqrt(T1),
                    mean(synth[time < t_int]) - tstar0 * 
                      sd(synth[time < t_int]) / sqrt(T0),
                    mean(synth[time >= t_int]) - tstar1 * att_se,
                    mean(att$Estimate[time >= t_int]) - tstar1 * att_se
                  ), 
                  "upper_bound" = c(
                    mean(X) + tstar0 * sd(X) / sqrt(T0),
                    mean(y) + tstar1 * sd(y) / sqrt(T1),
                    mean(synth[time < t_int]) + tstar0 * 
                      sd(synth[time < t_int]) / sqrt(T0),
                    mean(synth[time >= t_int]) + tstar1 * att_se,
                    mean(att$Estimate[time >= t_int]) + tstar1 * att_se
                  ))
  
  rownames(d) <- c("Treated, pre", "Treated, post",
                   "Synthetic, pre", "Synthetic, post", "ATT")
  
  colnames(d) <- c("Estimate", paste0(100 * conf.level, "% Lower Bound"),
                   paste0(100 * conf.level, "% Upper Bound"))
  
  class(d) <- c("averageOut", class(d))
  
  d
}

# Print method for results of averageOutcomes()
print.averageOut <- function(avo, digits = 3) {
  avo <- as.data.frame(avo)
  print(round(avo, digits))
}
