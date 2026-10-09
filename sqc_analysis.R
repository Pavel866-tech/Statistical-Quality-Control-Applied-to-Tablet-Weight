# Import the packages and the csv data
library(qcc)
library(ggplot2)
library(dplyr)
library(forecast)
library(moments)

tablet_weight <- read.csv("data/tablet_weight.csv")

# Check the data
head(tablet_weight)
tail(tablet_weight)

# Basic Data Audit
glimpse(tablet_weight)

cat("Dimensions (rows, columns):", dim(tablet_weight), "\n")

if (anyNA(tablet_weight)) {
  cat("⚠ Missing values detected in the dataset.\n")
} else {
  cat("✓ No missing values. Summary:\n\n")
  data.frame(
    col    = names(tablet_weight),
    type   = sapply(tablet_weight, class),
    unique = sapply(tablet_weight, \(x) length(unique(x)))
  )
}

# Get rid of X column and change the column names
tablet_weight <- tablet_weight |>
  select(-X) |>
  rename(
    target_value = Weight.target.value,
    weight_value = Weight.value
  )

glimpse(tablet_weight)

# Run Shewhart I-MR Chart
weight_value <- tablet_weight$weight_value
phase1 <- weight_value[1:100]
phase2 <- weight_value[101:length(weight_value)]

par(mar = c(9, 4, 4, 2))
imr <- qcc(
  data = phase1,
  newdata = phase2,
  type = "xbar.one",
  title = "Tablet Weight - I-MR Chart",
  xlab = "Observation",
  ylab = "Weight (mg)",
  add.stats = TRUE
)

summary(imr)

imr_beyond_limits <- imr$violations$beyond.limits
imr_violation_runs <- imr$violations$violating.runs

imr_violations <- data.frame(
  Type  = c("Beyond limits", "Violating runs", "Total (unique)"),
  Count = c(length(imr_beyond_limits),
            length(imr_violation_runs),
            length(unique(c(imr_beyond_limits,
                            imr_violation_runs))))
)
cat("I-MR Chart — Violations Summary\n\n")
imr_violations

# Histogram of Tablet Weight Values
optimal_bins <- round(1 + log(length(phase1), base = 2))
ggplot() +
  aes(x = phase1) +
  geom_histogram(bins = optimal_bins, fill = "steelblue", color = "black") +
  geom_density(fill = "blue", alpha = 0.5) +
  geom_vline(aes(xintercept = mean(weight_value),
                 color = "Mean"),
             linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = median(weight_value), 
                 color = "Median"),
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(name   = NULL,
                     values = c("Mean"   = "red",
                                "Median" = "darkgreen")) +
  labs(title = "Distribution of Tablet Weight",
       x = "Weight (mg)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "right")

# Box-plot of Tablet Weight Values
ggplot() +
  aes(x = phase1) +
  geom_boxplot(fill = "steelblue", outlier.color = "red", outlier.size = 2) +
  labs(title = "Tablet Weight — Boxplot", x = "Weight (mg)") +
  theme_minimal()

out_vals <- boxplot.stats(weight_value)$out
out_pos <- which(weight_value %in% out_vals)

data.frame(
  n_outliers = length(out_vals),
  Positions  = paste(out_pos, collapse = ", "),
  Values     = paste(round(out_vals, 2), collapse = ", "),
  Skewness   = round(skewness(phase1), 2)
)

# Autocorrelation Diagnosis
acf(phase1,
    lag.max = 20,
    type     = "correlation",
    plot     = TRUE,
    main = "ACF Plot of Tablet Weight"
    )

Box.test(phase1, lag = 20, type = "Ljung-Box")

# Normality Diagnosis
ggplot() +
  aes(sample = phase1) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot — Phase 1 Tablet Weight",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles (mg)") +
  theme_minimal()

shapiro.test(phase1)

# ARIMA Model
arima_model <- auto.arima(phase1)
arimaorder(arima_model)
residuals <- as.numeric(arima_model$residuals)

# ACF Plot for ARIMA Residuals
acf(residuals,
    lag.max = 20,
    type     = "correlation",
    plot     = TRUE,
    main = "ACF Plot of Residuals of ARIMA Model"
    )
Box.test(residuals, lag = 20, type = "Ljung-Box")

# Normality Test for ARIMA Residuals
ggplot() +
  aes(sample = residuals) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot — Residuals of ARIMA Model",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

shapiro.test(residuals)

# EWMA Plot
ewma_chart <- ewma(residuals,
                   lambda     = 0.1,
                   nsigmas    = 2.7,
                   add.stats = TRUE,
                   xlab = "Observation",
                   ylab = "Exponentially Weighted Residual"
                   )

summary(ewma_chart)

# EWMA Chart After Removing the Outlier
new_arima_model <- auto.arima(phase1[-out_pos])
arimaorder(new_arima_model)
new_residuals <- as.numeric(new_arima_model$residuals)
ewma_chart <- ewma(new_residuals,
                   lambda     = 0.1,
                   nsigmas    = 2.7,
                   add.stats = TRUE,
                   xlab = "Observation",
                   ylab = "Exponentially Weighted Residual"
                   )

# Capability Analysis
target <- unique(tablet_weight$target_value)
lsl <- target*(1-0.10)
usl <- target*(1+0.10)
capability <- process.capability(imr, spec.limits = c(lsl, usl))
capability$indices

# Process Classification
stable <- length(ewma_chart$violations) == 0

threshold <- 1.33
capable <- all(capability$indices[c("Cp", "Cp_k"), "Value"] >= threshold)

type <- if (stable && capable) "A" else
  if (!stable && capable) "B" else
    if (stable && !capable) "C" else "D"

# Type of Process
descriptions <- c(
  A = "Stable and capable — ideal process.",
  B = "Capable but unstable — meets specifications but unpredictable.",
  C = "Stable but incapable — predictable but produces nonconforming units.",
  D = "Unstable and incapable — worst state."
)

cat("Type:", type, "\n")
cat("Summary:", descriptions[type], "\n")