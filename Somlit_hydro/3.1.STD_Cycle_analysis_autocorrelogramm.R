##### Cycles analyses with autocorrelogramms #####

# Load the libraries
library(dplyr)
library(tidyverse)
library(ggplot2)

# Load the dataframes
table_reg_stl_filled_10 <- read_tsv(file="data/table_reg_stl_filled_10.tsv.gz")

###### Step 1 :  Remove the trends #####

#Pivot it
ts<- table_reg_stl_filled_10%>%select(name, target_date,deseason_filled)%>%pivot_wider(names_from="name",
                                                                                    values_from="deseason_filled")

# Start in 1992 
ts <- ts %>%filter(target_date > "1992-01-01")
#Remove the trends
residuals_table_corr <- ts%>%
  dplyr::select(target_date, CHLA, NO3,T, S, SIOH4) %>%
  pivot_longer(cols = c(CHLA, NO3, T, S, SIOH4),
               names_to = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  mutate(residual = resid(lm(value ~ target_date, data = cur_data()))) %>%
  ungroup()


r<- residuals_table_corr%>%select(-value)%>%
  pivot_wider(names_from="variable", values_from="residual")


##### Step 2 : Plot it to vizualize the cycles #####

# List of all variables
variables <- colnames(r)[colnames(r) != "target_date"]

# ACF plot for all the variables
for (var in variables) {
  ts_data <- r[[var]]
  cat("===> ACF for:", var, "\n")
  
  # max length = half of the time series 
  lag_max <- floor(length(ts_data) / 2)
  
  # Plot it 
  acf(ts_data, lag.max = lag_max, main = paste("Autocorrelogram of", var))
}


##### Step 3 : Calculate the length of the cycle #####

# Calculate ACF 
acf_res <- acf(ts$S, lag.max = length(ts$S)/2, plot = FALSE)

# Extract lags and ACF
lag_vals <- acf_res$lag[-1]
acf_vals <- acf_res$acf[-1]

# Moving average 
k <- 7
acf_vals_smoothed <- stats::filter(acf_vals, rep(1/k, k), sides = 2)

# First negative down
first_negative_index <- which(acf_vals_smoothed < 0)[1]

# Look at the first positive peak after this
sub_indices <- seq(first_negative_index + 1, length(acf_vals_smoothed))
peak_index_in_subset <- which.max(acf_vals_smoothed[sub_indices])
peak_index <- sub_indices[peak_index_in_subset]

# Find the lag
peak_lag <- lag_vals[peak_index]

cat("Peak is at lag =", peak_lag, "\n")

# Convert it in years
period_weeks <- peak_lag * 2
period_years <- period_weeks / 52
cat("=> Length of a cycle  ~", round(period_years, 2), "years\n")

# Look at it 
plot(lag_vals, acf_vals, type = "h",
     main = "Autocorrelogram", ylab = "ACF", xlab = "Lag")
lines(lag_vals, acf_vals_smoothed, col = "red", lwd = 2)
abline(v = peak_lag, col = "blue", lty = 2) 
