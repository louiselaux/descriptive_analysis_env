
##### Fourier analyses #####

# Load the libraries
library(dplyr)
library(tidyverse)
library(ggplot2)

# Load the dataframes
table_reg_stl_filled_10 <- read_tsv(file="data/table_reg_stl_filled_10.tsv.gz")

###### Step 1 :  Remove the trends #####

# Pivot it
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


###### Step 2 : Fourier transform analysis #####

# List of variables
variables <- setdiff(names(r), "target_date")

# Loop on the variables
for (var in variables) {
  cat("===> FFT for the variable:", var, "\n")
  
  ts_data <- r[[var]]
  
  # Remove NA
  ts_data <- ts_data[!is.na(ts_data)]
  
  # Should not be NA but just in case 
  n <- length(ts_data)
  if (n < 10) {
    cat("Not enough points for", var, "\n")
    next
  }
  
  # FFT
  T_fft <- fft(ts_data)
  amplitude <- Mod(T_fft)[1:(n/2)]
  freq <- (0:(n/2 - 1)) / n
  freq_annee <- freq * 26  # 26 step per year
  
  freq_annee <- freq_annee[-1]
  amplitude <- amplitude[-1]
  
  period_years <- 1 / freq_annee
  
  # Plot
  plot(period_years, amplitude, type = "h",
       main = paste("Frequency spectrum", var),
       xlab = "Frequence (cycles/an)", ylab = "Amplitude")
}



##### Step 3 : Fourier analysis differently #####

# variables 
ts_data <- r$T  

# Remove NA
ts_data <- ts_data[!is.na(ts_data)]

# Spectrum
spec <- spectrum(ts_data, spans = c(5,5), plot = FALSE, ci = 0.95)

# Conversion fequency into period
freq <- spec$freq
period_years <- 1 / (freq * 26)  # 26 per year
power <- spec$spec / sum(spec$spec) # to normalize

# Plot
plot(period_years, power, type = "l", log = "y",
     xlab = "Period (years)", ylab = "Spectrum", main = "Spectrum",
     xlim = c(0, 25))

##### On the whole time series######

table_reg_stl_10 <- dstl %>%
  group_by(name) %>%
  mutate(raw_filled = castr::interpolate(x = target_date,
                                              y = raw,
                                              xout = target_date))

ts<- table_reg_stl_10%>%select(name, target_date,raw_filled)%>%pivot_wider(names_from="name",
                                                                                       values_from="raw_filled")

ts <- ts %>% filter(target_date > "1997-01-01")
ts_data <- ts$T  

# Remove NA
ts_data <- ts_data[!is.na(ts_data)]

# Spectrum
spec <- spectrum(ts_data, spans = c(6,6), plot = TRUE, ci = 0.95)

# Convert to years
max(spec$spec)

m <- spec$freq[spec$spec==max(spec$spec)]
period_years <- 1 / (m * 26)

# Conversion fequency into period
freq <- spec$freq
period_years <- 1 / (freq * 26)  # 26 per year
power <- spec$spec / sum(spec$spec) # to normalize

# Plot
plot(period_years, power, type = "l", log = "y",
     xlab = "Period (years)", ylab = "Spectrum", main = "Spectrum",
     xlim = c(0, 25))

