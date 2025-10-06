##### Load the libraries #####

library(castr)
library(morphr)
library(stlplus)
library(readr)
library(dplyr)
library(tidyverse)
library(trend)
library(zoo)
library(missMDA)
library(broom)
library(nlme)

##### Load the data table#####

dstl <- read_tsv("data/dstl_interpolated_before.tsv")
head(dstl)
stats <- read_tsv("data/stats_interpolated_before.tsv")
head(stats)

##### Is summer becoming longer? ####
##### Raw data - gls regression #####

# Get the trend for each date
dstl <- dstl %>%
  left_join(stats %>% select(name, gls_intercept, gls_slope), by = "name") %>% group_by(name)%>%
  mutate(
    trend_gls = gls_intercept + gls_slope * as.numeric(target_date)
  )

# Remove the trend from the raw data

dstl <- dstl %>% group_by(name)%>%
  mutate(
    detrended_raw = raw - trend_gls
  )


# Interpolation at the day resolution
daily_grid <- tibble(date = seq(min(dstl$target_date, na.rm = TRUE),
                                max(dstl$target_date, na.rm = TRUE),
                                by = "1 day"))



# Daily interpolation

interp_daily <- dstl %>%
  group_by(name) %>%
  group_modify(~{
    xg <- .x %>%
      filter(!is.na(target_date), !is.na(detrended_raw)) %>%
      arrange(target_date)
    
    if (nrow(xg) < 2) {
      return(tibble(date = as.Date(character()), detrended = numeric()))
    }
    
    grid <- tibble(date = seq(min(xg$target_date), max(xg$target_date), by = "1 day"))
    
    # Interpolation
    y_interp <- castr::interpolate(
      x    = xg$target_date,
      y    = xg$detrended_raw,
      xout = grid$date
    )
    
    tibble(date = grid$date, detrended = y_interp)
  }) %>%
  ungroup()



head(interp_daily)

# Small test
target_date <- as.Date(c("2020-05-01", "2020-05-15"))
detrended_raw <- c(-3.0, 1.0) 

# Sequence of dates
daily_grid <- tibble(date = seq(min(target_date), max(target_date), by = "1 day"))

# Linear interpolation
interp_daily_test <- castr::interpolate(
  x    = target_date,
  y    = detrended_raw,
  xout = daily_grid$date
)

interp_daily_test <- tibble(date = daily_grid$date, detrended = interp_daily_test)

# Plot
ggplot() +
  #
  geom_point(data = interp_daily_test, aes(date, detrended), color = "purple") +
  # Points at the truth date
  geom_point(data = tibble(date = target_date, detrended = detrended_raw),
             aes(date, detrended), color = "red", size = 3, shape = 17) +
  labs(title = "Daily interpolation between two dates",
       y = "", x = "Date") +
  theme_bw()


# Test with one variable

# Choose the variable and the interval
var <- "T"
date_min <- as.Date("2000-01-04")
date_max <- as.Date("2000-01-20")

ggplot() +
  # Daily interpolation
  geom_point(
    data = interp_daily %>%
      filter(name == var, date >= date_min, date <= date_max),
    aes(date, detrended),
    color = "purple"
  ) +
  # Observed points
  geom_point(
    data = dstl %>%
      filter(name == var, target_date >= date_min, target_date <= date_max),
    aes(target_date, detrended_raw),
    color = "red", size = 2, shape = 17
  ) +
  labs(
    title = paste("Daily interpolation", var, "between", date_min, "&", date_max),
    x = "Date",
    y = "Detrended"
  ) +
  theme_bw()


##### Is summer becoming longer ? ######


###  Way 1 : Nb of fays above a threshold
head(interp_daily)

interp_daily <- interp_daily %>%filter(name=="T")%>%
  arrange(date) %>%
  mutate(
    detrended_smooth = rollmean(detrended, k = 20, fill = NA, align = "center") 
  )


# Check the smoothing
interp_daily %>% filter(name=="T")%>% mutate(year=year(date))%>%filter(year %in% c("2000"))%>%ggplot()+
  geom_line(aes(y = detrended, x= date), color = "darkgreen") +
  geom_line(aes(y = detrended_smooth, x=date), color = "blue", linewidth = 1) +
  geom_hline(yintercept = thresh_rel, color = "red", linetype = "dashed") +
  labs(title = "",
       x = "Date", y = "") +
  theme_bw()

# What is the threshold?
thresh_rel <- quantile(interp_daily$detrended_smooth, probs = 0.75, na.rm = TRUE)
thresh_rel

hist(interp_daily$detrended_smooth) 
abline (v= thresh_rel, col="red")


# Nb of days per year
stats_interp <- interp_daily %>% filter(name=="T")%>%mutate(year=year(date))%>%
  group_by(year)%>% filter(detrended_smooth>thresh_rel)%>%summarise(count=n())

stats_interp <- interp_daily %>% mutate(year=year(date))%>%
  group_by(year)%>% summarise(days_above=sum(detrended_smooth>thresh_rel, na.rm=TRUE))

# Statistical test
mean_days <- mean(stats_interp$days_above, na.rm = TRUE)
var_days  <- var(stats_interp$days_above,  na.rm = TRUE)

mean_days
var_days

model_stat <- glm(days_above ~ year,
                  data = stats_interp,
                  family = quasipoisson(link = "log"))

summary(model_stat)

# Plot one to see

year_test <- 1995

ggplot(interp_daily %>% filter(year(date) == year_test)) +
  geom_line(aes(date, detrended_smooth), color = "purple") +
  geom_hline(yintercept = thresh_rel, linetype = "dashed", color = "red") +
  geom_rug(data = interp_daily %>% 
             filter(year(date) == year_test, detrended_smooth > thresh_rel),
           aes(x = date), sides = "b", color = "red") +
  labs(title = paste("Year", year_test, "- days > threshold =", 
                     sum(interp_daily$detrended_smooth[year(interp_daily$date) == year_test] > thresh_rel, na.rm=TRUE)),
       y = "Year smoothed", x = "Date") +
  theme_bw()

# Plot the model
ggplot(stats_interp, aes(x = year, y = days_above)) +
  geom_point(color = "purple") +
  geom_smooth(
    method = "glm",
    method.args = list(family = quasipoisson(link = "log")),
    se = TRUE,   
    color = "red"
  ) +
  labs(
    title = "",
    x = "Year",
    y = "Days above threshold"
  ) +
  theme_bw()



###### On the raw data #####

interp_daily_raw <- dstl %>%
  group_by(name) %>%
  group_modify(~{
    xg <- .x %>%
      filter(!is.na(target_date), !is.na(raw)) %>%
      arrange(target_date)
    if (nrow(xg) < 2) return(tibble(date = as.Date(character()), raw = numeric()))
    grid <- tibble(date = seq(min(xg$target_date), max(xg$target_date), by = "1 day"))
    y_interp <- castr::interpolate(
      x    = xg$target_date,
      y    = xg$raw,
      xout = grid$date
    )
    tibble(date = grid$date, raw = y_interp)
  }) %>%
  ungroup()

head(interp_daily_raw)

# Focus on temperature, smooth and count

interp_T_raw <- interp_daily_raw %>%
  filter(name == "T") %>%
  arrange(date) %>%
  mutate(raw_smooth = rollmean(raw, k = 20, fill = NA, align = "center"))

thresh_rel_raw <- quantile(interp_T_raw$raw_smooth, probs = 0.75, na.rm = TRUE)

stats_interp_raw <- interp_T_raw %>%
  mutate(year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(days_above = sum(raw_smooth > thresh_rel_raw, na.rm = TRUE), .groups = "drop")

model_stat_raw <- glm(days_above ~ year,
                      data = stats_interp_raw,
                      family = quasipoisson(link = "log"))

summary(model_stat_raw)



# Compare the two
stats_det <- stats_interp %>%            
  mutate(approach = "Detrended (Q75)")

stats_raw2 <- stats_interp_raw %>%       
  mutate(approach = "Raw (Q75)")

common_years <- intersect(stats_det$year, stats_raw2$year)

stats_both <- bind_rows(
  stats_det  %>% filter(year %in% common_years),
  stats_raw2 %>% filter(year %in% common_years)
) %>%
  arrange(approach, year) %>%
  tidyr::drop_na(days_above)

ggplot(stats_both, aes(x = year, y = days_above, color = approach)) +
  geom_point(alpha = 0.8) +
  geom_smooth(
    method = "glm",
    method.args = list(family = quasipoisson(link = "log")),
    se = TRUE
  ) +
  labs(x = "", y = "",
       title = "") +
  theme_bw() +
  theme(legend.title = element_blank())
