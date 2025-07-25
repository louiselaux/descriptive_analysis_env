##### Trend analysis #####
#Load the libraries

#devtools::install_github("jiho/castr")
library(castr)
library(morphr)
library(stlplus)
library(readr)
library(dplyr)
library(tidyverse)
library(trend)
library(imputeTS)     # for na_seadec and na_interpolation functions
library(mar1s)        # for seasonal.smooth and create.fourier.basis functions

#Load the files
std <- read_tsv("data/std_all_depths.tsv")

#Pivot it here
std_pivot <- std %>%pivot_wider(names_from = "name", values_from="value" )
std_pivot<- as.data.frame(std_pivot)

#Plot relationships between variables
#env_mean_pivot%>%dplyr::select(-date)%>%ggpairs()

#First step: Regularization

# Define a sequence of dates 
ref <- tibble(
  target_date=seq(from=as.Date("1967-01-05"), to=as.Date("2022-12-31"), by=14),
  year=year(target_date)
)

#Start in 1992 if you want to be consistent with CCM analyses but you can also remove this line if you want to have it on the whole time series
#ref<- ref%>% filter(target_date>"1992-01-16")

# identify years in which the number of obs is larger than usual
pbs <- ref %>%
  count(year) %>%
  filter(n>26)

# ->this is often an extra date in very late decembre => just remove it
ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
ref %>%
  count(year) %>%
  filter(n>26)
# -> all OK

# Match data based on these reference dates
avail_dates <- unique(std_pivot$date)
ref <- ref %>%
  mutate(
    closest_date = castr::closest(ref$target_date, avail_dates),
    date_diff = abs(closest_date - target_date) %>% as.numeric()
  )

# Insert the data based on the matched dates
table_reg <- left_join(ref, std, by=c("closest_date"="date"), relationship="many-to-many")

# erase data for matches that are too far from the target
table_reg<- table_reg %>%
  mutate(value = if_else(date_diff > 6, NA, value))

ggplot(table_reg) + facet_wrap(~name, scales="free_y") +
  geom_point(aes(x=target_date, y=value), size=0.2) + theme_bw()+labs(x="date", y="value")
# -> OK

#Look at the distribution of variables to remove extreme values
t<-table_reg%>%pivot_wider(names_from = "name", values_from="value")
hist(t$CHLA)
hist(t$NO3)
hist(t$S)
hist(t$T)
hist(t$SIOH4)

#t<- t%>%mutate(CHLA=mask_extreme(CHLA, percent=c(0,0.5)),
#               S=mask_extreme(S, percent=c(0.2,0)))

t<- t%>%mutate(CHLA=log10(CHLA),
               S=mask_extreme(S, percent=c(0.2,0)))

table_reg_stl<- t %>%pivot_longer(cols=c("CHLA","NO3","S","SIOH4","T"))



#####Second step: STL decomposition ######
table_reg_stl<- table_reg_stl%>%select(-year)

##### At the depth you choose #####
x<-10  # Put the depth you want
table_reg_stl <-  table_reg_stl %>%filter(depth==x)

temp_1m <- table_reg_stl %>%
  filter(name == "S") %>%
  dplyr::select(target_date, closest_date, date_diff, value,name)


##### First : Change the values of the time windows #####

# Time windows you want to look at
t_windows <- c(25, 20, 15, 10, 5)

# STL 
stl_tests <- lapply(t_windows, function(tw) {
  dec <- stlplus(temp_1m$value, temp_1m$target_date, 
                 n.p = 26,                   # 26 per year
                 s.window = "periodic", 
                 t.window = tw * 26)         # How many years 
  
  data.frame(
    target_date = temp_1m$target_date,
    trend = dec$data$trend,
    seasonal = dec$data$seasonal,
    t_window = paste0("tw_", tw)
  )
})

# Put everything together
trends <- bind_rows(stl_tests) %>%
  pivot_longer(cols = c(trend, seasonal), names_to = "component", values_to = "value")

# Plot trends and seasonal variations
ggplot(trends, aes(x = target_date, y = value, color = t_window)) +
  geom_line() +
  facet_wrap(~component, ncol = 1, scales = "free_y") +
  labs(
    title = "",
    x = "Date", y = "Value",
    color = "T.window"
  ) +
  theme_bw(base_size = 14)



##### Second step : Change the seasonal values #####
s_windows <- c(1,3,5)  # in years

# Multiplication per 26 as 26 per year
s_window_pts <- s_windows * 26

# STL with fixed t window
stl_tests <- lapply(s_window_pts, function(sw) {
  dec <- stlplus(temp_1m$value, temp_1m$target_date,
                 n.p = 26,
                 s.window = sw,         # Seasonalisation changes
                 t.window = 10 * 26)    # Fixed 
  
  data.frame(
    target_date = temp_1m$target_date,
    trend = dec$data$trend,
    seasonal = dec$data$seasonal,
    s_window = paste0("sw_", sw/26, "years")
  )
})

# Bind all
trends <- bind_rows(stl_tests)

# Plot trends
ggplot(trends, aes(x = target_date, y = trend, color = s_window)) +
  geom_line(size = 1) +
  labs(
    title = "How does s.season impact the trend component",
    x = "Date", y = "STL trend",
    color = "s.window (years)"
  ) +
  theme_bw(base_size = 14)

# Plot season
ggplot(trends, aes(x = target_date, y = seasonal, color = s_window)) +
  geom_line(size = 1) +
  labs(
    title = "Impact of s window on seasonal component",
    x = "Date", y = "Seasonal component (°C)",
    color = "s.window (years)"
  ) +
  theme_bw(base_size = 14)


#### Small tests ####
dd<- stlplus(temp_1m$value, temp_1m$target_date,
             n.p = 26,
             s.window = 10*26,         # Seasonalisation changes
             t.window = 3* 26) 

plot(dd)


##### Comparisons of the two ways of doing it #####
# When it is periodic
dec_periodic <- stlplus(temp_1m$value, temp_1m$target_date,
                        n.p = 26, s.window = "periodic", t.window = 10 * 26)

# When the season is 5 years
dec_flexible <- stlplus(temp_1m$value, temp_1m$target_date,
                        n.p = 26, s.window = 5 * 26, t.window = 10 * 26)

# Comparison of the two
df_compare <- data.frame(
  date = temp_1m$target_date,
  periodic = dec_periodic$data$seasonal,
  flexible = dec_flexible$data$seasonal
)

ggplot(df_compare) +
  geom_line(aes(x = date, y = periodic, color = "Periodic")) +
  geom_line(aes(x = date, y = flexible, color = "Flexible 5 years")) +
  labs(title = "",
       x = "Date", y = "Seasonal component (°C)", color = "s.window") +
  theme_bw()

##### Is summer becoming longer? #####
dd<- stlplus(temp_1m$value, temp_1m$target_date,
             n.p = 26,
             s.window = 3*26,         # Seasonalisation changes
             t.window = 10* 26) 

plot(dd)

df_test <- data.frame(
  date = temp_1m$target_date,
  periodic = dd$data$seasonal
)


df_test%>% mutate(year=year(date))%>% filter(year%in% c(1977:1997))%>%
  ggplot()+ geom_line(aes(x=date,y=periodic))+theme_bw()


df_test%>% mutate(year=year(date))%>% filter(year%in% c(1997:2017))%>%
  ggplot()+ geom_line(aes(x=date,y=periodic))+theme_bw()



##### Now, kind of grid search with the two parameters that vary ####
# Parameters
t_windows <- c(25,10)     
s_windows <- c( 26*1, 26*2, 26*5) 

# All combinations possible
param_grid <- expand.grid(tw = t_windows, sw = s_windows, stringsAsFactors = FALSE)

# STL for each of them
stl_tests <- lapply(1:nrow(param_grid), function(i) {
  tw <- param_grid$tw[i]
  sw <- param_grid$sw[i]
  
  dec <- stlplus(temp_1m$value, temp_1m$target_date, 
                 n.p = 26,
                 s.window = sw, 
                 t.window = tw * 26)
  
  data.frame(
    target_date = temp_1m$target_date,
    trend = dec$data$trend,
    seasonal = dec$data$seasonal,
    t_window = paste0("tw_", tw),
    s_window = ifelse(sw == "periodic", "periodic", paste0("sw_", sw/26, "y"))
  )
})

# Put all the results together
trends <- bind_rows(stl_tests) %>%
  pivot_longer(cols = c(trend, seasonal),
               names_to = "component",
               values_to = "value")

# Plot trends and seasons
ggplot(trends, aes(x = target_date, y = value, color = t_window)) +
  geom_line() +
  facet_grid(component ~ s_window, scales = "free_y") +
  labs(
    title = "",
    x = "Date", y = "Value",
    color = "T.window"
  ) +
  theme_bw(base_size = 14)


###### Try with what Francesco sent #####
complete_data <- temp_1m   
value_col <- "value"
frequency <- 26  
#frequency = 365 # here you select in which time period the pattern is seen, for seasons, 365 days

# Assuming complete_data has a value column (like temp) and a date column (in date format, and without missing days)
# you can have missing values (you can impute them as below) but you can’t have missing days
# Create time series
x_ts <- ts(complete_data[[value_col]],
           frequency = frequency,
           start = min(complete_data$target_date))

# Seasonal decomposition and imputation
x_seas_imput <- imputeTS::na_seadec(x_ts, algorithm = 'mean')
x_seas <- rep(seasonal.smooth(x_seas_imput, basis = create.fourier.basis(nbasis = 10)),
              length.out = length(x_ts))
x_imput <- imputeTS::na_interpolation(x_ts - x_seas, maxgap = 14) + x_seas

# Create final dataset
complete_ts <- data.frame(
  date = complete_data$target_date,
  original = as.numeric(x_imput),
  seasonal_component = as.numeric(x_seas),
  deseasoned = as.numeric(x_imput - x_seas)
)
# the deseasoned time series is in complete_ts$deseasoned


# Plot
ggplot(complete_ts, aes(x = date)) +
  geom_line(aes(y = original, color = "Original")) +
  geom_line(aes(y = deseasoned, color = "Deseasoned")) +
  labs(
    title = "",
    x = "Date", y = "Temperature (°C)",
    color = "Série"
  ) +
  theme_bw()

ggplot(complete_ts, aes(x = date, y = seasonal_component)) +
  geom_line(color = "blue") +
  labs(title = "Seasonal component", x = "Date", y = "°C") +
  theme_bw()

##### Comparison of the 2 ways #####
# with STL
dec_stl_periodic <- stlplus(temp_1m$value, temp_1m$target_date,
                            n.p = 26,
                            s.window = "periodic",
                            t.window = 10 * 26)

seasonal_stl <- data.frame(
  date = temp_1m$target_date,
  seasonal_stl = dec_stl_periodic$data$seasonal
)

# With Fourier
value_col <- "value"
frequency <- 26  

x_ts <- ts(temp_1m[[value_col]],
           frequency = frequency,
           start = c(year(min(temp_1m$target_date)),
                     yday(min(temp_1m$target_date)) / 14))

x_seas_imput <- imputeTS::na_seadec(x_ts, algorithm = 'mean')
x_seas <- rep(seasonal.smooth(x_seas_imput, basis = create.fourier.basis(nbasis = 10)),
              length.out = length(x_ts))

seasonal_fourier <- data.frame(
  date = temp_1m$target_date,
  seasonal_fourier = as.numeric(x_seas)
)

# Put the two
df_compare <- seasonal_stl %>%
  left_join(seasonal_fourier, by = "date")

# Plot
ggplot(df_compare, aes(x = date)) +
  geom_line(aes(y = (seasonal_stl), color = "STL (periodic)")) +
  geom_line(aes(y = (seasonal_fourier), color = "Fourier"), linetype = "dashed") +
  labs(
    title = "",
    x = "Date",
    y = "Seasonal Component (°C)",
    color = "Method"
  ) +
  theme_bw()
