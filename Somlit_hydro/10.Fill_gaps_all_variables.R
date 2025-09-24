# Load the libraries
library(castr)
library(morphr)
library(stlplus)
library(readr)
library(dplyr)
library(tidyverse)
library(trend)
library(zoo)
library(missMDA)

#Load the files
std <- read_tsv("data/std_all_depths.tsv")

# Remove outliers
std <- std %>%
  mutate(value = ifelse(name == "S" & value < 37, NA, value))

# Pivot it here
std_pivot <- std %>%pivot_wider(names_from = "name", values_from="value" )
std_pivot<- as.data.frame(std_pivot)


##### Step 1: Interpolation through iterative PCA or linear interpolation depending on the length of the NA #####

head(std_pivot)


### Estimate the length of gaps ###
ref <- tibble(
  target_date = seq(from = as.Date("1967-01-05"), to = as.Date("2022-12-31"), by = 14),
  year = year(target_date)
)
pbs <- ref %>% count(year) %>% filter(n > 26)
ref <- ref %>% filter(!(year %in% pbs$year & yday(target_date) > 360)) %>% select(-year)

avail_dates <- unique(std_pivot$date)
ref <- ref %>%
  mutate(
    closest_date = castr::closest(target_date, avail_dates),
    date_diff = as.numeric(abs(closest_date - target_date)) # days 
  )

table_reg <- ref %>%
  left_join(std, by = c("closest_date" = "date"), relationship = "many-to-many")

table_reg <- table_reg %>% filter(name != "pH")

# Identification of gaps length

ts_all <- table_reg %>%
  #filter(depth == 10) %>%
  arrange(name, target_date, depth) %>%
  group_by(name, depth) %>%
  mutate(
    value_obs  = if_else(date_diff <= 3, value, NA_real_),
    is_na      = is.na(value_obs),
    is_na_next = lead(is_na, default = FALSE),
    int        = if_else(is_na, 1L, 0L),
    start_gap  = is_na & !lag(is_na, default = FALSE),
    nb_gap     = cumsum(start_gap),
    gap_id     = if_else(is_na, nb_gap, NA_integer_)
  ) %>%
  group_by(name, gap_id, depth) %>%
  mutate(size_gap = if_else(is.na(gap_id), 0L, n())) %>%
  ungroup()

# Check it seemed to have worked as expected
ts_all%>% select(size_gap, closest_date, name, value, target_date, depth)%>% filter(size_gap %in% c(2:20)) %>% filter(name == "T")

# Small gaps go to linear interpolation 

interp_linear <- ts_all %>%
  group_by(name, depth) %>%
  mutate(
    raw_value   = value,
    value_obs   = if_else(date_diff > 3, NA_real_, raw_value),
    value_interp = castr::interpolate(
      x    = closest_date[!is.na(raw_value)],
      y    = raw_value[!is.na(raw_value)],
      xout = target_date
    ),
    value_final = dplyr::coalesce(value_obs, value_interp),
    value_final = if_else(size_gap > 5, NA_real_, value_final)
  ) %>%
  ungroup()


# Plot to investigate

interp_linear %>% mutate(year=year(target_date))%>% filter(year==2009)%>% filter(depth =="50")%>%ggplot() + geom_point(aes(x=closest_date, y=raw_value), shape=15, colour="blue", size=2)+
  geom_point(aes(x=target_date,y=value_final), shape=12, color="red", size=2) + theme_bw() + facet_wrap ( ~name, scales = "free_y") 

# Taking into account the different depths 
interp_linear %>% mutate(year=year(target_date))%>% filter(year==2009)%>% filter(name=="T")%>%ggplot() + geom_point(aes(x=closest_date, y=raw_value), shape=15, colour="blue", size=2)+
  geom_point(aes(x=target_date,y=value_final), shape=12, color="red", size=2) + theme_bw() + facet_wrap ( ~depth, scales = "free_y") 


##### Iterative PCA for bigger gaps #####

interp_wide <- interp_linear %>%
  select(target_date, closest_date, name, value_final, depth) %>%
  pivot_wider(names_from = name, values_from = value_final,
              names_prefix = "lin_")

tablo_merged <- table_reg %>%
  pivot_wider(names_from = name, values_from = value) %>%
  #filter(depth == 10) %>%
  left_join(interp_wide, by = c("target_date","closest_date","depth"))


tablo_merged <- tablo_merged %>%
  mutate(
    T    = coalesce(T, lin_T),
    CHLA = coalesce(CHLA, lin_CHLA),
    NO3  = coalesce(NO3, lin_NO3),
    S    = coalesce(S, lin_S),
    O    = coalesce(O, lin_O),
    SIOH4= coalesce (SIOH4, lin_SIOH4),
    MES= coalesce(MES, lin_MES)
  ) %>%
  select(-starts_with("lin_"))


##### Problem regarding the seasonality so we need to take it into account #####
tab <- tablo_merged %>%
  transmute(
    target_date,
    depth,
    month = lubridate::month(target_date),
    T, CHLA, NO3, S, O, SIOH4, MES
  )


# Median per month per variable
clim <- tab %>% group_by(month, depth) %>% summarise(across(c(T, CHLA, NO3, S, O, SIOH4, MES), ~median(., na.rm = TRUE), .names = "{.col}_clim"), .groups = "drop") 


tabc <- tab %>% left_join(clim, by = c("month","depth")) 

#Anomalies= value - clim per month

tab_anom <- tabc %>% mutate(across(c(T, CHLA, NO3, S, O, SIOH4, MES), ~ . - get(paste0(cur_column(), "_clim")), .names = "{.col}_anom")) 


# PCA on anomalies

# We need to add the depth to the variables 
tab_wide <- tab_anom %>%
  select(target_date, depth, ends_with("_anom"), ends_with("_clim")) %>%
  pivot_longer(cols = ends_with("_anom"), names_to = "var", values_to = "anom") %>%
  left_join(
    tab_anom %>%
      select(target_date, depth, ends_with("_clim")) %>%
      pivot_longer(cols = ends_with("_clim"), names_to = "var_clim", values_to = "clim"),
    by = c("target_date","depth"),
    relationship = "many-to-many"
  ) %>%
  # keep the correspondance correct
  filter(gsub("_anom","",var) == gsub("_clim","",var_clim)) %>%
  mutate(var_depth = paste0(gsub("_anom","",var), "_", depth)) %>%
  select(target_date, var_depth, anom, clim)

# One column by var depth
X_anom <- tab_wide %>%
  select(target_date, var_depth, anom) %>%
  pivot_wider(names_from = var_depth, values_from = anom)


X_anom_pca <- X_anom %>% select(-target_date)

pca <- FactoMineR::PCA(X_anom_pca, axes =c(2,3))

# Estimate the number of components to reconstruct the PCA

ncp_est <- missMDA::estim_ncpPCA(X_anom_pca, scale = TRUE, method.cv=c("Kfold")) 

ncp_est$ncp # apparently it is 5

# Impute the missing values

res_anom <- missMDA::imputePCA(X_anom_pca, ncp = 5, scale = TRUE) 

imp_anom <- as.data.frame(res_anom$completeObs) 

names(imp_anom) <- c("T_anom","CHLA_anom","NO3_anom","S_anom","O_anom", "SIOH4_anom", "MES_anom") 

#  Recomposition = Imputed value = imputed anomaly + climatology (median per month)

vars <- c("T","CHLA","NO3","S","O","SIOH4","MES")

# Value_after
after <- map_dfc(vars, function(v){
  tibble(!!v := imp_anom[[paste0(v, "_anom")]] + tabc[[paste0(v, "_clim")]])
}) %>%
  bind_cols(tibble(target_date = tab$target_date), .)

# Data long: obs / recon / final / was_missing
df_long <- map_dfr(vars, function(v){
  tibble(
    target_date = tab$target_date,
    var         = v,
    obs         = tab[[v]],
    recon       = after[[v]]
  )
}) %>% mutate(was_missing = is.na(obs))

#  Ranges of each variable
ranges <- df_long %>%
  group_by(var) %>%
  summarise(
    first = suppressWarnings(min(target_date[!is.na(obs)], na.rm = TRUE)),
    last  = suppressWarnings(max(target_date[!is.na(obs)], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(has_obs = is.finite(first) & is.finite(last))

# Imputation only in the range of each variable
df_long <- df_long %>%
  left_join(ranges, by = "var") %>%
  mutate(
    in_range  = has_obs & target_date >= first & target_date <= last,
    # Last version of the series : in range or obs
    final     = ifelse(was_missing & in_range, recon, obs),
    show_red  = was_missing & in_range        # red : imputations
  )

# Plots
ggplot(df_long, aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = subset(df_long, !was_missing),
             aes(y = obs), color = "black", size = 1.1, alpha = 0.85) +
  geom_point(data = subset(df_long, show_red),
             aes(y = recon), color = "red", size = 1.3) +
  facet_wrap(~ var, scales = "free_y") +
  theme_bw() +
  labs(x = "Date", y = NULL,
       title = "Observed vs imputed")

# Final data table
final_wide <- df_long %>%
  select(target_date, var, final) %>%
  pivot_wider(names_from = var, values_from = final)


# Check this out
final_wide %>% ggplot() + geom_path(aes(x=target_date,y=T)) + theme_bw()

write_tsv(final_wide, file="data/final_wide.tsv")
