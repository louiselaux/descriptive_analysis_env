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

nrow(tabc) == nrow(tab)
anti_join(tab, clim, by = c("month","depth")) %>% nrow()

# Anomalies= value - clim per month

tab_anom <- tabc %>%
  mutate( T_anom= T-T_clim,
   CHLA_anom= CHLA- CHLA_clim,
   NO3_anom= NO3-NO3_clim,
   S_anom= S-S_clim,
   O_anom= O-O_clim,
   SIOH4_anom= SIOH4-SIOH4_clim,
   MES_anom= MES-MES_clim)

   


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

# PCA at 10 meters depth
pca10 <- FactoMineR::PCA(X_anom_pca %>% select(ends_with("_10")), scale.unit=TRUE)

# PCA at all depths 
pca <- FactoMineR::PCA(X_anom_pca, axes =c(1,2))

# Estimate the number of components to reconstruct the PCA

#ncp_est <- missMDA::estim_ncpPCA(X_anom_pca, scale = TRUE, method.cv=c("Kfold")) 

ncp_est$ncp # apparently it is 5 

# How many explained variance
pca$eig

# Impute the missing values

res_anom <- missMDA::imputePCA(X_anom_pca, ncp = 12, scale = TRUE) 


# Add the date column
imp_anom <- res_anom$completeObs
imp_anom <- dplyr::bind_cols(target_date = X_anom$target_date,
                             as.data.frame(imp_anom))


#####
# Change the format
clim_wide <- tab_wide %>%
  select(target_date, var_depth, clim) %>%
  tidyr::pivot_wider(names_from = var_depth, values_from = clim) %>%
  arrange(target_date)

# Recomposition : final value = anomaly imputed + climatology
common_cols <- setdiff(intersect(names(imp_anom), names(clim_wide)), "target_date")

final_wide <- imp_anom
for (col in common_cols) {
  final_wide[[col]] <- imp_anom[[col]] + clim_wide[[col]]
}
####

# Get everything back to a normal data frame
vars <- c("T","CHLA","NO3","S","O","SIOH4","MES")

obs_long <- tablo_merged %>%
  select(target_date, depth, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "var", values_to = "obs")


df_long <- final_wide %>%
  pivot_longer(cols = -target_date, names_to = "var_depth", values_to = "recon") %>%
  tidyr::separate(var_depth, into = c("var","depth"), sep = "_", convert = TRUE) %>%
  left_join(obs_long, by = c("target_date","depth","var")) %>%
  group_by(var, depth) %>%
  mutate(
    first   = suppressWarnings(min(target_date[!is.na(obs)], na.rm = TRUE)),
    last    = suppressWarnings(max(target_date[!is.na(obs)], na.rm = TRUE)),
    has_obs = is.finite(first) & is.finite(last),
    in_range = has_obs & target_date >= first & target_date <= last,
    was_missing = is.na(obs),
    final = ifelse(was_missing & in_range, recon, obs),
    show_red = was_missing & in_range
  ) %>%
  ungroup()

# Clipping the reconstructed values
bounds_final <- obs_long %>%
  filter(!is.na(obs)) %>%
  group_by(var, depth) %>%
  summarise(q01f = quantile(obs, 0.01, na.rm=TRUE),
            q99f = quantile(obs, 0.99, na.rm=TRUE),
            .groups="drop")

df_long <- df_long %>%
  left_join(bounds_final, by=c("var","depth")) %>%
  mutate(
    final = ifelse(
      was_missing & in_range,                      # only imputs
      pmin(pmax(final, q01f), q99f),               # clip
      final                                        
    )
  )

# See if worked
ggplot(df_long%>%filter(var=="T"), aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = subset(df_long, !was_missing),
             aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = subset(df_long, show_red),
             aes(y = recon), color = "red", size = 1.2) +
  facet_grid(var ~ depth, scales = "free_y") +
  theme_bw() +
  labs(x = "Date", y = NULL, title = "Observed vs imputed (by depth)")

ggplot(df_long %>% filter(var=="T"), aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = subset(df_long, !was_missing),
             aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = subset(df_long, show_red),
             aes(y = final), color = "red", size = 1.2) +   
  facet_grid(var ~ depth, scales = "free_y") +
  theme_bw()


ggplot(df_long %>% filter(var == "T") %>% mutate(year=year(target_date))%>%filter(year=="1995"), aes(x = target_date)) +geom_line(aes(y = final)) +
  geom_point(data = subset(df_long, var == "T" & !was_missing),
             aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = subset(df_long, var == "T" & show_red),
             aes(y = final), color = "red", size = 1.2) +
  facet_wrap(~ depth, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(x = "Date", y = "Temperature (°C)",
       title = "Observed vs imputed Temperature by depth")

# Focus year 1995
df_T_1995 <- df_long %>%
  filter(var == "T", year(target_date) %in% c(1993,1994,1995))

ggplot(df_T_1995, aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = df_T_1995 %>% filter(!was_missing),
             aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = df_T_1995 %>% filter(show_red),
             aes(y = final), color = "red", size = 1.2) +
  facet_wrap(~ depth, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(x = "Date", y = "Temperature (°C)",
       title = "Observed vs imputed Temperature by depth — 1995") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

# Final data frame
final_panel <- df_long %>%
  select(target_date, depth, var, final) %>%
  pivot_wider(names_from = var, values_from = final) %>%
  arrange(depth, target_date) 

# Save it 
readr::write_tsv(final_panel, "data/final_panel_by_depth.tsv")

write_tsv(final_wide, file="data/final_wide.tsv")


##### Comparison with raw data and not anomaly ####
vars <- c("T","CHLA","NO3","S","O","SIOH4","MES")

obs_long <- tablo_merged %>%
  select(target_date, depth, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to = "var", values_to = "obs")

X_raw <- tablo_merged %>%
  select(target_date, depth, all_of(vars)) %>%
  pivot_longer(cols = all_of(vars), names_to="var", values_to="val") %>%
  mutate(var_depth = paste0(var, "_", depth)) %>%
  select(target_date, var_depth, val) %>%
  pivot_wider(names_from = var_depth, values_from = val) %>%
  arrange(target_date)

X_raw_pca <- X_raw %>% select(-target_date)

res_raw <- missMDA::imputePCA(X_raw_pca, ncp = 8, scale = TRUE)

imp_raw <- dplyr::bind_cols(
  target_date = X_raw$target_date,
  as.data.frame(res_raw$completeObs)
)

df_long_raw <- imp_raw %>%
  pivot_longer(cols = -target_date, names_to = "var_depth", values_to = "recon") %>%
  tidyr::separate(var_depth, into = c("var","depth"), sep = "_", convert = TRUE) %>%
  left_join(obs_long, by = c("target_date","depth","var")) %>%
  group_by(var, depth) %>%
  mutate(
    first      = suppressWarnings(min(target_date[!is.na(obs)], na.rm = TRUE)),
    last       = suppressWarnings(max(target_date[!is.na(obs)], na.rm = TRUE)),
    has_obs    = is.finite(first) & is.finite(last),
    in_range   = has_obs & target_date >= first & target_date <= last,
    was_missing = is.na(obs),
    final      = ifelse(was_missing & in_range, recon, obs),  
    show_red   = was_missing & in_range
  ) %>%
  ungroup()

ggplot(df_long_raw %>% filter(var=="T"), aes(target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = subset(df_long_raw, !was_missing),
             aes(y = obs), size=1, alpha=.85) +
  geom_point(data = subset(df_long_raw, show_red),
             aes(y = final), color="red", size=1.2) +
  facet_grid(var ~ depth, scales = "free_y") +
  theme_bw()

# Focus on year 1995
df_T_1995 <- df_long_raw %>%
  filter(var == "T", year(target_date) %in% c(1994))

ggplot(df_T_1995, aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = df_T_1995 %>% filter(!was_missing),
             aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = df_T_1995 %>% filter(show_red),
             aes(y = final), color = "red", size = 1.2) +
  facet_wrap(~ depth, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(x = "Date", y = "Temperature (°C)",
       title = "Observed vs imputed Temperature by depth — 1995") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")



