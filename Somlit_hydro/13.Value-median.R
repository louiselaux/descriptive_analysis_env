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

# Calculate anomalies per variable per depth
tab <- table_reg %>% pivot_wider(names_from = "name", values_from="value")
head(tab)

# Median per month per variable
clim <- tab %>% mutate(month= month(target_date))%>% group_by(month, depth) %>% summarise(across(c(T, CHLA, NO3, S, O, SIOH4, MES), ~median(., na.rm = TRUE), .names = "{.col}_clim"), .groups = "drop") 

tab<- tab %>% mutate(month=month(target_date))
tabc <- tab %>% left_join(clim, by = c("month","depth")) 
head(tabc)
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



# Identification of gaps length

# Select anomalies = value - median per variable per depth per month
tabloo <- tab_anom %>% pivot_longer(cols=c("T_anom","CHLA_anom","NO3_anom","S_anom","O_anom","SIOH4_anom","MES_anom"), values_to="value", names_to="name")




ts_all <- tabloo %>%
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
ts_all%>% select(size_gap, closest_date, name, value, target_date, depth)%>% filter(size_gap %in% c(1:2)) %>% filter(name == "T_anom")

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
interp_linear %>% mutate(year=year(target_date))%>% filter(year==2009)%>% filter(name=="T_anom")%>%ggplot() + geom_point(aes(x=closest_date, y=raw_value), shape=15, colour="blue", size=2)+
  geom_point(aes(x=target_date,y=value_final), shape=12, color="red", size=2) + theme_bw() + facet_wrap ( ~depth, scales = "free_y") 


##### Iterative PCA for bigger gaps #####

interp_wide <- interp_linear %>%
  select(target_date, closest_date, name, value_final, depth) %>%
  pivot_wider(names_from = name, values_from = value_final,
              names_prefix = "lin_")


# Convert the table to perform an interative PCA 
anom_lin_long <- interp_wide %>%
  pivot_longer(
    cols = starts_with("lin_"),
    names_to = "col",
    values_to = "anom_lin"
  ) %>%
  mutate(
    var = col %>% sub("^lin_", "", .) %>% sub("_anom$", "", .)
  ) %>%
  select(target_date, depth, var, anom_lin)

# One column per depth_var
X_anom <- anom_lin_long %>%
  mutate(var_depth = paste0(var, "_", depth)) %>%
  select(target_date, var_depth, anom_lin) %>%
  pivot_wider(names_from = var_depth, values_from = anom_lin) %>%
  arrange(target_date)

X_mat <- X_anom %>% select(-target_date)  # for PCA

##### Imputation per PCA

 ncp_est <- missMDA::estim_ncpPCA(X_mat, scale = TRUE, method.cv = "Kfold")

# How much is it 
ncp <- ncp_est$ncp

res_imp <- missMDA::imputePCA(X_mat, ncp = 12, scale = TRUE)


FactoMineR::PCA(X_mat,scale.unit = TRUE, graph = TRUE)
pca10 <- FactoMineR::PCA(X_mat %>% select(ends_with("_10")), scale.unit=TRUE)

# Completed anomalies
anom_completed_wide <- bind_cols(
  target_date = X_anom$target_date,
  as.data.frame(res_imp$completeObs)
)

# Look at it 
anom_completed_long <- anom_completed_wide %>%
  pivot_longer(-target_date, names_to = "var_depth", values_to = "anom_imp") %>%
  mutate(
    var   = sub("_\\d+$", "", var_depth),
    depth = as.numeric(sub(".*_(\\d+)$", "\\1", var_depth))
  ) %>%
  select(target_date, depth, var, anom_imp)

##### Remove dates out of range ####
first_obs_tbl <- std %>%
  filter(name %in% c("T","CHLA","NO3","S","O","SIOH4","MES")) %>%
  group_by(name, depth) %>%
  summarise(first_obs = suppressWarnings(min(date[!is.na(value)], na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(key = paste0(name, "_", depth)) %>%
  select(key, first_obs)

# Mask dates out of range not to have values badly imputed
anom_completed_wide_masked <- anom_completed_wide %>%
  pivot_longer(-target_date, names_to = "key", values_to = "anom_imp") %>%
  left_join(first_obs_tbl %>% distinct(key, .keep_all = TRUE), by = "key") %>%
  mutate(
    anom_imp = if_else(!is.na(first_obs) & target_date < first_obs,
                       NA_real_, anom_imp)
  ) %>%
  select(target_date, key, anom_imp) %>%
  pivot_wider(names_from = key, values_from = anom_imp) %>%
  arrange(target_date)

#### Test it worked ####
anom_completed_long_masked <- anom_completed_wide_masked %>%
  pivot_longer(-target_date, names_to = "var_depth", values_to = "anom_imp") %>%
  mutate(
    var   = sub("_\\d+$", "", var_depth),
    depth = as.numeric(sub(".*_(\\d+)$", "\\1", var_depth)),
    month = lubridate::month(target_date)
  ) %>%
  left_join(clim %>%
              pivot_longer(cols = ends_with("_clim"),
                           names_to = "var_clim", values_to = "clim") %>%
              mutate(var = sub("_clim$", "", var_clim)) %>%
              select(-var_clim),
            by = c("var", "month", "depth"))

# Recomposition of final series
final_with_clim <- anom_completed_long_masked %>%
  mutate(final = anom_imp + clim)


df_T_1994 <- final_with_clim %>%
  filter(var == "T", lubridate::year(target_date) %in% c(1993, 1994, 1995))

# Graph by depth
ggplot(df_T_1994, aes(x = target_date)) +
  geom_line(aes(y = final), color = "black") +
  geom_point(aes(y = anom_imp + clim), color = "red", size = 1) +
  geom_line(aes(y = clim), color = "blue", linetype = "dashed") +
  facet_wrap(~ depth, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(
    x = "Date",
    y = "Temperature (°C)",
    title = "",
    subtitle = "black = final, Rouge = imputed values, dashed blue = clim"
  )


##### Trends on this #####
glance.gls <- function(m) {
  s <- summary(m)
  
  # r.squared
  f <- predict(m)
  mss <- sum((f - mean(f))^2)
  rss <- sum(residuals(m)^2)
  rsq <- mss / (mss + rss)
  
  # residuals
  shap <- shapiro.test(m$residuals)
  a <- pacf(residuals(m, type="normalized"), plot=FALSE)
  
  tibble(
    r.squared = rsq,
    
    statistic = s$tTable[2, "t-value"],
    p.value = s$tTable[2, "p-value"],
    
    intercept = m$coefficients[1],
    slope = m$coefficients[2],
    
    shapiro.p.value = shap$p.value,
    cor.struct=class(m$modelStruct$corStruct)[1],
    acf1 = a$acf[1],
    acf2 = a$acf[2]
  )
}

# compute all trends
dstl <- anom_completed_wide_masked%>% select(target_date, CHLA_10, T_10, S_10, O_10)%>%
  pivot_longer(cols=c("CHLA_10", "T_10", "S_10", "O_10"), names_to = "name", values_to="value")
  

stats <- dstl %>%
  group_by(name) %>%
  group_modify(~{
    x <- .x %>% filter(!is.na(value))
    
    # Remove what is before the truth of the beginning 
    if (nrow(x) < 5) return(tibble())
    
    # Mann-Kendall
    mkt <- trend::mk.test(x$value)
    
    # GLS
    m <- gls(value ~ target_date, data = x)
    a <- pacf(residuals(m, type="normalized"), plot = FALSE)
    
    if (abs(a$acf[1]) > 0.2) {
      m <- gls(value ~ target_date, data = x,
               correlation = corAR1(value = round(a$acf[1], 1)))
    }
    
    bind_cols(
      glance(mkt) |> select(mk_p.value = p.value),
      glance.gls(m) |>
        select(r.squared, p.value, intercept, slope, cor.struct, acf=acf1) |>
        rename_with(~paste0("gls_", .))
    )
  }) %>%
  ungroup()

# display the result
dstl<- dstl%>%mutate(year=year(target_date))
# add date range
start_stop <- dstl %>%
  group_by(name) %>%
  summarise(start=min(year), end=max(year))

# add significance stars
# mk signif + gls non-signif may mean a non linear trend
signif_stars <- function(x) {
  case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}
stats <- stats %>%
  mutate(
    gls_acf = abs(gls_acf),
    mk_signif = signif_stars(mk_p.value),
    gls_signif = signif_stars(gls_p.value)
  ) %>%
  left_join(start_stop)



ggplot() +
  facet_wrap(name~., scales="free_y", ncol=1) +
  geom_path(aes(target_date, value), data=dstl) +
  geom_abline(aes(slope=gls_slope, intercept=gls_intercept), data=stats%>%filter( gls_signif %in% c("*", "**", "***")), colour="red", size=1.5, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ylab("Value-median") +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    strip.background = element_rect(fill = "lightblue")
  )


