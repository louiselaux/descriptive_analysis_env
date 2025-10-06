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

final_wide <- read_tsv("data/final_wide.tsv")
head(final_wide)
final_wide <- as.data.frame(final_wide)


# Log the chla variable, log10+1 to better deal with 0 that will be -inf otherwise
final_wide <- final_wide %>%mutate(
  CHLA_10= log10(CHLA_10+1)
)


final_wide <- final_wide %>% select(target_date,T_10,CHLA_10,NO3_10,S_10,O_10,SIOH4_10,MES_10)
table_reg_stl <- final_wide %>% pivot_longer(cols=c("T_10","CHLA_10","NO3_10","S_10","O_10","SIOH4_10","MES_10"), names_to="name", values_to="value")
head(table_reg_stl)

##### STL decomposition #####

dstl <- table_reg_stl %>% filter(name %in% c("T_10","CHLA_10","NO3_10","S_10","O_10","SIOH4_10"))%>%
  group_by(name) %>%
  group_modify(.f=function(x,y) {
    # message(y)
    # if all is missing, do not do anything
    if ( all(is.na(x$value)) ) {
      out <- NULL
      # else perform stl
    } else {
      dec <- stlplus(x$value, x$target_date, n.p=26, s.window="periodic", t.window=10*26)
      out <- dec$data |> select(raw, seasonal, trend, remainder)
    }
    out <- bind_cols(select(x, target_date), out)
  }) |>
  ungroup() |>
  
  # cut the part before the variable becomes available for the first time
  group_by(name) |>
  group_modify(.f=function(x, y) {
    if (all(is.na(x$raw))) {
      out <- NULL
    } else {
      x <- arrange(x, target_date)
      start_idx <- min(which(!is.na(x$raw)))
      out <- x[start_idx:nrow(x),]
    }
    return(out)
  }) |>
  ungroup()

# Plot the result
dstl %>%
  pivot_longer(raw:remainder, names_to="component") |>
  mutate(component=factor(component, levels=c("raw", "trend", "seasonal", "remainder"))) |>
  ggplot() + facet_wrap(~interaction(name, component), scale="free", nrow=4) +
  geom_path(aes(x=target_date, y=value)) + theme_bw()


head(dstl)

write_tsv(dstl, file="data/dstl_interpolated_before.tsv")
##### Gls regression #####
#Second step : STL decomposition
## GLS regression ----

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
stats <- dstl %>%
  mutate(deseason = trend + remainder) %>%
  group_by(name) %>%
  group_modify(~{
    x <- .x %>% filter(!is.na(deseason))
    
    # Remove what is before the truth of the beginning 
    if (nrow(x) < 5) return(tibble())
    
    # Mann-Kendall
    mkt <- trend::mk.test(x$deseason)
    
    # GLS
    m <- gls(deseason ~ target_date, data = x)
    a <- pacf(residuals(m, type="normalized"), plot = FALSE)
    
    if (abs(a$acf[1]) > 0.2) {
      m <- gls(deseason ~ target_date, data = x,
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


# Save it 

write_tsv(stats, file="data/stats_interpolated_before.tsv")

#save it
dstl<- dstl%>%mutate(deseason=trend+remainder)


ggplot() +
  facet_wrap(name~., scales="free_y", ncol=1) +
  geom_path(aes(target_date, deseason), data=dstl%>%filter(name %in% c("CHLA_10","T_10","S_10","O_10","NO3_10"))) +
  geom_abline(aes(slope=gls_slope, intercept=gls_intercept), data=stats%>%filter( gls_signif %in% c("*", "**", "***")), colour="red", size=1.5, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ylab("Deseasonalized component") +
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




