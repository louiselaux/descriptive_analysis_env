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

t<- t%>%mutate(CHLA=mask_extreme(CHLA, percent=c(0,0.5)),
               S=mask_extreme(S, percent=c(0.2,0)))

table_reg_stl<- t %>%pivot_longer(cols=c("CHLA","NO3","S","SIOH4","T"))



#####Second step: STL decomposition ######
table_reg_stl<- table_reg_stl%>%select(-year)

##### At the depth you choose #####
table_reg_stl <-  table_reg_stl %>%filter(depth=="1")


dstl <- table_reg_stl %>%
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
    out <- bind_cols(select(x, target_date:date_diff), out)
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
  geom_path(aes(x=target_date, y=value))




## GLS regression ----

library("broom")
library("nlme")

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
  mutate(deseason = trend+remainder) |>
  group_by(name) |>
  group_modify(.f=function(x, y) {
    # message(y)
    if (all(is.na(x$raw))) {
      return(data.frame())
    }
    # 0. fill missing values through linear interpolation
    x$deseason_filled <- castr::interpolate(x=x$target_date, y=x$deseason, xout=x$target_date)
    
    # return(x)
    # 1. Mann-Kendall trend test
    mkt <- trend::mk.test(x$deseason_filled)
    
    # 2. GLS regression
    # simple model
    m <- gls(deseason_filled ~ target_date, data=x)
    a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    # if autocorrelation is too strong
    if (abs(a$acf[1]) > 0.2) {
      # add AR1 model on residuals
      m <- gls(deseason_filled ~ target_date, data=x, cor=corAR1(round(a$acf[1], 1)))
      a <- pacf(residuals(m, type="normalized"), plot=FALSE)
      
      
    }
    
    # extract diagnostic information for both approaches
    bind_cols(
      glance(mkt) |> select(mk_p.value=p.value),
      glance.gls(m) |> select(r.squared, p.value=p.value, intercept, slope, cor.struct, acf=acf1) |> rename_with(function(n) {str_c("gls_", n)})
    )
  }) |>
  ungroup()

#save it 
#save(stats, file="output/stats_env.Rdata")

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

#save it
dstl<- dstl%>%mutate(deseason=trend+remainder)

#save(dstl, file="output/dstl.Rdata")


#####Linear interpolation#####
#x$deseason_filled <- castr::interpolate(x=x$target_date, y=x$deseason, xout=x$target_date)


#table_reg_stl_filled <- dstl %>%
#  group_by(name) %>%
#  mutate(deseason_filled = castr::interpolate(x = target_date,
#                                              y = deseason,
#                                              xout = target_date)) %>%
# ungroup()
#save(table_reg_stl_filled, file="table_reg_stl_filled.Rdata")

#Figure 1: Inspection plots of the time series#

#Inspection plot when depth is 1 meter
std%>%filter(depth=="50")%>%filter(date>"1967-01-01")%>%filter(name%in%c("CHLA","T","S","O","NO2","NO3","Sigma","COP","NOP","PO4","SIOH4"))%>%ggplot() + facet_wrap(~name, scales="free_y") +
  geom_line(aes(x=date, y=value), size=0.2)+  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  ) +
  labs(x = "Date", y = "Value")


#Inspection plot when taking the mean of each variable
std%>%group_by(date,name)%>%summarise(value=mean(value, na.rm=TRUE))%>%filter(name%in%c("CHLA","T","S","O","NO2","NO3","Sigma","COP","NOP","PO4","SIOH4"))%>%ggplot() + facet_wrap(~name, scales="free_y") +
  geom_line(aes(x=date, y=value), size=0.2)+  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  ) +
  labs(x = "Date", y = "Valeur")



#Figure 2: Long-term trend with environment time series #

#The plot
ggplot() +
  facet_wrap(name~., scales="free_y", ncol=2) +
  geom_path(aes(target_date, deseason), data=dstl %>% filter(name %in% c("NO3","SIOH4","O")), colour="grey20") +
  geom_abline(aes(slope=gls_slope, intercept=gls_intercept), data=stats%>%filter(name %in% c("NO3","SIOH4","O") & gls_signif %in% c("*", "**", "***")), colour="red", size=0.75, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ggtitle("Long term trends") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  )

