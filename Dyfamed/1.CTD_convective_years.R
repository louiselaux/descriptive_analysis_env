# Library 
library(oce)
library(ncdf4)
library(castr)
library(FactoMineR)

# File directory
dir_mooring <- "/remote/complex/home/llaux/gdrive/shared/Proj_iMagine/series_JB/descriptive_analyses_env/data/dyfamed-vessel"


# Read them
library(dplyr)

read_ctd_vessel_profiles2 <- function(file) {
  lines  <- readLines(file, encoding = "latin1")
  starts <- grep("^\\*FI.*Data Type=H10", lines)
  ends   <- c(starts[-1] - 1, length(lines))
  
  lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    
    # How many parameters
    nb_line <- grep("^\\*NB PARAMETERS=", block, value = TRUE)[1]
    npar    <- as.integer(sub(".*\\*NB PARAMETERS=([0-9]+).*", "\\1", nb_line))
    
    # Name of those
    #    Remove what does not define a variable
    defs <- grep("^\\*[A-Z0-9]+ ", block, value = TRUE)
    defs <- defs[!grepl(
      "^\\*(FI|DATE|NB|GLOBAL|DC|DM|COMMENT|D[0-9]|SDN|EDMO|LOCAL_CDI|SURFACE)",
      defs
    )]
    pars <- sub("^\\*([^ ]+).*", "\\1", defs)[1:npar]
    
    # Read the table
    datalines <- grep("^[[:space:]]*-?[0-9]", block, value = TRUE)
    df <- read.table(text = datalines, header = FALSE)
    colnames(df) <- c(pars, "QC")
    
    # Date and time
    date_line <- grep("\\*DATE=", block, value = TRUE)[1]
    m <- regexec("\\*DATE=(\\d{2})(\\d{2})(\\d{4}) TIME=(\\d{2})(\\d{2})",
                 date_line)
    parts <- regmatches(date_line, m)[[1]][-1]
    dt <- as.POSIXct(sprintf("%s-%s-%s %s:%s",
                             parts[3], parts[2], parts[1],
                             parts[4], parts[5]),
                     format = "%Y-%m-%d %H:%M", tz = "UTC")
    
    df$time <- dt
    df$file <- basename(file)
    df
  })
}

# For all the files
ctd_files <- list.files(
  path       = "/remote/complex/home/llaux/gdrive/shared/Proj_iMagine/series_JB/descriptive_analyses_env/data/dyfamed-vessel",
  pattern    = "\\.ctd$",
  full.names = TRUE
)

all_profiles <- unlist(
  lapply(ctd_files, read_ctd_vessel_profiles2),
  recursive = FALSE
)

# Combine them all
ctd_all <- bind_rows(all_profiles)
head(ctd_all)

# Check something
ctd_all %>% mutate(year=year(time),
                   month=month(time))%>%filter(year==2000)%>%filter(PRES<5)

# Remove pressure problems
ctd_all <- ctd_all%>%
  filter(PRES != -999.9)

# A plot
ctd_all %>%
  filter(PRES != -999.9, TEMP != 99.999) %>%
  ggplot(aes(x = PSAL, y = PRES, group = time)) +
  geom_path(alpha = 0.7) +                           
  scale_y_reverse()  +          
  labs(
    title    = "temperatures profiles",
    x        = "Temperature (°C)",
    y        = "Pressure (dbar)"
  ) +
  theme_minimal() +
  theme(
    legend.position   = "bottom",
    legend.text        = element_text(size = 6),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

# Theta-S diagram

ctd_all %>% ggplot()+geom_path(aes(x=PSAL, y=TEMP)) # ??? what is that???

##### Despike the profiles #####

ctd_despiked <- ctd_all %>%  filter(PRES != -999.9, TEMP != 99.999) %>% 
  mutate(date  = as.Date(time),
    cast   = as.integer(factor(paste(file, time)))) %>% group_by(date, cast) %>%
  mutate(
    TEMP=castr::despike(TEMP, mult=2),
    PSAL=castr::despike(PSAL, mult=1.5)
  )


# Is the profile better?

ctd_despiked  %>%
  ggplot(aes(x = PSAL, y = PRES, group = time)) +
  geom_path(alpha = 0.7) +                           
  scale_y_reverse()
# Calculate the density

ctd_despikedd <- ctd_despiked %>%
  # Remove bad values
  filter(PRES != -999.9, PSAL != 99.999, TEMP != 99.999) %>%
  # Calculate density
  mutate(date     = as.Date(time),
         SIGMA = swRho(salinity = PSAL, 
                      temperature = TEMP, 
                      pressure    = PRES),
         sigma0   = swSigma0(salinity   = PSAL,
                             temperature= TEMP,
                             pressure   = PRES))

# Look at density profiles
ctd_despikedd %>%
  filter(date == as.Date("2000-01-30")) %>%
  ggplot(aes(x = sigma0, y = -PRES, group=time)) +
  geom_path() +
  scale_colour_datetime(name = "Heure (UTC)", date_labels = "%H:%M") +
  labs(x = expression(sigma[0]~"(kg/m"^3*")"),
       y = "Pressure (dbar)",
       title = "Profil DYFAMED du 30/01/2000") +
  theme_minimal()




# Check :
ctd_despikedd %>% distinct(cast, file, time) %>% head

ggplot(ctd_despikedd, aes(x = sigma0, y = PRES, group = cast)) +
  geom_path(size = 1, alpha = 0.8) +
  scale_y_reverse()  +
  labs(
    x     = expression(sigma[0]~"(kg/m"^3*")"),
    y     = "Pressure (dbar)",
    title = "Sigma profiles"
  ) +
  theme_minimal()

# Check it worked
ctd_despikedd%>% filter(cast=="2")%>%tail()


##### Calculate stratification index#####
strat_results <- ctd_despikedd %>%
  filter(!is.na(sigma0)) %>%
  group_by(cast,date) %>%
  summarise(
   # max depth for the profile
    bottom = max(PRES, na.rm = TRUE),
    # stratification using bottom
    strat = stratif(
      x          = sigma0,
      depth      = PRES,
      min.depths = c(0, 10),
      max.depths = c(bottom - 2, bottom),
      n.smooth   = 1,
      k          = 1,
      fun        = mean
    )
  ) %>%
  ungroup() %>%
  select(-bottom)%>%arrange(date)
head(strat_results)

# Check something
ctd_despikedd%>%group_by(date,time,cast)%>%summarize(dep=max(PRES))



##### Identify convective winters #####

##### First way : Define convective years based on mean stratification for winter months and take mean lower than a certain threshold

winter_strat2 <- strat_results %>%
  mutate(
    month = month(date),
    winter_year = if_else(month == 12, year(date) + 1, year(date))
  ) %>%
  filter(month %in% c(12,1,2,3)) %>%
  group_by(winter_year) %>%
  summarise(
    mean_strat = mean(strat, na.rm = TRUE)
  ) %>%
  ungroup()

#Check this out
print(winter_strat2)
  
strat_results %>% mutate(
  year=year(date),
  month=month(date))%>%arrange(date)
  
# Define a threshold to see what years are convective
qs <- winter_strat2 %>%
  summarise(
    q25 = quantile(mean_strat, 0.25, na.rm=TRUE),
    q50 = median(mean_strat,    na.rm=TRUE),
    q75 = quantile(mean_strat, 0.75, na.rm=TRUE)
  )

threshold <- qs$q25

winter_strat2<- winter_strat2%>%
  mutate(convective=mean_strat<=threshold)

winter_strat2%>%filter(convective==TRUE)

##### Second way : Perform a PCA and take years that contribute the most to PC1 #####
##### Try to do like Vandromme 1995 #####
profiles_env <- ctd_despikedd %>%
  filter(PRES <= 300) %>%
  group_by(cast, date) %>%
  summarise(
    T_int    = mean(TEMP,    na.rm = TRUE),  
    S_int    = mean(PSAL,    na.rm = TRUE),  
    D_int    = mean(sigma0,  na.rm = TRUE),  
    strat    = stratif(
      x          = sigma0,
      depth      = PRES,
      min.depths = c(0, 10),
      max.depths = c(max(PRES) - 2, max(PRES)),
      n.smooth   = 1,
      k          = 1,
      fun        = mean
    ),
    .groups = "drop"
  )

# Seasonal agregation
env_seas <- profiles_env %>%
  mutate(
    month       = month(date),
    winter_year = if_else(month == 12, year(date) + 1, year(date)),
    season      = case_when(
      month %in% 1:3 ~ "winter",
      month %in% 12 ~ "winter",
      month %in% 4:6 ~ "spring",
      month %in% 7:9 ~ "summer",
      TRUE           ~ "autumn"
    )
  ) %>%
  filter(season == "winter") %>%
  group_by(winter_year) %>%
  summarise(
    T_winter     = mean(T_int,  na.rm = TRUE),
    S_winter     = mean(S_int,  na.rm = TRUE),
    D_winter     = mean(D_int,  na.rm = TRUE),
    strat_winter = mean(strat,  na.rm = TRUE),
    .groups      = "drop"
  )

mat <- env_seas%>% select(-winter_year)

res.pca <- PCA(mat, scale.unit = TRUE, ncp = 5, graph = TRUE)

ind <- res.pca$ind$contrib
ind<- as.data.frame(ind)
ind$year<- env_seas$winter_year
head(ind)

ind%>% ggplot()+geom_point(aes(x=Dim.1, y=Dim.2, colour=as.factor(year)))
ind%>% arrange(Dim.1)


##### Third way : Like Dolan 2019 ####
casts_env <- ctd_despikedd %>%
  filter(PRES <= 300) %>%
  group_by(cast, date) %>%
  mutate(sigma0 = swSigma0(salinity   = PSAL,
                           temperature= TEMP,
                           pressure   = PRES)) %>%
  summarise(
    # Sigma 10 meters and 300 meters
    sigma_10   = sigma0[which.min(abs(PRES - 10))],
    sigma_300  = sigma0[which.min(abs(PRES - 300))],
    # Stratification difference 
    strat_diff = sigma_300 - sigma_10,
    # MLD
    MLD        = mld(x          = sigma0,
                     depth      = PRES,
                     ref.depths = 5:20,
                     criteria   = 0.03, # not really sure 
                     default.depth = 300,
                     n.smooth   = 1,
                     k          = 2),
    .groups = "drop"
  ) %>%
  mutate(
    nonstrat = strat_diff < 0.125    # TRUE if < threshold
  )


# Check something
casts_env%>%filter(nonstrat==TRUE)%>%mutate(year=year(date), month=month(date))%>% select(year,month)%>%
  distinct()%>% filter(month %in% c(1,2,3,12))%>%print(n=30)

casts_env%>%filter(nonstrat==TRUE)%>%mutate(year=year(date), month=month(date))%>% select(year,month, MLD)%>%
  distinct()%>%print(n=30)




# Agregate per winter
env_hiver <- casts_env %>%
  mutate(
    month       = month(date),
    winter_year = if_else(month == 12, year(date) + 1, year(date))
  ) %>%filter(month %in% c(12, 1, 2, 3))%>%
  group_by(winter_year) %>% 
  summarise(
    mean_strat_diff = mean(strat_diff, na.rm = TRUE),
    mean_MLD        = mean(MLD,         na.rm = TRUE),
    prop_nonstrat   = mean(strat_diff < 0.125, na.rm = TRUE),
    .groups = "drop"
  )

# Classification of convective winters based on stratif and MLD
env_hiver <- env_hiver %>%
  mutate(
    convectif = (mean_strat_diff < 0.125) & (mean_MLD > 80)
  )

# Result
print(env_hiver)
env_hiver%>%filter(convectif==TRUE)


##### Based on clustering#####
set.seed(42)
mat_clust <- env_hiver %>% select(mean_strat_diff, mean_MLD, prop_nonstrat) %>% scale()
km <- kmeans(mat_clust, centers = 2)
env_hiver$cluster <- km$cluster
# Convective cluster based on MLD:
conv_clust <- which.max(tapply(env_hiver$mean_MLD, km$cluster, mean))
env_hiver <- env_hiver %>%
  mutate(convectif = (cluster == conv_clust))
env_hiver%>%filter(cluster==2)
