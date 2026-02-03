############################################
# Boston Blue Bikes: iForest + proxy labels
# Columns:
# tripduration,starttime,stoptime,start_station_id,start_station_name,
# start_station_latitude,start_station_longitude,end_station_id,end_station_name,
# end_station_latitude,end_station_longitude,bikeid,usertype,birth year,gender
#
# Time format example: "1/1/2015 0:21"  (m/d/Y H:M)
############################################

# ---- 0) Packages ----
pkgs <- c("tidyverse", "lubridate", "geosphere", "isotree")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(tidyverse)
library(lubridate)
library(geosphere)  # distHaversine
library(isotree)    # isolation.forest

# ---- 1) Read data ----
path <- "test_data_set.csv"  # <-- change to your file path
df <- readr::read_csv(path, show_col_types = FALSE)

# ---- 2) Parse timestamps + clean basics ----
# Your times are m/d/Y H:M (no seconds). We'll parse that explicitly.
df2 <- df %>%
  mutate(
    tripduration = suppressWarnings(as.numeric(tripduration)),
    starttime = lubridate::parse_date_time(starttime, orders = c("mdY HM")),
    stoptime  = lubridate::parse_date_time(stoptime,  orders = c("mdY HM")),
    birth_year = suppressWarnings(as.integer(`birth year`)),
    gender = suppressWarnings(as.integer(gender))
  ) %>%
  # Hard filters for impossible/broken rows (not "anomalies"—just garbage)
  filter(
    !is.na(starttime),
    !is.na(stoptime),
    !is.na(tripduration),
    tripduration > 0,
    tripduration <= 24*60*60,  # <= 24 hours; adjust if desired
    !is.na(start_station_latitude), !is.na(start_station_longitude),
    !is.na(end_station_latitude), !is.na(end_station_longitude)
  ) %>%
  mutate(
    duration_min = tripduration / 60,
    start_hour   = hour(starttime),
    start_wday   = wday(starttime, week_start = 1),  # 1=Mon ... 7=Sun
    is_weekend   = as.integer(start_wday %in% c(6, 7)),
    start_month  = month(starttime),
    
    # Distance in meters (lon/lat order!)
    dist_m = geosphere::distHaversine(
      cbind(start_station_longitude, start_station_latitude),
      cbind(end_station_longitude,   end_station_latitude)
    ),
    dist_km = dist_m / 1000,
    speed_km_min = dist_km / pmax(duration_min, 1e-6),
    
    route_key = paste(start_station_id, end_station_id, sep = "->")
  ) %>%
  # optional: drop impossible coordinates/distances
  filter(is.finite(dist_m), dist_m >= 0, dist_m <= 50000)

# ---- 3) Build proxy "expected duration" baseline (route + hour-of-day) ----
# Robust stats: median + MAD
mad_robust <- function(x) mad(x, constant = 1.4826, na.rm = TRUE)

route_hour_stats <- df2 %>%
  group_by(route_key, start_hour) %>%
  summarise(
    n_route_hour = n(),
    med_dur = median(duration_min, na.rm = TRUE),
    mad_dur = mad_robust(duration_min),
    .groups = "drop"
  )

route_stats <- df2 %>%
  group_by(route_key) %>%
  summarise(
    n_route = n(),
    med_dur_route = median(duration_min, na.rm = TRUE),
    mad_dur_route = mad_robust(duration_min),
    .groups = "drop"
  )

global_med <- median(df2$duration_min, na.rm = TRUE)
global_mad <- mad_robust(df2$duration_min)

df3 <- df2 %>%
  left_join(route_hour_stats, by = c("route_key", "start_hour")) %>%
  left_join(route_stats,      by = "route_key") %>%
  mutate(
    # Use route+hour if we have enough samples, else route-only, else global
    exp_dur = case_when(
      !is.na(med_dur) & n_route_hour >= 20 ~ med_dur,
      !is.na(med_dur_route) & n_route >= 20 ~ med_dur_route,
      TRUE ~ global_med
    ),
    exp_mad = case_when(
      !is.na(mad_dur) & n_route_hour >= 20 & mad_dur > 0 ~ mad_dur,
      !is.na(mad_dur_route) & n_route >= 20 & mad_dur_route > 0 ~ mad_dur_route,
      TRUE ~ global_mad
    ),
    
    # Proxy deviation metrics
    z_robust = (duration_min - exp_dur) / pmax(exp_mad, 1e-6),
    ratio    = duration_min / pmax(exp_dur, 1e-6)
  )

# ---- 4) Create proxy anomaly labels (weak labels) ----
# Keep it rare (e.g., 1%). This is just for “relative labeling”.
p_anom <- 0.01  # 1% (try 0.005 to 0.02)

z_cut <- quantile(df3$z_robust, probs = 1 - p_anom, na.rm = TRUE)

df3 <- df3 %>%
  mutate(
    proxy_anom = as.integer(z_robust >= z_cut & duration_min >= 5),
    proxy_anom_hard = as.integer(z_robust >= 8 | ratio >= 4)
  )

# ---- 5) Prepare features for iForest (numeric + one-hot) ----
current_year <- year(Sys.Date())

df4 <- df3 %>%
  mutate(
    age = if_else(!is.na(birth_year) & birth_year >= 1900 & birth_year <= current_year,
                  current_year - birth_year,
                  NA_integer_),
    age = if_else(!is.na(age) & age >= 10 & age <= 90, age, NA_integer_),
    
    usertype = as.factor(usertype),
    gender_f = as.factor(gender)  # keep original numeric gender as factor for one-hot
  )

# numeric features
num_feats <- df4 %>%
  transmute(
    duration_min,
    dist_km,
    speed_km_min,
    start_hour,
    start_wday,
    is_weekend,
    start_month,
    age = if_else(is.na(age), median(age, na.rm = TRUE), as.numeric(age))
  )

# categorical one-hot
cat_mm <- model.matrix(~ usertype + gender_f - 1, data = df4)

X <- cbind(as.matrix(num_feats), cat_mm)

# drop rows with any non-finite (rare)
keep <- apply(X, 1, function(r) all(is.finite(r)))
X <- X[keep, , drop = FALSE]
df_model <- df4[keep, ]

# ---- 6) Train Isolation Forest (unsupervised) ----
set.seed(42)

iso <- isolation.forest(
  X,
  ntrees = 400,
  sample_size = min(256, nrow(X)),
  ndim = 1,
  nthreads = max(1, parallel::detectCores() - 1)
)

scores <- predict(iso, X, type = "score")

# Choose iForest anomalies as top p% by score (same rate as proxy labels by default)
score_cut <- quantile(scores, probs = 1 - p_anom, na.rm = TRUE)

df_out <- df_model %>%
  mutate(
    iforest_score = as.numeric(scores),
    iforest_anom  = as.integer(iforest_score >= score_cut)
  )

# ---- 7) Quick diagnostics ----
cat("\nOverlap between proxy labels and iForest flags:\n")
print(table(proxy = df_out$proxy_anom, iforest = df_out$iforest_anom))

cat("\nTop 25 iForest anomalies:\n")
top_anoms <- df_out %>%
  arrange(desc(iforest_score)) %>%
  select(
    iforest_score, iforest_anom,
    proxy_anom, proxy_anom_hard,
    tripduration, duration_min, dist_km, speed_km_min, z_robust, ratio,
    starttime, stoptime,
    start_station_name, end_station_name,
    usertype, `birth year`, gender
  ) %>%
  slice_head(n = 25)

print(top_anoms)

# ---- 8) Save outputs ----
readr::write_csv(df_out, "boston_bluebikes_iforest_scored.csv")
df_out %>% filter(iforest_anom == 1) %>%
  readr::write_csv("boston_bluebikes_iforest_anomalies.csv")

# ---- 9) Simple plots (optional) ----
ggplot(df_out, aes(iforest_score)) +
  geom_histogram(bins = 50) +
  labs(
    title = "Isolation Forest anomaly score distribution",
    x = "iForest score (higher = more anomalous)",
    y = "Count"
  )

ggplot(df_out, aes(z_robust, iforest_score)) +
  geom_point(alpha = 0.2) +
  labs(
    title = "Proxy deviation (robust z) vs iForest score",
    x = "Robust deviation vs route/hour expected duration",
    y = "iForest score"
  )

# Determine which streets these anomalies occur: 
anomaly_routes <- read.csv("boston_bluebikes_iforest_scored.csv")

# How many different street combinations are in the dataset?
num_combinations_anomaly <- anomaly_routes %>%
  filter(iforest_anom == 1) %>%
  distinct(start_station_name, end_station_name) %>%
  nrow()
num_combinations_anomaly

# What are the exact pairs?
start_end_anomaly <- anomaly_routes %>%
  filter(iforest_anom == 1) %>%
  distinct(start_station_name, end_station_name)

head(start_end_anomaly)

# What are the most common routes?
top_anomaly <- anomaly_routes %>%
  filter(iforest_anom == 1) %>%
  count(start_station_name, end_station_name, sort = TRUE)

head(top_anomaly, 10)

