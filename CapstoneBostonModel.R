# Install once if needed:
# install.packages(c("tidyverse", "lubridate", "igraph"))
install.packages(c("ggraph", "tidygraph", "ggplot2"))
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)

library(tidyverse)
library(lubridate)
library(igraph)

# Replace with your actual file path
trips_raw <- read_csv("test_data_set.csv")

# Check column names to confirm they match what we expect
names(trips_raw)
trips <- trips_raw %>%
  mutate(
    start_time = as.POSIXct(starttime, format = "%m/%d/%Y %I:%M:%S %p", tz = "America/New_York"),     # or started_at, etc.
    date       = as_date(start_time),
    hour_bin   = floor_date(start_time, unit = "hour")
  )

stations_start <- trips %>%
  select(
    station_id   = start_station_id,
    station_name = start_station_name,
    lat          = start_station_latitude,
    lng          = start_station_longitude
  )

stations_end <- trips %>%
  select(
    station_id   = end_station_id,
    station_name = end_station_name,
    lat          = end_station_latitude,
    lng          = end_station_longitude
  )

stations <- bind_rows(stations_start, stations_end) %>%
  filter(!is.na(station_id)) %>%
  distinct(station_id, .keep_all = TRUE)

# igraph expects a column called "name" for vertex IDs
vertices_df <- stations %>%
  mutate(
    station_id = as.character(station_id)
  ) %>%
  rename(name = station_id)

edges_df <- trips %>%
  filter(
    !is.na(start_station_id),
    !is.na(end_station_id),
    start_station_id != end_station_id  # optional: drop self-loops
  ) %>%
  count(
    from = start_station_id,
    to   = end_station_id,
    name = "weight"
  ) %>%
  mutate(
    from = as.character(from),
    to   = as.character(to)
  )

g_all <- graph_from_data_frame(
  d = edges_df,
  vertices = vertices_df,
  directed = TRUE
)

g_all

coords <- cbind(V(g_all)$lng, V(g_all)$lat)

plot(
  g_all,
  layout       = coords,
  vertex.size  = 3,
  vertex.label = NA,
  edge.arrow.size = 0.2,
  asp = 1                     # keep aspect ratio
)

## add on: new model

# --- STEP 1: Add stop_time and trip duration (in minutes) --------------------
trips = read_csv("test_data_set.csv")
summary(trips)
trips <- trips %>%
  mutate(
    # Parse start and stop times (24-hour clock, no seconds)
    start_time = as.POSIXct(starttime,
                            format = "%m/%d/%Y %H:%M",
                            tz = "America/New_York"),
    stop_time  = as.POSIXct(stoptime,
                            format = "%m/%d/%Y %H:%M",
                            tz = "America/New_York"),
    
    # Duration in minutes
    trip_duration_min = as.numeric(difftime(stop_time, start_time, units = "mins")),
    
    # Date and hour bin for later steps
    date     = as_date(start_time),
    hour_bin = floor_date(start_time, unit = "hour")
  ) %>%
  # Optional: drop weird trips
  filter(trip_duration_min > 0, trip_duration_min <= 180)

# --- STEP 2: Station × Hour key statistic: mean travel time -----------------

min_trips_per_bin <- 5  # require at least this many trips for stability

station_time_stats <- trips %>%
  mutate(
    station_id = as.character(start_station_id)  # ensure consistent type
  ) %>%
  group_by(station_id, hour_bin) %>%
  summarise(
    mean_travel_time = mean(trip_duration_min, na.rm = TRUE),
    trips_n          = n(),
    .groups = "drop"
  ) %>%
  # keep only bins with enough data
  filter(trips_n >= min_trips_per_bin)

# --- STEP 3: Build station neighborhoods (k nearest neighbors by distance) ---

# Use the original stations df; make sure station_id is character
stations_knn <- stations %>%
  mutate(station_id = as.character(station_id)) %>%
  select(station_id, lat, lng)

coords_mat <- as.matrix(stations_knn %>% select(lng, lat))  # lon/lat order

# Euclidean distance in lat/lng space is fine for city-scale analysis
dist_mat <- as.matrix(dist(coords_mat))
diag(dist_mat) <- Inf  # so we don't pick self as neighbor

k_neighbors <- 4  # fixed neighborhood cardinality like in the paper

neighbor_indices <- apply(dist_mat, 1, function(d) {
  order(d)[1:k_neighbors]  # indices of k nearest neighbors
})

station_neighbors <- tibble(
  station_id = rep(stations_knn$station_id, each = k_neighbors),
  neighbor_id = stations_knn$station_id[as.vector(neighbor_indices)]
)

head(station_neighbors)

# --- STEP 4: Join neighbors with station × hour stats and compute S(x,t) ----

# Focal station stats
st_focal <- station_time_stats %>%
  rename(
    focal_station   = station_id,
    focal_mean_time = mean_travel_time,
    focal_trips_n   = trips_n
  )

# Neighbor station stats
st_neighbor <- station_time_stats %>%
  rename(
    neighbor_station   = station_id,
    neighbor_mean_time = mean_travel_time,
    neighbor_trips_n   = trips_n
  )

# Map focal stations to their neighbors (spatial)
neighbor_pairs <- station_neighbors %>%
  rename(
    focal_station    = station_id,
    neighbor_station = neighbor_id
  )

# Combine: focal station × hour with neighbor travel times for same hour
neighbor_values <- st_focal %>%
  inner_join(neighbor_pairs, by = "focal_station") %>%
  inner_join(
    st_neighbor,
    by = c("neighbor_station", "hour_bin")
  )

# Now compute neighborhood-average travel time and S(x,t)
s_stats <- neighbor_values %>%
  group_by(focal_station, hour_bin) %>%
  summarise(
    mean_travel_time_focal = first(focal_mean_time),
    focal_trips_n          = first(focal_trips_n),
    neigh_mean_travel_time = mean(neighbor_mean_time, na.rm = TRUE),
    neigh_total_trips      = sum(neighbor_trips_n),
    .groups = "drop"
  ) %>%
  mutate(
    S = mean_travel_time_focal - neigh_mean_travel_time
  )

head(s_stats)

# --- STEP 5: Global distribution of S and outlier detection -----------------

s_dist <- s_stats %>%
  filter(is.finite(S), !is.na(S))

mu_S  <- mean(s_dist$S)
sd_S  <- sd(s_dist$S)

s_outliers <- s_stats %>%
  mutate(
    z_score   = (S - mu_S) / sd_S,
    is_outlier = abs(z_score) > 2  # 2-sigma rule; adjust if you want stricter/looser
  )

# How many outliers?
table(s_outliers$is_outlier)

# Look at the most extreme travel-time anomalies
top_anomalies <- s_outliers %>%
  arrange(desc(abs(z_score))) %>%
  slice_head(n = 50) %>%
  left_join(
    stations %>% mutate(station_id = as.character(station_id)),
    by = c("focal_station" = "station_id")
  )

head(top_anomalies)

