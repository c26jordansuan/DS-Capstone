library(data.table)
library(lubridate)
library(isotree)

dt <- fread("test_data_set.csv")

# Parse and clean duration/time
dt[, starttime := ymd_hms(starttime, quiet = TRUE)]
dt[, stoptime  := ymd_hms(stoptime,  quiet = TRUE)]
dt[, tripduration := as.numeric(tripduration)]
dt <- dt[is.finite(tripduration) & tripduration > 0]

# Haversine distance (meters)
haversine_m <- function(lat1, lon1, lat2, lon2) {
  R <- 6371000
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad; lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad; lon2 <- lon2 * to_rad
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

dt[, distance_m := haversine_m(
  start_station_latitude, start_station_longitude,
  end_station_latitude,   end_station_longitude
)]

dt[, speed_mps := distance_m / tripduration]
dt[, log_tripduration := log1p(tripduration)]
dt[, start_hour := hour(starttime)]

# Optional categorical -> numeric code
dt[, usertype_code := as.integer(as.factor(usertype))]

# Build feature table
features <- c("log_tripduration", "distance_m", "speed_mps", "start_hour", "usertype_code")
X <- dt[, ..features]

# Replace Inf/-Inf/NaN with NA (column-wise)
X[, (names(X)) := lapply(.SD, function(col) {
  if (is.numeric(col)) col[!is.finite(col)] <- NA_real_
  col
})]

# ---- IMPORTANT CHANGE ----
# Only REQUIRE the essential numeric trio; IMPUTE start_hour/usertype_code instead of dropping
essential <- c("log_tripduration", "distance_m", "speed_mps")
keep <- complete.cases(X[, ..essential])

dt_model <- dt[keep]
X <- X[keep]

# Impute start_hour if missing (from parse failures)
if (anyNA(X$start_hour)) {
  X$start_hour[is.na(X$start_hour)] <- as.numeric(stats::median(X$start_hour, na.rm = TRUE))
}

# Impute usertype_code if missing (rare)
if (anyNA(X$usertype_code)) {
  mode_val <- as.numeric(names(which.max(table(X$usertype_code, useNA = "no"))))
  X$usertype_code[is.na(X$usertype_code)] <- mode_val
}

cat("Rows used for training:", nrow(X), "\n")

# Guardrail: if still too few rows, stop with a clear message
if (nrow(X) < 10) {
  stop("Too few usable rows after cleaning (nrow(X) < 10). Check time parsing and missing lat/long.")
}

# Fit Isolation Forest
iso <- isolation.forest(
  as.matrix(X),                    # <- robust: give isotree a numeric matrix
  ntrees = 400,
  sample_size = min(256, nrow(X)), # must be <= nrow(X)
  seed = 42
)

dt_model[, iforest_anom_score := predict(iso, as.matrix(X), type = "score")]

contamination <- 0.02
cutoff <- quantile(dt_model$iforest_anom_score, probs = 1 - contamination, na.rm = TRUE)
dt_model[, is_anomaly := iforest_anom_score >= cutoff]

# Top anomalies
top_anoms <- dt_model[order(-iforest_anom_score)][1:25, .(
  starttime, tripduration, distance_m, speed_mps, start_hour, usertype,
  start_station_name, end_station_name, iforest_anom_score, is_anomaly
)]
print(top_anoms)

fwrite(dt_model, "bluebikes_isolation_forest_output.csv")
cat("\nSaved: bluebikes_isolation_forest_output.csv\n")

                       