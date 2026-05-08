library(myClim)
library(readxl)
library(lubridate)
library(here)
library(tidyverse)

# NOTE: sensor 181 was put in later than the rest

datpath <- here("data/TMS-logger-data/")

soil_data <- mc_read_files(datpath, dataformat_name = "TOMST", 
                           recursive = F, silent = T) %>%
  mc_prep_clean() 

tms_id <- read.csv("data/lagged-precip_TMS-ids.csv") %>%
  mutate(serial_number = as.character(serial_number)) %>%
  rename(bg = X)

trt.key <- read.csv("data/plot-layout/block-key.csv") %>%
  select(block:treatment) %>%
  mutate(treatment = ifelse(treatment == "dry", "drought", "ambient"))


# can use mc_prep_solar_tz() to calculate solar time rather than UTC

# View(mc_info(soil_data))

# crop the time series
start_date <- as.POSIXct("2024-01-27", tz = "UTC")
end_date <- as.POSIXct("2024-05-04", tz = "UTC")
soil_data <- mc_prep_crop(soil_data, start_date, end_date)


# aggregate to daily average
soil_daily <- soil_data %>%
  mc_agg(fun = "mean", period = "day", min_coverage = 0.95)


# Plotting everything
# mc_plot_line(soil_data)

# long format daily summary
tms_long <- soil_daily %>% 
  mc_reshape_long() %>%
  rename(sensor_name_combined = sensor_name) %>%
  separate_wider_delim(sensor_name_combined, "_", names = c(NA, "sensor_name", "metric")) %>%
  left_join(tms_id, by = join_by(locality_id == serial_number)) %>%
  left_join(trt.key, by = join_by(block))

tms_moisture <- tms_long %>% 
  filter(sensor_name == "moist") 


tms_control <- tms_moisture %>% filter(bg == "control")

drought.colors <- c("#0077CC", "#CD622E")

tms_control %>% 
  ggplot(aes(x = datetime, y = value, 
             group = locality_id, color = treatment)) +
  geom_line() +
  scale_color_manual(values = drought.colors)
