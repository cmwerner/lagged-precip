### Data reading and cleaning script for Stipa - Avena lagged precip study


library(here)
library(tidyverse)


datpath <- "data/"

# block and plot key
trt.key <- read.csv(paste0(datpath, "plot-layout/block-key.csv")) %>%
  select(block:treatment) %>%
  mutate(treatment = ifelse(treatment == "dry", "drought", "ambient"))

# Stipa fecundity 2023
stipa_fec_23 <- read.csv(paste0(datpath,"lagged-precip_stipa-fecundity_2023-04.csv")) %>%
  rename(stipa_type = Stipa..type,
         stipa_diam = Stipa.diam..cm.,
         weed_gram = X..weed.gram, 
         weed_forb = X..weed.forb) %>%
  rename_with(.,~ tolower(gsub( ".", "_", .x, fixed = TRUE))) %>%
  rename_with(.,~ gsub( "stipa", "sp", .x, fixed = TRUE)) %>%
  mutate(year = 2023) 

# Stipa fecundity 2024
stipa_fec_24 <- read.csv(paste0(datpath,"lagged-precip_stipa-fecundity_2024-04.csv")) %>%
  rename_with(.,~ tolower(gsub( ".", "_", .x, fixed = TRUE))) %>%
  mutate(year = 2024,
         sp_seeds_per_tiller = as.numeric(sp_seeds_per_tiller)) 

# Join 2023 and 2024 stipa fecundity
stipa_fec <- full_join(stipa_fec_23, stipa_fec_24, 
                       by = join_by(block, plot, bg, sp_type,
                                    sp_diam == sp_diam_1, sp_tillers, 
                                    seeds_tiller == sp_seeds_per_tiller, year,
                                    notes)) %>%
  left_join(trt.key, by = "block") %>%
  mutate(bg = case_when(  
    bg == "S_P low" | bg == "S_Plow" ~ "sp_low",
    bg == "S_P high" | bg == "S_Phigh" ~ "sp_high",
    bg == "S_s low" | bg == "S_slow" ~ "ss_low",
    bg == "S_s high" | bg == "S_shigh" ~ "ss_high",
    bg == "A low" | bg == "A_low" ~ "a_low",
    bg == "A high" | bg == "A_high" ~ "a_high",
    bg == "control" ~ "control",
    .default = "NA"
  )) # making plot background labels consistent



# Avena density and fecundity 
# NOTE: background plot labels are not consistent across years here -- need to update so they are

# Avena 2023
# data was entered in a weird wide format, took a little reshuffling to get it into a workable long format
avena_23 <-  read.csv(paste0(datpath,"lagged-precip_avena_2023-04.csv")) %>%
  rename_with(.,~ tolower(gsub( ".", "_", .x, fixed = TRUE))) %>%
  mutate(across(sub1_plants:full_im_tillers, as.numeric),
         year = 2023) %>%
  pivot_longer(cols = sub1_plants:full_im_tillers) %>%
  separate_wider_delim(cols = name, delim = "_", names = c("location", "info"), too_many = "merge") %>%
  pivot_wider(names_from = "info", values_from = "value") %>%
  pivot_longer(cols = seeds_1:im_5) %>%
  separate_wider_delim(cols = name, delim = "_", names = c("tiller_type", "rep")) %>%
  pivot_wider(names_from = "tiller_type", values_from ="value") %>%
  rename(im_seeds = im)

# Avena 2024
avena_24 <-  read.csv(paste0(datpath,"lagged-precip_avena_2024-05.csv")) %>%
  rename_with(.,~ tolower(gsub( ".", "_", .x, fixed = TRUE))) %>%
  mutate(year = 2024, 
         location = case_when(
           location == "" ~ "full",  # missing location info in one row
           .default = location))

avena <- full_join(avena_23, avena_24, by = join_by(block, plot, bg, location, 
                                                    plants, seeds, im_tillers, year)) %>%
  left_join(trt.key, by = "block") %>%
  mutate(bg = case_when(  
    bg == "S_P low" | bg == "S_Plow" ~ "sp_low",
    bg == "S_P high" | bg == "S_Phigh" ~ "sp_high",
    bg == "S_s low" | bg == "S_slow" ~ "ss_low",
    bg == "S_s high" | bg == "S_shigh" ~ "ss_high",
    bg == "A low" | bg == "A_low" ~ "a_low",
    bg == "A high" | bg == "A_high" ~ "a_high",
    bg == "control" ~ "control",
    .default = "NA"
  )) # making plot background labels consistent

# avena density calculation keeping subplots separate
# for comparing seed-added to non-seed-added side
avena <- avena %>%
  mutate(
    density = case_when(
      location == "full" ~ plants / 0.42, # 0.5 x 1m full plot, minus the two .2 x .2 subplots
      .default = plants * 25 # 0.2 x 0.2 m subplot
    )
  )


# plot-level information for overall avena densities per m2
avena_plot <- avena %>%
  group_by(block, plot, bg, year, treatment) %>%
  summarise(
    avena_count = sum(plants, na.rm = TRUE),
    area = case_when(
      is.na(sum(plants)) ~ 0.08,  # only subplots counted
      .default = 0.5  # whole plot counted
      ),
    avena_density = avena_count / area
  )

## survival/recruitment

# 2023, need to double-check what the sections of 1, 2, and 3 are. 2/3 of the plot was overseeded
# in the second year, so this is separating out first added vs newly added seed
# but need to check which side 1 vs 3 are

# this is also the dataset that has weed info for this time period
ss_23 <- read.csv(paste0(datpath,"lagged-precip_survival_2023-02.csv")) %>%
  mutate(ss_old_tot = ss.old.1 + ss.old.2 + ss.old.3, 
         ss_new_tot = ss.new.1 + ss.new.2 + ss.new.3) %>% 
  left_join(trt.key, by = "block") %>%
  mutate(bg = case_when(  
    bg == "S_P low" | bg == "S_Plow" ~ "sp_low",
    bg == "S_P high" | bg == "S_Phigh" ~ "sp_high",
    bg == "S_s low" | bg == "S_slow" ~ "ss_low",
    bg == "S_s high" | bg == "S_shigh" ~ "ss_high",
    bg == "A low" | bg == "A_low" ~ "a_low",
    bg == "A high" | bg == "A_high" ~ "a_high",
    bg == "control" ~ "control",
    .default = "NA"
  ))

# 2024 data needs cleaning of block/background names
# also not necessarily directly connected to 2023 data, need to think about how best to combine

# info for plot themes and colors
theme_cw <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      panel.background = element_blank(), 
      plot.background = element_blank(), 
      axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
      panel.grid.minor = element_blank(), 
      panel.grid.major.x = element_blank(),
      legend.background = element_blank(), 
      legend.key = element_blank(),
      strip.background = element_blank(), 
      # axis.text=element_text(size=12),
      # strip.text=element_text(size=12),
      complete = TRUE
    )
}

drought.colors <- c("#0077CC", "#CD622E")

