# Fish landings

# Load data

grid <- read.csv("data/fisheries_grid.csv", header = T)
bio <- read.csv("data/biomass_by_ground.csv", header = T)

library(dplyr)

bio <- bio %>% 
  rename(gridID = Row.Labels)

grid <- grid %>% 
  rename(gridID = Fishery_gr)

bio_grid <- left_join(grid, bio, by = "gridID")

bio_grid <- bio_grid %>% 
  mutate_all(~replace(., is.na(.), 0))

write.csv(bio_grid, "data/biomass_on_grid.csv")
