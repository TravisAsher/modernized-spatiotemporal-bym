# Modernized Spatio-Temporal BYM Model
# -----------------------------------------------------------
graphics.off(); rm(list=ls())

# 1. Environment and Package Initialization
# ------------------------------------------
required_packages <- c("INLA", "sf", "spdep", "lattice", "abind", "fields", "here", "tidyverse")
lapply(required_packages, require, character.only = TRUE)
inla.setOption(scale.model.default = FALSE)


# 2. Spatial Data Ingestion & Network Building
# --------------------------------------------
ohio_shape <- st_read(here("Ohio_data", "Ohio_data", "tl_2010_39_county00.shp"))
ohio_nb <- poly2nb(ohio_shape)

# Build an explicit INLA graph to guarantee system compatibility
nb2INLA("Ohio.graph", ohio_nb)
ohio_graph <- inla.read.graph("Ohio.graph")

# Label neighborhood list using native shapefile order
county_lookup <- ohio_shape |> 
  st_drop_geometry() |> 
  select(COUNTYFP00,NAME00)
names(ohio_nb) <- county_lookup$NAME00


# Neighborhood Adjacency Network Diagnostic Plot
plot(st_geometry(ohio_shape), border = "grey")
plot(ohio_nb, st_coordinates(st_centroid(st_geometry(ohio_shape))), 
     add = TRUE, col = "red", pch = 19, cex = 0.6)



# 3. Attribute Data Preparation & Precision Alignment
# ---------------------------------------------------
voronois_ohio <- read.csv(here("Ohio_data", "Ohio_data", "OhioRespMort.csv"))

# Map county names to their native spatial row numbers to prevent index scrambling
raw_ohio_countyfp <- ohio_shape |> 
  st_drop_geometry() |> 
  select(COUNTYFP00, NAME00) |> 
  mutate(spatial_id = row_number())

voronois_ohio_final <- voronois_ohio |> 
  left_join(raw_ohio_countyfp, by = c("NAME" = "NAME00")) |>
  mutate(spatial_id_interaction = spatial_id) # Twin column for space-time interaction


# 4a. Running the Static BYM Model
# ---------------------
bym_formula <- y ~ 1 + f(spatial_id, model = "bym", graph = ohio_graph)
model_bym_ohio <- inla(bym_formula, 
                       family = "poisson", 
                       data = voronois_ohio_final, 
                       E = E, 
                       control.compute = list(dic = TRUE))
summary(model_bym_ohio)




# 4b. Running the Bernardinelli Space-Time BYM Model
# -------------------------------------------------
spacetimebym_formula <- y ~ 1 + 
  f(spatial_id, model = "bym", graph = ohio_graph, constr = TRUE) +
  f(spatial_id_interaction, year, model = "iid", constr = TRUE) + 
  year

model_param <- inla(spacetimebym_formula, 
                    family = "poisson", 
                    data = voronois_ohio_final, 
                    E = E,
                    control.predictor = list(compute = TRUE),
                    control.compute = list(dic = TRUE, cpo = TRUE))

# View overall fixed effects trends
round(model_param$summary.fixed[, 1:5], 3)



# 5a. Plot Statewide Respiratory Motality Trend
# ----------------------------------------------------
years <- seq(1, 21)

# Extract parameters explicitly by name to ensure accuracy
beta_mean <- model_param$summary.fixed["year", "mean"]
beta_lower <- model_param$summary.fixed["year", "0.025quant"]
beta_upper <- model_param$summary.fixed["year", "0.975quant"]

# Plot the main estimated slope over time
plot(years, beta_mean * years, 
     type = "l", lwd = 2, xlab = "Year Index (t)", ylab = expression(beta * t), 
     ylim = c(-0.01, 0.1), main = "Statewide Respiratory Mortality Trend")

# Add the 95% Bayesian Credible Interval dashed lines
lines(years, beta_lower * years, lty = 2, col = "red")
lines(years, beta_upper * years, lty = 2, col = "red")


# 5b. Extract Risks & Generate Geographic Visualization
# ----------------------------------------------------
# Extract the baseline structured spatial random effects (first 88 rows)
spatial_risk <- model_param$summary.random$spatial_id[1:88, ] |> 
  mutate(relative_risk = exp(mean)) |> 
  select(ID, relative_risk)

# Merge results back to map geometries using synchronized indexing
ohio_risk_map <- ohio_shape |> 
  mutate(spatial_id = row_number()) |> 
  left_join(spatial_risk, by = c("spatial_id" = "ID"))

# Render map via ggplot2
ggplot(data = ohio_risk_map) +
  geom_sf(aes(fill = relative_risk)) +
  scale_fill_viridis_c(option = "plasma", name = "Relative Risk") +
  theme_minimal() +
  labs(title = "Baseline Spatial Relative Risk of Respiratory Mortality",
       subtitle = "Ohio Counties (BYM Model)")