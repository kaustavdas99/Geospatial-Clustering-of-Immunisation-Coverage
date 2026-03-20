# -----------------------------
# Install / load packages
# -----------------------------
library(sf)
library(dplyr)
library(readxl)
library(spdep)
library(tmap)
library(RColorBrewer)
library(terra)

# -----------------------------
# Read shapefile and filter West Bengal
# -----------------------------
wb_shp <- st_read("District_shape_West_Bengal.shp")

#wb_shp <- st_read("District_shape_West_Bengal.shp")%>%
  #filter(ST_NM == "West Bengal")

# Check districts
print(unique(wb_shp$NAME))

# -----------------------------
# Read immunisation index Excel
# -----------------------------
immun <- read_excel("District_Cumulative_Immunisation_Risk.xlsx")

# Rename column for consistency
immun <- immun %>% rename(Average_Index = `Average_Index`)

# -----------------------------
# District name mapping and category aggregation
# -----------------------------
district_map <- c(
  "Alipurduar" = "Alipurduar",
  "Bankura" = "Bankura",
  "Birbhum" = "Birbhum",
  "Coochbehar" = "Koch Bihar",
  "Dakshin Dinajpur" = "Dakshin Dinajpur",
  "Darjeeling" = "Darjiling",
  "Hooghly" = "Hugli",
  "Howrah" = "Haora",
  "Jalpaiguri" = "Jalpaiguri",
  "Jhargram" = "Jhargram",
  "Kalimpong" = "Kalimpong",
  "Kolkata" = "Kolkata",
  "Malda" = "Maldah",
  "Murshidabad" = "Murshidabad",
  "Nadia" = "Nadia",
  "North 24 Pargana" = "North 24 Parganas",
  "Paschim Bardhaman" = "Paschim Barddhaman",
  "Paschim Medinipur" = "Pashchim Medinipur",
  "Purba Bardhaman" = "Purba Barddhaman",
  "Purba Medinipur" = "Purba Medinipur",
  "Purulia" = "Puruliya",
  "South 24 Pargana" = "South 24 Parganas",
  "Uttar Dinajpur" = "Uttar Dinajpur"
)

# Apply mapping
immun <- immun %>%
  mutate(DISTRICT_2011 = district_map[District])

# -----------------------------
# Aggregate duplicates
# -----------------------------
immun_agg <- immun %>%
  group_by(DISTRICT_2011) %>%
  summarise(Immunisation_Index = mean(Average_Index, na.rm = TRUE))

# -----------------------------
# Merge aggregated index into shapefile
# -----------------------------
wb_data <- wb_shp %>%
  left_join(immun_agg, by = c("NAME" = "DISTRICT_2011"))

# Convert to numeric
wb_data$Immunisation_Index <- as.numeric(wb_data$Immunisation_Index)

# LISA and Gi* columns
wb_data$LISA_cluster <- NA
wb_data$GiZ <- NA

# Merge
wb_data %>% st_drop_geometry() %>% select(NAME, Immunisation_Index)

# -----------------------------
# Create neighbors
# -----------------------------
nb <- poly2nb(wb_data)
lw <- nb2listw(nb, style = "W")

# -----------------------------
# Global Moran's I
# -----------------------------
index_values <- wb_data$Immunisation_Index
global_moran <- moran.test(index_values, lw, na.action = na.omit)
print(global_moran)

# -----------------------------
# Local Moran's I (LISA)
# -----------------------------
non_na_idx <- !is.na(index_values)
wb_data_sub <- wb_data[non_na_idx, ]

nb_sub <- poly2nb(wb_data_sub)
lw_sub <- nb2listw(nb_sub, style = "W")

local_moran <- localmoran(wb_data_sub$Immunisation_Index, lw_sub)

# Results
wb_data_sub$Ii <- local_moran[,1]
wb_data_sub$Z_Ii <- local_moran[,4]
wb_data_sub$Ii_p <- local_moran[,5]

# Classify LISA clusters
wb_data_sub$LISA_cluster <- factor(
  ifelse(local_moran[,4] > 1.96 & wb_data_sub$Immunisation_Index > mean(wb_data_sub$Immunisation_Index), "High-High",
         ifelse(local_moran[,4] > 1.96 & wb_data_sub$Immunisation_Index < mean(wb_data_sub$Immunisation_Index), "Low-Low",
                ifelse(local_moran[,4] < -1.96 & wb_data_sub$Immunisation_Index > mean(wb_data_sub$Immunisation_Index), "High-Low",
                       "Low-High"))),
  levels = c("High-High","Low-Low","High-Low","Low-High")
)

# Fill shapefile
wb_data$LISA_cluster[non_na_idx] <- wb_data_sub$LISA_cluster

# -----------------------------
# Getis-Ord Gi* Hotspot/Coldspot
# -----------------------------
gi_star <- localG(wb_data_sub$Immunisation_Index, lw_sub)
wb_data_sub$GiZ <- as.numeric(gi_star)

wb_data$GiZ[non_na_idx] <- wb_data_sub$GiZ

#####################################################################################
#library(spdep)

# Subset districts with non-missing index
non_na_idx <- !is.na(wb_data$Immunisation_Index)
wb_data_sub <- wb_data[non_na_idx, ]

# Create neighbor list and spatial weights
nb_sub <- poly2nb(wb_data_sub)
lw_sub <- nb2listw(nb_sub, style = "W")

# Compute Getis-Ord Gi* Z-scores
gi_star <- localG(wb_data_sub$Immunisation_Index, lw_sub)

# Add to spatial data
wb_data_sub$GiZ <- as.numeric(gi_star)
wb_data$GiZ[non_na_idx] <- wb_data_sub$GiZ


########################################################################################
#GI_TABLE

library(dplyr)
library(writexl)

# Subset non-NA districts
non_na_idx <- !is.na(wb_data$Immunisation_Index)
wb_data_sub <- wb_data[non_na_idx, ]

# Create significance column based on Gi* Z-score
wb_data_sub <- wb_data_sub %>%
  mutate(GiSig = case_when(
    GiZ > 2.58  ~ "Hotspot (p<0.01)",
    GiZ > 1.96  ~ "Hotspot (p<0.05)",
    GiZ < -2.58 ~ "Coldspot (p<0.01)",
    GiZ < -1.96 ~ "Coldspot (p<0.05)",
    TRUE        ~ "Not significant"
  ))

# Select columns to include in the table
gi_table <- wb_data_sub %>%
  st_drop_geometry() %>% 
  select(NAME, Immunisation_Index, GiZ, GiSig) %>%
  arrange(desc(GiZ)) 
print(gi_table)

write_xlsx(gi_table, "WestBengal_Gi_Values_Table.xlsx")

##########################################################################################
#LISA TABLE

#library(spdep)
#library(dplyr)
#library(sf)
#library(writexl)

# ---------------------------------
# Subset non-missing districts
# ---------------------------------
non_na_idx <- !is.na(wb_data$Immunisation_Index)
wb_data_sub <- wb_data[non_na_idx, ]

# ---------------------------------
# Neighbours & weights
# ---------------------------------
nb_sub <- poly2nb(wb_data_sub)
lw_sub <- nb2listw(nb_sub, style = "W")

# ---------------------------------
# Local Moran’s I (LISA)
# ---------------------------------
local_moran <- localmoran(wb_data_sub$Immunisation_Index, lw_sub)

# ---------------------------------
# Add LISA results to spatial data
# ---------------------------------
wb_data_sub$Ii   <- local_moran[, 1] 
wb_data_sub$Z_Ii <- local_moran[, 4] 
wb_data_sub$Ii_p <- local_moran[, 5] 

# ---------------------------------
# LISA cluster classification
# ---------------------------------
mean_index <- mean(wb_data_sub$Immunisation_Index)

wb_data_sub$LISA_cluster <- factor(
  ifelse(wb_data_sub$Z_Ii >  1.96 & wb_data_sub$Immunisation_Index > mean_index, "High-High",
         ifelse(wb_data_sub$Z_Ii >  1.96 & wb_data_sub$Immunisation_Index < mean_index, "Low-Low",
                ifelse(wb_data_sub$Z_Ii < -1.96 & wb_data_sub$Immunisation_Index > mean_index, "High-Low",
                       "Low-High"))),
  levels = c("High-High", "Low-Low", "High-Low", "Low-High")
)

# ---------------------------------
# LISA significance labels
# ---------------------------------
wb_data_sub$LISA_Sig <- case_when(
  wb_data_sub$Z_Ii >  2.58 ~ "Significant (p<0.01)",
  wb_data_sub$Z_Ii >  1.96 ~ "Significant (p<0.05)",
  wb_data_sub$Z_Ii < -2.58 ~ "Significant (p<0.01)",
  wb_data_sub$Z_Ii < -1.96 ~ "Significant (p<0.05)",
  TRUE                     ~ "Not significant"
)

# ---------------------------------
# Create LISA table
# ---------------------------------
lisa_table <- wb_data_sub %>%
  st_drop_geometry() %>%
  select(
    NAME,
    Immunisation_Index,
    Ii,
    Z_Ii,
    Ii_p,
    LISA_cluster,
    LISA_Sig
  ) %>%
  arrange(desc(Z_Ii))

print(lisa_table)

write_xlsx(lisa_table, "WestBengal_LISA_Values_Table.xlsx")

##############################################################################################

#PLOTS

# -----------------------------
# Plot maps with tmap
# -----------------------------
tmap_mode("plot") 

tm_shape(wb_data) +
  tm_polygons("Immunisation_Index", palette = "YlGnBu",
              style = "quantile", title = "Immunisation Index") +
  tm_layout(main.title = "West Bengal Immunisation Index")

# LISA Cluster Map
tm_shape(wb_data) +
  tm_polygons("LISA_cluster", palette = c("red","blue","pink","lightblue"),
              title = "LISA Cluster") +
  tm_layout(main.title = "Local Moran's I Clusters")

# Gi* Hotspot / Coldspot
tm_shape(wb_data) +
  tm_polygons("GiZ", palette = "-RdBu", midpoint = 0,
              title = "Gi* Z-score") +
  tm_layout(main.title = "West Bengal Spatial Immunization Pattern Hotspot / Coldspot")


#####################################################################################

tmap_mode("plot")

tm_shape(wb_data) +
  tm_polygons(
    col = "Immunisation_Index",     
    palette = "YlGnBu",
    style = "quantile",              
    n = 5,                           
    title = "Immunisation Index",
    border.col = "grey40",
    lwd = 0.6
  ) +
  tm_text(
    "NAME",
    size = 0.6,
    col = "black"
  ) +
  tm_layout(
    main.title = "West Bengal: District-wise Immunisation Index",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

#####################################################################################

# -----------------------------
# Combined Map: Index + LISA + Gi* + Labels
# -----------------------------
tmap_mode("plot")  # static map

tm_shape(wb_data) +
  # Base polygons colored by Immunisation Index
  tm_polygons(
    col = "Immunisation_Index",
    palette = "YlGnBu",
    style = "quantile",
    title = "Immunisation Index",
    border.col = "grey50",   
    lwd = 0.5
  ) +
  # Overlay LISA clusters
  tm_borders(lwd = 0.5, col = "black") +  
  tm_shape(wb_data) +
  tm_polygons(
    col = "LISA_cluster",
    palette = c("red","blue","pink","lightblue"),
    alpha = 0.3,
    title = "LISA Cluster"
  ) +
  # Gi* hotspots/coldspots overlay
  tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ",
    palette = "-RdBu",
    style = "pretty",
    midpoint = 0,
    alpha = 0.4,
    title = "Gi* Z-score"
  ) +
  # District labels
  tm_text(
    text = "NAME",
    size = 0.6,
    col = "black",
    shadow = TRUE
  ) +
  # Layout and title
  tm_layout(
    main.title = "West Bengal: Immunisation Index with LISA & Gi* Hotspots",
    legend.outside = TRUE,
    frame = FALSE
  )

#############################################################################
#COLORING BASED ON IMMUNIZATION INDEX

library(tmap)

# Static plotting
tmap_mode("plot")

# Choropleth map
tm_shape(wb_data) +
  tm_polygons(
    col = "Immunisation_Index",    
    palette = "YlGnBu",            
    style = "quantile",             
    title = "Immunisation Index",
    border.col = "grey40",          
    lwd = 0.5
  ) +
  tm_text(
    "NAME", size = 0.6, col = "black"
  ) +
  tm_layout(
    main.title = "West Bengal: Immunisation Index by District",
    legend.outside = TRUE,
    frame = FALSE
  )


tm_shape(wb_data) +
  tm_polygons(
    col = "LISA_cluster",
    palette = c("red","blue","pink","lightblue"),
    title = "LISA Cluster",
    border.col = "grey50",
    lwd = 0.5
  ) +
  tm_text(
    text = "NAME",
    size = 0.6,
    col = "black",
    shadow = TRUE
  ) +
  tm_layout(
    main.title = "West Bengal: Spatial Immunization Pattern-(LISA) Clusters",
    main.title.size = 1.5,
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# -----------------------------
# Gi* Hotspot / Coldspot Map
# -----------------------------
tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ",
    palette = "-RdBu",
    midpoint = 0,
    style = "pretty",
    title = "Gi* Z-score",
    border.col = "grey50",
    lwd = 0.5
  ) +
  tm_text(
    text = "NAME",
    size = 0.6,
    col = "black",
    shadow = TRUE
  ) +
  tm_layout(
    main.title = "West Bengal: Spatial Immunization Pattern-(Hotspot / Coldspot)",
    main.title.size = 1.5,
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

##################################################################################################
# Separate Files for Rural Urban Public and Private

library(readxl)
library(writexl)
library(dplyr)

# Read
immun <- read_excel("6B.District_Average_Immunisation_Index_By_Category.xlsx")

# -----------------------------
# CREATE RURAL INDEX FILE
# -----------------------------
rural_data <- immun %>% 
  filter(Category == "Rural") %>% 
  select(District, Category, Index) %>% 
  arrange(District)

write_xlsx(rural_data, "Rural_Immunisation_Index.xlsx")
print(paste("Rural records:", nrow(rural_data)))

# -----------------------------
# CREATE URBAN INDEX FILE
# -----------------------------
urban_data <- immun %>% 
  filter(Category == "Urban") %>% 
  select(District, Category, Index) %>% 
  arrange(District)

write_xlsx(urban_data, "Urban_Immunisation_Index.xlsx")
print(paste("Urban records:", nrow(urban_data)))


# ===================================================
# CREATE SEPARATE EXCEL FILES: Public & Private 
# ===================================================

library(readxl)
library(writexl)
library(dplyr)

# Read
immun <- read_excel("6B.District_Average_Immunisation_Index_By_Category.xlsx")

# CREATE PUBLIC FACILITY FILE
public_data <- immun %>% 
  filter(Category == "Public_Facility") %>%  
  select(District, Category, Index) %>%
  arrange(District)

write_xlsx(public_data, "Public_Facility_Immunisation_Index.xlsx")
print(paste("Public records:", nrow(public_data)))

# CREATE PRIVATE FACILITY FILE
private_data <- immun %>% 
  filter(Category == "Private_Facility") %>%  # ← FIXED: underscore!
  select(District, Category, Index) %>%
  arrange(District)

write_xlsx(private_data, "Private_Facility_Immunisation_Index.xlsx")
print(paste("Private records:", nrow(private_data)))

###############################################################################
#IMMUNISATION INDEX MAPS Rural & Urban, Public & Private

# Libraries
library(sf)
library(dplyr)
library(readxl)
library(tmap)

# Read West Bengal district shapefile
wb_shp <- st_read("District_shape_West_Bengal.shp")

# Check district name field
names(wb_shp)

# Read immunisation index file (by category)
immun <- read_excel("6B.District_Average_Immunisation_Index_By_Category.xlsx")

# District name harmonisation
district_map <- c(
  "Alipurduar" = "Alipurduar",
  "Bankura" = "Bankura",
  "Birbhum" = "Birbhum",
  "Coochbehar" = "Koch Bihar",
  "Dakshin Dinajpur" = "Dakshin Dinajpur",
  "Darjeeling" = "Darjiling",
  "Hooghly" = "Hugli",
  "Howrah" = "Haora",
  "Jalpaiguri" = "Jalpaiguri",
  "Jhargram" = "Jhargram",
  "Kalimpong" = "Kalimpong",
  "Kolkata" = "Kolkata",
  "Malda" = "Maldah",
  "Murshidabad" = "Murshidabad",
  "Nadia" = "Nadia",
  "North 24 Pargana" = "North 24 Parganas",
  "Paschim Bardhaman" = "Paschim Barddhaman",
  "Paschim Medinipur" = "Pashchim Medinipur",
  "Purba Bardhaman" = "Purba Barddhaman",
  "Purba Medinipur" = "Purba Medinipur",
  "Purulia" = "Puruliya",
  "South 24 Pargana" = "South 24 Parganas",
  "Uttar Dinajpur" = "Uttar Dinajpur"
)

immun <- immun %>%
  mutate(DISTRICT_2011 = recode(District, !!!district_map)) %>%
  filter(!is.na(DISTRICT_2011))

# Aggregate indices by category
index_wide <- immun %>%
  group_by(DISTRICT_2011, Category) %>%
  summarise(Index = mean(Index, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Category,
    values_from = Index
  )

# Merge with shapefile
wb_data <- wb_shp %>%
  left_join(index_wide, by = c("NAME" = "DISTRICT_2011")) %>%
  mutate(
    Rural            = as.numeric(Rural),
    Urban            = as.numeric(Urban),
    Public_Facility  = as.numeric(Public_Facility),
    Private_Facility = as.numeric(Private_Facility)
  )

# 6. Static plotting
tmap_mode("plot")

# RURAL IMMUNISATION INDEX
tm_shape(wb_data) +
  tm_polygons(
    col = "Rural",
    palette = "YlGnBu",
    style = "quantile",
    title = "Rural Index"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Rural Immunisation Index (Birth Dose)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# URBAN IMMUNISATION INDEX
tm_shape(wb_data) +
  tm_polygons(
    col = "Urban",
    palette = "YlGnBu",
    style = "quantile",
    title = "Urban Index"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Urban Immunisation Index (Birth Dose)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# PUBLIC FACILITY IMMUNISATION INDEX
tm_shape(wb_data) +
  tm_polygons(
    col = "Public_Facility",
    palette = "YlGnBu",
    style = "quantile",
    title = "Public Facility Index"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Public Facility Immunisation Index (Birth Dose)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# PRIVATE FACILITY IMMUNISATION INDEX
tm_shape(wb_data) +
  tm_polygons(
    col = "Private_Facility",
    palette = "YlGnBu",
    style = "quantile",
    title = "Private Facility Index"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Private Facility Immunisation Index (Birth Dose)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )



##############################################################################

#GI* CALIBRATION — SEPARATE FOR EACH CATEGORY

###############################################################################
# Libraries
library(spdep)
library(sf)
library(dplyr)

# Compute Gi* for one variable
compute_gi_star <- function(sf_data, value_col) {
  
  # Subset only non-NA values
  sub_idx <- !is.na(sf_data[[value_col]])
  sf_sub  <- sf_data[sub_idx, ]
  
  # Neighbours and weights
  nb  <- poly2nb(sf_sub)
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Compute Gi*
  gi  <- localG(sf_sub[[value_col]], lw, zero.policy = TRUE)
  
  # Full-length vector with NA
  gi_full <- rep(NA_real_, nrow(sf_data))
  gi_full[sub_idx] <- as.numeric(gi)
  
  return(gi_full)
}


# Compute Gi* separately for each category
wb_data$GiZ_Rural   <- compute_gi_star(wb_data, "Rural")
wb_data$GiZ_Urban   <- compute_gi_star(wb_data, "Urban")
wb_data$GiZ_Public  <- compute_gi_star(wb_data, "Public_Facility")
wb_data$GiZ_Private <- compute_gi_star(wb_data, "Private_Facility")

# Check results
summary(wb_data$GiZ_Rural)
summary(wb_data$GiZ_Urban)
summary(wb_data$GiZ_Public)
summary(wb_data$GiZ_Private)

# Static plotting 
library(tmap)
tmap_mode("plot")

# GI* MAPS
###############################################################################

# -------------------------
# RURAL
# -------------------------
tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ_Rural",
    palette = "-RdBu",
    midpoint = 0,
    style = "pretty",
    title = "Gi* Z-score (Rural)"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Rural Immunisation Hotspots & Coldspots (Gi*)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# -------------------------
# URBAN
# -------------------------
tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ_Urban",
    palette = "-RdBu",
    midpoint = 0,
    style = "pretty",
    title = "Gi* Z-score (Urban)"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Urban Immunisation Hotspots & Coldspots (Gi*)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# -------------------------
# PUBLIC FACILITY
# -------------------------
tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ_Public",
    palette = "-RdBu",
    midpoint = 0,
    style = "pretty",
    title = "Gi* Z-score (Public Facility)"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Public Facility Immunisation Hotspots & Coldspots (Gi*)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )

# -------------------------
# PRIVATE FACILITY
# -------------------------
tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ_Private",
    palette = "-RdBu",
    midpoint = 0,
    style = "pretty",
    title = "Gi* Z-score (Private Facility)"
  ) +
  tm_text("NAME", size = 0.6) +
  tm_layout(
    main.title = "West Bengal: Private Facility Immunisation Hotspots & Coldspots (Gi*)",
    main.title.position = "center",
    legend.outside = TRUE,
    frame = FALSE
  )


###############################################################################
# 17. Save Gi* values


# Libraries
library(sf)
library(spdep)
library(dplyr)
library(openxlsx)

# Prepare spatial object

wb_data_sub <- wb_data %>%
  filter(!st_is_empty(geometry))

# 2. Neighbours & spatial weights
nb <- poly2nb(wb_data_sub)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Compute Gi*
compute_gi <- function(sf_obj, value_col, category_name) {
  
  # Extract variable
  x <- sf_obj[[value_col]]
  
  # Replace NA with 0
  x[is.na(x)] <- 0
  
  # Compute Gi*
  gi_z <- as.numeric(localG(x, lw, zero.policy = TRUE))
  
  # Significance classification
  gi_sig <- case_when(
    gi_z >=  1.96 ~ "Hotspot (p<0.05)",
    gi_z <= -1.96 ~ "Coldspot (p<0.05)",
    TRUE          ~ "Not significant"
  )
  
  # Final table
  out <- sf_obj %>%
    st_drop_geometry() %>%
    transmute(
      NAME,
      Immunisation_Index = x,
      GiZ = gi_z,
      GiSig = gi_sig
    )
  
  # Save
  write.xlsx(
    out,
    paste0("GiStar_", category_name, ".xlsx"),
    overwrite = TRUE
  )
  
  return(out)
}

# Compute Gi* for each category

# RURAL
gi_rural <- compute_gi(
  wb_data_sub,
  value_col   = "Rural",
  category_name = "Rural"
)

# URBAN
gi_urban <- compute_gi(
  wb_data_sub,
  value_col   = "Urban",
  category_name = "Urban"
)

# PUBLIC FACILITY
gi_public <- compute_gi(
  wb_data_sub,
  value_col   = "Public_Facility",
  category_name = "Public_Facility"
)

# PRIVATE FACILITY
gi_private <- compute_gi(
  wb_data_sub,
  value_col   = "Private_Facility",
  category_name = "Private_Facility"
)

cat("Gi* analysis completed and Excel files saved successfully.\n")
