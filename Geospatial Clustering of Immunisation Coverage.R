# ===================================================
# FULL CLEANED SCRIPT: WEST BENGAL IMMUNISATION INDEX
# ===================================================

# -----------------------------
# 0. Install / load packages
# -----------------------------
# install.packages(c("sf","dplyr","readxl","spdep","tmap","RColorBrewer","terra"))

library(sf)
library(dplyr)
library(readxl)
library(spdep)
library(tmap)
library(RColorBrewer)
library(terra)

# -----------------------------
# 1. Read shapefile and filter West Bengal
# -----------------------------
wb_shp <- st_read("District_shape_West_Bengal.shp")

#wb_shp <- st_read("District_shape_West_Bengal.shp")%>%
  #filter(ST_NM == "West Bengal")

# Check districts
print(unique(wb_shp$NAME))

# -----------------------------
# 2. Read immunisation index Excel
# -----------------------------
immun <- read_excel("District_Cumulative_Immunisation_Risk.xlsx")

# Rename column for consistency
immun <- immun %>% rename(Average_Index = `Average_Index`)

# -----------------------------
# 3. Map Excel districts to shapefile districts
# -----------------------------
# -----------------------------
# 2. District name mapping and category aggregation
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
# 4. Aggregate duplicates
# -----------------------------
immun_agg <- immun %>%
  group_by(DISTRICT_2011) %>%
  summarise(Immunisation_Index = mean(Average_Index, na.rm = TRUE))

# -----------------------------
# 5. Merge aggregated index into shapefile
# -----------------------------
wb_data <- wb_shp %>%
  left_join(immun_agg, by = c("NAME" = "DISTRICT_2011"))

# Convert to numeric
wb_data$Immunisation_Index <- as.numeric(wb_data$Immunisation_Index)

# Initialize LISA and Gi* columns
wb_data$LISA_cluster <- NA
wb_data$GiZ <- NA

# Check merge
wb_data %>% st_drop_geometry() %>% select(NAME, Immunisation_Index)

# -----------------------------
# 6. Create neighbors
# -----------------------------
nb <- poly2nb(wb_data)
lw <- nb2listw(nb, style = "W")

# -----------------------------
# 7. Global Moran's I
# -----------------------------
index_values <- wb_data$Immunisation_Index
global_moran <- moran.test(index_values, lw, na.action = na.omit)
print(global_moran)

# -----------------------------
# 8. Local Moran's I (LISA)
# -----------------------------
non_na_idx <- !is.na(index_values)
wb_data_sub <- wb_data[non_na_idx, ]

nb_sub <- poly2nb(wb_data_sub)
lw_sub <- nb2listw(nb_sub, style = "W")

local_moran <- localmoran(wb_data_sub$Immunisation_Index, lw_sub)

# Add results
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

# Fill into main shapefile
wb_data$LISA_cluster[non_na_idx] <- wb_data_sub$LISA_cluster

# -----------------------------
# 9. Getis-Ord Gi* Hotspot/Coldspot
# -----------------------------
gi_star <- localG(wb_data_sub$Immunisation_Index, lw_sub)
wb_data_sub$GiZ <- as.numeric(gi_star)

wb_data$GiZ[non_na_idx] <- wb_data_sub$GiZ


library(spdep)

# 1. Subset districts with non-missing index
non_na_idx <- !is.na(wb_data$Immunisation_Index)
wb_data_sub <- wb_data[non_na_idx, ]

# 2. Create neighbor list and spatial weights
nb_sub <- poly2nb(wb_data_sub)
lw_sub <- nb2listw(nb_sub, style = "W")

# 3. Compute Getis-Ord Gi* Z-scores
gi_star <- localG(wb_data_sub$Immunisation_Index, lw_sub)

# 4. Add to spatial data
wb_data_sub$GiZ <- as.numeric(gi_star)
wb_data$GiZ[non_na_idx] <- wb_data_sub$GiZ


###########################################################################################
#GI_TABLE

library(dplyr)
library(writexl)

# 1. Subset non-NA districts
non_na_idx <- !is.na(wb_data$Immunisation_Index)
wb_data_sub <- wb_data[non_na_idx, ]

# 2. Create significance column based on Gi* Z-score
wb_data_sub <- wb_data_sub %>%
  mutate(GiSig = case_when(
    GiZ > 2.58  ~ "Hotspot (p<0.01)",
    GiZ > 1.96  ~ "Hotspot (p<0.05)",
    GiZ < -2.58 ~ "Coldspot (p<0.01)",
    GiZ < -1.96 ~ "Coldspot (p<0.05)",
    TRUE        ~ "Not significant"
  ))

# 3. Select columns to include in the table
gi_table <- wb_data_sub %>%
  st_drop_geometry() %>%  # remove geometry for Excel
  select(NAME, Immunisation_Index, GiZ, GiSig) %>%
  arrange(desc(GiZ))  # optional: sort by Z-score

# 4. View the table
print(gi_table)

# 5. Save to Excel
write_xlsx(gi_table, "WestBengal_Gi_Values_Table.xlsx")

###############################################################################################
#LISA TABLE

library(spdep)
library(dplyr)
library(sf)
library(writexl)

# ---------------------------------
# 1. Subset non-missing districts
# ---------------------------------
non_na_idx <- !is.na(wb_data$Immunisation_Index)
wb_data_sub <- wb_data[non_na_idx, ]

# ---------------------------------
# 2. Neighbours & weights
# ---------------------------------
nb_sub <- poly2nb(wb_data_sub)
lw_sub <- nb2listw(nb_sub, style = "W")

# ---------------------------------
# 3. Local Moran’s I (LISA)
# ---------------------------------
local_moran <- localmoran(wb_data_sub$Immunisation_Index, lw_sub)

# ---------------------------------
# 4. Add LISA results to spatial data
# ---------------------------------
wb_data_sub$Ii   <- local_moran[, 1]   # Local Moran's I
wb_data_sub$Z_Ii <- local_moran[, 4]   # Z-score
wb_data_sub$Ii_p <- local_moran[, 5]   # p-value

# ---------------------------------
# 5. LISA cluster classification
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
# 6. LISA significance labels
# ---------------------------------
wb_data_sub$LISA_Sig <- case_when(
  wb_data_sub$Z_Ii >  2.58 ~ "Significant (p<0.01)",
  wb_data_sub$Z_Ii >  1.96 ~ "Significant (p<0.05)",
  wb_data_sub$Z_Ii < -2.58 ~ "Significant (p<0.01)",
  wb_data_sub$Z_Ii < -1.96 ~ "Significant (p<0.05)",
  TRUE                     ~ "Not significant"
)

# ---------------------------------
# 7. Create LISA table
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

# ---------------------------------
# 8. View & save
# ---------------------------------
print(lisa_table)

write_xlsx(lisa_table, "WestBengal_LISA_Values_Table.xlsx")

########################################################################################################

##PLOTS


# -----------------------------
# 10. Plot maps with tmap
# -----------------------------
tmap_mode("plot")  # static

# Choropleth
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
    col = "Immunisation_Index",     # ONLY index controls color
    palette = "YlGnBu",
    style = "quantile",              # or "pretty" / "fixed"
    n = 5,                            # number of classes
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
  # 1. Base polygons colored by Immunisation Index
  tm_polygons(
    col = "Immunisation_Index",
    palette = "YlGnBu",
    style = "quantile",
    title = "Immunisation Index",
    border.col = "grey50",   # district borders
    lwd = 0.5
  ) +
  # 2. Overlay LISA clusters
  tm_borders(lwd = 0.5, col = "black") +  # optional for clarity
  tm_shape(wb_data) +
  tm_polygons(
    col = "LISA_cluster",
    palette = c("red","blue","pink","lightblue"),
    alpha = 0.3,
    title = "LISA Cluster"
  ) +
  # 3. Gi* hotspots/coldspots overlay (semi-transparent)
  tm_shape(wb_data) +
  tm_polygons(
    col = "GiZ",
    palette = "-RdBu",
    style = "pretty",
    midpoint = 0,
    alpha = 0.4,
    title = "Gi* Z-score"
  ) +
  # 4. District labels
  tm_text(
    text = "NAME",
    size = 0.6,
    col = "black",
    shadow = TRUE
  ) +
  # 5. Layout and title
  tm_layout(
    main.title = "West Bengal: Immunisation Index with LISA & Gi* Hotspots",
    legend.outside = TRUE,
    frame = FALSE
  )

#############################################################################
#COLORING BASED ON IMMUNIZATION INDEX

library(tmap)

# Ensure we are in static plotting mode
tmap_mode("plot")

# Simple choropleth map
tm_shape(wb_data) +
  tm_polygons(
    col = "Immunisation_Index",    # variable to color districts
    palette = "YlGnBu",            # Yellow → Green → Blue
    style = "quantile",             # color by quantiles
    title = "Immunisation Index",
    border.col = "grey40",          # district borders
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
# 3. Gi* Hotspot / Coldspot Map
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

#######################################################################################################################
# Separate Files for Rural Urban Public and Private

library(readxl)
library(writexl)
library(dplyr)

# Read the original file
immun <- read_excel("6B.District_Average_Immunisation_Index_By_Category.xlsx")

# -----------------------------
# 1. CREATE RURAL INDEX FILE
# -----------------------------
rural_data <- immun %>% 
  filter(Category == "Rural") %>% 
  select(District, Category, Index) %>% 
  arrange(District)

write_xlsx(rural_data, "Rural_Immunisation_Index.xlsx")
print(paste("✅ Rural records:", nrow(rural_data)))

# -----------------------------
# 2. CREATE URBAN INDEX FILE
# -----------------------------
urban_data <- immun %>% 
  filter(Category == "Urban") %>% 
  select(District, Category, Index) %>% 
  arrange(District)

write_xlsx(urban_data, "Urban_Immunisation_Index.xlsx")
print(paste("✅ Urban records:", nrow(urban_data)))


# ===================================================
# CREATE SEPARATE EXCEL FILES: Public & Private 
# ===================================================

library(readxl)
library(writexl)
library(dplyr)

# Read the original file
immun <- read_excel("6B.District_Average_Immunisation_Index_By_Category.xlsx")

# 1. CREATE PUBLIC FACILITY FILE - CORRECT CATEGORY NAME
public_data <- immun %>% 
  filter(Category == "Public_Facility") %>%  # ← FIXED: underscore!
  select(District, Category, Index) %>%
  arrange(District)

write_xlsx(public_data, "Public_Facility_Immunisation_Index.xlsx")
print(paste("✅ Public records:", nrow(public_data)))

# 2. CREATE PRIVATE FACILITY FILE - CORRECT CATEGORY NAME  
private_data <- immun %>% 
  filter(Category == "Private_Facility") %>%  # ← FIXED: underscore!
  select(District, Category, Index) %>%
  arrange(District)

write_xlsx(private_data, "Private_Facility_Immunisation_Index.xlsx")
print(paste("✅ Private records:", nrow(private_data)))

###############################################################################
#IMMUNISATION INDEX MAPS Rural & Urban, Public & Private

###############################################################################
# 0. Libraries
###############################################################################
library(sf)
library(dplyr)
library(readxl)
library(tmap)

###############################################################################
# 1. Read West Bengal district shapefile
###############################################################################
wb_shp <- st_read("District_shape_West_Bengal.shp")

# Check district name field
names(wb_shp)
# EXPECTED: NAME

###############################################################################
# 2. Read immunisation index file (by category)
###############################################################################
immun <- read_excel("6B.District_Average_Immunisation_Index_By_Category.xlsx")

# EXPECTED COLUMNS:
# District | Category | Index

###############################################################################
# 3. District name harmonisation
###############################################################################
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

###############################################################################
# 4. Aggregate indices by category
###############################################################################
index_wide <- immun %>%
  group_by(DISTRICT_2011, Category) %>%
  summarise(Index = mean(Index, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Category,
    values_from = Index
  )

# EXPECTED columns after pivot:
# DISTRICT_2011 | Rural | Urban | Public_Facility | Private_Facility

###############################################################################
# 5. Merge with shapefile
###############################################################################
wb_data <- wb_shp %>%
  left_join(index_wide, by = c("NAME" = "DISTRICT_2011")) %>%
  mutate(
    Rural            = as.numeric(Rural),
    Urban            = as.numeric(Urban),
    Public_Facility  = as.numeric(Public_Facility),
    Private_Facility = as.numeric(Private_Facility)
  )

###############################################################################
# 6. Static plotting mode
###############################################################################
tmap_mode("plot")

###############################################################################
# 7. RURAL IMMUNISATION INDEX
###############################################################################
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

###############################################################################
# 8. URBAN IMMUNISATION INDEX
###############################################################################
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

###############################################################################
# 9. PUBLIC FACILITY IMMUNISATION INDEX
###############################################################################
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

###############################################################################
# 10. PRIVATE FACILITY IMMUNISATION INDEX
###############################################################################
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
# 11. Libraries
###############################################################################
library(spdep)
library(sf)
library(dplyr)

###############################################################################
# 12. Function to compute Gi* safely for one variable
###############################################################################
compute_gi_star <- function(sf_data, value_col) {
  
  # 1. Subset only non-NA values
  sub_idx <- !is.na(sf_data[[value_col]])
  sf_sub  <- sf_data[sub_idx, ]
  
  # 2. Neighbours and weights
  nb  <- poly2nb(sf_sub)
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # 3. Compute Gi*
  gi  <- localG(sf_sub[[value_col]], lw, zero.policy = TRUE)
  
  # 4. Create full-length vector with NA
  gi_full <- rep(NA_real_, nrow(sf_data))
  gi_full[sub_idx] <- as.numeric(gi)
  
  return(gi_full)
}

###############################################################################
# 13. Compute Gi* separately for each category (CORRECT WAY)
###############################################################################
wb_data$GiZ_Rural   <- compute_gi_star(wb_data, "Rural")
wb_data$GiZ_Urban   <- compute_gi_star(wb_data, "Urban")
wb_data$GiZ_Public  <- compute_gi_star(wb_data, "Public_Facility")
wb_data$GiZ_Private <- compute_gi_star(wb_data, "Private_Facility")

###############################################################################
# 14. Check results
###############################################################################
summary(wb_data$GiZ_Rural)
summary(wb_data$GiZ_Urban)
summary(wb_data$GiZ_Public)
summary(wb_data$GiZ_Private)

###############################################################################
# 15. Static plotting mode
###############################################################################
library(tmap)
tmap_mode("plot")

###############################################################################
# 16. GI* MAPS
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
# 17. Save Gi* values to Excel (Category-wise)

###############################################################################
# 0. Libraries
###############################################################################
library(sf)
library(spdep)
library(dplyr)
library(openxlsx)

###############################################################################
# 1. Prepare spatial object (remove empty geometries)
###############################################################################
wb_data_sub <- wb_data %>%
  filter(!st_is_empty(geometry))

###############################################################################
# 2. Neighbours & spatial weights (ONCE)
###############################################################################
nb <- poly2nb(wb_data_sub)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

###############################################################################
# 3. Helper function to compute Gi* safely
###############################################################################
compute_gi <- function(sf_obj, value_col, category_name) {
  
  # Extract variable
  x <- sf_obj[[value_col]]
  
  # Replace NA with 0 (MANDATORY for localG)
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
  
  # Save to Excel
  write.xlsx(
    out,
    paste0("GiStar_", category_name, ".xlsx"),
    overwrite = TRUE
  )
  
  return(out)
}

###############################################################################
# 4. Compute Gi* for each category
###############################################################################

# ---- RURAL ----
gi_rural <- compute_gi(
  wb_data_sub,
  value_col   = "Rural",
  category_name = "Rural"
)

# ---- URBAN ----
gi_urban <- compute_gi(
  wb_data_sub,
  value_col   = "Urban",
  category_name = "Urban"
)

# ---- PUBLIC FACILITY ----
gi_public <- compute_gi(
  wb_data_sub,
  value_col   = "Public_Facility",
  category_name = "Public_Facility"
)

# ---- PRIVATE FACILITY ----
gi_private <- compute_gi(
  wb_data_sub,
  value_col   = "Private_Facility",
  category_name = "Private_Facility"
)

###############################################################################
# 5. Done
###############################################################################
cat("Gi* analysis completed and Excel files saved successfully.\n")
