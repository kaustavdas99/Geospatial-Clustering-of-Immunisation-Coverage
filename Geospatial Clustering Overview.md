# Geospatial Clustering and Spatial Autocorrelation Analysis of Immunisation Coverage

## Overview

This repository presents a comprehensive geospatial analysis workflow for examining district-level variation in child immunisation coverage. The analysis integrates spatial data processing, clustering concepts, and spatial statistical techniques to identify patterns of similarity, clustering, and spatial dependence across regions.

The primary objective of this work is to understand how immunisation performance is distributed geographically and whether neighbouring districts exhibit similar or contrasting patterns. The workflow combines geospatial visualization with spatial autocorrelation statistics to generate actionable public health insights.

---

## What is Geospatial Clustering?

Geospatial clustering refers to the process of grouping geographic units (such as districts or regions) based on similarity in observed attributes while considering their spatial relationships.

Unlike traditional clustering, geospatial clustering explicitly accounts for the spatial arrangement of data points. In this context, districts are not independent entities—nearby districts may influence each other due to shared infrastructure, policy environments, or socio-economic conditions.

This project applies spatial clustering concepts to:
- Detect regions with consistently high or low immunisation coverage  
- Identify spatial inequalities  
- Reveal geographically structured patterns in public health data  

---

## Objectives

- To map district-level immunisation indices across West Bengal  
- To detect spatial clustering using statistical measures  
- To identify hotspots and coldspots of immunisation performance  
- To evaluate spatial dependence using global and local indicators  

---

## Data Description

The analysis uses:

- **District-level shapefile** of West Bengal  
- **Immunisation index dataset** derived from vaccination indicators  

The immunisation index represents aggregated vaccination coverage measures and is analyzed across:
- Overall population  
- Rural and Urban populations  
- Public and Private healthcare facilities  

---

## Methodological Framework

### 1. Data Integration

- Spatial boundaries (shapefile) are merged with immunisation data  
- District name harmonization ensures consistency across datasets  
- Aggregation is performed where multiple entries exist  

---

### 2. Spatial Neighbourhood Construction

Spatial relationships are defined using adjacency-based neighbours:

- Districts sharing boundaries are treated as neighbours  
- A spatial weights matrix is constructed using row-standardized weights  

This forms the basis for all spatial statistical computations.

---

## Spatial Autocorrelation Analysis

### Global Moran’s I

Global Moran’s I measures the overall spatial autocorrelation in the dataset.

- It evaluates whether similar values cluster together across the entire study region  
- Values range from:
  - **+1** → Strong clustering (similar values grouped together)  
  - **0** → Random spatial pattern  
  - **-1** → Dispersion (dissimilar values are neighbours)  

#### Interpretation:
A significant positive Moran’s I indicates that districts with similar immunisation levels are spatially clustered rather than randomly distributed.

---

### Local Moran’s I (LISA)

Local Indicators of Spatial Association (LISA) provide a localized version of Moran’s I.

Instead of a single global value, LISA identifies:
- Specific districts contributing to clustering  
- Local pockets of spatial association  

#### LISA Cluster Types:

- **High–High (HH):** High-value district surrounded by high-value neighbours  
- **Low–Low (LL):** Low-value district surrounded by low-value neighbours  
- **High–Low (HL):** High-value district surrounded by low-value neighbours (spatial outlier)  
- **Low–High (LH):** Low-value district surrounded by high-value neighbours (spatial outlier)  

#### Significance:
Statistical thresholds (Z-scores) are used to determine whether clusters are significant.

---

### Getis-Ord Gi\* Statistic (Hotspot Analysis)

The Getis-Ord Gi\* statistic identifies spatial concentrations of high or low values.

- It produces a **Z-score** for each district  
- Indicates whether a district is part of a statistically significant cluster  

#### Interpretation:

- **High positive Z-score:** Hotspot (high values clustered together)  
- **Low negative Z-score:** Coldspot (low values clustered together)  
- **Near zero:** No significant clustering  

#### Significance Levels:
- |Z| > 1.96 → Significant at p < 0.05  
- |Z| > 2.58 → Significant at p < 0.01  

---

## Outputs and Visualizations

The workflow generates multiple spatial outputs:

### 1. Choropleth Maps
- District-wise immunisation index visualization  

### 2. LISA Cluster Maps
- Identification of spatial clusters and outliers  

### 3. Gi\* Hotspot Maps
- Detection of high-risk and low-performance regions  

### 4. Combined Spatial Maps
- Overlay of index, LISA clusters, and hotspot analysis  

### 5. Tabular Outputs
- Exported Excel files containing:
  - Gi\* Z-scores and significance levels  
  - LISA statistics and cluster classifications  

---

## Category-wise Spatial Analysis

The workflow extends analysis across different population segments:

- Rural vs Urban  
- Public vs Private healthcare facilities  

Each category is independently analyzed using:
- Spatial clustering  
- Gi\* hotspot detection  

This enables comparative assessment of spatial inequalities across service delivery systems.

---

## Tools and Libraries

- **R Programming Language**
- Key packages:
  - `sf` – spatial data handling  
  - `spdep` – spatial autocorrelation analysis  
  - `tmap` – thematic mapping  
  - `dplyr` – data manipulation  
  - `readxl`, `writexl` – data input/output  

---

## Workflow Summary

1. Load spatial and immunisation datasets  
2. Harmonize and merge district-level data  
3. Construct spatial neighbour relationships  
4. Compute Global Moran’s I  
5. Perform Local Moran’s I (LISA) analysis  
6. Conduct Getis-Ord Gi\* hotspot analysis  
7. Generate maps and export analytical tables  

---

## Key Insights

- Spatial clustering reveals regional disparities in immunisation coverage  
- Hotspot analysis identifies priority intervention areas  
- Local spatial statistics provide district-level actionable insights  
- Category-wise analysis highlights structural inequalities in healthcare access  

---

## Future Directions

- Incorporation of spatial regression models  
- Integration with socio-economic and demographic variables  
- Application of Bayesian spatial models  
- Extension to temporal (spatio-temporal) analysis  

---
