# Spike at zero Measures
This repository contains R code, simulation output, and real-data analysis for 
zero-inflated extensions of three classical dependence measures:

- Spike at zero Kendall's $\tau^*_{ZI}$
- Spike at zero Spearman's $\rho^*_{ZI}$
- Spike at zero  Gini correlation $\Gamma^*_{ZI}$
  
# Here is a simple flow chart for Simulation study:


```mermaid
flowchart TD
    A[SparseDOSSA2 Data Generation: Stool / Vaginal / IBD] --> B[Simulation Setup: n=500, n=1000]
    B --> C[Compute Zero-Inflated Measures: tau*, rho*, Gini*]
    C --> D[Save Intermediate Outputs]
    D --> E[Summaries: RMSE, RMAD]
    E --> F[Create Plots: Box Plots]
    F --> G[Real Data Analysis: American Gut Project]
