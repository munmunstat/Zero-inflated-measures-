# Zero-Inflated Measures
This repository contains R code, simulation output, and real-data analysis for 
zero-inflated extensions of three classical dependence measures:

- Zero-inflated Kendall's $\tau^*_{ZI}$
- Zero-inflated Spearman's $\rho^*_{ZI}$
- Zero-inflated Gini correlation $\Gamma^*_{ZI}$
  
# Here is a simple flow chart for Simulation study:


```mermaid
flowchart TD
    A[SparseDOSSA2 Data Generation: Stool / Vaginal / IBD] --> B[Simulation Setup: n=500, n=1000]
    B --> C[Compute Zero-Inflated Measures: tau*, rho*, Gini*]
    C --> D[Save Intermediate Outputs]
    D --> E[Summaries: RMSE, RMAD]
    E --> F[Create Plots: Box Plots]
    F --> G[Real Data Analysis: American Gut Project]
