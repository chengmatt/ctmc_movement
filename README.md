# From Here to There and Back Again: Parsimonious Estimation of Environment- and Demography-Dependent Movement in Spatially-Stratified Stock Assessment Models

## Authors
Matthew LH. Cheng, James T. Thorson, Daniel R. Goethel, Curry J. Cunningham

### Abstract
Animal migrations and movement patterns are fundamental drivers that govern the dynamics of marine populations. Movement often varies across demographic attributes and through time in response to biology and environmental forcing, giving rise to a spectrum of connectivity patterns and resultant population structures. For commercially exploited marine fisheries, stock assessment models are widely used to estimate populations status and provide sustainable harvest advice, typically assuming a spatially homogeneous panmictic population and ignoring movement. Spatially-stratified stock assessment models can integrate population and spatial structure, as well as model movement, but their adoption is limited due to increased complexity and difficulties estimating movement. We introduce a parsimonious approach for representing movement using continuous-time Markov chains (CTMC). The theoretical basis for deriving diffusive and taxis-based movement from a CTMC process is provided, then a simulation-estimation experiment is undertaken that simulates age-, time-, and environmentally-driven movement rates to demonstrate the utility of the CTMC movement framework within spatially-stratified stock assessment models. A case-study based on a three-region model of Alaska sablefish (Anoplopoma fimbria) illustrates how the CTMC framework can be utilized in a real-world assessment and data context. Across simulation experiments, the CTMC movement model produced minimal bias and accurately represented complex movement. The sablefish CTMC application revealed age-varying movement dynamics, resulting in improved model fits relative to other movement parameterizations (i.e., an unstructured Markov approach with constant or age-varying movement). Overall, results support broader integration of CTMC-based movement models into assessment platforms, which could help reduce implementation barriers for spatially-stratified stock assessments.

### Files
- .png files represent plots generated from this study.
- ctmc_demo.R generates conceptual figure of the CTMC model
- movement_param_sablefish.R generates models rom the sablefish case-study
- simulate_movement.R simulates data needed for the study (need to run this first to run simulations)
- sim_demo.R plots simulation dynamics
- sim_plots.R plots simulation results
- setup_basic_spatial_em.R sets up a spatial EM for use in simulation
- run_movement_sim.R runs all simulations

#### Prerequisites

Ensure the following packages are installed:

```
install.packages("devtools")       # Development tools
install.packages("TMB")            # Template Model Builder
install.packages("RTMB")           # R interface to TMB
TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip") # Optional: multivariate OSA distributions
remotes::install_github("fishfollower/compResidual/compResidual") # Optional OSA residuals
```

#### Installing SPoRC
```
devtools::install_github(
  "chengmatt/SPoRC@dev-movement",
  dependencies = c("Depends", "Imports")
)
```
