# SequentialDesigns
Implementation and Simulations for new sequential and adaptive designs of clinical trials

These scripts shall provide an analysis of our newly developed method for sequential one-arm clinical trials with multivariate time-to event endpoints.\
The main scripts are "Comparison_InterimAnalyses.R" and "SimulationLoop.R".
The former simulates a fixed scenario, compares the results for two different choices of point of time for the interim analysis and provides figures to look at empirical distributions and correlations.
The latter executes simulations for a wide range of scenarios and computes empirical alpha-levels.\
The "setup_[...]" scripts are only used to setup the marginal and joint distributions, especially the conditional hazard functions.
These functions are essential for the test statistics developed in the manuscript
