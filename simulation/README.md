# Simulation study 

This folder implements the controlled simulation study of Appendix B.
It checks, under a known data-generating mechanism, that the
operator-matched estimators recover the truth when well specified and
stay bounded under operator misspecification, and compares them against
the additive model and the standard remedies for spatial confounding.

To reproduce:
1. Run run_simulation.R to run the simulation study (this takes time).
2. Run make_results.R to produce the figure and the summaries. 
