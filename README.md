# TSA_SGT
Transient Stability Assessment based on Small Gain Theorem

Description of Files:

compute_gain.m: MATLAB code, which contains building the 2nd order swing equation 
optimize_bound.jl: Julia code which use IPOPT to optimize bound and compute the perturbation bound (requires MATLAB and Julia with IPOPT and JuMP downloaded)
TSA_SGT_2bus.m: MATLAB code just for 2 bus system for plotting and analysis

TSA_SGT_Adaptive.m: Future work, which will be an extension of the current implementation (work in progress)
plot_bounding.m: used for the bounding plot in the paper
