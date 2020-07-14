# Arterial Signaling
MATLAB codes for a logic-based model of arterial wall signaling, which accompany the article:  Irons L., &amp; Humphrey, J. D. (2020). `Cell signaling model for arterial mechanobiology.' PLOS Computational Biology.

These codes are based on open source code available at: 
       github.com/saucermanlab/netflux
for simulating logic-based signaling networks with normalized Hill ODEs. The methods were described originally in: 
       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). `Modeling cardiac B-adrenergic signaling with normalized-Hill differential equations: comparison with a biochemical model.' BMC Systems Biology.

## Files include:

BaseScript_master.m - Script to simulate timecourses of our ODE system. Inputs and network parameters can be varied as desired. 

IO_Script.m - Script to simulate input-output relations for our optimal basal parameters and quantifies qualitative matches to experimental observations (Fig 2). The values of b, p, n, and EC50 were looped over to find this chosen/optimal combination, as described in S2 Appendix. 

CreateValidationMatrix.m - Function used in IO_Script to store qualitative experimental observations.

Knockdown_Script.m - Script to simulate node perturbations by reducing Ymax from 1 to 0.1 for each node in turn, whilst recording changing steady states of the other nodes (Fig 3).

DefineOrdering.m - Function used in Knockdown_Script and IO_Script to reorder the species in the output matrix.

Stress_Surfaces.m - Script to simulate the absolute and fold-change responses of a subset of species under varying magnitudes of basal inputs and stress perturbations (S3 Appendix)

StressAngII_doses.m - Script to simulate species responses to differing Stress and exogenous AngII combinations (Fig 4)

DoseResponseSurfaces_StressAngII.m - Script to simulate dose response surfaces as combinations of Stress and AngII are varied through [0,1] (Fig 5A)

Comparison_Wu_Script.m - Script to simulates Col1mRNA and Col3mRNA levels at 3 levels of Stress, with and without exogenous AngII, and compares model predictions to experimental data from Wu et al (2014), as described in our accompanying publication (Fig 6).

ODElist_final.mat - File containing the list of ODEs for our specific system, generated using the open source `Netflux' code, which converts a list of logic statements into a system of ODEs using normalized Hill functions. It is used in all of the above scripts.

reactions_final.mat - File containing the list of logic statements, as well as the list of species names, tau and ymax parameters. It was generated using the open source `Netflux' code and is used in all of the above scripts.
