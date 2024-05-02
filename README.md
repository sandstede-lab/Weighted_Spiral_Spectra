# Weighted_Spiral_Spectra

Matlab code to compute the spectra of spiral waves using and exponentially weighted preconditioner.


## Matlab Scripts
The main components of the code are:

``compute_spectra_weighted.m``: Computes the eigenvalues of a spiral wave using the exponentially weighted operator. 

``solve_spiral_wave.m``: Uses Newton's method to solve for a spiral wave solution in a polar co-rotating frame. 

``spiral_pseudospec_conditionNum.m``: Computes the condition number and minimum svd values of the linearized spiral wave operator across a grid of defined points in the complex plane.


### Auxiliary Functions
The following are helper functions that are called by the main scripts.

``Barkley_2D_rotating.m``: Defines the right-hand side of the Barkley model equation. Used to solve for spiral solutions in the co-rotating polar frame.

``Barkley_jacobian_weighted_operator``: Defines the weighted Jacobian matrix of the Barkley model.

``ComputeLinearOperator_shortGrid.m``: Defines the relevant 2D polar differentiation operators used for computing eigenvalues. 

``ComputeLinearOperator.m``: Defines the relevant 2D polar differentiation operators used for solving for spiral wave solutions.

``fourdif.m``: Defines 1D Fourier differentiation matrices

``plot_spiral.m``: Function to plot spiral wave solutions.


### Matlab Files
``Barkley_spiral_r25_h0p05_delta0p2_b0p001_a0p7_ep0p02.mat``: Spiral wave solution to the Barkley model with radius 25. Other parameters are define in the ``par`` structure and in the file name. Specifically, $\delta = 0.2$, $b = 0.001$, $a = 0.7$, and $\epsilon = 0.02$. 

``Barkley_spiral_r50_h0p05_delta0p2_b0p001_a0p7_ep0p02_positiveOmega.mat``: Spiral wave solution to the Barkley model with radius 50. Other parameters are define in the ``par`` structure and are equivalent to the radius 25 case.


