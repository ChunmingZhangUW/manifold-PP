readme_and_matlab_codes_for_PP-manifold.txt

This Readme description provides a summary of Matlab files for producing Figure 2 (left panel and right panel) in the paper. Inputs, such as sample size (n), data dimension (p), and types of data distribution used in codes for other figures, can be similarly and manually adjusted.

The implementations involve 3 steps:

Step 1: Download the 'Manopt' package, available at https://www.manopt.org/index.html. Follow instructions there to install the package.

Step 2: Download the Matlab functions (listed below) from GitHub at https://github.com/ChunmingZhangUW/manifold-PP:

- EDF_D_dim_func.m
- KS_stat_1_sample.m
- Least_squares_min_w_quad_equa_constraint.m
- Manopt_Z_crt_opt.m
- X_matrix_sim.m
- demo_Matlab_Figure_2.m
- est_pdf_from_histogram.m
- mean_cov_mixture_Gaussians_D_dim.m
- z_vector_critical_dep_on_X_G.m
- z_vector_feasible_dep_on_X_G.m
- z_vector_unit_Euclidean_norm.m

Store all Matlab codes from Step 1 and Step 2 in the same directory.

Step 3: Run the source code "demo_Matlab_Figure_2.m" from Step 2.
