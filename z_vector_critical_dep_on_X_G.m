
function [z_vector_crt, hat_LGMP, S_vector_crt, info_solver_Manopt] = ...
    z_vector_critical_dep_on_X_G(X_matrix_orig, a_vector_G, ...
    options_method_quad_min_over_sphere, own_options_Manopt, centering_mean)

%--------------------------------------------------------------------------
% Function: compute the vector of "critical" direction in the
%           feasibility results in BKN (2018), where p >= n
% Called  : Least_squares_min_w_quad_equa_constraint.m
%--------------------------------------------------------------------------
% <Inputs>:
% X_matrix_orig: original data matrix, p*n; X'X is n*n, we assume p >= n,
%           so X'X is invertible, if X is not mean-centered.
% a_vector_G: n*1 vector,
% options_method_quad_min_over_sphere:
% own_options_Manopt: input-information of the solver from "Manopt" toolbox
% centering_mean: 1 for centering the mean; 0 for not centering the mean
%--------------------------------------------------------------------------
% <Outputs>:
% z_vector_crt: p*1 vector, requiring \|z_vector\|_2 = 1
%   hat_LGMP  : selected Lagrange multiplier for method = 11, 12, 13;
%               [] for methods = 0 and 2
% S_vector_crt: n*1 vector, (S_1,...,S_n)'.
% info_solver_Manopt: output information of the solver from "Manopt" toolbox
%--------------------------------------------------------------------------

[~, num_obs] = size(X_matrix_orig); % X^o, matrix, p*n

%==================== Part A: Preprocessing ===============================
if     centering_mean == 0 % 0 for not centering the mean
    X_matrix_orig_ctrd = X_matrix_orig;
elseif centering_mean == 1 % 1 for centering the mean;
    vector_mean_X = mean(X_matrix_orig, 2); % vector, p*1

    X_matrix_orig_ctrd = X_matrix_orig ...
        - vector_mean_X * ones(num_obs, 1)'; % X_c^o, matrix, p*n
end

%==================== Part B: Extraction ==================================
[z_vector_crt, hat_LGMP, info_solver_Manopt] = ...
    Least_squares_min_w_quad_equa_constraint...
    (X_matrix_orig_ctrd', a_vector_G, 1, ...
    options_method_quad_min_over_sphere, own_options_Manopt);
% p*1 vector, on the unit sphere

%==================== Part C: Reconstruction ==============================
if     centering_mean == 0 % 0 for not centering the mean
    S_vector_crt = (z_vector_crt' * X_matrix_orig)';
elseif centering_mean == 1 % 1 for centering the mean;
    S_vector_crt = (z_vector_crt' * X_matrix_orig_ctrd)';
end  % (S_1,...,S_n)', n*1 vector

end
