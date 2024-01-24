
function [Z_matrix_crt, F_at_minimizer, info_solver_Manopt, options_Manopt] = ...
    Manopt_Z_crt_opt(X_matrix, A_matrix, own_options_Manopt)

%--------------------------------------------------------------------------
% Function: compute the optimizer with the "Manopt" toolbox：
%           Z_{opt} = \arg\min_{Z_matrix: Z_matrix^T Z_matrix = I_K}
%                     1/2 \|X_matrix^T Z_matrix - A_matrix\|_F^2.
%--------------------------------------------------------------------------
% <Inputs>:
% X_matrix: original data matrix, p*n
% A_matrix: n*K matrix of column vectors a_1,...,a_K,
% own_options_Manopt: my own options in "Manopt" toolbox
%--------------------------------------------------------------------------
% <Outputs>:
% Z_matrix_crt: = (z_vector_1,..., z_vector_K), p*K matrix,
%               implemented using the "Manopt" toolbox.
% F_at_minimizer: value of the cost function at the solution Z_matrix_crt
% info_solver_Manopt: output information of the solver from "Manopt" toolbox
% options_Manopt: options of the "Manopt" toolbox
%--------------------------------------------------------------------------

[dim_p, ~] = size(X_matrix); % p*n
num_K = size(A_matrix, 2); % the number of projection directions, e.g., 1, 2, 3.

%---------------------------- Part 1 --------------------------------------
% Create the problem structure.
% if     num_K == 1
%     manifold = spherefactory(dim_p); % faster
% elseif num_K >= 2
manifold = stiefelfactory(dim_p, num_K);
%end
problem.M = manifold;

%---------------------------- Part 2 --------------------------------------
% Define the problem cost function, its Euclidean gradient and Hessian.
% F(Z) = 1/2 \|X_matrix^T * Z - A_matrix\|_F^2: objective function
problem.cost = @(Z)  ...
    1/2 * norm(X_matrix' * Z - A_matrix, 'fro')^2;
% 1/2 * trace( (X_matrix' * Z - A_matrix) * (X_matrix' * Z - A_matrix)' );

problem.egrad = @(Z) ...
    X_matrix * (X_matrix' * Z - A_matrix);
% (X*X^T)*Z - X*A = X*(X^T*Z - A)
% notice the 'e' in 'egrad' for Euclidean

if     own_options_Manopt.choice_problem_ehess == 0
    %----- if problem.ehess in "Manopt" is NOT given
    warning('off', 'manopt:getHessian:approx')
elseif own_options_Manopt.choice_problem_ehess == 1
    %----- if problem.ehess in "Manopt" is YES given
    problem.ehess = @(Z, zdot) ...
        (X_matrix * X_matrix') * zdot;
end

%---------------------------- Part 3 --------------------------------------
% Numerically check gradient consistency (optional).
% checkgradient(problem);

%---------------------------- Part 4 --------------------------------------
% Solve
options_Manopt.verbosity = own_options_Manopt.verbose_Manopt;
% verbose_Manopt: Print detailed information. (A value < 0 means hiding
%      information in the optimization process, while a value > 0
%      will display detailed information.)

% options.subproblemsolver = @trs_tCG;
% options_Manopt.useRand = true;

z_0 = own_options_Manopt.z_0; % initial value for the solver
[Z_matrix_crt, F_at_minimizer, info_solver_Manopt, options_Manopt] = ...
    trustregions(problem, z_0, options_Manopt);

end
