
function [x_vector_opt, hat_LGMP, info_solver_Manopt] = ...
    Least_squares_min_w_quad_equa_constraint...
    (A_1_matrix, b_1_vector, a_circle, ...
    options_method_quad_min_over_sphere, own_options_Manopt)

%--------------------------------------------------------------------------
% Function: compute the optimizer
%           x_{opt} = \arg\min_{x_vector: \|x_vector\|_2 = a_circle}
%                     \|A_1_matrix x_vector - b_1_vector\|_2.
% Usage   : the "critical" direction in projection pursuit (Zhang, etal, 2023)
% Called: z_vector_unit_Euclidean_norm.m,
%         Manopt_Z_crt_opt.m.
%--------------------------------------------------------------------------
% <Inputs>:
% A_1_matrix: matrix, n*p
% b_1_vector: vector, n*1
% a_circle: > 0, scalar constraint in |x_vector|_2 = a_circle, e.g., a_circle = 1
%options_method_quad_min_over_sphere.method = 11; % Lagrangian, interval, fzero
%options_method_quad_min_over_sphere.method =  3; % "Manopt" toolbox
%----------------------------------------------------------
%  11: Lagrangian method, search for LGMP in an initial interval to solve a root
%                    Func_1(mu)=0 (with reliable result);
%   3: calling the Matlab toolbox "Manopt".
%----------------------------------------------------------
% own_options_Manopt: my own options in "Manopt" toolbox
%--------------------------------------------------------------------------
% <Outputs>:
% x_vector_opt: solution vector, p*1, on the sphere with radius a_circle,
%   hat_LGMP  : selected Lagrange multiplier for method = 11;
% info_solver_Manopt: output information of the solver from "Manopt" toolbox
%--------------------------------------------------------------------------

method = options_method_quad_min_over_sphere.method;

%--------------------------------------------------------------------------
hat_LGMP = []; info_solver_Manopt = [];
a_circle = abs(a_circle); % a > 0

if method == 11
    % for SVD method in
    [~, dim_p] = size(A_1_matrix); % A_1_matrix, n*p

    %_____________________ Part 1: b vector _________________________
    b_vector = A_1_matrix' * b_1_vector; % p*1 vector b_vector

    %_____________________ Part 2: SVD of matrix A_1 ___________________
    %svd    Singular value decomposition.
    %[U,S,V] = svd(X) produces a diagonal matrix S, of the same
    %dimension as X and with nonnegative diagonal elements in
    %decreasing order, and unitary matrices U and V so that X = U*S*V'.
    [~, Sigma_matrix, V_matrix] = svd(A_1_matrix); % [U,S,V] = svd(A_1)
    % U: n*n; S: n*p; V: p*p.

    diag_Sigma_1_matrix = diag(Sigma_matrix); % r*1 vector, diag(\Sigma_1)
    r = length(diag_Sigma_1_matrix);
    % diag(sigma_r,...,sigma_1), sigma_j>=0, decreasing

    Phi_matrix = V_matrix; % p*p

    %----------------- eigenvalues of matrix A = A_1^T A_1 ----------------
    vector_lambda_values = zeros(dim_p, 1); % p*1 vector, initialize
    vector_lambda_values(1:r) = diag_Sigma_1_matrix.^2;
    % first r values, diag(\Sigma_1^2)

    lambda_min = vector_lambda_values(dim_p);
    %min(vector_lambda_values); % minimum \lambda_j(A_1^T A_1)

    %_____________________ Part 3: 2 index sets _________________________
    E_set_one  = find(vector_lambda_values == lambda_min); % E_1
    E_set_plus = find(vector_lambda_values  > lambda_min); % E_+, could be []

    %_____________________ Part 4: beta_vector _________________________
    beta_vector      = Phi_matrix' * b_vector;  % beta_vector, p*1
    norm_beta_vector = sqrt( sum(beta_vector .^ 2) ); % |beta_vector|

    beta_vector_1      = beta_vector(E_set_one);       % beta_vector_1
    norm_beta_vector_1 = sqrt( sum(beta_vector_1 .^ 2) ); % |beta_vector_1|

    %_____________________ Part 5: solution c_vector ______________________
    c_vector = zeros(dim_p, 1); % p*1 vector, initialize

    %----- 5.1: define c_+, a^2-\|c_+\|_2^2 to separate 2 cases below -----
    c_vector_plus = [];    % c_vector_+
    d_plus_2 = a_circle^2; % scalar

    if length(E_set_plus) >= 1
        c_vector_plus = beta_vector(E_set_plus) ./ ...
            (vector_lambda_values(E_set_plus) - lambda_min); % c_vector_+
        d_plus_2 = a_circle^2 - sum(c_vector_plus .^ 2); % scalar
    end

    %------------ 5.2: determine Lagrange multiplier -----------
    if (norm_beta_vector_1 == 0) && (d_plus_2 >= 0) % Degenerate case
        hat_LGMP = -lambda_min; % selected Lagrange multiplier

        length_E_set_one = length(E_set_one); % # E_1

        c_vector(E_set_plus) = c_vector_plus; % c_+
        c_vector(E_set_one)  = sqrt(d_plus_2) * eye(length_E_set_one, 1); % c_1
        %c_vector(E_set_one) = sqrt(d_plus_2/length_E_set_one) ...
        %    * ones(length_E_set_one, 1);

    else % Nondegenerate case
        %--------------- search for Lagrange multiplier -----------------
        mu_Lower = norm_beta_vector_1 / a_circle - lambda_min;
        mu_Upper = norm_beta_vector   / a_circle - lambda_min;

        if method == 11
            Func_1 = @(mu) a_circle^2 - ...
                sum( (beta_vector ./ (vector_lambda_values + mu)).^2 );
            % F_1(mu) = a_circle^2 - sum_{j=1}^p [beta_j/(lambda_j+mu)]^2

            if     method == 11 % search mu in an interval
                mu_init = [mu_Lower, mu_Upper];

            end % if method == 11,

            [hat_LGMP, ~, EXITFLAG] = fzero(@(mu) Func_1(mu), mu_init);
            if EXITFLAG ~= 1
                disp([' EXITFLAG of fzero using Lagrangian = ', ...
                    num2str(EXITFLAG)])
            end

        end % if method == 11,

        if mu_Lower <= hat_LGMP && hat_LGMP <= mu_Upper
        else
            disp('   <-- !!!hat_mu fails!!! outside the true interval -->')
        end

        %------------------ 5.3: compute c_vector -------------------------
        c_vector = beta_vector ./ (vector_lambda_values + hat_LGMP); % p*1
    end % if Degenerate case and non-degenerate case

    %_____________________ Part 6: solution x_vector ______________________
    x_vector_opt = Phi_matrix * c_vector; % p*1 solution vector

elseif method == 3 % "Manopt" toolbox
    [x_vector_opt, ~, info_solver_Manopt, ~] = ...
        Manopt_Z_crt_opt(A_1_matrix', b_1_vector, own_options_Manopt);

    hat_LGMP = 0;
end % method == 11; 3.

x_vector_opt = z_vector_unit_Euclidean_norm(x_vector_opt);
% ensure \|x_vector_opt\|_2 = 1
end

