
function [X_matrix] = X_matrix_sim(choice_F, dim_p, num_obs, ...
    input_parameters_in_F_for_X)

%--------------------------------------------------------------------------
% Function: Generate a data matrix X = (X_{j,i}), consisting of n random
%           data vectors, where each data vector follows the
%           distribution F of a random vector X=(X_1,...,X_p)^T.
%--------------------------------------------------------------------------
% <Inputs>:
% choice_F: types of distribution of a random vector X,
%           0 for i.i.d. Gaussian N(mu, sigma) variables in the vector X;
%           1 for i.i.d. Uniform  U(A, B) variables in the vector X;
%           2 for i.i.d. Beta Betadist(A, B) variables in the vector X;
%           4 for multivariate skew-t distribution for the vector X;
%           5 for multivariate-Normal distribution for the vector X;
%           6 for multivariate mixture of 2 Gaussians for the vector X.
% dim_p   : dimension
% num_obs : # of data vectors,
% input_parameters_in_F_for_X: an array, include parameters in F
%--------------------------------------------------------------------------
% <Outputs>:
% X_matrix:  % p*n, data matrix
%--------------------------------------------------------------------------

%--------------------- generate data matrix ---------------------------
if     choice_F == 0 % i.i.d. Gaussian N(mu, sigma) variables in the vector X
    mu    = input_parameters_in_F_for_X.mu;    % \mu
    sigma = input_parameters_in_F_for_X.sigma; % \sigma

    mu_F = mu; %var_F = sigma^2;
    X_matrix = normrnd(mu, sigma, dim_p, num_obs) - mu_F;
    % ensure E(vector_X) = 0,

elseif choice_F == 1 % i.i.d. Uniform  U(A, B) variables in the vector X
    A = input_parameters_in_F_for_X.A;
    B = input_parameters_in_F_for_X.B;

    mu_F = (A+B)/2; %var_F = (B-A)^2/12;
    X_matrix = unifrnd(A, B, dim_p, num_obs) - mu_F;
    % ensure E(vector_X) = 0,

elseif choice_F == 2 % i.i.d. Beta Betadist(A, B) variables in the vector X
    A = input_parameters_in_F_for_X.A;
    B = input_parameters_in_F_for_X.B;

    mu_F = A/(A+B); %var_F = A*B/((A+B)^2 * (A+B+1));
    X_matrix = 6 * (betarnd(A, B, dim_p, num_obs) - mu_F);
    % ensure E(vector_X) = 0,

elseif choice_F == 4 % multivariate skew-t distribution for the vector X
    df_t = input_parameters_in_F_for_X.df_t;
    rho  = input_parameters_in_F_for_X.rho;
    gamma_1 = input_parameters_in_F_for_X.gamma_1;

    vector_gamma = gamma_1 * ones(dim_p, 1); % (\gamma_1,...,\gamma_p)'

    eye_d = eye(dim_p, dim_p);

    corr_matrix = rho * ones(dim_p, dim_p) + (1-rho) * eye_d;
    % dim_p * dim_p, correlation matrix R
    A_matrix = corr_matrix ^(1/2); % for matrix square-root
    % dim_p * dim_p, A is symmetric and A * A^T = R

    v    = df_t;
    mu_F_vector  = v/(v-2) * vector_gamma;
    %mu_F_centered = 0;

    type_depend = input_parameters_in_F_for_X.type_depend;
    if     type_depend == 1 % for unit variance & correlated
        cov_F_matrix = v/(v-2) * diag(diag(corr_matrix)) ...
            + 2*v^2 / ( (v-2)^2*(v-4) ) ...
            * diag(diag(vector_gamma * vector_gamma'));
    elseif type_depend == 0 % for unit variance & un-correlated
        cov_F_matrix = v/(v-2) * corr_matrix ...
            + 2*v^2 / ( (v-2)^2*(v-4) ) ...
            * (vector_gamma * vector_gamma');
    end
    cov_F_matrix_inv_half = cov_F_matrix^(-1/2);

    %---------------------
    W_vector_v = v ./ chi2rnd(v, 1, num_obs); % 1*n

    X_matrix = zeros(dim_p, num_obs);
    % dim_p * n, skewed t_d(v; R; (r_1, r_2, r_3))
    for i = 1:num_obs % for each column (each data vector)
        Z_i_vector = normrnd(0, 1, dim_p, 1);
        X_matrix(:, i) = vector_gamma * W_vector_v(i) ...
            + sqrt(W_vector_v(i)) * A_matrix * Z_i_vector;
        X_matrix(:, i) = ...
            cov_F_matrix_inv_half * (X_matrix(:, i) - mu_F_vector);
        % mu_F_centered = 0; unit variance
        % dim_p * n, t_d(v, 0, R, gamma_1, gamma_2)
    end
    % ensure E(vector_X) = 0,

elseif choice_F == 5 % multivariate-Normal distribution for the vector X
    rho = input_parameters_in_F_for_X.rho;

    eye_d = eye(dim_p, dim_p);

    corr_matrix = rho * ones(dim_p, dim_p) + (1-rho) * eye_d;
    % dim_p * dim_p, correlation matrix R
    A_matrix = corr_matrix^(1/2);
    % dim_p * dim_p, A is symmetric and A * A^T = R

    mu_F_vector = 0;
    %mu_F_centered = 0;

    type_depend = input_parameters_in_F_for_X.type_depend;
    if     type_depend == 1 % for unit variance & correlated
        cov_F_matrix = diag(diag(corr_matrix));
    elseif type_depend == 0 % for unit variance & un-correlated
        cov_F_matrix = corr_matrix;
    end
    cov_F_matrix_inv_half = cov_F_matrix^(-1/2);

    X_matrix = zeros(dim_p, num_obs);
    % dim_p * n, skewed t_d(v; R; (r_1, r_2, r_3))
    for i = 1:num_obs % for each column (each data vector)
        Z_i_vector = normrnd(0, 1, dim_p, 1);
        X_matrix(:, i) = A_matrix * Z_i_vector;
        X_matrix(:, i) = ...
            cov_F_matrix_inv_half * (X_matrix(:, i) - mu_F_vector);
    end
    % ensure E(vector_X) = 0,

elseif choice_F == 6 % multivariate mixture of 2 Gaussians for the vector X
    choice_prop_1 = input_parameters_in_F_for_X.choice_prop_1; % 2, used in F
    mu_1 = input_parameters_in_F_for_X.mu_1;

    if     choice_prop_1 == 1
        prop_1 = 1/2;
    elseif choice_prop_1 == 2
        prop_1 = 1/3;
    end
    mu_2 = -prop_1/(1-prop_1) * mu_1;

    X_matrix = mixture_Gaussians_D_dim_generator...
        ([prop_1; 1-prop_1], [mu_1*ones(dim_p, 1), mu_2*ones(dim_p, 1)], ...
        [ones(dim_p, 1), ones(dim_p, 1)], num_obs);
end
