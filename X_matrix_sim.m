
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
end
