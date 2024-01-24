
function [mu_vector, cov_matrix] = mean_cov_mixture_Gaussians_D_dim...
    (prop_vector, mu_matrix, sigma_matrix)

%-------------------------------------------------------------------------
% Function: compute the (population) mean vector, cov matrix of
%           a D-dim Gaussian mixture of K-components
%-------------------------------------------------------------------------
% <Input>:
% prop_vector : = (p_1,...,p_K)': K*1 vector,
%               prop. for K components in the Gaussian mixture
% mu_matrix   : = (mu_1_vector,...,mu_K_vector), D*K matrix,
%   where mu_1_vector = (mu_11,...,mu_D1)',...
% sigma_matrix: = (sigma_1_vector,...,sigma_K_vector): D*K matrix
%-------------------------------------------------------------------------
% <Output>:
% mu_vector : D*1, vector
% cov_matrix: D*D, matrix
%-------------------------------------------------------------------------

[D, ~] = size(mu_matrix);

mu_vector = zeros(D,1);
cov_matrix = zeros(D,D);
for d = 1:D
    %-----------------------------------------------------
    mu_vector(d) = sum( prop_vector' .* mu_matrix(d,:) );
    
    %-----------------------------------------------------
    var_d = sum( prop_vector' .* (mu_matrix(d,:).^2 + sigma_matrix(d,:).^2) ) ...
        - (mu_vector(d))^2;
    cov_matrix(d,d) = var_d;
end

for d = 1:D
    %-----------------------------------------------------
    for k = (d+1):D
        cov_matrix(k,d) = ...
            sum( prop_vector' .* mu_matrix(k,:) .* mu_matrix(d,:) ) ...
            - mu_vector(k) * mu_vector(d);
        cov_matrix(d, k) = cov_matrix(k,d);
    end
end