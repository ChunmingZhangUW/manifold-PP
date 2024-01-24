
function EDF_func = EDF_D_dim_func(X_matrix, x_vector)

%--------------------------------------------------------------------------
% Function: compute the EDF of n random vectors in X_matrix
%--------------------------------------------------------------------------
% <Inputs>:
% X_matrix: D*n matrix, (X_matrix(:,1),...,X_matrix(:,n))
% x_vector: D*1 vector, (x_1,...,x_D)'
%--------------------------------------------------------------------------
% <Outputs>:
% EDF_func: scalar
%--------------------------------------------------------------------------

[~, num_obs] = size(X_matrix);

vector_ind = zeros(num_obs, 1); % n*1, vector
for i = 1:num_obs
    vector_ind(i) = prod(X_matrix(:,i) <= x_vector); % 0 or 1
end

EDF_func = mean( vector_ind ); % scalar
