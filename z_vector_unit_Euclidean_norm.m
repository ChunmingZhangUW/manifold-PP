
function [z_vector] = z_vector_unit_Euclidean_norm( X_vector )

% Function: re-scale X_vector with unit Euclidean norm
%--------------------------------------------------------------------------
% <Input>:
% X_matrix: vector
%--------------------------------------------------------------------------
% <Output>:
% z_vector: vector, X_vector / |X_vector|_2
%--------------------------------------------------------------------------

z_vector = X_vector / sqrt(sum(X_vector.^2)); % fastest
%z_vector = X_vector / ( (sum(X_vector.^2))^(1/2) ); % faster
%z_vector =  X_vector / norm(X_vector, 2); % slower

% vector, on the unit sphere
