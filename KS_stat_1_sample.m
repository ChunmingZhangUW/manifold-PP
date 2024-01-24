
function [D_n] = KS_stat_1_sample(vector_X, vector_G_X)

%--------------------------------------------------------------------------
% Funcution: compute Kolmogorov-Smirnov one-sample goodness-of-fit statistics,
% D_n = sup_{x \in R^1} |EDF(x)-G(x)|, where
%    EDF(x) is the E.D.F. of (X_1,...,X_n), and $G$ is a specified CDF.
%--------------------------------------------------------------------------
% <Inputs>:
% vector_X  : n*1 vector, (X_1,...,X_n)'
% vector_G_X: n*1 vector, (G(X_1),..., G(X_n))'
%--------------------------------------------------------------------------
% <Outputs>:
% D_n: scalar
%--------------------------------------------------------------------------

n_obs = length(vector_X);

vector_ordered_G_X = sort(vector_G_X, 'ascend');
% sorted data, G(X_{(1)}) <= ... <= G(X_{(n)})

D_n_plus  = max( (1:n_obs)'/n_obs - vector_ordered_G_X );
D_n_minus = max( vector_ordered_G_X - ((1:n_obs)-1)'/n_obs );
D_n = max( [D_n_plus, D_n_minus, 0] );

end
