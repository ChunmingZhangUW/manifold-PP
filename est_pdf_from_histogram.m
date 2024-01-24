
function [bin_centers_vector, hat_f_vector] = ...
    est_pdf_from_histogram(X_vector, num_bins)

%--------------------------------------------------------------------------
% Function: density estimates from Histogram counts of X_vector
%--------------------------------------------------------------------------
% <Inputs>:
% X_vector: n*1 vector
% num_bins: the number B of bins, scalar
%--------------------------------------------------------------------------
% <Outputs>:
% bin_centers_vector: 1*B vector, 
% hat_f_vector:       1*B vector
%--------------------------------------------------------------------------

num_obs = length(X_vector); % n

[bin_counts_vector, bin_edges_vector] = histcounts(X_vector, num_bins);
% bin counts (n_1,...,n_B); bin edges (x_0, x_1, ..., x_B).
bin_centers_vector = ...
    ( bin_edges_vector(2:(num_bins+1)) + bin_edges_vector(1:num_bins))/2;
% bin centers (c_1,...,c_B)

rel_freq_vector = bin_counts_vector' / num_obs;
% (f_1,...,f_B) = (n_1,...,n_B)/n

hat_f_vector = rel_freq_vector / (bin_centers_vector(2) - bin_centers_vector(1));
% estimated pdfs at bin centers (\hat f(c_1), ..., \hat f(c_B))
