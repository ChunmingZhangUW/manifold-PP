
function [z_vector_feasible, tag_norm_z_vector_0, S_vector_feasible] = ...
    z_vector_feasible_dep_on_X_G(X_matrix_orig, a_vector_G, opt_display)

%--------------------------------------------------------------------------
% Function: compute the 'feasible' direction vector (in the "unit sphere")
%           constructed for plausible target distributions G in BKN (2018),
%           where p >= n is applicable under a condition on an initial z_0.
% Called  : z_vector_unit_Euclidean_norm.m
%--------------------------------------------------------------------------
% <Inputs>:
% X_matrix_orig: original data matrix, p*n; X'X is n*n, we assume p >= n,
%           so X'X is invertible, if X is not mean-centered.
% a_vector_G: n*1 vector,
% opt_display: 0 for "without" display; 1 for "with" display
%--------------------------------------------------------------------------
% <Outputs>:
% z_vector_feasible: p*1 vector, requiring \|z_vector_feasible\|_2 = 1
% tag_norm_z_vector_0: (norm_squared_z_vector_0 > 1);
%                      1 for fail: 0 for OK
% S_vector_feasible: n*1 vector, (S_1,...,S_n)'.
%--------------------------------------------------------------------------

z_vector_feasible   = [];
tag_norm_z_vector_0 = [];
S_vector_feasible   = [];

[dim_p, num_obs] = size(X_matrix_orig); % p*n

if     dim_p < num_obs  % p < n
    disp(['   dim_p < num_obs;' ...
        ' ''z_vector_feasible_dep_on_X_G'' exits and returns!!!']);
    % (X^T X) is not full-rank, so is not invertible.
    % The initial vector z_0 fails to apply.
    tag_norm_z_vector_0 = 1; % 1 for fail: 0 for OK

elseif dim_p >= num_obs % p >= n
    %==================== Part A: Preprocessing ===============================
    % Note: The data matrix X will NOT be mean-cnetered!!!

    %==================== Part B: Extraction ==================================
    %---------------------- obtain an initial vector z_vector_0 ---------------
    z_vector_0 = X_matrix_orig / (X_matrix_orig' * X_matrix_orig) * a_vector_G;
    if rank(X_matrix_orig) < num_obs % rank deficient
        disp([' Warning: rank(X_matrix_orig) < n, so ...' ...
            'X_matrix_orig'' * X_matrix_orig is not invertible!'])
    end
    norm_squared_z_vector_0 = sum(z_vector_0.^2); % \|z_vector_0\|_2^2

    tag_norm_z_vector_0 = (norm_squared_z_vector_0 > 1); % 1 for fail: 0 for OK

    if opt_display == 1
        disp([' \|z_vector_0\|_2^2 = ', num2str(norm_squared_z_vector_0)])
    end

    %---------------------- modify z_vector_0 ---------------------------------
    if     dim_p == num_obs % p = n
        if norm_squared_z_vector_0 == 1 % \|z_0\| = 1
            z_vector_feasible = z_vector_0;
            % p*1 vector, \|z_vector_feasible\|_2 = 1
        else
            disp(' norm_squared_z_vector_0 \ne 1; return!!!')
            return
        end

    elseif dim_p > num_obs % p >= n+1
        [~, ~, V_matrix] = svd(X_matrix_orig');
        %     [U,S,V] = svd(X) produces a diagonal matrix S, of the same
        %     dimension as X and with nonnegative diagonal elements in
        %     decreasing order, and unitary matrices U and V so that X = U*S*V'.
        % So, here, X_matrix_orig' = U S V^T, where X_matrix_orig' is n*p,
        % U is n*n, S is n*p, V is p*p.

        z_vector_feasible = z_vector_0 + ...
            sqrt( 1 - norm_squared_z_vector_0 ) * V_matrix(:, (num_obs+1));
        % p*1 vector, \|z_vector_feasible\|_2 = 1

        z_vector_feasible = z_vector_unit_Euclidean_norm(z_vector_feasible);
        % ensure \|z_vector_feasible\|_2 = 1
    end

    %==================== Part C: Reconstruction ==============================
    S_vector_feasible = (z_vector_feasible' * X_matrix_orig)';
    % (S_1,...,S_n)', n*1 vector
end % if dim_p < num_obs, dim_p >= num_obs

end
