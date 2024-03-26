%
% Name    : demo_Matlab_Figure_2.m
% Function: Simulation of projection persuit in BKN (2018) with
%           p >> n, K=1, for Theorem 1, and Theorem 2(i).
%           The target distribution G is plausible.
% Called  : mean_cov_mixture_Gaussians_D_dim.m,
%           X_matrix_sim.m,
%           mixture_Gaussians_D_dim_generator.m, z_vector_feasible_dep_on_X_G.m,
%           est_pdf_from_histogram.m, CDF_mixture_Gaussians_D_dim_func.m,
%           pdf_mixture_Gaussians_D_dim_func.m, KS_stat_1_sample.m,
%           z_vector_unit_Euclidean_norm.m, EDF_D_dim_func.m,
%           z_vector_critical_dep_on_X_G.m
%--------------------------------------------------------------------------

clear;
close all;

rng(12314357, 'twister');

centering_mean = 0; % 1 for mean-centered; 0 for not mean-centered; in z_crt

%=============== Illustrate Theorem 1, Theorem 2(i) =======================
%====================== Distribution F for Input Data =====================
choice_F = 0; % for data entries in the data matrix
% choice_F = input([' Input choice_F (0 for N(0,1); 1 for U(-3,3); \n',...
%     '   2 for centered Beta(1,2); 4 for multivariate skew-t; \n',...
%     '   5 for multivariate normal; \n',...
%     '   6 for multivariate mixture of 2 Gaussians) = ']);
% types of variables in data matrix X_matrix

%---------- Input parameters in the distribution F for data entries -------
if     choice_F == 0 % i.i.d. Gaussian N(mu, sigma) variables in the vector X
    input_parameters_in_F_for_X.mu    = 0;
    input_parameters_in_F_for_X.sigma = 1;
end % if choice_F == 0; 1; 2; 4; 5; 6

num_obs = 100;  % # of sampling units (data vetcors)
if     choice_F == 0
    dim_p =  550;
end
gamma = dim_p / num_obs;

num_sims = 100; % # of Monte-Carlo runs, for boxplots of |\hat G_z - G^*|_\infty

%====================== Target Distribution G =============================
choice_G = 1;  % mixture Gaussian
choice_CDF_in_target = 'G'; % target CDF = G

if     choice_G == 1 % mixture Gaussian
    num_K = 1; % for 1-dimensional projections

    prop_vector = [1/2; 1/2]; % prop. for 2 components in the Gaussian mixture

    mu = 2; sigma = 1/2;
    mu_matrix = [-mu  +mu]; sigma_matrix = [sigma  sigma];

    if strcmp(choice_CDF_in_target, 'G') == 1
        pdf_1D_in_target = @(x) pdf_mixture_Gaussians_D_dim_func...
            (prop_vector, mu_matrix, sigma_matrix, x); % pdf of target CDF = G

        CDF_1D_in_target = @(x) CDF_mixture_Gaussians_D_dim_func...
            (prop_vector, mu_matrix, sigma_matrix, x); % CDF of target CDF = G

        options_KS_randsearch.prop_vector  = prop_vector;
        options_KS_randsearch.mu_matrix    = mu_matrix;
        options_KS_randsearch.sigma_matrix = sigma_matrix;
    end

    [mu_G_vector, cov_G_matrix] = mean_cov_mixture_Gaussians_D_dim...
        (prop_vector, mu_matrix, sigma_matrix);
    disp(' mu_G_vector, cov_G_matrix are: ')
    disp( [mu_G_vector, cov_G_matrix] )

    mu_2_G = cov_G_matrix + mu_G_vector^2; % 2nd moment

    %----------------------------------------------------
    disp('===== For Theorem 1 =====');
    disp([' gamma = p/n = ', num2str(gamma), ...
        ', mu_2(G) = ', num2str(mu_2_G)])

    disp('===== For Theorem 2(i) =====');
end % if choice_G == 1 % mixture Gaussian

%====================== Projection Directions =============================
%--------------------------------------------------------------------------

num_randsearch_KS = 10; %4 * 200; % # of random search for z_KS
disp(' ');
choice_z_vector_KS = 0;

num_grids = 100; % # of grid points;

num_bins = 10; % # of bins in the empirical pdf estimator

choice_for_source_vector_a = 1;
%input(' Input choice_for_source_vector_a (1 for sample quantiles of G; 2 for i.i.d. samples from G = ');
% if num_K >= 2
%     choice_for_source_vector_a = 2;
% end
if choice_for_source_vector_a == 1; num_samplings = 10000; end

%======================== Simulations Studies =============================
D_n_fea     = zeros(num_sims, 1);
D_n_dep_1   = zeros(num_sims, 1);
D_n_indep_2 = zeros(num_sims, 1);
D_n_crt_11  = zeros(num_sims, 1);
D_n_crt_3   = zeros(num_sims, 1);

counts_generated_datasets = 0;
counts_failed_feas_directions = 0;

fig_name_common_1 = ['demo_Matlab_Figure_2',...
    '_G=',num2str(choice_G),'_F=',num2str(choice_F),...
    '_a=',num2str(choice_for_source_vector_a)];
tic
for sim = 1:num_sims
    disp([' sim = ', num2str(sim)])
    %--------------------- generate data matrix ---------------------------
    X_matrix = X_matrix_sim(choice_F, dim_p, num_obs, ...
        input_parameters_in_F_for_X);
    % dim_p * num_obs matrix

    counts_generated_datasets = counts_generated_datasets + 1;
    % # of generated datasets

    %--------------- get a_vector_G from G --------------------------------

    if choice_G == 1
        if strcmp(choice_CDF_in_target, 'G') == 1
            if     choice_for_source_vector_a == 1
                samples_G = mixture_Gaussians_D_dim_generator...
                    (prop_vector, mu_matrix, sigma_matrix, num_samplings);
                % 1*num_samplings, vector, random samples from G

                a_vector_G = prctile(samples_G, 100*(1:num_obs)'/(num_obs+1));
                % n*1, vector, sample quantiles of G
            end % if choice_for_source_vector_a == 1; 2
        end
    end

    %--------------- a data-dependent "feasible" direction ----------------
    if sim == 1;  opt_display = 1;  else;  opt_display = 0;  end
    [z_vector_data_fea, tag_norm_z_vector_0] = ...
        z_vector_feasible_dep_on_X_G...
        (X_matrix, a_vector_G, opt_display);  % p*1, vector
    if tag_norm_z_vector_0 == 1 % X fails for BKN condition
        disp('  initial feasible direction z_fea fails!')
    end

    while tag_norm_z_vector_0 == 1 % X fails for BKN condition
        counts_failed_feas_directions = counts_failed_feas_directions + 1;
        % # of failed feasible directions
        %_________________ re-generate data matrix ________________________

        %-------------------- generate data matrix ------------------------
        X_matrix = X_matrix_sim(choice_F, dim_p, num_obs, ...
            input_parameters_in_F_for_X);
        % dim_p * num_obs matrix

        counts_generated_datasets = counts_generated_datasets + 1;
        % # of generated datasets

        %--------------- a data-dependent direction -------------------
        [z_vector_data_fea, tag_norm_z_vector_0] = ...
            z_vector_feasible_dep_on_X_G...
            (X_matrix, a_vector_G, opt_display);  % p*1, vector

    end % while tag_norm_z_vector_0 == 1 % X fails for BKN condition

    %--------------- another data-dependent direction ---------------------
    z_vector_data_dep_1 = z_vector_unit_Euclidean_norm( ...
        median(X_matrix(:, (1:num_obs)), 2) );
    % p*1 vector, on the unit sphere

    %-------------- a data-independent direction --------------------------
    z_vector_data_indep_2 = z_vector_unit_Euclidean_norm...
        ( normrnd(0, 1, dim_p, 1) ); % N(0, 1)
    % p*1 vector, uniformly distributed on the unit sphere

    %--------------- a data-dependent "critical" direction ----------------

    for crt_algorithm = 1:2 % 1 for Lagrangian method; 2 for Stiefel manifold
        if crt_algorithm == 1
            %_____________________________ Method = 11 ____________________________
            method_for_z_crt_11 = 11;
            %input(' Input method (0; 11; 12; 13; 2; 3) for z_crt = ');
            options_method_quad_min_over_sphere.method = method_for_z_crt_11;

            own_options_Manopt = [];

            if  options_method_quad_min_over_sphere.method == 11 % for method = 11
                % search for LGMP in the interval
                options_method_quad_min_over_sphere.num_sims            = [];
                options_method_quad_min_over_sphere.num_sims_LGMP_Hager = [];
                options_method_quad_min_over_sphere.x_vector_initial    = [];
            end

            [z_vector_crt_11, ~, S_vector_crt_11] = z_vector_critical_dep_on_X_G...
                (X_matrix, a_vector_G, ...
                options_method_quad_min_over_sphere, own_options_Manopt, ...
                centering_mean);
            % p*1 vector, on the unit sphere

            %__________________________________________________________________

        elseif crt_algorithm == 2
            %_____________________________ Method = 3 ____________________________
            method_for_z_crt_3 = 3;
            %input(' Input method (0; 11; 12; 13; 2; 3) for z_crt = ');
            options_method_quad_min_over_sphere.method = method_for_z_crt_3;

            if method_for_z_crt_3 == 3
                own_options_Manopt.verbose_Manopt = -1;      % hiding output on screen
                own_options_Manopt.choice_problem_ehess = 1; % including ehess Hessian
                own_options_Manopt.z_0 = []; % initial value for the solver
            end

            if options_method_quad_min_over_sphere.method == 3 % for method = 3
                % "Manopt" toolbox
                options_method_quad_min_over_sphere.num_sims            = [];
                options_method_quad_min_over_sphere.num_sims_LGMP_Hager = [];
                options_method_quad_min_over_sphere.x_vector_initial    = [];
            end

            [z_vector_crt_3, ~, S_vector_crt_3] = z_vector_critical_dep_on_X_G...
                (X_matrix, a_vector_G, ...
                options_method_quad_min_over_sphere, own_options_Manopt, ...
                centering_mean);
            % p*1 vector, on the unit sphere

            %__________________________________________________________________
        end
    end  % for crt_algorithm = 1:2

    %--------------- a data-dependent KS-based direction ------------------
    if choice_z_vector_KS == 0
        choice_min_max = 0;
        num_randsearch_KS = [];
    end

    %------------------ vector S ------------------------------------------

    S_vector_fea     = (z_vector_data_fea'     * X_matrix)';
    S_vector_dep_1   = (z_vector_data_dep_1'   * X_matrix)';
    S_vector_indep_2 = (z_vector_data_indep_2' * X_matrix)';
    %S_vector_crt_11 = (z_vector_crt_11'       * X_matrix)';
    %S_vector_crt_3  = (z_vector_crt_3'        * X_matrix)';
    % (S_1,...,S_n)', n*1 vector

    %========== Part 1: plots of pdf and CDF via 1 simulated data =========
    if sim == 1
        %--------------- 1D plot -----------------------------------------
        x_grid = linspace(-4, 4, num_grids)';

        pdf_1D_in_target_grid = zeros(num_grids, 1);
        CDF_1D_in_target_grid = zeros(num_grids, 1);
        EDF_grid_data_fea     = zeros(num_grids, 1);
        EDF_grid_data_dep_1   = zeros(num_grids, 1);
        EDF_grid_data_indep_2 = zeros(num_grids, 1);
        EDF_grid_crt_11       = zeros(num_grids, 1);
        EDF_grid_crt_3        = zeros(num_grids, 1);
        for g = 1:num_grids
            pdf_1D_in_target_grid(g) = pdf_1D_in_target(x_grid(g));
            CDF_1D_in_target_grid(g) = CDF_1D_in_target(x_grid(g));

            EDF_grid_data_fea(g)   = EDF_D_dim_func(...
                S_vector_fea,   x_grid(g));
            EDF_grid_data_dep_1(g)   = EDF_D_dim_func(...
                S_vector_dep_1,   x_grid(g));
            EDF_grid_data_indep_2(g) = EDF_D_dim_func(...
                S_vector_indep_2, x_grid(g));
            EDF_grid_crt_11(g)       = EDF_D_dim_func(...
                S_vector_crt_11,  x_grid(g));
            EDF_grid_crt_3(g)        = EDF_D_dim_func(...
                S_vector_crt_3,   x_grid(g));
        end

        [x_data_fea,   Epdf_data_fea] = ...
            est_pdf_from_histogram(S_vector_fea,   num_bins);
        [x_data_dep_1,   Epdf_data_dep_1] = ...
            est_pdf_from_histogram(S_vector_dep_1,   num_bins);
        [x_data_indep_2, Epdf_data_indep_2] = ...
            est_pdf_from_histogram(S_vector_indep_2, num_bins);
        [x_data_crt_11,  Epdf_data_crt_11] = ...
            est_pdf_from_histogram(S_vector_crt_11,  num_bins);
        [x_data_crt_3,   Epdf_data_crt_3] = ...
            est_pdf_from_histogram(S_vector_crt_3,   num_bins);

        %----------------------- compare pdfs ---------
        h_1 = figure(1);
        subplot(2,2,1)
        plot(x_grid,         pdf_1D_in_target_grid, 'k-');
        hold on;
        plot(x_data_fea,     Epdf_data_fea,   'b:', 'LineWidth', 1.5);
        plot(x_data_dep_1,   Epdf_data_dep_1,   'r--', 'LineWidth', 1.0);
        plot(x_data_indep_2, Epdf_data_indep_2, 'magenta-.');
        plot(x_data_crt_11,  Epdf_data_crt_11,  'cyan-o', 'LineWidth', 0.8);
        plot(x_data_crt_3,   Epdf_data_crt_3,   'k+', 'LineWidth', 0.8);

        x_min_max = xlim; % get the x-axis limits for the current axes.

        xlabel('{\boldmath{$s$}}', 'interpreter', 'latex');
        ylabel('{\textbf{compare pdf}}', 'interpreter', 'latex');
        title(['{\boldmath$n = ', num2str(num_obs), ...
            ',  p = ', num2str(dim_p), '$}'], 'interpreter', 'latex')
        figure_name_1 = [fig_name_common_1,...
            '_Epdf_B=', num2str(num_bins), ...
            '_n=',num2str(num_obs),'_p=',num2str(dim_p),...
            '_crt=',num2str(method_for_z_crt_11),'-',num2str(method_for_z_crt_3),...
            '_zKS=',num2str(choice_z_vector_KS),num2str(choice_min_max),...
            '_CeM=',num2str(centering_mean)];
        print(h_1, '-depsc', [figure_name_1, '.eps']);

        %----------------------- compare KDEs ---------

        %----------------------- compare EDFs with true CDF ----------
    end % if sim == 1

    %-------------------- KS statistics ---------------------------------

    vector_G_S_fea     = zeros(num_obs, 1); % (G(S_1),...,G(S_n))', n*1 vector
    vector_G_S_dep_1   = zeros(num_obs, 1); % (G(S_1),...,G(S_n))', n*1 vector
    vector_G_S_indep_2 = zeros(num_obs, 1); % (G(S_1),...,G(S_n))', n*1 vector
    vector_G_S_crt_11  = zeros(num_obs, 1); % (G(S_1),...,G(S_n))', n*1 vector
    vector_G_S_crt_3   = zeros(num_obs, 1); % (G(S_1),...,G(S_n))', n*1 vector
    if     choice_G == 1 % 1-dim Gaussian mixture  of 2-components
        for i = 1:num_obs
            vector_G_S_fea(i)     = CDF_1D_in_target(S_vector_fea(i));
            vector_G_S_dep_1(i)   = CDF_1D_in_target(S_vector_dep_1(i));
            vector_G_S_indep_2(i) = CDF_1D_in_target(S_vector_indep_2(i));
            vector_G_S_crt_11(i)  = CDF_1D_in_target(S_vector_crt_11(i));
            vector_G_S_crt_3(i)   = CDF_1D_in_target(S_vector_crt_3(i));
        end
    end

    D_n_fea(sim)     = KS_stat_1_sample(S_vector_fea,     vector_G_S_fea);
    D_n_dep_1(sim)   = KS_stat_1_sample(S_vector_dep_1,   vector_G_S_dep_1);
    D_n_indep_2(sim) = KS_stat_1_sample(S_vector_indep_2, vector_G_S_indep_2);
    D_n_crt_11(sim)  = KS_stat_1_sample(S_vector_crt_11,  vector_G_S_crt_11);
    D_n_crt_3(sim)   = KS_stat_1_sample(S_vector_crt_3,   vector_G_S_crt_3);
end
disp('')
disp('======== Part 2: compare box plots of KS-statistics ======')
if choice_z_vector_KS == 0
    disp([...
        mean(D_n_fea)    ...
        mean(D_n_dep_1)  mean(D_n_indep_2) ...
        mean(D_n_crt_11) mean(D_n_crt_3)])
end

toc

disp(['feasible directions (total, failures, prop) = ', ...
    num2str([counts_generated_datasets,  counts_failed_feas_directions, ...
    counts_failed_feas_directions/counts_generated_datasets])])

h_3 = figure(3); %----- compare box plots of KS-statistics ------
subplot(2,2,1)
if choice_z_vector_KS == 0
    boxplot([...
        D_n_fea,    ...
        D_n_dep_1,  D_n_indep_2, ...
        D_n_crt_11, D_n_crt_3])
    set(gca, 'XTickLabel', {'v_fea', ...
        'v_1', 'v_2', ...
        'v_crt; I', 'v_crt; II'}, 'XTickLabelRotation', 45); %, 'FontSize', 8);
end

ylabel(['\textbf{compare} ', ...
    '{\boldmath{$\|\widehat{G}_{z}-G\|_{\infty}$}}'], 'interpreter', 'latex')
title(['{\boldmath$n = ', num2str(num_obs), ...
    ',  p = ', num2str(dim_p), '$}'], 'interpreter', 'latex');

if choice_z_vector_KS == 0
    max_D_n = max([...
        max(D_n_fea),    ...
        max(D_n_dep_1),  max(D_n_indep_2), ...
        max(D_n_crt_11), max(D_n_crt_3)]);
end
ylim([ 0-0.025 min((max(0.55, max_D_n) + 0.05), 1) ])

figure_name_3 = [fig_name_common_1,...
    '_KS',...
    '_sims=',num2str(num_sims), ...
    '_n=',num2str(num_obs),'_p=',num2str(dim_p),...
    '_crt=',num2str(method_for_z_crt_11),'-',num2str(method_for_z_crt_3),...
    '_zKS=',num2str(choice_z_vector_KS),num2str(choice_min_max),...
    '_CeM=',num2str(centering_mean)];
print(h_3, '-depsc', [figure_name_3, '.eps']);

%============ For Section 4.1 =============================================
% demo_Matlab_Figure_2
%  mu_G_vector, cov_G_matrix are:
%          0    4.2500
%
% ===== For Theorem 1 =====
%  gamma = p/n = 5.5, mu_2(G) = 4.25
% ===== For Theorem 2(i) =====
%
% ======== Part 2: compare box plots of KS-statistics ======
%     0.0144    0.4970    0.3957    0.0144    0.0144
%
% Elapsed time is 19.889336 seconds.
% feasible directions (total, failures, prop) = 116            16      0.137931
%
