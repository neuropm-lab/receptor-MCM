% =========================================================================
% Molecular MCM model fitting
% =========================================================================
%
% DESCRIPTION:
% Fit molecular-MCM models.
%
% =========================================================================
% INPUTS:
%   - data_for_MCM.mat: Path to MCM data structure from NeuroPM-box
%
% =========================================================================
% OUTPUTS:
%   - data_for_MCM.mat:
%
% =========================================================================
% AUTHOR: Ahmed Faraz Khan
% CONTACT: ahmed.faraz.khan@mai.mcgill.ca
% VERSION: 1.0.0
% =========================================================================
                                       
% =================Rename==========================================
ppmi_data,
% =================Rename==========================================
                        
cd('/export02/data/Work/re-MCM_PD')
addpath('/export02/data/Work/conn18b')
addpath(genpath('/export02/data/Work/MATLAB/'))
addpath(genpath('/export02/data/Work/MATLAB/LAR'))

addpath(genpath('/export02/data/Work/MATLAB/spm12'))
%addpath(genpath('/export02/data/Work/re-MCM_PD/nngarotte/'))
addpath(genpath('/export02/data/Work/lasso/'))
% Conflicts with stats.lasso
%addpath(genpath('/export02/data/Work/SPAMtoolbox/lar.m'))

load('output/PPMI_for_MCM.mat');
clear -regexp adni*
clear -regexp V*    
clear -regexp data_FUNC*
clear -regexp data_GM*
clear -regexp data_METB
%clear -regexp ppmi*
clear -regexp summ*
clear mod
close all

%DIAG_ID = 1;

isSnormalized = true;
isLongitudinal = true;
isEM = false;
isOneShot = true;

N_RAND_ITERS=1000;
N_nz_param = 50;
N_nz_rec = 3;
% Use 50th percentile of healthy subjects for baseline template
HEALTHY_PRCTILE = 50; 
rng(0);

%% 
% For PPMI paper: modelCodes 1,2,3,4,10
% Other model codes were experimental
for modelCode=0
    sprintf('Running model %d', modelCode)
    
    % Imputed regions flag
    % 1: use full NeuroPM atlas (Julich + Brodmann ROIs)
    % 0: use only Julich regions with receptor data
    % -1: use only imputed Brodmann regions from NeuroPM atlas
    useImputedROI = 1;
    if modelCode == 1
        % Only neuroimaging
        includeRec = false;
        includeInteract = false;
        isRandomZ = false;
    elseif modelCode == 2
        % Direct neuroimaging + receptor terms
        includeRec = true;
        includeInteract = false;
        isRandomZ = false;
    elseif modelCode == 3
        % Full model
        includeRec = true;
        includeInteract = true;
        isRandomZ = false;
    elseif modelCode == 4
        % Permutation analysis
        mkdir('output/Random_Maps');
        includeRec = true;
        includeInteract = true;
        isRandomZ = true;
    elseif modelCode == 5
        % Julich data regions only (no imputed receptor data)
        useImputedROI = 0;
        includeRec = true;
        includeInteract = true;
        isRandomZ = false;
    elseif modelCode == 6
        % Imputed, non-Julich data regions only
        useImputedROI = -1;
        includeRec = true;
        includeInteract = true;
        isRandomZ = false;
%     elseif modelCode == 7
%         % Receptor residuals: exclude one receptor at a time
%         includeRec = true;
%         includeInteract = true;
%         isRandomZ = false;
%     elseif modelCode == 8
%         % Cross-sectional model comparison
%         includeRec = true;
%         includeInteract = true;
%         isRandomZ = false;
%     elseif modelCode == 9
%         % Cross-sectional model comparison (null models)
%         includeRec = true;
%         includeInteract = true;
%         isRandomZ = false;
    elseif modelCode == 10
        % Personalized model, residuals with and without receptors and null
        % Receptor residuals
        includeRec = true;
        includeInteract = true;
        isRandomZ = false;
        N_iters = 1000;
        useImputedROI = 0;
        null_res_savedir = 'output/modelCode=10_null_residuals';
        mkdir(null_res_savedir)
        for fac=1:N_facs
            mkdir(sprintf('%s/fac_%d', null_res_savedir, fac));
        end
        for rec=1:N_recs
            for fac=1:N_facs
                mkdir(sprintf('%s/fac_%d/rec_%d', null_res_savedir, fac, rec));
            end
        end
    elseif modelCode == 11
        % Personalized model, residuals with and without receptors and null
        % Imaging residuals
        includeRec = false;
        includeInteract = false;
        isRandomZ = false;
        N_iters = 1000;
        useImputedROI = 0;
        null_res_savedir = 'output/modelCode=11_null_residuals';
        mkdir(null_res_savedir)
        for fac=1:N_facs
            mkdir(sprintf('%s/fac_%d', null_res_savedir, fac));
        end
        for rec=1:N_recs
            for fac=1:N_facs
                mkdir(sprintf('%s/fac_%d/rec_%d', null_res_savedir, fac, rec));
            end
        end
    end
    
    if ~isRandomZ
       N_map_iters=1; 
    else
       N_map_iters = N_RAND_ITERS;
    end

    r2s_random_maps = [];
    bs_random_maps = [];
    rss_random_maps = [];
    group_bs_random_maps = [];
    mcm_num_subj_samples = [];

    clear res_corrs % Correlation of residuals across time points

    for iter=1:N_map_iters
        if isRandomZ 
            Z0 = zeros(size(Z0_original,1), size(Z0_original,2));
            if ~mod(iter,100) 
                fprintf("Random Map iteration %d \n", iter)
            end
            for rec=1:size(Z0_original,1)
                Z0(rec,:) = Z0_original(rec,randperm(size(Z0,2)));
            end
        else 
            Z0 = Z0_original;
        end
        r2s_all_diags = [];
        rss_all_diags = [];
        ps_all_diags = [];
        r2s_sig_all_diags = [];
        rss_sig_all_diags = [];
        bs_all_diags = [];
        res_all_diags = [];
        group_bs_all_diags = [];
        residual_model_comp_all_diags = []; 

        if useImputedROI == 0
            ind_regs_data = find(regions_have_data);
        elseif useImputedROI == 1
            ind_regs_data = 1:N_regs;
        elseif useImputedROI == -1
            ind_regs_data = find(regions_have_data == 0);
        end

        N_regs_data = numel(ind_regs_data);

        for DIAG_ID=3%1:numel(diag_names)
            % Samples with diagnosis
            ind_diag_samples = find(long_diags == DIAG_ID);
            
            if ind_diag_samples == 0
                continue
            end
            diag_X0 = subjects_X0_adjusted(:,:, ind_diag_samples);
            diag_ages = long_ages(ind_diag_samples);

            % Subjects with diagnosis
            ind_diag_subj = find(subjects_diags == DIAG_ID);
            diag_num_subj_samples = num_subj_samples(ind_diag_subj);
            diag_samples_subj_ids = samples_subj_ids(ind_diag_samples);
            
            % Background/demographic data 
            diag_bg = subjects_background(ind_diag_subj,:);

            % Sort subjects by age
            %[sorted_ages, sorted_diag_ind_X0] = sort(diag_ages);
            %sorted_diag_X0 = diag_X0(:, sorted_diag_ind_X0);
            N_diag_samples = numel(ind_diag_samples);
            %sorted_diag_background = diag_backgrounds(sorted_diag_ind_X0, :);

            % Regression for receptor-factor interactions + factor diffusion

            if ~isLongitudinal

                % Treat 0s as missing data
                if isSnormalized
                    unnormed_X = reshape(sorted_diag_X0, N_regs, N_facs, numel(ind_diag_subjects));

                    % Normalize dims together
                    % [regions, factors, subjects] -> [factors, regions * subjects]
                    X_flat_norm = normalize(reshape(permute(unnormed_X, [2 1 3]), N_facs, N_regs * N_diag_samples), 2);
                    % [factors, regions, subjects] ([regions, factors, subjects])
                    diag_X_reg = permute(reshape(X_flat_norm, N_facs, N_regs, N_diag_samples), [2 1 3]);

                else
                    diag_X_reg = reshape(sorted_diag_X0, N_regs, N_facs, numel(ind_diag_subjects));
                end

                dS_dt = zeros(N_diag_samples - 1, N_regs, N_facs);
                dts = zeros(N_diag_samples-1,1);
                sums = zeros(N_diag_samples -1);

                for i=1:N_diag_samples - 1
                    dt = sorted_ages(i+1) - sorted_ages(i);
                    dts(i, 1) = dt;
                    %if dt~=0
                    dS_dt(i,:,:) = (squeeze(diag_X_reg(:,:,i+1)) - squeeze(diag_X_reg(:,:,i))) /dt;

                end


            % Longitudinal data, do not sort
            else
                if isSnormalized
                    %unnormed_X = reshape(diag_X0, N_regs, N_facs, N_diag_samples);

                    % Normalize dims together
                    % [regions, factors, subjects] -> [factors, regions * subjects]
                    X_flat_norm = normalize(reshape(permute(diag_X0, [2 1 3]), N_facs, N_regs * N_diag_samples), 2);
                    % [factors, regions, subjects] ([regions, factors, subjects])
                    diag_X_reg = permute(reshape(X_flat_norm, N_facs, N_regs, N_diag_samples), [2 1 3]);

                    % Youngest healthy samples as a baseline (normalize wrt all healthy subjects)
                    ind_baseline = find(long_diags == DIAG_ID);%find(strcmp("Control", diag_names)));
                    baseline_ages = long_ages(ind_baseline);
                    ind_youngest_baseline = find(baseline_ages < prctile(baseline_ages,HEALTHY_PRCTILE));
                    baseline_X0 = subjects_X0_adjusted(:,:,ind_baseline);

                    baseline_flat_norm = normalize(reshape(permute(baseline_X0, [2 1 3]), N_facs, N_regs * numel(ind_baseline)), 2);
                    baseline_X_reg = permute(reshape(baseline_flat_norm, N_facs, N_regs, numel(ind_baseline)), [2 1 3]);

                    yh_X_reg = squeeze(mean(baseline_X_reg(:, :, ind_youngest_baseline),3));

                else
                    diag_X_reg = reshape(diag_X0, N_regs, N_facs, N_diag_samples);
                end

                % wrt. baseline
                delta_S = zeros(N_regs, N_facs, sum(diag_num_subj_samples) - numel(diag_num_subj_samples));
                % wrt. previous sample
                dS_dt = zeros(sum(diag_num_subj_samples) - numel(diag_num_subj_samples), N_regs, N_facs);

                counter = 1; % all samples counter
                d_counter = 1; % baseline-subtracted samples counter
                for subj=1:numel(diag_num_subj_samples)
                    
                    for sample=1:diag_num_subj_samples(subj)
                        if sample > 1 % For the same subject after the baseline
                            dt = diag_ages(counter) - diag_ages(counter-1);
                            if dt == 0
                                sprintf("Warning: dt = 0")
                            end
                            dS_dt(d_counter,:,:) = (squeeze(diag_X_reg(:,:,counter)) - squeeze(diag_X_reg(:,:,counter-1))) /dt ;
                            delta_S(:,:,d_counter) = squeeze(diag_X_reg(:,:,counter)) - yh_X_reg;

                            d_counter = d_counter + 1;
                        end
                        counter = counter + 1;
                    end
                end

            end

            
            ind_remove = [];
            for i=1:numel(unique(diag_samples_subj_ids)) % For each subject with this diagnosis
                unique_ids = unique(diag_samples_subj_ids);
                subject = unique_ids(i); % Find subject id
                inds = find(diag_samples_subj_ids == subject); % Find sample indices for subject
                ind_remove = [ind_remove; inds(1)]; % Store the first sample index
                if numel(inds) < 3

                    sprintf("Warning: Subject %d has only %d samples", i, numel(inds))
                end
            end
            d_diag_samples_subj_ids = diag_samples_subj_ids;
            d_diag_samples_subj_ids(ind_remove) = [];
            d_diag_num_subj_samples = diag_num_subj_samples - 1; %ones(numel(diag_num_subj_samples), 1);
            mcm_num_subj_samples = cat(1, mcm_num_subj_samples, d_diag_num_subj_samples);

            % Normalize receptors
            Z0_norm = normalize(Z0, 2);

            %% Calculate regional spreading

            ks = 1e-5;%5e-3; for ridge regression

            % Receptor interactions, factor interactions, spreading, no inputs
            N_params = ((1 + N_recs) * N_facs ) + 1 + 1;

            % Factors, differential samples, regions
            net_spreads = zeros(N_facs,  size(dS_dt,1), N_regs);

            for reg=1:N_regs
                % fprintf("Region %d", reg);
                ind_other_regs = setdiff(1:N_regs, reg);
                conn_from = squeeze(sc(reg,ind_other_regs));
                conn_to = squeeze(sc(ind_other_regs, reg));

                for fac=1:N_facs
                    %y = squeeze(dS_dt(:, reg, fac));

                    % Spreading = S^sum( C_j->i S_j - C_i->j S_i )
                    %ind_other_facs = setdiff(1:N_facs, fac);
                    fac_in  = (conn_to .* squeeze(delta_S(ind_other_regs, fac, :)))'; % [subj, reg]
                    fac_out = (conn_from .*  squeeze(delta_S(reg, fac, :))); % [subj, reg]
                    net_spread = squeeze(sum(fac_in - fac_out, 2)); % Sum all regions
                    net_spreads(fac, :, reg) = net_spread;

                end
            end


            %% Regression on all patients, regions

            if isOneShot

                % By individual subject, find parameters that change

                %load('workspaces/before_combined_regress.mat');
                
                if modelCode == 1 || (modelCode==11)
                    N_params_combined = N_facs + 2;
                elseif modelCode == 2
                    N_params_combined =  N_facs + N_recs + 2;
                elseif (modelCode>=3) && (modelCode<=6) 
                    N_params_combined = (N_facs * N_recs) + N_facs + N_recs + 2;
                elseif (modelCode==7) || (modelCode==8) || (modelCode==9) || (modelCode==10) 
                    N_params_combined = (N_facs * N_recs) + N_facs + N_recs + 2;
                end
                
                method=1; % Linear regression
                
                %disp("Individual Models...")
                ks = 1e-5;%5e-3; for ridge regression

                delta_S(:, N_facs+1:N_facs+N_recs, :) = repmat(Z0_norm', 1,1,size(delta_S,3));
                delta_SZ_perm = permute(delta_S, [2 3 1]);

                % Skips subjects without at least 2 longitudinal samples
                unique_subject_ids = unique(d_diag_samples_subj_ids);

                bs_individual = zeros(numel(unique_subject_ids), N_facs, N_params_combined );
                r2s_individual = zeros(numel(unique_subject_ids), N_facs );
                RSS_individual = zeros(numel(unique_subject_ids), N_facs);
                ps_individual = zeros(numel(unique_subject_ids), N_facs, N_params_combined);
                r2s_sig_individual = zeros(numel(unique_subject_ids), N_facs );
                rss_sig_individual = zeros(numel(unique_subject_ids), N_facs );
                res_individual = zeros(numel(unique_subject_ids), N_facs, N_regs_data );
                
                % Mean squared residuals
                if modelCode == 10
                    msr_with = zeros(numel(unique_subject_ids), N_facs, N_recs, N_regs_data);
                    msr_without = zeros(numel(unique_subject_ids), N_facs, N_regs_data);
                    
                    r2_with = zeros(numel(unique_subject_ids), N_facs, N_recs);
                    r2_without = zeros(numel(unique_subject_ids), N_facs);
                    
                    res_with = zeros(numel(unique_subject_ids), N_facs, N_recs, N_regs_data * 3);
                    res_without = zeros(numel(unique_subject_ids), N_facs, N_regs_data * 3);
                    
                    leverage_with = zeros(numel(unique_subject_ids), N_facs, N_recs, N_regs_data * 3);
                    leverage_without = zeros(numel(unique_subject_ids), N_facs, N_regs_data * 3);
                    
                elseif modelCode == 11
                    msr_with = zeros(numel(unique_subject_ids), N_facs, N_facs, N_regs_data);
                    msr_without = zeros(numel(unique_subject_ids), N_facs, N_regs_data);
                    
                    r2_with = zeros(numel(unique_subject_ids), N_facs, N_facs);
                    r2_without = zeros(numel(unique_subject_ids), N_facs);
                    
                    res_with = zeros(numel(unique_subject_ids), N_facs, N_facs, N_regs_data * 3);
                    res_without = zeros(numel(unique_subject_ids), N_facs, N_regs_data * 3);
                    
                    leverage_with = zeros(numel(unique_subject_ids), N_facs, N_facs, N_regs_data * 3);
                    leverage_without = zeros(numel(unique_subject_ids), N_facs, N_regs_data * 3);
                    
                end
                
                    
                if modelCode ~= 8
                residual_model_comp = zeros(numel(unique_subject_ids), N_facs, N_regs_data, N_recs, 6);
                    
                for i=1:numel(unique_subject_ids)

                    subject = unique_subject_ids(i); % Find subject id
                    ind_subject = find(d_diag_samples_subj_ids == subject); % Find sample indices for subject
                    
                    % % [factors, regions, subjects] -> [factors, subjects, regions]
                    %delta_SZ = reshape(delta_SZ_perm(:,ind_subject,ind_regs_data), N_facs+N_recs,  numel(ind_subject) * N_regs_data);
                    
                    dS_dt_perm = permute(dS_dt(ind_subject,ind_regs_data,:), [3 1 2]);
                    % [N_facs, N_time, N_regs]
                    dS_dt_x = reshape(dS_dt_perm, N_facs, N_regs_data * numel(ind_subject));
                    subject_net_spreads = net_spreads(:, ind_subject, ind_regs_data);
                    
                    % Arrange neuroimaging data as [Timepoint1_reg123 tp2_reg123]
                    delta_SZ_tp = reshape(permute(delta_SZ_perm(:,ind_subject,ind_regs_data), [1 3 2]), ...
                        N_facs+N_recs,  numel(ind_subject) * N_regs_data);
                    dS_dt_perm_tp = permute(dS_dt(ind_subject,ind_regs_data,:), [3 2 1]);

                    %[b, r2, RSS, ps, r2_sig, rss_sig, cum_effect, res] = LassoPredict(dS_dt_perm, delta_SZ, subject_net_spreads, N_nz_param, N_regs_data, numel(ind_subject), includeRec, includeInteract, method, false);
                    [b, r2, RSS, ps, r2_sig, rss_sig, cum_effect, res] = LassoPredict(dS_dt_perm_tp, delta_SZ_tp, subject_net_spreads, N_nz_param, N_regs_data, numel(ind_subject), includeRec, includeInteract, method, false);
                    
                    % b = [method, fac, param]
                    bs_individual(i, :, :) = b;
                    r2s_individual(i, :) = r2;
                    RSS_individual(i,:) = RSS;
                    ps_individual(i,:,:) = ps;
                    r2s_sig_individual(i,:) = r2_sig;
                    rss_sig_individual(i,:) = rss_sig;
                    %cum_effect_individual(i,:,:,:) = cum_effect;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if modelCode == 10
                    
                    % Model performance with and without a specific receptor at a region
                    % Align timepoints and regions in 2D
                    % [fac, reg, timepoints]
                    
                    % Regression without receptors
                   
                    %inter_fac = squeeze(delta_SZ(1:N_facs, :))';
                    for fac=1:N_facs
                        %y = squeeze(reshape(dS_dt_perm(fac,:,:), 1, numel(ind_subject)*N_regs_data))';
                        y = dS_dt_perm_tp(fac,:)';
                        
                        net_spread = reshape(squeeze(subject_net_spreads(fac, :, :))', numel(y), 1);
                        X = [ones(numel(y),1) delta_SZ_tp(1:N_facs,:)' net_spread];
                        
                        X(:,2:end) = normalize(X(:,2:end));
                        y = normalize(y);
                        
                        leverage_without(i,fac,1:numel(y)) = leverage(X(:,2:end));
                        [~, ~, res, ~, stats] = regress(y,X);
                        res_without(i,fac,1:numel(res)) = res;
                        
                        r2_without(i,fac) = stats(1);
                        for reg=1:N_regs_data
                            %ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                            ind_curr_reg = [reg];
                            for d_ind=1:(numel(ind_subject)-1) % Add indices of all time points
                                ind_curr_reg = [ind_curr_reg reg+(N_regs_data*d_ind)];
                            end
                            msr_without(i,fac,reg) = sumsqr(res(ind_curr_reg))/numel(ind_subject) ;
                        end

                        
                    end
                    
                    for rec=1:N_recs
                        sprintf('[%d of %d] Model Code 10: Receptor %d', i, numel(unique_subject_ids), rec)

                        for fac=1:N_facs
                            sprintf('[%d of %d] Model Code 10: Receptor %d | Factor %d', i, numel(unique_subject_ids), rec, fac)
                            
                            sub_r2_nulls  = zeros(N_iters,1);
                            sub_msr_with_nulls = zeros(N_regs_data, N_iters);
                            sub_res_with_null = zeros(N_regs_data * 3, N_iters);
                            sub_leverage_with_null = zeros(N_regs_data * 3, N_iters);
                            
                            % Regression
                            %y = squeeze(dS_dt_perm(fac, :))';
                            y = dS_dt_perm_tp(fac,:)';
                            %inter_fac = squeeze(delta_SZ(1:N_facs, :))';
                            

                            rep_facs = repelem(delta_SZ_tp(1:N_facs,:), 1, 1);
                            rep_recs = repmat(delta_SZ_tp(N_facs+rec,:), N_facs, 1);
                            rec_fac = rep_facs .* rep_recs;
                            
                            % Net spreading (factor, subject, region)
                            net_spread = reshape(squeeze(subject_net_spreads(fac, :, :))', numel(y), 1);
                            
                            X = [ones(numel(y),1), rec_fac' delta_SZ_tp(1:N_facs,:)' delta_SZ_tp(N_facs+rec, :)' net_spread];
                            y = normalize(y);
                            X(:,2:end) = normalize(X(:,2:end));
                            
                            leverage_with(i,fac,rec,1:numel(y)) = leverage(X(:,2:end));
                            [b_temp, ~, res, ~, stats] = regress(y,X);
                            
                            res_with(i,fac,rec,1:numel(res)) = res;
                            
                            r2_with(i,fac,rec) = stats(1);
                            for reg=1:N_regs_data
                                %ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                            %ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                                ind_curr_reg = [reg];
                                for d_ind=1:(numel(ind_subject)-1) % Add indices of all time points
                                    ind_curr_reg = [ind_curr_reg reg+(N_regs_data*d_ind)];
                                end
                                msr_with(i,fac,rec,reg) = sumsqr(res(ind_curr_reg))/numel(ind_subject) ;
                            end
                            
                            % Random permutations
                            for r=1:N_iters
                                
                                for col=1:N_facs
                                    rec_fac(col,:) = rec_fac(col,randperm(size(rec_fac,2)));
                                end
                                
                                rec_shuffled = delta_SZ_tp(N_facs+rec, randperm(size(delta_SZ_tp,2)))';

                                X = [ones(numel(y),1), rec_fac' delta_SZ_tp(1:N_facs,:)' rec_shuffled net_spread];
                                %y = normalize(y);
                                X(:,2:end) = normalize(X(:,2:end));
                                [~, ~, res, ~, stats] = regress(y,X);
                                
                                sub_leverage_with_null(1:numel(y),r) = leverage(X(:,2:end));
                                
                                sub_r2_nulls(r) = stats(1);
                                
                                sub_res_with_null(1:numel(res), r) = res;
                                
                                for reg=1:N_regs_data
                                    ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                                    %sub_msr_with_nulls(fac,rec,reg,r) = sumsqr(res(ind_curr_reg))/numel(ind_subject);
                                    sub_msr_with_nulls(reg,r) = sumsqr(res(ind_curr_reg))/numel(ind_subject);
                                end
                            end
                            f_name = sprintf('%s/fac_%d/rec_%d/subject_%d.mat', null_res_savedir, fac,rec, i);
                            save(f_name, 'sub_r2_nulls', 'sub_res_with_null', 'sub_msr_with_nulls', 'sub_leverage_with_null')
                    
                            
                        end
                    end
                    
                    
                    % Save all sub stats
                    f_name = sprintf("output/modelCode=10.mat");
                    save(f_name, 'msr_with', 'msr_without', 'r2_with', 'r2_without', ...
                       'N_iters', 'res_with', 'res_without', 'leverage_with', 'leverage_without')%, ...
                       %'r2_nulls', 'res_with_null', 'msr_with_nulls', 'leverage_with_null');
                        
                    end
                    
                    if modelCode == 11
                        
                        % Residuals with and without imaging factors
                        
                        % Model performance with and without a specific receptor at a region
                        % Align timepoints and regions in 2D
                        % [fac, reg, timepoints]
                        
                        %inter_fac = squeeze(delta_SZ(1:N_facs, :))';
                        % Target imaging modality predicted only by itself
                        for tar_fac=1:N_facs
                            %y = squeeze(reshape(dS_dt_perm(tar_fac,:,:), 1, numel(ind_subject)*N_regs_data))';
                            y = squeeze(dS_dt_perm_tp(tar_fac,:))'; % tp1reg123 tp2reg123...
                            
                            net_spread = reshape(squeeze(subject_net_spreads(tar_fac, :, :))', numel(y), 1);
                            X = [ones(numel(y),1) delta_SZ_tp(tar_fac,:)' net_spread];
                            
                            X(:,2:end) = normalize(X(:,2:end));
                            y = normalize(y);
                            
                            leverage_without(i,tar_fac,1:numel(y)) = leverage(X(:,2:end));
                            [~, ~, res, ~, stats] = regress(y,X);
                            res_without(i,tar_fac,1:numel(res)) = res;
                            
                            r2_without(i,tar_fac) = stats(1);
                            for reg=1:N_regs_data
                                %ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                                ind_curr_reg = [reg];
                                for d_ind=1:(numel(ind_subject)-1) % Add indices of all time points
                                    ind_curr_reg = [ind_curr_reg reg+(N_regs_data*d_ind)];
                                end
                                msr_without(i,tar_fac,reg) = sumsqr(res(ind_curr_reg))/numel(ind_subject) ;
                            end
                        end
                        
                        
                        % Target imaging modality predicted by itself + 1
                        % other imaging modality at a time
                        for tar_fac=1:N_facs
                            sprintf('[%d of %d] Model Code 11: Tar fac %d', i, numel(unique_subject_ids), tar_fac)
                            mkdir(fullfile(null_res_savedir, sprintf('tar_fac_%d', tar_fac)));
                            for src_fac=1:N_facs
                                sprintf('[%d of %d] Model Code 11: Src fac %d | Tar fac %d', i, numel(unique_subject_ids), src_fac, tar_fac)
                                mkdir(fullfile(null_res_savedir, sprintf('tar_fac_%d', tar_fac), sprintf('src_fac_%d', src_fac)));
                                
                                sub_r2_nulls  = zeros(N_iters,1);
                                sub_msr_with_nulls = zeros(N_regs_data, N_iters);
                                sub_res_with_null = zeros(N_regs_data * 3, N_iters);
                                sub_leverage_with_null = zeros(N_regs_data * 3, N_iters);
                                
                                % Regression
                                %y = squeeze(dS_dt_perm(tar_fac, :))';
                                %inter_fac = squeeze(delta_SZ(1:N_facs, :))';
                                y = squeeze(dS_dt_perm_tp(tar_fac,:))';
                                
                                % Net spreading (factor, subject, region)
                                % Reg123_tp1... Reg123_tp2...
                                net_spread = reshape(squeeze(subject_net_spreads(tar_fac, :, :))', numel(y), 1);
                                
                                %X = [ones(numel(y),1), inter_fac(:, [tar_fac src_fac])  net_spread];
                                if tar_fac ~= src_fac
                                    X = [ones(numel(y),1), delta_SZ_tp([tar_fac src_fac], :)'  net_spread];
                                else
                                    X = [ones(numel(y),1), delta_SZ_tp(tar_fac,:)'  net_spread];
                                end
                                
                                y = normalize(y);
                                X(:,2:end) = normalize(X(:,2:end));
                                
                                leverage_with(i,tar_fac,src_fac,1:numel(y)) = leverage(X(:,2:end));
                                [b_temp, ~, res, ~, stats] = regress(y,X);
                                
                                res_with(i,tar_fac,src_fac,1:numel(res)) = res;
                                
                                r2_with(i,tar_fac,src_fac) = stats(1);
                                for reg=1:N_regs_data
                                    %ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                                    ind_curr_reg = [reg];
                                    for d_ind=1:(numel(ind_subject)-1) % Add indices of all time points
                                        ind_curr_reg = [ind_curr_reg reg+(N_regs_data*d_ind)];
                                    end
                                    msr_with(i,tar_fac,src_fac,reg) = sumsqr(res(ind_curr_reg))/numel(ind_subject) ;
                                end
                                
                                % Random permutations
                                for r=1:N_iters
                                    
                                    % Shuffle source but not target factor   
                                    if tar_fac == src_fac
                                        rand_inter_fac = delta_SZ_tp(tar_fac,:)';
                                    else
                                        rand_inter_fac = [delta_SZ_tp(tar_fac,:)' ...
                                            delta_SZ_tp(src_fac,randperm(size(delta_SZ_tp,2)))'];
                                    %rand_inter_fac = [inter_fac(:, tar_fac) ...
                                    %    inter_fac(randperm(size(inter_fac,1)), src_fac)];
                                    end
                                    
                                    X = [ones(numel(y),1), rand_inter_fac  net_spread];
                                    y = normalize(y);
                                    X(:,2:end) = normalize(X(:,2:end));
                                    [~, ~, res, ~, stats] = regress(y,X);
                                    
                                    sub_leverage_with_null(1:numel(y),r) = leverage(X(:,2:end));
                                    
                                    sub_r2_nulls(r) = stats(1);
                                    
                                    sub_res_with_null(1:numel(res), r) = res;
                                    
                                    for reg=1:N_regs_data
                                        ind_curr_reg = ((reg-1)*numel(ind_subject)) + 1: reg*numel(ind_subject);
                                        %sub_msr_with_nulls(fac,rec,reg,r) = sumsqr(res(ind_curr_reg))/numel(ind_subject);
                                        sub_msr_with_nulls(reg,r) = sumsqr(res(ind_curr_reg))/numel(ind_subject);
                                    end
                               end
                               f_name = sprintf('%s/tar_fac_%d/src_fac_%d/subject_%d.mat', null_res_savedir, tar_fac, src_fac, i);
                               % Save imaging for debugging
                               save(f_name, 'sub_r2_nulls', 'sub_res_with_null', 'sub_msr_with_nulls', 'sub_leverage_with_null');
                               
                            end
                        end
                        
                        % Save all sub stats
                        f_name = sprintf("output/modelCode=11.mat");
                        save(f_name, 'msr_with', 'msr_without', 'r2_with', 'r2_without', ...
                            'N_iters', 'res_with', 'res_without', 'leverage_with', 'leverage_without');

                    end
                    
                end
                
                r2s_all_diags = cat(1, r2s_all_diags, r2s_individual);%squeeze(max(r2s_individual,[])));
                rss_all_diags = cat(1, rss_all_diags, RSS_individual);%squeeze(max(RSS_individual,[])));
                ps_all_diags = cat(1, ps_all_diags, ps_individual);
                r2s_sig_all_diags = cat(1, r2s_sig_all_diags, r2s_sig_individual);
                rss_sig_all_diags = cat(1, rss_sig_all_diags, rss_sig_individual);
                bs_all_diags = cat(1, bs_all_diags, bs_individual);
                
                f_name = sprintf("output/individual_regression_DIAG=%d_modelCode=%d", DIAG_ID,modelCode);
                if modelCode == 7
                 residual_model_comp_all_diags = cat(1, residual_model_comp_all_diags, residual_model_comp);
                save(f_name, 'DIAG_ID', 'bs_individual', 'method', ...
                    'r2s_individual', 'RSS_individual', 'ps_individual',...
                    'res_individual', 'r2s_sig_all_diags', 'rss_sig_all_diags',...
                    'net_spreads', 'd_diag_samples_subj_ids', 'residual_model_comp_all_diags'); % 'cum_effect_individual',
                elseif modelCode < 7
                    save(f_name, 'DIAG_ID', 'bs_individual', 'method', ...
                    'r2s_individual', 'RSS_individual', 'ps_individual',...
                    'res_individual', 'r2s_sig_all_diags', 'rss_sig_all_diags',...
                    'net_spreads', 'd_diag_samples_subj_ids'); % 'cum_effect_individual'

                    
                end
                end
            end
        end
        r2s_random_maps = cat(3,r2s_random_maps, r2s_all_diags);
        rss_random_maps = cat(3,rss_random_maps, rss_all_diags);
        bs_random_maps = cat(4,  bs_random_maps, bs_all_diags);
        group_bs_random_maps = cat(4, group_bs_random_maps, group_bs_all_diags);

    %    if ~isRandomZ 
        p = size(bs_individual,3);
        %        if (includeInteract == true)

        %save('output/r2s_interact.mat', 'r2s_all_diags', 'rss_all_diags', 'p', 'ps_all_diags', 'bs_all_diags', 'group_bs_all_diags');
        %        else

        save(sprintf('output/r2s_modelCode=%d.mat',modelCode), 'r2s_all_diags',...
            'rss_all_diags', 'p', 'ps_all_diags', 'bs_all_diags', 'group_bs_all_diags', 'res_all_diags');
        %        end

    %    end

        % Save every 100 runs
        if isRandomZ %%&& ~(mod(iter, 100))
            save(sprintf('output/Random_Maps/r2s_random_maps_%d.mat', iter), ...
                'r2s_random_maps', 'rss_random_maps', 'bs_random_maps', 'group_bs_random_maps');
            r2s_random_maps = [];
            rss_random_maps = [];
            bs_random_maps = [];
            group_bs_random_maps = [];
        end

    end
end
sprintf("Done\n");



save('output/mcm_num_subj_samples.mat', 'mcm_num_subj_samples')



%% R2 of random maps

close all; 

mkdir('output/figures')
mkdir('output/figures/paper')

r2s_random = [];
%rss_random = [];
bs_random = [];
for i=1:N_RAND_ITERS
    fname = sprintf("output/Random_Maps/r2s_random_maps_%d.mat", i);
    if isfile(fname)
        load(fname, 'r2s_random_maps', 'bs_random_maps')
        r2s_random = cat(3, r2s_random, r2s_random_maps);
        %rss_random = cat(3, rss_random, rss_random_maps);
    end
end
% % % % % 
%load('output/Random_Maps/r2s_random_maps.mat', 'r2s_random_maps', 'bs_random_maps');
close all;

r2s_facs = zeros(N_facs,1);
for fac=1:N_facs
    r2s_facs = mean(mean(r2s_random,3),1);
end


load('output/r2s_modelCode=3.mat', 'r2s_all_diags', 'bs_all_diags');
r2s_all_diags_interact = r2s_all_diags;

r2s_mean_fac_diags = zeros(N_facs, numel(diag_names));
r2s_std_fac_diags = zeros(N_facs, numel(diag_names));
for diag=1:numel(diag_names)
    for fac=1:N_facs
        r2s_mean_fac_diags(fac,diag) = mean(r2s_all_diags_interact(subjects_diags==diag, fac));
        r2s_std_fac_diags(fac,diag) = std(r2s_all_diags_interact(subjects_diags==diag, fac));
        
    end
end

% Check how subject model fits distribution of R2s of null receptor models
r2s_zscores = zeros(size(r2s_random,1),N_facs);
r2s_random_means = mean(r2s_random, 3);
for subject=1:size(r2s_random,1)
    for fac=1:N_facs
        curr_random_r2s = r2s_random(subject,fac,:);
        r2s_zscores(subject, fac) = (r2s_all_diags_interact(subject, fac) - mean(curr_random_r2s))/std(curr_random_r2s);
    end
end

r2s_pvals = normpdf(r2s_zscores);

r2_rand_proportion = zeros(N_facs, numel(diag_names));

% Boxplots (interleave)

r2s_box = zeros(size(r2s_random,1), N_facs * 2);
labels = {};
labels_legend = {};
labels_interact = {};
labels_noninteract = {};

means = mean(r2s_pvals);
stds = std(r2s_pvals);
pval_counts = zeros(N_facs, 1);
diags_under = {};
diags_over = {};
for fac=1:N_facs
    if fac == 3
       labels{(2*fac)-1} = strcat('Activity', ' (Data)');
       labels{2*fac} = strcat('Activity', ' (Null)');
    else
       labels{(2*fac)-1} = strcat(fac_list{fac}, ' (Data)');
       labels{2*fac} = strcat(fac_list{fac}, ' (Null)');
    end
    
    r2s_box(:, (2*fac)-1) = r2s_all_diags_interact(:,fac);
    r2s_box(:, (2*fac)) = r2s_random_means(:,fac);
    
    ind_fac_under = find(r2s_pvals(:,fac)<0.05);
    pvals_counts(fac) = numel(ind_fac_under);
    diags_under{fac} = subjects_diags(ind_fac_under);
    diags_over{fac} = subjects_diags(setdiff(1:size(r2s_pvals,1), ind_fac_under));
    
    for diag=1:numel(diag_names)
        r2_rand_proportion(fac,diag) = nnz(diags_under{fac}==diag) / (nnz(diags_under{fac}==diag)+nnz(diags_over{fac}==diag));
    end
    
    
    labels_legend{fac} = strcat(string(fac_list(fac)), sprintf(" (%.1f%%)", 100*pvals_counts(fac)/size(r2s_pvals,1) ));
end



% Diagnosis for p<0.05 and p>0.05


for fac=1:N_facs
    figure;
    diff = numel(diags_under{fac}) - numel(diags_over{fac});
    if diff < 0
        % More over than under
        
        disp_hist = [diags_under{fac} NaN(1,-diff); diags_over{fac}];
    else
        if diff > 0
            disp_hist = [diags_under{fac}; diags_over{fac} NaN(1,diff)];
        else
            disp_hist = [diags_under{fac}; diags_over{fac}];
        end
    end 
    hist(disp_hist')
    ylabel("Frequency");
    xlabel("Diagnosis");
    xticks(1:numel(diag_names));
    
    
    xticklabels(string(diag_names));
    f_name=sprintf("output/figures/r2_diags_%d.svg", fac);
    %temp = [{["p<0.05"]}, {["p>0.05"]}];
    legend('p<0.05', 'p>0.05');%, 'Location','northeast');
    title(sprintf("%s", string(fac_list(fac))));
    saveas(gcf, f_name);
end
% % % 
close all;

figure('Renderer', 'painters', 'Position', [10 10 800 600])
bar(100*r2_rand_proportion(:,3))
ylabel("Proportion p<0.05 (%%)");
xlabel("Diagnosis");
%xticks(1:numel(diag_names));
%xticklabels(string(diag_names));
xticklabels(string(fac_list))

xlim([0.5 numel(fac_list)+0.5])
f_name=sprintf("output/figures/paper/r2_random_diag_dist.svg", fac);
%temp = [{["p<0.05"]}, {["p>0.05"]}];
%legend(string(fac_list), 'Color', 'none');%, 'Location','northeast') 
%legend boxoff
%set(lgnd, 'color', 'none');
title('Proportion of significant subjects by diagnosis');
saveas(gcf, f_name);


% % % 
figure
hist(r2s_pvals)
ylabel("Frequency");
xlabel("P Value");
f_name=sprintf("output/figures/paper/r2s_vs_null_pvals.svg");
legend(labels_legend, 'Location','northeast', 'Color', 'none');
title("p-values of individual model R^2");
saveas(gcf, f_name);


%% Fisher's method

combined_pvals = zeros(N_facs,1);
combined_chi_sq = zeros(N_facs,1);

for fac=1:N_facs
    
    % Effects: Fisher's method for combining p values
    curr_pvals = squeeze(r2s_pvals(:,fac));
    curr_chi = -2.*sum(log(curr_pvals));
    combined_chi_sq(fac) = curr_chi;
    combined_pvals(fac) = 1 - chi2cdf(sum(curr_chi), 2*length(curr_pvals));
end


%% Supplementary table: performance gain over random maps

    
 r2s_improvement_rand = zeros(N_facs,size(r2s_all_diags_interact,1));
 for fac=1:N_facs
     diff = r2s_all_diags_interact(:,fac) - r2s_random_means(:,fac);
     sprintf('%s: %.3f +/- %.3f (P<%.3f)', fac_list(fac), mean(diff)*100, std(diff)*100, combined_pvals(fac) )
 end


% Permutation test - coefficients

close all;

load('output/r2s_modelCode=3.mat', 'r2s_all_diags', 'bs_all_diags');

r2s_all_diags_interact = r2s_all_diags;

%load('output/r2s_random_maps.mat', 'r2s_random_maps', 'bs_random_maps');


% p value for each coefficient
N_params =  size(bs_all_diags,3);
N_subjects = size(bs_all_diags,1) ;
b_pval = zeros(N_facs, N_params, N_subjects);
pval_counts = zeros(N_facs,N_params);
for fac=1:N_facs
    bs_random = [];
    for i=1:N_RAND_ITERS
        fname = sprintf("output/Random_Maps/r2s_random_maps_%d.mat", i);
        if isfile(fname)
            load(fname, 'bs_random_maps');
            if numel(size(bs_random_maps)) == 4
                bs_random = cat(3,bs_random, squeeze(bs_random_maps(:,fac,:,:)));
            else
                bs_random = cat(3,bs_random, squeeze(bs_random_maps(:,fac,:)));
            end

        end
    end
    for p=1:N_params
       
        for subject=1:N_subjects
            z = (bs_all_diags(subject,fac,p) - mean(bs_random(subject,p,:))) / std(bs_random(subject,p,:));
            b_pval(fac,p,subject) = normpdf(z);
        end
    end
    
    p_thresh = 0.05/N_params;
    figure;
    pvals = squeeze(b_pval(fac,:,:));
    pvals_plot = zeros(N_params,N_subjects);
    pvals_plot(find(pvals > p_thresh)) = 1;
    imagesc(pvals);%_plot);
    for param=1:N_params
        pvals_counts(fac,param) = nnz(pvals_plot(param,:));
    end
    
    
    xlabel('Subject')
    ylabel('Parameter')
    title(sprintf("%s", string(fac_list(fac))))
    saveas(gcf, sprintf("output/figures/param_pval_%s.svg", string(fac_list(fac))));
end

% Coefficient p value confidence intervals
% (x+/- t*s/sqrt(n))
b_ci = zeros(N_facs, N_params, 2); % [lower, upper] bounds
b_mean = zeros(N_facs, N_params);
N_subjects = size(bs_all_diags,1);
t_crit = tinv(0.99, N_subjects-1); % Control for multiple comparisons?
param_ci_keep = zeros(N_facs,N_params);
for fac=1:N_facs
    for param=1:N_params
        b_mean(fac,param) = mean(bs_all_diags(:, fac,param));
        b_std = std(bs_all_diags(:,fac,param));
        b_ci(fac,param,1) = b_mean(fac,param) - (t_crit*b_std)/sqrt(N_subjects);
        b_ci(fac, param,2) = b_mean(fac,param) + (t_crit*b_std)/sqrt(N_subjects);
        if sign(b_ci(fac,param,1)) == sign(b_ci(fac,param,2))
            param_ci_keep(fac,param) = 1;
        end
    end
end
for fac=1:N_facs
    figure;
    plot(1:N_params, b_ci(fac,:,1), '--b')
    hold on
    plot(1:N_params, b_mean(fac,:), '-r')
    plot(1:N_params, b_ci(fac,:,2), '--b')
    hold off
    xlabel('Parameter')
    ylabel('Value')
    title(sprintf("%s", string(fac_list(fac))))
    saveas(gcf, sprintf("output/figures/param_ci_%s.svg", string(fac_list(fac))));
end

%%
% CI param names
fid_ci_params = fopen("output/params_ci_names.txt", 'w');
% Indices and mean values of parameters in 99% CI 
ci_param_inds_all = [];
b_mean_ci = [];
for fac=1:N_facs
    fprintf(fid_ci_params,"%s: ", fac_list(fac))
    curr_p = find(param_ci_keep(fac,:));
    for p=1:nnz(param_ci_keep(fac,:))
        b_mean_ci = [b_mean_ci b_mean(fac,curr_p(p))];
        fprintf(fid_ci_params, "%s, ", param_names(curr_p(p))) 
    end
    fprintf(fid_ci_params,"\n")
    
    param_offset = (fac-1)*numel(param_names);
    ci_param_inds_all = [ci_param_inds_all (curr_p + param_offset)];
end

%% 


%% If predictors are normalized before regression, parameters can be directly compared
%PLSClusters(sig_params_ind, sig_vars, param_ind, suffix, N_facs, N_recs,rec_type_codes,param_names,param_rec_code,param_fac_code, fac_list,save_dir)

% Repeated from ImagingCognitiveSVD.m
rec_type_codes = rec_type_codes(1:N_recs);
rec_type_codes(N_recs+1) = numel(unique(rec_type_codes))+1;

for i=1:numel(param_fac_code)
   if param_fac_code(i) == 0
       param_fac_code(i) = N_facs+1;
   end
end

for i=1:numel(param_rec_code)
   if param_rec_code(i) == 0
       param_rec_code(i) = N_recs+1;
   end
end

load(sprintf("output/individual_regression_DIAG=3_modelCode=3.mat"));
N_params = size(bs_individual,3);
param_ind_all = [];
for fac=1:N_facs
    param_ind_all = [param_ind_all; [repelem(fac,N_params)' (1:N_params)' rec_type_codes(param_rec_code(1:N_params))']];
end

PLSClusters(ci_param_inds_all, abs(b_mean_ci), param_ind_all, "CI_params", N_facs,...
    N_recs,rec_type_codes,param_names,param_rec_code,param_fac_code, fac_list, "ci_params")

save('output/param_ci_keep.mat', 'param_ci_keep', 'b_ci', 'b_mean', 'param_rec_code', 'param_fac_code', 'bs_all_diags');



% % R2


rmpath(genpath('/export02/data/Work/MATLAB/spm12'))
load('output/r2s_modelCode=3.mat', 'rss_all_diags', 'p', 'r2s_all_diags', 'ps_all_diags');
p_interact = p;
rss_interact = rss_all_diags;
r2s_interact = r2s_all_diags;
param_p_interact = ps_all_diags;
load('output/r2s_modelCode=2.mat', 'rss_all_diags', 'p', 'r2s_all_diags', 'ps_all_diags');
p_noninteract = p;
rss_noninteract = rss_all_diags;
r2s_noninteract = r2s_all_diags;
param_p_noninteract = ps_all_diags;
load('output/r2s_modelCode=1.mat', 'rss_all_diags', 'p', 'r2s_all_diags', 'ps_all_diags');
p_imaging = p;
rss_imaging = rss_all_diags;
r2s_imaging = r2s_all_diags;
param_p_imaging = ps_all_diags;
load('output/r2s_modelCode=5.mat', 'rss_all_diags', 'p', 'r2s_all_diags', 'ps_all_diags');
p_julich = p;
rss_julich = rss_all_diags;
r2s_julich = r2s_all_diags;
param_p_julich = ps_all_diags;
load('output/r2s_modelCode=6.mat', 'rss_all_diags', 'p', 'r2s_all_diags', 'ps_all_diags');
p_julich = p;
rss_imputed = rss_all_diags;
r2s_imputed = r2s_all_diags;
param_p_imputed = ps_all_diags;

% %  %
% % R2s by factor
mean_interact = mean(r2s_interact);
mean_noninteract = mean(r2s_noninteract);
mean_imaging = mean(r2s_imaging);
mean_julich = mean(r2s_julich);
mean_imputed = mean(r2s_imputed);
std_interact = std(r2s_interact);
std_noninteract = std(r2s_noninteract);
std_imaging = std(r2s_imaging);
std_julich = std(r2s_julich);
std_imputed = std(r2s_imputed);

labels_interact = {};
labels_noninteract = {};
labels_imaging = {};
labels_julich = {};
labels_imputed = {};

for fac=1:N_facs
   sprintf("%s %.3f %.3f | %.3f %.3f", string(fac_list(fac)), mean_noninteract(fac), std_noninteract(fac), mean_interact(fac), std_interact(fac) ) 
    
   labels_interact{fac} = strcat(string(fac_list(fac)), sprintf(" (%.2f ± %.2f)", mean_interact(fac), std_interact(fac) ));
   labels_noninteract{fac} = strcat(string(fac_list(fac)), sprintf(" (%.2f ± %.2f)", mean_noninteract(fac), std_noninteract(fac) ));
   labels_imaging{fac} = strcat(string(fac_list(fac)), sprintf(" (%.2f ± %.2f)", mean_imaging(fac), std_imaging(fac) ));
   labels_julich{fac} = strcat(string(fac_list(fac)), sprintf(" (%.2f ± %.2f)", mean_julich(fac), std_julich(fac) ));
   labels_imputed{fac} = strcat(string(fac_list(fac)), sprintf(" (%.2f ± %.2f)", mean_imputed(fac), std_imputed(fac) ));
end

% % %

figure('Renderer', 'painters', 'Position', [10 10 600 300])
hist(r2s_interact, 20)
xlabel("R^2")
ylabel("Frequency")
%xlim([0 1]);
%title(sprintf("%s, N=%d", string(diag_names(DIAG_ID)), size(r2s_all_diags,1)));
title(sprintf("Full re-MCM (Neuroimaging x Receptors)"));
f_name=sprintf("output/figures/paper/hist_interact.svg");
legend(labels_interact, 'Location','northwest');
saveas(gcf, f_name);

figure('Renderer', 'painters', 'Position', [10 10 600 300])
hist(r2s_noninteract, 20)
xlabel("R^2")
%xlim([0,1]);
ylabel("Frequency")
%title(sprintf("%s, N=%d", string(diag_names(DIAG_ID)), size(r2s_all_diags,1)));
title(sprintf("Restricted re-MCM (Neuroimaging + Receptors)"));
f_name=sprintf("output/figures/paper/hist_noninteract.svg");
legend(labels_noninteract, 'Location','northeast');
saveas(gcf, f_name);

figure('Renderer', 'painters', 'Position', [10 10 600 300])
hist(r2s_imaging, 20)
xlabel("R^2")
%xlim([0,1]);
ylabel("Frequency")
%title(sprintf("%s, N=%d", string(diag_names(DIAG_ID)), size(r2s_all_diags,1)));
title(sprintf("Restricted re-MCM (Neuroimaging)"));
f_name=sprintf("output/figures/paper/hist_imaging.svg");
legend(labels_imaging, 'Location','northeast');
saveas(gcf, f_name);

figure('Renderer', 'painters', 'Position', [10 10 600 300])
hist(r2s_julich, 20)
xlabel("R^2")
%xlim([0,1]);
ylabel("Frequency")
%title(sprintf("%s, N=%d", string(diag_names(DIAG_ID)), size(r2s_all_diags,1)));
title(sprintf("Full re-MCM (Julich regions)"));
f_name=sprintf("output/figures/paper/hist_julich_regions.svg");
legend(labels_julich, 'Location','northeast');
saveas(gcf, f_name);

figure('Renderer', 'painters', 'Position', [10 10 600 300])
hist(r2s_imputed, 20)
xlabel("R^2")
%xlim([0,1]);
ylabel("Frequency")
%title(sprintf("%s, N=%d", string(diag_names(DIAG_ID)), size(r2s_all_diags,1)));
title(sprintf("Full re-MCM (imputed regions)"));
f_name=sprintf("output/figures/paper/hist_imputed_regions.svg");
legend(labels_imputed, 'Location','northeast');
saveas(gcf, f_name);

%%
% % %
% % f test

rmpath(genpath('/export02/data/Work/MATLAB/spm12'))


load('output/mcm_num_subj_samples.mat')
mcm_subject_ids = [];
for DIAG_ID=3%1:numel(diag_names)
    load(sprintf("output/individual_regression_DIAG=%d_modelCode=3.mat", DIAG_ID));
    if numel(d_diag_samples_subj_ids) > 0
    mcm_subject_ids = cat(1,mcm_subject_ids, unique(d_diag_samples_subj_ids)');
    end
end

% REPLACE WHEN CRITERIA CHANGES
mcm_tab_facs = ones(numel(mcm_subject_ids), numel(fac_list));
%subjects_tab_facs(mcm_subject_ids,:);


%p2 = 113;
%p1 = 23;

close all

for fac=1:N_facs
    [H,Ptemp] = ttest2(squeeze(r2s_interact(:,fac)), squeeze(r2s_imaging(:,fac)));
    %[H,Ptemp] = ttest(squeeze(r2s_interact(fac,:) - r2s_imaging(fac,:)));
    r2_diff = r2s_interact(:,fac)-r2s_imaging(:,fac);
    %sprintf("%s, Mean R2 with =%.3f | without = %.3f", fac_list(fac), mean(r2s_interact(:,fac)), mean(r2s_imaging(:,fac)))
    %sprintf("%s, Mean gain in R2 = %.3f +/- .%3f | P=%.4f", fac_list(fac),mean(r2_diff), std(r2_diff), Ptemp)
    Ptemp
end

path_figure=sprintf("output/figures/fstat");

[~, fs_crit] = ModelFTest(rss_interact,rss_imaging, p_interact, p_imaging, 0.05, ...
    mcm_tab_facs, mcm_num_subj_samples * numel(ind_regs_data), param_p_interact, param_p_imaging, fac_list, path_figure, subjects_diags)

save('output/figures/paper/plot_r2s_vs_null.mat', 'r2s_zscores', 'r2s_pvals', 'r2s_all_diags_interact', 'r2s_random', 'rss_interact', 'rss_imaging',...
    'r2s_interact', 'r2s_imaging', 'mcm_num_subj_samples', 'fs_crit')

%%
clear ind* full_tab* ledd* data*
save('workspaces/after_model_tests.mat')

load('workspaces/after_model_tests.mat')
