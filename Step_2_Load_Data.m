% =========================================================================
% Load NeuroPM-box data
% =========================================================================
%
% DESCRIPTION:
% Prepare neuroimaging, molecular template and connectivity data for
% molecular-MCM analysis.
%
% This script takes as input preprocessed, harmonized and ROI-averaged multi-modal neuroimaging
% data organized using NeuroPM-box (https://www.neuropm-lab.com/neuropm-box-download.html).
%
% =========================================================================
% INPUTS:
%   - mcm_data.mat: Path to MCM data structure from NeuroPM-box
%   - connectivity.mat: Path to connectivity matrix .mat file (N_ROIs x N_ROIs)
%   - molecular_template.mat: Path to molecular template .mat file (N_ROIs x N_FEATURES)
%   - atlas.mat: Path to parcellation file (output/atlas_creation_subcortical)
%   - clincial.mat: Clinical data
%
% INPUT PARAMETERS:
%   - fac_list: List of neuroimaging modalities in the same order as the
%       NeuroPM-box preprocessing folder.
%       e.g., ["GM", "fALLF", "SPECT", "FA", "MD", "t1t2"];
%   - essential_facs: Required neuroimaging modalities to include a subject.
%       e.g., [1 3 6] to include GM, SPECT, t1t2
%   - imputed_facs: modalities that are okay to count if they have been imputed
%       e.g., [4 5 6] for FA, MD, t1t2
%
%
% =========================================================================
% OUTPUTS:
%   - data_for_MCM.mat:
%
% =========================================================================
% DEPENDENCIES:
%   - Functional connectivity toolbox (https://web.conn-toolbox.org/home)
%
% =========================================================================
% AUTHOR: Ahmed Faraz Khan
% CONTACT: ahmed.faraz.khan@mai.mcgill.ca
% VERSION: 1.0.0
% =========================================================================
                                       
% =================Rename==========================================
ppmi_data,
% =================Rename==========================================
                        
% Paths
addpath('/export02/data/Work/conn18b')
addpath(genpath('/export02/data/Work/MATLAB/'))

% Inputs
fac_list = ["GM", "fALLF", "SPECT", "FA", "MD", "t1t2"];
essential_facs = [1 3 6]; % GM, SPECT, t1/t2
imputed_facs = [4 5 6]; % FA, MD, t1/t2 (not relevant anymore)
ppmiDataPath = '/export02/data/Work/re-MCM_PD/data/MCM_Structures/Input_data_05-Nov-2021_11:26:01_MCM_imputed.mat';
ppmiImputedDataPath = '/export02/data/Work/re-MCM_PD/data/MCM_Structures//Input_data_05-Nov-2021_11:26:01_MCM_imputed.mat';

% Update imputed data and MCM data file names
load('/export02/data/Work/re-MCM_PD/data/imputed/GAIN_imputed_data_11-Jan-2022.mat'); % imputed_data
load('/export02/data/Work/re-MCM_PD/data/imputed/data_sub_inds_06-Jan-2022.mat'); % gain_sub_inds
load('/export02/data/Work/re-MCM_PD/data/imputed/data2impute_06-Jan-2022.mat'); %  data_2impute
load('/export02/data/Work/re-MCM_PD/data/imputed/mods_per_sub_05-Aug-2022.mat'); % is_mods_per_sub: modalities present per subject

load('output/atlas_creation_subcortical'); % V_, Z_, etc.
                                       
% Set to true to use GAN-imputed data (false for default MCM imputed data)
useGAIN = 0; 
useImputed = 1;

%% Combine multimodal data

% Levadopa-equivalent dose
%load('data/LEDDcomputed.mat'); % rowWithMotorPDmed...

%N_regs = size(Z0_Julich_scat, 1);

restoredefaultpath
rehash toolboxcache

ppmi_MCM = load(ppmiDataPath);
ppmi_MCM_imputed = load(ppmiImputedDataPath);
%ppmi_MCM_not_imputed = load(ppmiDataPath_not_imputed); 
ppmi_summ = load('data/summaryOfImagingData_February14_2021.mat');
ppmi_summ = ppmi_summ.summaryOfImagingData;

N_facs = ppmi_MCM.N_modalities;
N_regs = ppmi_MCM.N_im_regions;

ppmi_data = extractfield(ppmi_MCM, 'Data');
ppmi_data = ppmi_data{1};
ppmi_data_imputed = extractfield(ppmi_MCM_imputed, 'Data');
ppmi_data_imputed = ppmi_data_imputed{1};

MIN_SAMPLES = 3;
keep_ids = [];
keep_ppmi_ind = [];
modality_counts = [];
all_modality_count = 0;
subjects_tab_facs = zeros(ppmi_MCM.N_subjects, numel(fac_list));

% Check IDs and modalities
IDs_preprocessed = zeros(ppmi_MCM.N_subjects,1);
num_timepoints = zeros(ppmi_MCM.N_subjects,N_facs);

% Keep subjects with t1, SPECT, t1/t2
for subject=1:ppmi_MCM.N_subjects    
    
	n_samples = size(ppmi_data(subject).raw_data,2);
    % Times by itself doesn't guarantee there are that many data points...
    times = ppmi_data(subject).times;
    
    %%%%% Check IDs for some imaging modality count/time point criteria
    IDs_preprocessed(subject) = str2num(ppmi_data(subject).id);
    
    for fac=1:N_facs
        num_timepoints(subject,fac) = nnz(times(fac,:));
    end
    %%%%%%%
    
    % At keast 3 time points (make sure numbers of imaging and times match)
    modality_counts_subject = [];
    for row=1:numel(fac_list)
        modality_counts_subject = [modality_counts_subject nnz(times(row,:))];
        subjects_tab_facs(subject, row) = any(times(row,:));
    end
    if nnz(modality_counts_subject) == 3
        all_modality_count = all_modality_count +1;
    end
    modality_counts = [modality_counts; modality_counts_subject];

    % Check .S because outliers may have been excluded

    % All 3 essential modalities (GM, SPECT, t1/t2), at least 3 time points for 1 modality
    if (nnz(modality_counts_subject(essential_facs)) == 3 ) && (nnz(modality_counts_subject) > 3)  ... 
            && (n_samples >= MIN_SAMPLES) ...
            && (size(times,2) >= MIN_SAMPLES) ...
            && (size(ppmi_data(subject).S,2) >= (MIN_SAMPLES-1))
        times
        keep_ids = [keep_ids; round(str2double(ppmi_data(subject).id))];
        keep_ppmi_ind = [keep_ppmi_ind; subject];
    end
end

for fac=1:N_facs
    mean(num_timepoints(:,fac))
end

%% Healthy controls' modalities

summ_ids = cell2mat(ppmi_summ(2:end, 1));
summ_diag = ppmi_summ(2:end, 2);

temp_ind_control = find(strcmp(summ_diag, string('Control')));
%% Neuroimaging modalities per subject

for fac=1:N_facs
    sprintf("%s: %d subjects (%.0f%%)", fac_list(fac), nnz(subjects_tab_facs(keep_ppmi_ind,fac)),...
        100*nnz(subjects_tab_facs(keep_ppmi_ind,fac))/numel(subjects_tab_facs(keep_ppmi_ind,fac)))
end

%%
num_subj=zeros(6,1);
for mod=1:6
    num_subj(mod) = nnz(subjects_tab_facs(:,mod));
    
end

%
subjects_tab_facs = subjects_tab_facs(keep_ppmi_ind,:);
%% Get clinical data
gender_letters = 'MF';

keep_diags = [];
subjects_ages = [];
keep_gender = [];
keep_edus = [];
keep_hand = [];
keep_race = [];
keep_updrs_slope = [];
keep_moca_slope = [];

summ_ids = cell2mat(ppmi_summ(2:end, 1));
summ_diag = ppmi_summ(2:end, 2);
summ_ages = ppmi_summ(2:end, 3);
summ_gender = ppmi_summ(2:end, 4);
summ_race = ppmi_summ(2:end, 5);
summ_edus = ppmi_summ(2:end, 6);
summ_hand = ppmi_summ(2:end, 7);
summ_updrs_slope = ppmi_summ(2:end, 8);
summ_ledd_slope = ppmi_summ(2:end, 12);
summ_moca_slope = ppmi_summ(2:end, 13);

for i=1:numel(keep_ids)
    summ_ind = find(summ_ids == keep_ids(i));
    keep_diags = [keep_diags; summ_diag(summ_ind)];
    subjects_ages = [subjects_ages; cell2mat(summ_ages(summ_ind))];
    keep_gender = [keep_gender; find(gender_letters == cell2mat(summ_gender(summ_ind)))];
    keep_edus = [keep_edus; cell2mat(summ_edus(summ_ind))];
    keep_hand = [keep_hand; cell2mat(summ_hand(summ_ind))];
    keep_race = [keep_race; cell2mat(summ_race(summ_ind))];
    keep_updrs_slope = [keep_updrs_slope; cell2mat(summ_updrs_slope(summ_ind))];
    keep_moca_slope = [keep_moca_slope; cell2mat(summ_moca_slope(summ_ind))];
end

subjects_cog_slopes = [keep_updrs_slope, keep_ledd_slope, keep_moca_slope];
cog_names = {"UPDRS", "LEDD", "MOCA"};

diag_names = string(unique(keep_diags));
num_diags = zeros(numel(diag_names),1);
for i=1:numel(diag_names)
    num_diags(i)=numel(find(strcmp(keep_diags, string(diag_names(i)))));
    sprintf("%s: %d", string(diag_names(i)), num_diags(i))
end

figure('Renderer', 'painters', 'Position', [10 10 2500 900])
bar(num_diags)
xticklabels(diag_names)
ylabel("Number of subjects")
title(sprintf("Qualifying subjects by diagnosis (N=%d)", sum(num_diags)))
mkdir('output/figures');
f_name=sprintf("output/figures/diagnoses.png");
saveas(gcf, f_name);

subjects_diags = [];

% Calculate average LEDD (weighted by period of dosage)
clear subjects_ledd 
subjects_ledd_ave = zeros(numel(keep_ids),1);
subjects_ledd_cum = zeros(numel(keep_ids),1);
subjects_ledd_slope = zeros(numel(keep_ids),1);
subjects_ledd_int = zeros(numel(keep_ids),1);

ledd_patno = cell2mat(rowWithMotorPDmed_computedValues(2:end,2));
ledd = cell2mat(rowWithMotorPDmed_computedValues(2:end, 13));
ledd_drugtype = rowWithMotorPDmed_computedValues(2:end, 17);
ledd_startdates = cell2mat(rowWithMotorPDmed_computedValues(2:end, 8));
ledd_stopdates = cell2mat(rowWithMotorPDmed_computedValues(2:end, 9));

for i=1:numel(keep_ids)
    subjects_diags(i) = find(strcmp(diag_names, string(keep_diags(i))));
    
    % Some LEDD are NaN, and some stop dates are Inf
    ind_ledd = find(ledd_patno == keep_ids(i));
    sub_periods = ledd_stopdates(ind_ledd) - ledd_startdates(ind_ledd);
    sub_ledd = ledd(ind_ledd);
    ind_nonnan = ~isnan(sub_ledd) & isfinite(sub_periods);
    
    % Cumulative dose
    subjects_ledd_cum(i) = sum(sub_ledd(ind_nonnan) .* sub_periods(ind_nonnan));
    % Average dose
    subjects_ledd_ave(i) = subjects_ledd_cum(i) / sum(sub_periods(ind_nonnan));
    % Slope of dose
    mid_dates = (ledd_startdates(ind_ledd) + ledd_stopdates(ind_ledd)) / 2;
    mid_dates = mid_dates(ind_nonnan);
    if nnz(ind_nonnan)
        b = regress(sub_ledd(ind_nonnan), [ones(nnz(ind_nonnan),1) mid_dates]);
        subjects_ledd_int(i) = b(1);
        subjects_ledd_slope(i) = b(2);
    end
    subjects_ledd{i} = ledd(ind_ledd);
end

full_tab_facs =  extractfield(ppmi_MCM, 'table_factors');
full_tab_facs = full_tab_facs{1};
%full_tab_facs = reshape(full_tab_facs, numel(ppmi_MCM.N_subjects), N_facs);
keep_tab_facs = full_tab_facs(keep_ppmi_ind, 1:N_facs);


for fac=1:N_facs
    sprintf("Imputed %s: %d subjects (%.3f%%)", fac_list(fac), nnz(keep_tab_facs(:,fac)), ...
        100*nnz(keep_tab_facs(:,fac))/numel(keep_tab_facs(:,fac)))
end

%% Check healthy subjects' data

ind_healthy = [];
for i=1:numel(keep_diags)
    if strcmp(string(keep_diags(i)), 'Control')
        ind_healthy = [ind_healthy; i];
    end
end

subjects_tab_facs(ind_healthy,:)

%% Create dummy variables for covariate regression

subjects_male = keep_gender;
subjects_female = keep_gender - 1;
subjects_male(subjects_male > 1) = 0;

subjects_right_hand = keep_hand;
subjects_left_hand = ones(numel(keep_hand),1);
subjects_left_hand(subjects_right_hand > 0) = 0;

subjects_white = keep_race;
subjects_nonwhite = ones(numel(keep_race),1);
subjects_nonwhite(subjects_white > 0) = 0;

%% Extract longitudinal data
long_ages = [];
long_diags = [];
long_S = [];
long_background = [];  

subjects_background = [];

num_subj_samples = zeros(numel(keep_ids), 1);
num_imputed_subj_samples = zeros(numel(keep_ids), 1);
samples_subj_ids = [];
raw_mean_std = ppmi_MCM.raw_mean_std;
inds_exclude_rep = []; % Subjects to exclude because of repeated modalities

temp_repeated_age_counts = 0;

regions_mod_imputed = [];
for fac=1:N_facs
    if any(imputed_facs == fac)
        regions_mod_imputed = [regions_mod_imputed (fac-1)*N_regs+1:(fac-1)*N_regs+N_regs];
    end
end

untrained_imputation_counter = 0;
rawdat = zeros(size(data_2impute,1), size(data_2impute,2));
for j=1:numel(keep_ppmi_ind)
    subject = keep_ppmi_ind(j);
    
    % Imputation was optimized for FA,MD,t1/t2
    % Either re-impute the few remaining missing samples or  discard them 
    sub_imputed_data = imputed_data(gain_sub_inds == subject,:);
    sub_2impute_data = data_2impute(gain_sub_inds == subject,:);
    
    rawdat(gain_sub_inds==subject,:) = 1;i
    
    % Check if missing data pattern was unusual
    sprintf('Subject %d - %d timepoints, %d nan', j, size(sub_imputed_data,1),...
        nnz(isnan(sub_2impute_data(:))))

    for timepoint=1:size(sub_2impute_data,1)
    for fac=1:N_facs
        regions_mod = (fac-1)*N_regs+1:(fac-1)*N_regs+N_regs;
        
        % Missing non-imputed modalities
        if any(isnan(sub_2impute_data(timepoint,regions_mod))) && (~any(imputed_facs == fac)) 
            sprintf('[Non-imputed modality] Unexpected missing data for subject %d - %d (%s)', j, timepoint, fac_list(fac))
            untrained_imputation_counter = untrained_imputation_counter+1;
            
        % Missing imputed modality in non-trained pattern
        elseif any(isnan(sub_2impute_data(timepoint,regions_mod))) && (any(imputed_facs == fac))
            % Check if other imputed modalities are also NaN
            % If not all of the last 3 modalities were missing
            if nnz(isnan(sub_2impute_data(timepoint,regions_mod_imputed))) < numel(regions_mod_imputed)
                sprintf('[Imputed modality] Unexpected missing data for subject %d - %d (%s)',...
                    j, timepoint, fac_list(fac))
                untrained_imputation_counter = untrained_imputation_counter + 1;
            end
        end
        
    end
    end
end

%%
% List of data [sum(N_samples_per_subject), N_regs * N_facs] 
for j=1:numel(keep_ppmi_ind)
    subject = keep_ppmi_ind(j);
    
    % Imputation was optimized for FA,MD,t1/t2
    % Either re-impute the few remaining missing samples or  discard them 
    sub_imputed_data = imputed_data(gain_sub_inds == subject,:);
    sub_2impute_data = data_2impute(gain_sub_inds == subject,:);
    
    
    % Check if missing data pattern was unusual

    for timepoint=1:size(sub_2impute_data,1)
    for fac=1:N_facs
        regions_mod = (fac-1)*N_regs+1:(fac-1)*N_regs+N_regs;
        
        % Missing non-imputed modalities
        if any(isnan(sub_2impute_data(timepoint,regions_mod))) && (~any(imputed_facs == fac)) 
            sprintf('[Non-imputed modality] Unexpected missing data for subject %d - %d (%s)', j, timepoint, fac_list(fac))
            untrained_imputation_counter = untrained_imputation_counter+1;
            
        % Missing imputed modality in non-trained pattern
        elseif any(isnan(sub_2impute_data(timepoint,regions_mod))) && (any(imputed_facs == fac))
            % Check if other imputed modalities are also NaN
            % If not all of the last 3 modalities were missing
            if nnz(isnan(sub_2impute_data(timepoint,regions_mod_imputed))) < numel(regions_mod_imputed)
                sprintf('[Imputed modality] Unexpected missing data for subject %d - %d (%s)',...
                    j, timepoint, fac_list(fac))
                untrained_imputation_counter = untrained_imputation_counter + 1;
            end
        end
        
    end
    end
    
    
    subjects_background = [subjects_background; subjects_ages(j) keep_edus(j) subjects_male(j) subjects_female(j) ...
                    subjects_right_hand(j) subjects_left_hand(j) subjects_white(j) subjects_nonwhite(j)];
    
    % Taking imputed and standardized data
    %raw_mean_std = ppmi_data(subject).raw_mean_std;
    % Standardize baseline
    for fac = 1:N_facs
        regions_mod = (fac-1)*N_regs+1:(fac-1)*N_regs+N_regs;
        % Use toolbox imputation
        ppmi_data(subject).norm_data0(regions_mod,1) = (ppmi_data_imputed(subject).raw_data0(regions_mod,1)-raw_mean_std(fac,1))/raw_mean_std(fac,2);
    end    

    ppmi_data(subject).norm_data = [ppmi_data(subject).norm_data0 repmat(ppmi_data(subject).norm_data0,1,size(ppmi_data(subject).S,2)) + ppmi_data(subject).S(1:N_facs*N_regs,:)];
    ppmi_data(subject).times = min(nonzeros(ppmi_data(subject).times(:))) + [zeros(size(ppmi_data(subject).times_S,1),1) ppmi_data(subject).times_S];
    samples      = ppmi_data(subject).norm_data'; % normalized data across the whole population, including baseline
    sample_times = ppmi_data(subject).times;
    
    %samples = ppmi_data(subject).norm_data'; %ppmi_data(subject).S';
    num_subj_samples(j, 1) = size(samples, 1);
    %sample_times = ppmi_data(subject).times;%ppmi_data(subject).times_S;
    if num_subj_samples(j) < 3
        sprintf("Warning: Fewer than 3 samples for subject ind = %d (id = %d), norm_data has %d, S has %d", ...
            keep_ppmi_ind(j), keep_ids(j), size(ppmi_data(subject).norm_data,2), size(ppmi_data(subject).S,2) )
    end
    
    %sprintf("%d raw | %d S", size(ppmi_data(subject).raw_data,2), size(ppmi_data(subject).S,2))
    
    samples_to_exclude = zeros(size(samples,1),1);
    samples_to_exclude(1) = 0;
    for i=1:num_subj_samples(j,1)
        % If subject is known to require exclusion
        if any(inds_exclude_rep == j)
            samples_to_exclude(i) = 1;
            %num_subj_samples(j) = num_subj_samples(j) - 1;  
        else
        % Check if subject needs to be excluded
            if i>1
                sample_time = squeeze(sample_times(:, i));
                sample_age = min(sample_time(sample_time > 0));
                prev_sample_time = squeeze(sample_times(:,i-1));
                prev_sample_age = min(prev_sample_time(prev_sample_time > 0));
                for fac=1:N_facs
                    fac_inds = (1+((fac-1)*N_regs)):(fac*N_regs);
                    % Check if any factor or sample time stamp is identical
                    if (sum(samples(i,fac_inds) - samples(i-1, fac_inds)) == 0)...
                            ||  (sample_age==prev_sample_age)
                        samples_to_exclude(i) = 1;
                        num_subj_samples(j) = num_subj_samples(j) - 1;
                        if num_subj_samples(j) < 3
                            sprintf("Warning (repetition): Deleting subject ind = %d (id = %d), sample age diff = %.2f", ...
                                keep_ppmi_ind(j), keep_ids(j), sample_age-prev_sample_age )
                            inds_exclude_rep = [inds_exclude_rep j];
                        end
                        break
                    end
                end
            end
        end
    end
    if ~any(inds_exclude_rep == j)
        temp = [];
        
        for i=1:numel(samples_to_exclude)%num_subj_samples(j,1)
            if ~samples_to_exclude(i)
                long_S = [long_S; samples(i,:)];
                samples_subj_ids = [samples_subj_ids; subject];
                
                long_diags = [long_diags; subjects_diags(j)];
                sample_time = squeeze(sample_times(:, i));
                sample_age =  min(sample_time(sample_time > 0)); %+subjects_ages(subject);
                if (numel(long_ages)>0) && (sample_age == long_ages(numel(long_ages)))
                    sprintf("Warning sample age is identical")
                end
                long_ages = [long_ages; sample_age];
                temp = [temp sample_age];
                %sprintf("%.3f, ", sample_age)
                background = [sample_age keep_edus(j) subjects_male(j) subjects_female(j) ...
                    subjects_right_hand(j) subjects_left_hand(j) subjects_white(j) subjects_nonwhite(j)];
                
                long_background = [long_background; background];
            end
        end
        if numel(temp) > numel(unique(temp))
            sprintf("Warning: (This should not happen) Subject %d age repeated \n", j)
            temp_repeated_age_counts = temp_repeated_age_counts + 1;
        end
    end    
end


%% Exclude subjects that disqualified because of repeated data
keep_ids(inds_exclude_rep) = [];
keep_ppmi_ind(inds_exclude_rep) = [];
keep_diags(inds_exclude_rep) = [];
subjects_ages(inds_exclude_rep) = [];

subjects_background(inds_exclude_rep,:) = [];

keep_updrs_slope(inds_exclude_rep) = [];
keep_ledd_slope(inds_exclude_rep) = [];
keep_moca_slope(inds_exclude_rep) = [];
subjects_male(inds_exclude_rep) = [];
subjects_female(inds_exclude_rep) = [];
subjects_diags(inds_exclude_rep) = [];
subjects_tab_facs(inds_exclude_rep,:) = [];
subjects_ledd(inds_exclude_rep) = [];
subjects_ledd_ave(inds_exclude_rep) = [];
subjects_ledd_cum(inds_exclude_rep) = [];
subjects_ledd_slope(inds_exclude_rep) = [];
subjects_ledd_int(inds_exclude_rep) = [];
num_subj_samples(inds_exclude_rep) = [];

subjects_mods_presence = is_mod_per_sub(keep_ppmi_ind,:);

% After including imputation
sub_times = {};
for i=1:numel(keep_ppmi_ind)
   ind = keep_ppmi_ind(i);
   sub_times{i} = ppmi_data(ind).times; 
end

X0 = reshape(long_S, size(long_S,1), N_regs, N_facs);
X0 = permute(X0, [2 3 1]);

%%
close all
for i=1:size(X0,2)
    figure
    imagesc(squeeze(X0(:,i,:)))
    title(sprintf("%s", string(fac_list(i))))
    ylabel("Regions")
    xlabel("Samples from all subjects")
    colorbar
end

warning off;
X0_regrec = X0;

%% Keep or remove regions without data from Julich
N_background_covars = size(long_background, 2);
N_samples_total = size(long_ages,1);

% Exclude NeuroPM regions 60, 139 (Brodmann areas 29, 41)
% Todo: add checks
% Exclude imputed regions later
ind_regions_data = zeros(N_regs,1);

regions_have_data = zeros(N_regs,1);
regions_have_data(ind_regs_data) = 1;

if useImputed 
    ind_regions_data(1:N_regs) = 1;
else
    ind_regions_data(ind_regs_data) = 1;
end

%regs_remove = setdiff(1:N_regs, ind_regs_data);
%regs_remove = [139, 60];
regs_remove = sort([54 61 81 126 133 139], 'descend');
regs_remove = sort(ind_small_regs, 'descend');
                                
% Remove region names
reg_list_adjusted = reg_list_lr;
reg_list_adjusted(regs_remove) = [];


regs_keep = setdiff(1:N_regs, regs_remove);
N_recs = N_recs_J; % Only Julich receptors
Z0 = Z0(1:N_recs, regs_keep);
X0_regrec = X0_regrec(regs_keep, :, :);
sc = sc(regs_keep, regs_keep);
ind_regions_data = ind_regions_data(regs_keep);
regions_have_data = regions_have_data(regs_keep);
reg_ind_2_atlas = reg_ind_2_atlas(regs_keep);
ind_regs_data = find(regions_have_data);
N_regs = size(Z0, 2);

%X0_regrec = permute(X0_regrec, [2 3 1]);

% Remove covariates
warning off;
ind_healthy = find(long_diags == find(strcmp(diag_names, "Control")));
subjects_X0_adjusted = X0_regrec;

rec_list = ["AMPA", "NMDA", "kain", "musc", "flum", "cgp", "pire", ...
    "oxo", "damp", "epib", "praz", "rx", "dpat", "keta", "sch", "5-HT_{1A}", ...
    "5-HT_{1B}", "5-HT_{2A}", "5-HTT"];
rec_types(4) = "GABA";%"ACh (musc)";
rec_types(9) = "ACh (musc)";
rec_types(10) = "ACh (nico)";
rec_types(11) = "alpha 1 adrenergic";
rec_types(12) = "alpha 2 adrenergic";
rec_types(13) = "ser-5-HT1A";
rec_types(14) = "ser-5-HT2";% Ketansarin


rec_list = ["AMPA", "NMDA", "Kainate", "Muscimol", "Flumazenil", "CGP", "Pirenzepine", ...
    "Oxotremorine", "DAMP", "Epibatidine", "Prazosin", "RX781094", "8-OH-DPAT", "Ketanserin", "SCH23390", "5-HT1A", ...
    "5-HT1B", "5-HT2A", "5-HTT"];
rec_list_short = ["AMPA", "NMDA", "Kain.", "Muscimol", "Flum.", "CGP", "Pire.", ...
    "Oxo.", "DAMP", "Epib.", "Praz.", "RX", "DPAT", "Ketanserin", "SCH", "5-HT1A", ...
    "5-HT1B", "5-HT2A", "5-HTT"];
rec_list_short = ["AMPA", "NMDA", "Kainate", "GABA-A", "Bz site", "GABA-B", "M1", ...
    "M2", "M3", "α4β2", "$\alpha_1$", "$\alpha_2$", "5HT-1A", "5HT-2", "D1", "5-HT1A", ...
    "5-HT1B", "5-HT2A", "5-HTT"];
rec_list_short = ["AMPA", "NMDA", "Kainate", "GABA-A", "Bz site", "GABA-B", "M1", ...
    "M2", "M3", "α4β2", "α1", "α2", "5HT-1A", "5HT-2", "D1", "5-HT1A", ...
    "5-HT1B", "5-HT2A", "5-HTT"];

rec_types = ["Glutamate", "Glutamate", "Glutamate", "GABA", "GABA", "GABA", "ACh", "ACh", ...
    "ACh", "ACh", "Noradrenaline", "Noradrenaline", "Serotonin", "Serotonin", "Dopamine", ...
    "Serotonin", "Serotonin", "Serotonin", "Serotonin"];

N_recs = N_recs_J;
rec_list = rec_list(1:N_recs);
rec_types = rec_types(1:N_recs);
%N_recs = N_recs_J + N_recs_5ht;

%%

rec_type_names = ["Glutamate", "GABA", "Acetylcholine", "Noradrenaline", "Serotonin", "Dopamine", "Non-Receptor"];
rec_type_codes = [1 1 1 2 2 2 3 3 3 3 4 4 5 5 6 5 5 5 5];

param_names = ["offset"];
param_rec_code = [0]; % 0 for non rec
param_fac_code = [0]; % 0 for non fac
% [Fac1Rec1 Fac1Rec2 ... ]
for fac=1:N_facs
   for rec=1:N_recs
       param_names = [param_names; strcat(strcat(rec_list_short(rec), " x "), fac_list(fac))];
   end
   param_rec_code = [param_rec_code 1:N_recs];
   param_fac_code = [param_fac_code repelem(fac, N_recs)];
end
% for rec=1:N_recs
%    param_fac_code = [param_fac_code 1:N_facs]; 
% end
param_names = [param_names; string(fac_list)'; string(rec_list_short(1:N_recs))'; "spreading"];
param_rec_code = [param_rec_code repelem(0, N_facs) 1:N_recs 0];
param_fac_code = [param_fac_code 1:N_facs repelem(0, N_recs+1)];

% For permutation analysis
Z0_original = Z0;

clear -regexp adni*
clear -regexp V*
clear -regexp data_FUNC*
clear -regexp data_GM*
clear -regexp data_METB
clear -regexp summ*
clear mod
save('output/data_for_MCM.mat')
