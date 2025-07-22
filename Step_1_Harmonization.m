%% HarmCC
% Script for data harmonization on HPC cluster
% Author: Ahmed, April 2021
% Input: flattened neuroimaging data matrices with metadata
% - fac_all [n_samples, n_features]: flattened and masked neuroimaging data
% - ages_all [n_samples]
% - diags_all [n_samples]
% - gender_all [n_samples]
% - sites_all [n_samples]
% - edu_all [n_samples]
% - handed_all [n_samples]
% - subIDs_all [n_samples]

% ComBat 
% Harmonize separately for each modality using masked data
% (features x subjects matrix)
% (subjects) length site/scanner id
% (cov x subjects) covariates [age gender education handedness diagnosis]

%clearvars -except facSaveDirs DIMS   

addpath('/home/afkhan/combat/')
%addpath('/export02/data/Work/re-MCM_PD/combat/');
%addpath('/export02/data/Work/re-MCM_PD/data/flat_imaging/');
pathTempData = 'flat_constant_temp.mat';
flatFiles = dir('harm/flat_imaging/PPMI_flat_*.mat');

sprintf("Starting harmonization")

DIMS = [121 145 121];

for fac=1:numel(flatFiles)
    load(fullfile(flatFiles(fac).folder, flatFiles(fac).name));
    fac_name = erase(erase(flatFiles(fac).name, 'PPMI_flat_'), '.mat');
    
    sprintf('Harmonizing factor %d, %s', fac, fac_name)

    % 1) Remove single-subject sites or subjects with only 1 sample
    % Otherwise, combat does not converge
    singleSampleSubInds = []; % Sample ind
    singleSampleSubIDs = []; % Subject ind
    uniqueSubIDs = unique(subIDs_all);
    for indSub=1:numel(uniqueSubIDs)
        subID = find(subIDs_all == uniqueSubIDs(indSub));
        if numel(subID) == 1
            % Subjects with just one sample
            %sprintf('Single sample subject %d', indSub)
            singleSampleSubInds = [singleSampleSubInds subID];
            singleSampleSubIDs = [singleSampleSubIDs subIDs_all(subID)];
        end
    end
    
    uniqueSiteIDs = unique(sites_all);
    singleSubSiteInds = [];
    for i=1:numel(uniqueSiteIDs)
        code = uniqueSiteIDs(i);
        siteSubs = unique(subIDs_all(sites_all == code));
        
        if numel(siteSubs) == 1
            sprintf('Single subject site %d (site %d)', i, uniqueSiteIDs(i))
            
            site = find(sites_all == code);
            singleSubSiteInds = [singleSubSiteInds; site(:) ]; 
        else
            areSiteSubsSingle = zeros(numel(siteSubs),1);
            for indSub=1:numel(siteSubs)
                sub = siteSubs(indSub);
                areSiteSubsSingle(indSub) = nnz(ismember(singleSampleSubIDs,sub));
            end
            if nnz(areSiteSubsSingle) == (numel(areSiteSubsSingle) - 1)
                sprintf('All other subjects have a single sample for site %d', code)
                site = find(sites_all == code);
                singleSubSiteInds = [singleSubSiteInds; site(:) ];
            end
            
        end
    end

    single_sub_data = fac_all(singleSubSiteInds,:);
    single_sub_sites = sites_all(singleSubSiteInds);
    temp_ages = ages_all;
    
    fac_all(sort(singleSubSiteInds, 'descend'),:) = [];
    sites_all(sort(singleSubSiteInds, 'descend')) = [];
    diags_all(sort(singleSubSiteInds, 'descend')) = [];
    gender_all(sort(singleSubSiteInds, 'descend')) = [];
    edu_all(sort(singleSubSiteInds, 'descend')) = [];
    handed_all(sort(singleSubSiteInds, 'descend')) = [];
    ages_all(sort(singleSubSiteInds, 'descend')) = [];
    
    covariates = [diags_all gender_all edu_all handed_all ages_all'];
    
    % 2) Use atlas brain mask
%     masked_features = fac_all(:, indNonMask);
%     fac_all(:,indNonMask) = [];
    
    % 3) Remove constant rows and save indices
    [sds] = std(fac_all)';
    wh = find(sds==0);
    %[ns,ms] = size(wh);
    constant_features = fac_all(:,wh);
    fac_all(:,sort(wh, 'descend')) = [];

    %save(pathTempData, 'constant_features', 'masked_features', 'wh', 'single_sub_data', 'single_sub_sites', '-v7.3');
    save(pathTempData, 'constant_features', 'wh', 'single_sub_data', 'single_sub_sites', '-v7.3');
    clear constant_features
%    clear masked_features
    clear single_sub_data

    harmonized_fac = combat(fac_all', sites_all', covariates, 1)';
    clear fac_all
    
    % Add constant features
    load(pathTempData, 'constant_features', 'wh');
    if numel(constant_features)
        fac_all_constant = zeros(size(harmonized_fac,1), size(constant_features,2)+size(harmonized_fac,2));
        non_constant_ind = setdiff(1:size(constant_features,2)+size(harmonized_fac,2), wh);
        fac_all_constant(:,wh) = constant_features;
        fac_all_constant(:, non_constant_ind) = harmonized_fac;%harmonized_fac';
    else
        fac_all_constant = harmonized_fac;
    end
    clear harmonized_fac constant_features
    
    % Adjust means and add single subject sites/single sample subjects
    load(pathTempData, 'single_sub_data', 'single_sub_sites');
    groupMean = mean(fac_all_constant(:));
    subjectMean = mean(single_sub_data(:));
    if groupMean>subjectMean
        single_sub_data = single_sub_data + (groupMean-subjectMean);
    else
        single_sub_data = single_sub_data - (groupMean-subjectMean);
    end
    fac_all = zeros(size(fac_all_constant,1) + size(single_sub_data,1), size(fac_all_constant,2));
    fac_all(singleSubSiteInds,:) = single_sub_data;
    multiSubSiteInds = setdiff(1:numel(subIDs_all), singleSubSiteInds);
    fac_all(multiSubSiteInds,:) = fac_all_constant;
    clear fac_all_masked single_sub_data

    ages_all = temp_ages;

    save(sprintf('harmonized_%s', fac_name), 'fac_all', 'subIDs_all', 'ages_all', '-v7.3');
    sprintf('Saving %s', fac_name)
end
