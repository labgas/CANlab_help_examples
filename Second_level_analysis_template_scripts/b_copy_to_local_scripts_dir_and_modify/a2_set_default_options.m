%% a2_set_default_options.m
%
%
% USAGE
%
% Set default options used in various core secondlevel scripts included in
% CANlab_help_examples (LaBGAS fork). 
% This script is automatically called by many other scripts, including 
% a_set_up_paths_always_run_first, hence it does not often need to run
% standalone
% 
%
% NOTES
% 
% - Always make a study-specific copy of this script in your code subdataset, do NOT edit in the repo!
% - The below can be considered LaBGAS default options, but decisions will be study- or model-specific!
% - Various options added by @lukasvo76 spring 21, summer 22, and winter 23
% - If the title of the section below is capitalized, the scripts and their options have been revamped by @lukasvo76 already
%
%__________________________________________________________________________
%
% modified by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% a2_set_default_options.m         v4.2
% last modified: 2023/01/23


%% PREP_2_LOAD_IMAGE_DATA_AND_SAVE & PREP_3_CALC_UNIVARIATE_CONTRAST_MAPS_AND_SAVE
% --------------------------------------------------------------------------------

dofullplot = true;                                      % default true  Can set to false to save time
omit_histograms = false;                                % default false Histograms not useful for large samples 
dozipimages = false;                                    % default false to avoid load on data upload/download when re-running often, true is useful to save space when running final analyses; lukasvo76: changed from original CANlab default


%% PREP_3A_RUN_SECOND_LEVEL_REGRESSION_AND_SAVE
% ---------------------------------------------

maskname_glm = which('gray_matter_mask_sparse.img');    % default use of sparse gray matter mask
                                                            % model-specific maskdir defined in a_set_up_paths_always_run_first script
                                                            % if you do not want to mask, change to []
                                                            % if you want to use a custom mask, put it in maskdir and change name here
                                                            % used in this script only for visualization of unthresholded results
myscaling_glm = 'raw';                                  % 'raw', 'scaled', or 'scaledcontrasts'; 
                                                            % 'scaled': use z-scored condition images prior to computing contrasts
                                                            % 'scaled_contrasts': l2norm contrasts after computing them
design_matrix_type = 'onesample';                       % 'group', 'custom', or 'onesample'
                                                            % 'group': use DAT.BETWEENPERSON.group or 
                                                                % DAT.BETWEENPERSON.contrasts{c}.group;
                                                                % @lukasvo76: compare groups without controlling for covariates
                                                            % 'custom': use all columns of table object 
                                                                % DAT.BETWEENPERSON.contrasts{c};
                                                                % @lukasvo76: covariates to control for in addition to a group factor
                                                            % 'onesample': use intercept only
                                                                % @lukasvo76: no group factor, no covariates -
                                                                % one sample t-test with robust option
                                                                % (contrary to c_univariate_contrast_maps.m, which is being deprecated)
dorobust = true;                                        % robust statistics for voxel-based GLM [true, false] -- default true
dorobfit_parcelwise = false;                            % true runs robust parcelwise regression (CANlab's robfit_parcelwise() function) rather than voxel-based GLM (CANlab's regress() function)
    % robfit_parcelwise options
    csf_wm_covs = false;                                        % true adds global wm & csf regressors at second level
    remove_outliers = false;                                    % true removes outlier images/subjects based on mahalanobis distance 
doBayes = true;                                         % converts t-maps into Bayes Factor maps -- default true
domvpa_reg_cov = false;                                 % run MVPA regression model to predict covariate levels from (between-subject) brain data using CANlab's predict() function
    % mvpa_reg_covariate options
    algorithm_mvpa_reg_cov = 'cv_pcr';                          % default cv_pcr, will be passed into predict function (help fmri_data.predict for options)
    holdout_set_method_mvpa_reg_cov = 'no_group';               % 'no_group', or 'group'
                                                                    % 'group': use DAT.BETWEENPERSON.group or 
                                                                        % DAT.BETWEENPERSON.contrasts{c}.group;
                                                                        % @lukasvo76: balances holdout sets over groups
                                                                    % 'no_group'
                                                                        % @lukasvo76: no group factor, stratifies by
                                                                        % subject (i.e.leave whole subject out) since data is purely between-subject
    nfolds_mvpa_reg_cov = 5;                                    % default 5; number of cross-validation folds for kfold
    zscore_outcome_mvpa_reg_cov = false;                        % default false; zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models

    
%% C2A_SECOND_LEVEL_REGRESSION
% ----------------------------

% NOTE: the first option also applies to c_univariate_contrast_maps_ scripts, but the
% use of c2a_second_level_regression is preferred since it has more
% flexible options for scaling, design specification, and regression method

% GLM OPTIONS
% -----------
save_figures_glm = false;                                       % true saves .svg files of all figures generated by c2a_second_level_regression.m (slow, takes up space) % added by @lukasvo76 May 2022
q_threshold_glm = .05;                                          % threshold for FDR-corrected display items
p_threshold_glm = .005;                                         % threshold for uncorrected display items
k_threshold_glm = 10;                                           % extent threshold for both corrected and uncorrected display items
BF_threshold_glm = 10;                                          % threshold for Bayes Factor maps, |BF| > 10 indicates strong evidence in favour of H1 (positive value) or H0 (negative value) - see help.statistic_image.estimateBayesFactor for details

% MVPA OPTIONS
% ------------
dobootstrap_mvpa_reg_cov = false;                                % default false     bootstrapping; takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n is set to 10k samples
    % mvpa bootstrapping options
    boot_n_mvpa_reg_cov = 5000;                                      % default 5000      number of bootstrap samples, reduce number for quick results, increase to 10k for publication
    parallelstr_mvpa_reg_cov = 'noparallel';                         % parallel proc for boot.   'parallel' or 'noparallel'
    cons2boot = [];                                                  % default empty - bootstrap all conditions/contrasts in DAT.contrasts/conditions; specify indices for contrasts or conditions to bootstrap if you want to bootstrap a subset
    % mvpa thresholding options
    q_threshold_mvpa_reg_cov = .05;                                  % default .05       threshold for FDR-corrected display items
    k_threshold_mvpa_reg_cov = 10;                                   % default 10        extent threshold for FDR-corrected display items 

    
%% PREP_3C_RUN_SVMs_ON_CONTRASTS_MASKED 
% --------------------------------------------------------------------

ml_method_svm = 'predict';                                      % 'oofmridataobj', or 'predict'
                                                                    % 'oofmridataobj':
                                                                        % use @bogpetre's object-oriented method
                                                                        % https://github.com/canlab/ooFmriDataObjML
                                                                    % 'predict'
                                                                        % use CANlab's predict function
                                                                        % https://github.com/canlab/CanlabCore/blob/master/CanlabCore/%40fmri_data/predict.m
holdout_set_method_svm = 'onesample';                           % 'group', or 'onesample'
                                                                    % 'group': use DAT.BETWEENPERSON.group or 
                                                                        % DAT.BETWEENPERSON.contrasts{c}.group;
                                                                        % @lukasvo76: balances holdout sets over groups
                                                                    % 'onesample': use subject id only
                                                                        % @lukasvo76: no group factor, stratifies by
                                                                        % subject (i.e. leave whole subject out)
holdout_set_type_svm = 'kfold';                                 % default kfold     cross-validation method: 'kfold' or 'leave_one_subject_out' - the latter is not recommended
    nfolds_svm = 5;                                                 % default 5         number of cross-validation folds for kfold
maskname_svm = which('gray_matter_mask_sparse.img');            % default use of sparse gray matter mask; maskdir now defined in a_set_up_paths_always_run_first script; if you do not want to mask, change to []; if you want to use a custom mask, put it in maskdir and change name here.
myscaling_svm = 'raw';                                          % options are 'raw','subjectnorm','imagenorm','zscoreimages','zscorevoxels'
                                                                    % subjectnorm: normalize_each_subject_by_l2norm; normalizes images for each subject by L2 norm of Condition 1 image; can help with numerical scaling and inter-subject scaling diffs
                                                                    % imagenorm: normalize_images_by_l2norm; normalizes each image separately, not each subject/pair
                                                                    % zscoreimages: Z-score each input image, removing image mean and forcing std to 1. Removes overall effects of image intensity and scale. Can be useful across studies but also removes information. Use judiciously. lukasvo76: corresponds to 'scaled' in myscaling_glm option in prep_3a
                                                                    % zscorevoxels: Z-score each voxel across images
dosavesvmstats = true;                                          % default true      save statistics and weight map objects for SVM contrasts
dobootstrap_svm = false;                                        % default false     takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n is set to the default 10k samples
    boot_n_svm = 5000;                                              % default 5000      number of bootstrap samples, reduce number for quick results, increase to 10k for publication
parallelstr = 'parallel';                                       % parallel proc for boot.   'parallel' or 'noparallel'
dosearchlight_svm = false;                                          % default false     perform searchlight SVM analysis
    searchlight_radius_svm = 3;                                     % default 3         radius for searchlight sphere


%% C2_SVM_CONTRASTS_MASKED
% ------------------------

save_figures_svm = false;                                       % true saves .svg files of all figures generated by c2a_second_level_regression.m (slow, takes up space) % added by @lukasvo76 May 2022
q_threshold_svm = .05;                                          % threshold for FDR-corrected display items
p_threshold_svm = .005;                                         % threshold for uncorrected display items
k_threshold_svm = 10;                                           % extent threshold for both corrected and uncorrected display items


%% prep_3d_run_SVMs_betweenperson_contrasts options
% --------------------------------------------------------------------
% see prep_3b & prep_3c options above as well as the following:
myscaling_svm_between = 'raw'; % see above


%% PREP_3F_CREATE_FMRI_DATA_SINGLE_TRIAL_OBJECT
% ---------------------------------------------

cons2exclude_dat_st = {};                                       % cell array of condition names to exclude, separated by commas (or blanks)
behav_outcome_dat_st = 'rating';                                % name of outcome variable in DAT.BEHAVIOR.behavioral_data_table_st
subj_identifier_dat_st = 'participant_id';                      % name of subject identifier variable in same table
cond_identifier_dat_st = 'trial_type';                          % name of condition identifier variable in same table
% group_identifier_dat_st = 'group';                            % name of group identifier variable in same table; leave commented out if you don't have groups
vif_threshold_dat_st = 4;                                       % variance inflation threshold to exclude trials


%% C2F_RUN_MVPA_REGRESSION_SINGLE_TRIAL
% ----------------------------------------------------------------------------------------------------------------------

% GENERAL OPTIONS
%----------------
ml_method_mvpa_reg_st = 'predict';                              % 'oofmridataobj', or 'predict'
                                                                    % 'oofmridataobj':
                                                                    % use @bogpetre's object-oriented method
                                                                    % https://github.com/canlab/ooFmriDataObjML
                                                                    % 'predict'
                                                                    % use CANlab's predict function
                                                                    % https://github.com/canlab/CanlabCore/blob/master/CanlabCore/%40fmri_data/predict.m
algorithm_mvpa_reg_st = 'cv_pcr';                               % default cv_pcr, will be passed into predict function (help fmri_data.predict for options) if ml_method_mvpa_reg_st == 'predict' or adapted if 'oofmridataobj'
                                                                    % if ml_method_mvpa_reg_st = 'oofmridataobj', only cv_pls and cv_pcr are implemented in this script for now
holdout_set_method_mvpa_reg_st = 'onesample';                   % 'group', or 'onesample'
                                                                    % 'group': use DAT.BETWEENPERSON.group or 
                                                                    % DAT.BETWEENPERSON.contrasts{c}.group;
                                                                    % @lukasvo76: balances holdout sets over groups
                                                                    % 'onesample': use subject id only
                                                                    % @lukasvo76: no group factor, stratifies by
                                                                    % subject (i.e. leave whole subject out)
nfolds_mvpa_reg_st = 5;                                         % default 5; number of cross-validation folds for kfold
zscore_outcome_mvpa_reg_st = false;                             % default false; zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models
maskname_mvpa_reg_st = which('gray_matter_mask_sparse.img');    % see above
myscaling_mvpa_reg_st = 'raw';                                  % options are 'raw', 'centerimages', 'zscoreimages', 'l2normimages', 'zscorevoxels'

% STATISTICS AND RESULTS VISUALIZATION OPTIONS
%---------------------------------------------
dobootstrap_mvpa_reg_st = false;                                % default false     bootstrapping; takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n is set to 10k samples
    boot_n_mvpa_reg_st = 5000;                                      % default 5000      number of bootstrap samples, reduce number for quick results, increase to 10k for publication
    parallelstr_mvpa_reg_st = 'parallel';                           % parallel proc for boot.   'parallel' or 'noparallel'
doperm_mvpa_reg_st = false;                                     % default false     permutation testing; takes a very long time despite implementation of parfor loop
    perm_n_mvpa_reg_st = 5000;                                      % default 5000      number of permutations, reduce number for quick results, increase to 10k for publication
    perm_sidedness = 'both';                                        % default both      tails for permutation test, 'both','smaller', or 'larger'
dosourcerecon_mvpa_reg_st = false;                              % default false     source reconstruction/"structure coefficients", i.e. regressing each voxel's activity onto yhat - see Haufe et al NeuroImage 2014
    dosourcerecon_perm_mvpa_reg_st = false;                         % default false     permutation testing on source recon images; takes a very long time despite implementation of parfor loop
q_threshold_mvpa_reg_st = .05;                                  % default .05       threshold for FDR-corrected display items
k_threshold_mvpa_reg_st = 10;                                   % default 10        extent threshold for FDR-corrected display items 
dosavemvparegstats = true;                                      % see saving options above

domultilevel_mvpa_reg_st = false;                               % default false; fits multilevel mvpa models - WORK IN PROGRESS


%% C2G_RUN_MULTILEVEL_MEDIATION_SINGLE_TRIAL
% ----------------------------------------------------------------------------------------------------------------------

% GENERAL OPTIONS
%----------------

save_figures_pdm = false;                                       % default false; % true saves .svg files of all figures (slow, takes up space)
zscore_outcome_pdm = false;                                     % default false; zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models
maskname_pdm = which('gray_matter_mask_sparse.img');            % see above
myscaling_pdm = 'raw';                                          % options are 'raw', 'centerimages', 'zscoreimages', 'l2normimages', 'zscorevoxels'

% STATISTICS AND RESULTS VISUALIZATION OPTIONS
%---------------------------------------------
nPDM = 10;                                                      % default 10        number of PDMs to retain, chances are very low that meaningful variance is explained by PDM # > 10
dobootstrap_pdm = false;                                        % default false     bootstrapping; takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n is set to 10k samples
    boot_n_pdm = 5000;                                              % default 5000      number of bootstrap samples, reduce number for quick results, increase to 10k for publication
    k_threshold_pdm = 10;                                           % default 10        extent threshold for bootstrapped fdr-corrected results
dosourcerecon_pdm = false;                                      % default false     source reconstruction/"structure coefficients", i.e. regressing each voxel's activity onto yhat - see Haufe et al NeuroImage 2014
dosavepdmstats = true;                                          % see saving options above


%% PREP_4_APPLY_SIGNATURES_AND_SAVE 
% --------------------------------------------------------------------
myscaling_sigs = 'raw';                                         % default 'raw'             'raw', or 'scaled', see myscaling_glm above
similarity_metric_sigs = 'dotproduct';                          % default 'dotproduct'      other options 'cosine_similarity','correlation' - passed into apply_all_signatures
keyword_sigs = 'all';                                           % default 'all'             passed into load_image_set, help load_image_set for overview of many options


%% z_batch_publish_everything, z_batch_publish_analyses options 
% --------------------------------------------------------------------
% Enter string for which analyses to run, in any order
%
% 'contrasts'     : Coverage and univariate condition and contrast maps
% 'signatures'    : Pre-defined brain 'signature' responses from CANlab
% 'svm'           : Cross-validated Support Vector Machine analyses for each contrast
% 'bucknerlab'    : Decomposition of each condition and contrast into loadings on resting-state network maps
% 'meta_analysis' : Tests of "pattern of interest" analyses and ROIs derived from CANlab meta-analyses
%
% Default if you run z_batch_publish_analyses is to run all, in this order.
% Or run a custom set:
% batch_analyses_to_run = {'contrasts' 'signatures' 'svm' 'bucknerlab' 'meta_analysis'};
% z_batch_publish_analyses({'svm' 'bucknerlab'})