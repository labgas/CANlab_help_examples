%% a2_set_default_options.m
%
% Set options used in various core secondlevel scripts
% --------------------------------------------------------------------
% 
% LaBGAS notes
% 
% - Always make a study-specific copy of this script in your code subdataset, do NOT edit in the repo!
% - The below can be considered LaBGAS default options, but decisions will be study-specific!

% prep_2_load_image_data_and_save &
% prep_3_calc_univariate_contrast_maps_and_save options
% --------------------------------------------------------------------
dofullplot = true;         % default true  Can set to false to save time
omit_histograms = false;     % default false Histograms not useful for large samples 
dozipimages = false;        % default false to avoid load on data upload/download when re-running often, true is useful to save space when running final analyses; lukasvo76: changed from original CANlab default

% prep_3a_run_second_level_regression_and_save options 
% --------------------------------------------------------------------
dorobust = true;            % robust statistics [true, false] -- default true
myscaling_glm = 'raw';          % 'raw' or 'scaled'; @lukasvo76 edited: change to 'scaled' if you want to use z-scored condition images, change to 'scaled_contrasts' if you want to use l2norm scaled contrast images
design_matrix_type = 'custom';   % 'group' or 'custom'
                            % Group: use DAT.BETWEENPERSON.group or 
                            % DAT.BETWEENPERSON.contrasts{c}.group;
                            % lukasvo76: choose this option if you want to
                            % compare groups without controlling for
                            % covariates
                            % Custom: use all columns of table object 
                            % DAT.BETWEENPERSON.contrasts{c};
                            % lukasvo76: choose this option if you have
                            % covariates to control for in addition to a
                            % group factor
                            
% c_univariate_contrast_maps & c2a_second_level_regression_options
%---------------------------------------------------------------------
maskname_glm = which('gray_matter_mask_sparse.img'); %lukasvo76 edited: default use of sparse gray matter mask; maskdir now defined in a_set_up_paths_always_run_first script; if you do not want to mask, change to []; if you want to use a custom mask, put it in maskdir and change name here.

% c2e_second_level_robust_parcelwise_regression
%---------------------------------------------------------------------
% see prep_3a options above as well as the following
csf_wm_covs = false; % default false, true adds global wm & csf regressors
remove_outliers = false; % default false, true removes outlier images/subjects based on mahalanobis distance 

% prep_3b_run_SVMs_on_contrasts_and_save options 
% --------------------------------------------------------------------
dosubjectnorm = false;      % default false     normalize_each_subject_by_l2norm; can help with numerical scaling and inter-subject scaling diffs
dozscoreimages = false;     % default false     Z-score each input image, removing image mean and forcing std to 1. Removes overall effects of image intensity and scale. Can be useful across studies but also removes information. Use judiciously.
dosavesvmstats = true;      % default true      Save statistics and weight map objects for SVM contrasts
dobootstrap = false;        % default false     Takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n is set to the default 10k samples
boot_n = 10000;             % default number of bootstrap samples       Reduce number for quick results
parallelstr = 'parallel';   % parallel proc for boot.   'parallel' or 'noparallel'

% prep_3c_run_SVMs_on_contrasts_masked options 
% --------------------------------------------------------------------
% see prep_3b options above as well as the following:
maskname_svm = which('gray_matter_mask_sparse.img');

% prep_3d_run_SVMs_betweenperson_contrasts options
%---------------------------------------------------------------------
% see prep_3b & prep_3c options above as well as the following:
myscaling_svm_between = 'raw'; % see above

% prep_4_apply_signatures_and_save options 
% --------------------------------------------------------------------
use_scaled_images = false; % @lukasvo76: change to true to use z-scored images - see above

% emosymp_m1_s5_predict_symptom_ratings_lasso_pcr options
% --------------------------------------------------------------------
% see prep_3b options above as well as the following:
maskname_pcr = which('gray_matter_mask_sparse.img'); % see above
myscaling_pcr = 'raw'; % options are 'raw' or 'scaled'
dosavepcrstats = true; % see above
% lukasvo76: this refers to a study-specific script for the emosymp study,
% need to work on a more generic version

% z_batch_publish_everything, z_batch_publish_analyses options 
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