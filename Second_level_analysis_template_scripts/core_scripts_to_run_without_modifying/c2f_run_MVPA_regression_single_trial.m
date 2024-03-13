%% c2f_run_MVPA_regression_single_trial.m
%
%
% USAGE
% -----
%
% This script runs MVPA regression analysis on a continuous outcome Y
% (default pcr, but can easily be adapted to pls or
% other machine learning algorithms) on an fmri_data_st
% object created using prep_3f_create_fmri_data_single_trial_object.m
% That script should be run first, or the present script will load the data
% object if it is saved by the previous script.
%
% Options for this script are set in a2_set_default_options.m, see that
% script and below for more info. 
% Many of these options get passed into CANlab's predict function:
% help fmri_data.predict in Matlab command window for more info
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('c2f_run_MVPA_regression_single_trial','outputDir',htmlsavedir)
%
%
% TUTORIALS AND DOCUMENTATION
% ---------------------------
%
% A. CANLAB'S PREDICT FUNCTION
%
% This is the classic CANlab method of running ML models, and can be chosen
% by setting the ml_method_mvpa_reg_st option in a2_set_default_options.m
% to 'predict'
%
% This script is based on the extremely helpful tutorials on 
% single trial MVPA analysis by @bogpetre @CANlab.
%
% Here are Bogdan's walkthroughs:
% https://canlab.github.io/_pages/canlab_single_trials_demo/demo_norming_comparison.html
% https://canlab.github.io/_pages/mlpcr_demo/mlpcr_demo.html (WiP)
%
% Here are two scripts @lukasvo76 adapted from these walkthroughs
% https://www.dropbox.com/sh/e17nl3ew1db1twk/AACO9QAEt6Sy3TejH-n-tbdEa?dl=0
% https://www.dropbox.com/sh/bm0at2dr81isk70/AABD67D_bF8A0NFa4gtt2dHNa?dl=0
% 
% Another highly helpful resource in this context is this Nature Methods
% paper by Tor and Wani Woo
% https://www.nature.com/articles/s41596-019-0289-5
%
% @lukasvo76's version of the script for this paper can be found here
% https://www.dropbox.com/sh/v2nsgoqmbi0cqnk/AAD6I1Gn5KUM6aViom4TLeVJa?dl=0
%
% B. BOGDAN'S MACHINE LEARNING TOOLKIT FOR FMRI_DATA OBJECTS
%
% This is a newer method inspired by Python's scikit-learn, including more
% flexible options for algorithm and feature selection, 
% hyperparameter optimization, nested cross-validation, etc. However, it
% does require more advanced programming skills and understanding the logic
% of the method
%
% Dependency: https://github.com/canlab/ooFmriDataObjML
%
% Tutorial: https://canlab.github.io/_pages/canlab_pipelines_walkthrough/estimateBestRegionPerformance.html
% Example script: https://github.com/labgas/LaBGAScore/blob/main/secondlevel/LaBGAScore_secondlevel_ooFmriDataObjML_example.m
% 
%
% OPTIONS
% -------
%
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options, 
% you can make a copy of this script with a letter index (e.g. _s6a_) 
% and change the default options below
%
% GENERAL OPTIONS
%
% ml_method_mvpa_reg_st: 'oofmridataobj', or 'predict'
%       'oofmridataobj':
%           use @bogpetre's object-oriented method
%           https://github.com/canlab/ooFmriDataObjML
%       'predict'
%           use CANlab's predict function
%           https://github.com/canlab/CanlabCore/blob/master/CanlabCore/%40fmri_data/predict.m
% algorithm_mvpa_reg_st: default cv_pcr
%       will be passed into predict function (help predict for options) if ml_method_mvpa_reg_st == 'predict' or adapted correctly if 'oofmridataobj'
%       if ml_method_mvpa_reg_st = 'oofmridataobj', only cv_pls and cv_pcr are implemented in this script for now
% holdout_set_method_mvpa_reg_st: 'group', or 'onesample'
%       'group': use DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.contrasts{c}.group;
%           balances holdout sets over groups
%       'onesample': use subject id only
%           no group factor, stratifies by subject (i.e. leave whole subject out)
% nfolds_mvpa_reg_st: default 5; number of cross-validation folds for kfold
% zscore_outcome_mvpa_reg_st: default false; true zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models
% maskname_mvpa_reg_st: default which('gray_matter_mask_sparse.img');
%       - default use of sparse gray matter mask
%       - maskdir now defined in a_set_up_paths_always_run_first script
%       - if you do not want to mask, change to []
%       - if you want to use a custom mask, put it in maskdir and change name here.
% myscaling_mvpa_reg_st: default 'raw'; options are 'raw', 'centerimages', 'zscoreimages', 'l2normimages', 'zscorevoxels'
%
% STATISTICS AND RESULTS VISUALIZATION OPTIONS
%
% dobootstrap_mvpa_reg_st: default false; true bootstraps weights - takes AN AWFUL LOT OF TIME, hence only use true for final analysis
%    boot_n_mvpa_reg_st: default 5000; number of bootstrap samples, reduce number for quick results
%    parallelstr_mvpa_reg_st: default 'parallel'; parallel proc for bootstrapping
% doperm_mvpa_reg_st: default false; true performs permutation testing - takes AN AFWUL LOT OF TIME, hence only use true for final analysis
%    perm_n_mvpa_reg_st: default 5000; number of permutations, reduce number for quick results
%    perm_sidedness: default 'both'; tails for permutation test, 'both','smaller', or 'larger'
% dosourcerecon_mvpa_reg_st: default false; true performs source reconstruction/"structure coefficients", i.e. regressing each voxel's activity onto yhat - see Haufe et al NeuroImage 2014
%    dosourcerecon_perm_mvpa_reg_st: default false; true performs
%    permutation testing on source recon images; takes an AWFUL LOT OF TIME
% q_threshold_mvpa_reg_st: default .05; threshold for FDR-corrected display items
% k_threshold_mvpa_reg_st: default = true;                                      % see saving options above 10; extent threshold for FDR-corrected display items 
% dosavemvparegstats: default true; Save statistics and weight map objects
% domultilevel_mvpa_reg_st: default false; fits multilevel mvpa models - WORK IN PROGRESS
%
%
% IMPORTANT NOTE
% --------------
%
% This script is work in progress, particularly the permutation and multilevel
% options are still undergoing improvement and full testing
% lukasvo76 note to self: see Bogdan's hyp_opt_and_mlpcr_demo script in
% this repo to implement the latter!
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be, bogpetre@gmail.com
% date:   April, 2021
%__________________________________________________________________________
% @(#)% c2f_run_MVPA_regression_single_trial     v6.0        
% last modified: 2024/02/28


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% COPY MANDATORY OPTION FROM CORRESPONDING PREP_3f_ SCRIPT IF

results_suffix = ''; % suffix of your choice added to .mat file with saved results

% COPY OPTIONS FROM CORRESPONDING PREP_3f_ SCRIPT

% cons2exclude_dat_st = {}; % cell array of condition names to exclude, separated by commas (or blanks)
% behav_outcome_dat_st = 'rating'; % name of outcome variable in DAT.BEHAVIOR.behavioral_data_table_st
% subj_identifier_dat_st = 'participant_id'; % name of subject identifier variable in same table
% group_identifier_dat_st = 'group'; % name of group identifier variable in same table; leave commented out if you don't have groups

% SET CUSTOM OPTIONS FOR CURRENT SCRIPT
 
% See documentation above and a2_set_default_options.m for list of options

% NOTE: the latter two option categories only need to be specified if you want to run a second version of your model 
% with different options than the defaults you set in your model-specific version of a2_set_default_options.m


%% LOAD FMRI_DATA_ST OBJECT AND OTHER NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('LOADING DATA');
fprintf('\n\n');

if ~exist('DSGN','var') || ~exist('DAT','var')

    load(fullfile(resultsdir,'image_names_and_setup.mat'));

    if ~isfield(DAT,'BEHAVIOR')
        fprintf('\n');
        error('Behavioral data not yet added to DAT structure - run prep_1b script first');
    end

end

if ~exist('fmri_dat','var')
    
    if ~isempty(cons2exclude_dat_st)
        load(fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '.mat']));
        
    else
        load(fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome_dat_st, '_', results_suffix, '.mat']));

    end
end


%% DEFINE SUBJECT IDENTIFIERS
% -------------------------------------------------------------------------

subject_id = fmri_dat.metadata_table.(subj_identifier_dat_st);
[uniq_subject_id, ~, subject_id] = unique(subject_id,'stable');
% fmri_dat.metadata_table.subject_id = subject_id; % not needed but does not harm
n_subj = size(uniq_subject_id,1);


%% SCALE AND/OR MASK IMAGES AND BEHAVIORAL OUTCOME ACCORDING TO OPTIONS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('MASKING AND SCALING IMAGES IF REQUESTED IN OPTIONS');
fprintf('\n\n');

% MASKING IMAGES
%---------------

    if exist('maskname_mvpa_reg_st','var') && ~isempty(maskname_mvpa_reg_st) && exist(maskname_mvpa_reg_st, 'file')

        [~,maskname_short] = fileparts(maskname_mvpa_reg_st);
        fprintf('\nMasking data with %s\n',maskname_short);
        mask_string = sprintf('masked with %s', maskname_short);

        mvpamask = fmri_mask_image(maskname_mvpa_reg_st);
        mvpamask = resample_space(mvpamask,fmri_dat); % resample mask to fmri_data obj space
        fmri_dat = fmri_dat.apply_mask(mvpamask);
        fmri_dat.mask_descrip = maskname_mvpa_reg_st;

    else

        fprintf('\nNo mask found; using full original image data\n\n');
        mask_string = sprintf('without_masking');

    end % if loop mask


% SCALING IMAGES
%---------------

    switch myscaling_mvpa_reg_st

        case 'raw'

            fprintf('\nNo scaling of input images\n\n');

        case 'centerimages'

            fmri_dat = fmri_dat.rescale('centerimages'); 
            fprintf('\nCentering input images\n\n');

        case 'l2normimages'

            fmri_dat = fmri_dat.rescale('l2norm_images');
            fprintf('\nNormalizing input images by l2norm\n\n');

        case 'zscoreimages'

            fmri_dat = fmri_dat.rescale('zscoreimages');
            fprintf('\nZ-scoring input images\n\n');

        case 'zscorevoxels'

            fmri_dat = fmri_dat.rescale('zscorevoxels');
            fprintf('\nZ-scoring voxels across input images\n\n');

        otherwise 

            error('\ninvalid scaling option %s specified in myscaling_mvpa_reg_st variable defined in a2_set_default_options script, please correct\n\n', myscaling_mvpa_reg_st);

    end % switch scaling
   

% ZSCORE BEHAVIORAL OUTCOME
%--------------------------

% NOTE: useful for more interpretable values of prediction MSE

    if zscore_outcome_mvpa_reg_st

        fmri_dat.Y = zscore(fmri_dat.Y);
        fprintf('\nZ-scoring outcome fmri_dat.Y across subjects\n\n');

    end


%% DATA VISUALISATION PRIOR TO MODEL BUILDING
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('PLOTTING DATA');
fprintf('\n\n');

% BETA IMAGES
%------------

h1=figure;

    for sub = 1:n_subj
        subj_idx = sub == subject_id;
        this_subj_dat = fmri_dat.dat(:,subj_idx);
        q(sub,:) = quantile(this_subj_dat(:),[0.025,0.5,0.975]);
        mu = mean(mean(this_subj_dat(:)));
        sd = std(this_subj_dat(:));
        h1 = plot([mu-sd, mu+sd],[sub,sub],'-');
        hold on;
        h2 = plot(mu,sub,'o');
        h2.Color = h1.Color;
    end

box off
title('Distribution of beta weights');
xlabel('\beta');
ylabel('Subject');
hold off

p = get(gcf,'Position');
set(gcf,'Position',[p(1:2),1024,2048],'WindowState','Maximized');
drawnow, snapnow;

clear sub

% BEHAVIORAL OUTCOME
%-------------------

% OVER SUBJECTS

b1=figure;

hold off;
b1=histogram(fmri_dat.Y);
box off
title(['Histogram of single trial ' behav_outcome_dat_st]);
xlabel(behav_outcome_dat_st);
ylabel('n(observations)');
set(gcf,'WindowState','Maximized');
drawnow, snapnow;

% PER SUBJECT

b2=figure;

    for sub = 1:n_subj
        this_idx_Y = find(sub == subject_id);
        this_Y = fmri_dat.Y(this_idx_Y);

        subplot(ceil(sqrt(n_subj)), ceil(n_subj/ceil(sqrt(n_subj))), sub);
        hold off
        b2 = histogram(this_Y);
        box off
        title(uniq_subject_id{sub});
        xlabel(behav_outcome_dat_st);
        ylabel('n(obs)');
    end

set(gcf,'WindowState','Maximized');
drawnow, snapnow;

clear sub

        
%% FIT SINGLE-LEVEL MVPA MODELS
% -------------------------------------------------------------------------

% NOTE: algorithm choice set in a2_default_options.m

fprintf('\n\n');
printhdr(['RUNNING ', upper(algorithm_mvpa_reg_st)]);
fprintf('\n\n');

% RUN PREDICTIVE REGRESSION MODEL
%--------------------------------
    switch ml_method_mvpa_reg_st

    % CANLAB's PREDICT FUNCTION

        case 'predict'

            % CROSS-VALIDATION FOLD SELECTION

            % NOTE: balancing over groups, stratifying over subjects (i.e. leave whole
            % subject out)

            switch holdout_set_method_mvpa_reg_st

                case 'group'

                    group = fmri_dat.metadata_table.(group_identifier_dat_st);
                    cv = cvpartition2(group, 'Group',subject_id, 'GroupKFold', nfolds_mvpa_reg_st);
                        fold_labels = zeros(size(fmri_dat.dat,2),1);
                        for sub = 1:cv.NumTestSets
                            fold_labels(cv.test(sub)) = sub;
                        end

                case 'onesample'

                    cv = cvpartition2(size(fmri_dat.dat,2),'Group',subject_id, 'GroupKFold', nfolds_mvpa_reg_st);
                        fold_labels = zeros(size(fmri_dat.dat,2),1);
                        for sub = 1:cv.NumTestSets
                            fold_labels(cv.test(sub)) = sub;
                        end

                otherwise

                    error('\ninvalid option "%s" defined in holdout_set_method_mvpa_reg_st variable, choose between "group" and "onesample"\n\n',holdout_set_method_mvpa_reg_st);

            end % switch holdout_set_method

            % FIT MODEL

            t0 = tic;

            [cverr, stats, optout] = predict(fmri_dat, 'algorithm_name', algorithm_mvpa_reg_st, ...
                        'nfolds', fold_labels, 'error_type', 'mse', parallelstr_mvpa_reg_st, 'verbose', 0);

            t_end = toc(t0);        

            
    % BOGDAN'S OOFMRIDATAOBJ METHOD

        case 'oofmridataobj'

            % DEFINE ALGORITHM

            switch algorithm_mvpa_reg_st

                case 'cv_pls'
                    alg = plsRegressor(); % intiate alg as a plsRegressor estimator object, other estimators in Github repo/estimators; ; numcomponents is arbitrary, but typically low for pls - we will optimize this hyperparm later

                case 'cv_pcr'
                    alg = pcrRegressor();

                otherwise
                    error('\nchoice of algorithm "%s" in algorithm_mvpa_reg_st option variable is not compatible with choice of machine learning method "%s" in ml_method_mvpa_reg_st option variable,\n Either chance method to "predict" or change algorithm to "cv_pcr" or "cv_pls"\n\n', algorithm_mvpa_reg_st, ml_method_mvpa_reg_st);

            end

            alg.fit(fmri_dat.dat', fmri_dat.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
            % NOTE: fit is not strictly necessary at this stage, but a good test
            alg_params = alg.get_params; % get to know the hyperparams for this algorithm, which we want to optimize

            % DEFINE FEATURE EXTRACTOR
            
            featConstructor_han = @(X)([]);
            extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han); % initiate extractVxl as an empty fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
            extractVxl.fit(fmri_dat); % transformer takes fmri_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)
            % NOTE: fit is not strictly necessary at this stage, but a good test

            % DEFINE PIPELINE

            fmri_pipeline = pipeline({{'featExt',extractVxl},{'alg',alg}}); % define fmri_pcr as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
            fmri_pipeline.fit(fmri_dat,fmri_dat.Y);
            % NOTE: fit is not strictly necessary at this stage, but a good test

            % INNER CROSS-VALIDATION FUNCTION

            switch holdout_set_method_mvpa_reg_st

                case 'group'
                    innercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier_dat_st), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st));

                case 'onesample'
                    innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

            end
            % NOTE: we use metadata_table here, since the input to bo.fit is the
            % fmri_data_st object fmri_dat, which has a metadata_table field

            cv_folds_mvpa_reg_st = innercv(fmri_dat,fmri_dat.Y); 
            % NOTE: get cross-validation folds, not strictly necessary, but a good test

            % DEFINE BAYESIAN OPTIMIZATION

            switch algorithm_mvpa_reg_st

                case 'cv_pls'
                    dims = optimizableVariable('alg__numcomponents',[1,30],'Type','integer','Transform','log');

                case 'cv_pcr'
                    dims = optimizableVariable('alg__numcomponents',[1,floor(rank(fmri_dat.dat)*(nfolds_mvpa_reg_st-1)/nfolds_mvpa_reg_st)],'Type','integer');

            end

            % NOTE: Type and number of hyperparams to optimize depends on algorithm (check alg.get_params above), as well as other settings

            bayesOptParams = {dims, 'AcquisitionFunctionName','expected-improvement-plus',...
                'MaxObjectiveEvaluations',30, 'UseParallel', false, 'verbose',1, 'PlotFcn', {}};

            bo = bayesOptCV(fmri_pipeline,innercv,@get_mse,bayesOptParams);
            bo.fit(fmri_dat,fmri_dat.Y);
            bo_numcomponents = bo.estimator.estimator.numcomponents;

            % OUTER CROSS-VALIDATION FUNCTION

            switch holdout_set_method_mvpa_reg_st

                case 'group'
                    outercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier_dat_st), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st));

                case 'onesample'
                    outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

            end
            % NOTE: we use metadata_table here, since the input to cvGS.do is the
            % fmri_data_st object fmri_dat, which has a metadata_table field

            % ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE

            cvGS = crossValScore(bo, outercv, @get_mse, 'n_parallel', nfolds_mvpa_reg_st, 'verbose', true);
            % NOTE: Bogdan advises not parallizing too much for the purpose of Bayesian
            % model optimization, since each step learns from the previous one, so we
            % only parallelize the outer cv loop with 1 core per outer cv fold,
            % resulting in very acceptable runtimes (on LaBGAS server with 128 GB RAM)

            cvGS.do(fmri_dat, fmri_dat.Y);
            cvGS.do_null(); % fits null model - intercept only
            fold_labels = cvGS.fold_lbls;

            % CREATE AN FMRI_DATA OBJECT WITH THE BETAS FOR VISUALIZATION PURPOSES
            
            weight_obj = bo.estimator.transformers{1}.brainModel; % empty .dat at this stage
            weight_obj.dat = bo.estimator.estimator.B(:); % fills mdl.dat with betas

        otherwise

            error('\ninvalid option "%s" defined in ml_method_mvpa_reg_st variable, choose between "oofmridataobj" and "predict"\n\n',ml_method_mvpa_reg_st);

    end % switch machine learning method


% VISUALIZE RESULTS
%------------------

fprintf('\n\n');
printhdr('VISUALIZING RESULTS');
fprintf('\n\n');

% PLOT OBSERVED VERSUS PREDICTED

fprintf('\n\n');
printhdr('Plotting observed versus predicted');
fprintf('\n\n');

    switch ml_method_mvpa_reg_st

        case 'predict'

            fprintf('\n%s r = %0.3f\n\n', algorithm_mvpa_reg_st, corr(stats.yfit, fmri_dat.Y));

            figure
            
            nr_subjs = max(unique(subject_id));
            
            colors = distinguishable_colors(nr_subjs);
            colors_cell = mat2cell(colors,ones(1,nr_subjs));
            [han,~,~,slope_stats] = line_plot_multisubject(fmri_dat.Y, stats.yfit, 'subjid', subject_id, 'group_avg_ref_line', 'MarkerTypes','d', 'colors', colors_cell);
            xlabel({['Observed ' behav_outcome_dat_st],'(across conditions)'},'FontSize',14); 
            ylabel({['Estimated ' behav_outcome_dat_st],'(cross validated)'},'FontSize',14);

            for i = 1:size(han.line_handles,2)
                han.line_handles(1,i).LineWidth = 2.5;
                han.point_handles(1,i).MarkerSize = 5;
            end
            
            han.grpline_handle.LineWidth = 7.5;
            han.grpline_handle.LineStyle = '-';
            han.grpline_handle.Color = 'r';
            
            ax = gca;
            ax.FontWeight = 'bold';

            set(gcf,'WindowState','Maximized','Color','w');
            drawnow, snapnow;


        case 'oofmridataobj'

            for y = 1:size(cvGS.Y,2)
                r(y,1) = corr(cvGS.yfit{y},cvGS.Y{y});
            end

            cv_r = mean(r); % average r(yfit,y) over folds = cv correlation

            cv_mse = mean(cvGS.scores); % average MSE over folds = cv model performance

            fprintf('\n%s cross-validated r = %0.3f\n\n', algorithm_mvpa_reg_st, cv_r);
            fprintf('\n%s cross-validated mse = %0.3f\n\n', algorithm_mvpa_reg_st, cv_mse);

            f1 = cvGS.plot; % plots predicted versus observed

            set(gcf,'WindowState','Maximized');
            drawnow, snapnow;

    end


% PLOT MONTAGE OF UNTHRESHOLDED WEIGHTS

fprintf('\n\n');
printhdr('Plotting unthresholded weight maps');
fprintf('\n\n');

    whmontage = 5;

    fprintf ('\nSHOWING UNTHRESHOLDED %s RESULTS, %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

    figure

    o2 = canlab_results_fmridisplay([], 'compact', 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        switch ml_method_mvpa_reg_st

            case 'predict'
                w = region(stats.weight_obj);

            case 'oofmridataobj'
                w = region(weight_obj);

        end

    o2 = addblobs(o2, w, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' unthresholded ' mask_string]);

    figtitle = sprintf('%s_unthresholded_montage_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;

    clear w, clear o2, clear figtitle


%% PERFORM BOOTSTRAPPING, AND/OR PERMUTATION TESTING ON SINGLE-LEVEL MVPA RESULTS IF REQUESTED
% --------------------------------------------------------------------------------------------

% BOOTSTRAP IF REQUESTED
%-----------------------

    if dobootstrap_mvpa_reg_st
        
        delete(gcp('nocreate'));
        c = parcluster('local'); % determine local number of cores, and initiate parallel pool with 80% of them
        nw = c.NumWorkers;
        parpool(round(0.8*nw));
        
        fprintf('\n\n');
        printhdr('BOOTSTRAPPING WEIGHT MAPS');
        fprintf('\n\n');

        t0_boot = tic;

        switch ml_method_mvpa_reg_st
    
            case 'predict'
                [~ , bs_stats] = predict(fmri_dat, 'algorithm_name', algorithm_mvpa_reg_st,...
                    'bootsamples', boot_n_mvpa_reg_st, 'nfolds', 1, 'error_type', 'mse', ...
                    parallelstr_mvpa_reg_st, 'verbose', 0);
                
            case 'oofmridataobj'
                [~ , bs_stats] = predict(fmri_dat, 'algorithm_name', algorithm_mvpa_reg_st,...
                    'bootsamples', boot_n_mvpa_reg_st, 'nfolds', 1, 'numcomponents', bo_numcomponents, ...
                    'error_type', 'mse', parallelstr_mvpa_reg_st, 'verbose', 0);
                
        end

        t_end_boot = toc(t0_boot);

    elseif ~dobootstrap_mvpa_reg_st && doperm_mvpa_reg_st % bootstrap with few samples just to create a statistic_image weight object for use in permutation code below

        t0_boot = tic;

        [~ , bs_stats] = predict(fmri_dat, 'algorithm_name', algorithm_mvpa_reg_st,...
            'bootsamples', 20, 'nfolds', 1 , 'error_type', 'mse', ...
            parallelstr_mvpa_reg_st, 'verbose', 0);

        t_end_boot = toc(t0_boot);

    end

    
% PERFORM PERMUTATION TESTING IF REQUESTED
%----------------------------------------

% NOTE: code adapted from @pkragel s1_predict_anxiety_ratings_lassopcr.m -
% NEEDS CHECK BY PHIL!

    if doperm_mvpa_reg_st
        
        fprintf('\n\n');
        printhdr('PERFORMING PERMUTATION TESTING ON WEIGHT MAPS');
        fprintf('\n\n');

        t0_perm = tic;

        % OBTAIN PERMUTATION WEIGHTS

        null_beta = zeros(perm_n_mvpa_reg_st,size(fmri_dat.dat,1)); % number of permutations, number of voxels
        null_weights = zeros(perm_n_mvpa_reg_st,size(fmri_dat.dat,1));
        size_dat = size(fmri_dat.dat,1); % more efficient in parfor loop
        size_Y = size(fmri_dat.Y,1);

            parfor perm = 1:perm_n_mvpa_reg_st

                regress_stats_null = cell(max(fold_labels),1);
                betas_null = zeros(max(fold_labels),size_dat);

                random_inds=randperm(size_Y); % number of images/ratings over subjects
                temp_dat=fmri_dat;
                temp_dat.Y=temp_dat.Y(random_inds);

                switch ml_method_mvpa_reg_st

                    case 'predict'
                        [~, stats_null] = predict(temp_dat, 'algorithm_name', algorithm_mvpa_reg_st, 'nfolds', fold_labels, ...
                            'error_type', 'mse', parallelstr_mvpa_reg_st, 'verbose', 0);

                    case 'oofmridataobj'
                        [~, stats_null] = predict(temp_dat, 'algorithm_name', algorithm_mvpa_reg_st, 'nfolds', fold_labels, ...
                            'numcomponents', bo_numcomponents, 'error_type', 'mse', parallelstr_mvpa_reg_st, 'verbose', 0);

                end

                for k = 1:max(fold_labels) % number of cv folds
                    regress_data_null = fmri_dat;
                    regress_data_null.X = stats_null.yfit(fold_labels==k);
                    regress_data_null.dat = regress_data_null.dat(:,fold_labels==k);
                    regress_stats_null{k} = regress(regress_data_null,'nodisplay','noverbose');
            %                         tv = replace_empty(regress_stats(k).b);
                    betas_null(k,:) = regress_stats_null{k}.b.dat(:,1);
                end

                null_beta(perm,:) = mean(betas_null);
                null_weights(perm,:) = stats_null.weight_obj.dat;

                perm/perm_n_mvpa_reg_st

            end % parfor loop over number of permutations

        % GET PROBABILITY FOR OBSERVED WEIGHTS FROM PERMUTATON WEIGHTS

            switch ml_method_mvpa_reg_st

                case 'predict'
                    mean_weight = stats.weight_obj.dat;

                case 'oofmridataobj'
                    mean_weight = weight_obj.dat;

            end

            for weight = 1:size(null_weights,2)

                phat(weight) = 1-sum(mean_weight(weight) > null_weights(:,weight)) / (1+size(null_weights,1));

                switch perm_sidedness
                    case 'both'
                        phat(weight) = (length(find(abs(null_weights(:,weight)) > abs(mean_weight(weight))))+1) / (perm_n_mvpa_reg_st+1); 
                    case 'smaller'
                        phat(weight) = (length(find(null_weights(:,weight) < mean_weight(weight)))+1) / (perm_n_mvpa_reg_st+1);
                    case 'larger'
                        phat(weight) = (length(find(null_weights(:,weight) > mean_weight(weight)))+1) / (perm_n_mvpa_reg_st+1);
                end

            end

        clear weight;

        perm_stats_obj = bs_stats.weight_obj; % initiate perm_stats_obj as a statistic_image object identical to bs_stats.weight_obj - is easiest from existing object, but will not work if bootstrapping not (yet) performed
        perm_stats_obj.dat = mean_weight; % Store weights and p-values in statistic_image object
        perm_stats_obj.p = phat';

        t_end_perm = toc(t0_perm);
    
    end % if loop permutation


% VISUALIZE BOOTSTRAPPING, AND/OR PERMUTATION RESULTS
%----------------------------------------------------

    if dobootstrap_mvpa_reg_st || doperm_mvpa_reg_st
        
        fprintf('\n\n');
        printhdr('VISUALIZING BOOTSTRAPPING AND/OR PERMUTATION RESULTS');
        fprintf('\n\n');
        
    end

    % PLOT MONTAGE OF THRESHOLDED WEIGHTS AFTER BOOTSTRAPPING AND/OR PERMUTATION TESTING

    if dobootstrap_mvpa_reg_st

        fprintf ('\nSHOWING BOOTSTRAPPED %s RESULTS, %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

        figure

        o2 = canlab_results_fmridisplay([], 'compact', 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        t = bs_stats.weight_obj;
        t = threshold(t, q_threshold_mvpa_reg_st, 'fdr', 'k', k_threshold_mvpa_reg_st); 
        r = region(t);

        o2 = addblobs(o2, r, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' bootstrapped ' mask_string]);

        figtitle = sprintf('%s_%1.4f_bootstrap_FDR_montage_%s_%s', algorithm_mvpa_reg_st, q_threshold_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

        clear w, clear o2, clear figtitle

    end % if loop bootstrap
    
    if doperm_mvpa_reg_st

        fprintf ('\nSHOWING PERMUTATION %s RESULTS, %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

        figure

        o2 = canlab_results_fmridisplay([], 'compact', 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        t = perm_stats_obj;
        t = threshold(t, q_threshold_mvpa_reg_st, 'fdr', 'k', k_threshold_mvpa_reg_st); 
        r = region(t);

        o2 = addblobs(o2, r, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' permutation ' mask_string]);

        figtitle = sprintf('%s_%1.4f_permutation_FDR_montage_%s_%s', algorithm_mvpa_reg_st, q_threshold_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

        clear w, clear o2, clear figtitle

    end % if loop permutation
    

%% PERFORM SOURCE RECONSTRUCTION AKA STRUCTURE COEFFICIENTS ON SINGLE LEVEL MVPA MODELS IF REQUESTED
% --------------------------------------------------------------------------------------------------

    if dosourcerecon_mvpa_reg_st
        
        fprintf('\n\n');
        printhdr('SOURCE RECONSTRUCTION');
        fprintf('\n\n');
        
        fprintf('\n\n');
        printhdr('Calculating source reconstruction weights');
        fprintf('\n\n');

    % OBTAIN SOURCE RECONSTRUCTION WEIGHTS
    %-------------------------------------

        % NOTE: regress each voxel's activity onto yhat
        % Ref: Haufe et al, NeuroImage 2014

        regress_stats = cell(max(fold_labels),1);
        betas = zeros(max(fold_labels),size(fmri_dat.dat,1));

        for k = 1:max(fold_labels) % number of cv folds

            regress_data = fmri_dat;
            
            switch ml_method_mvpa_reg_st
    
                case 'predict'
                    regress_data.X = stats.yfit(fold_labels==k);
                    
                case 'oofmridataobj'
                    regress_data.X = cvGS.yfit{k};
                    
            end
            
            regress_data.dat = regress_data.dat(:,fold_labels==k);
            regress_stats{k} = regress(regress_data,'nodisplay','noverbose');
%                         tv = replace_empty(regress_stats(k).b);
            betas(k,:) = regress_stats{k}.b.dat(:,1);

        end % for loop over folds

        clear k;

        mean_beta = mean(betas);
        
        source_recon_data_obj = stats.weight_obj; % initialize as fmri_data object with correct volInfo
        source_recon_data_obj.dat = mean_beta';
        

    % VISUALIZE UNTHRESHOLDED SOURCE RECONSTRUCTION RESULTS
    % -----------------------------------------------------
    
        fprintf('\n\n');
        printhdr('Visualizing unthresholded source reconstruction maps');
        fprintf('\n\n');
        
        fprintf ('\nSHOWING UNTHRESHOLDED SOURCE RECONSTRUCTION %s RESULTS, %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

        figure

        o2 = canlab_results_fmridisplay([], 'compact', 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        w = region(source_recon_data_obj);

        o2 = addblobs(o2, w, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' unthresholded source reconstruction' mask_string]);

        figtitle = sprintf('%s_unthresholded_source_reconstruction_montage_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

        clear w, clear o2, clear figtitle
    
    end % if loop source reconstruction
    
    
    if dosourcerecon_perm_mvpa_reg_st
        
    % PERFORM PERMUTATION ON SOURCE RECONSTRUCTION WEIGHTS
    % ----------------------------------------------------
    
        fprintf('\n\n');
        printhdr('Permutation tests on source reconstruction maps');
        fprintf('\n\n');

        % GET PROBABILITY FOR OBSERVED SOURCE RECONSTRUCTION WEIGHTS FROM PERMUTATON WEIGHTS

        for beta = 1:size(null_beta,2)

            phat_srec(beta) = 1-sum(mean_beta(beta) > null_beta(:,beta)) / (1+size(null_beta,1));

            switch perm_sidedness
                case 'both'
                    phat_srec(beta) = (length(find(abs(null_beta(:,beta)) > abs(mean_beta(beta))))+1) / (perm_n_mvpa_reg_st+1); 
                case 'smaller'
                    phat_srec(beta) = (length(find(null_beta(:,beta) < mean_beta(beta)))+1) / (perm_n_mvpa_reg_st+1);
                case 'larger'
                    phat_srec(beta) = (length(find(null_beta(:,beta) > mean_beta(beta)))+1) / (perm_n_mvpa_reg_st+1);
            end

        end

        source_recon_perm_stats_obj = bs_stats.weight_obj;
        source_recon_perm_stats_obj.dat = mean_beta';
        source_recon_perm_stats_obj.p = phat_srec';
        
        % VISUALIZE THRESHOLDED SOURCE RECONSTRUCTION RESULTS AFTER PERMUTATION
        
        fprintf('\n\n');
        printhdr('Visualizing thresholded source reconstruction maps after permutation tests');
        fprintf('\n\n');

        fprintf ('\nSHOWING PERMUTATION SOURCE RECONSTRUCTION %s RESULTS, %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

        figure

        o2 = canlab_results_fmridisplay([], 'compact', 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        t = source_recon_perm_stats_obj;
        t = threshold(t, q_threshold_mvpa_reg_st, 'fdr', 'k', k_threshold_mvpa_reg_st); 
        r = region(t);

        o2 = addblobs(o2, r, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' permutation source reconstruction ' mask_string]);

        figtitle = sprintf('%s_%1.4f_permutation_source_reconstruction_FDR_montage_%s_%s', algorithm_mvpa_reg_st, q_threshold_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

        clear w, clear o2, clear figtitle
       
    end % if loop permutation testing on source reconstruction
    
    
%% SAVE STATS FOR SINGLE LEVEL MODELS
% -------------------------------------------------------------------------

if dosavemvparegstats
    
    fprintf('\n\n');
    printhdr('SAVING SINGLE LEVEL RESULTS');
    fprintf('\n\n');
    
    if exist('maskname_short', 'var')
        
        if ~isempty(cons2exclude_dat_st)
            savefilename = fullfile(resultsdir, ['single_trial_', myscaling_mvpa_reg_st, '_', maskname_short, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '_results.mat']);
            
        else
            savefilename = fullfile(resultsdir, ['single_trial_', myscaling_mvpa_reg_st, '_', maskname_short, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_', results_suffix, '_results.mat']);
        
        end
        
    else
        
        if ~isempty(cons2exclude_dat_st)
                        savefilename = fullfile(resultsdir, ['single_trial_', myscaling_mvpa_reg_st, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '_results.mat']);
            
        else
            savefilename = fullfile(resultsdir, ['single_trial_', myscaling_mvpa_reg_st, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_', results_suffix, '_results.mat']);
        
        end
    end
    
    
    switch ml_method_mvpa_reg_st
        
        case 'predict'
    
            save(savefilename, 'stats','slope_stats', '-v7.3');
            
        case 'oofmridataobj'
            
            save(savefilename, 'stats', '-v7.3');
            
    end
    
    
    if dobootstrap_mvpa_reg_st
        save(savefilename,'bs_stats', '-append');
    end

    if doperm_mvpa_reg_st
        save(savefilename,'perm_stats_obj', '-append');
    end

    if dosourcerecon_mvpa_reg_st
        save(savefilename,'source_recon_data_obj', '-append');
    end 

    if dosourcerecon_perm_mvpa_reg_st
        save(savefilename,'source_recon_perm_stats_obj', '-append');
    end 
        
end
        

%% FIT MULTILEVEL MVPA MODEL
% -------------------------------------------------------------------------

if domultilevel_mvpa_reg_st

    % NEW WIP CODE
    %-------------
    
    fmri_dat.metadata_table.(subj_identifier_dat_st) = subject_id;
    
    % DEFINE ALGORITHM
    
    clear alg

    alg = mlpcrRegressor('randInt',true,'cpca',false,'batch_id_funhan',@(X)(X.metadata));

%             alg.fit(fmri_dat.dat', fmri_dat.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
    % NOTE: fit is not strictly necessary at this stage, but a good test
    alg_params = alg.get_params; % get to know the hyperparams for this algorithm, which we want to optimize

    % NOTE: hyperparms for mlpcrRegressor

%             bt_dim = Inf;
%             wi_dim = Inf;
%             cpca = false; % balances number of trials per subject, but we
%             keep it to false now to be consistent with single level
%             scenario - it could be implemented there too
%             randInt = false;
%             randSlope = false; % probably not really needed, increases computation time

    % DEFINE FEATURE EXTRACTOR

    featConstructor_han = @(X)(X.metadata_table.(subj_identifier_dat_st));
    extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han); % initiate extractVxl as an empty fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
    extractVxl.fit(fmri_dat); % transformer takes fmri_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)
    % NOTE: fit is not strictly necessary at this stage, but a good test

    % DEFINE PIPELINE

    fmri_pipeline = pipeline({{'featExt',extractVxl},{'alg',alg}}); % define fmri_pcr as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
    fmri_pipeline.fit(fmri_dat,fmri_dat.Y);
    % NOTE: fit is not strictly necessary at this stage, but a good test

    % INNER CROSS-VALIDATION FUNCTION

    switch holdout_set_method_mvpa_reg_st

        case 'group'
            innercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier_dat_st), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st));

        case 'onesample'
            innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

    end
    % NOTE: we use metadata_table here, since the input to bo.fit is the
    % fmri_data_st object fmri_dat, which has a metadata_table field

    cv_folds_mvpa_reg_st = innercv(fmri_dat,fmri_dat.Y); 
    % NOTE: get cross-validation folds, not strictly necessary, but a good test

    % DEFINE BAYESIAN OPTIMIZATION

    scaling = 0.8;
    wi_dims = optimizableVariable('alg__wi_dim',[1,floor(scaling*(rank(fmri_dat.dat)-n_subj)*(nfolds_mvpa_reg_st-1)/nfolds_mvpa_reg_st)],'Type','integer');
    bt_dims = optimizableVariable('alg__bt_dim',[0,floor((nfolds_mvpa_reg_st-1)/nfolds_mvpa_reg_st*n_subj)-1],'Type','integer');

    % NOTE: Type and number of hyperparams to optimize depends on algorithm (check alg.get_params above), as well as other settings

    bayesOptParams = {[wi_dims,bt_dims], 'AcquisitionFunctionName','expected-improvement-plus',...
        'MaxObjectiveEvaluations',45, 'UseParallel', false, 'verbose',1};

    bo = bayesOptCV(fmri_pipeline,innercv,@get_mse,bayesOptParams); % consider l2norm images scaling option above to work with @get_mse - look at loss landscape in plots; alternative is to get rid of optimization and use max number of components
    bo.fit(fmri_dat,fmri_dat.Y);
    bo_numcomponents = [bo.estimator.estimator.wi_dim,bo.estimator.estimator.bt_dim];

    % OUTER CROSS-VALIDATION FUNCTION

    switch holdout_set_method_mvpa_reg_st

        case 'group'
            outercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier_dat_st), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st));

        case 'onesample'
            outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.(subj_identifier_dat_st)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

    end
    % NOTE: we use metadata_table here, since the input to cvGS.do is the
    % fmri_data_st object fmri_dat, which has a metadata_table field

    % ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE

    cvGS = crossValScore(bo, outercv, @get_mse, 'n_parallel', nfolds_mvpa_reg_st, 'verbose', true);
    % NOTE: Bogdan advises not parallizing too much for the purpose of Bayesian
    % model optimization, since each step learns from the previous one, so we
    % only parallelize the outer cv loop with 1 core per outer cv fold,
    % resulting in very acceptable runtimes (on LaBGAS server with 128 GB RAM)

    cvGS.do(fmri_dat, fmri_dat.Y);
    cvGS.do_null(); % fits null model - intercept only
    fold_labels = cvGS.fold_lbls;

    % CREATE FMRI_DATA OBJECTS WITH THE BETAS FOR VISUALIZATION PURPOSES

    weight_obj = bo.estimator.transformers{1}.brainModel; % empty .dat at this stage
    weight_obj.dat = bo.estimator.estimator.B(:); % fills mdl.dat with betas of the sum

    wi_weight_obj = bo.estimator.transformers{1}.brainModel; % empty .dat at this stage
    wi_weight_obj.dat = bo.estimator.estimator.Bw(:); % fills mdl.dat with betas

    bt_weight_obj = bo.estimator.transformers{1}.brainModel; % empty .dat at this stage
    bt_weight_obj.dat = bo.estimator.estimator.Bb(:); % fills mdl.dat with betas
    
    % PLOT OBSERVED VERSUS PREDICTED
    
    fprintf('\n\n');
    printhdr('Plotting observed versus predicted');
    fprintf('\n\n');

%     switch ml_method_mvpa_reg_st
% 
%         case 'predict'
% 
%             fprintf('\n%s r = %0.3f\n\n', algorithm_mvpa_reg_st, corr(stats.yfit, fmri_dat.Y));
% 
%             figure
%             
%             nr_subjs = max(unique(subject_id));
%             
%             colors = distinguishable_colors(nr_subjs);
%             colors_cell = mat2cell(colors,ones(1,nr_subjs));
%             [han,~,~,slope_stats] = line_plot_multisubject(fmri_dat.Y, stats.yfit, 'subjid', subject_id, 'group_avg_ref_line', 'MarkerTypes','d', 'colors', colors_cell);
%             xlabel({['Observed ' behav_outcome_dat_st],'(across conditions)'},'FontSize',14); 
%             ylabel({['Estimated ' behav_outcome_dat_st],'(cross validated)'},'FontSize',14);
% 
%             for i = 1:size(han.line_handles,2)
%                 han.line_handles(1,i).LineWidth = 2.5;
%                 han.point_handles(1,i).MarkerSize = 5;
%             end
%             
%             han.grpline_handle.LineWidth = 7.5;
%             han.grpline_handle.LineStyle = '-';
%             han.grpline_handle.Color = 'r';
%             
%             ax = gca;
%             ax.FontWeight = 'bold';
% 
%             set(gcf,'WindowState','Maximized','Color','w');
%             drawnow, snapnow;
% 
% 
%         case 'oofmridataobj'

            for y = 1:size(cvGS.Y,2)
                r(y,1) = corr(cvGS.yfit{y},cvGS.Y{y});
            end

            cv_r = mean(r); % average r(yfit,y) over folds = cv correlation

            cv_mse = mean(cvGS.scores); % average MSE over folds = cv model performance

            fprintf('\n%s cross-validated r = %0.3f\n\n', algorithm_mvpa_reg_st, cv_r);
            fprintf('\n%s cross-validated mse = %0.3f\n\n', algorithm_mvpa_reg_st, cv_mse);

            f1 = cvGS.plot; % plots predicted versus observed

            set(gcf,'WindowState','Maximized');
            drawnow, snapnow;

%     end


    % PLOT MONTAGE OF UNTHRESHOLDED WEIGHTS

    fprintf('\n\n');
    printhdr('Plotting unthresholded weight maps');
    fprintf('\n\n');

    whmontage = 5;
    
    % total

    fprintf ('\nSHOWING UNTHRESHOLDED MULTILEVEL %s RESULTS, %s, TOTAL EFFECT, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

    figure;

    o2 = canlab_results_fmridisplay([], 'compact');

%         switch ml_method_mvpa_reg_st
% 
%             case 'predict'
%                 w = region(stats.weight_obj);

%             case 'oofmridataobj'
                w = region(weight_obj);

%         end

    o2 = addblobs(o2, w, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    o2 = title_montage(o2, whmontage, ['multilevel ' algorithm_mvpa_reg_st ' total effect, unthresholded ' mask_string]);

    figtitle = sprintf('%s_unthresholded_montage_multilevel_total_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;

    clear w, clear o2, clear figtitle
    
    % within
    
    fprintf ('\nSHOWING UNTHRESHOLDED MULTILEVEL %s RESULTS, %s, WITHIN-SUBJECT EFFECT, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

    figure;

    o2 = canlab_results_fmridisplay([], 'compact');

%         switch ml_method_mvpa_reg_st
% 
%             case 'predict'
%                 w = region(stats.weight_obj);

%             case 'oofmridataobj'
                w = region(wi_weight_obj);

%         end

    o2 = addblobs(o2, w, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    o2 = title_montage(o2, whmontage, ['multilevel ' algorithm_mvpa_reg_st ' within-subject effect, unthresholded ' mask_string]);

    figtitle = sprintf('%s_unthresholded_montage_multilevel_within_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;

    clear w, clear o2, clear figtitle
    
    % between
    
    fprintf ('\nSHOWING UNTHRESHOLDED MULTILEVEL %s RESULTS, %s, BETWEEN-SUBJECT EFFECT, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

    figure;

    o2 = canlab_results_fmridisplay([], 'compact');

%         switch ml_method_mvpa_reg_st
% 
%             case 'predict'
%                 w = region(stats.weight_obj);

%             case 'oofmridataobj'
                w = region(bt_weight_obj);

%         end

    o2 = addblobs(o2, w, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    o2 = title_montage(o2, whmontage, ['multilevel ' algorithm_mvpa_reg_st ' between-subject effect, unthresholded ' mask_string]);

    figtitle = sprintf('%s_unthresholded_montage_multilevel_between_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;

    clear w, clear o2, clear figtitle
            
            
    % OLD PREDICT CODE, COULD BE REINTEGRATED
    %---------------------------------------- 
    
%     % DETERMINE MAXIMUM NUMBER OF COMPONENTS
%     %---------------------------------------
% 
%     max_comp = floor(size(fmri_dat.dat,2).*0.75 - n_subj);
% 
%     % NOTE: we specify the maximum number of components as < the number of columns in
%     % fmri_dat.dat (n_subjects*n_conditions(in every fold)) to avoid overfitting in multilevel models, 
%     % where we need to leave df for the random intercepts (upper bound 1 df per random intercept hence subject)
% 
%     % FIT SINGLE LEVEL MODEL
%     %-----------------------
% 
%     [sl_cverr, sl_stats,sl_optout] = fmri_dat.predict('algorithm_name','cv_pcr',...
%         'nfolds',fold_labels, 'numcomponents',max_comp);
%     fprintf('PCR r = %0.3f\n', corr(sl_stats.yfit,fmri_dat.Y));
% 
%     figure
% 
%     line_plot_multisubject(fmri_dat.Y, sl_stats.yfit, 'subjid', subject_id);
%     xlabel({['Observed ' behav_outcome_dat_st],'(average over conditions)'}); ylabel({['PCR Estimated ' behav_outcome_dat_st],'(cross validated)'})
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     figure
% 
%     sl_stats.weight_obj.montage;
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
% 
%     % FIT MULTILEVEL MODEL W/FIXED PARAMETER ESTIMATION
%     %--------------------------------------------------
% 
%     % split maximum amount of components in between and within
% 
%     n_bt_comp = floor(0.75*n_subj);
%     n_wi_comp = max_comp - n_bt_comp;
% 
%     % NOTE: max between = n_subj IN EVERY FOLD (hence n_subj - 20% in 5-fold CV), 
%     % and you want to put more money on within since this typically explains
%     % more variance
% 
%     % overall model prediction
% 
%     [ml_cverr, ml_stats, ml_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr',...
%         'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id);
%     fprintf('multilevel PCR r = %0.3f\n',corr(ml_stats.yfit, fmri_dat.Y));
% 
%     % NOTE: algorithm option 'cv_mlpcr' requires subject identifier,
%     % which makes sense since this is a multilevel/mixed model
%     % note that fold labels are the same, since they respect subject membership
% 
%     figure
% 
%     line_plot_multisubject(fmri_dat.Y, ml_stats.yfit, 'subjid', subject_id);
%     xlabel({['Observed ' behav_outcome_dat_st],'(average over conditions)'}); ylabel({['MLPCR Estimated ' behav_outcome_dat_st],'(cross validated)'})
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     figure
% 
%     ml_stats.weight_obj.montage;
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     figure
% 
%     subplot(1,2,1)
%     line_plot_multisubject(sl_stats.yfit, ml_stats.yfit, 'subjid', subject_id);
%     xlabel({'PCR model prediction'}); ylabel('Multilevel PCR model prediction');
%     axis square
%     subplot(1,2,2);
%     plot(sl_optout{1}(:),ml_optout{1}(:),'.');
%     lsline;
%     xlabel('PCR model weights'); ylabel('Multilevel PCR model weights');
%     axis square
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     % NOTE: contrary to @bogpetre's walkthrough, the
%     % pcr and the multilevel pcr models are not exactly equivalent anymore
%     % since I have been specifying the number of components
% 
%     % get the variance explained by the between and within component
% 
%     % NOTE: These functions call the same thing under the hood, 
%     % but simply perform cross validation using ONLY between or within
%     % subject models.
% 
%     [ml_bt_cverr, ml_bt_stats] = fmri_dat.predict('algorithm_name','cv_mlpcr_bt',...
%         'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
%     pred_bt = ml_bt_stats.yfit;
% 
%     [ml_wi_cverr, ml_wi_stats] = fmri_dat.predict('algorithm_name','cv_mlpcr_wi',...
%         'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
%     pred_wi = ml_wi_stats.yfit;
% 
%     % NOTE: algorithm options are created by bogpetre
% 
%     fprintf('Between subject PCR components r = %0.3f\n', corr(ml_bt_stats.yfit, fmri_dat.Y));
%     fprintf('Within subject PCR components r = %0.3f\n', corr(ml_wi_stats.yfit, fmri_dat.Y));
% 
%     figure
% 
%     subplot(1,2,1)
%     line_plot_multisubject(fmri_dat.Y, pred_bt, 'subjid', subject_id);
%     xlabel({['Observed ' behav_outcome_dat_st]}); ylabel('Between subject components'' prediction');
%     axis square
%     subplot(1,2,2)
%     line_plot_multisubject(fmri_dat.Y, pred_wi, 'subjid', subject_id);
%     xlabel({['Observed ' behav_outcome_dat_st]}); ylabel('Within subject components'' prediction');
%     axis square
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
% 
%     % FIT MULTiLEVEL MODEL W/ MIXED PARAMETER ESTIMATION
%     %---------------------------------------------------
% 
%     % NOTE: this is not part of the walkthrough yet, but @bogpetre pushed
%     % mlpcr3 function to CanlabCore
% 
%     % NOTE bogpetre:
%     % main function is mlpcr3, so help mlpcr3 for usage options.
%     % basically it's the same as cv_mlpcr, except there's a randInt, randSlope and fitlmeOpts now
%     % fitlmeOpts get passed on to fitlme. It picks some sensible defaults
%     % randSlope makes things much slower. randInt is roughly the sae order of magnitude as running the fixed effects version
%     % the function defaults to fixed effects by default, so it's a drop in replacement
%     % for mlpcr2.m (aka cv_mlpcr)
% 
%     % overall model prediction including random intercept only
% 
%     [ml3_cverr, ml3_stats, ml3_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr3',...
%         'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1);
% 
%     fprintf('multilevel PCR r = %0.3f\n',corr(ml3_stats.yfit, fmri_dat.Y));
% 
%     % NOTE: compare with code in previous section and note change of algorithm_name option to cv_mlpcr3 and addition of randInt
%     % option - see help mlpcr3 for more details
% 
%     figure
% 
%     line_plot_multisubject(fmri_dat.Y, ml3_stats.yfit, 'subjid', subject_id);
%     xlabel({['Observed ' behav_outcome_dat_st],'(average over conditions)'}); ylabel({['MLPCR3 Estimated ' behav_outcome_dat_st],'(cross validated)'})
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     figure
% 
%     ml3_stats.weight_obj.montage;
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     figure
% 
%     subplot(1,2,1)
%     line_plot_multisubject(ml_stats.yfit, ml3_stats.yfit, 'subjid', subject_id);
%     xlabel({'Multilevel PCR model prediction'}); ylabel('Multilevel PCR model random int model prediction');
%     axis square
%     subplot(1,2,2);
%     plot(ml_optout{1}(:),ml3_optout{1}(:),'.');
%     lsline;
%     xlabel('Multilevel PCR model weights'); ylabel('Multilevel PCR model random int model weights');
%     axis square
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     % NOTE: models are very similar but not exactly equivalent
% 
%     % overall model prediction including random intercept and random slope
% 
%     [ml3rs_cverr, ml3rs_stats, ml3rs_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr3',...
%         'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1, 'randSlope', 1);
% 
%     fprintf('multilevel PCR r = %0.3f\n',corr(ml3rs_stats.yfit, fmri_dat.Y));
% 
%     figure
% 
%     line_plot_multisubject(fmri_dat.Y, ml3rs_stats.yfit, 'subjid', subject_id);
%     xlabel({['Observed ' behav_outcome_dat_st],'(average over conditions)'}); ylabel({['MLPCR3 Estimated ' behav_outcome_dat_st],'(cross validated)'})
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;
% 
%     figure
% 
%     ml3rs_stats.weight_obj.montage;
% 
%     set(gcf, 'WindowState','maximized');
%     drawnow, snapnow;


    %% SAVE STATS FOR MULTILEVEL MODELS
    %--------------------------------------------------------------------------
    
    if dosavemvparegstats
        
        fprintf('\n\n');
        printhdr('SAVING MULTILEVEL RESULTS');
        fprintf('\n\n');
    
        if exist('maskname_short', 'var')
        
            if ~isempty(cons2exclude_dat_st)
                savefilename = fullfile(resultsdir, ['single_trial_multilevel_', myscaling_mvpa_reg_st, '_', maskname_short, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '_results.mat']);

            else
                savefilename = fullfile(resultsdir, ['single_trial_multilevel_', myscaling_mvpa_reg_st, '_', maskname_short, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_', results_suffix, '_results.mat']);

            end
        
        else
        
            if ~isempty(cons2exclude_dat_st)
                            savefilename = fullfile(resultsdir, ['single_trial_multilevel_', myscaling_mvpa_reg_st, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '_results.mat']);

            else
                savefilename = fullfile(resultsdir, ['single_trial_multilevel_', myscaling_mvpa_reg_st, '_', algorithm_mvpa_reg_st, '_', behav_outcome_dat_st, '_', results_suffix, '_results.mat']);

            end
        end
        
        save(savefilename, 'weight_obj','wi_weight_obj','bt_weight_obj', '-v7.3');
        
    end
    
end % if loop multilevel

