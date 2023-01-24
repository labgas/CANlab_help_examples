%% prep_3c_run_SVMs_on_contrasts_masked.m
%
%
% USAGE
%
% This script 
% 1) runs second⁻level (i.e. across subjects) support vector machines
% for each within-subject CONTRAST registered in DAT.contrasts
% 2) plots montages of the uncorrected results
% 3) saves the results using standard naming and location
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('prep_3c_run_SVMs_on_contrasts_masked','outputDir',htmlsavedir)
%
% To get results reports after bootstrapping, publish c2_SVM_contrasts_masked
%
%
% OPTIONS
%
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options (for example
% comparing different machine learning methods), you can make a copy of this script with
% a letter index (e.g. _s6a_) and change the default option here
%
% - maskname_svm: 
%       - default use of sparse gray matter mask
%       - maskdir now defined in a_set_up_paths_always_run_first script
%       - if you do not want to mask, change to []
%       - if you want to use a custom mask, put it in maskdir and change name here.
% - myscaling_svm: options are 'raw','subjectnorm','imagenorm','zscoreimages','zscorevoxels'
%                       subjectnorm: normalize_each_subject_by_l2norm; normalizes images for each subject by L2 norm of Condition 1 image; can help with numerical scaling and inter-subject scaling diffs
%                       imagenorm: normalize_images_by_l2norm; normalizes each image separately, not each subject/pair
%                       zscoreimages: Z-score each input image, removing image mean and forcing std to 1. Removes overall effects of image intensity and scale. Can be useful across studies but also removes information. Use judiciously. lukasvo76: corresponds to 'scaled' in myscaling_glm option in prep_3a
%                       zscorevoxels: Z-score each voxel across images
% - dosavesvmstats: default true; Save statistics and weight map objects for SVM contrasts
% - dobootstrap_svm: default false; Takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n_svm is set to the default 10k samples
%       boot_n_svm: default 5000; number of bootstrap samples       Reduce number for quick results
% - parallelstr_svm: default 'parallel'; parallel processing for bootstrapping.   'parallel' or 'noparallel'
% - holdout_set_method_svm: 'group', or 'onesample'
%                       Group: use DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.contrasts{c}.group for balancing holdout set over groups;
%                       Onesample: use subject identifier only for stratifying by subject (i.e. leave whole subject out)
%
%       'group' option 
%           Assuming that groups are concatenated in contrast image lists, and
%           regressor values of 1 or -1 will specify the group identity for each image. 
%           Requires DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.contrasts{c}.group 
%           field specifying group membership for each image, calling
%           plugin_get_holdout_sets_balanced_groups to achieve balanced
%           holdout sets over groups
%
%       'onesample' option:
%           Uses subject identifier only, calling plugin_get_holdout_sets
%           to stratify by subject (i.e. leave both images of same subject
%           out in holdout set)
%
% - holdout_set_type_svm: 'kfold', or 'leave_one_subject_out': default 'kfold'; choose between kfold or leave one subject out cross-validation - the latter is not recommended
%       nfolds_svm: default 5, number of folds for kfold CV
% - dosearchlight_svm: perform searchlight SVM analysis 
%       searchlight_radius_svm: radius for searchlight sphere
%       IMPORTANT NOTE: searchlight option is currently not working because of
%                       weird memory issues in fmri_data.searchlight - DO NOT USE UNTIL FIXED
% - ml_method_svm: 'oofmridataobj', or 'predict'
%       'oofmridataobj' option:
%           use @bogpetre's object-oriented method
%           https://github.com/canlab/ooFmriDataObjML
%       'predict' option:
%           use CANlab's predict function
%           https://github.com/canlab/CanlabCore/blob/master/CanlabCore/%40fmri_data/predict.m
%
% MANDATORY OPTIONS TO BE SPECIFIED IN THIS SCRIPT
%
% - results_suffix: name to add to results file to specify in case of multiple version of model, e.g. 'oofmridataobj'
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% prep_3c_run_SVMs_on_contrasts_masked.m         v3.2
% last modified: 2022/08/25


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% SET MANDATORY OPTIONS

results_suffix = ''; % adds a suffix of your choice to .mat file with results that will be saved
% NOTE: do NOT delete this option, leave empty if not needed
% NOTE: do NOT use to add a suffix specifying the scaling or masking option, this will be added automatically

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;
% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT!
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% GET DEFAULT OPTIONS IF NOT SET IN A2_SET_DEFAULT_OPTIONS

options_needed = {'dosavesvmstats', 'dobootstrap_svm', 'boot_n_svm', 'holdout_set_method_svm', 'holdout_set_type_svm', 'nfolds_svm', 'ml_method_svm'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {true false 5000 'onesample' 'kfold' 5 'predict'};

plugin_get_options_for_analysis_script;

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% ml_method_svm: 'oofmridataobj'/'predict';
% holdout_set_method_svm = 'onesample'/'group';
% holdout_set_type_svm = 'kfold'/'leave_one_subject_out';
%    nfolds_svm = x;
% maskname_svm = []/which(maskname);
% myscaling_svm = 'raw'/'subjectnorm'/'imagenorm'/'zscoreimages'/'zscorevoxels'
% dosavesvmstats = true/false;
% dobootstrap_svm = true/false;
%    boot_n_svm = yyyy;
% parallelstr_svm = 'parallel'/'noparallel';
% dosearchlight_svm = true/false;
%    searchlight_radius_svm = z;


%% LOAD NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

if ~exist('DATA_OBJ','var')
    
    load(fullfile(resultsdir,'data_objects.mat'));
    
end


%% CHECK DEPENDENCIES
% -------------------------------------------------------------------------

switch ml_method_svm
    
    case 'predict'
        spath = which('use_spider.m');
        if isempty(spath)
            error('Spider toolbox not found on Matlab path, clone CanlabCore repo and add to Matlab path to prevent prediction from breaking')
        end
        
    case 'oofmridataobj'
        opath = which('bayesOptCV.m');
        if isempty(opath)
            error('ooFmriDataObjML repo not found on Matlab path, clone repo and add to Matlab path to prevent prediction from breaking')
        end
        
    otherwise
        error('\ninvalid option "%s" defined in ml_method_svm variable, choose between "oofmridataobj" and "predict"\n',ml_method_svm);
        
end
            

%% GET MASK
% -------------------------------------------------------------------------

if exist('maskname_svm', 'var') && ~isempty(maskname_svm) 
    svmmask = fmri_mask_image(maskname_svm, 'noverbose');
    [~,maskname_short] = fileparts(maskname_svm);
    mask_string = sprintf('masked with %s', maskname_short);

end


%% RUN SUPPORT VECTOR MACHINES FOR EACH CONTRAST
% -------------------------------------------------------------------------

kc = size(DAT.contrasts, 1);

svm_stats_results = cell(1, kc);

if dosearchlight_svm
    searchlight_svm_stats = cell(1, kc);
    searchlight_svm_objs = cell(1, kc);
    
    for cond = 1:size(DATA_OBJ,2) % need to convert our fmri_data_st objects to fmri_data to prevent searchlight function from breaking on some removed voxel issue
        DATA_OBJ{1,cond} = fmri_data(DATA_OBJ{1,cond});
    
    end
    
end

if dobootstrap_svm
   bootstrap_svm_stats = cell(1, kc); 
end

for c = 1:kc
    
    analysisname = DAT.contrastnames{c};
    
    fprintf('\n\n');
    printhdr(['CONTRAST #', num2str(c), ': ', upper(analysisname)]);
    fprintf('\n\n');
    
    mycontrast = DAT.contrasts(c, :);
    wh = find(mycontrast); % wh is which conditions have non-zero contrast weights
    
    % CREATE COMBINED DATA OBJECT WITH INPUT IMAGES FOR BOTH CONDITIONS
    % ---------------------------------------------------------------------
    
    [cat_obj, condition_codes] = cat(DATA_OBJ{wh}); 
    cat_obj = enforce_variable_types(cat_obj); % @lukasvo76 added as the fmri_data.cat function includes replace_empty on the objects, which causes problems later on with the stats_object output of predict
    
    % APPLY MASK IF SPECIFIED IN OPTIONS
    %----------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('Masking and scaling images if requested in options');
    fprintf('\n\n');
    
    if exist('svmmask', 'var')
        fprintf('\nMasking data with %s\n\n',maskname_short);
        cat_obj = apply_mask(cat_obj, svmmask);
        cat_obj.mask_descrip = maskname_svm;
        
    else
        fprintf('\nNo mask found; using full existing image data\n\n');
        
    end
    
    % NORMALIZE IF SPECIFIED IN OPTIONS
    % ---------------------------------------------------------------------
    
    % NORMALIZE BY L2NORM
    
    % possibly normalize_each_subject_by_l2norm; can help with numerical scaling and inter-subject scaling diffs
    % Sometimes condition differences are very small relative to baseline
    % values and SVM is numerically unstable. If so, re-normalizing each
    % subject can help.
    
    switch myscaling_svm
        
        case 'raw'
    
            scaling_string = 'no_scaling';
            
        case 'subjectnorm'
    
            cat_obj = normalize_each_subject_by_l2norm(cat_obj, condition_codes);
            scaling_string = 'scaling_l2norm_subjects';
            fprintf('\nNormalizing condition images for each subject by L2 norm of Condition 1 image before SVM\n\n');
            
        case 'imagenorm'
         
            cat_obj = normalize_images_by_l2norm(cat_obj);
            scaling_string = 'scaling_l2norm_conditions';
            fprintf('\nNormalizing each condition image for each subject by L2 norm before SVM\n\n');
    
    % Z-SCORE
    
    % Z-score each input image, removing image mean and forcing std to 1.
    % Removes overall effects of image intensity and scale. Can be useful
    % across studies but also removes information. Use judiciously.
    
        case 'zscoreimages'
            fprintf('\nZ-scoring each condition image for each subject before SVM\n\n');
            scaling_string = 'scaling_z_score_conditions';
            cat_obj = rescale(cat_obj, 'zscoreimages');
            
        case 'zscorevoxels'
            fprintf('\nZ-scoring each condition image for each subject before SVM\n\n');
            scaling_string = 'scaling_z_score_conditions';
            cat_obj = rescale(cat_obj, 'zscoreimages');
            
        otherwise
            error('incorrect scaling option %s specified in myscaling_svm option in a2_set_default_options.\nChoose between "raw", "subjectnorm", "imagenorm", "zscoreimages", or "zscorevoxels"\n\n', myscaling_svm);
        
    end
    
    % DEFINE HOLDOUT SETS
    % -------------------
    % NOTE: use plugin scripts according to holdout_set_method_svm

        switch holdout_set_method_svm

            case 'group' 
                plugin_get_holdout_sets_balanced_groups; % @lukasvo76: built in this option to manually define your holdout sets balancing for group variable, wrote a new plugin based on @bogpetre's walkthrough code for single-trials and between-within MVPA

            case 'onesample'
                plugin_get_holdout_sets; % @lukasvo76: this is the original CANlab code which works fine if you do not have to balance your holdout sets for a group variable

            otherwise
                error('\ninvalid option "%s" defined in holdout_set_method_svm variable, choose between "group" and "onesample"\n\n',holdout_set_method_svm);

        end
    
    % FORMAT AND ATTACH OUTCOME
    % -------------------------
    % NOTES:
    % a. 1, -1 for the two conditions in the contrast
    % b. assume that subjects are in same position in each input file!
            
    cat_obj.Y = outcome_value;
    cat_obj.metadata_table.subject_id = [[1:(size(cat_obj.Y,1)/2)]';[1:(size(cat_obj.Y,1)/2)]'];
        if strcmp(holdout_set_method_svm,'group')
            cat_obj.metadata_table.group_id = [group';group'];      
        end
    
    % SANITY CHECK ON OUTCOME
    
    if all(cat_obj.Y > 0) || all(cat_obj.Y < 0)
        % Only positive or negative weights - nothing to compare
        fprintf('\n');
        warning(' Only positive or negative weights - nothing to compare');
        fprintf('\n');
        
        continue    
        
    end
    
    % RUN PREDICTIVE SVM MODEL
    % --------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('Running cross-validated SVM');
    fprintf('\n\n');
    
    switch ml_method_svm
    
        case 'predict'

            % RUN MODEL USING CANLAB'S PREDICT FUNCTION

            [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, ...
                'error_type', 'mcr', parallelstr_svm, 'verbose' ,0);

        case 'oofmridataobj'

            % DEFINE ALGORITHM

            alg = linearSvmClf();
            alg.fit(cat_obj.dat', cat_obj.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
            % NOTE: fit is not strictly necessary at this stage, but a good test
            alg_params = alg.get_params; % get to know the hyperparams for this algorithm, which we want to optimize
            
            % DEFINE FEATURE EXTRACTOR
            
            featConstructor_han = @(X)([]);
            extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han); % initiate extractVxl as an empty fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
            extractVxl.fit(cat_obj); % transformer takes fmri_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)
            % NOTE: fit is not strictly necessary at this stage, but a good test
            
            % DEFINE PIPELINE

            fmri_pipeline = pipeline({{'featExt',extractVxl},{'alg',alg}}); % define fmri_pcr as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
            fmri_pipeline.fit(cat_obj,cat_obj.Y);
            % NOTE: fit is not strictly necessary at this stage, but a good test
            
            % INNER CROSS-VALIDATION FUNCTION

            switch holdout_set_method_svm

                case 'group'
                    innercv = @(X,Y) cvpartition2(X.metadata_table.group_id, 'GroupKFold', nfolds_svm, 'Group', X.metadata_table.subject_id);

                case 'onesample'
                    innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_svm, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

            end
            % NOTE: we use metadata_table here, since the input to bo.fit is the
            % fmri_data_st object cat_obj, which has a metadata_table field
            
            cv_folds_svm = innercv(cat_obj,cat_obj.Y); 
            % NOTE: get cross-validation folds, not strictly necessary, but a good test
            
            % DEFINE BAYESIAN OPTIMIZATION

            dims = optimizableVariable('alg__C',[0.1,100]);

            % NOTE: Type and number of hyperparams to optimize depends on algorithm (check alg.get_params above), as well as other settings

            bayesOptParams = {dims, 'AcquisitionFunctionName','expected-improvement-plus',...
                'MaxObjectiveEvaluations',30, 'UseParallel', false, 'verbose',1, 'PlotFcn', {}};

            bo = bayesOptCV(fmri_pipeline,innercv,@get_hinge_loss,bayesOptParams);
            bo.fit(cat_obj,cat_obj.Y);
            bo_C = bo.estimator.estimator.C;
            
            % OUTER CROSS-VALIDATION FUNCTION

            switch holdout_set_method_svm

                case 'group'
                    outercv = @(X,Y) cvpartition2(X.metadata_table.group_id, 'GroupKFold', nfolds_svm, 'Group', X.metadata_table.subject_id);

                case 'onesample'
                    outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_svm, 'Group', X.metadata_table.subject_id);

            end
            % NOTE: we use metadata_table here, since the input to cvGS.do is the
            % fmri_data_st object fmri_dat, which has a metadata_table field
            
            % ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE

            cvGS = crossValScore(bo, outercv, @get_f1_macro, 'n_parallel', nfolds_svm, 'verbose', true);
            % NOTE: Bogdan advises not parallizing too much for the purpose of Bayesian
            % model optimization, since each step learns from the previous one, so we
            % only parallelize the outer cv loop with 1 core per outer cv fold,
            % resulting in very acceptable runtimes (on LaBGAS server with 128 GB RAM)

            cvGS.do(cat_obj, cat_obj.Y);
            cvGS.do_null(); % fits null model - intercept only
            fold_labels = cvGS.fold_lbls;
            
            % CREATE AN FMRI_DATA OBJECT WITH THE BETAS FOR VISUALIZATION PURPOSES
            
            weight_obj = bo.estimator.transformers{1}.brainModel; % empty .dat at this stage
            weight_obj.dat = bo.estimator.estimator.B(:); % fills mdl.dat with betas
            
            % KEEP IMPORTANT VARIABLES/OBJECTS IN STRUCTURE FOR SAVING
            
            stats = struct('bo',bo,'cvGS',cvGS,'weight_obj',weight_obj);

        otherwise

            error('\ninvalid option "%s" defined in ml_method_mvpa_svm variable, choose between "oofmridataobj" and "predict"\n\n',ml_method_svm);

    end % switch machine learning method
    
    
    % RUN SEARCHLIGHT SVM MODEL IF REQUESTED IN OPTIONS
    %----------------------------------------------------------------------
    
    if dosearchlight_svm
        
        fprintf('\n\n');
        printhdr('Running cross-validated searchlight SVM');
        fprintf('\n\n');
        
        delete(gcp('nocreate'));
        clust = parcluster('local'); % determine local number of cores, and initiate parallel pool with 80% of them
        nw = clust.NumWorkers;
        parpool(round(0.8*nw));
        
        [searchlight_obj, searchlight_stats, ~] = searchlightLukas(cat_obj, 'algorithm_name', 'cv_svm', ...
            'r', searchlight_radius_svm, 'holdout_set', holdout_set, 'do_online', 'no_weights');
    
    end
    
    % PLOT MONTAGE OF UNTHRESHOLDED SVM RESULTS
    % ---------------------------------------------------------------------
    
    whmontage = 5;
    
    fprintf('\n\n');
    printhdr('Visualizing cross-validated SVM results');
    fprintf('\n\n');

    fprintf ('\nMONTAGE UNTHRESHOLDED SVM RESULTS, CONTRAST: %s, %s, SCALING: %s\n\n', analysisname, mask_string, scaling_string);
    
        switch ml_method_svm

            case 'predict'
                r = region(stats.weight_obj);
                
            case 'oofmridataobj'
                r = region(weight_obj);
        end
    
    o2 = montage(r, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
    o2 = title_montage(o2, whmontage, [analysisname ' unthresholded ' mask_string ' ' scaling_string]);

    figtitle = sprintf('%s_unthresholded_montage_%s_%s', analysisname, mask_string, scaling_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
        if save_figures_svm
            plugin_save_figure;
        end
    clear o2, clear figtitle
    
    % PLOT MONTAGE OF UNTHRESHOLDED SEARCHLIGHT SVM RESULTS IF REQUESTED
    % ---------------------------------------------------------------------
    
    if dosearchlight_svm

        fprintf ('\nMONTAGE UNTHRESHOLDED SEARCHLIGHT SVM RESULTS, CONTRAST: %s, %s, SCALING: %s\n\n', analysisname, mask_string, scaling_string);

        r = region(stats.weight_obj);

        o2 = montage(r, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [analysisname ' unthresholded searchlight ' mask_string ' ' scaling_string]);

        figtitle = sprintf('%s_unthresholded_searchlight_montage_%s_%s', analysisname, mask_string, scaling_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;
            if save_figures_svm
                plugin_save_figure;
            end
        clear o2, clear figtitle
        
    end % if loop bootstrap
    
    
    % BOOTSTRAP IF REQUESTED
    % --------------------------------------------------------------------
    
    if dobootstrap_svm
        
        fprintf('\n\n');
        printhdr('Bootstrapping SVM weights');
        fprintf('\n\n');
        
        delete(gcp('nocreate'));
        c = parcluster('local'); % determine local number of cores, and initiate parallel pool with 80% of them
        nw = c.NumWorkers;
        parpool(round(0.8*nw));
        
        t0_boot = tic;

        switch ml_method_svm
    
            case 'predict'
                [~, bs_stats] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', 1, ...
                    'bootsamples', boot_n_svm, 'error_type', 'mcr', parallelstr_svm, 'verbose', 0);
                
            case 'oofmridataobj'
                [~ , bs_stats] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', 1, ...
                    'bootsamples', boot_n_svm,  'C', bo_C, ...
                    'error_type', 'mse', parallelstr_svm, 'verbose', 0);
                
        end

        t_end_boot = toc(t0_boot);
        disp('Cumulative run time:');
        toc(t_end_boot); 
        
        
    end % if loop bootstrap
    
    
    % STORE STATISTIC OBJECTS IN CELL ARRAY
    % --------------------------------------------------------------------
    
%     stats.weight_obj = enforce_variable_types(stats.weight_obj);
    svm_stats_results{c} = stats;

    
    if dosearchlight_svm
        searchlight_svm_stats{c} = searchlight_stats;
        searchlight_svm_objs{c} = searchlight_svm_obj;
    end
    
    if dobootstrap_svm
        bootstrap_svm_stats{c} = bs_stats;
    end
        
    if exist('svmmask', 'var')
%         svm_stats_results{c}.mask = svmmask;
        svm_stats_results{c}.maskname = maskname_svm;
        if dosearchlight_svm
            searchlight_svm_stats{c}.maskname = maskname_svm;
        end
        if dobootstrap_svm
            bootstrap_svm_stats{c}.maskname = maskname_svm;
        end
            
    end

end  % loop over contrasts


%% SAVE RESULTS
% -------------------------------------------------------------------------

if dosavesvmstats
    
    fprintf('\n\n');
    printhdr('Saving SVM results');
    fprintf('\n\n');
    
    if exist('maskname_short', 'var')
        savefilenamedata = fullfile(resultsdir, ['svm_stats_results_contrasts_', scaling_string, '_', maskname_short, '_', results_suffix,'.mat']);
    else
        savefilenamedata = fullfile(resultsdir, ['svm_stats_results_contrasts_', scaling_string, '_', results_suffix,'.mat']);
    end
    save(savefilenamedata, 'svm_stats_results', '-v7.3');
    fprintf('\nSaved svm_stats_results_contrasts\n');
    
    if dosearchlight_svm
        save(savefilenamedata, 'searchlight_svm_stats','searchlight_svm_objs','-append');
        fprintf('\nAdded searchlight results to saved svm_stats_results_contrasts\n');
    end
    
    if dobootstrap_svm
        save(savefilenamedata, 'bootstrap_svm_stats','-append');
        fprintf('\nAdded bootstrapped results to saved svm_stats_results_contrasts\n');
    end
    
end


%% SUBFUNCTIONS FOR NORMALIZING
% -------------------------------------------------------------------------

function cat_obj = normalize_each_subject_by_l2norm(cat_obj, condition_codes)
% normalize_each_subject_by_l2norm; can help with numerical scaling and inter-subject scaling diffs
% Sometimes condition differences are very small relative to baseline
% values and SVM is numerically unstable. If so, re-normalizing each
% subject can help.

disp('Normalizing images for each subject by L2 norm of Condition 1 image');

wh = find(condition_codes == 1);

wh2 = find(condition_codes == 2);

% nv: normalization values, to be determined from condition 1 and applied
% to conditions 1 and 2.  This keeps same scaling applied to both
% conditions, for each participant

nv = zeros(size(wh));

for i = 1:length(wh)
    
    nv(i) = norm(cat_obj.dat(:, wh(i)));

    % do normalization
    cat_obj.dat(:, wh(i)) = cat_obj.dat(:, wh(i)) ./ nv(i);
    
    cat_obj.dat(:, wh2(i)) = cat_obj.dat(:, wh2(i)) ./ nv(i);
     
end

end


function cat_obj = normalize_images_by_l2norm(cat_obj)
% normalize_images_by_l2norm; can help with numerical scaling and inter-subject scaling diffs
% Sometimes condition differences are very small relative to baseline
% values and SVM is numerically unstable. If so, re-normalizing each
% subject can help.
%
% This version normalizes each image separately, not each subject/pair

disp('Normalizing images for each image by L2 norm');
cat_obj = rescale(cat_obj, 'l2norm_images');


end
