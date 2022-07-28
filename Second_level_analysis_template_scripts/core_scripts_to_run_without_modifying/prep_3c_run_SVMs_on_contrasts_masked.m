%%% prep_3c_run_SVMs_on_contrasts_masked.m

% USAGE
%
% This script 
% 1) runs second⁻level (i.e. across subjects) support vector machines
% for each within-subject CONTRAST registered in DAT.contrasts
% 2) plots montages of the uncorrected results
% 3) saves the results using standard naming and location
%
% - To specify analysis options, run a2_set_default_options
% - To get results reports, run c2a_second_level_regression
%
% OPTIONS SPECIFIED IN a2_set_default_options
%
% - maskname_svm: default use of sparse gray matter mask; maskdir now defined in a_set_up_paths_always_run_first script; if you do not want to mask, change to []; if you want to use a custom mask, put it in maskdir and change name here.
% - dosubjectnorm: default false; normalize_each_subject_by_l2norm; normalizes images for each subject by L2 norm of Condition 1 image; can help with numerical scaling and inter-subject scaling diffs
% - doimagenorm: default false; normalize_images_by_l2norm; normalizes each image separately, not each subject/pair
% - dozscoreimages: default false; Z-score each input image, removing image mean and forcing std to 1. Removes overall effects of image intensity and scale. Can be useful across studies but also removes information. Use judiciously. lukasvo76: corresponds to 'scaled' in myscaling_glm option in prep_3a
% - dosavesvmstats: default true; Save statistics and weight map objects for SVM contrasts
% - dobootstrap: default false; Takes a lot of time, hence only use true for final analysis, since this takes a lot of time, especially if boot_n is set to the default 10k samples
% - boot_n: default 5000; number of bootstrap samples       Reduce number for quick results
% - parallelstr: default 'parallel'; parallel processing for bootstrapping.   'parallel' or 'noparallel'
% - holdout_set_method: 'group', or 'onesample'
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
% - holdout_set_type: 'kfold', or 'leave_one_subject_out': default 'kfold'; choose between kfold or leave one subject out cross-validation - the latter is not recommended
% - nfolds: default 5, number of folds for kfold CV
%
% LAST THREE OPTIONS ADDED BY @LUKASVO76 JULY 2022
%
% 
% OPTIONS TO BE SPECIFIED IN THIS SCRIPT
%
% -results_suffix: name to add to results file to specify model in case of multiple models, e.g. 'masked_gray_matter'
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% prep_3c_run_SVMs_on_contrasts_masked.m         v2.2
% last modified: 2022/07/28


%% SETTINGS
%--------------------------------------------------------------------------

% options to be specified here
 
results_suffix = ''; % do not delete, leave empty if not needed

% options set in a2_set_default_options

options_needed = {'dosavesvmstats', 'dobootstrap', 'boot_n', 'holdout_set_method', 'holdout_set_type', 'nfolds'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {true false 5000 'onesample' 'kfold' 5};          % defaults if we cannot find info in a2_set_default_options at all 

plugin_get_options_for_analysis_script

% specify which montage to add title to

whmontage = 5; % see region.montage docs


%% CHECK SPIDER TOOLBOX AND START CLOCK
% -------------------------------------------------------------------------

spath = which('use_spider.m');
if isempty(spath)
    error('Spider toolbox not found on Matlab path, clone CanlabCore repo and add to Matlab path to prevent prediction from breaking')
end

if dobootstrap
    svmtime = tic; 
end


%% GET MASK AND MASKNAME
% -------------------------------------------------------------------------

if exist('maskname_svm', 'var') && ~isempty(maskname_svm) 
    svmmask = fmri_mask_image(maskname_svm, 'noverbose');
    [~,maskname_short] = fileparts(maskname_svm);
    mask_string = sprintf('within mask %s', maskname_short);
    
else
    fprintf('\nNo mask found; using full original image data\n');

end


%% RUN SUPPORT VECTOR MACHINES FOR EACH CONTRAST
% -------------------------------------------------------------------------

kc = size(DAT.contrasts, 1);

svm_stats_results = cell(1, kc);

for c = 1:kc
    
    analysisname = DAT.contrastnames{c};
    printhdr(analysisname);
    
    mycontrast = DAT.contrasts(c, :);
    wh = find(mycontrast); % wh is which conditions have non-zero contrast weights
    
    % CREATE COMBINED DATA OBJECT WITH INPUT IMAGES FOR BOTH CONDITIONS
    % ---------------------------------------------------------------------
    [cat_obj, condition_codes] = cat(DATA_OBJ{wh}); 
    cat_obj = remove_empty(cat_obj); % @lukasvo76: added this as the cat function includes replace_empty on the objects, which causes problems later on with the stats_object output of predict; I also guess we do not want to run the SVMs on these "artificial" zeroes?
    
    % NORMALIZE IF SPECIFIED IN OPTIONS
    % ---------------------------------------------------------------------
    
    % NORMALIZE BY L2NORM
    
    % possibly normalize_each_subject_by_l2norm; can help with numerical scaling and inter-subject scaling diffs
    % Sometimes condition differences are very small relative to baseline
    % values and SVM is numerically unstable. If so, re-normalizing each
    % subject can help.
    
    scaling_string = 'no_scaling';
    
    if exist('dosubjectnorm', 'var') && dosubjectnorm
        cat_obj = normalize_each_subject_by_l2norm(cat_obj, condition_codes);
        scaling_string = 'scaling_l2norm_subjects';
        fprintf('\nNormalizing condition images for each subject by L2 norm of Condition 1 image before SVM\n')
        
    elseif exist('doimagenorm','var') && doimagenorm  
        cat_obj = normalize_images_by_l2norm(cat_obj);
        scaling_string = 'scaling_l2norm_conditions';
        fprintf('\nNormalizing each condition image for each subject by L2 norm before SVM\n')
        
    end
    
    % Z-SCORE
    
    % Z-score each input image, removing image mean and forcing std to 1.
    % Removes overall effects of image intensity and scale. Can be useful
    % across studies but also removes information. Use judiciously.
    
    if exist('dozscoreimages', 'var') && dozscoreimages
        fprintf('\nZ-scoring each condition image for each subject before SVM\n');
        scaling_string = 'scaling_z_score_conditions';
        cat_obj = rescale(cat_obj, 'zscoreimages');
        
    end
    
    % APPLY MASK IF SPECIFIED IN OPTIONS
    %----------------------------------------------------------------------
    if exist('svmmask', 'var')
        fprintf('\nMasking data with %s\n',maskname_short);
        cat_obj = apply_mask(cat_obj, svmmask);
        cat_obj.mask_descrip = maskname_svm;
        
    else
        fprintf('\nNo mask found; using full existing image data\n');
        
    end
    
    % FORMAT AND ATTACH OUTCOME: 
    % 1, -1 for the two conditions in the contrast
    % DEFINE HOLDOUT SETS:
    % use plugin scripts according to holdout_set_method
    % --------------------------------------------------------------------
    
    % NOTE:
    % assume that subjects are in same position in each input file!
    
    switch holdout_set_method
        
        case 'group' 
            plugin_get_holdout_sets_balanced_groups; % @lukasvo76: built in this option to manually define your holdout sets balancing for group variable, wrote a new plugin based on @bogpetre's walkthrough code for single-trials and between-within MVPA
            cat_obj.Y = outcome_value;
        
        case 'onesample'
            plugin_get_holdout_sets; % @lukasvo76: this is the original CANlab code which works fine if you do not have to balance your holdout sets for a group variable
            cat_obj.Y = outcome_value;
            
        otherwise
            error('\ninvalid option "%s" defined in holdout_set_method variable, choose between "group" and "onesample"\n',holdout_set_method);
    
    end
    
    % SANITY CHECK ON OUTCOME
    % --------------------------------------------------------------------
    
    if all(cat_obj.Y > 0) || all(cat_obj.Y < 0)
        % Only positive or negative weights - nothing to compare
        printhdr(' Only positive or negative weights - nothing to compare');
        
        continue    
        
    end
    
    % RUN PREDICTIVE SVM MODEL
    % --------------------------------------------------------------------
    if dobootstrap
        [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, ...
            'bootsamples', boot_n, 'error_type', 'mcr', parallelstr);
        % Threshold, if possible - can re-threshold later with threshold() method
%         stats.weight_obj = threshold(stats.weight_obj, .05, 'unc'); %
%         @lukasvo76: commented out since we want to threshold flexibly at
%         a later stage in c2_SVM_contrasts_masked
        
    else
        [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, ...
            'error_type', 'mcr', parallelstr);
    
    end
    
    % PLOT MONTAGE OF UNTHRESHOLDED RESULTS
    % ---------------------------------------------------------------------

    fprintf ('\nShowing unthresholded SVM results, : %s\nEffect: %s\n\n', analysisname, mask_string);
    
    r = region(stats.weight_obj);
    
    o2 = montage(r, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
    o2 = title_montage(o2, whmontage, [analysisname ' unthresholded ' mask_string]);

    figtitle = sprintf('%s_unthresholded_montage_%s_%s', analysisname, scaling_string, mask_string);
    set(gcf, 'Tag', figtitle, 'WindowState','maximized');
    drawnow, snapnow;
        if save_figures
            plugin_save_figure;
        end
    clear o2, clear figtitle
    
    
    % STORE STATISTIC OBJECTS IN CELL ARRAY
    % --------------------------------------------------------------------
    
%     stats.weight_obj = enforce_variable_types(stats.weight_obj);
    svm_stats_results{c} = stats;
        
    if exist('svmmask', 'var')
        svm_stats_results{c}.mask = svmmask;
        svm_stats_results{c}.maskname = maskname_svm;
    end
    
    if dobootstrap
        disp('Cumulative run time:');
        toc(svmtime); 
    
    end

end  % loop over contrasts


%% SAVE RESULTS
%--------------------------------------------------------------------------

if dosavesvmstats
    savefilenamedata = fullfile(resultsdir, ['svm_stats_results_contrasts_', scaling_string, '_', results_suffix,'.mat']);
    save(savefilenamedata, 'svm_stats_results', '-v7.3');
    fprintf('\nSaved svm_stats_results_contrasts\n');
    
end


%% SUBFUNCTIONS FOR NORMALIZING
%--------------------------------------------------------------------------

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
