%% prep_3a_run_second_level_regression_and_save.m
%
%
% USAGE
%
% This script
%
% 1) runs second⁻level (i.e. across subjects) regression analyses
%   for each within-subject CONTRAST or CONDITION registered in the DAT
%   structure, either
%
%       a) voxel-wise, calling CANlab's regress() function under the hood,
%       including its robust regression option if specified
%
%       b) parcel-wise, calling CANlab's robfit_parcelwise() function under the
%       hood, which is robust by default
%
%   the option to convert t-maps into BayesFactor maps using CANlab's
%   estimateBayesFactor() function is built in- see walkthrough 
%   https://canlab.github.io/_pages/EmoReg_BayesFactor_walkthrough/EmoReg_BayesFactor_walkthrough.html
%
% 2) runs cross-validated MVPA regression models predicting continuous
%   covariates if desired using CANlab's predict() function
%
% 3) saves the results using standard naming and location
% 
% Run this script with Matlab's publish function to generate html report of results:
% publish('prep_3a_run_second_level_regression_and_save','outputDir',htmlsavedir)
%
% To get results reports after thresholding, publish
% c2a_second_level_regression
%
%
% OPTIONS
%
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options (for example
% voxel- and parcelwise regression), you can make a copy of this script with
% a letter index (e.g. _s6a_) and change the default option here
%
% - dorobust : robust regression or OLS (true/false)
% - dorobfit_parcelwise: voxel- or parcelwise regression (true/false) - % OPTION ADDED BY @LUKASVO76 MAY 2022
%       - csf_wm_covs: true adds global wm & csf regressors at second level
%       - remove_outliers: true removes outlier images/subjects based on mahalanobis distance 
% - myscaling_glm: 'raw', 'scaled', or 'scaled_contrasts' (defined in a2_set_..., image scaling done in prep_2_... and prep_3_... data load)
% - maskname_glm: 
%       - default use of sparse gray matter mask
%       - model-specific maskdir defined in a_set_up_paths_always_run_first script
%       - if you do not want to mask, change to []
%       - if you want to use a custom mask, put it in maskdir and change name here
%       - only used for visualization of uncorrected results in this script
% - design_matrix_type: 'group', 'custom', or 'onesample'
%                       Group: use DAT.BETWEENPERSON.group
%                       Custom: use all columns of table object DAT.BETWEENPERSON.(mygroupnamefield){c}
%                               NOTE: you can flexibly use the columns as
%                                       covariates by editing line 308 in the code below
%                               EXAMPLE: if you only want to use the first
%                                       column, you can change to
%                                   table_obj = DAT.BETWEENPERSON.(mygroupnamefield){c}(:,1);
%                       Onesample: use constant (i.e. intercept) only
%
%       - 'group' option 
%           Assuming that groups are concatenated in contrast image lists, and
%           regressor values of 1 or -1 will specify the group identity for each
%           image. Requires DAT.BETWEENPERSON.group field specifying group membership for
%           each image.
%
%       - 'custom' option: 
%           Can enter a multi-column design matrix for each contrast
%           Design matrix can be different for each contrast
%
%       - 'onesample' option:
%           Only adds intercept, hence performs a one-sample t-test on contrast
%           images across all subjects, similarly to c_univariate_contrast_maps_
%           scripts, but with more flexible options including scaling and robustfit
%   
%       OPTION ADDED BY @LUKASVO76 MAY 2022
%
%       NOTE: To set up group and custom variables, see prep_1b_prep_behavioral_data
%
% - doBayes: convert t-maps into Bayes Factors 
%
%       OPTION ADDED BY @LUKASVO76 JANUARY 2023
%
% - domvpa_reg_cov: run MVPA regression model to predict covariate levels from (between-subject) brain data using CANlab's predict() function
%
%   NOTE: THIS OPTION ONLY APPLIES WHEN DESIGN_MATRIX_TYPE = 'CUSTOM' SINCE OTHERWISE THERE IS NO CONTINUOUS OUTCOME TO PREDICT!
%         TO CLASSIFY GROUPS USING MVPA MODELS, USE SVM SCRIPTS PREP_3C AND C2
%     
%     mvpa_reg_covariate options
%       algorithm_mvpa_reg_cov = 'cv_pcr': will be passed into predict function (help fmri_data.predict for options)
%       holdout_set_method_mvpa_reg_cov: 'group', or 'no_group'
%                                        group: use DAT.BETWEENPERSON.group to balance holdout sets over groups
%                                        no_group: no group factor, stratifies by subject (i.e.leave whole subject out) since data is purely between-subject
%       nfolds_mvpa_reg_cov: number of cross-validation folds for kfold
%       zscore_outcome_mvpa_reg_cov: zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models
%
%       OPTION ADDED BY @LUKASVO76 JANUARY 2023
%
%
% MANDATORY OPTIONS TO BE SPECIFIED IN THIS SCRIPT
%
% - mygroupfieldname: 'contrasts' or 'conditions'
% - results_suffix: name to add to results file to specify in case of multiple versions of model, e.g. 'covariate_rating'
% 
% MANDATORY OPTIONS TO BE SPECIFIED IN THIS SCRIPT ONLY IF DOMVPA_REG_COV = true
%
% - behav_outcome: variable name in DAT.BETWEENPERSON.(mygroupnamefield){:}
% - group_identifier = name of group identifier variable in same table 
%
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% prep_3a_run_second_level_regression_and_save.m         v4.0
% last modified: 2023/01/12


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% SET MANDATORY OPTIONS

mygroupnamefield = 'contrasts'; 
results_suffix = ''; % adds a suffix of your choice to .mat file with results that will be saved

% NOTE: do NOT delete the latter option, leave empty if not needed
% NOTE: do NOT use to add a suffix specifying the regressors, scaling or masking option, this will be added automatically

% MANDATORY OPTIONS IF DOMVPA_COV_REG = true

behav_outcome_mvpa_reg_cov = {'delta_intensity','delta_rating'}; % needs to correspond to variable name(s) in DAT.BETWEENPERSON.(mygroupnamefield){:} AND TO THE ORDER IN WHICH THEY APPEAR THERE

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% GET DEFAULT OPTIONS IF NOT SET IN A2_SET_DEFAULT_OPTIONS

options_needed = {'dorobust', 'dorobfit_parcelwise', 'myscaling_glm', 'design_matrix_type', 'maskname_glm'};
options_exist = cellfun(@exist, options_needed); 

option_default_values = {false, false, 'raw', 'onesample', which('gray_matter_mask_sparse.img')};

plugin_get_options_for_analysis_script;

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% dorobust = true/false;
% dorobfit_parcelwise = true/false;
%   csf_wm_covs = true/false;
%   remove_outliers = true/false;
% myscaling_glm = 'raw'/'scaled'/'scaled_contrasts';
% design_matrix_type = 'onesample'/'group'/'custom';
% doBayes = true/false;
% domvpa_reg_cov = true/false;
%   algorithm_mvpa_reg_cov = ;
%   holdout_set_mvpa_reg_cov = 'group'/'no_group';
%   nfolds_mvpa_reg_cov = x;
%   zscore_outcome_mvpa_reg_cov = true/false;

% SANITY CHECK

if ~strcmpi(design_matrix_type,'custom') && domvpa_reg_cov
    error('\noption "%s" defined in design_matrix_type not compatible with do_mvpa_reg_cov, change design_matrix_type to "custom" or turn off do_mvpa_reg_cov\n', design_matrix_type);
end
    


%% LOAD NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

if ~exist('DATA_OBJ','var') || ~exist('DATA_OBJsc','var')
    
    load(fullfile(resultsdir,'data_objects.mat'));
    load(fullfile(resultsdir,'data_objects_scaled.mat'));
    
end

if ~exist('DATA_OBJ_CON','var') || ~exist('DATA_OBJ_CONsc','var') || ~exist('DATA_OBJ_CONscc','var')
    
    load(fullfile(resultsdir,'contrast_data_objects.mat'));
    
end


%% CHECK REQUIRED DAT FIELDS
% -------------------------------------------------------------------------

% List required fields in DAT, in cell array

if ~strcmpi(design_matrix_type,'onesample')
    
    required_fields = {'BETWEENPERSON', 'contrastnames', 'contrasts' 'contrastcolors', 'conditions', 'colors'};

    ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
    if ~ok_to_run
        return
    end
    
else
    
    required_fields = {'contrastnames', 'contrasts' 'contrastcolors', 'conditions', 'colors'};

    ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
    if ~ok_to_run
        return
    end
    
end


%% MASKING
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('MASKING IMAGES IF REQUESTED IN OPTIONS');
fprintf('\n\n');

if exist('maskname_glm', 'var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')
    [~,maskname_short] = fileparts(maskname_glm);
    mask_string = sprintf('masked with %s', maskname_short);
    glmmask = fmri_mask_image(maskname_glm, 'noverbose'); 
    fprintf('\nMasking results visualization with %s\n\n', maskname_short);
else
    mask_string = sprintf('without masking');
    fprintf('\nShowing results without masking\n\n');
end  


%% RUN SECOND LEVEL REGRESSION FOR EACH CONTRAST
% -------------------------------------------------------------------------

switch mygroupnamefield
    
    case 'contrasts'

        kc = size(DAT.contrasts, 1);
       
        fprintf('\nRUNNING SECOND LEVEL REGRESSIONS ON FIRST LEVEL CONTRASTS\n\n');
        
    case 'conditions'
        
        kc = size(DAT.conditions, 2);
       
        fprintf('\nRUNNING SECOND LEVEL REGRESSIONS ON FIRST LEVEL CONDITIONS\n\n');
        
    otherwise
        
        error('\ninvalid option "%s" defined in mygroupnamefield variable, choose between "contrasts" and "conditions"\n\n',mygroupnamefield)

end

if ~dorobfit_parcelwise
    regression_stats_results = cell(1, kc);
else
    parcelwise_stats_results = cell(1,kc);
end

for c = 1:kc
    
    % GET DESIGN MATRIX FOR THIS CONTRAST OR CONDITION
    % ---------------------------------------------------------------------
    
    switch mygroupnamefield
        
        case 'contrasts'
            fprintf('\n\n');
            printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.contrastnames{c})]);
            fprintf('\n\n');
            
        case 'conditions'
            fprintf('\n\n');
            printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.conditions{c})]);
            fprintf('\n\n');
    
    end
      
    fprintf('\n\n');
    printhdr('Building design matrix');
    fprintf('\n\n');
    
    groupnames_string = 'intercept';
    
    switch design_matrix_type
        
        case 'custom'
            
            % Define design matrix X "design_matrix"
            % Use custom matrix for each condition/contrast
            table_obj = DAT.BETWEENPERSON.(mygroupnamefield){c};
            groupnames = table_obj.Properties.VariableNames;
            X = table2array(table_obj);
            idx_nan = ~isnan(X);
            idx_nan = ~(sum(idx_nan,2) < size(idx_nan,2)); % at least one column of X contains NaN
            imgs_nan = 1:size(X,1);
            imgs_nan = imgs_nan(idx_nan');
            X = X(idx_nan,:);
            
            for name = 1:size(groupnames,2)
                groupnames_string = [groupnames_string, ' ', groupnames{name}];
            end
            
        case 'group'
            
            % Use 'groups' single regressor
            if ~isempty(DAT.BETWEENPERSON.group)
                group = DAT.BETWEENPERSON.group;
                groupnames = {'group'};
                X = group;
                imgs_nan = [];
                groupnames_string = [groupnames_string, ' ', groupnames];
            else
                error('\nGroup not defined in DAT.BETWEENPERSON.group, which is required for option "%s" defined in design_matrix_type\n', design_matrix_type);
            end

        case 'onesample'
            
                % Use intercept only
                switch mygroupnamefield
                    case 'conditions'
                        X = ones((size(DAT.imgs{c},1)),1);
                    case 'contrasts'
                        X = ones((size(DAT.gray_white_csf_contrasts{c},1)),1);
                end
                groupnames = {'intercept'};
                imgs_nan = [];
            
        otherwise
            
            error('\ninvalid option "%s" defined in design_matrix_type variable, choose between "group", "custom", or "onesample"\n\n', design_matrix_type);
            
    end
    
    fprintf('\nREGRESSOR(S): %s\n\n', groupnames_string);
    
    % SELECT DATA FOR THIS CONTRAST/CONDITION
    % ---------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('Scaling data if requested in options');
    fprintf('\n\n');
    
    switch mygroupnamefield
        
        case 'contrasts'
            
            switch myscaling_glm

                case 'raw'
                    fprintf('\nContrast calculated on raw (unscaled) condition images used in second-level GLM\n\n');
                    scaling_string = 'no_scaling';
                    cat_obj = DATA_OBJ_CON{c};
                    if imgs_nan
                        cat_obj = cat_obj.get_wh_image(imgs_nan);
                    end

                case 'scaled'
                    fprintf('\nContrast calculated on z-scored condition images used in second-level GLM\n\n');
                    scaling_string = 'scaling_z_score_conditions';
                    cat_obj = DATA_OBJ_CONsc{c};
                    if imgs_nan
                        cat_obj = cat_obj.get_wh_image(imgs_nan);
                    end

                case 'scaled_contrasts'
                    fprintf('\nl2norm scaled contrast images used in second-level GLM\n\n');
                    scaling_string = 'scaling_l2norm_contrasts';
                    cat_obj = DATA_OBJ_CONscc{c};
                    if imgs_nan
                        cat_obj = cat_obj.get_wh_image(imgs_nan);
                    end

                otherwise
                    error('\nInvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw", "scaled", or "scaled_constrast" given option "%s" defined in mygroupnamefield variable\n\n', myscaling_glm, mygroupnamefield);

            end
            
        case 'conditions'
            
            switch myscaling_glm

                case 'raw'
                    fprintf('\nRaw (unscaled) condition images used in second-level GLM\n\n');
                    scaling_string = 'no_scaling';
                    cat_obj = DATA_OBJ{c};
                    if imgs_nan
                        cat_obj = cat_obj.get_wh_image(imgs_nan);
                    end

                case 'scaled'
                    fprintf('\nZ-scored condition images used in second-level GLM\n\n');
                    scaling_string = 'scaling_z_score_conditions';
                    cat_obj = DATA_OBJsc{c};
                    if imgs_nan
                        cat_obj = cat_obj.get_wh_image(imgs_nan);
                    end

                case 'scaled_contrasts'
                    error('\nInvalid combination of option "%s" defined in myscaling_glm_variable in a2_set_default_options script and option "%s" defined in mygroupnamefield variable, choose between "raw" and "scaled" options\n\n',myscaling_glm,mygroupnamefield);

                otherwise
                    error('\nInvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw",  and "scaled", given option "%s" defined in mygroupnamefield variable\n\n', myscaling_glm, mygroupnamefield);

            end
            
    end % switch mygroupnamefield - contrasts or conditions
    
  
    % FORMAT AND ATTACH DESIGN MATRIX
    % ---------------------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('Checking design matrix');
    fprintf('\n\n');
    
    if ~strcmpi(design_matrix_type,'onesample')
        
        % Confirm design_matrix is 1, -1, or mean-centered
        meancentered = ~(abs(mean(X)) > 1000 * eps);
        effectscoded = all(X == 1 | X == -1 | X == 0, 1);
        isconstant = all(X == mean(X, 1), 1);
        vifs = getvif(X);

        if any(isconstant)
            
            fprintf('\n');
            warning('An intercept appears to be added manually. Do not include an intercept - it will be added automatically.');
            warning('Skipping this contrast.');
            fprintf('\n');
            
            continue
        end

        % Report
        design_table = table;
        design_table.Mean = mean(X)';
        design_table.Var = var(X)';
        design_table.EffectsCode = effectscoded';
        design_table.VIF = vifs';
        design_table.Properties.RowNames = groupnames';
        disp(design_table)
        disp(' ');

        if any(~meancentered & ~effectscoded)
            fprintf('\n');
            warning('Some columns are not mean-centered or effects coded. Intercept may not be interpretable');
            fprintf('\nColumns: ')
            fprintf('%d \n', find(~meancentered & ~effectscoded));
        else
            fprintf('\nChecked OK: All columns mean-centered or are effects-coded [1 -1 0]\n\n');
        end

        if any(vifs > 2)
            fprintf('\n');
            warning('Some regressors have high variance inflation factors. Parameters might be poorly estimated or uninterpretable.');
            fprintf('\n');
        else
            fprintf('\nChecked OK: VIFs for all columns are < 2\n\n');
        end
    
    else
        design_table = table;
        design_table.Mean = mean(X)';
        design_table.Var = var(X)';
        
    end % if loop design_matrix_type
    
    cat_obj.X = X;
    
    % SANITY CHECK ON REGRESSORS, SKIP CONTRAST IF NEEDED
    % ---------------------------------------------------------------------
    
     if ~strcmpi(design_matrix_type,'onesample')
        
        for col = 1:size(cat_obj.X,2)
    
            if all(cat_obj.X(:,col) > 0) || all(cat_obj.X(:, col) < 0)
                % Only positive or negative weights - nothing to compare

                fprintf('\n');
                warning('Only positive or negative regressor values - bad design, please check');
                fprintf('\n');

                continue
            end
            
        end
        
    end
    
    % ADD COVARIATE TO .Y FIELD OF CAT_OBJ FOR MVPA IF REQUESTED
    % ---------------------------------------------------------------------
    
    if domvpa_reg_cov
        
        for covar = 1:size(behav_outcome_mvpa_reg_cov,2)
            
            if ~ismember(behav_outcome_mvpa_reg_cov{covar},table_obj.Properties.VariableNames)
                error('\nCovariate "%s" defined in behav_outcome not present in DAT.BETWEENPERSON.(mygroupfield){%d}, please correct before proceeding\n',behav_outcome_mvpa_reg_cov{covar},covar);
            end
            
            if size(behav_outcome_mvpa_reg_cov,2) > 1
            
                if ~strcmp(behav_outcome_mvpa_reg_cov{covar},table_obj.Properties.VariableNames{covar})
                    error('\norder of covariates in behav_outcome not consistent with DAT.BETWEENPERSON.(mygroupfield){%d}, please correct before proceeding\n');
                end
                
                mvpa_data_objects{covar} = cat_obj;
                mvpa_data_objects{covar}.Y = table2array(table_obj(:,covar));
                mvpa_data_objects{covar}.Y_names = behav_outcome_mvpa_reg_cov{covar};
   
                
            else 
               
                idx= ismember(table_obj.Properties.VariableNames,behav_outcome_mvpa_reg_cov{covar});
                [~,covar_idx] = find(idx);
                
                mvpa_data_objects{covar} = cat_obj;
                mvpa_data_objects{covar}.Y = table2array(table_obj(:,covar_idx));
                mvpa_data_objects{covar}.Y_names = behav_outcome_mvpa_reg_cov{covar};
                
            end
            
        end
        
    end
        
        
    % RUN GLM MODEL
    % ---------------------------------------------------------------------
    
    % VOXEL-WISE
    
    if ~dorobfit_parcelwise
        
        if dorobust
            robuststring = 'robust';
            regresstime = tic;
        else
            robuststring = 'norobust';
        end
        
        fprintf('\n\n');
        printhdr(['Running voxel-wise ', robuststring ' regression']);
        fprintf('\n\n');

        if ~strcmpi(design_matrix_type,'onesample')
            % out.t has t maps for all regressors, intercept is last
            switch mygroupnamefield
                case 'contrasts'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.contrastnames{c}, 'variable_names', groupnames, 'nodisplay'); % trick to get unthresholded maps
                case 'conditions'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.conditions{c}, 'variable_names', groupnames, 'nodisplay');
            end
        else
            % out.t has t maps for intercept only
            switch mygroupnamefield
                case 'contrasts'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.contrastnames{c}, 'variable_names', groupnames, 'nointercept', 'nodisplay');
                case 'conditions'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.conditions{c}, 'variable_names', groupnames, 'nointercept', 'nodisplay');
            end
        end

        % make sure variable types are right data formats
        regression_stats.design_table = design_table;
        regression_stats.t = enforce_variable_types(regression_stats.t);
        regression_stats.b = enforce_variable_types(regression_stats.b);
        regression_stats.df = enforce_variable_types(regression_stats.df);
        regression_stats.sigma = enforce_variable_types(regression_stats.sigma);
        
        % calculate Bayes Factors from t-maps if requested and add to results
        if doBayes
            
            fprintf('\n\n');
            printhdr('Calculating Bayes Factor maps');
            fprintf('\n\n');
            
            N = single(sum(~(isnan(cat_obj.dat') | cat_obj.dat' == 0) , 1)); % code from fmri_data.ttest to correct N in regressions_stat.t to make estimateBayesFactor function work on regress() output
            t_for_Bayes = regression_stats.t;
            t_for_Bayes.N = N';
            regression_stats.BF = estimateBayesFactor(t_for_Bayes,'t');
        end

        % add analysis name, regressor names and other meta-data
        switch mygroupnamefield
            case 'contrasts'
                regression_stats.contrastname = DAT.contrastnames{c};
                regression_stats.contrast = DAT.contrasts(c, :);
                regression_stats.analysis_name = DAT.contrastnames{c};
            case 'conditions'
                regression_stats.contrastname = DAT.conditions{c};
                regression_stats.contrast = 1;
                regression_stats.analysis_name = DAT.conditions{c};
        end

        % add names for variables 
        if ~strcmpi(design_matrix_type,'onesample')
            regression_stats.variable_names = [groupnames {'intercept'}];
        else
            regression_stats.variable_names = groupnames;
        end

        % PLOT ORTHVIEWS (MASKED IF SPECIFIED IN MASKNAME_GLM OPTION)
        % --------------------------------------------------------------------
        
        fprintf ('\nORTHVIEWS GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR(S): %s, %s, SCALING: %s\n\n', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
        
        t = threshold(regression_stats.t,.05,'unc');
            if maskname_short
                t = apply_mask(t,glmmask);
            end
        orthviews(t);
            for kk = 1:length(regression_stats.variable_names)
                switch mygroupnamefield
                    case 'contrasts'
                        spm_orthviews_name_axis([regression_stats.variable_names{kk},' ',DAT.contrastnames{c}], kk);
                    case 'conditions'
                        spm_orthviews_name_axis([regression_stats.variable_names{kk},' ',DAT.conditions{c}], kk);
                end
            end
        drawnow;snapnow;
        
        if doBayes
            fprintf ('\nORTHVIEWS BAYESIAN GLM RESULTS AT |BF| > 3, EFFECT: %s, REGRESSOR(S): %s, %s, SCALING: %s\n\n', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
        
            BF = threshold(regression_stats.BF,[-2.1972 2.1972],'raw-outside');
                if maskname_short
                    BF = apply_mask(BF,glmmask);
                end
            orthviews(BF);
                for kk = 1:length(regression_stats.variable_names)
                    switch mygroupnamefield
                        case 'contrasts'
                            spm_orthviews_name_axis([regression_stats.variable_names{kk},' ',DAT.contrastnames{c}], kk);
                        case 'conditions'
                            spm_orthviews_name_axis([regression_stats.variable_names{kk},' ',DAT.conditions{c}], kk);
                    end
                end
        drawnow;snapnow;
            
        end

        % KEEP RESULTS OBJECTS IN CELL ARRAY FOR SAVING
        % ---------------------------------------------------------------------

        regression_stats_results{c} = regression_stats;
        
        if exist(maskname_glm,'file')
            regression_stats_results{c}.maskname = maskname_glm;
        end
        
        if dorobust
            fprintf('\nCumulative run time:\n'), toc(regresstime); 
        end
        
    % PARCEL-WISE    
        
    else
        
        fprintf('\n\n');
        printhdr('Running parcel-wise robust regression');
        fprintf('\n\n');
        
        if csf_wm_covs && remove_outliers
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',true,'doplot',false);
        elseif csf_wm_covs && ~remove_outliers
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',false,'doplot',false);
        elseif ~csf_wm_covs && remove_outliers
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',false,'remove_outliers',true,'doplot',false);
        else
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'doplot',false);
        end
        
        % calculate Bayes Factors from t-maps if requested and add to results
        if doBayes
            
            fprintf('\n\n');
            printhdr('Calculating Bayes Factor maps');
            fprintf('\n\n');
            
            N = single(size(cat_obj.dat,2).*(ones(size(parcelwise_stats.t_obj.dat,1),1))); % code from fmri_data.ttest to correct N in regressions_stat.t to make estimateBayesFactor function work on regress() output
            t_for_Bayes = parcelwise_stats.t_obj;
            t_for_Bayes.N = N;
            parcelwise_stats.BF = estimateBayesFactor(t_for_Bayes,'t');
            
        end
        
        % add design table
        parcelwise_stats.design_table = design_table;

        % add analysis name, regressor names, and other meta-data
        switch mygroupnamefield
            case 'contrasts'
                parcelwise_stats.contrastname = DAT.contrastnames{c};
                parcelwise_stats.contrast = DAT.contrasts(c, :);
                parcelwise_stats.analysis_name = DAT.contrastnames{c};
            case 'conditions'
                parcelwise_stats.contrastname = DAT.conditions{c};
                parcelwise_stats.contrast = 1;
                parcelwise_stats.analysis_name = DAT.conditions{c};
        end

        % add names for variables 
        if ~strcmpi(design_matrix_type,'onesample')
            parcelwise_stats.variable_names = [groupnames {'Intercept'}];
        else
            parcelwise_stats.variable_names = groupnames;
        end
        
        % PLOT PARCELWISE SPECIFIC WEIGHTS AND DIAGNOSTICS
        % --------------------------------------------------------------------
        fprintf('\nPlotting parcel weights and diagnostics\n\n');
        
        create_figure('parcelwise weights and metrics', 2, 2);
        set(gcf, 'WindowState','maximized');
        xlabel('Image'); ylabel('Weights');
        errorbar(mean(parcelwise_stats.weights), std(parcelwise_stats.weights), 'bo', 'MarkerFaceColor', [0 0 .5])
        title('Mean weights across parcels (s.d. error bars) per image');
        axis tight; 

        subplot(2, 2, 2);
        imagesc(parcelwise_stats.weights);
        xlabel('Image'); ylabel('Parcel');
        title('Weights by parcel');
        colorbar;
        axis tight; set(gca, 'YDir', 'Reverse');

        subplot(2, 2, 3);
        xlabel('Image'); ylabel('Z(Weights)');
        errorbar(zscore(mean(parcelwise_stats.weights)), ste(parcelwise_stats.weights), 'bo-', 'MarkerFaceColor', [0 0 .5], 'LineWidth', 2)
        title('Mean weights (s.e. error bars) and quality metrics');
        plot(zscore(parcelwise_stats.individual_metrics.gm_L1norm), 'LineWidth', 2);
        plot(zscore(parcelwise_stats.individual_metrics.csf_L1norm), 'LineWidth', 2);
        plot(zscore(parcelwise_stats.ind_quality_dat.Mahal_corr), 'LineWidth', 2);
        plot(zscore(parcelwise_stats.ind_quality_dat.Mahal_cov), 'LineWidth', 2);
        legend({'Z(Weights)' 'Z(GM L1 norm)' 'Z(CSF L1 norm)' 'Mahal corr dist' 'Mahal cov dist'});
        axis tight; 

        % mark off who are outliers
        wh_out = find(parcelwise_stats.outliers_uncorr);
            for i = 1:length(wh_out)

                hh = plot_vertical_line(wh_out(i));
                set(hh, 'Color', 'r', 'LineStyle', '--');

                if i == 1
                        legend({'Z(Weights)' 'Z(GM L1 norm)' 'Z(CSF L1 norm)' 'Mahal corr dist' 'Mahal cov dist' 'Mah. outliers p<.05 uncor'});
                end
            end

        subplot(2, 2, 4)
        plot_correlation_matrix(parcelwise_stats.datmatrix, 'dofigure', false);
        title('inter-parcel correlations across images');
        drawnow, snapnow;
        

        % PLOT MONTAGE (MASKED IF SPECIFIED IN MASKNAME_GLM OPTION)
        % ---------------------------------------------------------------------
        
        fprintf ('\nMONTAGE PARCELWISE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR(S): %s, %s, SCALING: %s\n\n', parcelwise_stats.analysis_name, groupnames_string, mask_string, scaling_string);
        
        num_effects = size(parcelwise_stats.t_obj.dat, 2); % number of regressors
        o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        for j = 1:num_effects

            tj = get_wh_image(parcelwise_stats.t_obj, j);
                if maskname_short
                    tj = apply_mask(tj, glmmask);
                end
            tj = threshold(tj, .05, 'unc'); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [parcelwise_stats.analysis_name ' ' parcelwise_stats.variable_names{j} ' ' mask_string ' ' scaling_string]);

        end

        figtitle = sprintf('%s_05_unc_montage_%s_%s_%s', parcelwise_stats.analysis_name, groupnames_string, mask_string, scaling_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;
            if save_figures
                plugin_save_figure;
            end
        clear o2, clear figtitle, clear j, clear tj
        
        if doBayes
           
            fprintf ('\nMONTAGE BAYESIAN PARCELWISE GLM RESULTS AT |BF| > 3, EFFECT: %s, REGRESSOR(S): %s, %s, SCALING: %s\n\n', parcelwise_stats.analysis_name, groupnames_string, mask_string, scaling_string);
        
            num_effects = size(parcelwise_stats.BF.dat, 2); % number of regressors
            o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

            for j = 1:num_effects

                BFj = get_wh_image(parcelwise_stats.BF, j);
                    if maskname_short
                        BFj = apply_mask(BFj, glmmask);
                    end
                BFj = threshold(BFj,[-2.1972 2.1972],'raw-outside'); 

                o2 = addblobs(o2, region(BFj), 'wh_montages', (2*j)-1:2*j);
                o2 = title_montage(o2, 2*j, [parcelwise_stats.analysis_name ' ' parcelwise_stats.variable_names{j} ' ' mask_string ' ' scaling_string]);

            end

            figtitle = sprintf('%s_BF_3_montage_%s_%s_%s', parcelwise_stats.analysis_name, groupnames_string, mask_string, scaling_string);
            set(gcf, 'Tag', figtitle, 'WindowState','maximized');
            drawnow, snapnow;
                if save_figures
                    plugin_save_figure;
                end
            clear o2, clear figtitle, clear j, clear BFj
            
        end

        % KEEP RESULTS OBJECTS IN CELL ARRAY FOR SAVING
        % ---------------------------------------------------------------------

        parcelwise_stats_results{c} = parcelwise_stats;
        
        if exist(maskname_glm,'file')
            parcelwise_stats_results{c}.maskname = maskname_glm;
        end
        
    end % if loop voxel- versus parcelwise
    
    
    % RUN MVPA MODEL IF REQUESTED IN OPTIONS
    % ---------------------------------------------------------------------
    
    if domvpa_reg_cov
        
        for covar = 1:size(mvpa_data_objects,2)
            
            mvpa_dat = mvpa_data_objects{covar};
            
            % DATA VISUALIZATION
            
            fprintf('\n\n');
            printhdr('PLOTTING DATA');
            fprintf('\n\n');

                % con images

                h1=figure;

                    for subj = 1:size(mvpa_dat.dat,2)
                        this_subj_dat = mvpa_dat.dat(:,subj);
                        q(subj,:) = quantile(this_subj_dat(:),[0.025,0.5,0.975]);
                        mu = mean(mean(this_subj_dat(:)));
                        sd = std(this_subj_dat(:));
                        h1 = plot([mu-sd, mu+sd],[subj,subj],'-');
                        hold on;
                        h2 = plot(mu,subj,'o');
                        h2.Color = h1.Color;
                    end

                box off
                title('Distribution of con weights');
                xlabel('\beta');
                ylabel('Subject');
                hold off

                p = get(gcf,'Position');
                set(gcf,'Position',[p(1:2),1024,2048],'WindowState','Maximized');
                drawnow, snapnow;

                clear subj

                % behavioral outcome

                b1=figure;

                hold off;
                b1=histogram(mvpa_dat.Y);
                box off
                title(['Histogram of ' behav_outcome_mvpa_reg_cov{covar}]);
                xlabel(behav_outcome_mvpa_reg_cov{covar});
                ylabel('n(observations)');
                set(gcf,'WindowState','Maximized');
                drawnow, snapnow;
            
            % RUN MODEL
            
                % cross-validation fold selection

                switch holdout_set_method_mvpa_reg_cov

                    case 'no_group'

                        if ~isempty(DAT.BETWEENPERSON.group)
                            fprintf('\n');
                            warning('DAT.BETWEENPERSON.group defines a grouping factor, please change holdout_set_method_mvpa_reg_cov to "group" for correctly stratified CV fold selection.');
                            fprintf('\n');
                        end

                        cv=cvpartition(size(mvpa_dat.dat,2),'KFold',nfolds_mvpa_reg_cov);
                        fold_labels = zeros(size(mvpa_dat.dat,2),1);
                            for subj = 1:cv.NumTestSets
                                fold_labels(cv.test(subj)) = subj;
                            end
                        clear subj

                    case 'group'

                        if ~isempty(DAT.BETWEENPERSON.group)
                            group = DAT.BETWEENPERSON.group;
                        else
                            error('\nGroup not defined in DAT.BETWEENPERSON.group, which is required for option "%s" chosen in holdout_set_method_mvpa_reg_cov\n', holdout_set_method_mvpa_reg_cov);
                        end

                        cv = cvpartition(group, 'KFold',nfolds_mvpa_reg_cov);
                            fold_labels = zeros(size(mvpa_dat.dat,2),1);
                            for subj = 1:cv.NumTestSets
                                fold_labels(cv.test(subj)) = subj;
                            end
                        clear subj

                end % switch holdout set method

                % fit model

                t0 = tic;

                [mvpa_cverr, mvpa_stats, mvpa_optout] = predict(mvpa_dat, 'algorithm_name', algorithm_mvpa_reg_cov, ...
                            'nfolds', fold_labels, 'error_type', 'mse', 'parallel', 'verbose', 0);

                t_end = toc(t0);  
            
            % VISUALIZE UNTHRESHOLDED RESULTS
            
            fprintf('\n\n');
            printhdr('VISUALIZING RESULTS');
            fprintf('\n\n');

                % plot observed versus predicted

                fprintf('\n\n');
                printhdr('Plotting observed versus predicted');
                fprintf('\n\n');

                fprintf('\n%s r = %0.3f\n\n', algorithm_mvpa_reg_cov, corr(mvpa_stats.yfit, mvpa_dat.Y));
                
                observed = mvpa_dat.Y;
                predicted = mvpa_stats.yfit;
                tbl = table(observed, predicted);
                mdl = fitlm(tbl,'predicted ~ observed');
                
                figure
                
                plot(mdl);
                xlabel({['Observed ' behav_outcome_mvpa_reg_cov{covar}]}); ylabel({['Estimated ' behav_outcome_mvpa_reg_cov{covar}],'(cross validated)'})

                set(gcf,'WindowState','Maximized');
                drawnow, snapnow;

                % plot montage of unthresholded weights

                fprintf('\n\n');
                printhdr('Plotting unthresholded weight maps');
                fprintf('\n\n');

                whmontage = 5;

                fprintf ('\nSHOWING UNTHRESHOLDED %s RESULTS, %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_st), mask_string, myscaling_mvpa_reg_st);

                figure

                o2 = canlab_results_fmridisplay([], 'compact', 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

                w = region(mvpa_stats.weight_obj);

                o2 = addblobs(o2, w);
                o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' unthresholded ' mask_string]);

                figtitle = sprintf('%s_unthresholded_montage_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;

                clear w, clear o2, clear figtitle
                    
            % KEEP RESULTS IN CELL ARRAY FOR SAVING

            mvpa_stats_and_maps{c,covar} = mvpa_stats;

        end % for loop over covariates

    end % if loop mvpa option

end  % for loop over contrasts or conditions


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING GLM RESULTS');
fprintf('\n\n');

if ~dorobfit_parcelwise
        savefilenamedata = fullfile(resultsdir, ['regression_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
        save(savefilenamedata, 'regression_stats_results', '-v7.3');
        fprintf('\nSaved regression_stats_results for %s\n', mygroupnamefield);

else
        savefilenamedata = fullfile(resultsdir, ['parcelwise_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
        save(savefilenamedata, 'parcelwise_stats_results', '-v7.3');
        fprintf('\nSaved parcelwise_stats_results for %s\n', mygroupnamefield);
end

fprintf('\nFilename: %s\n', savefilenamedata);

fprintf('\n\n');
printhdr('SAVING MVPA RESULTS');
fprintf('\n\n');

savefilenamedata_mvpa = fullfile(resultsdir, ['mvpa_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
save(savefilenamedata_mvpa, 'mvpa_stats_results', '-v7.3');
fprintf('\nSaved mvpa_stats_results for %s\n', mygroupnamefield);

fprintf('\nFilename: %s\n', savefilenamedata_mvpa);

