%%% prep_3a_run_second_level_regression_and_save.m

% USAGE
%
% This script 
% 1) runs second⁻level (i.e. across subjects) regression analyses
% for each within-subject CONTRAST or CONDITION registered in the DAT
% structure, either
%   a) voxel-wise, calling CANlab's regress() function under the hood,
%   including its robust regression option if specified
%   b) parcel-wise, calling CANlab's robfit_parcelwise() function under the
%   hood, which is robust by default
% 2) saves the results using standard naming and location
% 
% - To specify analysis options, run a2_set_default_options
% - To choose between conditions and contrasts, set option below
% - To get results reports, run c2a_second_level_regression
%
% OPTIONS SPECIFIED IN a2_set_default_options
%
% - dorobust : robust regression or OLS (true/false)
% - dorobfit_parcelwise: voxel- or parcelwise regression (true/false) - % OPTION ADDED BY @LUKASVO76 MAY 2022
% - myscaling: 'raw', 'scaled', or 'scaled_contrasts' (defined in a2_set_..., image scaling done in prep_2_... and prep_3_... data load)
% - design_matrix_type: 'group', 'custom', or 'onesample'
%                       Group: use DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.contrasts{c}.group;
%                       Custom: use all columns of table object DAT.BETWEENPERSON.contrasts{c};
%                       Onesample: use constant (i.e. intercept) only
%
% 'group' option 
% Assuming that groups are concatenated in contrast image lists, and
% regressor values of 1 or -1 will specify the group identity for each
% image. Requires DAT.BETWEENPERSON.group field specifying group membership for
% each image.
%
% 'custom' option: 
% Can enter a multi-column design matrix for each contrast
% Design matrix can be different for each contrast
%
% 'onesample' option:
% Only adds intercept, hence performs a one-sample t-test on contrast
% images across all subjects, similarly to c_univariate_contrast_maps_
% scripts, but with more flexible options including scaling and robustfit
% OPTION ADDED BY @LUKASVO76 MAY 2022
%
% To set up group and custom variables, see prep_1b_prep_behavioral_data
%
% OPTIONS TO BE SPECIFIED IN THIS SCRIPT
%
% -mygroupfieldname: 'contrasts' or 'conditions'
% -results_suffix: name to add to results file to specify model in case of multiple models, e.g. 'covariate_rating'
%
%__________________________________________________________________________
%
% revamped by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% prep_3a_run_second_level_regression_and_save.m         v2.3
% last modified: 2022/07/28


%% SETTINGS
%--------------------------------------------------------------------------

% options to be specified here

mygroupnamefield = 'contrasts'; 
results_suffix = ''; % do not delete, leave empty if not needed

% options set in a2_set_default_options

options_needed = {'dorobust', 'myscaling_glm', 'design_matrix_type'};
options_exist = cellfun(@exist, options_needed); 

option_default_values = {true, 'raw', 'onesample'}; % defaults if we cannot find info in a2_set_default_options at all ; @lukasvo76: changed the defaults to align with a2_set_default_options

plugin_get_options_for_analysis_script


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
%--------------------------------------------------------------------------

if exist('maskname_glm', 'var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')
    [~,maskname_short] = fileparts(maskname_glm);
    mask_string = sprintf('within mask %s', maskname_short);
    glmmask = fmri_mask_image(maskname_glm, 'noverbose'); 
else
    mask_string = sprintf('without masking');
end  


%% RUN SECOND LEVEL REGRESSION FOR EACH CONTRAST
% -------------------------------------------------------------------------

switch mygroupnamefield
    
    case 'contrasts'

        kc = size(DAT.contrasts, 1);
        
        printhdr('running second-level regressions on first-level contrasts');
        
    case 'conditions'
        
        kc = size(DAT.conditions, 2);
        
        printhdr('running second-level regressions on first-level conditions');
        
    otherwise
        
        error('\ninvalid option "%s" defined in mygroupnamefield variable, choose between "contrasts" and "conditions"\n',mygroupnamefield)

end

if ~dorobfit_parcelwise
    regression_stats_results = cell(1, kc);
else
    parcelwise_stats_results = cell(1,kc);
end

for c = 1:kc
    
    % GET DESIGN MATRIX FOR THIS CONTRAST OR CONDITION
    % ---------------------------------------------------------------------
    
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
            
        case 'group'
            
            % Use 'groups' single regressor
            if ~isempty(DAT.BETWEENPERSON.group)
                group = DAT.BETWEENPERSON.group;
                groupnames = {'group'};
                X = group;
                imgs_nan = [];
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
            
            error('\ninvalid option "%s" defined in design_matrix_type variable, choose between "group", "custom", or "onesample"\n', design_matrix_type);
            
    end

    switch mygroupnamefield
        
        case 'contrasts'
            printstr(DAT.contrastnames{c});
            printstr(dashes)
        case 'conditions'
            printstr(DAT.conditions{c});
            printstr(dashes)
    
    end

    
    % SELECT DATA FOR THIS CONTRAST/CONDITION
    % ---------------------------------------------------------------------
    switch mygroupnamefield
        
        case 'contrasts'
            
            switch myscaling_glm

                case 'raw'
                    fprintf('\ncontrast calculated on raw (unscaled) condition images used in second-level GLM\n\n');
                    scaling_string = 'no_scaling';
                    cat_obj = DATA_OBJ_CON{c};
                    if imgs_nan
                        cat_obj = cat_obj.get_wh_image(imgs_nan);
                    end

                case 'scaled'
                    fprintf('\ncontrast calculated on z-scored condition images used in second-level GLM\n\n');
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
                    error('\ninvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw", "scaled", or "scaled_constrast" given option "%s" defined in mygroupnamefield variable\n', myscaling_glm, mygroupnamefield);

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
                    error('\ninvalid combination of option "%s" defined in myscaling_glm_variable in a2_set_default_options script and option "%s" defined in mygroupnamefield variable, choose between "raw" and "scaled" options\n',myscaling_glm,mygroupnamefield);

                otherwise
                    error('\ninvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw",  and "scaled", given option "%s" defined in mygroupnamefield variable\n', myscaling_glm, mygroupnamefield);

            end
            
    end % switch mygroupnamefield - contrasts or conditions
    
  
    % FORMAT AND ATTACH DESIGN MATRIX
    % ---------------------------------------------------------------------
    
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
            warning('Some columns are not mean-centered or effects coded. \nIntercept may not be interpretable');
            fprintf('\nColumns: ')
            fprintf('%d \n', find(~meancentered & ~effectscoded));
        else
            fprintf('\nChecked OK: All columns mean-centered or are effects-coded [1 -1 0]\n\n');
        end

        if any(vifs > 2)
            fprintf('\n');
            warning('Some regressors have high variance inflation factors: Parameters might be poorly estimated or uninterpretable.');
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
    
        if all(cat_obj.X > 0) || all(cat_obj.X < 0)
            % Only positive or negative weights - nothing to compare

            fprintf('\n');
            warning('Only positive or negative regressor values - bad design, please check');
            fprintf('\n');
            
            continue
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

        if ~strcmpi(design_matrix_type,'onesample')
            % out.t has t maps for all regressors, intercept is last
            switch mygroupnamefield
                case 'contrasts'
                    regression_stats = regress(cat_obj, .05, 'unc', robuststring, 'analysis_name', DAT.contrastnames{c}, 'variable_names', groupnames, 'nodisplay');
                case 'conditions'
                    regression_stats = regress(cat_obj, .05, 'unc', robuststring, 'analysis_name', DAT.conditions{c}, 'variable_names', groupnames, 'nodisplay');
            end
        else
            % out.t has t maps for intercept only
            switch mygroupnamefield
                case 'contrasts'
                    regression_stats = regress(cat_obj, .05, 'unc', robuststring, 'analysis_name', DAT.contrastnames{c}, 'variable_names', groupnames, 'nointercept', 'nodisplay');
                case 'conditions'
                    regression_stats = regress(cat_obj, .05, 'unc', robuststring, 'analysis_name', DAT.conditions{c}, 'variable_names', groupnames, 'nointercept', 'nodisplay');
            end
        end

        % Make sure variable types are right data formats
        regression_stats.design_table = design_table;
        regression_stats.t = enforce_variable_types(regression_stats.t);
        regression_stats.b = enforce_variable_types(regression_stats.b);
        regression_stats.df = enforce_variable_types(regression_stats.df);
        regression_stats.sigma = enforce_variable_types(regression_stats.sigma);

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
        
        fprintf ('\nShowing results at p uncor < 0.05: %s\nEffect: %s\n\n', regression_stats.analysis_name, mask_string);
        
        t = regression_stats.t;
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

        % KEEP RESULTS OBJECTS IN CELL ARRAY FOR SAVING
        % ---------------------------------------------------------------------

        regression_stats_results{c} = regression_stats;
        
        if exist(maskname_glm,'file')
            regression_stats_results{c}.maskname = maskname_glm;
        end
        
        if dorobust
            fprintf('\nCumulative run time:'), toc(regresstime); 
        end
        
    % PARCEL-WISE    
        
    else
        
        if csf_wm_covs && remove_outliers
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',true,'doplot',false);
        elseif csf_wm_covs && ~remove_outliers
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',false,'doplot',false);
        elseif ~csf_wm_covs && remove_outliers
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',false,'remove_outliers',true,'doplot',false);
        else
            parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'doplot',false);
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
        num_effects = size(parcelwise_stats.t_obj.dat, 2); % number of regressors
        o2 = canlab_results_fmridisplay([], 'multirow', num_effects, 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        for j = 1:num_effects

            fprintf ('\nShowing GLM results at p uncor < 0.05: %s\nEffect: %s, %s\n\n', parcelwise_stats.analysis_name, parcelwise_stats.variable_names{j}, mask_string);

            tj = get_wh_image(parcelwise_stats.t_obj, j);
                if maskname_short
                    tj = apply_mask(tj, glmmask);
                end
            tj = threshold(tj, .05, 'unc'); 

            o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
            o2 = title_montage(o2, 2*j, [parcelwise_stats.analysis_name ' ' parcelwise_stats.variable_names{j}]);

        end

        figtitle = sprintf('%s_05_unc_montage_%s_%s', parcelwise_stats.analysis_name, scaling_string, mask_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;
            if save_figures
                plugin_save_figure;
            end
        clear o2, clear figtitle, clear j, clear tj

        % KEEP RESULTS OBJECTS IN CELL ARRAY FOR SAVING
        % ---------------------------------------------------------------------

        parcelwise_stats_results{c} = parcelwise_stats;
        
        if exist(maskname_glm,'file')
            parcelwise_stats_results{c}.maskname = maskname_glm;
        end
        
    end % if loop voxel- versus parcelwise
    
end  % for loop over contrasts or conditions


%% SAVE RESULTS
% -------------------------------------------------------------------------
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

