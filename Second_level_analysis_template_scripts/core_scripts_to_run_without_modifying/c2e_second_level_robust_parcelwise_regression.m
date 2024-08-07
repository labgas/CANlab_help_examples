% NOTE: the use of this script is being deprecated, as the option to run
% robust parcelwise second-level regression on conditions is now built into
% prep_3a_run_second_level_regression_and_save.m; please use that script
% instead of this one
% @lukasvo76, May 2022

% THIS SCRIPT RUNS BETWEEN-PERSON (2nd-level) Robust Parcelwise Regression 
% analyses for each within-person CONTRAST registered in the analysis, 
% using Tor's new robfit_parcelwise method for fmri_data objects
% 
% To specify analysis options, run a2_set_default_options
%
% Analysis options include:
% - myscaling: 'raw' or 'scaled' (image scaling done in prep_2_... data load)
% - design_matrix_type: 'group' or 'custom'
%                       Group: use DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.contrasts{c}.group;
%                       Custom: use all columns of table object DAT.BETWEENPERSON.contrasts{c};
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
% To set up group and custom variables, see prep_1b_prep_behavioral_data

% --------------------------------------------------------------------

% USER OPTIONS
% This is a standard block of code that can be used in multiple scripts.
% Each script will have its own options needed and default values for
% these.
% The code: 
% (1) Checks whether the option variables exist
% (2) Runs a2_set_default_options if any are missing
% (3) Checks again and uses the default options if they are still missing
% (e.g., not specified in an older/incomplete copy of a2_set_default_options)


%% Settings
%--------------------------------------------------------------------------

% Now set in a2_set_default_options
options_needed = {'myscaling_glm', 'design_matrix_type', 'csf_wm_covs', 'remove_outliers'}; % options we are looking for. Set in a2_set_default_options % @lukasvo76: or in study-specific version of that script
options_exist = cellfun(@exist, options_needed); 

option_default_values = {'raw', 'group', false, false}; % defaults if we cannot find info in a2_set_default_options at all ; @lukasvo76: changed the defaults to align with a2_emosymp_m1_s1_set_default_options

plugin_get_options_for_analysis_script


%% Check for required DAT fields
% -------------------------------------------------------------------------

% List required fields in DAT, in cell array:
required_fields = {'BETWEENPERSON', 'contrastnames', 'contrasts' 'contrastcolors'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end


%% Run between-person robust parcelwise second-level regression for each contrast
% -------------------------------------------------------------------------------

printhdr('Second-level robust parcelwise regressions discriminate between-person contrasts');

kc = size(DAT.contrasts, 1);

robfit_parcel_stats_results = cell(1, kc);

for c = 1:kc
    
    % Get design matrix for this contrast
    % ---------------------------------------------------------------------
    
    mygroupnamefield = 'contrasts';  % 'conditions' or 'contrasts'
    
    switch design_matrix_type
        case 'custom'
            
            % Define design matrix X "design_matrix"
            % Use custom matrix for each condition/contrast
            table_obj = DAT.BETWEENPERSON.(mygroupnamefield){c};
            groupnames = table_obj.Properties.VariableNames;
            X = table2array(table_obj);
            idx_nan = ~isnan(X);
            idx_nan = ~(sum(idx_nan,2) < size(idx_nan,2)); % at least one column of X contains NaN
            imgs_nan = [1:size(X,1)];
            imgs_nan = imgs_nan(idx_nan');
            X = X(idx_nan,:);
            
        case 'group'
            
            % Use 'groups' single regressor
            group = DAT.BETWEENPERSON.(mygroupnamefield){c}.group;
            groupnames = {'group'};
            X = group;
            imgs_nan = [];
            
            if isempty(group)
                fprintf('Group not defined for contrast %s. Skipping.\n', DAT.contrastnames{c});
                continue
            end
            
        otherwise error('Incorrect option specified for design_matrix_type');
    end
    
    printstr(DAT.contrastnames{c});
    printstr(dashes)
    
    mycontrast = DAT.contrasts(c, :);
    wh = find(mycontrast);
    
    % Select data for this contrast
    % ---------------------------------------------------------------------
    
    switch myscaling_glm
        case 'raw'
            printstr('Raw (unscaled) images used in between-person robust parcelwise GLM');
            scaling_string = 'no_scaling';
            cat_obj = DATA_OBJ_CON{c};
            if imgs_nan
                cat_obj = cat_obj.get_wh_image(imgs_nan);
            end
            
        case 'scaled'
            printstr('Z-scored images used in robust parcelwise between-person GLM');
            scaling_string = 'scaling_z_score_conditions';
            cat_obj = DATA_OBJ_CONsc{c};
            if imgs_nan
                cat_obj = cat_obj.get_wh_image(imgs_nan);
            end
            
        case 'scaled_contrasts'
            printstr('l2norm scaled contrast images used in robust parcelwise between-person GLM');
            scaling_string = 'scaling_l2norm_contrasts';
            cat_obj = DATA_OBJ_CONscc{c};
            if imgs_nan
                cat_obj = cat_obj.get_wh_image(imgs_nan);
            end
            
        otherwise
            error('myscaling must be ''raw'' or ''scaled'' or ''scaled_contrasts''');
    end
    
  
    
    % Format and attach regressors
    % ---------------------------------------------------------------------
    
    % Confirm design_matrix is 1, -1, or mean-centered
    meancentered = ~(abs(mean(X)) > 1000 * eps);
    effectscoded = all(X == 1 | X == -1 | X == 0, 1);
    isconstant = all(X == mean(X, 1), 1);
    vifs = getvif(X);
    
    if any(isconstant)
        disp('An intercept appears to be added manually. Do not include an intercept - it will be added automatically.');
        disp('Skipping this contrast.')
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
        disp('Warning: some columns are not mean-centered or effects coded. \nIntercept may not be interpretable.\n');
        fprintf('Columns: ')
        fprintf('%d ', find(~meancentered & ~effectscoded));
        fprintf('\n');
    else
        disp('Checked OK: All columns mean-centered or are effects-coded [1 -1 0]');
    end
    
    if any(vifs > 2)
        disp('Some regressors have high variance inflation factors: Parameters might be poorly estimated or uninterpretable.');
    else
        disp('Checked OK: VIFs for all columns are < 2');
    end

    cat_obj.X = X;
    
    % Skip if necessary
    % ---------------------------------------------------------------------
    
    if all(cat_obj.X > 0) | all(cat_obj.X < 0)
        % Only positive or negative weights - nothing to compare
        
        printhdr(' Only positive or negative regressor values - bad design');
        
        continue
    end
    
    % Run robust parcelwise regression model
    % ---------------------------------------------------------------------
    
    % out.t has t maps for all regressors, intercept is last
    if csf_wm_covs && remove_outliers
        robfit_parcel_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',true);
    elseif csf_wm_covs && ~remove_outliers
        robfit_parcel_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',false);
    elseif ~csf_wm_covs && remove_outliers
        robfit_parcel_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',false,'remove_outliers',true);
    else
        robfit_parcel_stats = robfit_parcelwise(cat_obj,'names', groupnames);
    end

    % add design table
    robfit_parcel_stats.design_table = design_table;

    % add regressor names and other meta-data
    robfit_parcel_stats.contrastname = DAT.contrastnames{c};
    robfit_parcel_stats.contrast = DAT.contrasts(c, :);
    
    % add names for analyses and variables: 
    robfit_parcel_stats.analysis_name = DAT.contrastnames{c};
    robfit_parcel_stats.variable_names = [groupnames {'Intercept'}];
    
    % Save stats objects for results later
    % ---------------------------------------------------------------------

    robfit_parcel_stats_results{c} = robfit_parcel_stats;
    
end  % between-person contrast


%% Save
% -------------------------------------------------------------------------
savefilenamedata = fullfile(resultsdir, ['robfit_parcel_stats_and_maps_',scaling_string,'.mat']);

save(savefilenamedata, 'robfit_parcel_stats_results', '-v7.3');
printhdr('Saved robfit_parcel_stats_results for contrasts');
fprintf('Filename: %s\n', savefilenamedata);


%% Subfunction
% -------------------------------------------------------------------------
% @lukasvo76: this subfunction may be redundant as there is a standalone
% CANlab function with the same name

function [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i)

group = []; groupnames = []; groupcolors = [];

if isfield(DAT, 'BETWEENPERSON') && ...
        isfield(DAT.BETWEENPERSON, mygroupnamefield) && ...
        iscell(DAT.BETWEENPERSON.(mygroupnamefield)) && ...
        length(DAT.BETWEENPERSON.(mygroupnamefield)) >= i && ...
        ~isempty(DAT.BETWEENPERSON.(mygroupnamefield){i})
    
    group = DAT.BETWEENPERSON.(mygroupnamefield){i};
    
elseif isfield(DAT, 'BETWEENPERSON') && ...
        isfield(DAT.BETWEENPERSON, 'group') && ...
        ~isempty(DAT.BETWEENPERSON.group)
    
    group = DAT.BETWEENPERSON.group;

end

if isfield(DAT, 'BETWEENPERSON') && isfield(DAT.BETWEENPERSON, 'groupnames')
    groupnames = DAT.BETWEENPERSON.groupnames;
elseif istable(group)
    groupnames = group.Properties.VariableNames(1);
else
    groupnames = {'Group-Pos' 'Group-neg'};
end

if isfield(DAT, 'BETWEENPERSON') && isfield(DAT.BETWEENPERSON, 'groupcolors')
    groupcolors = DAT.BETWEENPERSON.groupcolors;
else
    groupcolors = seaborn_colors(2);
end

if istable(group), group = table2array(group); end

end