% THIS SCRIPT DISPLAYS BETWEEN-PERSON CONTRASTS as 
% 1) The effect of DAT.BETWEENPERSON.group if design_matrix_type in a2_
% script == 'group'
%      - effect of the regressor with group membership (1,-1)
% OR the effect of DAT.BETWEENPERSON.contrasts{i} if design_matrix_type in a2_
% script == 'custom'
%      - effect of the regressor with group membership (1,-1) as above
%      - effect of any other covariates entered in
%      DAT.BETWEENPERSON.contrasts{i} in prep_1b script
% 2) Average activation controlling for the covariates
%      = 'Intercept' 
%
% Requires DAT.BETWEENPERSON.group field specifying group membership for
% each image, or DAT.BETWEENPERSON.contrasts{i}.xxx fields specifying group
% membership and potentially other covariates for each image
%
% Marta Ceko 2017 (based on existing scripts) 
% adapted by @lukasvo76 on LaBGAS fork 2021

%% LOAD REGRESSION RESULTS (if needed)

switch myscaling_glm
    case 'raw'
        printstr('Raw (unscaled) images used in between-person GLM');
        scaling_string = 'no_scaling';
        
    case 'scaled'
        printstr('Z-scored images used in between-person GLM');
        scaling_string = 'scaling_z_score_conditions';
        
    case 'scaled_contrasts'
        printstr('l2norm scaled contrast images used in between-person GLM');
        scaling_string = 'scaling_l2norm_contrasts';
        
    otherwise
        error('myscaling_glm must be ''raw'' or ''scaled'' or ''scaled_contrasts''');
        
end       
        
savefilenamedata = fullfile(resultsdir, ['regression_stats_and_maps_',scaling_string,'.mat']);

if ~exist(savefilenamedata, 'file')
    disp('Run prep_3a_run_second_level_regression_and_save.m to get regression results.'); 
    disp('No saved results file.  Skipping this analysis.')
    return
end

fprintf('\nLoading regression results and maps from %s\n\n', savefilenamedata);
load(savefilenamedata, 'regression_stats_results');


%% MASKING

if exist(maskname_glm, 'file')
    apply_mask_before_fdr = true;
    maskname_short = split(maskname_glm,"\");
    maskname_short = maskname_short{size(maskname_short,1)};
    maskname_short = maskname_short(:,1:end-4);
    mask_string = sprintf('within mask %s', maskname_short);
    mask = fmri_data_st(maskname_glm, 'noverbose'); 
else
    apply_mask_before_fdr = false;
    mask_string = sprintf('without masking');
end  


%% UNIVARIATE CONTRASTS WHOLE BRAIN

ncontrasts = size (regression_stats_results, 2);

for c = 1:ncontrasts

    analysisname = regression_stats_results{c}.analysis_name;
    names = regression_stats_results{c}.variable_names;
    t = regression_stats_results{c}.t;
    
    printhdr(analysisname)
    disp('Regressors: ')
    disp(names)
    
    if isfield(regression_stats_results{c}, 'design_table')
        disp(regression_stats_results{c}.design_table);
    end
    
    num_effects = size(t.dat, 2); % number of regressors
    
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: montage at 0.05 FDR corrected
    % ---------------------------------------------------------------------
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects);
    
    for j = 1:num_effects
        
        fprintf ('\nShowing results at FDR q < 0.05: %s\nEffect: %s, %s\n\n', analysisname, names{j}, mask_string);
        
        tj = get_wh_image(t, j);
            if apply_mask_before_fdr
                tj = apply_mask(tj, mask);
            end
        tj = threshold(tj, .05, 'fdr'); 
        
        o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
        o2 = title_montage(o2, 2*j, [analysisname ' ' names{j}]);
    end
    
    figtitle = sprintf('Regression results %s_05_FDR_%s_%s', analysisname, scaling_string, mask_string);
    set(gcf, 'Tag', figtitle);
    plugin_save_figure;
    clear o2, clear figtitle
    
    for j = 1:num_effects
        
        fprintf('\n\nTable of results for clusters >= 5 contiguous voxels.');
        
        tj = get_wh_image(t, j);
            if apply_mask_before_fdr
                tj = apply_mask(tj, mask);
            end
        tj = threshold(tj, .05, 'fdr'); 
        
        r = region(tj, 'noverbose');
        r(cat(1, r.numVox) < 5) = [];                   % r = extent_threshold(r);
        [rpos, rneg] = table(r);       % add labels
        r = [rpos rneg];               % re-concatenate labeled regions

        % Montage of regions in table (plot and save)
        if ~isempty(r)
            o3 = montage(r, 'colormap', 'regioncenters');

            % Activate, name, and save figure - then close
            figtitle = sprintf('%s_05_FDR_regions_%s_%s_%s', analysisname, names{j}, scaling_string, mask_string);
            region_fig_han = activate_figures(o3);
            if ~isempty(region_fig_han)
                set(region_fig_han{1}, 'Tag', figtitle);
                plugin_save_figure;
                close(region_fig_han{1}), clear o3, clear figtitle
            else
                disp('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
            end

        end % end conditional montage plot if there are regions to show
    end

    
    % BETWEEN-SUBJECT REGRESSORS & INTERCEPT: montage at 0.01 uncorrected
    % ---------------------------------------------------------------------    
    o2 = canlab_results_fmridisplay([], 'multirow', num_effects);
    
    for j = 1:num_effects
        
        fprintf ('\nShowing results at uncorrected p < 0.01: %s\nEffect: %s, %s\n\n', analysisname, names{j}, mask_string);
        
        tj = get_wh_image(t, j);
            if apply_mask_before_fdr
                tj = apply_mask(tj, mask);
            end
        tj = threshold(tj, .01, 'unc'); 
        
        o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j);
        o2 = title_montage(o2, 2*j, [analysisname ' ' names{j}]);
    end
    
    figtitle = sprintf('Regression results %s_01_unc %s_%s', analysisname, scaling_string, mask_string);
    set(gcf, 'Tag', figtitle);
    plugin_save_figure;
    clear o2, clear figtitle
        
    for j = 1:num_effects
        
        fprintf('\n\nTable of results for clusters >= 10 contiguous voxels.');
        
        tj = get_wh_image(t, j);
            if apply_mask_before_fdr
                tj = apply_mask(tj, mask);
            end
        tj = threshold(tj, .01, 'unc'); 
        
        r = region(tj, 'noverbose');
        r(cat(1, r.numVox) < 10) = [];                   % r = extent_threshold(r);
        [rpos, rneg] = table(r);       % add labels
        r = [rpos rneg];               % re-concatenate labeled regions

        % Montage of regions in table (plot and save)
        if ~isempty(r)
            o3 = montage(r, 'colormap', 'regioncenters');

            % Activate, name, and save figure - then close
            figtitle = sprintf('%s_01_unc_regions_%s_%s_%s', analysisname, names{j}, scaling_string, mask_string);
            region_fig_han = activate_figures(o3);
            if ~isempty(region_fig_han)
                set(region_fig_han{1}, 'Tag', figtitle);
                plugin_save_figure;
                close(region_fig_han{1}), clear o3, clear figtitle
            else
                disp('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
            end

        end % end conditional montage plot if there are regions to show
    end
   
end % c contrasts

