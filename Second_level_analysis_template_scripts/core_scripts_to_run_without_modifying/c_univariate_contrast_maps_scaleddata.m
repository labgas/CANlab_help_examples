%% Check whether contrasts are defined
% ------------------------------------------------------------------------
if ~isfield(DAT, 'contrasts') || isempty(DAT.contrasts)
    
    % skip
    disp('No contrasts specified. Skipping contrast estimation');

    return
end


%% Specify which montage to add title to. This is fixed for a given slice display
% ------------------------------------------------------------------------
whmontage = 5; 
plugin_check_or_create_slice_display; % script, checks for o2 and uses whmontage

scaling_string = 'scaling_z_score_conditions';


%% Masking
% ------------------------------------------------------------------------
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
 

%% T-test on each contrast image
% ------------------------------------------------------------------------
docompact2 = 0;  % 0 for default, 1 for compact2 version

printhdr('Contrast maps - OLS t-tests');

k = size(DAT.contrasts, 1);
contrast_t_fdr = {};

if docompact2
    o2 = canlab_results_fmridisplay([], 'compact2', 'noverbose');
    whmontage = 1;
else
    create_figure('fmridisplay'); axis off
    o2 = canlab_results_fmridisplay([], 'noverbose');
    whmontage = 5;
end

for i = 1:k
    
    figtitle = sprintf('%s_%s_05_FDR_%s', DAT.contrastnames{i}, scaling_string, mask_string);
    figstr = format_strings_for_legend(figtitle); 
    figstr = figstr{1};
    
    disp(' ')
    printhdr(figstr);

    contrast_t_fdr{i} = ttest( DATA_OBJ_CONsc{i}, .05, 'fdr');
    
    if apply_mask_before_fdr
        contrast_t_fdr{i} = apply_mask(contrast_t_fdr{i}, mask);
        contrast_t_fdr{i} = threshold(contrast_t_fdr{i}, .05, 'fdr');
    end
    
    fprintf('\nShowing results at FDR q < .05, %s\n', mask_string);    
    
    % 1st plot at 0.05 FDR
    % -----------------------------------------------
    o2 = removeblobs(o2);
    r = region(contrast_t_fdr{i}, 'noverbose');
    o2 = addblobs(o2, r, 'splitcolor', 'noverbose'); %  To change colors: 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 
    o2 = legend(o2);
    
    %axes(o2.montage{whmontage}.axis_handles(5));
    [o2, title_handle] = title_montage(o2, whmontage, figstr);
    set(title_handle, 'FontSize', 18);

    % Activate, name, and save figure
    fighan = activate_figures(o2); % Find and activate figure associated with existing fmridisplay object o2
    if ~isempty(fighan)
        set(fighan{1}, 'Tag', figtitle);
        plugin_save_figure;             % re-activate the slice montage figure for saving, name, and save
    else
        disp('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
    end
    
    % Table of results (5 vox or greater)
    % @lukasvo76: changed the default from 3 to 5 voxels
    fprintf('Table of results for clusters >= 5 contiguous voxels.');
    r(cat(1, r.numVox) < 5) = [];                   % r = extent_threshold(r);
    [rpos, rneg] = table(r);       % add labels
    r = [rpos rneg];               % re-concatenate labeled regions
    
    % Montage of regions in table (plot and save)
    if ~isempty(r)
        o3 = montage(r, 'colormap', 'regioncenters');
        
        % Activate, name, and save figure - then close
        figtitle = sprintf('%s_%s_05_FDR_regions_%s', DAT.contrastnames{i}, scaling_string, mask_string);
        region_fig_han = activate_figures(o3);
        if ~isempty(region_fig_han)
            set(region_fig_han{1}, 'Tag', figtitle);
            plugin_save_figure;
            close(region_fig_han{1}), clear o3
        else
            disp('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
        end
        
    end % end conditional montage plot if there are regions to show
    
    
    % 2nd plot at 0.01 uncorrected
    % -----------------------------------------------
    figtitle = sprintf('%s_%s_01_unc_%s', DAT.contrastnames{i}, scaling_string, mask_string);
    figstr = format_strings_for_legend(figtitle);
    figstr = figstr{1};
    
    disp(' ')
    printhdr(figstr);
    
    o2 = removeblobs(o2);
    contrast_t_unc{i} = threshold(contrast_t_fdr{i}, .01, 'unc');
    
    r = region(contrast_t_unc{i}, 'noverbose');
    o2 = addblobs(o2, r, 'splitcolor', {[0 0 1] [0 1 1] [1 .5 0] [1 1 0]}, 'noverbose');
    o2 = legend(o2);

    [o2, title_handle] = title_montage(o2, whmontage, figstr);
    set(title_handle, 'FontSize', 18);

    % Activate, name, and save figure
    fighan = activate_figures(o2); % Find and activate figure associated with existing fmridisplay object o2
    if ~isempty(fighan)
        set(fighan{1}, 'Tag', figtitle);
        plugin_save_figure;             % re-activate the slice montage figure for saving, name, and save
    else
        disp('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
    end
    
    % Table of results (10 vox or greater)
    fprintf('Table of results for clusters >= 10 contiguous voxels.');
    r(cat(1, r.numVox) < 10) = [];    
    [rpos, rneg] = table(r);       % add labels
    r = [rpos rneg];               % re-concatenate labeled regions
    
    
    % Montage of regions in table (plot and save)
    if ~isempty(r)
        o3 = montage(r, 'colormap', 'regioncenters');
        
        % Activate, name, and save figure - then close
        figtitle = sprintf('%s_%s_01_unc_regions_%s', DAT.contrastnames{i}, scaling_string, mask_string);
        region_fig_han = activate_figures(o3);
        if ~isempty(region_fig_han)
            set(region_fig_han{1}, 'Tag', figtitle);
            plugin_save_figure;
            close(region_fig_han{1}), clear o3
        else
            disp('Cannot find figure - Tag field was not set or figure was closed. Skipping save operation.');
        end
        
    end % end conditional montage plot if there are regions to show
     
end