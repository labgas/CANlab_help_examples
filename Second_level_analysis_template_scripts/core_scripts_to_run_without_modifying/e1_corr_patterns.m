%% ery_4a_secondlevel_m6m_s8_e1_corr_patterns_conds.m
% 
% 
% USAGE
% 
% This script displays calculates, thresholds, and plots correlation maps
% between all pairwise combinations of condition or contrast images using
% the searchlight approach, calling the CANlab function
% searchlight_correlation() under the hood
% 
% Run this script with Matlab's publish function to generate html report of results:
% publish('e1_corr_patterns_conds','outputDir',htmlsavedir)
% 
% 
% MANDATORY OPTIONS
% 
% * mygroupnamefield = 'conditions'/'contrasts';      calculate correlations between conditions or contrasts defined in DAT  
% * results_suffix = '';                              adds a suffix of your choice to .mat file with results that will be saved
% 
% NOTE: do NOT delete the results_suffix option, leave empty if not needed
% NOTE: do NOT use to add a suffix specifying the scaling or masking option, this will be added automatically
% 
% CUSTOM OPTIONS FOR THIS SCRIPT SET IN A2 SCRIPT 
% 
% * r_threshold_corr = 0.xx;                          r threshold for correlation coeffs 
% * corr_type = 'Pearson'/'Spearman'/'Kendall';
% 
% COPY OPTIONS FROM PREP_3a_ SCRIPT IF DIFFERENT FROM DEFAULTS SET IN A2 SCRIPT
% 
% * atlasname_glm = 'atlas_name';
% * maskname_glm = 'mask_name';
% * myscaling_glm = 'raw/scaled';
% 
% -------------------------------------------------------------------------
% 
% author: Lukas Van Oudenhove
% date:   KU Leuven, February, 2023
% 
% -------------------------------------------------------------------------
%
% e1_corr_patterns_conds.m            v1.0
% last modified: 2023/02/15
% 
% 
%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6m_s0_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% SET MANDATORY OPTIONS

mygroupnamefield = 'conditions'; 
results_suffix = '';                            % adds a suffix of your choice to .mat file with results that will be saved

% CUSTOM OPTIONS FOR THIS SCRIPT SET IN A2 SCRIPT 

r_threshold_corr = 0.60;                      % r threshold for correlation coeffs 
%  corr_type = 'Pearson'/'Spearman'/'Kendall';

% NOTE: do NOT delete the results_suffix option, leave empty if not needed
% NOTE: do NOT use to add a suffix specifying the scaling or masking option, this will be added automatically

% COPY OPTIONS FROM PREP_3a_ SCRIPT IF DIFFERENT FROM DEFAULTS SET IN A2 SCRIPT

% atlasname_glm = 'atlas_name';
% maskname_glm = 'mask_name';
% myscaling_glm = 'raw/scaled';


%% GET MASKING OPTION
% -------------------------------------------------------------------------

if exist('atlasname_glm','var') && ~isempty(atlasname_glm) && exist(atlasname_glm, 'file')

    if contains(atlasname_glm,'.mat')
        
        load(atlasname_glm);
        clear mask
        
        [~,maskname_short] = fileparts(atlasname_glm);
        mask_string = sprintf('masked with %s', maskname_short);
        mask = fmri_mask_image(combined_atlas,'noverbose');
        
        fprintf('\nShowing correlations in custom-made atlas/mask %s\n\n', maskname_short);

    elseif ischar(atlasname_glm)

        combined_atlas = load_atlas(atlasname_glm);
        
        [~,maskname_short] = fileparts(atlasname_glm);
        mask_string = sprintf('in atlas %s',atlasname_glm);
        mask = fmri_mask_image(combined_atlas,'noverbose');
        
        fprintf('\nShowing correlations in atlas/mask %s\n\n', maskname_short);

    else

         error('\ninvalid option "%s" defined in atlasname_glm variable, should be a keyword for load_atlas.m or a .mat file containing an atlas object, check docs"\n\n',atlasname_glm)

    end

else

    if exist('maskname_glm', 'var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')

        [~,maskname_short] = fileparts(maskname_glm);
        mask_string = sprintf('masked with %s', maskname_short);
        mask = fmri_mask_image(maskname_glm, 'noverbose'); 
        
        fprintf('\nShowing correlations in mask %s\n\n', maskname_short);

    else

        mask_string = sprintf('without masking');
        
        fprintf('\nShowing correlations without masking\n\n');

    end  

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


%% PREP WORK
% -------------------------------------------------------------------------

switch mygroupnamefield
    
    case 'contrasts'

        kc = size(DAT.contrasts, 1);
        con_names = DAT.contrastnames;
       
        fprintf('\nCALCULATING CORRELATION MAPS ON FIRST LEVEL CONTRASTS\n\n');
        
    case 'conditions'
        
        kc = size(DAT.conditions, 2);
        con_names = DAT.conditions;
       
        fprintf('\nCALCULATING CORRELATION MAPS ON FIRST LEVEL CONDITIONS\n\n');
        
    otherwise
        
        error('\ninvalid option "%s" defined in mygroupnamefield variable, choose between "contrasts" and "conditions"\n\n',mygroupnamefield)

end


fprintf('\n\n');
printhdr('SCALING DATA IF REQUESTED IN OPTIONS');
fprintf('\n\n');

cat_obj = cell(1,kc);

for c = 1:kc
    
    switch mygroupnamefield
        
        case 'contrasts'
            
            switch myscaling_glm

                case 'raw'
                    if c ==1
                        fprintf('\nContrast calculated on raw (unscaled) condition images\n\n');
                        scaling_string = 'no_scaling';
                    end
                    cat_obj{c} = DATA_OBJ_CON{c};

                case 'scaled'
                    if c == 1
                        fprintf('\nContrast calculated on z-scored condition images\n\n');
                        scaling_string = 'scaling_z_score_conditions';
                    end
                    cat_obj{c} = DATA_OBJ_CONsc{c};

                case 'scaled_contrasts'
                    if c == 1
                        fprintf('\nl2norm scaled contrast images used\n\n');
                        scaling_string = 'scaling_l2norm_contrasts';
                    end
                    cat_obj{c} = DATA_OBJ_CONscc{c};

                otherwise
                    error('\nInvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw", "scaled", or "scaled_constrast" given option "%s" defined in mygroupnamefield variable\n\n', myscaling_glm, mygroupnamefield);

            end
            
        case 'conditions'
            
            switch myscaling_glm

                case 'raw'
                    if c == 1
                        fprintf('\nRaw (unscaled) condition images used\n\n');
                        scaling_string = 'no_scaling';
                    end
                    cat_obj{c} = DATA_OBJ{c};

                case 'scaled'
                    if c == 1
                        fprintf('\nZ-scored condition images used\n\n');
                        scaling_string = 'scaling_z_score_conditions';
                    end
                    cat_obj{c} = DATA_OBJsc{c};

                case 'scaled_contrasts'
                    error('\nInvalid combination of option "%s" defined in myscaling_glm_variable in a2_set_default_options script and option "%s" defined in mygroupnamefield variable, choose between "raw" and "scaled" options\n\n',myscaling_glm,mygroupnamefield);

                otherwise
                    error('\nInvalid option "%s" defined in myscaling_glm variable in a2_set_default_options script, choose between "raw",  and "scaled", given option "%s" defined in mygroupnamefield variable\n\n', myscaling_glm, mygroupnamefield);

            end
            
    end % switch mygroupnamefield - contrasts or conditions
    
    if c == 1 % we only need to do this once
        mask = resample_space(mask,cat_obj{c});
    end
    
    cat_obj{c} = apply_mask(cat_obj{c},mask);
    
end


%% CALCULATE CORRELATION MAPS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('CALCULATING CORRELATION MAPS');
fprintf('\n\n');

j = 1;
idx = 0;

corr_maps = cell(1,size(combntns([1:kc],2),1));
sl_size = cell(1,size(combntns([1:kc],2),1));
con_title = cell(1,size(combntns([1:kc],2),1));


while j < kc
    
    idx = idx +1;

    dat1 = cat_obj{j};
    dat2 = cat_obj{j+1};

    [~, corr_maps{idx}, sl_size{idx}] = searchlight_correlation(dat1,dat2,'type', corr_type);       
    con_title{idx} = [con_names{j} ' ' con_names{j+1}];
            
    
        for k = 2:(kc - j)
            
            idx = idx +1;
            
            dat1 = cat_obj{j};
            dat2 = cat_obj{j+k};

            [~, corr_maps{idx}, sl_size{idx}] = searchlight_correlation(dat1,dat2);
            con_title{idx} = [con_names{j} ' ' con_names{j+k}];
            
        end
        
    j = j+1;
    
end


%% PLOT CORRELATION MAPS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('PLOTTING CORRELATION MAPS');
fprintf('\n\n');

o2 = canlab_results_fmridisplay([], 'multirow', idx, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

n_voxels_sig = cell(1,idx);

for l = 1:idx
    
    corr_maps{l} = threshold(corr_maps{l},[r_threshold_corr 1],'raw-between');
    n_voxels_sig{l} = sum(corr_maps{l}.sig);
    o2 = addblobs(o2, region(corr_maps{l}), 'wh_montages', (2*l)-1:2*l, 'colormap', 'maxcolor', [0.94 0.98 0.13], 'mincolor', [0.47 0.11 0.43], 'cmaprange', [r_threshold_corr 1]);
    o2 = title_montage(o2, 2*l, [con_title{l} ' r > ' num2str(r_threshold_corr) ' ' mask_string ' ' scaling_string]);
    
end

h = findobj('Type','figure');

for f = 1:size(h,1)
    set(h(f),'WindowState','maximized');
end


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING CORRELATION RESULTS');
fprintf('\n\n');

savefilenamedata = fullfile(resultsdir, ['correlation_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
save(savefilenamedata, 'corr_maps', 'sl_size', 'n_voxels_sig', '-v7.3');
fprintf('\nSaved correlation_maps for %s\n', mygroupnamefield);

fprintf('\nFilename: %s\n', savefilenamedata);

