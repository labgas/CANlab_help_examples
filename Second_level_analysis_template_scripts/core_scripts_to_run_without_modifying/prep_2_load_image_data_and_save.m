%% prep_2_load_image_data_and_save.m
%
%
% *USAGE*
%
% This prep script 
%
% # loads first-level beta/con images into CANlab's fmri_data_st objects 
% # performs quality control, including plots if requested in a2 script
% # z-scores images and then repeats 1) and 2)
% # saves the relevant resulting objects/variables in a .mat file to resultsdir
% # publishes an html report (run using Matlab's publish function)
%
%
% *OPTIONS*
%
% * dofullplot              default true, can set to false to save time, but not recommended for quality control purposes
%
% * omit_histograms         default false, can set to true to save time, especially in case of large samples but not recommended for quality control purposes
%
% * dozipimages             default false, to avoid load on data upload/download when re-running often, true is useful to save space when running final analyses
%
% * maskname_brain          path to brainmask
%
% * subjs2exclude_data      default empty, subjects to be excluded from data objects, for example because of missing session not allowing all contrasts to be calculated, example {'sub-010' 'sub-018'}
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   Dartmouth, May, 2022
%
% -------------------------------------------------------------------------
%
% prep_2_load_image_data_and_save.m         v2.2
%
% last modified: 2024/12/02
%
%
%% RUN SCRIPT A_SET_UP_PATHS_ALWAYS_RUN_FIRST AND LOAD/CREATE DAT IF NEEDED
% -------------------------------------------------------------------------

a_set_up_paths_always_run_first;

if ~exist('DAT','var')
    
    try
    
        load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
    catch
        
        prep_1_set_conditions_contrasts_colors;
        prep_1b_prep_behavioral_data;
        
    end
    
end


%% SET DEFAULT OPTIONS IF NEEDED
% -------------------------------------------------------------------------

% This is a standard block of code that can be used in multiple scripts.
% Each script will have its own options needed and default values for
% these.
% The code: 
% (1) Checks whether the option variables exist
% (2) Runs a2_set_default_options if any are missing
% (3) Checks again and uses the default options if they are still missing
% (e.g., not specified in an older/incomplete copy of a2_set_default_options)

options_needed = {'dofullplot', 'omit_histograms', 'dozipimages', 'maskname_brain'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed);        % initializing this means a2_set_defaults_options will never run

option_default_values = {true, false, false, which('brain_mask_fmriprep20_template_1000.nii')};          % defaults if we cannot find info in a2_set_default_options at all; @lukasvo76: changed the default for zipping images

plugin_get_options_for_analysis_script


%% PREP AND CHECK IMAGES NAMES
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('PREP WORK');
fprintf('\n\n');

clear imgs cimgs

if ~isempty(subjs2exclude_data) % we have subjects to exclude
    idx_include = ~contains(firstsubjdirs,subjs2exclude_data);
end

for i = 1:size(DAT.conditions,2)
    
    % @lukasvo76: adapted to LaBGAS/BIDS conventional directory structure,
    % if should return 1 since we use wildcards for subject subfolders on
    % Linux OS (see prep_1 script)
    
    if ~isempty(DAT.subfolders) && ~isempty(DAT.subfolders{i})
        
        str = fullfile(datadir, DAT.subfolders{i}, DAT.functional_wildcard{i});
        
%         % Unzip if needed - not needed in LaBGAS case since we typically do
%         not have zipped con images, although that could be implemented
%         % note, Matlab's gunzip() does not remove .gz images, so use eval( ) version.
%         % note, replace spaces with '\ ' 
%         
%         try eval(['!gunzip ' strrep(str, ' ', '\ ') '.gz']), catch, end     % gunzip([str '.gz'])
%         cimgs{i} = filenames(str, 'absolute');
        
        cimgs{i} = plugin_unzip_images_if_needed(str);
        
            if ~isempty(subjs2exclude_data) % we have subjects to exclude
                if size(cimgs{i},1) == size(idx_include,1) % to be excluded subjects are not missing condition i
                    cimgs{i} = cimgs{i}(idx_include);
                end
            end
    
    % @lukasvo76: this is the fallback option for Windows OS (which does
    % not accept wildcards before the last separator in the path)
    % it requires different definitions of subfolder & functional wildcard in
    % DAT structure set up in prep_1 script, see example there
    % spm_select uses regular expressions as filter . is wildcard, not *!    
    
    else 
        
        str = spm_select('ExtFPListRec',datadir, DAT.functional_wildcard{i}, Inf); 
        
%         % Unzip if needed - not needed in LaBGAS case since we typically do
%         not have zipped con images, although that could be implemented
%         
%         try eval(['!gunzip ' strrep(str, ' ', '\ ') '.gz']), catch, end
%         cimgs{i} = filenames(str, 'absolute');
        
        cimgs{i} = cellstr(str);
        
            for j = 1:size(cimgs{i},1)
                cimgs{i}{j} = cimgs{i}{j}(1,1:end-2); % lukasvo76: gets rid of the ',1' added by spm_select at the end of the filename (first volume, but con images only have one volume)
            end
        
            if ~isempty(subjs2exclude_data) % we have subjects to exclude
                if size(cimgs{i},1) == size(idx_include,1) % to be excluded subjects are not missing condition i
                    cimgs{i} = cimgs{i}(idx_include);
                end
            end
        
    end
    
    %  check whether files exist
    if isempty(cimgs{i}), fprintf('\nLooking in: %s\n', str)
        error('CANNOT FIND IMAGES. Check path names and wildcards.'); 
    end
    
    cimgs{i} = cellfun(@check_valid_imagename, cimgs{i}, repmat({1}, size(cimgs{i}, 1), 1), 'UniformOutput', false);
    
    DAT.imgs{i} = cimgs{i};

end

brainmask = fmri_mask_image(maskname_brain,'noverbose');


%% LOAD FULL OBJECTS AND QC
% -------------------------------------------------------------------------

%%
% *PREP SAMPLING*

% Determine whether we want to sample to the mask (2 x 2 x 2 mm) or native
% space, whichever is more space-efficient

test_image = fmri_data(deblank(DAT.imgs{1}(1, :)), 'noverbose');
voxelsize = diag(test_image.volInfo.mat(1:3, 1:3))';

if prod(abs(voxelsize)) < 8
    sample_type_string = 'sample2mask'; 
    fprintf('\nLoading images into canonical mask space (2 x 2 x 2 mm)\n\n');
else
    sample_type_string = 'native_image_space'; 
    fprintf('\nLoading images in native space (%3.2f x %3.2f x %3.2f mm)\n\n', voxelsize);

end

%%
% *LOAD CONDITION IMAGES INTO FMRI_DATA_ST OBJECT, PERFORM QC, AND PLOT*

fprintf('\n\n');
printhdr('LOADING RAW IMAGES INTO FMRI_DATA_ST OBJECTS');
fprintf('\n\n');

for i = 1:size(DAT.conditions,2)
    
    fprintf('\n\n');
    printhdr(sprintf('Loading raw images: condition #%d, %s', i, DAT.conditions{i}));
    fprintf('\n\n');
    
    DATA_OBJ{i} = fmri_data_st(DAT.imgs{i}, maskname_brain, sample_type_string); % @lukasvo76: changed to @bogpetre's improved data_st object class, and changed to more sparse brainmask
    
    % make sure we are using right variable types (space-saving)
    % NOTE CANlab (old): this is new and could be a source of errors - beta testing!
    % NOTE lukasvo76: this also includes removing empty voxels using the fmri_data.remove_empty function!
    DATA_OBJ{i} = enforce_variable_types(DATA_OBJ{i});
     
    if dozipimages
        % zip original files to save space and delete the unzipped images (we are done using them now).
        for j=1:size(DAT.imgs{i},1)
            gzip(DAT.imgs{i}{j});
            delete(DAT.imgs{i}{j});
        end 
    end
    
    % QUALITY CONTROL METRICS
    
    fprintf('\n\n');
    printhdr(sprintf('QC metrics for images: condition #%d, %s', i, DAT.conditions{i}));
    fprintf('\n\n');
    
    [group_metrics,individual_metrics,values,gwcsf,gwcsfmean,gwcsfl2norm] = qc_metrics_second_level(DATA_OBJ{i});
    
    DAT.quality_metrics_by_condition{i} = group_metrics;
    DAT.gray_white_csf{i} = values;
    
    fprintf('\nSaving quality control metrics in DAT.quality_metrics_by_condition\n');
    fprintf('\nSaving gray, white, CSF means in DAT.gray_white_csf\n\n');
    
    drawnow; snapnow
    
    % PLOT (OPTIONAL)
    
    if dofullplot
        if ischar(DAT.functional_wildcard{i})
            fprintf('\n');
            fprintf('%s\nPlot of raw images: %s\n%s\n', dashes, DAT.functional_wildcard{i}, dashes);  % This fails when trying to pass in a cell array of wildcards - Michael Sun 10/22/2021
            fprintf('\n');
        elseif iscellstr(DAT.functional_wildcard{i}) || isstring(DAT.functional_wildcard{i})
            fprintf('\n');
            fprintf('%s\nPlot of raw images: %s\n%s\n', dashes, DAT.conditions{i}, dashes);
            fprintf('\n');
        end
        
        disp(DATA_OBJ{i}.fullpath)

        plot(DATA_OBJ{i},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run
        
        drawnow; snapnow
        
        if ~omit_histograms
            
            create_figure('histogram');
            set(gcf,'WindowState','maximized');
            hist_han = histogram(DATA_OBJ{i}, 'byimage', 'by_tissue_type');
            
            drawnow; snapnow
            
        end
        
    end
    
    % DERIVED MEASURES
    
    DATA_OBJ{i} = remove_empty(DATA_OBJ{i});
    DAT.globalmeans{i} = mean(DATA_OBJ{i}.dat)';
    DAT.globalstd{i} = std(DATA_OBJ{i}.dat)';
    
    drawnow; snapnow

end


%% Z-SCORE IMAGES, LOAD INTO OBJECTS, AND QC
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('LOADING Z-SCORED IMAGES INTO FMRI_DATA_ST OBJECTS');
fprintf('\n\n');

for i=1:size(DAT.conditions,2)
    
    % Z-SCORING
    
    fprintf('\n\n');
    printhdr(sprintf('Z-scoring images: condition %d, %s', i, DAT.conditions{i}));
    fprintf('\n\n');

    DATA_OBJsc{i} = rescale(DATA_OBJ{i}, 'zscoreimages');

    DATA_OBJsc{i} = enforce_variable_types(DATA_OBJsc{i});

    % QUALITY CONTROL METRICS

    printhdr(sprintf('QC metrics for z-scored images: condition %3.0f, %s', i, DAT.conditions{i}));
    
    [group_metrics,individual_metrics,values,gwcsf,gwcsfmean,gwcsfl2norm] = qc_metrics_second_level(DATA_OBJsc{i});
    
    DAT.sc_quality_metrics_by_condition{i} = group_metrics;
    DAT.sc_gray_white_csf{i} = values;
    
    fprintf('\nSaving quality control metrics in DAT.sc_quality_metrics_by_condition\n');
    fprintf('\nSaving gray, white, CSF means in DAT.sc_gray_white_csf\n\n');
    
    drawnow; snapnow
    
    % PLOT (OPTIONAL)
    
    if dofullplot
        if ischar(DAT.functional_wildcard{i})
            fprintf('\n');
            fprintf('%s\nPlot of z-scored images: %s\n%s\n', dashes, DAT.functional_wildcard{i}, dashes);  % This fails when trying to pass in a cell array of wildcards - Michael Sun 10/22/2021
            fprintf('\n');
        elseif iscellstr(DAT.functional_wildcard{i}) || isstring(DAT.functional_wildcard{i})
            fprintf('\n');
            fprintf('%s\nPlot of z-scored images: %s\n%s\n', dashes, DAT.conditions{i}, dashes);
            fprintf('\n');
        end

        disp(DATA_OBJsc{i}.fullpath)
        
        plot(DATA_OBJsc{i},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run; 
        
        drawnow; snapnow
        
        if ~omit_histograms
            
            create_figure('histogram');
            set(gcf,'WindowState','maximized');
            hist_han = histogram(DATA_OBJsc{i}, 'byimage', 'by_tissue_type');
            drawnow; snapnow
            
        end
        
    end
    
    % DERIVED MEASURES
    
    DATA_OBJsc{i} = remove_empty(DATA_OBJsc{i});
    DAT.sc_globalmeans{i} = mean(DATA_OBJsc{i}.dat)';
    DAT.sc_globalstd{i} = std(DATA_OBJsc{i}.dat)';
    
    drawnow, snapnow

end


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVE UPDATED DAT STRUCTURE IN images_names_and_setup.mat, AND CONDITION DATA OBJECTS IN data_objects(_scaled).mat ');
fprintf('\n\n');

cd(resultsdir); % unannex image_names_and_setup.mat file if already datalad saved to prevent write permission problems
! git annex unannex image_names_and_setup.mat
cd(rootdir);

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, '-append', 'DAT');

savefilenamedata = fullfile(resultsdir, 'data_objects.mat');
save(savefilenamedata, 'DATA_OBJ', '-v7.3');                 % Note: 6/7/17 Tor switched to -v7.3 format by default

savefilenamedata = fullfile(resultsdir, 'data_objects_scaled.mat');
save(savefilenamedata, 'DATA_OBJsc', '-v7.3');
