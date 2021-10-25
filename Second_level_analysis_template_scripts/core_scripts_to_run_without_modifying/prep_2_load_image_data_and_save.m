% load brain
%  get filenames
%  load data
%  apply NPS
% --------------------------------------------------------

% USER OPTIONS
% This is a standard block of code that can be used in multiple scripts.
% Each script will have its own options needed and default values for
% these.
% The code: 
% (1) Checks whether the option variables exist
% (2) Runs a2_set_default_options if any are missing
% (3) Checks again and uses the default options if they are still missing
% (e.g., not specified in an older/incomplete copy of a2_set_default_options)

% Now set in a2_set_default_options
options_needed = {'dofullplot', 'omit_histograms' 'dozipimages'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {true false false};          % defaults if we cannot find info in a2_set_default_options at all; lukasvo76: changed the default for zipping images

plugin_get_options_for_analysis_script


%% Prep and check image names
% -------------------------------------------------------------------------

clear imgs cimgs

for i = 1:length(DAT.conditions)
    
    printhdr(sprintf('Raw data, condition %3.0f, %s', i, DAT.conditions{i}));
    
    % This will vary based on your naming conventions
    % This version assumes that FOLDERS have names of CONDITIONS and images
    % are in each folder
    %
    % @lukasvo76: LaBGAS convention is to have subfolders by subject
    % throughout analysis, so we leave DAT.subfolders blank in prep_1, and
    % I adapted the else part of the if loop below to our conventional
    % folder structure
    
    if ~isempty(DAT.subfolders) && ~isempty(DAT.subfolders{i})  % if we have subfolders
        
        str = fullfile(datadir, DAT.subfolders{i}, DAT.functional_wildcard{i});
        
        % Unzip if needed - 
        % note, Matlab's gunzip() does not remove .gz images, so use eval( ) version.
        % note, replace spaces with '\ ' 
        
        try eval(['!gunzip ' strrep(str, ' ', '\ ') '.gz']), catch, end     % gunzip([str '.gz'])
        cimgs{i} = filenames(str, 'absolute');
        
        cimgs{i} = plugin_unzip_images_if_needed(str);
        
    else
        
        str = spm_select('ExtFPListRec',datadir, DAT.functional_wildcard{i}, Inf); % lukasvo76: spm_select uses regular expressions as filter . is wildcard, not *!
        
%         % Unzip if needed - not needed in LaBGAS case since we typically do
%         not have zipped con images, although that could be implemented
%         
%         try eval(['!gunzip ' strrep(str, ' ', '\ ') '.gz']), catch, end
%         cimgs{i} = filenames(str, 'absolute');
        
        cimgs{i} = cellstr(str);
        
        for j = 1:size(cimgs{i},1)
            cimgs{i}{j} = cimgs{i}{j}(1,1:end-2); % lukasvo76: gets rid of the ',1' added by spm_select at the end of the filename (first volume, but con images only have one volume)
        end
        
    end
    
    %  CHECK that files exist
    if isempty(cimgs{i}), fprintf('Looking in: %s\n', str), error('CANNOT FIND IMAGES. Check path names and wildcards.'); end
    cimgs{i} = cellfun(@check_valid_imagename, cimgs{i}, repmat({1}, size(cimgs{i}, 1), 1), 'UniformOutput', false);
    
    DAT.imgs{i} = cimgs{i};

end


%% Load full objects
% -------------------------------------------------------------------------

% Determine whether we want to sample to the mask (2 x 2 x 2 mm) or native
% space, whichever is more space-efficient

test_image = fmri_data(deblank(DAT.imgs{1}(1, :)), 'noverbose');
voxelsize = diag(test_image.volInfo.mat(1:3, 1:3))';
if prod(abs(voxelsize)) < 8
    sample_type_string = 'sample2mask'; 
    disp('Loading images into canonical mask space (2 x 2 x 2 mm)');
else
    sample_type_string = 'native_image_space'; 
    fprintf('Loading images in native space (%3.2f x %3.2f x %3.2f mm)\n', voxelsize);

end

printhdr('LOADING IMAGES INTO FMRI_DATA_ST OBJECTS');

for i = 1:length(DAT.conditions)
    
    printhdr(sprintf('Loading images: condition %3.0f, %s', i, DAT.conditions{i}));
    printstr(dashes);
    
    % If images are less than 2 mm res, sample in native space:
    DATA_OBJ{i} = fmri_data_st(DAT.imgs{i}); % @lukasvo76: changed to @bogpetre's improved data_st object class
    
    % If images are very large/high-res, you may want to sample to the mask space instead:
%     DATA_OBJ{i} = fmri_data_st(DAT.imgs{i}, which('brainmask.nii'), sample_type_string, 'noverbose');
    
    % make sure we are using right variable types (space-saving)
    % this is new and could be a source of errors - beta testing!
    DATA_OBJ{i} = enforce_variable_types(DATA_OBJ{i});
     
    if dozipimages
        % zip original files to save space and delete the unzipped images (we are done using them now).
        for j=1:size(DAT.imgs{i},1)
            gzip(DAT.imgs{i}{j});
            delete(DAT.imgs{i}{j});
        end 
    end
    
    % QUALITY CONTROL METRICS
    % ---------------------------------------------------------------------

    printhdr(sprintf('QC metrics for images: condition %3.0f, %s', i, DAT.conditions{i}));
    printstr(dashes);
    
    [group_metrics individual_metrics values gwcsf gwcsfmean gwcsfl2norm] = qc_metrics_second_level(DATA_OBJ{i});
    
    DAT.quality_metrics_by_condition{i} = group_metrics;
    DAT.gray_white_csf{i} = values;
    
    disp('Saving quality control metrics in DAT.quality_metrics_by_condition');
    disp('Saving gray, white, CSF means in DAT.gray_white_csf');
    
    drawnow; snapnow
    
    % optional: plot
    % ---------------------------------------------------------------------
    
    if dofullplot
        fprintf('%s\nPlot of images: %s\n%s\n', dashes, DAT.functional_wildcard{i}, dashes);
        disp(DATA_OBJ{i}.fullpath)
        
        plot(DATA_OBJ{i}); drawnow; snapnow
        
        if ~omit_histograms
            
            hist_han = histogram(DATA_OBJ{i}, 'byimage', 'singleaxis');
            title([DAT.conditions{i} ' histograms for each image']);
            drawnow; snapnow
            
            hist_han = histogram(DATA_OBJ{i}, 'byimage', 'by_tissue_type');
            drawnow; snapnow
            
        end
        
    end
    
    % derived measures
    
    DATA_OBJ{i} = remove_empty(DATA_OBJ{i});
    DAT.globalmeans{i} = mean(DATA_OBJ{i}.dat)';
    DAT.globalstd{i} = std(DATA_OBJ{i}.dat)';
    
    drawnow, snapnow
end


%% Z-score images
% -------------------------------------------------------------------------

printhdr('Z-SCORING IMAGES');

for i=1:length(DAT.conditions)
    
    % Z-scoring
    % ---------------------------------------------------------------------
    
    printhdr(sprintf('Z-scoring images: condition %3.0f, %s', i, DAT.conditions{i}));
    printstr(dashes);

    DATA_OBJsc{i} = rescale(DATA_OBJ{i}, 'zscoreimages');

    DATA_OBJsc{i} = enforce_variable_types(DATA_OBJsc{i});

    % QUALITY CONTROL METRICS
    % ---------------------------------------------------------------------

    printhdr(sprintf('QC metrics for z-scored images: condition %3.0f, %s', i, DAT.conditions{i}));
    printstr(dashes);
    
    [group_metrics individual_metrics values gwcsf gwcsfmean gwcsfl2norm] = qc_metrics_second_level(DATA_OBJsc{i});
    
    DAT.sc_quality_metrics_by_condition{i} = group_metrics;
    DAT.sc_gray_white_csf{i} = values;
    
    disp('Saving quality control metrics in DAT.sc_quality_metrics_by_condition');
    disp('Saving gray, white, CSF means in DAT.sc_gray_white_csf');
    
    drawnow; snapnow
    
    % optional: plot
    % ---------------------------------------------------------------------
    
    if dofullplot
        fprintf('%s\nPlot of z-scored images: %s\n%s\n', dashes, DAT.functional_wildcard{i}, dashes);
        disp(DATA_OBJsc{i}.fullpath)
        
        plot(DATA_OBJsc{i}); drawnow; snapnow
        
        if ~omit_histograms
            
            hist_han = histogram(DATA_OBJsc{i}, 'byimage', 'singleaxis');
            title([DAT.conditions{i} ' histograms for each image']);
            drawnow; snapnow
            
            hist_han = histogram(DATA_OBJsc{i}, 'byimage', 'by_tissue_type');
            drawnow; snapnow
            
        end
        
    end
    
    % derived measures
    
    DATA_OBJsc{i} = remove_empty(DATA_OBJsc{i});
    DAT.sc_globalmeans{i} = mean(DATA_OBJsc{i}.dat)';
    DAT.sc_globalstd{i} = std(DATA_OBJsc{i}.dat)';
    
    drawnow, snapnow

end


%% SAVE
% -------------------------------------------------------------------------

printhdr('Save results');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'DAT', 'basedir', 'datadir', 'resultsdir', 'scriptsdir', 'figsavedir');

savefilenamedata = fullfile(resultsdir, 'data_objects.mat');
save(savefilenamedata, 'DATA_OBJ', '-v7.3');                 % Note: 6/7/17 Tor switched to -v7.3 format by default

savefilenamedata = fullfile(resultsdir, 'data_objects_scaled.mat');
save(savefilenamedata, 'DATA_OBJsc', '-v7.3');
