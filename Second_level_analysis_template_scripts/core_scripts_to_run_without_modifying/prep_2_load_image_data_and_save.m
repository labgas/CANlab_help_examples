%% prep_2_load_image_data_and_save.m
%
%
% *USAGE*
%
% This prep script
%
% # calls a_set_up_paths_always_run_first, then loads DAT from
%   image_names_and_setup.mat if not already in the workspace (falling back
%   to prep_1_set_conditions_contrasts_colors and prep_1b_prep_behavioral_data
%   if that file doesn't exist yet)
% # loads first-level beta/con images into CANlab's fmri_data_st objects
% # performs quality control, including plots if requested in a2 script
% # z-scores images and then repeats steps 2-3
% # saves DAT (appended to image_names_and_setup.mat), raw condition images
%   (data_objects.mat), and z-scored condition images
%   (data_objects_scaled.mat) to resultsdir
%
% * the quality-control plots and metrics produced in steps 2-4 should be
%   inspected in the resulting html report before proceeding to
%   prep_3_calc_univariate_contrast_maps_and_save.m - this is the routine
%   check that images loaded correctly
%
% Run this script headless from the Linux command line (the default), which
% publishes the html report and fails loudly if the script errors:
%
%   labgascore_run_headless.sh -d /data/proj_xxx \
%       -s <proj>_secondlevel_m<M>_s0_a_set_up_paths_always_run_first \
%       <proj>_secondlevel_m<M>_s<N>_prep_2_load_image_data_and_save
%
% Or, interactively from the Matlab terminal (use this when you want
% higher-resolution figures, or are debugging):
%
%   LaBGAScore_prov_publish('prep_2_load_image_data_and_save', htmlsavedir)
%
% NOTE: publish() catches a script error into the html and returns normally, so
% a crashed run looks exactly like a successful one. Prefer the routes above,
% which read the report back and check for a caught error, over a bare
% publish('prep_2_load_image_data_and_save','outputDir',htmlsavedir).
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
%                           THIS WILL ALSO REMOVE THOSE SUBJECTS IN
%                           DAT.BEHAVIOR AND DAT.BETWEENPERSON
%
% * docombat                default false, run ComBat harmonization on the RAW condition images before scaling and before contrasts are formed
%
% * combat_batch            batch/site labels, either the name of a column in DAT.BETWEENPERSON.conditions{i} (e.g. 'center') or an n x 1 vector; required if docombat is true
%
% * combat_mod              default {}, cell array of column names in DAT.BETWEENPERSON.conditions{i} whose effects are PRESERVED, e.g. {'group'}; see the warning about classifiers in the ComBat section below
%
% * combat_parametric       default true, parametric (true) or non-parametric (false) empirical Bayes adjustment
%
% * combat_ref_batch        default empty, label of the batch to harmonize towards; empty harmonizes to the grand mean rather than any one site's distribution
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   Dartmouth, May, 2022
%
% -------------------------------------------------------------------------
%
% prep_2_load_image_data_and_save.m         v2.5
%
% last modified: 2026/08/14
%
%
%% RUN SCRIPT A_SET_UP_PATHS_ALWAYS_RUN_FIRST AND LOAD/CREATE DAT IF NEEDED
% -------------------------------------------------------------------------
% Remember where the study's own setup put the results, so the call below can be
% checked against it (see the guard immediately after).
resultsdir_before_setup = '';
if exist('resultsdir','var'), resultsdir_before_setup = resultsdir; end


a_set_up_paths_always_run_first;


% GUARD: did the path setup just move the output directory?
%
% This line is meant to be replaced, in a study's copy, by that study's own
% s0 (e.g. mystudy_secondlevel_m2a_s0_a_set_up_paths_always_run_first). Left
% as the generic call, it RE-DERIVES resultsdir - typically from the
% FIRST-LEVEL model name - and silently overwrites whatever the study's setup
% had already set. Every result then lands in a different model's directory
% while the published report still goes to the right one, so the split is easy
% to miss. This has happened three times: proj_discoverie's SVM wrote into
% secondlevel/model_2_basic, proj_moodbugs wrote into secondlevel/model_3_basic,
% and all seven core scripts of a new discoverie model were about to do the same.
if ~isempty(resultsdir_before_setup) && ~strcmp(resultsdir_before_setup, resultsdir)
    error(['\nPATH SETUP MOVED THE RESULTS DIRECTORY.\n\n' ...
           '  before: %s\n  after : %s\n\n' ...
           'The generic a_set_up_paths_always_run_first re-derived resultsdir and\n' ...
           'discarded the one your study setup had set. In your copy of this script,\n' ...
           'replace that call with your study''s own s0 path script.\n'], ...
           resultsdir_before_setup, resultsdir);
end

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
        
    if isfield(DAT,'BEHAVIOR')
        DAT.BEHAVIOR.behavioral_data_table = DAT.BEHAVIOR.behavioral_data_table(idx_include,:);
    end
    
    if isfield(DAT,'BETWEENPERSON')
        if isfield(DAT.BETWEENPERSON,'group')
            DAT.BETWEENPERSON.group = DAT.BETWEENPERSON.group(idx_include,:);
        end
    end
    
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

        % capture existing figures first: these CANlab calls open more than one
        % (plot(fmri_data) opens canlab_orthviews AND the data-matrix figure), and
        % sizing only gcf leaves the others at their created size. keepaspect stops
        % the wide orthviews panel being stretched to the default 16:10.
        fh_before = findobj('Type','figure');
        plot(DATA_OBJ{i},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run
        % Size every figure plot() produced, found by TAG as well as by newness.
        % Newness alone is not enough: plot() builds its panels with create_figure,
        % which REUSES a figure carrying the same tag rather than opening a new one.
        % From the second loop iteration on, the 6-panel 'fmri data matrix' figure is
        % therefore not new, setdiff missed it, and it never got sized or title-scaled
        % - the titles came out shrunk on the first contrast and full size on every
        % one after it.
        fh_plot = findobj('Type','figure','Tag','fmri data matrix');
        fh_plot = [fh_plot; findobj('Type','figure','Tag','means by condition (unique Y values)')];
        fh_plot = [fh_plot; findobj('Type','figure','Tag','Slice_montage')];
        fh_plot = unique([fh_plot; setdiff(findobj('Type','figure'), fh_before)]);
        if ~isempty(fh_plot)
            plugin_set_figure_size('fig', fh_plot, 'keepaspect', true, ...
                'titlescale', 0.5);   % the 6-panel data-matrix figure's 15.4 pt titles overlap even at the 2/3 default
        end
        
        drawnow; snapnow
        
        if ~omit_histograms
            
            create_figure('histogram');
            hist_han = histogram(DATA_OBJ{i}, 'byimage', 'by_tissue_type');
            fh_dens = findobj('Type','figure','Tag','histogram');
            fh_rel  = findobj('Type','figure','Tag','relationships');
            % The 'histogram' figure is a grid of ONE density plot per subject, so with
            % 64 or 158 subjects each panel becomes unreadable at a fixed canvas size.
            % minpanel grows the canvas until every panel is at least 1.2 x 1.0 inches,
            % reading the actual layout rather than assuming one. The sibling
            % 'relationships' figure (mean/SD by tissue type) keeps the normal sizing -
            % it is a 1x3 layout, and minpanel on three panels asks for a canvas
            % hundreds of inches wide.
            %
            % Both figures are found by TAG across all open figures, deliberately, not
            % by diffing against a snapshot taken before the call. create_figure REUSES
            % an existing figure with the same tag rather than opening a new one, so on
            % every loop iteration after the first the density grid is not a NEW figure.
            % A snapshot-based test therefore skipped it silently, and in one variant
            % sized the 1x3 'relationships' figure instead - which is what produced the
            % 480x420 grids and the 2880x205 ribbon in the first proj_cfs reports.
            if ~isempty(fh_dens)
                plugin_set_figure_size('fig', fh_dens(1), 'minpanel', [1.2 1.0]);
            end
            if ~isempty(fh_rel)
                plugin_set_figure_size('fig', fh_rel(1), 'keepaspect', true);
            end
            
            drawnow; snapnow
            
        end
        
    end
    
    % DERIVED MEASURES
    
    DAT.globalmeans{i} = mean(DATA_OBJ{i}.dat)';
    DAT.globalstd{i} = std(DATA_OBJ{i}.dat)';
    
    drawnow; snapnow

end


%% COMBAT HARMONIZATION OF RAW IMAGES (OPTIONAL)
% -------------------------------------------------------------------------
% Harmonizes multi-site (or multi-scanner) additive and multiplicative
% differences out of the RAW condition images - before any scaling, and before
% contrasts are formed - so that everything downstream inherits harmonized
% data. Runs only if docombat is true.
%
% Requires ComBatHarmonization on the path (Jfortin1/ComBatHarmonization,
% Matlab/scripts/combat.m).
%
% *combat_mod AND SUBSEQUENT DECODING - READ THIS*
% Variables named in combat_mod are PRESERVED: their effects are protected
% from removal as site effects. Include the biological effect of interest
% (typically group) when the harmonized data feed a GLM. Do NOT include it
% when the harmonized data feed a classifier trained on that same variable -
% that leaks label information into the features and inflates accuracy.
% combat.m's own help gives the same warning. For decoding, harmonize inside
% the cross-validation loop instead.
%
% *IDENTIFIABILITY*
% A site whose subjects are all one group carries no within-site contrast, so
% its site effect and the group effect are separable only through the other
% sites. ComBat does not repair such a design: it removes that site's offset
% and leaves the group effect estimated from the sites that do vary, so those
% subjects contribute little independent evidence. The per-batch counts
% printed below show whether this applies.

if ~exist('docombat','var') || isempty(docombat), docombat = false; end

if docombat

    fprintf('\n\n');
    printhdr('COMBAT HARMONIZATION OF RAW IMAGES');
    fprintf('\n\n');

    if isempty(which('combat'))
        error('docombat is true but combat.m is not on the path. Add ComBatHarmonization/Matlab/scripts.');
    end
    if ~exist('combat_batch','var') || isempty(combat_batch)
        error('docombat is true but combat_batch is empty. Set it to a column name in DAT.BETWEENPERSON.conditions{1}, or an n x 1 vector of site labels.');
    end
    if ~exist('combat_mod','var'), combat_mod = {}; end
    if ~exist('combat_parametric','var') || isempty(combat_parametric), combat_parametric = true; end
    if ~exist('combat_ref_batch','var'), combat_ref_batch = []; end
    if ischar(combat_mod) || isstring(combat_mod), combat_mod = cellstr(combat_mod); end

    DAT.combat = struct();
    DAT.combat.mod_vars   = combat_mod;
    DAT.combat.parametric = logical(combat_parametric);
    DAT.combat.ref_batch  = combat_ref_batch;

    for i = 1:size(DAT.conditions,2)

        fprintf('\n\n');
        printhdr(sprintf('ComBat: condition #%d, %s', i, DAT.conditions{i}));
        fprintf('\n\n');

        n_i = size(DATA_OBJ{i}.dat, 2);
        T   = DAT.BETWEENPERSON.conditions{i};

        % --- resolve the batch vector -------------------------------------
        if ischar(combat_batch) || isstring(combat_batch)
            bname = char(combat_batch);
            % Look in the per-condition covariate table first, then in
            % DAT.BETWEENPERSON itself. Site labels usually should NOT be a
            % column of the covariate table: that table becomes the second-level
            % design matrix, and a categorical label column would break it.
            % Keeping the raw labels in DAT.BETWEENPERSON.<name> lets ComBat use
            % them while the design carries dummy columns, or none at all.
            if istable(T) && ismember(bname, T.Properties.VariableNames)
                batch_raw = T.(bname);
            elseif isstruct(DAT.BETWEENPERSON) && isfield(DAT.BETWEENPERSON, bname)
                batch_raw = DAT.BETWEENPERSON.(bname);
            else
                error(['combat_batch ''%s'' was found neither as a column of ' ...
                       'DAT.BETWEENPERSON.conditions{%d} nor as a field of DAT.BETWEENPERSON.'], bname, i);
            end
            DAT.combat.batch_var = bname;
        else
            batch_raw = combat_batch;
            DAT.combat.batch_var = '<supplied as vector>';
        end
        batch_raw = batch_raw(:);
        if numel(batch_raw) ~= n_i
            error('combat_batch has %d entries but condition #%d has %d images.', numel(batch_raw), i, n_i);
        end

        % combat.m compares batch entries with == (find(batch == uniq_batch(i))),
        % so the batch vector MUST be numeric - a cellstr of site names errors
        % there. Map labels to integer codes, keeping the labels for reporting.
        batch_labels_all = cellstr(string(batch_raw));
        [batch_labels, ~, batch] = unique(batch_labels_all, 'stable');
        batch = double(batch);

        % --- resolve the reference batch ----------------------------------
        ref_code = [];
        if ~isempty(combat_ref_batch)
            ref_code = find(strcmp(batch_labels, char(string(combat_ref_batch))));
            if isempty(ref_code)
                error('combat_ref_batch ''%s'' is not one of the batches present (%s).', ...
                    char(string(combat_ref_batch)), strjoin(batch_labels(:)', ', '));
            end
        end

        % --- build the preserved-effects design ---------------------------
        mod = [];
        for k = 1:numel(combat_mod)
            if ~istable(T) || ~ismember(combat_mod{k}, T.Properties.VariableNames)
                error('combat_mod variable ''%s'' is not in DAT.BETWEENPERSON.conditions{%d}.', combat_mod{k}, i);
            end
            mod = [mod double(T.(combat_mod{k})(:))]; %#ok<AGROW>
        end

        % --- report the design --------------------------------------------
        fprintf('  batch variable : %s\n', DAT.combat.batch_var);
        if isempty(combat_mod)
            fprintf('  preserved      : none\n');
        else
            fprintf('  preserved      : %s\n', strjoin(combat_mod, ', '));
        end
        if isempty(ref_code)
            fprintf('  reference      : none (harmonized to grand mean)\n');
        else
            fprintf('  reference      : %s\n', batch_labels{ref_code});
        end
        if logical(combat_parametric)
            fprintf('  adjustment     : parametric empirical Bayes\n\n');
        else
            fprintf('  adjustment     : non-parametric empirical Bayes\n\n');
        end
        for b = 1:numel(batch_labels)
            sel = batch == b;
            lvlstr = '';
            for k = 1:size(mod,2)
                u = unique(mod(sel,k));
                lvlstr = [lvlstr sprintf('   %s: %s', combat_mod{k}, mat2str(u(:)'))]; %#ok<AGROW>
            end
            fprintf('    %-12s n = %3d%s\n', batch_labels{b}, sum(sel), lvlstr);
        end
        if ~isempty(mod)
            fprintf(['\n    a batch showing a single level of a preserved variable carries no\n' ...
                     '    within-batch contrast for it - see the identifiability note above\n']);
        end

        % --- harmonize ------------------------------------------------------
        dat = DATA_OBJ{i}.dat;

        % ComBat standardizes by the within-batch variance, so a voxel that is
        % constant within ANY batch yields Inf/NaN. Harmonize only voxels that
        % vary in every batch and pass the rest through untouched.
        ok = true(size(dat,1),1);
        for b = 1:numel(batch_labels)
            ok = ok & (std(double(dat(:, batch == b)), 0, 2) > 0);
        end
        if ~all(ok)
            fprintf('\n  %d of %d voxels are constant within at least one batch; left unharmonized\n', ...
                sum(~ok), numel(ok));
        end

        cb_args = {double(dat(ok,:)), batch, mod, double(logical(combat_parametric))};
        if ~isempty(ref_code), cb_args = [cb_args {'ref', ref_code}]; end %#ok<AGROW>

        fprintf('\n');
        [harmonized, gamma_star, delta_star, gamma_hat, delta_hat] = combat(cb_args{:});

        dat(ok,:) = harmonized;
        DATA_OBJ{i}.dat = dat;
        DATA_OBJ{i} = enforce_variable_types(DATA_OBJ{i});

        % combat.m's help: a dramatic empirical -> posterior shift means the
        % priors are driving the fit rather than the data. Print both so that
        % is visible in the report rather than buried in the returned structs.
        fprintf('\n  empirical -> posterior batch parameters (mean over voxels):\n');
        for b = 1:numel(batch_labels)
            fprintf('    %-12s gamma %+8.4f -> %+8.4f    delta %8.4f -> %8.4f\n', ...
                batch_labels{b}, mean(gamma_hat(b,:)), mean(gamma_star(b,:)), ...
                mean(delta_hat(b,:)), mean(delta_star(b,:)));
        end

        DAT.combat.batch_labels      = batch_labels;
        DAT.combat.n_per_batch       = accumarray(batch, 1)';
        DAT.combat.n_voxels_adjusted(i) = sum(ok);
        DAT.combat.gamma_hat_mean{i}    = mean(gamma_hat, 2)';
        DAT.combat.gamma_star_mean{i}   = mean(gamma_star, 2)';
        DAT.combat.delta_hat_mean{i}    = mean(delta_hat, 2)';
        DAT.combat.delta_star_mean{i}   = mean(delta_star, 2)';

        drawnow; snapnow

    end

    DAT.combat.applied = true;
    fprintf('\n\nComBat applied to raw condition images; all downstream scaling and contrasts inherit harmonized data\n\n');

else

    DAT.combat.applied = false;

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
        
        % capture existing figures first: these CANlab calls open more than one
        % (plot(fmri_data) opens canlab_orthviews AND the data-matrix figure), and
        % sizing only gcf leaves the others at their created size. keepaspect stops
        % the wide orthviews panel being stretched to the default 16:10.
        fh_before = findobj('Type','figure');
        plot(DATA_OBJsc{i},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run; 
        % Size every figure plot() produced, found by TAG as well as by newness.
        % Newness alone is not enough: plot() builds its panels with create_figure,
        % which REUSES a figure carrying the same tag rather than opening a new one.
        % From the second loop iteration on, the 6-panel 'fmri data matrix' figure is
        % therefore not new, setdiff missed it, and it never got sized or title-scaled
        % - the titles came out shrunk on the first contrast and full size on every
        % one after it.
        fh_plot = findobj('Type','figure','Tag','fmri data matrix');
        fh_plot = [fh_plot; findobj('Type','figure','Tag','means by condition (unique Y values)')];
        fh_plot = [fh_plot; findobj('Type','figure','Tag','Slice_montage')];
        fh_plot = unique([fh_plot; setdiff(findobj('Type','figure'), fh_before)]);
        if ~isempty(fh_plot)
            plugin_set_figure_size('fig', fh_plot, 'keepaspect', true, ...
                'titlescale', 0.5);   % the 6-panel data-matrix figure's 15.4 pt titles overlap even at the 2/3 default
        end
        
        drawnow; snapnow
        
        if ~omit_histograms
            
            create_figure('histogram');
            hist_han = histogram(DATA_OBJsc{i}, 'byimage', 'by_tissue_type');
            fh_dens = findobj('Type','figure','Tag','histogram');
            fh_rel  = findobj('Type','figure','Tag','relationships');
            % The 'histogram' figure is a grid of ONE density plot per subject, so with
            % 64 or 158 subjects each panel becomes unreadable at a fixed canvas size.
            % minpanel grows the canvas until every panel is at least 1.2 x 1.0 inches,
            % reading the actual layout rather than assuming one. The sibling
            % 'relationships' figure (mean/SD by tissue type) keeps the normal sizing -
            % it is a 1x3 layout, and minpanel on three panels asks for a canvas
            % hundreds of inches wide.
            %
            % Both figures are found by TAG across all open figures, deliberately, not
            % by diffing against a snapshot taken before the call. create_figure REUSES
            % an existing figure with the same tag rather than opening a new one, so on
            % every loop iteration after the first the density grid is not a NEW figure.
            % A snapshot-based test therefore skipped it silently, and in one variant
            % sized the 1x3 'relationships' figure instead - which is what produced the
            % 480x420 grids and the 2880x205 ribbon in the first proj_cfs reports.
            if ~isempty(fh_dens)
                plugin_set_figure_size('fig', fh_dens(1), 'minpanel', [1.2 1.0]);
            end
            if ~isempty(fh_rel)
                plugin_set_figure_size('fig', fh_rel(1), 'keepaspect', true);
            end
            drawnow; snapnow
            
        end
        
    end
    
    % DERIVED MEASURES
    
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
