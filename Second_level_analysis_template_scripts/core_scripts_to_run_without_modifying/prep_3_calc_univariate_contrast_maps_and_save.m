%% prep_3_calc_univariate_contrast_maps_and_save.m
%
%
% *USAGE*
%
% This prep script
%
% # calls a_set_up_paths_always_run_first, then loads DAT/DATA_OBJ/DATA_OBJsc
%   from image_names_and_setup.mat/data_objects.mat/data_objects_scaled.mat
%   if not already in the workspace (falling back to prep_1/prep_1b/prep_2
%   if those files don't exist yet)
% # calculates contrast images from prep_2's raw condition images (DATA_OBJ)
%   and z-scored condition images (DATA_OBJsc), storing them as CANlab's
%   fmri_data_st objects in DATA_OBJ_CON and DATA_OBJ_CONsc respectively
% # l2norm-rescales the raw contrast images into DATA_OBJ_CONscc
% # performs quality control, including plots if requested in a2 script, on
%   all three contrast-object variants
% # saves DATA_OBJ_CON/DATA_OBJ_CONsc/DATA_OBJ_CONscc to
%   contrast_data_objects.mat, and appends DAT (including the new
%   DAT.gray_white_csf_contrasts field) to image_names_and_setup.mat
%
% * the quality-control plots and metrics produced in step 4 should be
%   inspected in the resulting html report before proceeding to
%   prep_3a_run_second_level_regression_and_save.m (or other downstream
%   scripts) - this is the routine check that images loaded and contrasts
%   were computed correctly
%
% Run this script headless from the Linux command line (the default), which
% publishes the html report and fails loudly if the script errors:
%
%   labgascore_run_headless.sh -d /data/proj_xxx \
%       -s <proj>_secondlevel_m<M>_s0_a_set_up_paths_always_run_first \
%       <proj>_secondlevel_m<M>_s<N>_prep_3_calc_univariate_contrast_maps_and_save
%
% Or, interactively from the Matlab terminal (use this when you want
% higher-resolution figures, or are debugging):
%
%   LaBGAScore_prov_publish('prep_3_calc_univariate_contrast_maps_and_save', htmlsavedir)
%
% NOTE: publish() catches a script error into the html and returns normally, so
% a crashed run looks exactly like a successful one. Prefer the routes above,
% which read the report back and check for a caught error, over a bare
% publish('prep_3_calc_univariate_contrast_maps_and_save','outputDir',htmlsavedir).
%
% *NOTE*
%   We can include image sets with different numbers of images, as occurs with between-person designs, as
%   long as the contrast weights are zero for all elements with different numbers of images.
%
%
% *OPTIONS*
%
% * docombat_contrasts      default false, run ComBat on the CONTRAST images after they are formed.
%                           Independent of docombat in prep_2; shares combat_batch / combat_mod /
%                           combat_parametric / combat_ref_batch. Harmonizing conditions does NOT
%                           harmonize contrasts, because contrast variance depends on the
%                           between-condition covariance that condition-level ComBat leaves alone.
%
% * dofullplot              default true, can set to false to save time, but not recommended for quality control purposes
%
% * omit_histograms         default false, can set to true to save time, especially in case of large samples but not recommended for quality control purposes
%
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   Dartmouth, May, 2022
%
% -------------------------------------------------------------------------
% prep_3_calc_univariate_contrast_maps_and_save.m         v2.1
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

if ~exist('DATA_OBJ','var')
    
    try
    
        load(fullfile(resultsdir,'data_objects.mat'));
    
    catch
        
        prep_2_load_image_data_and_save;
        
    end
    
end

if ~exist('DATA_OBJsc','var')
    
    try
    
        load(fullfile(resultsdir,'data_objects_scaled.mat'));
    
    catch
        
        prep_2_load_image_data_and_save;
        
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

options_needed = {'dofullplot', 'omit_histograms'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed);        % initializing this means a2_set_defaults_options will never run

option_default_values = {true false};          % defaults if we cannot find info in a2_set_default_options at all; @lukasvo76: changed the default for zipping images

plugin_get_options_for_analysis_script


%% CONDITION-LEVEL HARMONIZATION (OPTIONAL) - BEFORE ANY CONTRAST IS FORMED
% -------------------------------------------------------------------------
% Two independent options, both acting on the CONDITION objects (DATA_OBJ and
% DATA_OBJsc) before the contrast loops below:
%
%   prescale_conditions   per-site prescaling of the condition images
%   docombat_conditions   ComBat on the condition images
%
% They are independent of prescale/ComBat at the CONTRAST level further down,
% so the full grid of "where is each operation applied" is reachable:
%
%   ComBat on conditions, contrasts formed after   -> docombat_conditions
%   ComBat on contrasts                            -> docombat_contrasts
%   prescale conditions, ComBat on contrasts       -> prescale_conditions + docombat_contrasts
%   prescale contrasts,  ComBat on contrasts       -> combat_prescale_sites + docombat_contrasts
%
% Both objects are harmonized, so myscaling_glm then selects whether the GLM
% sees contrasts built from RAW or from Z-SCORED conditions - exactly the
% arrangement used at contrast level, which keeps the two levels comparable.
%
% NOTE harmonizing conditions does NOT harmonize contrasts: contrast variance is
% var(A)+var(B)-2cov(A,B) and neither option here touches the between-condition
% covariance. Measured on proj_discoverie, condition-level ComBat left UM at
% 0.708 of KUL's within-case contrast variance while contrast-level reached
% 0.872. That is a reason to compare the levels, not to prefer one blindly.

if ~exist('prescale_conditions','var') || isempty(prescale_conditions), prescale_conditions = false; end
if ~exist('docombat_conditions','var') || isempty(docombat_conditions), docombat_conditions = false; end

if prescale_conditions || docombat_conditions

    fprintf('\n\n');
    printhdr('CONDITION-LEVEL HARMONIZATION');
    fprintf('\n\n');

    if docombat_conditions && isempty(which('combat'))
        error('docombat_conditions is true but combat.m is not on the path.');
    end
    if ~exist('combat_batch','var') || isempty(combat_batch)
        error('condition-level harmonization requires combat_batch.');
    end
    if ~exist('combat_mod','var'), combat_mod = {}; end
    if ~exist('combat_parametric','var') || isempty(combat_parametric), combat_parametric = true; end
    if ~exist('combat_ref_batch','var'), combat_ref_batch = []; end
    if ischar(combat_mod) || isstring(combat_mod), combat_mod = cellstr(combat_mod); end

    DAT.combat_conditions = struct('prescaled', prescale_conditions, ...
        'combat', docombat_conditions, 'mod_vars', {combat_mod}, 'ref_batch', combat_ref_batch);

    for i = 1:size(DAT.conditions,2)

        fprintf('\n');
        printhdr(sprintf('condition #%d, %s', i, DAT.conditions{i}));

        Ti = DAT.BETWEENPERSON.conditions{i};
        if ischar(combat_batch) || isstring(combat_batch)
            bname = char(combat_batch);
            if istable(Ti) && ismember(bname, Ti.Properties.VariableNames)
                batch_raw = Ti.(bname);
            elseif isstruct(DAT.BETWEENPERSON) && isfield(DAT.BETWEENPERSON, bname)
                batch_raw = DAT.BETWEENPERSON.(bname);
            else
                error('combat_batch ''%s'' not found for condition %d.', bname, i);
            end
        else
            batch_raw = combat_batch;
        end
        [batch_labels, ~, batch] = unique(cellstr(string(batch_raw(:))), 'stable');
        batch = double(batch);

        ref_code = [];
        if ~isempty(combat_ref_batch)
            ref_code = find(strcmp(batch_labels, char(string(combat_ref_batch))));
            if isempty(ref_code)
                error('combat_ref_batch ''%s'' is not among the batches present.', char(string(combat_ref_batch)));
            end
        end

        mod = [];
        for k = 1:numel(combat_mod)
            if ~istable(Ti) || ~ismember(combat_mod{k}, Ti.Properties.VariableNames)
                error('combat_mod variable ''%s'' is not in DAT.BETWEENPERSON.conditions{%d}.', combat_mod{k}, i);
            end
            mod = [mod double(Ti.(combat_mod{k})(:))]; %#ok<AGROW>
        end

        for o = 1:2
            if o == 1
                if ~exist('DATA_OBJ','var') || numel(DATA_OBJ) < i, continue; end
                Xc = double(DATA_OBJ{i}.dat);  oname = 'DATA_OBJ  ';
            else
                if ~exist('DATA_OBJsc','var') || numel(DATA_OBJsc) < i, continue; end
                Xc = double(DATA_OBJsc{i}.dat); oname = 'DATA_OBJsc';
            end
            if size(Xc,2) ~= numel(batch)
                error('%s{%d} has %d images but batch has %d entries.', oname, i, size(Xc,2), numel(batch));
            end

            if prescale_conditions
                if isempty(ref_code)
                    error('prescale_conditions requires combat_ref_batch.');
                end
                subj_rms = sqrt(mean(Xc.^2, 1));
                site_mag = zeros(1, numel(batch_labels));
                for b = 1:numel(batch_labels), site_mag(b) = median(subj_rms(batch == b)); end
                fac = site_mag(ref_code) ./ site_mag;
                for b = 1:numel(batch_labels)
                    Xc(:, batch == b) = Xc(:, batch == b) * fac(b);
                    fprintf('  %s prescale %-10s factor %6.3f\n', oname, batch_labels{b}, fac(b));
                end
            end

            if docombat_conditions
                ok = true(size(Xc,1),1);
                for b = 1:numel(batch_labels)
                    ok = ok & (std(Xc(:, batch == b), 0, 2) > 0);
                end
                cb = {Xc(ok,:), batch, mod, double(logical(combat_parametric))};
                if ~isempty(ref_code), cb = [cb {'ref', ref_code}]; end %#ok<AGROW>
                fprintf('  %s ComBat on %d of %d voxels\n', oname, sum(ok), numel(ok));
                Xc(ok,:) = combat(cb{:});
            end

            if o == 1
                DATA_OBJ{i}.dat = Xc;   DATA_OBJ{i} = enforce_variable_types(DATA_OBJ{i});
            else
                DATA_OBJsc{i}.dat = Xc; DATA_OBJsc{i} = enforce_variable_types(DATA_OBJsc{i});
            end
        end
    end

    DAT.combat_conditions.applied = true;
    fprintf('\n\nCondition-level harmonization done; contrasts below are formed from these images\n\n');

else
    DAT.combat_conditions.applied = false;
end


%% RAW AND L2NORM-RESCALED CONTRAST IMAGES FROM RAW CONDITION IMAGES
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('CALCULATING CONTRAST IMAGES FROM RAW CONDITION IMAGES AND CONVERTING TO FMRI_DATA_ST OBJECTS');
fprintf('\n\n');

if ~isfield(DAT, 'contrasts') || isempty(DAT.contrasts)
    % skip
    return
end

k = size(DAT.conditions,2);

%%
% *GET SIZES OF DATA_OBJ*

clear sz

for i = 1:k
    sz(i, :) = size(DATA_OBJ{i}.dat); 
end

sz = sz(:, 2);

for i = 1:k
    DATA_OBJ{i} = replace_empty(DATA_OBJ{i},'voxels');
end


%%
% CREATE DATA_OBJ_CON, RESCALE, AND QC

for c = 1:size(DAT.contrasts, 1)
    
    fprintf('\n');
    fprintf('%s\nCONTRAST: %s\n%s\n', dashes, upper(DAT.contrastnames{c}), dashes);
    fprintf('\n');
    
    % PREP WORK
    % ---------
    
    % initialize : shell object, keep same space/volume info
    wh = find(DAT.contrasts(c, :));
    
    my_size = sz(wh(1));  
    
    % check sizes and make sure they are the same
    if ~all(sz(wh) == sz(wh(1)))
        fprintf('\nNot all image set sizes are the same for contrast %d\n\n', c);
    end
    
    % CREATE CONTRAST OBJECTS & RESCALE BY L2NORM
    % -------------------------------------------
    
    fprintf('\n');
    fprintf('%s\nCreating fmri_data_st object for raw contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
    fprintf('\n');
    
    DATA_OBJ_CON{c} = DATA_OBJ{wh(1)};
    [DATA_OBJ_CON{c}.image_names, DATA_OBJ_CON{c}.fullpath] = deal([]);
        
    DATA_OBJ_CON{c}.dat = zeros(size(DATA_OBJ{wh(1)}.dat));
    
    for i = 1:k
        
        % add data * contrast weight
        condat = DATA_OBJ{i}.dat .* DAT.contrasts(c, i);
        
        if DAT.contrasts(c, i) == 0
            % Skip.  This allows us to include image sets with different
            % numbers of images, as occurs with between-person designs, as
            % long as the contrast weights are zero for all elements with
            % different numbers of images.
            continue
        end
        
        if size(condat, 2) ~= my_size
            fprintf('\nCondition %d : number of images does not match. Check DATA_OBJ images and contrasts\n', i)
            error('exiting...')
        end
        
        DATA_OBJ_CON{c}.dat = DATA_OBJ_CON{c}.dat + condat;
        
    end
    
    DATA_OBJ_CON{c}.image_names = DAT.contrastnames;
    DATA_OBJ_CON{c}.source_notes = DAT.contrastnames;
    
    % rescale contrast objects by l2norm    % added by @lukasvo76 01/03/21
    fprintf('\n');
    fprintf('%s\nRescaling fmri_data_st object by l2norm for raw contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
    fprintf('\n');
    
    DATA_OBJ_CONscc{c} = rescale(DATA_OBJ_CON{c}, 'l2norm_images');
    
    % enforce variable types in objects to save space
    DATA_OBJ_CON{c} = enforce_variable_types(DATA_OBJ_CON{c}); 
    DATA_OBJ_CONscc{c} = enforce_variable_types(DATA_OBJ_CONscc{c}); 
    
    
    % QUALITY CONTROL METRICS & PLOT (OPTIONAL)
    % -----------------------------------------
    
    % RAW CONTRAST OBJECTS
    
    fprintf('\n');
    fprintf('%s\nQC metrics for raw contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
    fprintf('\n');
    
    % qc
    [group_metrics, individual_metrics, gwcsf, gwcsfmean] = qc_metrics_second_level(DATA_OBJ_CON{c});
    drawnow; snapnow
    
    % plot
    if dofullplot
        fprintf('\n');
        fprintf('%s\nPlot of raw contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
        fprintf('\n');
        
        disp(DATA_OBJ_CON{c}.fullpath)
        
        % capture existing figures first: these CANlab calls open more than one
        % (plot(fmri_data) opens canlab_orthviews AND the data-matrix figure), so
        % sizing only gcf leaves the others as created. keepaspect stops the wide
        % orthviews panel being stretched to the default 16:10.
        fh_before = findobj('Type','figure');
        plot(DATA_OBJ_CON{c},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run
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
            
            % The density-plot figure is a grid of ONE panel per subject, so at a fixed
            % canvas each panel becomes unreadable as n grows. minpanel grows the canvas
            % until every panel is at least 1.2 x 1.0 inches, reading the actual layout
            % rather than assuming one. Sibling figures ('relationships', 45 panels) keep
            % the normal sizing - matching their per-panel size would make them absurd.
            create_figure('histogram');
            hist_han = histogram(DATA_OBJ_CON{c}, 'byimage', 'by_tissue_type');
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
    
    % RESCALED CONTRAST OBJECTS
    
    fprintf('\n');
    fprintf('%s\nQC metrics for l2norm-rescaled contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
    fprintf('\n');
    
    % qc
    [group_metrics, individual_metrics, gwcsf, gwcsfmean] = qc_metrics_second_level(DATA_OBJ_CONscc{c});
    drawnow; snapnow
    
    % plot
    if dofullplot
        fprintf('\n');
        fprintf('%s\nPlot of l2norm-rescaled contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
        fprintf('\n');
        
        disp(DATA_OBJ_CONscc{c}.fullpath)
        
        % capture existing figures first: these CANlab calls open more than one
        % (plot(fmri_data) opens canlab_orthviews AND the data-matrix figure), so
        % sizing only gcf leaves the others as created. keepaspect stops the wide
        % orthviews panel being stretched to the default 16:10.
        fh_before = findobj('Type','figure');
        plot(DATA_OBJ_CONscc{c},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run
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
            
            % The density-plot figure is a grid of ONE panel per subject, so at a fixed
            % canvas each panel becomes unreadable as n grows. minpanel grows the canvas
            % until every panel is at least 1.2 x 1.0 inches, reading the actual layout
            % rather than assuming one. Sibling figures ('relationships', 45 panels) keep
            % the normal sizing - matching their per-panel size would make them absurd.
            create_figure('histogram_l2norm');
            hist_han_l2norm = histogram(DATA_OBJ_CONscc{c}, 'byimage', 'by_tissue_type');
            fh_dens = findobj('Type','figure','Tag','histogram_l2norm');
            fh_rel  = findobj('Type','figure','Tag','relationships');
            % minpanel grows the canvas until every panel of the per-subject density grid
            % is at least 1.2 x 1.0 inches, reading the actual layout rather than assuming
            % one. The sibling 'relationships' figure (mean/SD by tissue type) keeps the
            % normal sizing - it is a 1x3 layout, and minpanel on three panels asks for a
            % canvas hundreds of inches wide.
            %
            % Both are found by TAG across all open figures, deliberately, not by diffing
            % against a snapshot taken before the call. create_figure REUSES an existing
            % figure with the same tag rather than opening a new one, so on every loop
            % iteration after the first the density grid is not a NEW figure. A snapshot
            % test therefore skipped it silently, and in one variant sized the 1x3
            % 'relationships' figure instead - which is what produced the 480x420 grids
            % and the 2880x205 ribbon in the first proj_cfs reports.
            if ~isempty(fh_dens)
                plugin_set_figure_size('fig', fh_dens(1), 'minpanel', [1.2 1.0]);
            end
            if ~isempty(fh_rel)
                plugin_set_figure_size('fig', fh_rel(1), 'keepaspect', true);
            end

            drawnow; snapnow
            
        end
        
    end
    
end


%% CONTRAST IMAGES FROM Z-SCORED CONDITION IMAGES
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('CALCULATING CONTRAST IMAGES FROM Z-SCORED CONDITION IMAGES AND CONVERTING TO FMRI_DATA_ST OBJECTS');
fprintf('\n\n');

for i = 1:k
    DATA_OBJsc{i} = replace_empty(DATA_OBJsc{i});
end


%%
% CREATE DATA_OBJ_CONsc, AND QC

for c = 1:size(DAT.contrasts, 1)
    
    fprintf('\n');
    fprintf('%s\nCONTRAST: %s\n%s\n', dashes, upper(DAT.contrastnames{c}), dashes);
    fprintf('\n');

    % PREP
    % ----
    
    fprintf('\n');
    fprintf('%s\nCreating fmri_data_st object for z-scored contrast: %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
    fprintf('\n');
    
    % initialize : shell object, keep same space/volume info
    wh = find(DAT.contrasts(c, :));
    
    my_size = sz(wh(1));  
    
    % check sizes and make sure they are the same
    if ~all(sz(wh) == sz(wh(1)))
        fprintf('\nNot all image set sizes are the same for contrast %d\n\n', c);
    end
    
    % CREATE CONTRAST OBJECTS
    
    DATA_OBJ_CONsc{c} = DATA_OBJsc{wh(1)};
    [DATA_OBJ_CONsc{c}.image_names, DATA_OBJ_CONsc{c}.fullpath] = deal([]);
    DATA_OBJ_CONsc{c}.dat = zeros(size(DATA_OBJsc{wh(1)}.dat));
    
    
    for i = 1:k
        
        % add data * contrast weight
        condat = DATA_OBJsc{i}.dat .* DAT.contrasts(c, i);
        
        if DAT.contrasts(c, i) == 0
            % Skip.  This allows us to include image sets with different
            % numbers of images, as occurs with between-person designs, as
            % long as the contrast weights are zero for all elements with
            % different numbers of images.
            continue
        end
        
        if size(condat, 2) ~= my_size
            fprintf('Condition %3.0f : number of images does not match. Check DATA_OBJsc images and contrasts.', i);
            error('exiting.')
        end
        
        DATA_OBJ_CONsc{c}.dat = DATA_OBJ_CONsc{c}.dat + condat;
        
    end
    
    DATA_OBJ_CONsc{c}.image_names = DAT.contrastnames;
    DATA_OBJ_CONsc{c}.source_notes = DAT.contrastnames;
    
    % enforce variable types in objects to save space
    DATA_OBJ_CONsc{c} = enforce_variable_types(DATA_OBJ_CONsc{c}); 

    
    % QUALITY CONTROL METRICS & PLOT (OPTIONAL)
    
    fprintf('\n');
    fprintf('%s\nQC metrics for contrast (from z-scored condition images): %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
    fprintf('\n');
    
    % qc
    [group_metrics, individual_metrics, gwcsf, gwcsfmean] = qc_metrics_second_level(DATA_OBJ_CONsc{c});
    drawnow; snapnow
    
    % plot
    if dofullplot
        fprintf('\n');
        fprintf('%s\nPlot of contrast (from z-scored condition images): %s\n%s\n', dashes, DAT.contrastnames{c}, dashes);
        disp(DATA_OBJ_CONsc{c}.fullpath)
        
        % capture existing figures first: these CANlab calls open more than one
        % (plot(fmri_data) opens canlab_orthviews AND the data-matrix figure), so
        % sizing only gcf leaves the others as created. keepaspect stops the wide
        % orthviews panel being stretched to the default 16:10.
        fh_before = findobj('Type','figure');
        plot(DATA_OBJ_CONsc{c},'norunmontages'); % @lukasvo76 turned run montages off, since second level con images are most often not per run
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
            
            % The density-plot figure is a grid of ONE panel per subject, so at a fixed
            % canvas each panel becomes unreadable as n grows. minpanel grows the canvas
            % until every panel is at least 1.2 x 1.0 inches, reading the actual layout
            % rather than assuming one. Sibling figures ('relationships', 45 panels) keep
            % the normal sizing - matching their per-panel size would make them absurd.
            create_figure('histogram_zscore');
            hist_han_zscore = histogram(DATA_OBJ_CONsc{c}, 'byimage', 'by_tissue_type');
            fh_dens = findobj('Type','figure','Tag','histogram_zscore');
            fh_rel  = findobj('Type','figure','Tag','relationships');
            % minpanel grows the canvas until every panel of the per-subject density grid
            % is at least 1.2 x 1.0 inches, reading the actual layout rather than assuming
            % one. The sibling 'relationships' figure (mean/SD by tissue type) keeps the
            % normal sizing - it is a 1x3 layout, and minpanel on three panels asks for a
            % canvas hundreds of inches wide.
            %
            % Both are found by TAG across all open figures, deliberately, not by diffing
            % against a snapshot taken before the call. create_figure REUSES an existing
            % figure with the same tag rather than opening a new one, so on every loop
            % iteration after the first the density grid is not a NEW figure. A snapshot
            % test therefore skipped it silently, and in one variant sized the 1x3
            % 'relationships' figure instead - which is what produced the 480x420 grids
            % and the 2880x205 ribbon in the first proj_cfs reports.
            if ~isempty(fh_dens)
                plugin_set_figure_size('fig', fh_dens(1), 'minpanel', [1.2 1.0]);
            end
            if ~isempty(fh_rel)
                plugin_set_figure_size('fig', fh_rel(1), 'keepaspect', true);
            end

            drawnow; snapnow
            
        end
        
    end
    
end


%% COMBAT HARMONIZATION OF CONTRAST IMAGES (OPTIONAL)
% -------------------------------------------------------------------------
% Harmonizes site differences at the CONTRAST level, after the contrasts have
% been formed. Runs only if docombat_contrasts is true. Independent of
% docombat in prep_2: a model may harmonize conditions, contrasts, both, or
% neither. Batch, mod, parametric and reference settings are shared with
% prep_2 (combat_batch / combat_mod / combat_parametric / combat_ref_batch).
%
% *WHY THIS EXISTS SEPARATELY FROM prep_2*
% Harmonizing the CONDITIONS does not harmonize the CONTRASTS. Contrast
% variance is var(A) + var(B) - 2cov(A,B), and ComBat on conditions equalizes
% the two marginal variances while leaving the between-condition covariance
% untouched. Measured on proj_discoverie after condition-level ComBat: UM's
% marginal variances came to 0.84 and 0.87 of KUL's, but its stress-control
% correlation is 0.801 against KUL's 0.758, and because the contrast leverages
% 2(1-r) that gap pulled UM's contrast variance down to 0.706. Harmonizing
% here acts directly on what the second-level GLM actually consumes.
%
% The same caveats as prep_2 apply: combat_mod PRESERVES the named effects, so
% do not name a variable you intend to decode, and a site whose subjects are
% all one level of a preserved variable has no within-site contrast for it.

if ~exist('docombat_contrasts','var') || isempty(docombat_contrasts), docombat_contrasts = false; end

if docombat_contrasts

    fprintf('\n\n');
    printhdr('COMBAT HARMONIZATION OF CONTRAST IMAGES');
    fprintf('\n\n');

    if isempty(which('combat'))
        error('docombat_contrasts is true but combat.m is not on the path. Add ComBatHarmonization/Matlab/scripts.');
    end
    if ~exist('combat_batch','var') || isempty(combat_batch)
        error('docombat_contrasts is true but combat_batch is empty.');
    end
    if ~exist('combat_mod','var'), combat_mod = {}; end
    if ~exist('combat_parametric','var') || isempty(combat_parametric), combat_parametric = true; end
    if ~exist('combat_ref_batch','var'), combat_ref_batch = []; end
    if ~exist('combat_prescale_sites','var') || isempty(combat_prescale_sites), combat_prescale_sites = false; end
    if ischar(combat_mod) || isstring(combat_mod), combat_mod = cellstr(combat_mod); end

    DAT.combat_contrasts = struct();
    DAT.combat_contrasts.mod_vars   = combat_mod;
    DAT.combat_contrasts.parametric = logical(combat_parametric);
    DAT.combat_contrasts.ref_batch  = combat_ref_batch;

    for c = 1:size(DAT.contrasts,1)

        fprintf('\n\n');
        printhdr(sprintf('ComBat on contrast #%d, %s', c, DAT.contrastnames{c}));
        fprintf('\n\n');

        Tc = DAT.BETWEENPERSON.contrasts{c};

        % --- batch vector, same resolution rule as prep_2 -------------------
        if ischar(combat_batch) || isstring(combat_batch)
            bname = char(combat_batch);
            if istable(Tc) && ismember(bname, Tc.Properties.VariableNames)
                batch_raw = Tc.(bname);
            elseif isstruct(DAT.BETWEENPERSON) && isfield(DAT.BETWEENPERSON, bname)
                batch_raw = DAT.BETWEENPERSON.(bname);
            else
                error(['combat_batch ''%s'' was found neither as a column of ' ...
                       'DAT.BETWEENPERSON.contrasts{%d} nor as a field of DAT.BETWEENPERSON.'], bname, c);
            end
            DAT.combat_contrasts.batch_var = bname;
        else
            batch_raw = combat_batch;
            DAT.combat_contrasts.batch_var = '<supplied as vector>';
        end
        batch_raw = batch_raw(:);

        [batch_labels, ~, batch] = unique(cellstr(string(batch_raw)), 'stable');
        batch = double(batch);

        ref_code = [];
        if ~isempty(combat_ref_batch)
            ref_code = find(strcmp(batch_labels, char(string(combat_ref_batch))));
            if isempty(ref_code)
                error('combat_ref_batch ''%s'' is not among the batches present (%s).', ...
                    char(string(combat_ref_batch)), strjoin(batch_labels(:)', ', '));
            end
        end

        mod = [];
        for k = 1:numel(combat_mod)
            if ~istable(Tc) || ~ismember(combat_mod{k}, Tc.Properties.VariableNames)
                error('combat_mod variable ''%s'' is not in DAT.BETWEENPERSON.contrasts{%d}.', combat_mod{k}, c);
            end
            mod = [mod double(Tc.(combat_mod{k})(:))]; %#ok<AGROW>
        end

        fprintf('  batch variable : %s\n', DAT.combat_contrasts.batch_var);
        if isempty(combat_mod), fprintf('  preserved      : none\n');
        else,                   fprintf('  preserved      : %s\n', strjoin(combat_mod, ', ')); end
        if isempty(ref_code),   fprintf('  reference      : none (grand mean)\n');
        else,                   fprintf('  reference      : %s\n', batch_labels{ref_code}); end
        for b = 1:numel(batch_labels)
            fprintf('    %-12s n = %3d\n', batch_labels{b}, sum(batch == b));
        end

        % --- harmonize the unscaled and the scaled contrast objects ---------
        % Written out explicitly for both objects rather than looped over
        % variable NAMES: a loop would need eval/assignin, and assignin('base')
        % targets the base workspace, which is not the script's workspace when
        % the script is run under publish().
        % --- OPTIONAL per-site prescaling, BEFORE ComBat ---------------------
        % Recommended for fMRI by the ENIGMA-lineage implementation
        % (combat.enigma, Radua et al.): "prescaling is a good option for fMRI,
        % where different devices can have varying units of measurement", and is
        % "especially beneficial when the sites use different scales".
        %
        % ONE factor per site, deliberately. That is the whole difference from
        % per-image z-scoring: a per-SUBJECT normalization also removes genuine
        % between-subject differences in effect magnitude - part of the group
        % effect being tested - whereas a per-SITE factor rescales units while
        % leaving within-site between-subject variation intact.
        %
        % The factor is the MEDIAN over that site's subjects of each subject's
        % RMS across voxels, i.e. a robust measure of typical image magnitude at
        % that site, expressed relative to the reference site so the reference
        % passes through unchanged. Median rather than mean so one extreme
        % subject cannot set a whole site's scale. Note this choice of factor is
        % a judgement, not something the source prescribes: any monotone measure
        % of site magnitude would serve, and the arms should be compared rather
        % than one assumed correct.
        if combat_prescale_sites
            if isempty(ref_code)
                error(['combat_prescale_sites requires combat_ref_batch: without a reference ' ...
                       'site there is no scale to harmonize the other sites towards.']);
            end
            fprintf('\n  PER-SITE PRESCALING (before ComBat), reference %s\n', batch_labels{ref_code});
            for o = 1:2
                if o == 1
                    if ~exist('DATA_OBJ_CON','var') || numel(DATA_OBJ_CON) < c, continue; end
                    Xps = double(DATA_OBJ_CON{c}.dat);
                else
                    if ~exist('DATA_OBJ_CONsc','var') || numel(DATA_OBJ_CONsc) < c, continue; end
                    Xps = double(DATA_OBJ_CONsc{c}.dat);
                end
                subj_rms = sqrt(mean(Xps.^2, 1));                 % 1 x n, per subject
                site_mag = zeros(1, numel(batch_labels));
                for b = 1:numel(batch_labels)
                    site_mag(b) = median(subj_rms(batch == b));
                end
                fac = site_mag(ref_code) ./ site_mag;             % reference gets 1
                for b = 1:numel(batch_labels)
                    Xps(:, batch == b) = Xps(:, batch == b) * fac(b);
                end
                if o == 1
                    DATA_OBJ_CON{c}.dat = Xps;
                    for b = 1:numel(batch_labels)
                        fprintf('    CON    %-10s magnitude %8.4f  factor %6.3f\n', batch_labels{b}, site_mag(b), fac(b));
                    end
                else
                    DATA_OBJ_CONsc{c}.dat = Xps;
                    for b = 1:numel(batch_labels)
                        fprintf('    CONsc  %-10s magnitude %8.4f  factor %6.3f\n', batch_labels{b}, site_mag(b), fac(b));
                    end
                end
            end
            DAT.combat_contrasts.prescaled = true;
        else
            DAT.combat_contrasts.prescaled = false;
        end

        cb_common = {batch, mod, double(logical(combat_parametric))};
        if ~isempty(ref_code), cb_common = [cb_common {'ref', ref_code}]; end %#ok<AGROW>

        if exist('DATA_OBJ_CON','var') && numel(DATA_OBJ_CON) >= c && ~isempty(DATA_OBJ_CON{c})
            dat = DATA_OBJ_CON{c}.dat;
            if size(dat,2) ~= numel(batch)
                error('DATA_OBJ_CON{%d} has %d images but the batch vector has %d entries.', ...
                    c, size(dat,2), numel(batch));
            end
            ok = true(size(dat,1),1);
            for b = 1:numel(batch_labels)
                ok = ok & (std(double(dat(:, batch == b)), 0, 2) > 0);
            end
            fprintf('\n  DATA_OBJ_CON: harmonizing %d of %d voxels\n', sum(ok), numel(ok));
            dat(ok,:) = combat(double(dat(ok,:)), cb_common{:});
            DATA_OBJ_CON{c}.dat = dat;
            DATA_OBJ_CON{c} = enforce_variable_types(DATA_OBJ_CON{c});
        end

        if exist('DATA_OBJ_CONsc','var') && numel(DATA_OBJ_CONsc) >= c && ~isempty(DATA_OBJ_CONsc{c})
            dat = DATA_OBJ_CONsc{c}.dat;
            if size(dat,2) ~= numel(batch)
                error('DATA_OBJ_CONsc{%d} has %d images but the batch vector has %d entries.', ...
                    c, size(dat,2), numel(batch));
            end
            ok = true(size(dat,1),1);
            for b = 1:numel(batch_labels)
                ok = ok & (std(double(dat(:, batch == b)), 0, 2) > 0);
            end
            fprintf('  DATA_OBJ_CONsc: harmonizing %d of %d voxels\n', sum(ok), numel(ok));
            dat(ok,:) = combat(double(dat(ok,:)), cb_common{:});
            DATA_OBJ_CONsc{c}.dat = dat;
            DATA_OBJ_CONsc{c} = enforce_variable_types(DATA_OBJ_CONsc{c});
        end

        % DATA_OBJ_CONscc is DERIVED from DATA_OBJ_CON (l2norm), so it must be
        % rebuilt from the harmonized version or it would silently stay stale.
        if exist('DATA_OBJ_CONscc','var') && numel(DATA_OBJ_CONscc) >= c
            DATA_OBJ_CONscc{c} = rescale(DATA_OBJ_CON{c}, 'l2norm_images');
            DATA_OBJ_CONscc{c} = enforce_variable_types(DATA_OBJ_CONscc{c});
            fprintf('  DATA_OBJ_CONscc rebuilt by l2norm from the harmonized DATA_OBJ_CON\n');
        end

        drawnow; snapnow

    end

    DAT.combat_contrasts.applied = true;
    fprintf('\n\nComBat applied to contrast images\n\n');

else

    DAT.combat_contrasts.applied = false;

end


%% SAVE RESULTS
% ------------------------------------------------------------------------

fprintf('\n\n');
% contrast_objects_tag lets one model hold more than one set of contrast objects,
% which is what a second, differently-harmonised path needs: prep_3 can be run twice
% over the same conditions with different combat_mod and the two results no longer
% collide on one filename. Empty (the default) reproduces the original behaviour
% exactly, so every existing study is unaffected.
%
% DEFINED HERE, above the first use below - an earlier version defined it further
% down, beside savefilenamedata, and the printhdr line that also uses it ran first.
if ~exist('contrast_objects_tag','var') || isempty(contrast_objects_tag)
    contrast_objects_tag = '';
end
printhdr(['SAVE CONTRAST DATA OBJECTS IN contrast_data_objects' contrast_objects_tag '.mat']);
fprintf('\n\n');

savefilenamedata = fullfile(resultsdir, ['contrast_data_objects' contrast_objects_tag '.mat']);   % both unscaled and two versions of scaled
if isempty(contrast_objects_tag)
    save(savefilenamedata, 'DATA_OBJ_CON*', '-v7.3');
else
    % A tagged run is a SECOND harmonisation path over the same conditions. Its
    % provenance travels with its own objects rather than in the shared DAT,
    % which belongs to the untagged path.
    combat_record = DAT.combat_conditions; %#ok<NASGU>
    save(savefilenamedata, 'DATA_OBJ_CON*', 'combat_record', '-v7.3');
end                       % Note: 6/7/17 Tor switched to -v7.3 format by default 


%% GET CONTRASTS IN GLOBAL GRAY, WHITE, CSF VALUES
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('CALCULATE CONTRASTS IN GRAY/WHITE/CSF VALUES');
fprintf('\n\n');

DAT.gray_white_csf_contrasts = {};

for c = 1:size(DAT.contrasts, 1)
    
    wh = find(DAT.contrasts(c, :));
    
    DAT.gray_white_csf_contrasts{c} = zeros(size(DAT.gray_white_csf{wh(1)}));
    
    for i = 1:k
        
        if DAT.contrasts(c, i) == 0
            % Skip.  This allows us to include image sets with different
            % numbers of images, as occurs with between-person designs, as
            % long as the contrast weights are zero for all elements with
            % different numbers of images.
            continue
        end
        
        % add data * contrast weight
        DAT.gray_white_csf_contrasts{c} = DAT.gray_white_csf_contrasts{c} + DAT.gray_white_csf{i} .* DAT.contrasts(c, i);
        
    end
    
end


%% ADD TO PREVIOUSLY SAVED RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('ADDED CONTRAST GRAY/WHITE/CSF TO DAT in image_names_and_setup.mat');
fprintf('\n\n');

cd(resultsdir); % unannex image_names_and_setup.mat file if already datalad saved to prevent write permission problems
! git annex unannex image_names_and_setup.mat
cd(rootdir);

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');

% GUARD: does the stored DAT carry signature results this one would discard?
%
% prep_4 appends DAT.SIG_conditions / DAT.SIG_contrasts to this same file. The
% '-append' below replaces the whole DAT variable, so re-running this script
% after prep_4 destroys them silently - no error, nothing in the report. That is
% exactly how proj_moodbugs model_1b lost its signature analysis: prep_4 ran on
% 2026-09-08, this script was re-run on 2026-09-09, and the SIG fields have been
% absent ever since.
%
% They are deliberately NOT carried over: the image objects have just been
% rebuilt, so any signature response computed from the previous ones is stale,
% and keeping it quietly would be worse than losing it. Re-run prep_4.
if exist(savefilename, 'file')
    stored = whos('-file', savefilename);
    if ismember('DAT', {stored.name})
        prior = load(savefilename, 'DAT');
        lost = intersect({'SIG_conditions','SIG_contrasts'}, fieldnames(prior.DAT));
        lost = lost(~ismember(lost, fieldnames(DAT)));
        if ~isempty(lost)
            warning(['\n\n*** SIGNATURE RESULTS OVERWRITTEN ***\n' ...
                     'The stored DAT carried %s, which this save discards.\n' ...
                     'They were computed from the image objects that have just been rebuilt,\n' ...
                     'so they are stale and are NOT being kept.\n' ...
                     'RE-RUN prep_4_apply_signatures_and_save (and any h_/d_ script after it)\n' ...
                     'before reporting any signature result from this model.\n\n'], ...
                     strjoin(lost, ' and '));
        end
    end
end

% Only the UNTAGGED path owns the shared DAT. A tagged run rebuilds the same
% contrasts under a different harmonisation, so appending its DAT here would
% overwrite combat_conditions and gray_white_csf_contrasts with values that do
% not describe the objects the GLM actually uses - a silent metadata/contrast
% mismatch. Its own record is saved beside its contrast objects instead.
if isempty(contrast_objects_tag)
    save(savefilename, '-append', 'DAT');
else
    fprintf(['\nDAT NOT appended to %s: this is a tagged run (%s), and the shared\n' ...
             'DAT belongs to the untagged harmonisation path. The harmonisation record\n' ...
             'for this path is saved as combat_record in %s.\n'], ...
            savefilename, contrast_objects_tag, savefilenamedata);
end
