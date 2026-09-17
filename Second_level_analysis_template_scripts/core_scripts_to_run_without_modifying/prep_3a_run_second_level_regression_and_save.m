%% prep_3a_run_second_level_regression_and_save.m
%
%
% *USAGE*
%
% This script
%
% 1. Runs second-level (i.e. across subjects) regression analyses
%   for each within-subject CONTRAST or CONDITION registered in the DAT
%   structure, either
%
%       * voxel-wise, calling CANlab's regress() function under the hood,
%       including its robust regression option if specified
%
%       * parcel-wise, calling CANlab's robfit_parcelwise() function under the
%       hood, which is robust by default
%
% 2. The option to convert t-maps into BayesFactor maps using CANlab's
%   estimateBayesFactor() function is built in - see walkthrough 
%   https://canlab.github.io/_pages/EmoReg_BayesFactor_walkthrough/EmoReg_BayesFactor_walkthrough.html
%
% 3. Runs cross-validated MVPA regression models predicting continuous
%   covariates if desired using CANlab's predict() function
%
% 4. Saves the results using standard naming and location
% 
% Run this script headless from the Linux command line (the default), which
% publishes the html report and fails loudly if the script errors:
%
%   labgascore_run_headless.sh -d /data/proj_xxx \
%       -s <proj>_secondlevel_m<M>_s0_a_set_up_paths_always_run_first \
%       <proj>_secondlevel_m<M>_s<N>_prep_3a_run_second_level_regression_and_save
%
% Or, interactively from the Matlab terminal (use this when you want
% higher-resolution figures, or are debugging):
%
%   LaBGAScore_prov_publish('prep_3a_run_second_level_regression_and_save', htmlsavedir)
%
% NOTE: publish() catches a script error into the html and returns normally, so
% a crashed run looks exactly like a successful one. Prefer the routes above,
% which read the report back and check for a caught error, over a bare
% publish('prep_3a_run_second_level_regression_and_save','outputDir',htmlsavedir).
%
% To get results reports after thresholding, publish
% c2a_second_level_regression
%
%
%
% *SETTING results_suffix FROM A CALLING SCRIPT - READ THIS FIRST*
%
% This script re-declares results_suffix = '' in its own option block below.
% An override placed BEFORE that point - for instance in the study's
% a_set_up_paths or a2_set_default_options - is therefore silently wiped, and
% the results are saved with an empty suffix, overwriting whatever a previous
% model wrote to the same filename. Set results_suffix AFTER the script's own
% default, not before it. Learned the hard way in proj_cfs model_2a.
%
% *OPTIONS*
%
% * NOTE 
%       defaults are specified in a2_set_default_options for any given model,
%       but if you want to run the same model with different options (for example
%       voxel- and parcelwise regression), you can make a copy of this script with
%       a letter index (e.g. _s6a_) and change the default option here
%
% * dorobust                    robust regression or OLS (true/false)
%
% * dorobfit_parcelwise         voxel- or parcelwise regression (true/false)
%                              NOTE ON TFCE: the TFCE block sits inside
%                              "if ~dorobfit_parcelwise" BY DESIGN. TFCE is a voxel-level method -
%                              it integrates over cluster-forming thresholds using spatial
%                              contiguity between neighbouring voxels - and a parcelwise fit has
%                              discrete parcels rather than a spatial field, so there is no cluster
%                              extent for it to operate on. doTFCE is therefore ignored in a
%                              parcelwise run, correctly. Set it false explicitly anyway, so the
%                              option block states what actually happens.
%
%       * csf_wm_covs               true adds global wm & csf regressors at second level
%       * remove_outliers           true removes outlier images/subjects based on mahalanobis distance
%
% * myscaling_glm               'raw', 'scaled', or 'scaled_contrasts' (defined in a2_set_..., image scaling done in prep_2_... and prep_3_... data load)
%
% * maskname_glm
%
%       * default use of sparse gray matter mask
%       * model-specific maskdir defined in a_set_up_paths_always_run_first script
%       * if you do not want to mask, change to []
%       * if you want to use a custom mask, put it in maskdir and change name here
%       * only used for visualization of uncorrected results in this script
%
% * atlasname_glm               atlas object used for 
%                                   1. defining parcels and masking in parcelwise analysis
%                                       in this case make sure it is an atlas object corresponding to the binary mask defined in maskname_glm
%                                   2. labeling regions in both voxelwise and parcelwise analyses
%
%                                   option a - atlas name from load_atlas.m for different atlas than
%                                       default canlab_2018 (used if you do not specify this option)
%                                   option b - which('xxx.mat') name of .mat file in maskdir, generated by
%                                       https://github.com/labgas/LaBGAScore/blob/main/atlas_mask_tools/LaBGAScore_atlas_binary_mask_from_atlas.m
%                                       this file should contain a 'combined_atlas' var containing an atlas object
%                                       USAGE: use option b if you want to restrict the parcels to the
%                                               ones included in maskname_glm, or any atlas or a subset thereof
%
%       * atlas_granularity         level of granularity of canlab2023 atlas for 
%                                       1. defining parcels in parcelwise analysis
%                                       2. labeling regions in both voxel- and parcelwise analysis
%
%                                       options: 1 = fine (595 parcels), 2 = intermediate (525 parcels), 3 = coarse (264 parcels)
%
% * design_matrix_type
%
%       1. 'group' 
%           Assuming that groups are concatenated in contrast image lists, and
%           regressor values of 1 or -1 will specify the group identity for each image. 
%           Requires DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.(mygroupfieldname){c}.groupfield specifying group membership for
%           each image.
%
%       2. 'custom'
%           Uses all columns of table object DAT.BETWEENPERSON.(mygroupnamefield){c}
%           NOTE: you can flexibly use one or more of these columns as
%                                       covariates by specifying the covs2use option below
%           Can enter a multi-column design matrix for each contrast
%           Design matrix can be different for each contrast
%
%       3. 'onesample' option:
%           Only adds intercept, hence performs a one-sample t-test on contrast
%           images across all subjects, similarly to c_univariate_contrast_maps_
%           scripts, but with more flexible options including scaling and robustfit
%
%       NOTE: To set up group and custom variables, see prep_1b_prep_behavioral_data
%
% * doBayes                     convert t-maps into Bayes Factors 
%
% * doTFCE                      calculate TFCE maps from fmri_data_object
%
%       _TFCE analysis options_
%
%         TFCE here is classic threshold-free cluster enhancement (Smith &
%         Nichols 2009), computed by LaBGAScore's group_tfce_from_subject_maps
%         with a sign-flip (one-sample) or label-exchange (two-sample)
%         permutation null. Height and extent exponents and connectivity are
%         left at their defaults (H = 2, E = 0.5, conn = 26); pass them through
%         group_tfce_from_subject_maps directly if you need to change them.
%     cons2tfce:
%         vector of contrast indices to run TFCE on, if you only want it for a
%         subset. Empty (default) runs TFCE on every contrast. TFCE is by far
%         the most expensive step in this script - perm_n_tfce permutations per
%         contrast - so restricting it to the contrast(s) of interest is often
%         the difference between an overnight job and a coffee break.
%
%         NOTE: results produced before the 2026 TFCE overhaul of LaBGAScore
%         are not comparable. That work replaced a mis-parameterised pTFCE call
%         and repaired both permutation schemes; the two-sample null had been
%         degenerate, and one-sample with covariates could never reject.
%
%
%         * perm_n_tfce         number of permutations for TFCE-based stats
%         * tfce_sidedness      'one' versus 'two'-tailed test for TFCE-based stats
%         * tfce_tail           'pos' or 'neg' if tfce_sidedness = 'one'
%         * tfce_seed           base RNG seed for the permutation null. Optional:
%                               if unset a seed is drawn and printed, and is also
%                               returned in tfce_info.seed. Set it to reproduce a
%                               previous run exactly.
%
% * doroi_analysis              extract roi averages from condition (beta) or contrast (con) images using atlas objects created by LaBGAScore_atlas_binary_mask_from_atlas.m and written in secondlevel/modeldir/masks as input
%
%       PREREQUISITE: the roi masks must already exist. They are NOT created
%       here. Generate them first with LaBGAScore_atlas_rois_from_atlas.m
%       (LaBGAScore atlas_mask_tools/), run from the root of your superdataset,
%       which writes them into the model's maskdir as
%
%           <maskdir>/<roi_modelname>_rois_<roi_set_name>.mat
%
%       This script loads exactly that path. A missing file used to surface as a
%       bare "Unable to read file" from load(), partway through a long run; it
%       is now checked before the regression starts.
%
%       NOTE roi_modelname is only a FILENAME PREFIX. The file is loaded from
%       the maskdir of the model you are running now, whatever roi_modelname
%       says. To reuse an roi set generated for another model, COPY its .mat
%       into this model's maskdir - pointing roi_modelname at that model is not
%       enough, and is an easy way to hit the error above.
%
%       _roi analysis options_
%
%         * roi_names            cell array of names corresponding to roiname variables in LaBGAScore_atlas_rois_from_atlas.m which writes atlas objects for each roi in secondlevel/modeldir/masks
%         * roi_modelname        from same script; also the prefix of the .mat filename
%         * roi_set_name         from same script; also part of the .mat filename
%
%         * doroi_glm            true runs inference on the roi averages, false (default) keeps
%                                the earlier behaviour, where barplot_columns plots covariate-
%                                adjusted roi means but no test is run, so an roi effect could
%                                not be called significant without refitting by hand.
%
%                                Which covariates are of interest and which are nuisance comes
%                                from nuisance_covs - NOT from a new option, and NOT from
%                                covs2use, which SUBSETS the design matrix (dropping every
%                                covariate not listed) rather than labelling its columns.
%                                Everything in the design that nuisance_covs does not name is an
%                                effect of interest.
%
%                                Two levels are reported, in this order:
%
%                                  1. MANOVA across the whole roi set (roi_manova_stats), one
%                                     omnibus test per effect of interest. Wilks' Lambda from the
%                                     full model against a reduced model without that effect, so
%                                     nuisance covariates are adjusted for, converted to Rao's F.
%                                     manova1 cannot do this: it is one-way and takes no
%                                     covariates. Skipped with a message if there are fewer than
%                                     2 rois, or if the error df do not exceed the number of rois
%                                     (the residual covariance would be singular).
%                                  2. A GLM per roi (roi_glm_stats): one model per roi holding
%                                     every covariate in the design, no interaction,
%                                         roi_mean ~ effect(s) of interest + nuisance...
%                                     Each effect is reported separately with BOTH q_BH and
%                                     q_Storey. Storey needs many tests to estimate pi0 and falls
%                                     back to BH when it cannot (an roi set is usually far too
%                                     small for it), so with a handful of rois expect the two
%                                     columns to agree. Effect size is Cohen's d for a two-level
%                                     effect, otherwise the partial correlation, and 'estimate'
%                                     for a two-level effect is the DIFFERENCE between levels,
%                                     not the raw beta.
%
%                                Both tables are printed into the published report and saved in
%                                roi_stats_*.mat. Read them together: individual rois surviving
%                                FDR under a null omnibus test should be treated cautiously.
%                                Ignored when design_matrix_type is 'onesample' (no covariates).
%
% * doneurotransmitter_maps     calculate spatial similarity with neurotransmitter maps from Hansen et al Nat Neurosci 2022 for each contrast/condition
%
%       _neurotransmitter map options_
%
%         * neurotransmitter_maps_metric      'correlation' (default) or 'cosine_similarity'
%
% * domvpa_reg_cov              run MVPA regression model to predict covariate levels from (between-subject) brain data using CANlab's predict() function
%
%       NOTE: THIS OPTION ONLY APPLIES WHEN DESIGN_MATRIX_TYPE = 'CUSTOM' SINCE OTHERWISE THERE IS NO CONTINUOUS OUTCOME TO PREDICT!
%         TO CLASSIFY GROUPS USING MVPA MODELS, USE SVM SCRIPTS PREP_3C AND C2
%     
%     _mvpa_reg_covariate options_
%
%       * algorithm_mvpa_reg_cov                e.g. 'cv_pcr', or other option passed into predict function (help fmri_data.predict for options)
%
%       * holdout_set_method_mvpa_reg_cov
%
%           1. group: use DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.(mygroupfieldname){c}.group to balance holdout sets over groups
%                                        
%
%           2. no_group: no group factor, stratifies by subject (i.e.leave whole subject out) since data is purely between-subject
%
%       * nfolds_mvpa_reg_cov                   number of cross-validation folds for kfold
%
%       * zscore_outcome_mvpa_reg_cov           zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models
%
%
% *MANDATORY OPTIONS TO BE SPECIFIED IN THIS SCRIPT*
%
% * mygroupfieldname            'contrasts' or 'conditions'
%
% * results_suffix              name to add to results file to specify in case of multiple versions of model, e.g. 'covariate_rating'
%
% 
% *OPTIONS TO BE SPECIFIED IN THIS SCRIPT IF DESIGN_MATRIX_TYPE = CUSTOM*
%
% * covs2use                    variable name(s) in DAT.BETWEENPERSON.(mygroupnamefield){:} to be used as covariates in GLM and, if domvpa_reg_cov = true, outcome in MVPA regression
% * nuisance_covs               variable name(s) in DAT.BETWEENPERSON.(mygroupnamefield){:} to be used as nuisance covariate rather than covariate of interest in GLM
%
%       NOTE: only use the first option if you don't want to use all variables in the above table as covariate, otherwise delete or comment out below
%
%
% *OPTIONS TO BE SPECIFIED IN THIS SCRIPT IF DESIGN_MATRIX_TYPE = GROUP*
%
% * group_id                    name of group identifier variable in same table
%
%       NOTE: not needed if DAT.BETWEENPERSON.group contains group identifier, in that case delete or comment out
%
%
% -------------------------------------------------------------------------
%
% revamped by: Lukas Van Oudenhove
%
% date:   Dartmouth, May, 2022
%
% -------------------------------------------------------------------------
%
% prep_3a_run_second_level_regression_and_save.m         v9.2
%
% last modified: 2026/08/14
%
%
%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS
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

% NOTES 
%   1. CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
%   2. THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% s0 MUST RUN BEFORE THE OPTION BLOCK BELOW. It calls a2_set_default_options,
% so every option a2 defines is (re)assigned at this point - anything set above
% this line is silently discarded. Placing it first means the script-specific
% options below reliably override the a2 defaults. c2a already orders it this
% way; prep_3a did not, which silently reverted dorobfit_parcelwise in the
% parcelwise variants and made them run voxelwise.

% SET MANDATORY OPTIONS

mygroupnamefield = 'contrasts'; 
results_suffix = ''; % adds a suffix of your choice to .mat file with results that will be saved

% NOTES 
%   1. do NOT delete the latter option, leave empty if not needed
%   2. do NOT use to add a suffix specifying the regressors, scaling or masking option, this will be added automatically


% OPTIONS IF DESIGN_MATRIX_TYPE = CUSTOM

% covs2use = {'delta_wanting'};      % needs to correspond to variable name(s) in DAT.BETWEENPERSON.(mygroupnamefield){:} AND THE ORDER IN WHICH THEY APPEAR THERE
% nuisance_covs = {'center'};

% NOTE: if you want to use all variables in DAT.BETWEENPERSON.(mygroupnamefield){:} as covariates, comment this option out


% OPTIONS IF DESIGN_MATRIX_TYPE = GROUP

% group_id = {'group'};             % needs to correspond to variable name(s) in DAT.BETWEENPERSON.(mygroupnamefield){:} AND THE ORDER IN WHICH THEY APPEAR THERE

% NOTE: if DAT.BETWEENPERSON.group contains group identifier, you can comment this option out


% GET DEFAULT OPTIONS IF NOT SET IN A2_SET_DEFAULT_OPTIONS

options_needed = {'dorobust', 'dorobfit_parcelwise', 'myscaling_glm', 'design_matrix_type', 'maskname_glm'};
options_exist = cellfun(@exist, options_needed); 

option_default_values = {false, false, 'raw', 'onesample', which('gm_mask_canlab2023_coarse_fmriprep20_0_20.nii')};

plugin_get_options_for_analysis_script;


% SET CUSTOM OPTIONS

% NOTE
%   only specify if you want to run multiple versions of your model with different options
%   than the defaults you set in your model-specific version of a2_set_default_options.m

% maskname_glm = 'mask_name';
% atlasname_glm = 'atlas_name';
%   atlas_granularity = 1/2/3;
% dorobust = true/false;
% dorobfit_parcelwise = true/false;
%   csf_wm_covs = true/false;
%   remove_outliers = true/false;
% myscaling_glm = 'raw'/'scaled'/'scaled_contrasts';
% design_matrix_type = 'custom';
% doBayes = true/false;
% doTFCE = true/false;
%     perm_n_tfce = [number];                                                    
%     tfce_sidedness = 'two'/'one';                                                
%     tfce_tail = 'pos/neg';
%     tfce_seed = [integer];                                               
% doroi_analysis = true/false;
%   roi_names = {'x','y','z'};
%   roi_modelname = 'modelname';
%   roi_set_name = 'setname';
% doneurotransmitter_maps = true/false;
%   neurotransmitter_maps_metric = 'cosine_similarity'/'correlation';
% domvpa_reg_cov = true/false;
%   algorithm_mvpa_reg_cov = 'cv_pcr';
%   holdout_set_method_mvpa_reg_cov = 'no_group'/'group';
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


%% RESTRICT THE SAMPLE (OPTIONAL)
% -------------------------------------------------------------------------
%
% subject_filter subsets the analysis to a subgroup, so a study can run the
% same model on nested samples - "tiers" - without building a separate model
% directory for each. It is the GLM counterpart of the option of the same name
% in LaBGAScore_decoding_SVM_between_subjects, so the two analyses can be run
% on identical samples and compared.
%
% Format: a cell of {column, values-to-keep} pairs, ANDed, naming columns of
% DAT.BETWEENPERSON.(mygroupnamefield){:}:
%
%   subject_filter = { {'center_UM', 0}, {'center_UGOT', 0} };   % KUL only
%   subject_filter = { {'center_UM', 0} };                       % drop UM
%
% Everything is subset together - the image objects, the covariate tables and
% the group vector - so the design cannot fall out of step with the data.
% ALWAYS set results_suffix as well: without it each tier overwrites the last.
if exist('subject_filter','var') && ~isempty(subject_filter)

    % ALWAYS reload from disk before filtering.
    %
    % publish() runs a batch of scripts in ONE workspace, so objects left by an
    % earlier script are still resident. prep_3a's loader is guarded by
    % "if ~exist('DATA_OBJ_CON','var')", which means a second filtered run in
    % the same batch would inherit the FIRST run's already-subset data and
    % filter it again - "keeping 101 of 101" instead of 124 of 158, with no
    % error and a plausible-looking report. Reloading makes each filtered run
    % independent of whatever ran before it.
    clear DATA_OBJ_CON DATA_OBJ_CONsc DATA_OBJ_CONscc DATA_OBJ DATA_OBJsc
    load(fullfile(resultsdir,'contrast_data_objects.mat'));
    tmp_sf_dat = load(fullfile(resultsdir,'image_names_and_setup.mat'),'DAT');
    DAT = tmp_sf_dat.DAT; clear tmp_sf_dat
    fprintf('\nreloaded unfiltered data and DAT before applying subject_filter\n');

    keep_sf = [];
    for sf = 1:numel(subject_filter)
        col = subject_filter{sf}{1};
        val = subject_filter{sf}{2};
        tbl_sf = DAT.BETWEENPERSON.(mygroupnamefield){1};
        if ~ismember(col, tbl_sf.Properties.VariableNames)
            error(['\nsubject_filter names ''%s'', which is not a column of ' ...
                   'DAT.BETWEENPERSON.%s{1}.\nAvailable: %s\n'], ...
                   col, mygroupnamefield, strjoin(tbl_sf.Properties.VariableNames, ', '));
        end
        k = ismember(double(tbl_sf.(col)), val);
        if isempty(keep_sf), keep_sf = k; else, keep_sf = keep_sf & k; end
    end

    if ~any(keep_sf)
        error('\nsubject_filter left 0 subjects.\n');
    end

    fprintf('\n=== SUBJECT FILTER: keeping %d of %d subjects ===\n', sum(keep_sf), numel(keep_sf));
    for sf = 1:numel(subject_filter)
        fprintf('    %s in %s\n', subject_filter{sf}{1}, mat2str(subject_filter{sf}{2}));
    end

    idx_sf = find(keep_sf);

    for cc_sf = 1:numel(DAT.BETWEENPERSON.conditions)
        if ~isempty(DAT.BETWEENPERSON.conditions{cc_sf})
            DAT.BETWEENPERSON.conditions{cc_sf} = DAT.BETWEENPERSON.conditions{cc_sf}(keep_sf,:);
        end
    end
    for cc_sf = 1:numel(DAT.BETWEENPERSON.contrasts)
        if ~isempty(DAT.BETWEENPERSON.contrasts{cc_sf})
            DAT.BETWEENPERSON.contrasts{cc_sf} = DAT.BETWEENPERSON.contrasts{cc_sf}(keep_sf,:);
        end
    end
    if isfield(DAT.BETWEENPERSON,'group') && ~isempty(DAT.BETWEENPERSON.group)
        DAT.BETWEENPERSON.group = DAT.BETWEENPERSON.group(keep_sf);
    end

    for objname_sf = {'DATA_OBJ_CON','DATA_OBJ_CONsc','DATA_OBJ_CONscc','DATA_OBJ','DATA_OBJsc'}
        if exist(objname_sf{1},'var')
            tmp_sf = eval(objname_sf{1});
            if iscell(tmp_sf)
                for cc_sf = 1:numel(tmp_sf)
                    if ~isempty(tmp_sf{cc_sf})
                        tmp_sf{cc_sf} = get_wh_image(tmp_sf{cc_sf}, idx_sf);
                    end
                end
                eval([objname_sf{1} ' = tmp_sf;']);
            end
        end
    end
    clear tmp_sf idx_sf keep_sf k col val tbl_sf cc_sf sf objname_sf

    fprintf('    image objects, covariate tables and group vector all subset together\n\n');

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

if ~dorobfit_parcelwise

    if exist('maskname_glm', 'var') && ~isempty(maskname_glm) && exist(maskname_glm, 'file')
        
        [~,maskname_short] = fileparts(maskname_glm);
            if contains(maskname_short,'nii')
                [~,maskname_short] = fileparts(maskname_short);
            end
        mask_string = sprintf('masked with %s', maskname_short);
        glmmask = fmri_mask_image(maskname_glm, 'noverbose'); 
           
             if any(unique(glmmask.dat) ~= 1)
                glmmask.dat(glmmask.dat > 0) = 1; % binarize mask if needed
             end
            
        fprintf('\nMasking voxelwise results visualization with %s\n\n', maskname_short);
        
        % MONTAGE OF MASK
        
        % canlab_results_fmridisplay's 'compact' layout only calls axes('Position',...);
        % unlike 'multirow' it never opens a figure of its own, so it draws into whatever
        % figure is current - the previous block's montage, or the last figure left open by
        % the previous script in the same session. Open a fresh one. ('multirow' does
        % create its own figure, so those call sites are deliberately left alone: adding
        % figure; there would leave an empty figure behind for every montage.)
        figure;
        o2 = canlab_results_fmridisplay([], 'compact');
        o2 = addblobs(o2, glmmask,'onecolor','color',[0.4 0.2 0.6],'trans','transvalue',0.50);
        o2 = title_montage(o2, 5, ['voxel-wise analysis masked with: ' maskname_short]);
        plugin_set_figure_size;
        drawnow,snapnow;
        
        clear o2
        
    else
        
        mask_string = sprintf('without masking');
        fprintf('\nShowing voxelwise results without masking\n\n');
        
    end 

end
    
if exist('atlasname_glm','var') && ~isempty(atlasname_glm)

    if contains(atlasname_glm,'.mat')
        
        load(atlasname_glm);
        [~,atlasname_short] = fileparts(atlasname_glm);
        clear mask
        
        if logical(exist('atlas_granularity','var')) && atlas_granularity ~= 1
            combined_atlas = combined_atlas.downsample_parcellation(['labels_' num2str(atlas_granularity)]);
        end
        
        if dorobfit_parcelwise
            
            maskname_short = atlasname_short;
            mask_string = sprintf('masked with %s', maskname_short);
            fprintf('\nRunning parcelwise analysis in custom-made atlas %s at granularity level labels_%d\n\n', atlasname_short, atlas_granularity);
            
        end
        
        fprintf('\nLabeling regions using custom-made atlas %s at granularity level labels_%d\n\n', atlasname_short, atlas_granularity);
        
        % MONTAGE OF ATLAS
        
        cmap = colormap('colorcube');
        close gcf;
                
        figure;
        o2 = canlab_results_fmridisplay([], 'compact');
        o2 = addblobs(o2, atlas2region(combined_atlas),'indexmap',cmap,'interp','nearest');
        if dorobfit_parcelwise
            o2 = title_montage(o2, 5, ['parcel-wise analysis in atlas: ' atlasname_short ' at granularity level labels_' num2str(atlas_granularity)]);
        else
            o2 = title_montage(o2, 5, ['voxel-wise analysis labeled with atlas: ' atlasname_short ' at granularity level labels_' num2str(atlas_granularity)]);
        end
        plugin_set_figure_size;
        drawnow,snapnow;
        
        clear o2

    elseif ischar(atlasname_glm)

        combined_atlas = load_atlas(atlasname_glm);
        
        if logical(exist('atlas_granularity','var')) && atlas_granularity ~= 1
            combined_atlas = combined_atlas.downsample_parcellation(['labels_' num2str(atlas_granularity)]);
        end
        
        if contains(atlasname_glm,'canlab2023') || contains(atlasname_glm,'canlab2024')
            combined_atlas = combined_atlas.threshold(0.20); % only keep probability values > 0.20 in probabistic canlab2023/4 atlas
        end
        
        if dorobfit_parcelwise
            mask_string = sprintf('in atlas %s',atlasname_glm);
            
            fprintf('\nRunning parcelwise analysis in custom-made atlas %s at granularity level labels_%d\n\n', atlasname_glm, atlas_granularity);
            
        end
        
        fprintf('\nLabeling regions using custom-made atlas %s at granularity level labels_%d\n\n', atlasname_glm, atlas_granularity);
        
        
        % MONTAGE OF ATLAS
        
        cmap = colormap('colorcube');
        close gcf;
        
        figure;
        o2 = canlab_results_fmridisplay([], 'compact');
        o2 = addblobs(o2, atlas2region(combined_atlas),'indexmap',cmap,'interp','nearest');
        if dorobfit_parcelwise
            o2 = title_montage(o2, 5, ['parcel-wise analysis in atlas: ' atlasname_glm ' at granularity level labels_' num2str(atlas_granularity)]);
        else
            o2 = title_montage(o2, 5, ['voxel-wise analysis labeled with atlas: ' atlasname_glm ' at granularity level labels_' num2str(atlas_granularity)]);
        end
        plugin_set_figure_size;
        drawnow,snapnow;
        
        clear o2
        
    else

         error('\ninvalid option "%s" defined in atlasname_glm variable, should be a keyword for load_atlas.m or a .mat file containing an atlas object, check docs"\n\n',atlasname_glm)

    end
    
else

    if dorobfit_parcelwise
        mask_string = sprintf('without masking');
        
        fprintf('\nShowing parcelwise results without masking in 489 parcels of canlab_2018 atlas, which is now deprecated\n\n');
        
    end

end

brainmask = fmri_mask_image(maskname_brain,'noverbose');


%% MERGE ROI ATLAS OBJECTS
% -------------------------------------------------------------------------

if doroi_analysis
    
    % Fail with something actionable rather than a bare load() error. These
    % masks are generated by LaBGAScore_atlas_rois_from_atlas, not here, so a
    % model whose maskdir was never populated otherwise dies partway through a
    % long run on "Unable to read file".
    roifile = fullfile(maskdir,[roi_modelname '_rois_' roi_set_name '.mat']);
    if exist(roifile,'file') ~= 2
        error(['\nroi analysis requested (doroi_analysis = true) but the roi masks are missing:\n  %s\n\n' ...
               'Generate them with LaBGAScore_atlas_rois_from_atlas.m (LaBGAScore atlas_mask_tools/),\n' ...
               'run from the root of your superdataset with the SAME roi_modelname (''%s'') and\n' ...
               'roi_set_name (''%s''), which writes them to this model''s maskdir.\n\n' ...
               'If the set already exists for another model, COPY the .mat into\n  %s\n' ...
               'roi_modelname is only a filename prefix - the file is always loaded from the\n' ...
               'maskdir of the model being run.\n'], ...
               roifile, roi_modelname, roi_set_name, maskdir);
    end
    load(roifile);
    
    if logical(exist('roi_names','var')) && ~isempty(roi_names)
        
        all_rois = false;
    
        roi_idx = zeros(1,size(roi_atlases_flat,2));
        roi_idx = logical(roi_idx);

        for r = 1:size(roi_atlases_flat,2)

            roi_idx(r) = contains(roi_atlases_flat{r}.atlas_name,roi_names);

        end

        clear r

        % GUARD: roi_names must actually match something.
        %
        % a2 ships EXAMPLE roi_names (the bit_rew reward set). A study that turns
        % doroi_analysis on without replacing them selects nothing here, and
        % roi_atlases_flat silently becomes empty - the roi analysis then either
        % dies much later with an opaque error or quietly produces nothing. Fail
        % here instead, naming what is actually available.
        if ~any(roi_idx)
            avail = cellfun(@(x) x.atlas_name, roi_atlases_flat, 'UniformOutput', false);
            error(['\nNone of the roi_names matched the roi set in this model''s maskdir.\n\n' ...
                   '  roi_names  : %s\n  available  : %s\n\n' ...
                   'a2_set_default_options ships EXAMPLE roi_names - replace them with the\n' ...
                   'rois of the set you generated for this model.\n'], ...
                   strjoin(roi_names, ', '), strjoin(avail, ', '));
        end
        if sum(roi_idx) < numel(roi_names)
            avail = cellfun(@(x) x.atlas_name, roi_atlases_flat, 'UniformOutput', false);
            warning('CANlab:prep_3a:roiNamesPartial', ...
                ['only %d of %d roi_names matched. Unmatched names are ignored, so the roi\n' ...
                 'analysis will silently cover fewer rois than intended.\n  available: %s'], ...
                 sum(roi_idx), numel(roi_names), strjoin(avail, ', '));
        end

        roi_atlases_flat = roi_atlases_flat(roi_idx);
        
    else
        
        all_rois = true;
        
        roi_names = cell(1,size(roi_atlases_flat,2));
        
    end
        
    roi_atlas = roi_atlases_flat{1};
    roi_atlas.atlas_name = 'combined_roi_atlas';
    roi_atlas.labels{1} = roi_atlases_flat{1}.atlas_name;
    
    if all_rois
        
        roi_names{1} = roi_atlases_flat{1}.atlas_name;
        
    end
    
    roi = 2;
    
    while roi < (size(roi_atlases_flat,2) + 1)
        
        roi_atlas = merge_atlases(roi_atlas, roi_atlases_flat{roi},'noreplace');
%           LVO: whether or not 'noreplace' option is chosen or not should
%           not matter if all regions are from same atlas

        roi_atlas.labels{roi} = roi_atlases_flat{roi}.atlas_name;

            if all_rois

                roi_names{roi} = roi_atlases_flat{roi}.atlas_name;

            end
        
        roi = roi+1;
        
    end
    
    clear roi
    
    for roi = 1:size(roi_atlases_flat,2)
        roi_atlas.labels{roi} = roi_atlases_flat{roi}.atlas_name;
    end
    
    % MONTAGE OF ROI ATLAS OBJECT
        
        cmap2cell = scn_standard_colors(size(roi_atlas.labels,2));
        cmap2 = zeros(size(roi_atlas.labels,2),3);
        
        for color = 1:size(cmap2cell,2)
            cmap2(color,:) = cmap2cell{color};
        end
    
        figure;
        o2 = canlab_results_fmridisplay([], 'compact');
        o2 = addblobs(o2, atlas2region(roi_atlas),'indexmap',cmap2,'interp','nearest');
        o2 = title_montage(o2, 5, 'atlas used for extraction of roi averages');
        plugin_set_figure_size;
        drawnow,snapnow;
        
        clear o2
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
    
    if doBayes
        bayesian_regression_stats_results = cell(1, kc);
    end
    
    if doTFCE
        tfce_regression_stats_results = cell(1, kc);
    end
    
else
    
    parcelwise_stats_results = cell(1,kc);
    
end

if doroi_analysis
    
    roi_means = cell(1,kc);
    roi_means_table = cell(1,kc);
    roi_adjusted_means = cell(1,kc);
    roi_glm_stats = cell(1,kc);       % per-roi GLM, filled when doroi_glm is true
    roi_manova_stats = cell(1,kc);    % omnibus MANOVA across the roi set, same condition
    
end

if doneurotransmitter_maps
    
    neurotransmitter_stats = cell(1,kc);
    
    if isequal(design_matrix_type,'group') || (isequal(design_matrix_type,'custom') && ~isempty(DAT.BETWEENPERSON.group))
        
        neurotransmitter_group_stats = cell(1,kc);
        neurotransmitter_group_tables = cell(1,kc);
        neurotransmitter_multcomp_group = cell(1,kc);
        
    end
    
end


for c = 1:kc
    
    %%
    % *PREP WORK*
    
    % GET DESIGN MATRIX FOR THIS CONTRAST OR CONDITION
    % ------------------------------------------------
    
    switch mygroupnamefield
        
        case 'contrasts'
            fprintf('\n\n');
            printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.contrastnames{c})]);
            fprintf('\n\n');
            
        case 'conditions'
            fprintf('\n\n');
            printhdr(['CONDITION #', num2str(c), ': ', upper(DAT.conditions{c})]);
            fprintf('\n\n');
    
    end
      
    fprintf('\n\n');
    printhdr('BUILDING DESIGN MATRIX');
    fprintf('\n\n');
    
    groupnames_string = 'intercept';
    
    switch design_matrix_type
        
        case 'custom'
            
            % Define design matrix X "design_matrix"
            % Use custom matrix for each condition/contrast
            table_obj = DAT.BETWEENPERSON.(mygroupnamefield){c};
            groupnames = table_obj.Properties.VariableNames;
            
                if exist('covs2use','var')
            
                    idx_covar = ismember(groupnames,covs2use);

                        if sum(idx_covar) == 0
                            error('\nOne or more covariates defined in covs2use not present in DAT.BETWEENPERSON.%s{%d}, please correct before proceeding\n',mygroupnamefield,c);
                        end

                    table_obj = table_obj(:,idx_covar);
                    groupnames = groupnames(idx_covar);

                end
                
                % An EMPTY nuisance_covs is legitimate: a single-site tier, or any
                % model with no nuisance regressor at all, has nothing to name here.
                % The old guard fired on empty as well as on unmatched, so a model
                % with nuisance_covs = {} could not run. The message also said
                % covs2use while testing nuisance_covs, which sent debugging the
                % wrong way.
                if exist('nuisance_covs','var') && ~isempty(nuisance_covs)
                    idx_nuisance = ismember(groupnames,nuisance_covs);

                        if sum(idx_nuisance) == 0
                            error(['\nOne or more covariates named in NUISANCE_COVS are not present in ' ...
                                   'DAT.BETWEENPERSON.%s{%d}.\n  nuisance_covs: %s\n  available    : %s\n'], ...
                                   mygroupnamefield, c, strjoin(nuisance_covs, ', '), strjoin(groupnames, ', '));
                        end

                end

            % DUMMY-CODE UNORDERED FACTORS
            %
            % A k-level factor needs k-1 columns. Squeezing one into a single
            % numeric column treats its levels as ordered and spends one degree
            % of freedom where k-1 are needed, so it removes only part of the
            % between-level variance. proj_discoverie's model_3a did exactly
            % this: three centres (UGOT, KUL, UM) in one -1/0/1 column, which
            % both imposes UGOT < KUL < UM and under-adjusts for site.
            %
            % Name such columns in categorical_covs and they are expanded here,
            % before X is built, with the nuisance index expanded to match so a
            % dummy-coded nuisance stays nuisance in regression_stats.
            if exist('categorical_covs','var') && ~isempty(categorical_covs)

                cat_here = categorical_covs(ismember(categorical_covs, groupnames));
                if isempty(cat_here)
                    error(['\ncategorical_covs names none of the columns of ' ...
                           'DAT.BETWEENPERSON.%s{%d}.\nAvailable: %s\n'], ...
                           mygroupnamefield, c, strjoin(groupnames, ', '));
                end

                new_tbl = table();
                new_names = {};
                new_nuis = false(1,0);
                for gi = 1:numel(groupnames)
                    gn = groupnames{gi};
                    is_nuis_gi = exist('idx_nuisance','var') && idx_nuisance(gi);
                    if ismember(gn, cat_here)
                        [Xd, lev] = LaBGAScore_dummy_code(table_obj.(gn));
                        for di = 1:size(Xd,2)
                            dn = matlab.lang.makeValidName(sprintf('%s_%s', gn, lev{di+1}));
                            new_tbl.(dn) = Xd(:,di);
                            new_names{end+1} = dn; %#ok<AGROW>
                            new_nuis(end+1) = is_nuis_gi; %#ok<AGROW>
                        end
                        fprintf('\ndummy-coded ''%s'': %d levels (%s) -> %d column(s), reference ''%s''\n', ...
                            gn, numel(lev), strjoin(lev', ', '), size(Xd,2), lev{1});
                    else
                        new_tbl.(gn) = table_obj.(gn);
                        new_names{end+1} = gn; %#ok<AGROW>
                        new_nuis(end+1) = is_nuis_gi; %#ok<AGROW>
                    end
                end
                table_obj  = new_tbl;
                groupnames = new_names;
                if exist('idx_nuisance','var'), idx_nuisance = new_nuis; end

            end

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
            elseif ismember(DAT.BETWEENPERSON.(mygroupnamefield){c}.Properties.VariableNames,group_id{1})
                group = DAT.BETWEENPERSON.(mygroupnamefield){c}.(group_id{1});
            else
                error('\nGroup not defined in DAT.BETWEENPERSON.group nor DAT.BETWEENPERSON.%s.{%d}, which is required for option "%s" defined in design_matrix_type\n', mygroupnamefield,c,design_matrix_type);
            end
            
            groupnames = {'group'};
                X = group;
                imgs_nan = [];
                groupnames_string = [groupnames_string, ' ', groupnames{1}];

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
    % ---------------------------------------
    
    fprintf('\n\n');
    printhdr('SCALING DATA IF REQUESTED IN OPTIONS');
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
    
    % RESAMPLE MASK SPACE TO IMAGE SPACE
    % ----------------------------------
    
    voxelsize_cat_obj = abs(diag(cat_obj.volInfo.mat(1:3, 1:3)))';
    
    if ~dorobfit_parcelwise
        
        if exist('maskname_short','var')
            voxelsize_glmmask = abs(diag(glmmask.volInfo.mat(1:3, 1:3)))';
            if ~isequal(voxelsize_glmmask,voxelsize_cat_obj)
                glmmask = resample_space(glmmask,cat_obj);
                glmmask.dat(glmmask.dat < 1) = 0; % re-binarize mask, resample_space causes non-zero non-one values at inner cortical boundaries
                
                if c == 1
                
                    % MONTAGE OF RESAMPLED MASK

                    figure;
                    o2 = canlab_results_fmridisplay([], 'compact');
                    o2 = addblobs(o2, glmmask);
                    o2 = title_montage(o2, 5, ['resampled ' maskname_short]);
                    plugin_set_figure_size;
                    drawnow,snapnow;

                    clear o2
                    
                end
                
            end
            
        end
        
    else 
        
        if exist('combined_atlas','var')
            
            voxelsize_atlas = abs(diag(combined_atlas.volInfo.mat(1:3, 1:3)))';
            if ~isequal(voxelsize_atlas,voxelsize_cat_obj)
                combined_atlas = resample_space(combined_atlas,cat_obj);
                
                if c == 1
                
                    % MONTAGE OF RESAMPLED ATLAS

                    figure;
                    o2 = canlab_results_fmridisplay([], 'compact');
                    o2 = addblobs(o2, atlas2region(combined_atlas),'indexmap',cmap,'interp','nearest');
                    if exist('atlasname_short','var')
                        o2 = title_montage(o2, 5, ['resampled ' atlasname_short]);
                    else
                        o2 = title_montage(o2, 5, ['resampled ' atlasname_glm]);
                    end
                    plugin_set_figure_size;
                    drawnow,snapnow;

                    clear o2
                    
                end
                
            end
            
            if doneurotransmitter_maps
                    
                glmmask = fmri_mask_image(combined_atlas);
                voxelsize_glmmask = abs(diag(glmmask.volInfo.mat(1:3, 1:3)))';
                voxelsize_cat_obj = abs(diag(cat_obj.volInfo.mat(1:3, 1:3)))';
                
                if ~isequal(voxelsize_glmmask,voxelsize_cat_obj)
                    glmmask = resample_space(glmmask,cat_obj);
                    glmmask.dat(glmmask.dat < 1) = 0; % re-binarize mask, resample_space causes non-zero non-one values at inner cortical boundaries
                    
                end
                    
            end
  
        end
        
    end
    
    % FORMAT AND ATTACH DESIGN MATRIX
    % -------------------------------
    
    fprintf('\n\n');
    printhdr('CHECKING DESIGN MATRIX');
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
        disp(design_table)
        disp(' ');
        
    end % if loop design_matrix_type
    
    cat_obj.X = X;
    
    % SANITY CHECK ON REGRESSORS, SKIP CONTRAST IF NEEDED
    % ---------------------------------------------------
    
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
    % ----------------------------------------------------------
    
    if domvpa_reg_cov
        
        mvpa_data_objects = cell(size(cat_obj.X,2),1);
        
        for covar = 1:size(cat_obj.X,2)
            
            mvpa_data_objects{covar} = cat_obj;
            mvpa_data_objects{covar}.Y = cat_obj.X(:,covar);
            mvpa_data_objects{covar}.Y_names = groupnames{covar};
            
        end
        
    end
    
    
    %%
    % *EXTRACT ROI AVERAGES*
    
    if doroi_analysis
            
        roi_means{c} = apply_parcellation(cat_obj,roi_atlas);
        roi_means_table{c} = array2table(roi_means{c},'VariableNames',roi_names');
        roi_colors = mat2cell(cmap2,ones(1,size(roi_atlas.labels,2)));
            
        switch mygroupnamefield
        
            case 'contrasts'
                
                fprintf('\n\n');
                printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.contrastnames{c})]);
                fprintf('\n\n');

                switch design_matrix_type

                    case 'custom'
                        roi_means_table{c} = [roi_means_table{c} table_obj];
                        [~, roi_adjusted_means{c}, ~] = barplot_columns(roi_means_table{c}(:,1:end-size(table_obj,2)),'covs',table2array(table_obj),'title',['ROI means, EFFECT: ' DAT.contrastnames{c} ', COVARIATE(S): ' groupnames{:}],'color',roi_colors');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                    case 'group'
                        group_table = array2table(group,'VariableNames',groupnames);
                        roi_means_table{c} = [roi_means_table{c} group_table];
                        [~, roi_adjusted_means{c}, ~] = barplot_columns(roi_means_table{c}(:,1:end-size(group_table,2)),'covs',table2array(group_table),'title',['ROI means, EFFECT: ' DAT.contrastnames{c} ', COVARIATE(S): ' groupnames{:}],'color',roi_colors');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                    case 'onesample'
                        [~, roi_adjusted_means{c}, ~] = barplot_columns(roi_means_table{c},'title',['ROI means, EFFECT: ' DAT.contrastnames{c}],'color',roi_colors');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                end
        
            case 'conditions'
                
                fprintf('\n\n');
                printhdr(['CONDITION #', num2str(c), ': ', upper(DAT.conditions{c})]);
                fprintf('\n\n');
                
                switch design_matrix_type

                    case 'custom'
                        roi_means_table{c} = [roi_means_table{c} table_obj];
                        [~, roi_adjusted_means{c}, ~] = barplot_columns(roi_means_table{c}(:,1:end-size(table_obj,2)),'covs',table2array(table_obj),'title',['ROI means, EFFECT: ' DAT.conditions{c} ', COVARIATE(S): ' groupnames{:}],'color',roi_colors');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                    case 'group'
                        group_table = array2table(group,'VariableNames',groupnames);
                        roi_means_table{c} = [roi_means_table{c} group_table];
                        [~, roi_adjusted_means{c}, ~] = barplot_columns(roi_means_table{c}(:,1:end-size(group_table,2)),'covs',table2array(group_table),'title',['ROI means, EFFECT: ' DAT.conditions{c} ', COVARIATE(S): ' groupnames{:}],'color',roi_colors');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                    case 'onesample'
                        [~, roi_adjusted_means{c}, ~] = barplot_columns(roi_means_table{c},'title',['ROI means, EFFECT: ' DAT.conditions{c}],'color',roi_colors');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                end
    
        end
               

        %%
        % *GLM ON ROI MEANS*
        %
        % Formal inference on the roi averages. Until this existed the roi
        % analysis produced a picture and no statistics: barplot_columns above
        % adjusts the plotted means for covariates but reports no per-roi test,
        % so an roi effect could not be called significant without refitting
        % the model by hand.
        %
        % The split between effects of interest and nuisance is NOT a new
        % option: it reuses nuisance_covs, the same variable the voxelwise GLM
        % uses to set regression_stats.nuisance_columns and to build the
        % Freedman-Lane permutations in the TFCE branch. Anything in the design
        % that is not named there is an effect of interest. Note this is
        % deliberately NOT covs2use, which SUBSETS the design matrix (dropping
        % every covariate not listed) rather than labelling its columns.
        %
        % Two levels of inference are reported:
        %
        %   1. MANOVA across the whole roi set, adjusted for nuisance. One
        %      omnibus test of whether the roi profile differs, which does not
        %      spend multiple comparisons and is the univariate analogue of a
        %      multivariate classifier trained on the same rois.
        %   2. A GLM per roi, FDR-corrected across the set.
        %
        % Reported in that order deliberately: if the omnibus test is null,
        % individual rois surviving FDR should be read with that in mind.
        if exist('doroi_glm','var') && doroi_glm

            if isequal(design_matrix_type,'onesample')

                fprintf('\nroi GLM skipped: design_matrix_type is ''onesample'', so there are no covariates to model.\n');

            else

                switch design_matrix_type
                    case 'custom'
                        cov_tbl_glm = table_obj;
                    case 'group'
                        cov_tbl_glm = group_table;
                end

                covnames_glm = cov_tbl_glm.Properties.VariableNames;
                roinames_glm = roi_means_table{c}.Properties.VariableNames;
                roinames_glm = roinames_glm(~ismember(roinames_glm, covnames_glm));

                % nuisance from the same option the voxelwise GLM uses
                if exist('nuisance_covs','var') && ~isempty(nuisance_covs)
                    nuis_glm = covnames_glm(ismember(covnames_glm, nuisance_covs));
                else
                    nuis_glm = {};
                end
                eoi_glm = covnames_glm(~ismember(covnames_glm, nuis_glm));

                if isempty(eoi_glm)
                    error(['\nroi GLM: every covariate in the design is named in nuisance_covs, ' ...
                           'so there is no effect of interest left to test.\n']);
                end

                Xglm = table2array(cov_tbl_glm(:, [eoi_glm nuis_glm]));
                if ~isnumeric(Xglm)
                    error('\nroi GLM predictors must be numeric; %s are not.\n', ...
                        strjoin([eoi_glm nuis_glm], ', '));
                end
                Yglm = table2array(roi_means_table{c}(:, roinames_glm));

                nroi_glm = numel(roinames_glm);
                ncov_glm = numel(eoi_glm);
                nnui_glm = numel(nuis_glm);

                fprintf('\n\n');
                if isempty(nuis_glm)
                    printhdr(['GLM ON ROI MEANS: ' upper(strjoin(eoi_glm, ', '))]);
                else
                    printhdr(['GLM ON ROI MEANS: ' upper(strjoin(eoi_glm, ', ')) ...
                              ' controlling for ' upper(strjoin(nuis_glm, ', '))]);
                end
                fprintf('\nmodel: roi_mean ~ %s\n', strjoin([eoi_glm nuis_glm], ' + '));
                fprintf('n = %d, %d roi(s)\n', size(Xglm,1), nroi_glm);


                %%
                % *MANOVA ACROSS ROIS*
                %
                % Wilks' Lambda from the full model against a reduced model
                % without the effect being tested, so nuisance covariates are
                % adjusted for rather than ignored. manova1 cannot do this - it
                % is one-way and takes no covariates - so the SSCP matrices are
                % formed directly and converted to Rao's F.
                roi_manova_stats{c} = table();

                if nroi_glm < 2

                    fprintf('\nMANOVA skipped: needs 2 or more rois, %d present.\n', nroi_glm);

                else

                    Xfull_man = [ones(size(Xglm,1),1) Xglm];
                    dfe_man   = size(Xfull_man,1) - rank(Xfull_man);

                    if dfe_man <= nroi_glm
                        fprintf(['\nMANOVA skipped: %d roi(s) but only %d error df. ' ...
                                 'The residual covariance is singular, so Wilks'' Lambda is undefined.\n'], ...
                                 nroi_glm, dfe_man);
                    else

                        Rfull_man = Yglm - Xfull_man*(Xfull_man\Yglm);
                        E_man     = Rfull_man' * Rfull_man;

                        man_eff = {}; man_lam = []; man_F = []; man_df1 = []; man_df2 = []; man_p = [];

                        for ee = 1:ncov_glm
                            % reduced model: drop just this effect, keep the rest
                            keep_man = true(1, size(Xglm,2));
                            keep_man(ee) = false;
                            Xred_man = [ones(size(Xglm,1),1) Xglm(:,keep_man)];
                            Rred_man = Yglm - Xred_man*(Xred_man\Yglm);
                            E0_man   = Rred_man' * Rred_man;

                            lambda_man = det(E_man) / det(E0_man);   % = det(E)/det(E+H)

                            pp_man = nroi_glm;                        % dependent variables
                            vh_man = 1;                               % hypothesis df (one column)
                            ve_man = dfe_man;

                            denom_man = pp_man^2 + vh_man^2 - 5;
                            if denom_man > 0
                                t_man = sqrt((pp_man^2*vh_man^2 - 4) / denom_man);
                            else
                                t_man = 1;
                            end
                            df1_man = pp_man * vh_man;
                            df2_man = t_man*(ve_man - (pp_man - vh_man + 1)/2) - (pp_man*vh_man - 2)/2;
                            lam_t   = lambda_man^(1/t_man);
                            F_man   = ((1 - lam_t)/lam_t) * (df2_man/df1_man);
                            p_man   = 1 - fcdf(F_man, df1_man, df2_man);

                            man_eff{end+1,1} = eoi_glm{ee};
                            man_lam(end+1,1) = lambda_man;
                            man_F(end+1,1)   = F_man;
                            man_df1(end+1,1) = df1_man;
                            man_df2(end+1,1) = df2_man;
                            man_p(end+1,1)   = p_man;
                        end

                        roi_manova_stats{c} = table(man_eff, man_lam, man_F, man_df1, man_df2, man_p, ...
                            'VariableNames', {'effect','wilks_lambda','F','df1','df2','p'});

                        fprintf('\n');
                        printhdr(['MANOVA ACROSS ' num2str(nroi_glm) ' ROIS (Wilks'' Lambda, Rao''s F)']);
                        fprintf('\n');
                        if nnui_glm > 0
                            fprintf('each effect tested against a reduced model, adjusted for %s\n\n', ...
                                strjoin(nuis_glm, ', '));
                        else
                            fprintf('each effect tested against a reduced model (no nuisance covariates)\n\n');
                        end
                        disp(roi_manova_stats{c});

                        for ee = 1:height(roi_manova_stats{c})
                            if roi_manova_stats{c}.p(ee) < 0.05
                                fprintf('%s: roi profile DIFFERS, F(%.0f,%.1f) = %.3f, p = %.4f\n', ...
                                    roi_manova_stats{c}.effect{ee}, roi_manova_stats{c}.df1(ee), ...
                                    roi_manova_stats{c}.df2(ee), roi_manova_stats{c}.F(ee), ...
                                    roi_manova_stats{c}.p(ee));
                            else
                                fprintf('%s: no omnibus difference across rois, F(%.0f,%.1f) = %.3f, p = %.4f\n', ...
                                    roi_manova_stats{c}.effect{ee}, roi_manova_stats{c}.df1(ee), ...
                                    roi_manova_stats{c}.df2(ee), roi_manova_stats{c}.F(ee), ...
                                    roi_manova_stats{c}.p(ee));
                            end
                        end

                    end % enough df

                end % enough rois


                %%
                % *GLM PER ROI*

                B = nan(nroi_glm, ncov_glm); SE = B; TT = B; PP = B; DF = B; ES = B;

                for rr = 1:nroi_glm
                    y_glm = double(roi_means_table{c}.(roinames_glm{rr}));
                    mdl_glm = fitlm(Xglm, y_glm, 'VarNames', ...
                        [eoi_glm nuis_glm {roinames_glm{rr}}]);
                    for ee = 1:ncov_glm
                        k = ee + 1;   % +1 for the intercept
                        B(rr,ee)  = mdl_glm.Coefficients.Estimate(k);
                        SE(rr,ee) = mdl_glm.Coefficients.SE(k);
                        TT(rr,ee) = mdl_glm.Coefficients.tStat(k);
                        PP(rr,ee) = mdl_glm.Coefficients.pValue(k);
                        DF(rr,ee) = mdl_glm.DFE;
                        lv = unique(Xglm(~isnan(Xglm(:,ee)), ee));
                        if numel(lv) == 2
                            ES(rr,ee) = B(rr,ee) * diff(lv) / mdl_glm.RMSE;
                        else
                            ES(rr,ee) = sign(TT(rr,ee)) * ...
                                sqrt(TT(rr,ee)^2 / (TT(rr,ee)^2 + mdl_glm.DFE));
                        end
                    end
                end

                eff_col = {}; roi_col = {}; est_col = []; se_col = []; t_col = [];
                df_col = []; p_col = []; qbh_col = []; qst_col = []; es_col = []; esname_col = {};
                strel_col = []; pi0_col = [];   % Storey reliability verdict and pi0, per effect

                for ee = 1:ncov_glm
                    p_ee   = PP(:,ee);
                    % FDR through the lab's canonical implementation, so prep_3a and the
                    % decoding scripts cannot drift apart. LaBGAScore_Storey_FDR defaults to
                    % SAS PROC MULTTEST's PFDR (spline, falling back to the Storey &
                    % Tibshirani bootstrap on SAS's own trigger), estimates pi0, judges
                    % whether pi0 is identifiable at all, and returns Benjamini-Hochberg when
                    % it is not. It prints its own diagnostics into the report.
                    %
                    % This replaced an inline copy that took mafdr's spline pi0 and guarded
                    % only aprioriprob > 0.99. That catches the conservative failure but not
                    % pi0 -> 0: on these very 8 roi p-values the spline returned pi0 = 0.012,
                    % every q fell below its own p, and the q >= p floor turned the column
                    % back into the RAW p-values under the heading q_Storey.
                    %
                    % Note roi means from the same subjects are strongly correlated, which
                    % violates Storey's independence assumption; BH holds under positive
                    % regression dependency. When the verdict is unreliable, q_Storey below
                    % simply equals q_BH.
                    fprintf('\n  %s:', eoi_glm{ee});
                    [qst_ee, pi0_ee, storey_info_ee] = LaBGAScore_Storey_FDR(p_ee);

                    qbh_ee             = storey_info_ee.q_BH(:);
                    storey_reliable_ee = storey_info_ee.reliable;
                    qst_ee             = qst_ee(:);

                    lv = unique(Xglm(~isnan(Xglm(:,ee)), ee));
                    if numel(lv) == 2
                        est_ee = B(:,ee) * diff(lv);      % difference between the two levels
                        esn = 'cohens_d';
                    else
                        est_ee = B(:,ee);                 % slope
                        esn = 'partial_r';
                    end
                    eff_col    = [eff_col;    repmat(eoi_glm(ee), nroi_glm, 1)];
                    roi_col    = [roi_col;    roinames_glm'];
                    est_col    = [est_col;    est_ee];
                    se_col     = [se_col;     SE(:,ee)];
                    t_col      = [t_col;      TT(:,ee)];
                    df_col     = [df_col;     DF(:,ee)];
                    p_col      = [p_col;      p_ee];
                    qbh_col    = [qbh_col;    qbh_ee];
                    qst_col    = [qst_col;    qst_ee];
                    es_col     = [es_col;     ES(:,ee)];
                    esname_col = [esname_col; repmat({esn}, nroi_glm, 1)];
                    strel_col  = [strel_col;  repmat(storey_reliable_ee, nroi_glm, 1)];
                    pi0_col    = [pi0_col;    repmat(pi0_ee, nroi_glm, 1)];
                end

                roi_glm_stats{c} = table(eff_col, roi_col, est_col, se_col, t_col, ...
                    df_col, p_col, qbh_col, qst_col, logical(strel_col), pi0_col, ...
                    es_col, esname_col, 'VariableNames', ...
                    {'effect','roi','estimate','se','t','df','p','q_BH','q_Storey', ...
                     'storey_reliable','pi0','effect_size','effect_size_type'});

                fprintf('\n');
                printhdr(['GLM PER ROI, FDR ACROSS ' num2str(nroi_glm) ' ROIS']);
                fprintf('\n');
                disp(roi_glm_stats{c});

                for ee = 1:ncov_glm
                    lv  = unique(Xglm(~isnan(Xglm(:,ee)), ee));
                    sel = strcmp(roi_glm_stats{c}.effect, eoi_glm{ee});
                    qb  = roi_glm_stats{c}.q_BH(sel);
                    qs  = roi_glm_stats{c}.q_Storey(sel);
                    if numel(lv) == 2
                        fprintf(['\n%s: ''estimate'' is the difference between %s = %g and %s = %g ' ...
                                 '(positive = higher at %g).\n'], ...
                            eoi_glm{ee}, eoi_glm{ee}, lv(2), eoi_glm{ee}, lv(1), lv(2));
                    end
                    if any(qb < 0.05)
                        fprintf('%s: %d/%d roi(s) at q_BH < .05: %s\n', eoi_glm{ee}, ...
                            sum(qb < 0.05), nroi_glm, strjoin(roinames_glm(qb < 0.05), ', '));
                    else
                        fprintf('%s: no roi at q_BH < .05 (smallest q_BH = %.4f)\n', ...
                            eoi_glm{ee}, min(qb));
                    end
                    rel_ee = roi_glm_stats{c}.storey_reliable(find(sel,1));
                    if any(qs < 0.05)
                        fprintf('%s: %d/%d roi(s) at q_Storey < .05: %s\n', eoi_glm{ee}, ...
                            sum(qs < 0.05), nroi_glm, strjoin(roinames_glm(qs < 0.05), ', '));
                    else
                        fprintf('%s: no roi at q_Storey < .05 (smallest q_Storey = %.4f)\n', ...
                            eoi_glm{ee}, min(qs));
                    end
                    if ~rel_ee
                        fprintf(['%s: the q_Storey column above is reported for completeness only - ' ...
                                 'its pi0 is not trustworthy for this set, so read q_BH.\n'], eoi_glm{ee});
                    end
                end

            end % onesample check

        end % if roi glm requested

    end % if loop roi analysis
    
    
    %%
    % *CALCULATE SIMILARITY WITH NEUROTRANSMITTER MAPS*
    
    if doneurotransmitter_maps
            
        switch mygroupnamefield
        
            case 'contrasts'
                
                fprintf('\n\n');
                printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.contrastnames{c})]);
                fprintf('\n\n');
                
                switch neurotransmitter_maps_metric
                    
                    case 'correlation'
                        
                        if exist('glmmask','var')               
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'doAverage','mask',glmmask);
                        else
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'doAverage');
                        end
                            
                        title(DAT.contrastnames{c},'Interpreter','none');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                        if isequal(design_matrix_type,'group') || (isequal(design_matrix_type,'custom') && ~isempty(DAT.BETWEENPERSON.group))
                            
                            if exist('glmmask','var')
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'doAverage','compareGroups',DAT.BETWEENPERSON.group, 'mask',glmmask);
                            else
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'doAverage','compareGroups',DAT.BETWEENPERSON.group);
                            end
                            
                            title(DAT.contrastnames{c},'Interpreter','none');
                            plugin_set_figure_size;
                            drawnow,snapnow;

                        end
                        
                    case 'cosine_similarity'
                        
                        if exist('glmmask','var') 
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage','mask',glmmask);
                        else
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage');
                        end
                        
                        title(DAT.contrastnames{c},'Interpreter','none');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                        if isequal(design_matrix_type,'group') || (isequal(design_matrix_type,'custom') && ~isempty(DAT.BETWEENPERSON.group))
                            
                            if exist('glmmask','var')
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage','compareGroups',DAT.BETWEENPERSON.group, 'mask', glmmask);
                            else
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage','compareGroups',DAT.BETWEENPERSON.group);
                            end
                            
                            title(DAT.contrastnames{c} ,'Interpreter','none');
                            plugin_set_figure_size;
                            drawnow,snapnow;

                        end
                        
                end  % switch similarity metric
        
            case 'conditions'
                
                fprintf('\n\n');
                printhdr(['CONDITION #', num2str(c), ': ', upper(DAT.conditions{c})]);
                fprintf('\n\n');
                
                switch neurotransmitter_maps_metric
                    
                    case 'correlation'
                        
                        if exist('glmmask','var')               
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'doAverage','mask',glmmask);
                        else
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'doAverage');
                        end

                        title(DAT.conditions{c},'Interpreter','none');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                        if isequal(design_matrix_type,'group') || (isequal(design_matrix_type,'custom') && ~isempty(DAT.BETWEENPERSON.group))

                            if exist('glmmask','var')
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'doAverage','compareGroups',DAT.BETWEENPERSON.group, 'mask',glmmask);
                            else
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'doAverage','compareGroups',DAT.BETWEENPERSON.group);
                            end
                            
                            title(DAT.conditions{c} ,'Interpreter','none');
                            plugin_set_figure_size;
                            drawnow,snapnow;

                        end
                        
                    case 'cosine_similarity'
                        
                        if exist('glmmask','var') 
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage','mask',glmmask);
                        else
                            neurotransmitter_stats{c} = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage');
                        end
                        
                        title(DAT.conditions{c},'Interpreter','none');
                        plugin_set_figure_size;
                        drawnow,snapnow;

                        if isequal(design_matrix_type,'group') || (isequal(design_matrix_type,'custom') && ~isempty(DAT.BETWEENPERSON.group))

                            if exist('glmmask','var')
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage','compareGroups',DAT.BETWEENPERSON.group, 'mask', glmmask);
                            else
                                [neurotransmitter_group_stats{c},~,~,~,neurotransmitter_group_tables{c},neurotransmitter_multcomp_group{c}] = hansen_neurotransmitter_maps(cat_obj,'cosine_similarity','doAverage','compareGroups',DAT.BETWEENPERSON.group);
                            end
                            
                            title(DAT.conditions{c} ,'Interpreter','none');
                            plugin_set_figure_size;
                            drawnow,snapnow;

                        end
                        
                end % switch similarity metric
    
        end % switch conditions or contrasts
               
    end % if loop neurotransmitter maps
    
   
    %%
    % *RUN GLM MODEL*
    
    % VOXEL-WISE
    % ----------
    
    switch mygroupnamefield
        
        case 'contrasts'
            fprintf('\n\n');
            printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.contrastnames{c})]);
            fprintf('\n\n');
            
        case 'conditions'
            fprintf('\n\n');
            printhdr(['CONDITION #', num2str(c), ': ', upper(DAT.conditions{c})]);
            fprintf('\n\n');
    
    end
    
    if ~dorobfit_parcelwise
        
        if dorobust
            robuststring = 'robust';
            regresstime = tic;
        else
            robuststring = 'norobust';
        end
        
        fprintf('\n\n');
        printhdr(['RUNNING VOXEL-WISE ', upper(robuststring) ' REGRESSION']);
        fprintf('\n\n');

        if ~strcmpi(design_matrix_type,'onesample')
            % regression_stats.t has t maps for all regressors, intercept is last
            switch mygroupnamefield
                case 'contrasts'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.contrastnames{c}, 'variable_names', groupnames, 'nodisplay','residual'); % trick to get unthresholded maps
                case 'conditions'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.conditions{c}, 'variable_names', groupnames, 'nodisplay','residual');
            end
        else
            % regression_stats.t has t maps for intercept only
            switch mygroupnamefield
                case 'contrasts'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.contrastnames{c}, 'variable_names', groupnames, 'nointercept', 'nodisplay','residual');
                case 'conditions'
                    regression_stats = regress(cat_obj, 1, 'unc', robuststring, 'analysis_name', DAT.conditions{c}, 'variable_names', groupnames, 'nointercept', 'nodisplay','residual');
            end
        end

        
        % RUN DIAGNOSTICS ON FITTED MODEL AND SUMMARIZE
        
        % idx_nuisance is only created inside "if exist('nuisance_covs','var')"
        % further up, so a design with no nuisance covariates - the common case
        % for a plain group comparison - reaches here with it undefined and the
        % whole script dies AFTER the regression has been computed but BEFORE
        % anything is saved. Default it to "no nuisance columns" instead.
        if ~exist('idx_nuisance','var')
            idx_nuisance = false(size(groupnames));
        end
        
        regression_stats.nuisance_columns = find(idx_nuisance);
        
        regression_stats = validate_object(regression_stats);
        regression_stats = run_diagnostics(regression_stats);
        summary(regression_stats);
        
        
        if doBayes
            
            % CALCULATE BAYES FACTORS FROM T-MAPŜ AND SAVE TO SEPARATE
            % RESULTS STRUCT
            
            fprintf('\n\n');
            printhdr('Calculating voxel-wise Bayes Factor maps');
            fprintf('\n\n');
            
            N = single(sum(~(isnan(cat_obj.dat') | cat_obj.dat' == 0) , 1)); % code from fmri_data.ttest to correct N in regressions_stat.t to make estimateBayesFactor function work on regress() output
            
            for reg = 1:size(regression_stats.t.dat,2)
                
                t_for_Bayes = get_wh_image(regression_stats.t, reg);
                t_for_Bayes.N = N';
                bayesian_regression_stats.BF(reg) = estimateBayesFactor(t_for_Bayes,'t');
            
            end
            
        end
        
        LaBGAScore_smart_parallel_pool_setup;      
        
        % TFCE can be restricted to a subset of contrasts: it is the dominant cost
        % of this script, and most designs have one contrast of interest.
        if ~exist('cons2tfce','var')
            cons2tfce = [];
        end
        
        if doTFCE && (isempty(cons2tfce) || ismember(c,cons2tfce))
            
            % CALCULATE TFCE STATS FROM DATA OBJECT
            
            % Resolve ONE seed for this run and pass it to every call below, so the
            % permutation null is reproducible. group_tfce_from_subject_maps permutes
            % inside a parfor, whose workers a client-side rng() never reaches, so
            % without this the TFCE maps differ from run to run - which matters most
            % where it is least visible, near threshold and at the p floor.
            % Set tfce_seed in a2_set_default_options to fix it across runs; left
            % unset, a seed is drawn here and reported, so the run stays random but
            % can be reproduced exactly afterwards.
            if ~exist('tfce_seed','var') || isempty(tfce_seed)
                tfce_seed = randi(2^31-1);
            end
            fprintf('\nTFCE permutation seed for this run: %d\n', tfce_seed);
            
            % MASK THE TFCE INPUT
            % -------------------------------------------------------------
            % TFCE is the one statistic prep_3a produces that is already
            % CORRECTED: group_tfce_from_subject_maps builds its own
            % max-statistic permutation null and returns FWE p-values. Every
            % other map here is uncorrected, with masking and correction left
            % to c2a - so TFCE is the exception, and its correction has to see
            % the mask, or it is computed over the whole image extent while the
            % FDR beside it is computed within grey matter.
            %
            % That is not a small difference. The null is a MAXIMUM over
            % voxels, so including non-grey-matter and edge voxels inflates it,
            % and every in-mask voxel is then judged against competitors that
            % are not part of the analysis.
            %
            % Only a COPY is masked. cat_obj itself stays whole-brain, so the
            % regression, its t-map and its p-values remain uncorrected and
            % unmasked exactly as before, and c2a keeps full freedom to mask
            % them at reporting time.
            %
            % Set mask_tfce_input = false to reproduce the previous behaviour.
            if ~exist('mask_tfce_input','var') || isempty(mask_tfce_input)
                mask_tfce_input = true;
            end
            cat_obj_tfce = cat_obj;
            if mask_tfce_input && exist('glmmask','var') && ~isempty(glmmask)
                glmmask_tfce = glmmask;
                vs_m = abs(diag(glmmask_tfce.volInfo.mat(1:3,1:3)))';
                vs_d = abs(diag(cat_obj.volInfo.mat(1:3,1:3)))';
                if ~isequal(vs_m, vs_d)
                    glmmask_tfce = resample_space(glmmask_tfce, cat_obj);
                    glmmask_tfce.dat(glmmask_tfce.dat < 1) = 0;
                end
                n_before_tfce = size(cat_obj_tfce.dat,1);
                cat_obj_tfce = apply_mask(cat_obj_tfce, glmmask_tfce);
                fprintf(['\nTFCE input masked with %s: %d of %d voxels (%.1f%%).\n' ...
                         'The permutation null is therefore a maximum over the MASKED volume,\n' ...
                         'matching the FDR that c2a computes after masking. cat_obj itself is\n' ...
                         'untouched, so the regression stays whole-brain and uncorrected.\n\n'], ...
                         maskname_short, size(cat_obj_tfce.dat,1), n_before_tfce, ...
                         100*size(cat_obj_tfce.dat,1)/n_before_tfce);
            elseif mask_tfce_input
                fprintf('\nmask_tfce_input is true but no glmmask exists; TFCE runs on the full extent\n\n');
            end

            fprintf('\n\n');
            printhdr('Calculating voxel-wise TFCE maps');
            fprintf('\n\n');
            
            switch design_matrix_type
                
                case 'onesample'
                    
                    switch tfce_sidedness
                        
                        case 'two'
                            
                            [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'onesample',[],[],perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness);
                        
                        case 'one'
                            
                            [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'onesample',[],[],perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness,'tail',tfce_tail);
                            
                    end
                    
                case 'group'
                    
                    switch tfce_sidedness
                        
                        case 'two'
                            
                            [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'twosample',DAT.BETWEENPERSON.group,[],perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness);
                        
                        case 'one'
                            
                            [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'twosample',DAT.BETWEENPERSON.group,[],perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness,'tail',tfce_tail);
                            
                    end
                    
                case 'custom'
                    
                    if ~isempty(DAT.BETWEENPERSON.group)
                              
                           
                       if exist('nuisance_covs','var') && ~isempty(nuisance_covs)
                               
                           switch tfce_sidedness

                                case 'two'

                                    [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'twosample',DAT.BETWEENPERSON.group,regression_stats.X(:,regression_stats.wh_nuisance),perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness);

                                case 'one'

                                    [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'twosample',DAT.BETWEENPERSON.group,regression_stats.X(:,regression_stats.wh_nuisance),perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness,'tail',tfce_tail);
                                    
                           end
                           
                       else
                           
                           switch tfce_sidedness
                        
                                case 'two'

                                    [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'twosample',DAT.BETWEENPERSON.group,[],perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness);

                                case 'one'

                                    [tfce_dat,tfce_stat_img,tfce_info] = group_tfce_from_subject_maps(cat_obj_tfce,'twosample',DAT.BETWEENPERSON.group,[],perm_n_tfce,'seed',tfce_seed,'sidedness',tfce_sidedness,'tail',tfce_tail);
                            
                            end

                       end
                       
                    else
                        
                        warning('TFCE stats not implemented yet for continuous regressors, skipping TFCE for this contrast')
                        
                    end
                    
            end % switch design_matrix_type
            
            tfce_regression_stats.tfce_dat = tfce_dat;
            tfce_regression_stats.tfce_stat_img = tfce_stat_img;
            tfce_regression_stats.tfce_info = tfce_info;
            
            fprintf('\nMaximum real TFCE = %g\n',tfce_info.TFCE_real_max);
            fprintf('\nMaximum null TFCE = %g\n',max(tfce_info.TFCE_null_max));
            fprintf('\nGlobal TFCE p-value = %g\n',tfce_info.p_TFCE_global);
                    
            
        end % if doTFCE

        
        % PLOT MONTAGE (MASKED IF SPECIFIED IN MASKNAME_GLM)

        fprintf('\n\n');
        printhdr('Plotting voxel-wise GLM results');
        fprintf('\n\n');
        
        t = regression_stats.t;
        
        t = apply_mask(t,brainmask); % re-apply brainmask just to be sure
        
        if exist('maskname_short','var')
            t = apply_mask(t,glmmask);
        end
        
        t = threshold(t,.05,'unc');
        
        fprintf ('\nMONTAGE VOXELWISE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR(S): %s, MASK: %s, SCALING: %s\n\n', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
                
        num_effects = size(t.dat, 2); % number of regressors
        o2 = canlab_results_fmridisplay([], 'multirow', num_effects);

        for j = 1:num_effects

            tj = get_wh_image(t, j);
            tj = threshold(tj, .05, 'unc'); 
            
            datsig = tj.dat(logical(tj.sig));
            datsigneg = datsig(datsig<0);
            datsigpos = datsig(datsig>0);

                if isempty(datsigneg) && ~isempty(datsigpos)

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'mincolor',[.9 .4 0], 'maxcolor', [1 1 0]);%, 'cmaprange', [min(datsigpos) max(datsigpos)]);

                elseif isempty(datsigpos) && ~isempty(datsigneg)

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'mincolor',[.1 .8 .8], 'maxcolor', [.1 .1 .8]);%, 'cmaprange', [min(datsigneg) max(datsigneg)]);

                else

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});%, 'cmaprange', [min(datsigneg) max(datsigneg) min(datsigpos) max(datsigpos)]);

                end
                
                if num_effects < 4
                    o2 = legend(o2);
                end
            
            o2 = title_montage(o2, 2*j, [regression_stats.analysis_name ' ' regression_stats.variable_names{j} ' ' mask_string ' ' scaling_string]);

        end

        figtitle = sprintf('%s_05_unc_montage_%s_%s_%s', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
        set(gcf, 'Tag', figtitle); plugin_set_figure_size;
        drawnow, snapnow;
            if save_figures_glm
                plugin_save_figure;
            end
        clear o2, clear figtitle, clear j, clear tj
        
        if doBayes
            
            fprintf('\n\n');
            printhdr('Plotting voxel-wise Bayesian GLM results');
            fprintf('\n\n');
            
            BF = bayesian_regression_stats.BF;
            
            fprintf ('\nMONTAGE VOXELWISE BAYESIAN GLM RESULTS AT |BF| > 3, EFFECT: %s, REGRESSOR(S): %s, MASK: %s, SCALING: %s\n\n', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
            
            o2 = canlab_results_fmridisplay([], 'multirow', num_effects);
            
            for img = 1:size(BF,2)
                
                BF(img) = apply_mask(BF(img),brainmask);
                
                    if exist('maskname_short','var')
                        BF(img) = apply_mask(BF(img),glmmask);
                    end
                
                BF(img) = threshold(BF(img),[-2.1972 2.1972],'raw-outside');
                    
                datsig = BF(img).dat(logical(BF(img).sig));
                datsigneg = datsig(datsig<0);
                datsigpos = datsig(datsig>0);
                
                    if isempty(datsigneg) && ~isempty(datsigpos)

                        o2 = addblobs(o2, region(BF(img)), 'wh_montages', (2*img)-1:2*img, 'mincolor',[0 0.25 0], 'maxcolor', [0 1 0]);%, 'cmaprange', [min(datsigpos) max(datsigpos)]);

                    elseif isempty(datsigpos) && ~isempty(datsigneg)

                        o2 = addblobs(o2, region(BF(img)), 'wh_montages', (2*img)-1:2*img, 'mincolor',[.25 0 0], 'maxcolor', [1 0 0]);%, 'cmaprange', [min(datsigneg) max(datsigneg)]);

                    else

                        o2 = addblobs(o2, region(BF(img)), 'wh_montages', (2*img)-1:2*img, 'splitcolor',{[.25 0 0] [1 0 0] [0 0.25 0] [0 1 0]});%, 'cmaprange', [min(datsigneg) max(datsigneg) min(datsigpos) max(datsigpos)]); % red in favor of H0, green in favor of H1 for BF maps

                    end

                    if size(BF,2) < 4
                        o2 = legend(o2);
                    end
                
                o2 = title_montage(o2, 2*img, [regression_stats.analysis_name ' ' regression_stats.variable_names{img} ' ' mask_string ' ' scaling_string]);
            
            end

            figtitle = sprintf('%s_BF_3_montage_%s_%s_%s', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
            set(gcf, 'Tag', figtitle); plugin_set_figure_size;
            drawnow, snapnow;
                if save_figures_glm
                    plugin_save_figure;
                end
            clear o2, clear figtitle, clear img, clear BF
            
        end
        
        if doTFCE && (isempty(cons2tfce) || ismember(c,cons2tfce))
            
            fprintf('\n\n');
            printhdr('Plotting voxel-wise TFCE GLM results');
            fprintf('\n\n');
            
            tfce_stat_img_thr_unc_05 = threshold(tfce_stat_img,0.05,'unc');
            tfce_dat_thr_unc_05 = thresholded_fmri_data_from_statistic_image(tfce_stat_img_thr_unc_05,tfce_dat.dat,combined_atlas,0.05,'tfce','unc',0);
            
            tfce_regression_stats.tfce_stat_img_thr_unc_05 = tfce_stat_img_thr_unc_05;
            tfce_regression_stats.tfce_dat_thr_unc_05 = tfce_dat_thr_unc_05;
            
                if exist('maskname_short','var')
                    tfce_dat_thr_unc_05 = apply_mask(tfce_dat_thr_unc_05, glmmask);
                end
                
                switch design_matrix_type
                    
                    case 'onesample'
            
                        fprintf ('\nMONTAGE VOXELWISE TFCE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
                        figtitle = sprintf('%s_TFCE_05_unc_montage_%s_%s_%s', regression_stats.analysis_name, groupnames_string, mask_string, scaling_string);
                        
                    case 'group'
                        
                        fprintf ('\nMONTAGE VOXELWISE TFCE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', regression_stats.analysis_name, strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), ', '), mask_string, scaling_string);
                        figtitle = sprintf('%s_TFCE_05_unc_montage_%s_%s_%s', regression_stats.analysis_name, strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), '_'), mask_string, scaling_string);
                        
                    case 'custom'
                        
                        if exist('nuisance_covs','var') && ~isempty(nuisance_covs)
                            
                            fprintf ('\nMONTAGE VOXELWISE TFCE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR: %s, NUISANCE COVARIATE(S): %s, MASK: %s, SCALING: %s\n\n', regression_stats.analysis_name, strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), ', '), strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_nuisance)), ', '), mask_string, scaling_string);
                            figtitle = sprintf('%s_TFCE_05_unc_montage_%s_nuisance_%s_%s_%s', regression_stats.analysis_name, strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), '_'), strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_nuisance)), '_'), mask_string, scaling_string);
                            
                        else
                        
                            fprintf ('\nMONTAGE VOXELWISE TFCE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR: %s, MASK: %s, SCALING: %s\n\n', regression_stats.analysis_name, strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), ', '), mask_string, scaling_string);
                            figtitle = sprintf('%s_TFCE_05_unc_montage_%s_%s_%s', regression_stats.analysis_name, strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), '_'), mask_string, scaling_string);
                            
                        end
                        
                end
            
            % montage() on a DATA object routes through canlab_results_fmridisplay WITHOUT
            % create_figure, so it draws into whatever figure is current - which is the
            % previous block's montage. Verified: after canlab_results_fmridisplay the
            % figure count stays 1 and the axes accumulate, so the TFCE blobs landed on top
            % of the Bayes montage and both were captured in one snapnow. Open a fresh
            % figure first. (region montages are fine - they go through create_figure.)
            figure;

            o2 = montage(tfce_dat_thr_unc_05,'mincolor',[0.47 0.11 0.43], 'maxcolor', [0.94 0.98 0.13]);
            o2 = title_montage(o2, 5, ['tfce ' regression_stats.analysis_name ' ' strjoin(cellstr(regression_stats.variable_names(regression_stats.wh_interest)), ', ') ' ' mask_string ' ' scaling_string]);
            set(gcf, 'Tag', figtitle); plugin_set_figure_size;
            drawnow, snapnow;
                if save_figures_glm
                    plugin_save_figure;
                end
            clear o2, clear figtitle
            
        end

        % KEEP RESULTS OBJECTS IN CELL ARRAY FOR SAVING

        regression_stats_results{c} = regression_stats;
        
        if doBayes
            bayesian_regression_stats_results{c} = bayesian_regression_stats;
        end
        
        if doTFCE && (isempty(cons2tfce) || ismember(c,cons2tfce))
            tfce_regression_stats_results{c} = tfce_regression_stats;
        end

        
    % PARCEL-WISE
    % -----------
        
    else
        
        fprintf('\n\n');
        printhdr('RUNNING PARCEL-WISE ROBUST REGRESSION');
        fprintf('\n\n');
        
        if exist('combined_atlas','var')
        
            if csf_wm_covs && remove_outliers
                parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',true,'doplot',false,'mask',combined_atlas);
            elseif csf_wm_covs && ~remove_outliers
                parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',true,'remove_outliers',false,'doplot',false,'mask',combined_atlas);
            elseif ~csf_wm_covs && remove_outliers
                parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'csf_wm_covs',false,'remove_outliers',true,'doplot',false,'mask',combined_atlas);
            else
                parcelwise_stats = robfit_parcelwise(cat_obj,'names', groupnames,'doplot',false,'mask',combined_atlas);
            end
            
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
            
        end
        
        if doBayes
            
            % CALCULATE BAYES FACTORS FROM T-MAPŜ AND ADD TO RESULTS
            
            fprintf('\n\n');
            printhdr('Calculating parcel-wise Bayes Factor maps');
            fprintf('\n\n');
           
            N = single(size(cat_obj.dat,2).*(ones(size(parcelwise_stats.t_obj.dat,1),1))); % code from fmri_data.ttest to correct N in regressions_stat.t to make estimateBayesFactor function work on regress() output
            
            for reg = 1:size(parcelwise_stats.t_obj.dat,2)
                
                t_for_Bayes = get_wh_image(parcelwise_stats.t_obj, reg);
                t_for_Bayes.N = N;
                parcelwise_stats.BF(reg) = estimateBayesFactor(t_for_Bayes,'t');
            
            end
            
        end
        
        % ADD DESIGN TABLE
        
        parcelwise_stats.design_table = design_table;

        % ADD CONTRASTNAMES, REGRESSORS, AND OTHER METADATA
        
        switch mygroupnamefield
            case 'contrasts'
                parcelwise_stats.contrastname = DAT.contrastnames{c};
                parcelwise_stats.contrast = DAT.contrasts(c, :);
            case 'conditions'
                parcelwise_stats.contrastname = DAT.conditions{c};
                parcelwise_stats.contrast = 1;
        end

        % ADD VARIABLE NAMES
        
        if ~strcmpi(design_matrix_type,'onesample')
            parcelwise_stats.variable_names = [groupnames {'Intercept'}];
        else
            parcelwise_stats.variable_names = groupnames;
        end
        
        % PLOT PARCELWISE SPECIFIC WEIGHTS AND DIAGNOSTICS
        
        fprintf('\n\n');
        printhdr('Plotting parcel-wise weights and diagnostics');
        fprintf('\n\n');
        
        create_figure('parcelwise weights and metrics', 2, 2);
        plugin_set_figure_size;
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
        
        % PLOT MONTAGE (MASKING ALREADY DONE AS PART OF ROBFIT_PARCELWISE RUN)
        
        fprintf('\n\n');
        printhdr('Plotting parcel-wise GLM results');
        fprintf('\n\n');
        
        fprintf ('\nMONTAGE PARCELWISE GLM RESULTS AT UNCORRECTED p < 0.05, EFFECT: %s, REGRESSOR(S): %s, MASK: %s, SCALING: %s\n\n', parcelwise_stats.contrastname, groupnames_string, mask_string, scaling_string);
        
        num_effects = size(parcelwise_stats.t_obj.dat, 2); % number of regressors
        o2 = canlab_results_fmridisplay([], 'multirow', num_effects);

        for j = 1:num_effects

            tj = get_wh_image(parcelwise_stats.t_obj, j);
            tj = threshold(tj, .05, 'unc'); 

            datsig = tj.dat(logical(tj.sig));
            datsigneg = datsig(datsig<0);
            datsigpos = datsig(datsig>0);

                if isempty(datsigneg) && ~isempty(datsigpos)

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'mincolor',[.9 .4 0], 'maxcolor', [1 1 0]);%, 'cmaprange', [min(datsigpos) max(datsigpos)]);

                elseif isempty(datsigpos) && ~isempty(datsigneg)

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'mincolor',[.1 .8 .8], 'maxcolor', [.1 .1 .8]);%, 'cmaprange', [min(datsigneg) max(datsigneg)]);

                else

                    o2 = addblobs(o2, region(tj), 'wh_montages', (2*j)-1:2*j, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});%, 'cmaprange', [min(datsigneg) max(datsigneg) min(datsigpos) max(datsigpos)]);

                end
                
            if num_effects < 4 % if too many rows, legend gets messy
                o2 = legend(o2);
            end 
            
            o2 = title_montage(o2, 2*j, [parcelwise_stats.contrastname ' ' parcelwise_stats.variable_names{j} ' ' mask_string ' ' scaling_string]);

        end

        figtitle = sprintf('%s_05_unc_montage_%s_%s_%s', parcelwise_stats.contrastname, groupnames_string, mask_string, scaling_string);
        set(gcf, 'Tag', figtitle); plugin_set_figure_size;
        drawnow, snapnow;
            if save_figures_glm
                plugin_save_figure;
            end
        clear o2, clear figtitle, clear j, clear tj
        
        if doBayes
            
            fprintf('\n\n');
            printhdr('Plotting parcel-wise Bayesian GLM results');
            fprintf('\n\n');
           
            fprintf ('\nMONTAGE BAYESIAN PARCELWISE GLM RESULTS AT |BF| > 3, EFFECT: %s, REGRESSOR(S): %s, MASK: %s, SCALING: %s\n\n', parcelwise_stats.contrastname, groupnames_string, mask_string, scaling_string);
        
            o2 = canlab_results_fmridisplay([], 'multirow', num_effects);
            
            for img = 1:size(parcelwise_stats.BF,2)
                
                BF = threshold(parcelwise_stats.BF(1,img),[-2.1972 2.1972],'raw-outside');
                    
                datsig = BF.dat(logical(BF.sig));
                datsigneg = datsig(datsig<0);
                datsigpos = datsig(datsig>0);
                
                    if isempty(datsigneg) && ~isempty(datsigpos)

                        o2 = addblobs(o2, region(BF), 'wh_montages', (2*img)-1:2*img, 'mincolor',[0 0.25 0], 'maxcolor', [0 1 0]);%, 'cmaprange', [min(datsigpos) max(datsigpos)]);

                    elseif isempty(datsigpos) && ~isempty(datsigneg)

                        o2 = addblobs(o2, region(BF), 'wh_montages', (2*img)-1:2*img, 'mincolor',[.25 0 0], 'maxcolor', [1 0 0]);%, 'cmaprange', [min(datsigneg) max(datsigneg)]);

                    else

                        o2 = addblobs(o2, region(BF), 'wh_montages', (2*img)-1:2*img, 'splitcolor',{[.25 0 0] [1 0 0] [0 0.25 0] [0 1 0]});%, 'cmaprange', [min(datsigneg) max(datsigneg) min(datsigpos) max(datsigpos)]); % red in favor of H0, green in favor of H1 for BF maps

                    end

                    if size(BF,2) < 4
                        o2 = legend(o2);
                    end
                
                o2 = title_montage(o2, 2*img, [parcelwise_stats.contrastname ' ' parcelwise_stats.variable_names{img} ' ' mask_string ' ' scaling_string]);
            
            end

            figtitle = sprintf('%s_BF_3_montage_%s_%s_%s', parcelwise_stats.contrastname, groupnames_string, mask_string, scaling_string);
            set(gcf, 'Tag', figtitle); plugin_set_figure_size;
            drawnow, snapnow;
                if save_figures_glm
                    plugin_save_figure;
                end
            clear o2, clear figtitle, clear img, clear BF
            
        end

        % KEEP RESULTS OBJECTS IN CELL ARRAY FOR SAVING

        % robfit_parcelwise takes no analysis_name, so it leaves the CANlab default
        % 'Regression analysis' on every contrast - which is what the parcelwise c2a
        % report then prints as its heading for all of them. The voxelwise branch
        % passes DAT.contrastnames{c} to regress(); do the same here.
        switch mygroupnamefield
            case 'contrasts'
                parcelwise_stats.analysis_name = DAT.contrastnames{c};
            case 'conditions'
                parcelwise_stats.analysis_name = DAT.conditions{c};
        end

        parcelwise_stats_results{c} = parcelwise_stats;

        
    end % if loop voxel- versus parcelwise
    
    
    %%
    % *RUN MVPA MODEL IF REQUESTED IN OPTIONS*
    
    if domvpa_reg_cov
        
        switch mygroupnamefield

            case 'contrasts'
                fprintf('\n\n');
                printhdr(['CONTRAST #', num2str(c), ': ', upper(DAT.contrastnames{c})]);
                fprintf('\n\n');

            case 'conditions'
                fprintf('\n\n');
                printhdr(['CONDITION #', num2str(c), ': ', upper(DAT.conditions{c})]);
                fprintf('\n\n');

        end
        
        fprintf('\n\n');
        printhdr('RUNNING VOXEL-WISE MVPA REGRESSION ANALYSIS');
        fprintf('\n\n');
        
        mvpa_stats_results = cell(kc,size(mvpa_data_objects,2));
        mvpa_dats = cell(kc,size(mvpa_data_objects,2));
        
        for covar = 1:size(mvpa_data_objects,2)
            
            mvpa_dat = mvpa_data_objects{covar};
            
            fprintf('\n\n');
            printhdr(['COVARIATE #', num2str(covar), ': ', upper(mvpa_dat.Y_names)]);
            fprintf('\n\n');
            
            % DATA VISUALIZATION
            % ------------------
            
            fprintf('\n\n');
            printhdr('Plotting X (brain) and Y (behavioural outcome) data');
            fprintf('\n\n');

                % CON IMAGES

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
                title(['Distribution of con weights for ' groupnames{covar}]);
                xlabel('\beta');
                ylabel('Subject');
                hold off

                p = get(gcf,'Position');
                plugin_set_figure_size('width', 6, 'height', 12);   % portrait, was 1024x2048 px
                drawnow, snapnow;

                clear subj

                % BEHAVIORAL OUTCOME

                b1=figure;

                hold off;
                b1=histogram(mvpa_dat.Y);
                box off
                title(['Histogram of ' groupnames{covar}]);
                xlabel(groupnames{covar});
                ylabel('n(observations)');
                plugin_set_figure_size;
                drawnow, snapnow;
            
            % RUN MODEL
            % ---------
            
                % CROSS-VALIDATION FOLD SELECTION
                
                fprintf('\n\n');
                printhdr('Cross-validation fold selection');
                fprintf('\n\n');

                switch holdout_set_method_mvpa_reg_cov

                    case 'no_group'

                        if ~isempty(DAT.BETWEENPERSON.group)
                            fprintf('\n');
                            warning('DAT.BETWEENPERSON.group defines a grouping factor, please change holdout_set_method_mvpa_reg_cov to "group" for correctly stratified CV fold selection.');
                            fprintf('\n');
                        end

                        cv = cvpartition(size(mvpa_dat.dat,2),'KFold',nfolds_mvpa_reg_cov);
                        fold_labels = zeros(size(mvpa_dat.dat,2),1);
                            for subj = 1:cv.NumTestSets
                                fold_labels(cv.test(subj)) = subj;
                            end
                        clear subj

                    case 'group'
                        
                        if ~isempty(DAT.BETWEENPERSON.group)
                            group = DAT.BETWEENPERSON.group;
                        elseif ismember(DAT.BETWEENPERSON.(mygroupnamefield){c}.Properties.VariableNames,group_id{1})
                            group = DAT.BETWEENPERSON.(mygroupnamefield){c}.(group_id{1});
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

                % FIT MODEL
                
                fprintf('\n\n');
                printhdr('Fit MVPA regression model');
                fprintf('\n\n');

                t0 = tic;
                
                switch algorithm_mvpa_reg_cov
                    
                    case 'cv_lassopcr'

                        [mvpa_cverr, mvpa_stats, mvpa_optout, pm] = predict(mvpa_dat, 'algorithm_name', algorithm_mvpa_reg_cov, ...
                            'nfolds', fold_labels, 'error_type', 'mse', 'estimateparams', 'parallel', 'verbose', 0, 'newapi');
                        
                    case 'cv_lassopcrmatlab'

                        [mvpa_cverr, mvpa_stats, mvpa_optout] = predict(mvpa_dat, 'algorithm_name', algorithm_mvpa_reg_cov, ...
                            'nfolds', fold_labels, 'error_type', 'mse', 'EstimateParams', 'parallel', 'verbose', 0);
                                
                    otherwise
                        
                        [mvpa_cverr, mvpa_stats, mvpa_optout] = predict(mvpa_dat, 'algorithm_name', algorithm_mvpa_reg_cov, ...
                            'nfolds', fold_labels, 'error_type', 'mse', 'parallel', 'verbose', 0);
                        
                end

                t_end = toc(t0); 
                
                mvpa_stats.Y_names = mvpa_dat.Y_names;
                mvpa_stats.contrastname = cat_obj.image_names{c};
            
            % VISUALIZE UNTHRESHOLDED RESULTS
            % -------------------------------
            
            fprintf('\n\n');
            printhdr('Plotting MVPA regression results');
            fprintf('\n\n');

                % PLOT OBSERVED VERSUS PREDICTED

                fprintf('\nPLOTTING OBSERVED VERSUS PREDICTED\n');

                fprintf('\n%s r = %0.3f\n\n', algorithm_mvpa_reg_cov, corr(mvpa_stats.yfit, mvpa_dat.Y));
                
                observed = mvpa_dat.Y;
                predicted = mvpa_stats.yfit;
                tbl = table(observed, predicted);
                mdl = fitlm(tbl,'predicted ~ observed','RobustOpts','on');
                
                figure
                
                plot(mdl);
                xlabel({['Observed ' groupnames{covar}]}); ylabel({['Estimated ' groupnames{covar}],'(cross validated)'})

                plugin_set_figure_size;
                drawnow, snapnow;

                % PLOT MONTAGE OF UNTHRESHOLDED WEIGHTS

                fprintf('\nPLOTTTING UNTHRESHOLDED WEIGHT MAPS\n');

                whmontage = 5;

                fprintf ('\nSHOWING UNTHRESHOLDED %s RESULTS, EFFECT: %s, MASK: %s, SCALING: %s\n\n', upper(algorithm_mvpa_reg_cov), mvpa_stats.Y_names, mask_string, myscaling_glm);

                figure

                % canlab_results_fmridisplay's 'compact' layout only calls axes('Position',...);
                % unlike 'multirow' it never opens a figure of its own, so it draws into whatever
                % figure is current - the previous block's montage, or the last figure left open by
                % the previous script in the same session. Open a fresh one. ('multirow' does
                % create its own figure, so those call sites are deliberately left alone: adding
                % figure; there would leave an empty figure behind for every montage.)
                figure;
                o2 = canlab_results_fmridisplay([], 'compact');
                w = mvpa_stats.weight_obj;
                
                w = apply_mask(w,brainmask);
                
                    if exist('maskname_short','var')
                        w = apply_mask(w,glmmask);
                    end
                    
                w = region(w);
                
                o2 = addblobs(o2, w, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
                o2 = legend(o2);
                o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_cov ' unthresholded ' mvpa_stats.Y_names ' ' mask_string ' ' myscaling_glm]);

                figtitle = sprintf('%s_unthresholded_montage_%s_%s', algorithm_mvpa_reg_cov, myscaling_glm, mask_string);
                set(gcf, 'Tag', figtitle); plugin_set_figure_size;
                drawnow, snapnow;

                clear w, clear o2, clear figtitle
                    
            % KEEP RESULTS IN CELL ARRAY FOR SAVING

            mvpa_stats_results{c,covar} = mvpa_stats;
            mvpa_dats{c,covar} = mvpa_dat;

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
        if doBayes
            save(savefilenamedata, 'bayesian_regression_stats_results', '-append');
        end
        if doTFCE
            save(savefilenamedata, 'tfce_regression_stats_results', '-append');
        end
        fprintf('\nSaved regression_stats_results for %s\n', mygroupnamefield);

else
        savefilenamedata = fullfile(resultsdir, ['parcelwise_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
        save(savefilenamedata, 'parcelwise_stats_results', '-v7.3');

        % The parcelwise branch computes Bayes Factors too (see the doBayes block
        % above, which assigns bayesian_regression_stats_results), but this save
        % used to write only parcelwise_stats_results - so the Bayes maps were
        % computed and then silently discarded, and c2a's parcelwise report had no
        % Bayesian section at all. Append them, as the voxelwise save does.
        if doBayes && exist('bayesian_regression_stats_results','var')
            save(savefilenamedata, 'bayesian_regression_stats_results', '-append');
        end

        fprintf('\nSaved parcelwise_stats_results for %s\n', mygroupnamefield);
end

fprintf('\nFilename: %s\n', savefilenamedata);


if doroi_analysis
    
    fprintf('\n\n');
    printhdr('SAVING ROI RESULTS');
    fprintf('\n\n');
    
    savefilenamedata_roi = fullfile(resultsdir, ['roi_stats_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
    save(savefilenamedata_roi, 'roi_means_table', 'roi_adjusted_means', 'roi_glm_stats', 'roi_manova_stats', '-v7.3');
    fprintf('\nSaved roi_stats for %s\n', mygroupnamefield);
    fprintf('\nFilename: %s\n', savefilenamedata_roi);
    
end


if doneurotransmitter_maps
    
    fprintf('\n\n');
    printhdr('SAVING NEUROTRANSMITTER MAP RESULTS');
    fprintf('\n\n');
    
    % FDR ACROSS NEUROTRANSMITTER MAPS
    % ---------------------------------------------------------------------
    % hansen_neurotransmitter_maps tests each map separately and returns one
    % ANOVA table per map. There are ~30 of them, so reading those p-values
    % uncorrected treats thirty chances as one. Correct within contrast,
    % across maps, with the same LaBGAScore_Storey_FDR used by the roi GLM
    % above and by the decoding scripts.
    %
    % Unlike the roi set (8 tests), ~30 maps is enough that Storey can
    % sometimes estimate pi0; when it cannot it falls back to pi0 = 1 and
    % q_Storey equals q_BH, and the verdict below says so.
    %
    % The p-value sits at {2,6} of each ANOVA cell (row 'Groups', column
    % 'Prob>F'); a map whose table is missing or malformed yields NaN and is
    % excluded rather than silently scored.
    neurotransmitter_group_fdr = cell(1, numel(neurotransmitter_group_tables));
    if exist('neurotransmitter_group_tables','var') && ~all(cellfun(@isempty, neurotransmitter_group_tables))
        fprintf('\n\n');
        printhdr('FDR CORRECTION ACROSS NEUROTRANSMITTER MAPS');
        fprintf('\n\n');
        for cnt = 1:numel(neurotransmitter_group_tables)
            Tnt = neurotransmitter_group_tables{cnt};
            if isempty(Tnt), continue, end
            pnt = nan(1, numel(Tnt));
            for mm_i = 1:numel(Tnt)
                x = Tnt{mm_i};
                if iscell(x) && size(x,1) >= 2 && size(x,2) >= 6 && isnumeric(x{2,6}) && isscalar(x{2,6})
                    pnt(mm_i) = x{2,6};
                end
            end
            ok_nt = ~isnan(pnt);
            if ~any(ok_nt)
                fprintf('contrast %d: no usable p-values from the neurotransmitter tables\n', cnt);
                continue
            end
            [q_st_nt, pi0_nt, info_nt] = LaBGAScore_Storey_FDR(pnt(ok_nt));
            q_bh_nt = info_nt.q_BH(:)';
            nmz = {};
            if ~isempty(neurotransmitter_stats) && numel(neurotransmitter_stats) >= cnt && ~isempty(neurotransmitter_stats{cnt})
                if isfield(neurotransmitter_stats{cnt},'networknames') && ~isempty(neurotransmitter_stats{cnt}.networknames)
                    nmz = neurotransmitter_stats{cnt}.networknames;
                end
            end
            idx_nt = find(ok_nt);
            neurotransmitter_group_fdr{cnt} = struct('map_index',idx_nt,'p',pnt(ok_nt), ...
                'q_BH',q_bh_nt,'q_Storey',q_st_nt,'pi0',pi0_nt,'storey_reliable',info_nt.reliable);
            [~, ord_nt] = sort(pnt(ok_nt));
            fprintf('\n%s: %d map(s) tested\n', mygroupnamefield, numel(idx_nt));
            fprintf('  %-24s %10s %9s %10s\n','map','p','q_BH','q_Storey');
            for jj = 1:min(10, numel(ord_nt))
                k = ord_nt(jj); lbl = sprintf('map %d', idx_nt(k));
                if numel(nmz) >= idx_nt(k), lbl = strtrim(char(string(nmz{idx_nt(k)}))); end
                fprintf('  %-24s %10.4f %9.4f %10.4f\n', lbl, pnt(idx_nt(k)), q_bh_nt(k), q_st_nt(k));
            end
            fprintf('\n  %d/%d at q_BH < .05, %d at q_Storey < .05 (pi0 = %.3f)\n', ...
                sum(q_bh_nt < .05), numel(q_bh_nt), sum(q_st_nt < .05), pi0_nt);
            if ~info_nt.reliable
                fprintf('  pi0 not identifiable for this set, so q_Storey equals q_BH - read q_BH\n');
            end
        end
    end

    savefilenamedata_nt = fullfile(resultsdir, ['neurotransmitter_stats_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
    
    if isequal(design_matrix_type,'group') || (isequal(design_matrix_type,'custom') && ~isempty(DAT.BETWEENPERSON.group))
        save(savefilenamedata_nt, 'neurotransmitter_stats', 'neurotransmitter_group_stats', 'neurotransmitter_group_tables', 'neurotransmitter_multcomp_group', 'neurotransmitter_group_fdr', '-v7.3');
    else
        save(savefilenamedata_nt, 'neurotransmitter_stats', '-v7.3');
    end
        
    fprintf('\nSaved neurotransmitter_stats for %s\n', mygroupnamefield);
    fprintf('\nFilename: %s\n', savefilenamedata_nt);
    
end


if domvpa_reg_cov

    fprintf('\n\n');
    printhdr('SAVING MVPA RESULTS');
    fprintf('\n\n');

    savefilenamedata_mvpa = fullfile(resultsdir, ['mvpa_stats_and_maps_', mygroupnamefield, '_', scaling_string, '_', results_suffix, '.mat']);
    save(savefilenamedata_mvpa, 'mvpa_stats_results', 'mvpa_dats','-v7.3');
    fprintf('\nSaved mvpa_stats_results for %s\n', mygroupnamefield);

    fprintf('\nFilename: %s\n', savefilenamedata_mvpa);
    
end