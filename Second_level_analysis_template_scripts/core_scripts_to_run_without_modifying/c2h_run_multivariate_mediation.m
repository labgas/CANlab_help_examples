%% c2h_run_multivariate_mediation.m
%
%
% *USAGE*
%
% Single-level multivariate mediation (PDM) on second-level CONTRAST images:
% does a multivariate brain pattern mediate the effect of a between-person
% treatment X on a between-person outcome Y?
%
%   X   group (one value per subject), from DAT.BETWEENPERSON.group
%   M   that subject's contrast image
%   Y   a between-person behavioural outcome, one value per subject
%
% This is the single-level counterpart of c2g_run_multivariate_mediation_single_trial.m.
% That script needs prep_3f to build a single-trial object first, because
% single-trial data is not part of the standard pipeline. Contrast images ARE,
% so this script needs no prep step of its own: it reads the objects prep_3
% already saved and the design prep_1b already built.
%
% Run headless from the Linux command line (the default), which publishes the
% html report and fails loudly if the script errors:
%
%   labgascore_run_headless.sh -d /data/proj_xxx \
%       -s <proj>_secondlevel_m<M>_s0_a_set_up_paths_always_run_first \
%       <proj>_secondlevel_m<M>_s<N>_c2h_run_multivariate_mediation
%
% NOTE: publish() catches a script error into the html and returns normally, so
% a crashed run looks exactly like a successful one. Prefer the route above,
% which reads the report back and checks for a caught error.
%
%
% *DOCUMENTATION*
%
% * PDM toolbox in CANlab's MediationToolbox: https://github.com/canlab/MediationToolbox/tree/master/PDM_toolbox
%       - example: https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
%       - README with papers: https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/README.md
% * multilevel single-trial version: c2g_run_multivariate_mediation_single_trial.m
%
%
% *OPTIONS*
%
% Defaults live in a2_set_default_options for a given model. To run the same
% model with different options, copy this script with a letter index and change
% the options below - AFTER the plugin call, or a re-run of a2 will reset them.
%
% * behav_outcome_pdm       REQUIRED. Name of the outcome variable Y. Looked up in order:
%                           a column of DAT.BETWEENPERSON.contrasts{c}; a FIELD of
%                           DAT.BETWEENPERSON; a column of BIDS/phenotype.csv. The first two
%                           are aligned to image order by construction; only the third relies
%                           on row order (see the alignment note below).
% * myscaling_pdm           default 'scaled'; 'raw', 'scaled' or 'scaled_contrasts'. Selects
%                           DATA_OBJ_CON / DATA_OBJ_CONsc / DATA_OBJ_CONscc, exactly as
%                           myscaling_glm does in prep_3a, so the mediation runs on the same
%                           images as the GLM it accompanies.
% * maskname_pdm            default: maskname_glm if set, else no mask. The mediator is masked
%                           BEFORE the PDM, because here the mask is an analysis decision: it
%                           determines the feature space the components are estimated in.
% * pdm_covs                default {}; names of DAT.BETWEENPERSON.contrasts{c} columns to
%                           residualise out of M (and of Y if pdm_resid_outcome is true)
%                           before the mediation. See the covariate note below.
% * pdm_resid_outcome       default false; also residualise Y on pdm_covs
% * nPDM_svd                default 3; number of PDMs to compute
% * nPDM_B                  default 20; number of PVD components retained before the PDM step
% * dobootstrap_pdm         default true; bootstrap the PDMs
% * boot_n_pdm              default 5000; bootstrap samples (use >=5000 for publication)
% * contrasts2include_pdm   default []; empty = every contrast in DAT.contrasts
% * dosavepdmstats          default true
%
%
% *NOTES*
%
% COVARIATES. multivariateMediation takes only X, Y and M - the PDM framework
% has no covariate facility. pdm_covs therefore residualises the covariate out
% of the mediator (and optionally the outcome) beforehand. That is NOT the same
% as modelling it: it removes variance shared between covariate and mediator,
% including any that the treatment also explains, so an X collinear with the
% covariate will lose real signal. Where a site or scanner term is confounded
% with group by design, no amount of residualising repairs it - restrict the
% sample instead.
%
% OUTCOME ALIGNMENT. Y must be in the same subject order as the columns of the
% contrast object. A column of DAT.BETWEENPERSON.contrasts{c} is aligned by
% construction, so that route is preferred and needs no assumption. The
% phenotype.csv fallback assumes the file's row order matches the image order,
% which is what prep_1b relies on; the script checks the counts and refuses to
% proceed if they differ, but it cannot detect a re-ordered file. If you can,
% add the outcome to prep_1b so it travels in DAT.BETWEENPERSON.
%
% MISSING OUTCOMES. Subjects with NaN Y are dropped, and the number dropped is
% reported. The mediation sample is then smaller than the GLM sample, so the
% two are not describing quite the same subjects - say so when reporting.
%
% -------------------------------------------------------------------------
% Author: Lukas Van Oudenhove
% Date:   Leuven, September, 2026
% -------------------------------------------------------------------------
% c2h_run_multivariate_mediation.m         v1.0
% -------------------------------------------------------------------------


%% RUN SCRIPT A_SET_UP_PATHS_ALWAYS_RUN_FIRST AND LOAD/CREATE DAT IF NEEDED
% -------------------------------------------------------------------------

resultsdir_before_setup = '';
if exist('resultsdir','var'), resultsdir_before_setup = resultsdir; end

a_set_up_paths_always_run_first;

if ~isempty(resultsdir_before_setup) && ~isequal(resultsdir_before_setup, resultsdir)
    error(['\nThe generic a_set_up_paths_always_run_first has overwritten resultsdir:\n' ...
           '  before: %s\n  after : %s\n' ...
           'Run this model''s own s0 script instead.\n'], resultsdir_before_setup, resultsdir);
end

if ~exist('DAT','var')
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
end


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

options_needed = {'myscaling_glm','maskname_glm'};
options_exist  = cellfun(@exist, options_needed);
option_default_values = {'scaled', which('gm_mask_canlab2023_coarse_fmriprep20_0_20.nii')};
plugin_get_options_for_analysis_script;

if ~exist('behav_outcome_pdm','var') || isempty(behav_outcome_pdm)
    error('behav_outcome_pdm is not set: name the outcome variable Y.');
end
if ~exist('myscaling_pdm','var') || isempty(myscaling_pdm),   myscaling_pdm = 'scaled'; end
if ~exist('maskname_pdm','var'),                              maskname_pdm = maskname_glm; end
if ~exist('pdm_covs','var'),                                  pdm_covs = {}; end
if ~exist('pdm_resid_outcome','var') || isempty(pdm_resid_outcome), pdm_resid_outcome = false; end
if ~exist('dobootsamples_pdm','var') || isempty(dobootsamples_pdm), dobootsamples_pdm = false; end
if ~exist('nPDM_svd','var') || isempty(nPDM_svd),             nPDM_svd = 3; end
if ~exist('nPDM_B','var') || isempty(nPDM_B),                 nPDM_B = 20; end
if ~exist('dobootstrap_pdm','var') || isempty(dobootstrap_pdm), dobootstrap_pdm = true; end
if ~exist('boot_n_pdm','var') || isempty(boot_n_pdm),         boot_n_pdm = 5000; end
if ~exist('contrasts2include_pdm','var'),                     contrasts2include_pdm = []; end
if ~exist('dosavepdmstats','var') || isempty(dosavepdmstats), dosavepdmstats = true; end
if ischar(pdm_covs) || isstring(pdm_covs), pdm_covs = cellstr(pdm_covs); end

if isempty(which('multivariateMediation'))
    error('multivariateMediation not found. Add MediationToolbox/PDM_toolbox to the path.');
end

fprintf('\n\n');
printhdr('SINGLE-LEVEL MULTIVARIATE MEDIATION (PDM)');
fprintf('\n\n');
fprintf('  outcome (Y)     : %s\n', behav_outcome_pdm);
fprintf('  mediator scaling: %s\n', myscaling_pdm);
fprintf('  PDMs            : %d (from %d PVD components)\n', nPDM_svd, nPDM_B);
if dobootstrap_pdm, fprintf('  bootstrap       : %d samples\n', boot_n_pdm);
else,               fprintf('  bootstrap       : off\n'); end
if isempty(pdm_covs), fprintf('  residualised on : none\n');
else,                 fprintf('  residualised on : %s%s\n', strjoin(pdm_covs,', '), ...
                              string(pdm_resid_outcome).replace("true"," (mediator and outcome)").replace("false"," (mediator only)")); end


%% LOAD CONTRAST OBJECTS
% -------------------------------------------------------------------------

if ~exist('DATA_OBJ_CON','var')
    load(fullfile(resultsdir,'contrast_data_objects.mat'));
end

switch myscaling_pdm
    case 'raw'
        OBJ = DATA_OBJ_CON;      scaling_string_pdm = 'no_scaling';
    case 'scaled'
        OBJ = DATA_OBJ_CONsc;    scaling_string_pdm = 'scaling_z_score_conditions';
    case 'scaled_contrasts'
        OBJ = DATA_OBJ_CONscc;   scaling_string_pdm = 'scaling_l2norm_contrasts';
    otherwise
        error('invalid myscaling_pdm "%s": choose raw, scaled or scaled_contrasts', myscaling_pdm);
end

cons = contrasts2include_pdm;
if isempty(cons), cons = 1:size(DAT.contrasts,1); end

mediationresultsdir = fullfile(resultsdir,'mediation_analysis','pdm');
if ~exist(mediationresultsdir,'dir'), mkdir(mediationresultsdir); end


%% RUN THE MEDIATION, ONE CONTRAST AT A TIME
% -------------------------------------------------------------------------

pdm_results = cell(1, size(DAT.contrasts,1));

for c = cons

    fprintf('\n\n');
    printhdr(sprintf('CONTRAST #%d: %s', c, upper(DAT.contrastnames{c})));
    fprintf('\n\n');

    obj = OBJ{c};
    n_img = size(obj.dat,2);

    % ---- X: treatment ----------------------------------------------------
    X_all = DAT.BETWEENPERSON.group(:);
    if isempty(X_all)
        error('DAT.BETWEENPERSON.group is empty; this script mediates a GROUP effect.');
    end
    if numel(X_all) ~= n_img
        error('group has %d entries but contrast %d has %d images.', numel(X_all), c, n_img);
    end

    % ---- Y: outcome ------------------------------------------------------
    T = [];
    if isfield(DAT.BETWEENPERSON,'contrasts') && numel(DAT.BETWEENPERSON.contrasts) >= c
        T = DAT.BETWEENPERSON.contrasts{c};
    end
    if istable(T) && ismember(behav_outcome_pdm, T.Properties.VariableNames)
        Y_all = double(T.(behav_outcome_pdm)(:));
        fprintf('  Y taken from DAT.BETWEENPERSON.contrasts{%d} (aligned by construction)\n', c);
    elseif isstruct(DAT.BETWEENPERSON) && isfield(DAT.BETWEENPERSON, behav_outcome_pdm)
        % Carried by prep_1b as its own field rather than a covs column: the covs
        % table becomes the second-level design matrix, so an outcome placed
        % there would silently enter the GLM. A field is aligned to image order
        % by the same construction, without that side effect - so this route is
        % as safe as the table one and avoids the phenotype.csv fallback.
        Y_all = double(DAT.BETWEENPERSON.(behav_outcome_pdm)(:));
        fprintf('  Y taken from DAT.BETWEENPERSON.%s (aligned by construction)\n', behav_outcome_pdm);
    else
        phenofile = fullfile(BIDSdir,'phenotype.csv');
        if exist(phenofile,'file') ~= 2
            error('"%s" is not in the design table and %s does not exist.', behav_outcome_pdm, phenofile);
        end
        P = readtable(phenofile,'FileType','text','Delimiter',',');
        if ~ismember(behav_outcome_pdm, P.Properties.VariableNames)
            error('"%s" is in neither the design table nor %s.', behav_outcome_pdm, phenofile);
        end
        if height(P) ~= n_img
            error(['%s has %d rows but contrast %d has %d images, so the two cannot be\n' ...
                   'matched by order. Add %s to prep_1b instead.'], ...
                   phenofile, height(P), c, n_img, behav_outcome_pdm);
        end
        Y_all = double(P.(behav_outcome_pdm));
        fprintf(['  Y read from %s by ROW ORDER - this assumes the file is in the same order\n' ...
                 '  as the images, which is what prep_1b relies on but cannot be verified here\n'], phenofile);
    end

    % ---- drop subjects without an outcome --------------------------------
    ok = ~isnan(Y_all) & ~isnan(X_all);
    if ~all(ok)
        fprintf('  %d of %d subject(s) dropped for a missing outcome; n = %d for the mediation\n', ...
            sum(~ok), numel(ok), sum(ok));
    end
    if sum(ok) < 20
        error('only %d subjects have both X and Y; too few for a stable PDM.', sum(ok));
    end

    % Subset with get_wh_image, not by slicing .dat: the object carries
    % per-image bookkeeping (removed_images, image_names) that a raw slice
    % leaves at the original length, and apply_mask then fails inside
    % remove_empty with "Arrays have incompatible sizes".
    obj_c = get_wh_image(obj, find(ok));
    X = X_all(ok);  Y = Y_all(ok);

    % ---- mask the mediator ----------------------------------------------
    if ~isempty(maskname_pdm) && exist(maskname_pdm,'file')
        mk = fmri_mask_image(maskname_pdm,'noverbose');
        if any(unique(mk.dat) ~= 1), mk.dat(mk.dat > 0) = 1; end
        if ~isequal(abs(diag(mk.volInfo.mat(1:3,1:3)))', abs(diag(obj_c.volInfo.mat(1:3,1:3)))')
            mk = resample_space(mk, obj_c);
            mk.dat(mk.dat < 1) = 0;
        end
        nv_before = size(obj_c.dat,1);
        obj_c = apply_mask(obj_c, mk);
        [~,mask_short] = fileparts(maskname_pdm);
        fprintf('  mediator masked with %s: %d of %d voxels\n', mask_short, size(obj_c.dat,1), nv_before);
    else
        mask_short = 'none';
        fprintf('  mediator not masked\n');
    end

    % ---- optional residualisation of covariates --------------------------
    if ~isempty(pdm_covs)
        if ~istable(T)
            error('pdm_covs is set but DAT.BETWEENPERSON.contrasts{%d} is not a table.', c);
        end
        miss = pdm_covs(~ismember(pdm_covs, T.Properties.VariableNames));
        if ~isempty(miss)
            error('pdm_covs names %s, not in the design table.', strjoin(miss,', '));
        end
        Cov = [];
        for k = 1:numel(pdm_covs)
            Cov = [Cov double(T.(pdm_covs{k})(:))]; %#ok<AGROW>
        end
        Cov = Cov(ok,:);
        Xd  = [ones(size(Cov,1),1) Cov];
        Mdat = double(obj_c.dat)';                 % subjects x voxels
        var_before = mean(var(Mdat,0,1));
        Mdat = Mdat - Xd*(Xd\Mdat);
        obj_c.dat = Mdat';
        fprintf('  mediator residualised on %s: mean voxel variance %.4g -> %.4g\n', ...
            strjoin(pdm_covs,', '), var_before, mean(var(Mdat,0,1)));
        if pdm_resid_outcome
            Y = Y - Xd*(Xd\Y);
            fprintf('  outcome also residualised on the same covariates\n');
        end
        fprintf(['  NOTE residualising is not the same as modelling the covariate: variance\n' ...
                 '  shared with the treatment goes with it, so a treatment collinear with the\n' ...
                 '  covariate loses real signal.\n']);
    end

    % ---- assemble in the cell format multivariateMediation expects -------
    % One cell per subject, each holding that subject's scalar x, scalar y and
    % 1 x nvox mediator. This is the same layout the multilevel version uses,
    % with a single observation per subject.
    n = numel(X);
    xx = num2cell(X);
    yy = num2cell(Y);
    mm = cell(n,1);
    Md = double(obj_c.dat);
    for i = 1:n, mm{i} = Md(:,i); end

    fprintf('\n  running multivariateMediation: n = %d, %d voxels\n', n, size(Md,1));

    args = {xx, yy, mm, 'B', nPDM_B, 'svd', nPDM_svd};
    if dobootstrap_pdm
        args = [args {'bootPDM', 1:nPDM_svd, 'bootJPDM', 'Bsamp', boot_n_pdm}]; %#ok<AGROW>
        if dobootsamples_pdm
            args = [args {'returnbootsamples'}]; %#ok<AGROW>
        end
    end
    pdm = multivariateMediation(args{:});
    pdm_results{c} = pdm;

    % PATH-COEFFICIENT INFERENCE IS DELIBERATELY NOT REPORTED.
    % An earlier version bootstrapped CIs and p-values for c, c', a, b and ab
    % from pdm.boot.SamplesTheta. That was removed: Chen et al. (2018)
    % Biostatistics 19(2):121 base inference on the voxel weights |w_k| alone,
    % via a pseudo-null and half-normal fit, precisely because the signs of the
    % DMs are unidentifiable - a sign flip in w_k is offset by sign flips in
    % both alpha_k and beta_k. The paper reports no CIs or p-values for
    % alpha_k*beta_k and no joint test across DMs.
    %
    % In practice the implementation resolves the sign indeterminacy with a
    % positive-alpha convention, so every bootstrap sample returned a > 0 and
    % the a-path p-value pinned to its floor, 2/(Bsamp+1), in every PDM. That
    % is a property of the software, not evidence. The product alpha*beta is
    % invariant to the sign flip and so is better behaved in principle, but it
    % is still exposed to direction correspondence across bootstrap samples
    % (the k-th DM in a resample need not be the k-th DM of the original fit),
    % which is what makes a point estimate fall outside its own CI.
    %
    % SamplesTheta is still retained when dobootsamples_pdm is true - it costs
    % ~0.5 MB against SamplesW's many GB, and having it is what allowed the
    % above to be checked rather than assumed. It is simply not reported.

    % ---- report and write out -------------------------------------------
    dat_template = obj_c;
    dat_template.dat = zeros(size(obj_c.dat,1),1);
    try
        [~, figh] = plotPDM(pdm, dat_template); %#ok<ASGLU>
        drawnow, snapnow
    catch ME
        fprintf('  plotPDM failed (%s); continuing\n', ME.message);
    end

    % The covariate set is part of the analysis identity for the written maps
    % exactly as it is for the .mat: without it a residualised run overwrites
    % the unadjusted PDM*.nii of the same contrast, silently, and the two
    % analyses end up sharing one set of images. covtag_pdm is built here
    % rather than reused from the save block because that block runs later.
    if isempty(pdm_covs)
        covtag_dir = '';
    else
        covtag_dir = ['_adj_' strjoin(pdm_covs, '_')];
        if pdm_resid_outcome, covtag_dir = [covtag_dir '_yadj']; end
    end
    condir = fullfile(mediationresultsdir, ...
        [matlab.lang.makeValidName(DAT.contrastnames{c}) covtag_dir]);
    if ~exist(condir,'dir'), mkdir(condir); end

    if isfield(pdm,'boot') && isfield(pdm.boot,'p')
        for k = 1:numel(pdm.boot.p)
            d = dat_template;
            d.dat = pdm.Wfull{k} .* (pdm.boot.p{k} < pdm.pThreshold(k));
            nsig = sum(d.dat ~= 0);
            fprintf('  PDM%d: %d voxel(s) below the bootstrap threshold p < %.4g\n', k, nsig, pdm.pThreshold(k));
            write(d, 'fname', fullfile(condir, sprintf('PDM%d.nii', k)), 'overwrite');
        end
    else
        fprintf('  no bootstrap performed, so no thresholded PDM images written\n');
    end

end

%% SAVE
% -------------------------------------------------------------------------

if dosavepdmstats
    % The covariate set is part of the analysis identity: without it a
    % residualised run silently overwrites the unadjusted one of the same
    % outcome and scaling.
    if isempty(pdm_covs)
        covtag_pdm = '';
    else
        covtag_pdm = ['_adj_' strjoin(pdm_covs, '_')];
        if pdm_resid_outcome, covtag_pdm = [covtag_pdm '_yadj']; end
    end
    savefilename_pdm = fullfile(mediationresultsdir, ...
        ['pdm_mediation_', behav_outcome_pdm, '_', scaling_string_pdm, covtag_pdm, '.mat']);
    save(savefilename_pdm, 'pdm_results', 'behav_outcome_pdm', 'myscaling_pdm', ...
         'pdm_covs', 'pdm_resid_outcome', 'nPDM_svd', 'nPDM_B', 'boot_n_pdm', '-v7.3');
    fprintf('\nSaved: %s\n', savefilename_pdm);
end
