function pm = mvpa_reg_cov_predictive_model(mvpa_dat, fold_labels, varargin)
% mvpa_reg_cov_predictive_model  MVPA regression on a covariate via @predictive_model.
%
% Drop-in replacement for the hand-rolled predict()/permutation/bootstrap code
% in prep_3a's domvpa_reg_cov branch and c2a's MVPA bootstrap branch. Everything
% this function does was previously spread across ~200 lines in two scripts.
%
% :Usage:
% ::
%     pm = mvpa_reg_cov_predictive_model(mvpa_dat, fold_labels);
%     pm = mvpa_reg_cov_predictive_model(mvpa_dat, fold_labels, ...
%              'algorithm', 'pcr', 'nperm', 5000, 'nboot', 5000, 'seed', 20260923);
%
% :Inputs:
%   **mvpa_dat:**   fmri_data object; .dat is voxels x subjects, .Y the outcome.
%   **fold_labels:** n x 1 integer fold ids, one per subject. Pass the vector the
%                   caller already builds for centre stratification - see NOTE.
%
% :Optional Inputs:
%   **'algorithm':**     'lassopcr' (default), 'pcr', 'linear_svr', 'ridge', 'svr'.
%                      NOTE 'lassopcr' with NO modeloptions reduces to plain
%                      PCR - the shrinkage only exists with 'estimateparam'
%                      (nested-CV penalty selection) or {'lasso_num', k}. The
%                      default below therefore passes 'estimateparam'.
%   **'numcomponents':** component count for pcr/lassopcr; [] = registry default.
%   **'nperm':**         permutations for the null (0 = skip).
%   **'nboot':**         bootstrap samples for the weight map (0 = skip).
%   **'nstab':**         resamples for STABILITY SELECTION (0 = skip). This is
%                      an alternative to the bootstrap z/p, not a duplicate of
%                      it - see NOTE ON WHICH PATTERN INFERENCE below.
%   **'stab_k':**        top-k features by |w| per resample. [] = the method's
%                      own default (10% of features, floor 10).
%   **'stab_threshold':** selection frequency at/above which a feature counts
%                      as stable. Default 0.6, per Meinshausen & Buhlmann.
%   **'seed':**          RNG seed; set it or nothing is reproducible.
%   **'use_parallel':**  true (default).
%
% :Outputs:
%   **pm:** fitted @predictive_model, with .weights.weight_obj attached so
%           montage(pm) / surface(pm) work without further arguments.
%
% NOTE ON FOLD STRATIFICATION
% cv_splitter.stratified_kfold stratifies on Y, which is meaningless for a
% CONTINUOUS outcome - it is built for class labels. Centre stratification is
% therefore supplied as a custom partition: the caller computes fold ids from
% num_center exactly as the legacy 'strata' branch does, and they are wrapped
% with cv_splitter.custom_partition. This keeps the one piece of the legacy
% code that encodes a study decision, and delegates everything else.
%
% NOTE ON WHICH PATTERN INFERENCE
% bootstrap() and stability_selection() answer DIFFERENT questions and both are
% offered; neither is a stricter version of the other.
%
%   bootstrap            "is this voxel's weight reliably non-zero?" -> z, p,
%                        FDR-thresholded weight map.
%   stability_selection  "is this voxel reliably among the top-k |weights|?"
%                        -> selection frequency in [0,1], and a stable mask.
%
% The second is the appropriate one for a high-dimensional REGULARISED linear
% model, which is what the default lassopcr + estimateparam is. Regularisation
% pins the weights to nearly the same solution on every resample, so the
% bootstrap z/p collapses and reports implausibly sharp inference; asking
% instead whether a voxel keeps its RANK is what stays informative in that
% regime (Meinshausen & Buhlmann, JRSS-B 2010).
%
% Conversely, for unregularised 'pcr' the bootstrap is the natural choice and
% stability selection is the less motivated one - its rationale is about
% regularised solutions. Running both on the same fit is a reasonable way to
% see whether they agree; disagreement is informative, not an error.
%
% The selection frequencies are mapped into voxel space the same way the
% weights are, so the result is a brain map in [0,1] rather than a bare vector.
%
% NOTE ON STATUS
% Written against the @predictive_model API as it exists in CanlabCore. It has
% NOT yet been run end to end - see README_predictive_model_port.md for the
% verification checklist that must pass before this replaces the legacy path.
%
% ..
%     Author and copyright information:
%     Copyright (C) 2026 Lukas Van Oudenhove
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% ..

% -------------------------------------------------------------------------
% PARSE
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('algorithm',     'lassopcr', @(x) ischar(x) || isstring(x));
p.addParameter('modeloptions',  {'estimateparam'}, @iscell);
p.addParameter('numcomponents', [],     @(x) isempty(x) || isscalar(x));
p.addParameter('nperm',         0,      @isscalar);
p.addParameter('nboot',         0,      @isscalar);
p.addParameter('nstab',         0,      @isscalar);
p.addParameter('stab_k',        [],     @(x) isempty(x) || isscalar(x));
p.addParameter('stab_threshold',0.6,    @isscalar);
p.addParameter('stab_reuse_boot', true, @islogical);   % see the stability block
p.addParameter('seed',          [],     @(x) isempty(x) || isscalar(x));
p.addParameter('use_parallel',  true,   @islogical);
p.parse(varargin{:});
o = p.Results;

% -------------------------------------------------------------------------
% SHAPE THE DATA
% -------------------------------------------------------------------------
% @predictive_model takes numeric X (observations x features), not an image
% object, so transpose out of CANlab's voxels x subjects convention.

X = double(mvpa_dat.dat)';
Y = mvpa_dat.Y(:);

if size(X,1) ~= numel(Y)
    error('mvpa_dat has %d image(s) but %d outcome value(s).', size(X,1), numel(Y));
end
if numel(fold_labels) ~= numel(Y)
    error('fold_labels has %d entries but there are %d subject(s).', ...
          numel(fold_labels), numel(Y));
end
if isempty(o.seed)
    warning(['seed is empty, so the fold split and the permutation null are ' ...
             'NOT reproducible. Set it before reporting anything.']);
end

% -------------------------------------------------------------------------
% BUILD AND FIT
% -------------------------------------------------------------------------

args = {'algorithm', char(o.algorithm), 'task', 'regression', ...
        'use_parallel', o.use_parallel, 'scorer', cv_scorer.pearson_r()};
if ~isempty(o.seed),          args = [args, {'random_state', o.seed}]; end
mo = o.modeloptions;
if ~isempty(o.numcomponents), mo = [mo, {'numcomponents', o.numcomponents}]; end
if ~isempty(mo), args = [args, {'modeloptions', mo}]; end

pm = predictive_model(args{:});

% Centre-stratified folds, computed by the caller - see NOTE above.
cv = cv_splitter.custom_partition(fold_labels(:));

pm = crossval(pm, X, Y, 'cv', cv);

% Inference on the PREDICTION.
if o.nperm > 0
    pm = permutation_test(pm, X, Y, 'nperm', o.nperm);
end

% Inference on the PATTERN. Different question from the permutation test;
% neither substitutes for the other.
if o.nboot > 0
    pm = bootstrap(pm, X, Y, 'nboot', o.nboot);
end

% Stability selection: inference on the RANK of a voxel's weight rather than
% its magnitude. See NOTE ON WHICH PATTERN INFERENCE above.
if o.nstab > 0
    % SPEED: DERIVE STABILITY FROM THE BOOTSTRAP WE ALREADY PAID FOR.
    %
    % @predictive_model/stability_selection refits a model per resample in a
    % SERIAL `for b = 1:nboot` loop. @predictive_model/bootstrap does the same
    % resampling - the two blocks are line-for-line identical, randi(n,[n,1]) or
    % whole-group sampling, then clone+fit - and it RETAINS every weight vector
    % in pm.weights.boot_w as [p x nboot].
    %
    % CAUTION ON THE PARALLELISM. On this machine bootstrap.m runs its loop under
    % `parfor`, but that is an UNCOMMITTED working-tree edit in CanlabCore
    % (for -> parfor, made 2026-07-01; permutation_test.m carries the same edit
    % from 2026-07-13; stability_selection.m never got it). Upstream both are
    % serial. Reuse therefore helps EITHER WAY, and more so upstream: with a
    % parallel bootstrap it converts nstab serial fits into a few sorts, and
    % with a serial bootstrap it still removes nstab fits entirely. Do not rely
    % on the parfor being present - check, or commit it.
    %
    % Stability selection is then just "how often is each feature in the top-k
    % by |w|", which is one sort per column of a matrix already in memory. So
    % when a bootstrap has already run with at least nstab samples, computing it
    % from boot_w is arithmetically the same estimator (same scheme, same fits)
    % at a cost of seconds instead of nstab serial model fits.
    %
    % Set stab_reuse_boot = false to force the method's own refitting path, e.g.
    % to validate this against it.
    did_reuse = false;
    if o.stab_reuse_boot && o.nboot > 0 && isfield(pm.weights,'boot_w') ...
            && ~isempty(pm.weights.boot_w) && size(pm.weights.boot_w,2) >= o.nstab
        bw = pm.weights.boot_w(:, 1:o.nstab);
        okb = ~all(isnan(bw), 1);                 % same validity test the method uses
        bw = bw(:, okb);
        kk = o.stab_k; if isempty(kk), kk = min(2000, size(bw,1)); end
        cnt = zeros(size(bw,1), 1);
        for b = 1:size(bw,2)
            [~, ord] = sort(abs(bw(:,b)), 'descend');
            cnt(ord(1:kk)) = cnt(ord(1:kk)) + 1;
        end
        ss = struct();
        ss.selection_count = cnt;
        ss.valid_boots     = size(bw,2);
        ss.selection_freq  = cnt / max(ss.valid_boots,1);
        ss.stable          = ss.selection_freq >= o.stab_threshold;
        ss.n_stable        = sum(ss.stable);
        ss.k               = kk;
        ss.threshold       = o.stab_threshold;
        ss.derived_from_bootstrap = true;
        pm.diagnostics.stability_selection = ss;
        did_reuse = true;
        fprintf(['\nstability selection derived from the existing %d bootstrap ' ...
                 'weight vectors (no refits)\n'], ss.valid_boots);
    end
    if ~did_reuse
        stabargs = {'nboot', o.nstab, 'threshold', o.stab_threshold};
        if ~isempty(o.stab_k), stabargs = [stabargs, {'k', o.stab_k}]; end
        pm = stability_selection(pm, X, Y, stabargs{:});
    end
end

% Attach the brain map so montage/surface work downstream.
pm = weight_map_object(pm, mvpa_dat);

% Map the selection frequencies into voxel space too, by the route the method's
% own documentation recommends: stash them as a weight vector and re-run
% weight_map_object. Done on a COPY so pm.weights.w keeps the real weights -
% overwriting them here would silently corrupt every downstream weight map.
if o.nstab > 0 && isfield(pm.diagnostics, 'stability_selection')
    freq = pm.diagnostics.stability_selection.selection_freq;
    tmp  = pm;
    tmp.weights.w = freq(:);
    tmp  = weight_map_object(tmp, mvpa_dat);
    pm.diagnostics.stability_selection.freq_obj = tmp.weights.weight_obj;
    fprintf('\nstability selection: %d of %d feature(s) stable at freq >= %.2f\n', ...
        pm.diagnostics.stability_selection.n_stable, numel(freq), o.stab_threshold);
end

end
