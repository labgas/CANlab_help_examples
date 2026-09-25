function out = mvpa_reg_cov_oofmri(mvpa_dat, strata, varargin)
% mvpa_reg_cov_oofmri  MVPA-on-covariate via ooFmriDataObjML's nested CV.
%
% Third engine for the domvpa_reg_cov analysis, alongside the legacy predict()
% path and the @predictive_model path. Exists because ooFmriDataObjML solves
% the nesting problem architecturally rather than by convention.
%
% :Usage:
% ::
%     out = mvpa_reg_cov_oofmri(mvpa_dat, strata);
%     out = mvpa_reg_cov_oofmri(mvpa_dat, strata, 'outer_k',5, 'inner_k',4, ...
%               'grid', (1:12)', 'seed', 20260923);
%
% :Inputs:
%   **mvpa_dat:** fmri_data; .dat voxels x subjects, .Y the outcome.
%   **strata:**   n x 1 stratification key (e.g. num_center), attached to the
%                 object's metadata_table so the partitioners can re-derive
%                 folds from any subset - see THE ARCHITECTURAL POINT.
%
% :Optional Inputs:
%   **'outer_k' / 'inner_k':** fold counts, default 5 and 4. Independent on
%                 purpose: the legacy predict() path welds them together.
%   **'grid':**   column vector of numcomponents values to search.
%   **'seed':**   RNG seed.
%
% :Outputs:
%   **out:** struct with .scores (per outer fold), .r, .cvGS (the fitted
%            crossValScore object), .estimator.
%
% THE ARCHITECTURAL POINT
% -----------------------
% The partitioners here are FUNCTION HANDLES, not fold vectors:
%
%     innercv = @(X,Y) cvpartition2(X.metadata_table.strata, 'KFold', inner_k);
%
% so when the outer loop hands a training SUBSET down, the inner partitioner
% re-derives its folds from that subset's own metadata. Stratification (and
% grouping, via 'Group') therefore propagates to every nesting level by
% construction, with nothing sliced by hand.
%
% That is the difference from the other two engines:
%
%   predict()           inner folds = the OUTER partition reused. Structure
%                       aware, but inner k is welded to outer k and the inner
%                       partitions are near-identical across outer folds.
%   @predictive_model   'estimateparam' uses a round-robin over ROW INDEX,
%                       which ignores structure entirely; grid_search does not
%                       nest itself, so the caller must slice labels by hand.
%   this                partitioner re-derived per level, inner k free.
%
% A second advantage is composability: gridSearchCV(estimator, grid, innercv)
% RETURNS AN ESTIMATOR, so crossValScore(gs, outercv, scorer) nests by
% composition and the outer loop is library code. In @predictive_model,
% grid_search returns a model with fixed options and there is no tuned-estimator
% object, so permutation_test and bootstrap - which call crossval internally -
% cannot wrap a tuned model at all.
%
% CAVEAT: ooFmriDataObjML is UNMAINTAINED (last commit 2024-08-23). It is used
% by the fork's prep_3c and c2f, so it is already a dependency, but it is not
% actively developed. Treat this engine as a reference implementation and a
% cross-check on the other two, not as the default.
%
% ..
%     Copyright (C) 2026 Lukas Van Oudenhove. GPLv3.
% ..

p = inputParser;
p.addParameter('outer_k', 5,   @isscalar);
p.addParameter('inner_k', 4,   @isscalar);
p.addParameter('grid',    (1:12)', @isnumeric);
p.addParameter('seed',    [],  @(x) isempty(x) || isscalar(x));
p.addParameter('n_parallel', 5, @isscalar);
p.parse(varargin{:});
o = p.Results;

% PATH ORDER MATTERS HERE. CanlabCore defines its OWN @pipeline class, and
% addpath PREPENDS, so whichever of the two is added LAST wins. ooFmriDataObjML
% must come after CanlabCore or `pipeline` resolves to CanlabCore's and fails
% with an unrelated error inside normalize_step.
pw = which('pipeline');
if isempty(strfind(pw, 'ooFmriDataObjML'))
    error(['`pipeline` resolves to %s, not ooFmriDataObjML''s.\n' ...
           'addpath(genpath(''/data/master_github_repos/ooFmriDataObjML'')) ' ...
           'AFTER CanlabCore.'], pw);
end

if isempty(which('crossValScore'))
    error(['ooFmriDataObjML is not on the path. Add it with\n' ...
           '  addpath(genpath(''/data/master_github_repos/ooFmriDataObjML''))\n']);
end
if ~isempty(o.seed), rng(o.seed, 'twister'); end

% Attach the stratification key to the object's metadata so the partitioners
% can re-derive folds from any subset. The library's contract is that any
% metadata_table field referenced by a partitioner must survive cat() and
% get_wh_image(), which a plain column does.
dat = mvpa_dat;
if isempty(dat.metadata_table) || height(dat.metadata_table) ~= numel(dat.Y)
    dat.metadata_table = table();
end
dat.metadata_table.strata = strata(:);

% estimator: voxels -> features -> PCR, with numcomponents as the tunable
est = pipeline({{'featurizer', fmri2VxlFeatTransformer()}, ...
                {'model',      pcrRegressor()}});

% The grid column must be named as pipeline.get_params() reports it:
% STEP-PREFIXED with a double underscore, i.e. 'model__numcomponents' for the
% step named 'model'. A bare 'numcomponents' fails inside gridSearchCV with
%   optimizableVariable names must match pipeline.get_params()
% Verified by calling get_params() on the constructed pipeline rather than
% guessing; rename the step and this name changes with it.
gridname = 'model__numcomponents';
gp = est.get_params();
if ~ismember(gridname, gp)
    error(['grid name %s is not in pipeline.get_params(): %s\n' ...
           'Name the grid column exactly as get_params() reports it.'], ...
           gridname, strjoin(gp, ', '));
end
grid = table(o.grid(:), 'VariableNames', {gridname});

innercv = @(X,Y) cvpartition2(X.metadata_table.strata, 'KFold', o.inner_k);
outercv = @(X,Y) cvpartition2(X.metadata_table.strata, 'KFold', o.outer_k);

gs   = gridSearchCV(est, grid, innercv, @get_r);
cvGS = crossValScore(gs, outercv, @get_r, 'n_parallel', o.n_parallel, 'verbose', true);
cvGS = cvGS.do(dat, dat.Y);

out = struct();
out.cvGS      = cvGS;
out.estimator = gs;
out.scores    = cvGS.scores;
out.r         = mean(cvGS.scores, 'omitnan');

fprintf('\nooFmriDataObjML nested CV: mean r over %d outer fold(s) = %+.4f\n', ...
        numel(out.scores), out.r);
fprintf('per-fold r: %s\n', mat2str(round(out.scores(:)', 4)));

end
