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
%   **out:** struct with
%              .r_pooled     correlation of ALL held-out predictions with Y.
%                            THIS is the number comparable to the legacy and
%                            @predictive_model engines, which pool the same way.
%              .r_per_fold   correlation within each outer fold. Very noisy at
%                            n/outer_k ~ 19: SE(r) ~ 0.25 against ~0.11 pooled.
%                            Read it as a diagnostic of fold stability, not as a
%                            result.
%              .r_mean_fold  plain mean of .r_per_fold. Reported only because
%                            crossValScore would give this; prefer .r_fisher.
%              .r_fisher     (n_k - 3)-weighted mean in Fisher z, the defensible
%                            way to average correlations.
%              .yfit         the held-out predictions themselves.
%              .cvP          the fitted crossValPredict object (carries .cvpart).
%              .estimator    the gridSearchCV estimator.
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

if isempty(which('crossValPredict'))
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

gs = gridSearchCV(est, grid, innercv, @get_r);

% ONE PASS, NOT TWO. An earlier version ran crossValScore (for the per-fold
% losses) AND crossValPredict (for the pooled predictions) with the same
% estimator and partitioner, i.e. it executed the entire nested CV TWICE - 24.7
% minutes serially where half that was enough. crossValPredict exposes both
% .yfit and .cvpart, so every quantity below is derivable from it alone.
%
% NOTE on get_r, which the earlier version also misread: it returns 1 - r, a
% LOSS to be MINIMISED, which is correct for gridSearchCV. Reading cvGS.scores
% as correlations gave an impossible "+0.9148" with per-fold values of 1.1888
% and 1.0661. Nothing below depends on that convention any more - the
% correlations here are computed directly from the held-out predictions.
cvP = crossValPredict(gs, outercv, 'n_parallel', o.n_parallel, 'verbose', true);
cvP = cvP.do(dat, dat.Y);

Y    = dat.Y(:);
yfit = cvP.yfit(:);
okp  = ~isnan(yfit) & ~isnan(Y);

out = struct();
out.cvP       = cvP;
out.estimator = gs;
out.yfit      = yfit;

% POOLED r - the statistic comparable to the legacy and @predictive_model
% engines, which correlate all held-out predictions at once.
out.r_pooled = corr(yfit(okp), Y(okp));

% PER-FOLD r, from the same predictions, plus the true fold sizes for the
% Fisher weights below.
nfold = cvP.cvpart.NumTestSets;
out.r_per_fold = nan(1, nfold);
nk             = nan(1, nfold);
for k = 1:nfold
    te    = cvP.cvpart.test(k) & okp;
    nk(k) = sum(te);
    if nk(k) > 2 && std(Y(te)) > 0 && std(yfit(te)) > 0
        out.r_per_fold(k) = corr(yfit(te), Y(te));
    end
end
out.r_mean_fold = mean(out.r_per_fold, 'omitnan');

% CORRELATION IS NOT AN ADDITIVE LOSS, so "average over folds" and "compute over
% the pooled predictions" are genuinely different numbers - for MSE or for R2
% against a global mean the two are algebraically identical once folds are
% weighted by n_k, and no such choice arises. Where an average IS wanted, average
% in Fisher z, never in raw r: r is biased for rho (bias ~ -rho(1-rho^2)/2n) and
% skewed, whereas z = atanh(r) is near-normal with variance ~ 1/(n-3)
% INDEPENDENT OF RHO, so (n_k - 3) are the right weights. nk is now measured
% rather than assumed equal.
out.r_fisher = NaN;
z  = atanh(min(max(out.r_per_fold, -0.999999), 0.999999));
w  = max(nk - 3, 0);
ok = isfinite(z) & w > 0;
if any(ok), out.r_fisher = tanh(sum(w(ok) .* z(ok)) / sum(w(ok))); end

fprintf('\nooFmriDataObjML nested CV\n');
fprintf('  per-fold r      : %s\n', mat2str(round(out.r_per_fold, 4)));
fprintf('  mean per-fold r : %+.4f   (NOT comparable to the pooled r of other engines)\n', ...
        out.r_mean_fold);
fprintf('  Fisher-z mean r : %+.4f   (the defensible version of the line above)\n', out.r_fisher);
fprintf('  POOLED r        : %+.4f   <- use this for engine comparison\n', out.r_pooled);
fprintf(['  NOTE: pooled r mixes predictions from %d different models. It is only\n' ...
         '  honest if fold membership does not itself predict Y - check the\n' ...
         '  stratification before reading it.\n'], numel(out.r_per_fold));

end
