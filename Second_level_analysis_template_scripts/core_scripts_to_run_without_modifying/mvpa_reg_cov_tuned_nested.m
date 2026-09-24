function out = mvpa_reg_cov_tuned_nested(mvpa_dat, fold_labels, strata, varargin)
% mvpa_reg_cov_tuned_nested  Nested CV for MVPA-on-covariate, structure-aware.
%
% Implements the nesting pattern CanlabCore's own tutorials teach
% (docs/markdown_tutorials/.../part3, "Nested cross-validation (the honest
% estimate)"), which @predictive_model does NOT do for you:
%
%   outer folds  -> for each: grid_search on the TRAINING rows only, using an
%                   inner splitter built under the SAME structural constraints
%                   -> fit at the chosen hyperparameter -> predict the untouched
%                   outer test fold.
%
% :Usage:
% ::
%     out = mvpa_reg_cov_tuned_nested(mvpa_dat, fold_labels, strata);
%     out = mvpa_reg_cov_tuned_nested(mvpa_dat, fold_labels, strata, ...
%               'algorithm','lassopcr', 'grid', struct('lasso_num',1:12), ...
%               'inner_k', 4, 'nperm', 1000, 'seed', 20260923);
%
% :Inputs:
%   **mvpa_dat:**    fmri_data; .dat voxels x subjects, .Y the outcome.
%   **fold_labels:** n x 1 OUTER fold ids.
%   **strata:**      n x 1 stratification key (e.g. num_center). The inner
%                    folds are rebuilt from the TRAINING SUBSET of this, which
%                    is the whole point - see WHY THIS EXISTS.
%
% :Optional Inputs:
%   **'algorithm':** default 'lassopcr'.
%   **'grid':**      struct of hyperparameter vectors. Default
%                    struct('lasso_num', 1:12) - the PATH STEP, which is what
%                    the tutorials tune, NOT 'estimateparam'.
%   **'inner_k':**   inner fold count, default 4. Chosen INDEPENDENTLY of the
%                    outer count on purpose.
%   **'nperm':**     permutations of the whole nested procedure (0 = skip).
%   **'seed':**      RNG seed.
%
% :Outputs:
%   **out:** struct with .yfit, .r, .r2, .rmse, .chosen (per outer fold),
%            .perm (.p, .null_r, .n), .pm_full (model refit on all data).
%
% WHY THIS EXISTS
% ---------------
% @predictive_model can tune lasso-PCR two ways, and both are unsatisfactory
% for a structured design:
%
%   'estimateparam'  runs an INTERNAL nested CV whose folds are a deterministic
%                    round-robin over ROW INDEX (fit_lassopcr:
%                    cv_assignment = mod(0:n-1,5)'+1). That ignores grouping and
%                    stratification entirely. Under exchangeable rows it is
%                    valid; with grouped data it splits a dependent cluster
%                    across inner folds and leaks; with stratified outer folds
%                    it selects lambda under a different sampling model than the
%                    one being estimated. The function accepts a cv_assignment
%                    argument that would fix this, but NOTHING ever supplies it
%                    - fit.m calls it with three arguments.
%
%   grid_search      does not nest itself; its own documentation says nested CV
%                    is "not yet automated". Wrapping it is left to the caller,
%                    which is what this function does.
%
% Legacy fmri_data/predict takes a third route: it reuses the OUTER partition
% as the inner one (predict.m: options.cv_assignment = cv_assignment). That is
% structure-aware by construction but welds inner k to outer k - degenerate at
% outer k = 2 - and makes the inner partitions near-identical across outer
% folds, so the selected lambdas are strongly dependent.
%
% This function takes the third option: rebuild the inner folds under the SAME
% constraints as the outer ones, with inner k free. Strata are re-derived from
% the training subset each time, so nothing has to be sliced by hand.
%
% NOT USED HERE, DELIBERATELY: select_features. It is applied to the FULL data
% and carried via omitted_features, so calling it before a cross-validation
% leaks the outcome into feature selection. Any feature selection must happen
% inside the fold loop below.
%
% ..
%     Copyright (C) 2026 Lukas Van Oudenhove. GPLv3.
% ..

p = inputParser;
p.addParameter('algorithm', 'lassopcr', @(x) ischar(x) || isstring(x));
p.addParameter('grid',      struct('lasso_num', 1:12), @isstruct);
p.addParameter('inner_k',   4,  @isscalar);
p.addParameter('nperm',     0,  @isscalar);
p.addParameter('seed',      [], @(x) isempty(x) || isscalar(x));
p.addParameter('verbose',   true, @islogical);
p.parse(varargin{:});
o = p.Results;

X = double(mvpa_dat.dat)';
Y = mvpa_dat.Y(:);
n = numel(Y);
if numel(fold_labels) ~= n, error('fold_labels has %d entries, Y has %d.', numel(fold_labels), n); end
if numel(strata)      ~= n, error('strata has %d entries, Y has %d.',      numel(strata), n);      end
if ~isempty(o.seed), rng(o.seed, 'twister'); end

out = local_nested(X, Y, fold_labels, strata, o);

if o.verbose
    fprintf('\nnested CV (%s, tuned on %s): r = %+.4f, R2 = %+.4f, RMSE = %.4f\n', ...
        char(o.algorithm), strjoin(fieldnames(o.grid)', ', '), out.r, out.r2, out.rmse);
    fprintf('hyperparameter chosen per outer fold: %s\n', mat2str(out.chosen(:)'));
    if numel(unique(out.chosen)) > 1
        fprintf(['  NOTE it varies across folds, so "the best value" is itself\n' ...
                 '  uncertain - the tutorials make the same point.\n']);
    end
end

% ---- permutation null over the WHOLE nested procedure -------------------
% Permuting inside the nesting is the only honest null: the tuning is part of
% the procedure being tested, so it has to be redone on every permutation.
% This is why @predictive_model's permutation_test cannot be used here - it
% calls crossval internally and cannot wrap a tuned model.
out.perm = struct('n', 0, 'p', NaN, 'null_r', []);
if o.nperm > 0
    if o.verbose
        fprintf('\npermuting the whole nested procedure %d time(s)\n', o.nperm);
    end
    nullr = nan(o.nperm, 1);
    permidx = zeros(n, o.nperm);
    for pp = 1:o.nperm, permidx(:,pp) = randperm(n)'; end
    parfor pp = 1:o.nperm
        try
            op = o; op.nperm = 0; op.verbose = false;
            rp = local_nested(X, Y(permidx(:,pp)), fold_labels, strata, op);
            nullr(pp) = rp.r;
        catch
            nullr(pp) = NaN;
        end
    end
    nv = nullr(~isnan(nullr));
    out.perm = struct('n', numel(nv), 'p', (sum(nv >= out.r) + 1) / (numel(nv) + 1), ...
                      'null_r', nv, 'n_failed', sum(isnan(nullr)));
    if o.verbose
        fprintf('null: n = %d, mean %+.4f, sd %.4f ; p = %.4f\n', ...
            numel(nv), mean(nv), std(nv), out.perm.p);
        if abs(mean(nv)) > 0.10
            fprintf(['WARNING: null not centred near zero (%+.4f) - check the fold\n' ...
                     '         structure before trusting this p.\n'], mean(nv));
        end
    end
end

% full-data refit at the modal hyperparameter, for the weight map
best = mode(out.chosen);
pmf  = predictive_model('algorithm', char(o.algorithm), 'task', 'regression', ...
                        'modeloptions', local_opts(o.grid, best));
out.pm_full = weight_map_object(fit(pmf, X, Y), mvpa_dat);
out.best_overall = best;

end % main


function res = local_nested(X, Y, fold_labels, strata, o)
% One full outer CV, with an inner grid search per outer training set.
ufold  = unique(fold_labels(:))';
yfit   = nan(numel(Y), 1);
chosen = nan(numel(ufold), 1);
gname  = fieldnames(o.grid); gname = gname{1};
gvals  = o.grid.(gname);

for k = 1:numel(ufold)
    te = fold_labels(:) == ufold(k);
    tr = ~te;

    % INNER folds rebuilt from the TRAINING subset's own strata, under the same
    % constraint as the outer split. This is the step the class's round-robin
    % skips and the reason this function exists.
    inner = local_stratified_folds(strata(tr), o.inner_k);

    % grid search on training rows only
    sc = nan(numel(gvals), 1);
    Xtr = X(tr,:); Ytr = Y(tr);
    for g = 1:numel(gvals)
        yf_in = nan(numel(Ytr), 1);
        for j = unique(inner(:))'
            ite = inner == j; itr = ~ite;
            m = predictive_model('algorithm', char(o.algorithm), 'task', 'regression', ...
                                 'modeloptions', local_opts(o.grid, gvals(g)));
            m = fit(m, Xtr(itr,:), Ytr(itr));
            yf_in(ite) = predict(m, Xtr(ite,:));
        end
        sc(g) = corr(yf_in, Ytr);
    end
    [~, bi]   = max(sc);
    chosen(k) = gvals(bi);

    m = predictive_model('algorithm', char(o.algorithm), 'task', 'regression', ...
                         'modeloptions', local_opts(o.grid, gvals(bi)));
    m = fit(m, Xtr, Ytr);
    yfit(te) = predict(m, X(te,:));
end

sse = sum((Y - yfit).^2); sst = sum((Y - mean(Y)).^2);
res = struct('yfit', yfit, 'chosen', chosen, 'r', corr(yfit, Y), ...
             'r2', 1 - sse/sst, 'rmse', sqrt(sse/numel(Y)));
end


function f = local_stratified_folds(strata, k)
% Stratified k-fold over an arbitrary key: each stratum is dealt round-robin
% across folds, so every fold gets a proportional share of every level.
s = string(strata(:));
f = zeros(numel(s), 1);
for lev = unique(s)'
    idx = find(s == lev);
    idx = idx(randperm(numel(idx)));
    f(idx) = mod(0:numel(idx)-1, k) + 1;
end
end


function c = local_opts(grid, val)
gname = fieldnames(grid); gname = gname{1};
c = {gname, val};
end
