function T = mvpa_reg_cov_benchmark_algorithms(mvpa_dat, fold_labels, varargin)
% mvpa_reg_cov_benchmark_algorithms  Compare @predictive_model algorithms on one dataset.
%
% Re-runs the algorithm comparison under @predictive_model. The earlier
% comparison (cv_pcr 0.73, cv_pls 0.73, cv_svr 0.68) was measured with the
% LEGACY predict() path and a Spider-based SVR, so those numbers do not carry
% over: 'svr' here is @fitrsvm and 'linear_svr' is @fitrlinear.
%
% Every algorithm is given the SAME folds, so differences are the estimator and
% not the split.
%
% :Usage:
% ::
%     T = mvpa_reg_cov_benchmark_algorithms(mvpa_dat, fold_labels);
%     T = mvpa_reg_cov_benchmark_algorithms(mvpa_dat, fold_labels, ...
%             'algorithms', {'pcr','lassopcr_estimateparam','ridge'}, 'seed', 20260923);
%
% :Inputs:
%   **mvpa_dat:**    fmri_data object; .dat voxels x subjects, .Y the outcome.
%   **fold_labels:** n x 1 fold ids - pass the centre-stratified vector, so the
%                    comparison is made under the split actually used.
%
% :Optional Inputs:
%   **'algorithms':** cellstr of keys from the table below.
%   **'seed':**       RNG seed.
%   **'legacy_r':**   the legacy pred_outcome_r for this dataset, if known.
%                     Printed alongside as a regression check.
%
% :Output:
%   **T:** table, one row per algorithm: pearson_r, r2, rmse, fit seconds.
%
% ALGORITHM KEYS
%   'pcr'                      PCA + OLS on components. Unregularised.
%                              NOTE this is also what CANlab's DEFAULT
%                              cv_lassopcr computes - fit_lassopcr with no
%                              modeloptions reduces to PCR.
%   'lassopcr_estimateparam'   PCA + LASSO with the penalty chosen by NESTED
%                              CV + relaxed-OLS refit. The genuinely
%                              regularised CANlab method.
%   'lassopcr_num'             as above, fixed path step via {'lasso_num', k}.
%   'linear_svr'               @fitrlinear. Scales to p >> n.
%   'svr'                      @fitrsvm, linear kernel. Does NOT scale
%                              comfortably to ~150k features - expect it to be
%                              the slow row, and consider dropping it.
%   'ridge', 'lasso'           @fitrlinear with the named regularisation.
%
% HYPERPARAMETERS ARE NOT TUNED HERE, except for
% 'lassopcr_estimateparam', whose nesting is internal. grid_search exists but
% its own docs say nested CV is "not yet automated", so tuning the others on
% this sample would bias their scores upward relative to the nested one. Read
% the untuned rows as lower bounds, not as verdicts.
%
% NOTE ON STATUS: not yet run. See README_predictive_model_port.md.
%
% ..
%     Copyright (C) 2026 Lukas Van Oudenhove. GPLv3.
% ..

p = inputParser;
p.addParameter('algorithms', {'pcr','lassopcr_estimateparam','linear_svr','ridge'}, @iscell);
p.addParameter('seed',       [],  @(x) isempty(x) || isscalar(x));
p.addParameter('legacy_r',   [],  @(x) isempty(x) || isscalar(x));
p.parse(varargin{:});
o = p.Results;

X = double(mvpa_dat.dat)';
Y = mvpa_dat.Y(:);
if numel(fold_labels) ~= numel(Y)
    error('fold_labels has %d entries but %d subject(s).', numel(fold_labels), numel(Y));
end

% One splitter, built once and reused, so every algorithm sees identical folds.
cv = cv_splitter.custom_partition(fold_labels(:));

name = {}; rr = []; r2 = []; rmse = []; secs = [];

for a = 1:numel(o.algorithms)

    key = o.algorithms{a};
    switch key
        case 'pcr',                    alg = 'pcr';        mo = {};
        case 'lassopcr_estimateparam', alg = 'lassopcr';   mo = {'estimateparam'};
        case 'lassopcr_num',           alg = 'lassopcr';   mo = {'lasso_num', 5};
        otherwise,                     alg = key;          mo = {};
    end

    args = {'algorithm', alg, 'task', 'regression', 'use_parallel', true, ...
            'scorer', cv_scorer.pearson_r()};
    if ~isempty(o.seed), args = [args, {'random_state', o.seed}]; end %#ok<AGROW>
    if ~isempty(mo),     args = [args, {'modeloptions', mo}];      end %#ok<AGROW>

    fprintf('\n--- %s ---\n', key);
    t0 = tic;
    try
        pm = predictive_model(args{:});
        pm = crossval(pm, X, Y, 'cv', cv);
        el = toc(t0);

        yfit = pm.fitted_values.yfit(:);
        thisr = corr(yfit, Y);
        sse   = sum((Y - yfit).^2);
        sst   = sum((Y - mean(Y)).^2);

        name{end+1,1} = key;              %#ok<AGROW>
        rr(end+1,1)   = thisr;            %#ok<AGROW>
        r2(end+1,1)   = 1 - sse/sst;      %#ok<AGROW>
        rmse(end+1,1) = sqrt(sse/numel(Y)); %#ok<AGROW>
        secs(end+1,1) = el;               %#ok<AGROW>

        fprintf('  r = %+.4f   R2 = %+.4f   RMSE = %.4f   (%.1f s)\n', ...
                thisr, 1 - sse/sst, sqrt(sse/numel(Y)), el);
    catch ME
        % One algorithm failing must not lose the rest of the comparison.
        fprintf('  FAILED: %s\n', ME.message);
        name{end+1,1} = key;   rr(end+1,1) = NaN;   %#ok<AGROW>
        r2(end+1,1) = NaN;     rmse(end+1,1) = NaN; %#ok<AGROW>
        secs(end+1,1) = toc(t0);                    %#ok<AGROW>
    end
end

T = table(name, rr, r2, rmse, secs, ...
    'VariableNames', {'algorithm','pearson_r','R2','RMSE','seconds'});
T = sortrows(T, 'pearson_r', 'descend', 'MissingPlacement', 'last');

fprintf('\n');
disp(T);

% R2 < 0 means the model predicts worse than the sample mean - the thing that
% matters for a null result, and invisible if only r is reported.
if any(T.R2 < 0)
    fprintf(['\nNOTE: %d algorithm(s) have R2 < 0, i.e. worse than predicting\n' ...
             '      the mean. A positive r with a negative R2 is not prediction.\n'], ...
             sum(T.R2 < 0));
end

if ~isempty(o.legacy_r)
    [~, ib] = max(T.pearson_r);
    fprintf('\nlegacy pred_outcome_r = %+.4f; best here = %+.4f (%s), difference %+.4f\n', ...
        o.legacy_r, T.pearson_r(ib), T.algorithm{ib}, T.pearson_r(ib) - o.legacy_r);
end

end
