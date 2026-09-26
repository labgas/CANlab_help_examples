function pm = mvpa_reg_cov_stability_from_boot(pm, X, Y, mvpa_dat, varargin)
% mvpa_reg_cov_stability_from_boot  Stability selection from existing bootstrap weights.
%
% :Usage:
% ::
%     pm = mvpa_reg_cov_stability_from_boot(pm, X, Y, mvpa_dat, ...
%              'nstab', 5000, 'k', [], 'threshold', 0.9);
%
% For a predictive_model that has ALREADY been bootstrapped, derive stability
% selection from pm.weights.boot_w instead of refitting.
%
% WHY THIS EXISTS. @predictive_model/stability_selection resamples exactly as
% @predictive_model/bootstrap does - the two blocks are line-for-line identical,
% randi(n,[n,1]) or whole-group sampling, then clone+fit - but it does so in a
% SERIAL `for b = 1:nboot` loop and discards each weight vector after ranking it.
% bootstrap RETAINS every one in pm.weights.boot_w as [p x nboot]. Stability
% selection is then just "how often is each feature in the top-k by |w|", i.e. one
% sort per column of a matrix already in memory. Measured on 93 x 149154 at
% 1.02 s per fit, that is the difference between seconds and ~85 minutes at
% nstab = 5000.
%
% Same estimator, same resampling scheme, same fits - only the bookkeeping differs.
% Falls back to stability_selection() when boot_w is absent or has fewer than
% nstab columns, so it is always safe to call.
%
% :Inputs:
%   **pm:**        predictive_model, already bootstrapped (or not - see fallback)
%   **X, Y:**      the design and outcome the model was fitted on
%   **mvpa_dat:**  fmri_data used to map selection frequencies back into voxels
%
% :Optional Inputs:
%   **'nstab':**     resamples to use; default = all columns of boot_w
%   **'k':**         top-k by |w| per resample; default min(2000, p). See the
%                    c2a header on choosing k WITH the threshold - the
%                    Meinshausen-Buhlmann bound couples them
%   **'threshold':** 'stable' if selected in >= this fraction; default 0.9
%
% :Outputs:
%   **pm:** with pm.diagnostics.stability_selection holding .selection_count,
%           .selection_freq, .stable, .n_stable, .valid_boots, .k, .threshold,
%           .derived_from_bootstrap and .freq_obj (the frequency map)
%
% ..
%     Copyright (C) 2026 Lukas Van Oudenhove. GPLv3.
% ..

p = inputParser;
p.addParameter('nstab',     [],  @(x) isempty(x) || isscalar(x));
p.addParameter('k',         [],  @(x) isempty(x) || isscalar(x));
p.addParameter('threshold', 0.9, @isscalar);
p.parse(varargin{:});
o = p.Results;

if o.threshold <= 0.5
    error('mvpa_reg_cov_stability_from_boot:threshold', ...
        ['threshold must exceed 0.5 - the Meinshausen-Buhlmann error bound is ' ...
         'undefined at or below it (got %.2f).'], o.threshold);
end

p_feat = size(X, 2);
if isempty(o.k), o.k = min(2000, p_feat); end
o.k = max(1, min(o.k, p_feat));

have_boot = isfield(pm.weights,'boot_w') && ~isempty(pm.weights.boot_w);
nstab = o.nstab;
if have_boot
    if isempty(nstab), nstab = size(pm.weights.boot_w, 2); end
    nstab = min(nstab, size(pm.weights.boot_w, 2));
end

if have_boot && nstab > 0

    bw = pm.weights.boot_w(:, 1:nstab);
    bw = bw(:, ~all(isnan(bw), 1));          % the validity test the method uses
    cnt = zeros(p_feat, 1);
    for b = 1:size(bw, 2)
        [~, ord] = sort(abs(bw(:, b)), 'descend');
        cnt(ord(1:o.k)) = cnt(ord(1:o.k)) + 1;
    end

    ss = struct();
    ss.selection_count = cnt;
    ss.valid_boots     = size(bw, 2);
    ss.selection_freq  = cnt / max(ss.valid_boots, 1);
    ss.stable          = ss.selection_freq >= o.threshold;
    ss.n_stable        = sum(ss.stable);
    ss.k               = o.k;
    ss.threshold       = o.threshold;
    ss.derived_from_bootstrap = true;
    pm.diagnostics.stability_selection = ss;

    fprintf(['stability selection derived from %d existing bootstrap weight ' ...
             'vector(s), no refits\n'], ss.valid_boots);

else

    if isempty(nstab) || nstab <= 0
        error('mvpa_reg_cov_stability_from_boot:nstab', ...
            'no bootstrap weights on pm and no usable nstab, so nothing to do.');
    end
    fprintf(['no bootstrap weights on pm - falling back to ' ...
             'stability_selection(), which REFITS %d time(s) serially\n'], nstab);
    pm = stability_selection(pm, X, Y, 'nboot', nstab, 'k', o.k, ...
                             'threshold', o.threshold);

end

% Map the frequencies into voxel space, by the route the method's own docs
% recommend. Done on a COPY so pm.weights.w keeps the real weights - overwriting
% them here would silently corrupt every downstream weight map.
if isfield(pm.diagnostics, 'stability_selection')
    tmp = pm;
    tmp.weights.w = pm.diagnostics.stability_selection.selection_freq(:);
    tmp = weight_map_object(tmp, mvpa_dat);
    pm.diagnostics.stability_selection.freq_obj = tmp.weights.weight_obj;
    ev = o.k^2 / ((2*o.threshold - 1) * p_feat);
    fprintf('  k = %d, pi = %.2f, p = %d -> E(V) <= %.2f (Meinshausen & Buhlmann 2010)\n', ...
            o.k, o.threshold, p_feat, ev);
end

end
