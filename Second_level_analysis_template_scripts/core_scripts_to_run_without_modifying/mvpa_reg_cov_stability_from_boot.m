function ss = mvpa_reg_cov_stability_from_boot(pm, X, Y, mvpa_dat, varargin)
% mvpa_reg_cov_stability_from_boot  Stability selection from existing bootstrap weights.
%
% :Usage:
% ::
%     ss = mvpa_reg_cov_stability_from_boot(pm, X, Y, mvpa_dat, ...
%              'nstab', 5000, 'k', [], 'threshold', 0.9);
%
% RETURNS A STRUCT, AND DOES NOT MUTATE pm. predictive_model declares
% diagnostics, weights and most other fitted state under
% properties (SetAccess = protected), so only class methods may write them - an
% external function assigning pm.diagnostics.stability_selection fails with
% "Unable to set the 'diagnostics' property ... because it is read-only".
% The frequency map is therefore built by COPYING pm.weights.weight_obj, which is
% already a statistic_image in the right space, and replacing its .dat.
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
%   **ss:** struct with .selection_count, .selection_freq, .stable, .n_stable,
%           .valid_boots, .k, .threshold, .derived_from_bootstrap and .freq_obj
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

    fprintf(['stability selection derived from %d existing bootstrap weight ' ...
             'vector(s), no refits\n'], ss.valid_boots);

else

    if isempty(nstab) || nstab <= 0
        error('mvpa_reg_cov_stability_from_boot:nstab', ...
            'no bootstrap weights on pm and no usable nstab, so nothing to do.');
    end
    fprintf(['no bootstrap weights on pm - falling back to ' ...
             'stability_selection(), which REFITS %d time(s) serially\n'], nstab);
    % The class METHOD may write diagnostics; read the result back out.
    pm = stability_selection(pm, X, Y, 'nboot', nstab, 'k', o.k, ...
                             'threshold', o.threshold);
    ss = pm.diagnostics.stability_selection;
    ss.derived_from_bootstrap = false;

end

% Frequency map in voxel space. Copy the weight image and swap its .dat rather
% than round-tripping through weight_map_object on a mutated pm - that path
% needs to write pm.weights, which is protected. mvpa_dat is accepted for the
% fallback case where no weight_obj exists yet.
if isfield(pm.weights,'weight_obj') && ~isempty(pm.weights.weight_obj)
    fobj = pm.weights.weight_obj;
else
    fobj = weight_map_object(pm, mvpa_dat);
    fobj = fobj.weights.weight_obj;
end
if numel(fobj.dat) == numel(ss.selection_freq)
    fobj.dat = ss.selection_freq(:);
    if isprop(fobj,'p')   || isfield(struct(fobj),'p'),   fobj.p   = []; end
    if isprop(fobj,'sig') || isfield(struct(fobj),'sig'), fobj.sig = []; end
    ss.freq_obj = fobj;
else
    warning('mvpa_reg_cov_stability_from_boot:sizeMismatch', ...
        ['weight image has %d voxel(s) but there are %d selection frequencies; ' ...
         'no freq_obj returned.'], numel(fobj.dat), numel(ss.selection_freq));
end

ev = o.k^2 / ((2*o.threshold - 1) * p_feat);
fprintf('  k = %d, pi = %.2f, p = %d -> E(V) <= %.2f (Meinshausen & Buhlmann 2010)\n', ...
        o.k, o.threshold, p_feat, ev);

end
