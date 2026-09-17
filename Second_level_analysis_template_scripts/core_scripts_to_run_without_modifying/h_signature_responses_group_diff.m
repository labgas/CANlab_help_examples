%% h_signature_responses_group_diff.m
%
%
% *USAGE*
%
% This script runs a two-sample t-test for a set of signature responses,
% for each contrast, plotting group differences in signature response and
% printing between-group test statistics for each. Calls
% plugin_get_group_names_colors, barplot_columns, and ttest2_printout
% under the hood.
%
% NOTE: NPS-subregion group-difference code exists at the bottom of the
% script but is currently commented out/inactive.
%
%
% *OPTIONS*
%
% NOTE: this script has no dedicated section of its own in
% a2_set_default_options.m. mysignature/scalenames/simnames below are set
% directly from that file's PREP_4_APPLY_SIGNATURES_AND_SAVE &
% H_SIGNATURE_RESPONSES_GROUP_DIFF section (keyword_sigs, myscaling_sigs,
% similarity_metric_sigs respectively) - the same options prep_4 uses.
%
% mygroupnamefield = 'contrasts';   hardcoded below (not an a2 option);
%                                   change to 'conditions' directly in
%                                   this script if needed
%
%
% *NOTES*
%
% Group membership comes from DAT.BETWEENPERSON.group, or the
% condition/contrast-specific DAT.BETWEENPERSON.conditions/.contrasts
% fields if set, both defined in prep_1b_prep_behavioral_data.m. Continuous
% grouping variables are binarized via median split.
%
% Unlike other Group 2 scripts, this script does not call
% a_set_up_paths_always_run_first or reload DAT/DATA_OBJ from saved .mat
% files itself; it assumes these are already in the workspace from a
% previous script run earlier in the same MATLAB session (e.g.
% prep_4_apply_signatures_and_save.m).
%
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove
%
% date:   Leuven, December, 2024
%
% -------------------------------------------------------------------------
%
% h_signature_responses_group_diff.m       v1.1
%
% last modified: 2026/08/13


%% USER OPTIONS
% -------------------------------------------------------------------------

% Covariates to adjust the between-group test for. Empty = the unadjusted
% t-test only, which is what this script used to do.
%
% Names must be columns of DAT.BETWEENPERSON.contrasts{i} (or .conditions{i}),
% i.e. whatever prep_1b put in the design. An UNORDERED FACTOR must already be
% dummy-coded there, with one column per non-reference level, and every column
% named here: a three-centre study passes {'center_UGOT','center_UM'}, not
% {'center'}. Squeezing k levels into one numeric column treats them as ordered
% and spends one df where k-1 are needed.
%
% When set, each contrast additionally gets a covariate-adjusted group test
% (signature ~ group + covariates) and a companion barplot of the adjusted
% responses, so the adjusted and unadjusted results sit side by side.
adjust_for_covs = {};


% Now set in a2 script
mysignature =   keyword_sigs;                           % 'NPS' 'NPSpos' 'NPSneg' 'SIIPS' etc.  See load_image_set('npsplus')
scalenames =    {myscaling_sigs};                       % or scaled
simnames =      {similarity_metric_sigs};               % or 'cosine_sim' 'dotproduct'
mygroupnamefield = 'contrasts';                         % 'conditions' or 'contrasts'


%% DEFINE GROUPS IN PREP_1b_PREP_BEHAVIORAL_DATA
% -------------------------------------------------------------------------

% There are two ways to define groups. The other is condition- and
% contrast-specific.  These are entered in DAT.BETWEENPERSON.conditions and
% DAT.BETWEENPERSON.contrasts, in cells.  1 -1 codes work best. 
% These are set up in prep_1b_prep_behavioral_data.m
% If they are missing, using DAT.BETWEENPERSON.group will be used as a
% generic option.
% -------------------------------------------------------------------------

group = DAT.BETWEENPERSON.group;


%% LOOP THROUGH SIGNATURES, TEST GROUP DIFFERENCE, CREATE ONE PLOT PER CONTRAST
% -------------------------------------------------------------------------
% Collect every signature's group test so the set can be FDR-corrected below.
% Signatures are tested one at a time above; without this the report gives a
% dozen or more uncorrected p-values and no way to read them as a family.
sig_fdr = struct('name',{},'contrast',{},'p_unadj',{},'p_adj',{});

for s = 1:length(mysignature)
    
    % Get data
    % -------------------------------------------------------------------------
    if isstruct(DAT.SIG_contrasts.(scalenames{1}).(simnames{1}).(mysignature{s}))    % this is a group of signatures rather than an individual one
        
        siggroup = DAT.SIG_contrasts.(scalenames{1}).(simnames{1}).(mysignature{s});
        
        for sig = 1:size(siggroup.signaturenames,2)
            
            signature = siggroup.signaturenames{1,sig};
            
            contrastdata = table2array(siggroup.(signature));
            
            kc = size(contrastdata, 2);

            % Plot
            % -------------------------------------------------------------------------
            printhdr(sprintf('%s responses: Scale = %s Metric = %s', signature, scalenames{1}, simnames{1}));

            figtitle = sprintf('%s group diffs %s %s', signature, scalenames{1}, simnames{1});
            fh_sig = create_figure(figtitle, 1, kc);

            for i = 1:kc

                % Load group variable
                [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i);

                if isempty(group), continue, end % skip this condition/contrast - no groups

                if size(group, 2) > 1            % we have multiple variables
                    disp('Warning: Group has > 1 column. Using first column only.')
                    group = group(:, 1);
                end

                if length(unique(group)) > 2  % this is a continuous variable
                    disp('Binarizing continuous grouping variable via median split.')
                    group = mediansplit(group);  
                end

                % In case we forgot to assign colors or groupnames in prep 1b script
                if length(groupcolors) < 2, groupcolors = seaborn_colors(2); end
                if length(groupnames) < 2, groupnames = {'High' 'Low'}; end

                % Must code data with pos or neg values
                y = {contrastdata(group > 0, i) contrastdata(group < 0, i)};

                subplot(1, kc, i)

                printstr(' ');
                printstr(sprintf('Group differences: %s, %s', signature, DAT.contrastnames{i}));
                printstr(dashes)

                barplot_columns(y, 'nofig', 'colors', groupcolors, 'names', groupnames);

                title(DAT.contrastnames{i})
                xlabel('Group');
                ylabel(sprintf('%s Response', signature));

                printstr('Between-groups test:');

                [H,p,ci,stats] = ttest2_printout(y{1}, y{2});
            p_unadj_this = p;
                p_unadj_this = p;

                printstr(dashes)

            % covariate-adjusted test + companion plot, when requested
            if ~isempty(adjust_for_covs)
                [y_adj_all{i}, adj_stats_this] = h_adjusted_group_test(y, group, DAT, i, adjust_for_covs);
                sig_fdr(end+1) = struct('name',signature,'contrast',DAT.contrastnames{i}, ...
                    'p_unadj',p_unadj_this,'p_adj',adj_stats_this.p); %#ok<SAGROW>
            end

            end % panels

            % Size the panel figure before it is captured. The template never did,
            % so these barplots came out at MATLAB's 480x420 default while every other
            % figure in the reports is sized - most visibly next to a covariate-adjusted
            % companion plot, which is sized and so looked like a different report.
            % Sized by HANDLE, not gcf: barplot_columns can leave another figure
            % current, and create_figure reuses a figure carrying the same tag rather
            % than opening a new one, so neither gcf nor a new-figure test is reliable.
            if exist('fh_sig','var') && all(isgraphics(fh_sig))
                            plugin_set_figure_size('fig', fh_sig);
            end

            if ~isempty(adjust_for_covs)
                h_plot_adjusted(y_adj_all, DAT, kc, groupcolors, groupnames, signature, adjust_for_covs);
            end

            drawnow, snapnow
            
            clear signature contrastdata
            
        end
        
    else
    
        contrastdata = table2array(DAT.SIG_contrasts.(scalenames{1}).(simnames{1}).(mysignature{s}));
    
        kc = size(contrastdata, 2);

        % Plot
        % -------------------------------------------------------------------------
        printhdr(sprintf('%s responses: Scale = %s Metric = %s', mysignature{s}, scalenames{1}, simnames{1}));

        figtitle = sprintf('%s group diffs %s %s', mysignature{s}, scalenames{1}, simnames{1});
        fh_sig = create_figure(figtitle, 1, kc);

        for i = 1:kc

            % Load group variable
            [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i);

            if isempty(group), continue, end % skip this condition/contrast - no groups

            if size(group, 2) > 1            % we have multiple variables
                disp('Warning: Group has > 1 column. Using first column only.')
                group = group(:, 1);
            end

            if length(unique(group)) > 2  % this is a continuous variable
                disp('Binarizing continuous grouping variable via median split.')
                group = mediansplit(group);  
            end

            % In case we forgot to assign colors or groupnames in prep 1b script
            if length(groupcolors) < 2, groupcolors = seaborn_colors(2); end
            if length(groupnames) < 2, groupnames = {'High' 'Low'}; end

            % Must code data with pos or neg values
            y = {contrastdata(group > 0, i) contrastdata(group < 0, i)};

            subplot(1, kc, i)

            printstr(' ');
            printstr(sprintf('Group differences: %s, %s', mysignature{s}, DAT.contrastnames{i}));
            printstr(dashes)

            barplot_columns(y, 'nofig', 'colors', groupcolors, 'names', groupnames);

            title(DAT.contrastnames{i})
            xlabel('Group');
            ylabel(sprintf('%s Response', mysignature{s}));

            printstr('Between-groups test:');

            [H,p,ci,stats] = ttest2_printout(y{1}, y{2});

            printstr(dashes)

            if ~isempty(adjust_for_covs)
                [y_adj_all{i}, adj_stats_this] = h_adjusted_group_test(y, group, DAT, i, adjust_for_covs);
                sig_fdr(end+1) = struct('name',mysignature{s},'contrast',DAT.contrastnames{i}, ...
                    'p_unadj',p_unadj_this,'p_adj',adj_stats_this.p); %#ok<SAGROW>
            end

        end % panels

        % Size the panel figure before it is captured. The template never did,
        % so these barplots came out at MATLAB's 480x420 default while every other
        % figure in the reports is sized - most visibly next to a covariate-adjusted
        % companion plot, which is sized and so looked like a different report.
        % Sized by HANDLE, not gcf: barplot_columns can leave another figure
        % current, and create_figure reuses a figure carrying the same tag rather
        % than opening a new one, so neither gcf nor a new-figure test is reliable.
        if exist('fh_sig','var') && all(isgraphics(fh_sig))
                    plugin_set_figure_size('fig', fh_sig);
        end

        if ~isempty(adjust_for_covs)
            h_plot_adjusted(y_adj_all, DAT, kc, groupcolors, groupnames, mysignature{s}, adjust_for_covs);
        end

        drawnow, snapnow
    
    end % if loop group or individual signature
    
end % signature


%% FDR CORRECTION ACROSS SIGNATURES
% -------------------------------------------------------------------------
% Each signature above is tested on its own. A panel of a dozen or more
% signatures is a family, and reading those p-values uncorrected overstates
% the evidence - which is exactly how an apparently strong result can turn out
% to be one of fourteen chances.
%
% Correction is applied WITHIN contrast, across signatures, since that is the
% family actually being read. Both q_BH and q_Storey are reported, using the
% same LaBGAScore_Storey_FDR the roi GLM in prep_3a and the decoding scripts
% use, so the three cannot drift apart.
%
% Storey needs enough tests to estimate pi0. With a handful of signatures it
% usually cannot, declares the estimate unreliable and falls back to pi0 = 1,
% at which point q_Storey is identical to q_BH. That is the intended
% behaviour, not a failure - read q_BH whenever the note below says so.
%
% The ADJUSTED p is corrected when adjust_for_covs is set, because that is the
% test being interpreted; the unadjusted column is corrected too, for
% comparison.

if ~isempty(sig_fdr)

    fprintf('\n\n');
    printhdr('FDR CORRECTION ACROSS SIGNATURES');
    fprintf('\n\n');

    contrasts_done = unique({sig_fdr.contrast}, 'stable');

    for cc = 1:numel(contrasts_done)

        sel = strcmp({sig_fdr.contrast}, contrasts_done{cc});
        nms = {sig_fdr(sel).name};
        pu  = [sig_fdr(sel).p_unadj];
        pa  = [sig_fdr(sel).p_adj];

        fprintf('\n%s  (%d signatures)\n', contrasts_done{cc}, numel(nms));

        use_adj = ~isempty(adjust_for_covs) && any(~isnan(pa));
        if use_adj
            [q_st_a, pi0_a, info_a] = LaBGAScore_Storey_FDR(pa(~isnan(pa)));
            q_bh_a = info_a.q_BH(:)';
        end
        [q_st_u, pi0_u, info_u] = LaBGAScore_Storey_FDR(pu(~isnan(pu)));
        q_bh_u = info_u.q_BH(:)';

        if use_adj
            fprintf('  %-16s %10s %9s %10s | %10s %9s %10s\n', ...
                'signature','p_unadj','q_BH','q_Storey','p_adj','q_BH','q_Storey');
        else
            fprintf('  %-16s %10s %9s %10s\n','signature','p_unadj','q_BH','q_Storey');
        end

        ku = 0; ka = 0;
        for k = 1:numel(nms)
            if ~isnan(pu(k)), ku = ku + 1; qbu = q_bh_u(ku); qsu = q_st_u(ku); else, qbu = NaN; qsu = NaN; end
            if use_adj && ~isnan(pa(k)), ka = ka + 1; qba = q_bh_a(ka); qsa = q_st_a(ka); else, qba = NaN; qsa = NaN; end
            if use_adj
                fprintf('  %-16s %10.4f %9.4f %10.4f | %10.4f %9.4f %10.4f\n', nms{k}, pu(k), qbu, qsu, pa(k), qba, qsa);
            else
                fprintf('  %-16s %10.4f %9.4f %10.4f\n', nms{k}, pu(k), qbu, qsu);
            end
        end

        fprintf('\n  unadjusted: %d/%d at q_BH < .05, %d at q_Storey < .05 (pi0 = %.3f)\n', ...
            sum(q_bh_u < .05), numel(q_bh_u), sum(q_st_u < .05), pi0_u);
        if ~info_u.reliable
            fprintf('    pi0 not identifiable for this set, so q_Storey equals q_BH - read q_BH\n');
        end
        if use_adj
            fprintf('  adjusted  : %d/%d at q_BH < .05, %d at q_Storey < .05 (pi0 = %.3f)\n', ...
                sum(q_bh_a < .05), numel(q_bh_a), sum(q_st_a < .05), pi0_a);
            if ~info_a.reliable
                fprintf('    pi0 not identifiable for this set, so q_Storey equals q_BH - read q_BH\n');
            end
        end

    end

end



% %% NPS SUBREGIONS, TEST GROUP DIFFERENCE, CREATE ONE PLOT PER CONTRAST
% % -------------------------------------------------------------------------
% 
% % POSITIVE
% % --------
% 
% % which variables to use
% mysubrfield = 'npspos_by_region_contrasts'; % 'npspos_by_region_cosinesim';     %'npspos_by_regionsc';
% mysubrfieldneg = 'npsneg_by_region_contrasts'; % 'npsneg_by_region_cosinesim';  % 'npsneg_by_regionsc';
% 
% posnames = DAT.NPSsubregions.posnames;
% negnames = DAT.NPSsubregions.negnames;
% 
% clear means p T
% 
% for i = 1:kc  % for each contrast
%     
%     [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i);
%     if isempty(group), continue, end % skip this condition/contrast - no groups
%     
%     if size(group, 2) > 1            % we have multiple variables
%         disp('Warning: Group has > 1 column. Using first column only.')
%         group = group(:, 1);
%     end
%     
%     if length(unique(group)) > 2  % this is a continuous variable
%         disp('Binarizing continuous grouping variable via median split.')
%         group = mediansplit(group);
%     end
%     
%     % In case we forgot to assign colors or groupnames in prep 1b script
%     if length(groupcolors) < 2, groupcolors = seaborn_colors(2); end
%     if length(groupnames) < 2, groupnames = {'High' 'Low'}; end
%         
%     mydat = DAT.NPSsubregions.(mysubrfield){i};
%     k = size(mydat, 2);
%     
%     create_figure(sprintf('NPS subregions by group %s', DAT.contrastnames{i}), 1, k);
%     pos = get(gcf, 'Position');
%     pos(4) = pos(4) .* 2.5;
%     set(gcf, 'Position', pos);
%     
%     clear means p T
%     
%     for j = 1:k  % for each subregion
%         
%         subplot(1, k, j);
%         
%         y = {mydat(group == 1, j) mydat(group == -1, j)};
%         
%         printhdr(posnames{j});
%         
%         barplot_columns(y, 'nofig', 'colors', groupcolors, 'noviolin', 'noind', 'names', groupnames );
%         
%         title(posnames{j})
%         xlabel('Group');
%         if j == 1, ylabel('NPS Response'); end
%         
%         printstr('Between-groups test:');
%         [H,p(j, 1),ci,stats] = ttest2_printout(y{1}, y{2});
%         
%         means(j, :) = stats.means;
%         T(j, 1) = stats.tstat;
%         
%     end
%     
%     drawnow, snapnow
%     
%     % Print between-subject Table
%     printhdr('Between-group tests');
%     Region = posnames';
%     regionmeans = table(Region, means, T, p);
%     
%     disp(regionmeans);
%     
% end % panels
% 
% % NEGATIVE
% % --------
% 
% for i = 1:kc
%     
%     [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i);
%     
%     if size(group, 2) > 1            % we have multiple variables
%         disp('Warning: Group has > 1 column. Using first column only.')
%         group = group(:, 1);
%     end
%     
%     if length(unique(group)) > 2  % this is a continuous variable
%         disp('Binarizing continuous grouping variable via median split.')
%         group = mediansplit(group);
%     end
%     
%     % In case we forgot to assign colors or groupnames in prep 1b script
%     if length(groupcolors) < 2, groupcolors = seaborn_colors(2); end
%     if length(groupnames) < 2, groupnames = {'High' 'Low'}; end
%     
%         
%     if isempty(group), continue, end % skip this condition/contrast - no groups
%     
%     mydat = DAT.NPSsubregions.(mysubrfieldneg){i};
%     k = size(mydat, 2);
%     
%     create_figure(sprintf('NPS neg subregions by group %s', DAT.contrastnames{i}), 1, k);
%     pos = get(gcf, 'Position');
%     pos(4) = pos(4) .* 2.5;
%     set(gcf, 'Position', pos);
%     
%     clear means p T
%     
%     for j = 1:k
%         
%         subplot(1, k, j);
%         
%         y = {mydat(group == 1, j) mydat(group == -1, j)};
%         
%         printhdr(negnames{j});
%         
%         barplot_columns(y, 'nofig', 'colors', groupcolors, 'noviolin', 'noind', 'names', groupnames );
%         
%         title(negnames{j})
%         xlabel('Group');
%         if j == 1, ylabel('NPS Response'); end
%         
%         printstr('Between-groups test:');
%         [H,p(j, 1),ci,stats] = ttest2_printout(y{1}, y{2});
%         
%         means(j, :) = stats.means;
%         T(j, 1) = stats.tstat;
%         
%     end
%     
%     drawnow, snapnow
%     
%     % Print between-subject Table
%     printhdr('Between-group tests');
%     Region = negnames';
%     regionmeans = table(Region, means, T, p);
%     
%     disp(regionmeans);
%     
% end


% =========================================================================

function [y_adj, adj_stats] = h_adjusted_group_test(y, group_i, DAT, i, adjust_for_covs)
% Print the group difference ADJUSTED for the named covariate(s), and return
% the covariate-adjusted values split by group so they can be plotted.
%
% Returns [] when the adjustment cannot be made, so the caller can skip the
% panel without special-casing.
%
% The unadjusted ttest2 printed by the caller is left exactly as it was, so
% the two are directly comparable and nothing existing changes meaning.

y_adj = [];
% Second output carries the adjusted test so the caller can correct across
% signatures. Initialised before every early return, so a skipped adjustment
% yields NaN rather than an undefined variable.
adj_stats = struct('p',NaN,'t',NaN,'df',NaN,'diff',NaN,'covs',{{}});
if isempty(adjust_for_covs), return, end

if ~isfield(DAT,'BETWEENPERSON') || ~isfield(DAT.BETWEENPERSON,'contrasts') ...
        || numel(DAT.BETWEENPERSON.contrasts) < i || isempty(DAT.BETWEENPERSON.contrasts{i})
    fprintf('\nno covariate table for this contrast; adjusted test skipped\n');
    return
end

T = DAT.BETWEENPERSON.contrasts{i};
have = adjust_for_covs(ismember(adjust_for_covs, T.Properties.VariableNames));
if isempty(have)
    fprintf('\ncovariate(s) %s not in the design; adjusted test skipped\n', ...
        strjoin(adjust_for_covs, ', '));
    return
end

% Rebuild the full-length response in the original subject order: barplot
% wants it split by group, the model wants it whole.
group_i = group_i(:);
yy = nan(numel(group_i), 1);
yy(group_i > 0) = y{1};
yy(group_i < 0) = y{2};

ok = ~isnan(yy);
X  = table();
X.signature = yy(ok);
X.group     = group_i(ok);
for k = 1:numel(have)
    X.(have{k}) = double(T.(have{k})(ok));
end

mdl = fitlm(X, ['signature ~ group + ' strjoin(have, ' + ')]);

r  = mdl.Coefficients;
b  = r.Estimate('group');  se = r.SE('group');
tv = r.tStat('group');     pv = r.pValue('group');
ci = b + [-1 1] * se * tinv(0.975, mdl.DFE);

fprintf('%s\n', sprintf('Between-groups test ADJUSTED for %s:', strjoin(have, ', ')));
% group is coded -1/1, so the difference BETWEEN groups is 2*beta
fprintf('adjusted group difference = %3.4f, 95%% CI [%3.4f %3.4f], t(%d) = %3.2f, p = %3.6f\n', ...
    2*b, 2*ci(1), 2*ci(2), mdl.DFE, tv, pv);
fprintf('model R-squared = %3.4f (adjusted %3.4f), n = %d\n', ...
    mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted, mdl.NumObservations);

adj_stats = struct('p',pv,'t',tv,'df',mdl.DFE,'diff',2*b,'covs',{have});
for k = 1:numel(have)
    fprintf('   covariate %-10s beta = %3.4f, t(%d) = %3.2f, p = %3.6f\n', ...
        have{k}, r.Estimate(have{k}), mdl.DFE, r.tStat(have{k}), r.pValue(have{k}));
end

% Adjusted values for plotting: residualise on the COVARIATES ONLY and add
% the grand mean back, so the bars keep the original units and differ from
% the raw plot only by what the covariate explains. The group term is
% deliberately not removed - that is the effect being shown.
Xc = X;  Xc.group = [];
mdl_cov = fitlm(Xc, ['signature ~ ' strjoin(have, ' + ')]);
adj = mdl_cov.Residuals.Raw + mean(X.signature);

g = X.group;
y_adj = {adj(g > 0) adj(g < 0)};

end


function h_plot_adjusted(y_adj_all, DAT, kc, groupcolors, groupnames, siglabel, adjust_for_covs)
% One figure of adjusted panels, mirroring the raw figure the caller drew.

if isempty(adjust_for_covs) || all(cellfun(@isempty, y_adj_all)), return, end

figtitle = sprintf('%s group diffs adjusted for %s', siglabel, strjoin(adjust_for_covs, ', '));
fh_adj = create_figure(figtitle, 1, kc);

for i = 1:kc
    if isempty(y_adj_all{i}), continue, end
    subplot(1, kc, i)
    barplot_columns(y_adj_all{i}, 'nofig', 'colors', groupcolors, 'names', groupnames);
    title(DAT.contrastnames{i})
    xlabel('Group');
    ylabel(sprintf('%s (adj. %s)', siglabel, strjoin(adjust_for_covs, ', ')));
end

set(fh_adj, 'Tag', figtitle);
% Size by HANDLE, not gcf: barplot_columns can leave another figure current,
% so plugin_set_figure_size() with no 'fig' may size the wrong one and leave
% this figure at the 480x420 default.
plugin_set_figure_size('fig', fh_adj);
drawnow, snapnow;

end
