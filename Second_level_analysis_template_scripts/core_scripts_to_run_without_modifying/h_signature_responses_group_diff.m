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


% NOTE ON ORDER: the path/data block MUST come before the user options below.
% It calls a2_set_default_options, which (re)assigns every option a2 defines -
% so anything set above it is silently discarded, and options that READ a2
% values (mysignature = keyword_sigs, and the scaling/metric names) would be
% referencing variables that do not exist yet. The other core scripts already
% order it this way; this one did not.

%% LOAD PATHS AND DATA IF NEEDED
% -------------------------------------------------------------------------
% These scripts were written to run straight after their s10 counterpart, in
% the same MATLAB session, and so relied on DAT already being in the
% workspace. Run on their own they failed at the first DAT reference with
% "Unable to resolve the name DAT.BETWEENPERSON.group". The cfs h1 script has
% always carried this guard; adding it here makes the two studies behave the
% same and lets these be re-run without redoing prep_4.

% Remember where the study's own setup put the results, so the call below can be
% checked against it (see the guard immediately after).
resultsdir_before_setup = '';
if exist('resultsdir','var'), resultsdir_before_setup = resultsdir; end


a_set_up_paths_always_run_first;   % NOTE: replace with your study-specific s0 script

% GUARD: did the path setup just move the output directory?
%
% The call above is meant to be replaced, in a study's copy, by that study's own
% s0 (e.g. mystudy_secondlevel_m2a_s0_a_set_up_paths_always_run_first). Left as
% the generic call, it RE-DERIVES resultsdir - typically from the FIRST-LEVEL
% model name - and silently overwrites whatever the study's setup had already
% set. Every result then lands in a different model's directory while the
% published report still goes to the right one, so the split is easy to miss.
if ~isempty(resultsdir_before_setup) && ~strcmp(resultsdir_before_setup, resultsdir)
    error(['\nPATH SETUP MOVED THE RESULTS DIRECTORY.\n\n' ...
           '  before: %s\n  after : %s\n\n' ...
           'The generic a_set_up_paths_always_run_first re-derived resultsdir and\n' ...
           'discarded the one your study setup had set. In your copy of this script,\n' ...
           'replace that call with your study''s own s0 path script.\n'], ...
           resultsdir_before_setup, resultsdir);
end


if ~exist('DAT','var') || ~isfield(DAT,'SIG_contrasts')
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
end

% The reload above is not proof the data arrived: if prep_2 or prep_3 was re-run
% after prep_4, the stored DAT was rebuilt without SIG_* and the file has none
% either. Without this check the script carries on and fails later somewhere
% less obvious - or appears to run on nothing. proj_moodbugs model_1b sat in
% exactly this state, with a published prep_4 report and no SIG fields on disk.
if ~isfield(DAT,'SIG_contrasts') || isempty(DAT.SIG_contrasts)
    error(['\nNO SIGNATURE RESULTS IN DAT.\n\n' ...
           'DAT.SIG_contrasts is missing from %s.\n\n' ...
           'prep_4_apply_signatures_and_save appends it, but prep_2 and prep_3 ' ...
           'rebuild DAT\nand overwrite the stored copy, so re-running either of ' ...
           'them after prep_4\ndiscards the signature results.\n\n' ...
           'Re-run prep_4_apply_signatures_and_save, then this script.\n'], ...
           fullfile(resultsdir,'image_names_and_setup.mat'));
end



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
adjust_for_covs = {};   % names of DAT.BETWEENPERSON.contrasts{i} columns to adjust for; {} for none.

% TEST TYPE: 'group' (default) or 'continuous'.
%
% 'group'      two-sample t-test of the signature response between the two
%              levels of DAT.BETWEENPERSON.group. Unchanged behaviour - every
%              existing study copy keeps working without edits.
% 'continuous' linear regression of the signature response on a CONTINUOUS
%              predictor named by sig_covariate_name. Reports beta, t and
%              partial r instead of a mean difference and Cohen's d.
%
% Everything downstream - the correction family across signatures, the summary
% table, the subregion pass - is shared. Only the test and the plot differ.
% prep_4 is untouched: it computes the signature responses either way.
sig_test_type = 'group';

% For sig_test_type = 'continuous': the column of
% DAT.BETWEENPERSON.contrasts{i} (or .conditions{i}) holding the predictor.
sig_covariate_name = '';

% Validate the test type before any analysis runs, so a typo or a missing
% predictor fails immediately rather than part way through the signature loop.
sig_test_type = lower(char(sig_test_type));
if ~ismember(sig_test_type, {'group','continuous'})
    error('sig_test_type must be ''group'' or ''continuous'', not ''%s''.', sig_test_type);
end
if isequal(sig_test_type,'continuous')
    if isempty(sig_covariate_name)
        error(['sig_test_type = ''continuous'' requires sig_covariate_name to name ' ...
               'a column of DAT.BETWEENPERSON.%s{:}.'], mygroupnamefield);
    end
    fprintf('\nTEST TYPE: continuous - regression on %s\n', sig_covariate_name);
else
    fprintf('\nTEST TYPE: group - two-sample t-test\n');
end

                        % ComBat already removed centre from the condition images, and the
                        % corresponding GLM arm (nocov) does not covary centre either. Adjusting
                        % here would correct twice for the same thing and would make the
                        % signature results describe a different model from the GLM they sit
                        % beside. Set to {'center_UGOT','center_UM'} only against a cov_center arm.


% Now set in a2 script
mysignature =   keyword_sigs;                           % 'NPS' 'NPSpos' 'NPSneg' 'SIIPS' etc.  See load_image_set('npsplus')
scalenames =    {myscaling_sigs};                       % or scaled
simnames =      {similarity_metric_sigs};               % or 'cosine_sim' 'dotproduct'
mygroupnamefield = 'contrasts';                         % 'conditions' or 'contrasts'

% Signature selection, matched to proj_cfs so the two studies report the same
% set. Empty = report every signature in the group.
% FIRST PASS: NPSpos and NPSneg are LEFT OUT. They are the positive- and
% negative-weight halves of NPS, not independent signatures, so including all
% three tests the same pattern three times and inflates the family the FDR
% correction is applied over. Test NPS itself first; only if NPS is significant
% (unadjusted) is it worth decomposing it into NPSpos/NPSneg and the NPS
% subregions, for which prep_4 already stores DAT.NPSsubregions and there is
% commented-out group-difference code at the bottom of this script.
subsets_i_want = {'NPS','SIIPS','PINES','GSR','Heart','FM_pain'};



% WHICH MULTIPLE-COMPARISON CORRECTIONS TO REPORT
%
% 'BH'          Benjamini-Hochberg linear step-up                    FDR
% 'Storey'      Storey q-values via LaBGAScore_Storey_FDR ('sas')     FDR
% 'adaptiveFDR' Benjamini & Hochberg (2000) adaptive step-up, m0 by
%               lowest slope (SAS ADAPTIVEFDR default)               FDR
% 'BKY'         Benjamini, Krieger & Yekutieli (2006) two-stage      FDR
% 'holmSidak'   Holm step-down with Sidak multiplier                 FWER
%
% BH and Storey only by default: those are the two the pipeline has always
% reported, and a table with every method in it is harder to read, not easier.
% Add the others when the question is specifically whether a result depends on
% the choice of correction - at small m it often does.
%
% Note holmSidak is reported as p_holmSidak, not q_: it is an FWER-adjusted
% p-value, a different quantity from the FDR q-values, not a stricter version
% of one.
if ~exist('corrections_i_want','var') || isempty(corrections_i_want)
    corrections_i_want = {'BH','Storey'};
end

% Output tag. This script writes its summary table to a name built from the
% metric, the scaling and the contrast - none of which change when you run it a
% second time over a DIFFERENT subsets_i_want. Two runs on the same contrast
% therefore collide, and the second silently overwrites the first: running an
% NPS decomposition after the main panel replaced the main panel's table with a
% two-row one. Set sig_results_tag in the second script to keep them apart.
if ~exist('sig_results_tag','var') || isempty(sig_results_tag), sig_results_tag = ''; end

if ~exist('corrections_i_want_subregions','var') || isempty(corrections_i_want_subregions)
    corrections_i_want_subregions = {'BH','Storey','holmSidak'};
end

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
% The field list must match the records appended below EXACTLY, including the
% continuous-only fields, or the first append fails with "Subscripted
% assignment between dissimilar structures". They are populated for both test
% types - NaN under 'group' - so one record shape serves both.
sig_fdr = struct('name',{},'contrast',{},'p_unadj',{},'p_adj',{}, ...
                 'y1',{},'y2',{},'t',{},'df',{},'d',{},'diff_adj',{},'t_adj',{}, ...
                 'y1_adj',{},'y2_adj',{},'beta',{},'r_partial',{},'n_cont',{});

for s = 1:length(mysignature)
    
    % Get data
    % -------------------------------------------------------------------------
    if isstruct(DAT.SIG_contrasts.(scalenames{1}).(simnames{1}).(mysignature{s}))    % this is a group of signatures rather than an individual one
        
        siggroup = DAT.SIG_contrasts.(scalenames{1}).(simnames{1}).(mysignature{s});
        
        for sig = 1:size(siggroup.signaturenames,2)
            
            signature = siggroup.signaturenames{1,sig};

            % Report only the signatures named in subsets_i_want, so this study
            % and proj_cfs cover the same set. Empty = report all of them.
            if ~isempty(subsets_i_want) && ~ismember(signature, subsets_i_want)
                continue
            end
            
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

                if isequal(sig_test_type,'continuous')
                
                    % CONTINUOUS BRANCH. y is left empty on purpose: the group summary
                    % columns (n_g1, mean_g1, ...) are undefined for a regression, and the
                    % table block below branches on sig_test_type rather than trying to
                    % synthesise them.
                    xv = h_sig_get_covariate(DAT, mygroupnamefield, i, sig_covariate_name);
                    [cm, cn] = h_sig_get_covmat(DAT, mygroupnamefield, i, adjust_for_covs);
                    cont_this = h_sig_continuous_test(contrastdata(:,i), xv, cm, cn, ...
                                    signature, DAT.contrastnames{i}, kc, i);
                    y = {[] []};
                    p_unadj_this = cont_this.p;
                    stats = struct('tstat', cont_this.t, 'df', cont_this.df);
                
                else
                
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
                    cont_this = [];
                
                end

                printstr(dashes)

            % Covariate-adjusted test, when requested. The sig_fdr record is made
            % UNCONDITIONALLY: it used to sit inside this guard, so turning the
            % adjustment off emptied sig_fdr and silently skipped both the FDR
            % correction and the summary table - the unadjusted analysis lost its
            % multiple-comparison correction precisely when it was the only one being
            % reported. The adjusted fields are NaN/[] when no adjustment was made.
            if ~isempty(adjust_for_covs)
                [y_adj_all{i}, adj_stats_this] = h_adjusted_group_test(y, group, DAT, i, adjust_for_covs);
                if iscell(y_adj_all{i}) && numel(y_adj_all{i}) == 2
                    y1a = y_adj_all{i}{1}; y2a = y_adj_all{i}{2};
                else
                    y1a = []; y2a = [];
                end
            else
                y_adj_all{i} = [];
                adj_stats_this = struct('p',NaN,'t',NaN,'df',NaN,'diff',NaN,'covs',{{}});
                y1a = []; y2a = [];
            end
            sig_fdr(end+1) = struct('name',signature,'contrast',DAT.contrastnames{i}, ...
                    'p_unadj',p_unadj_this,'p_adj',adj_stats_this.p, ...
                    'y1',{y{1}},'y2',{y{2}},'t',stats.tstat,'df',stats.df, ...
                    'd',stats.tstat*sqrt(1/numel(y{1})+1/numel(y{2})), ...
                    'diff_adj',adj_stats_this.diff,'t_adj',adj_stats_this.t, ...
                    'y1_adj',{y1a},'y2_adj',{y2a}, ...
                    'beta',h_sig_field(cont_this,'beta'), ...
                    'r_partial',h_sig_field(cont_this,'r_partial'), ...
                    'n_cont',h_sig_field(cont_this,'n')); %#ok<SAGROW>

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

            if isequal(sig_test_type,'continuous')
            
                % CONTINUOUS BRANCH. y is left empty on purpose: the group summary
                % columns (n_g1, mean_g1, ...) are undefined for a regression, and the
                % table block below branches on sig_test_type rather than trying to
                % synthesise them.
                xv = h_sig_get_covariate(DAT, mygroupnamefield, i, sig_covariate_name);
                [cm, cn] = h_sig_get_covmat(DAT, mygroupnamefield, i, adjust_for_covs);
                cont_this = h_sig_continuous_test(contrastdata(:,i), xv, cm, cn, ...
                                mysignature{s}, DAT.contrastnames{i}, kc, i);
                y = {[] []};
                p_unadj_this = cont_this.p;
                stats = struct('tstat', cont_this.t, 'df', cont_this.df);
            
            else
            
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
                p_unadj_this = p;
                cont_this = [];
            
            end

            printstr(dashes)

            % The sig_fdr record is made UNCONDITIONALLY: it used to sit inside this
            % guard, so turning the adjustment off emptied sig_fdr and silently
            % skipped both the FDR correction and the summary table. Adjusted fields
            % are NaN/[] when no adjustment was made.
            if ~isempty(adjust_for_covs)
                [y_adj_all{i}, adj_stats_this] = h_adjusted_group_test(y, group, DAT, i, adjust_for_covs);
                if iscell(y_adj_all{i}) && numel(y_adj_all{i}) == 2
                    y1a = y_adj_all{i}{1}; y2a = y_adj_all{i}{2};
                else
                    y1a = []; y2a = [];
                end
            else
                y_adj_all{i} = [];
                adj_stats_this = struct('p',NaN,'t',NaN,'df',NaN,'diff',NaN,'covs',{{}});
                y1a = []; y2a = [];
            end
            sig_fdr(end+1) = struct('name',mysignature{s},'contrast',DAT.contrastnames{i}, ...
                    'p_unadj',p_unadj_this,'p_adj',adj_stats_this.p, ...
                    'y1',{y{1}},'y2',{y{2}},'t',stats.tstat,'df',stats.df, ...
                    'd',stats.tstat*sqrt(1/numel(y{1})+1/numel(y{2})), ...
                    'diff_adj',adj_stats_this.diff,'t_adj',adj_stats_this.t, ...
                    'y1_adj',{y1a},'y2_adj',{y2a}, ...
                    'beta',h_sig_field(cont_this,'beta'), ...
                    'r_partial',h_sig_field(cont_this,'r_partial'), ...
                    'n_cont',h_sig_field(cont_this,'n')); %#ok<SAGROW>

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



%% SUMMARY TABLE AND VIOLIN PLOT ACROSS SIGNATURES
% -------------------------------------------------------------------------
% One row per signature: the unadjusted and covariate-adjusted group
% differences with their raw p-values, and both FDR corrections.
%
% q_Storey comes from LaBGAScore_Storey_FDR, the same function the ROI GLM and
% the decoding scripts use. With a handful of signatures Storey usually cannot
% identify pi0, falls back to pi0 = 1, and q_Storey then EQUALS q_BH. That is
% the documented behaviour, not a bug - the printout says so when it happens,
% and q_BH is the column to read in that case.
%
% The figure plots the UNADJUSTED cosine similarities, since those are the
% observed data; the adjusted statistics live in the table. Cosine is bounded
% and unitless, so all signatures can share one axis - which is exactly what
% dot product does not allow.

SIGSUM = table();

if ~isempty(sig_fdr)

    fprintf('\n\n');
    printhdr('SUMMARY TABLE: GROUP DIFFERENCES IN COSINE SIMILARITY');
    fprintf('\n\n');

    contrasts_done = unique({sig_fdr.contrast}, 'stable');

    for cc = 1:numel(contrasts_done)

        sel = find(strcmp({sig_fdr.contrast}, contrasts_done{cc}));
        nms = {sig_fdr(sel).name}';
        pu  = [sig_fdr(sel).p_unadj]';
        pa  = [sig_fdr(sel).p_adj]';

        % same correction as the section above, recomputed here so the table is
        % self-contained and cannot silently disagree with the printout
        % Corrections requested in corrections_i_want. Each is computed by
        % LaBGAScore_Storey_FDR under a different 'method', so the table cannot
        % drift from the function: there is one implementation, selected here.
        corr_defs = { 'BH',          'bh',             'q_BH'
                      'Storey',      'sas',            'q_Storey'
                      'adaptiveFDR', 'adaptivefdr',    'q_adaptiveFDR'
                      'BKY',         'bky',            'q_BKY'
                      'holmSidak',   'stepdown_sidak', 'p_holmSidak' };
        unknown_corr = setdiff(corrections_i_want, corr_defs(:,1));
        if ~isempty(unknown_corr)
            error('unknown entry in corrections_i_want: %s. Valid: %s', ...
                strjoin(unknown_corr, ', '), strjoin(corr_defs(:,1)', ', '));
        end
        wh_corr = find(ismember(corr_defs(:,1), corrections_i_want));

        ok_u   = ~isnan(pu);
        corr_u = cell(numel(wh_corr),1);
        for z = 1:numel(wh_corr)
            v = nan(size(pu));
            v(ok_u) = LaBGAScore_Storey_FDR(pu(ok_u), 'method', corr_defs{wh_corr(z),2}, 'verbose', false);
            corr_u{z} = v(:);
        end

        % pi0 and its verdict are Storey-specific; only meaningful if asked for
        if ismember('Storey', corrections_i_want)
            [~, pi0_u, info_u] = LaBGAScore_Storey_FDR(pu(ok_u));
        else
            pi0_u = NaN; info_u = struct('reliable', true);
        end

        % Additional corrections, so the table shows what the choice of method
        % is actually worth rather than resting on one of them. At the panel
        % sizes here Storey's pi0 is estimated from a handful of tests and is
        % biased low, which makes q_Storey the most permissive column by some
        % margin; the two adaptive FDR procedures are biased conservative, and
        % Holm-Sidak controls FWER rather than FDR so it is stricter again and
        % is NOT comparable with the q columns.
        qau = nan(size(pu)); qku = nan(size(pu)); qhu = nan(size(pu));
        qau(~isnan(pu)) = LaBGAScore_Storey_FDR(pu(~isnan(pu)), 'method', 'adaptivefdr',    'verbose', false);
        qku(~isnan(pu)) = LaBGAScore_Storey_FDR(pu(~isnan(pu)), 'method', 'bky',            'verbose', false);
        qhu(~isnan(pu)) = LaBGAScore_Storey_FDR(pu(~isnan(pu)), 'method', 'stepdown_sidak', 'verbose', false);

        use_adj = ~isempty(adjust_for_covs) && any(~isnan(pa));
        corr_a  = cell(numel(wh_corr),1);
        pi0_a   = NaN; info_a = struct('reliable', true);
        if use_adj
            ok_a = ~isnan(pa);
            for z = 1:numel(wh_corr)
                v = nan(size(pa));
                v(ok_a) = LaBGAScore_Storey_FDR(pa(ok_a), 'method', corr_defs{wh_corr(z),2}, 'verbose', false);
                corr_a{z} = v(:);
            end
            if ismember('Storey', corrections_i_want)
                [~, pi0_a, info_a] = LaBGAScore_Storey_FDR(pa(ok_a));
            end
        end

        n = numel(sel);
        m1 = nan(n,1); m2 = nan(n,1); sd1 = nan(n,1); sd2 = nan(n,1);
        n1 = nan(n,1); n2 = nan(n,1); dmean = nan(n,1);
        for k = 1:n
            y1 = sig_fdr(sel(k)).y1(:); y2 = sig_fdr(sel(k)).y2(:);
            if isempty(y1) && isempty(y2), continue, end   % continuous branch: no groups
            m1(k) = mean(y1,'omitnan');  sd1(k) = std(y1,'omitnan');  n1(k) = sum(~isnan(y1));
            m2(k) = mean(y2,'omitnan');  sd2(k) = std(y2,'omitnan');  n2(k) = sum(~isnan(y2));
            dmean(k) = m1(k) - m2(k);
        end

        % Holm-Sidak is named p_, not q_, on purpose: it is an FWER-adjusted
        % P-VALUE (probability of ANY false positive), while the q_ columns are
        % FDR q-values (expected PROPORTION of false positives among those
        % called). Different quantities, not different strengths of the same
        % one - the prefix is there so the table cannot be read as if they were
        % interchangeable.
        % The group summary columns (n_g1, mean_g1, cohens_d, ...) have no meaning
        % for a regression on a continuous predictor, so the continuous branch
        % gets its own columns rather than NaN-filled group ones. Everything
        % after this point - the correction columns, the adjusted block - is
        % shared, because it only ever touches p-values and t-statistics.
        if isequal(sig_test_type,'continuous')
            T = table(nms, [sig_fdr(sel).n_cont]', [sig_fdr(sel).beta]', ...
                      [sig_fdr(sel).t]', [sig_fdr(sel).r_partial]', pu, ...
                'VariableNames', {'signature','n','beta_unadj','t_unadj', ...
                                  'partial_r','p_unadj'});
        else
            T = table(nms, n1, n2, m1, sd1, m2, sd2, dmean, ...
                      [sig_fdr(sel).t]', [sig_fdr(sel).d]', pu, ...
                'VariableNames', {'signature','n_g1','n_g2','mean_g1','sd_g1','mean_g2','sd_g2', ...
                                  'diff_unadj','t_unadj','cohens_d','p_unadj'});
        end
        corr_cols_u = cell(numel(wh_corr),1);
        for z = 1:numel(wh_corr)
            corr_cols_u{z} = [corr_defs{wh_corr(z),3} '_unadj'];
            T.(corr_cols_u{z}) = corr_u{z};
        end
        T.diff_adj = [sig_fdr(sel).diff_adj]';
        T.t_adj    = [sig_fdr(sel).t_adj]';
        T.p_adj    = pa;
        corr_cols_a = {};
        if use_adj
            corr_cols_a = cell(numel(wh_corr),1);
            for z = 1:numel(wh_corr)
                corr_cols_a{z} = [corr_defs{wh_corr(z),3} '_adj'];
                T.(corr_cols_a{z}) = corr_a{z};
            end
        end
        T = sortrows(T, 'p_unadj');
        T.Properties.UserData = struct('contrast', contrasts_done{cc}, ...
            'metric', simnames{1}, 'scaling', scalenames{1}, ...
            'groups', {groupnames}, 'adjusted_for', {adjust_for_covs}, ...
            'pi0_unadj', pi0_u, 'storey_reliable_unadj', info_u.reliable);

        fprintf('\ncontrast: %s   |   metric: %s   |   scaling: %s\n', ...
            contrasts_done{cc}, simnames{1}, scalenames{1});
        % The group framing is wrong under a continuous predictor, and n1/n2 are
        % NaN there because the continuous branch leaves y1/y2 empty on purpose.
        % Printing "high X (n=NaN) vs low X (n=NaN)" above a regression table is
        % exactly the kind of stale label that survives into a manuscript.
        if isequal(sig_test_type,'continuous')
            fprintf('predictor: %s (continuous, n=%d)', sig_covariate_name, ...
                    max([sig_fdr(sel).n_cont]));
        else
            fprintf('groups: %s (n=%d) vs %s (n=%d)', groupnames{1}, n1(1), groupnames{2}, n2(1));
        end
        if use_adj, fprintf('   |   adjusted for: %s', strjoin(adjust_for_covs, ', ')); end
        fprintf('\n\n');
        % full table is wide; print the requested corrections compactly too
        disp(T);
        fprintf('\n  corrections side by side (unadjusted), %d method(s):\n', numel(wh_corr));
        disp(T(:, [{'signature','p_unadj'}, corr_cols_u(:)']));

        if ~info_u.reliable
            fprintf('  NOTE pi0 not identifiable (unadjusted): q_Storey == q_BH, read q_BH\n');
        end
        if use_adj && ~info_a.reliable
            fprintf('  NOTE pi0 not identifiable (adjusted): q_Storey == q_BH, read q_BH\n');
        end

        csvname = fullfile(resultsdir, sprintf('signature_group_diff_%s_%s_%s%s.csv', ...
            simnames{1}, scalenames{1}, matlab.lang.makeValidName(contrasts_done{cc}), ...
            sig_results_tag));
        writetable(T, csvname);
        fprintf('\n  saved: %s\n', csvname);

        SIGSUM = [SIGSUM; T]; %#ok<AGROW>

        % ---- violin plots: one figure unadjusted, one adjusted --------------
        % Two figures rather than one, because the adjusted values are on the
        % same scale but are not the observed data: overlaying them would
        % invite reading the covariate-removed spread as raw variability.
        variants = {'unadjusted'};
        if use_adj && ~isempty(sig_fdr(sel(1)).y1_adj), variants{end+1} = 'adjusted'; end

        for vv = 1:numel(variants)

            isadj = strcmp(variants{vv}, 'adjusted');

            Y1 = cell(1,n); Y2 = cell(1,n); keep = true(1,n);
            for k = 1:n
                if isadj
                    Y1{k} = sig_fdr(sel(k)).y1_adj(:); Y2{k} = sig_fdr(sel(k)).y2_adj(:);
                else
                    Y1{k} = sig_fdr(sel(k)).y1(:);     Y2{k} = sig_fdr(sel(k)).y2(:);
                end
                keep(k) = ~isempty(Y1{k}) && ~isempty(Y2{k});
            end
            if ~any(keep), continue, end
            Y1 = Y1(keep); Y2 = Y2(keep);
            labs = T.signature(keep);
            nk = numel(Y1);

            figttl = sprintf('signature group diffs %s %s %s', simnames{1}, variants{vv}, contrasts_done{cc});
            fh_v = create_figure(figttl);

            xg1 = (1:3:3*nk); xg2 = xg1 + 1;

            violinplot(Y1, 'x', xg1, 'facecolor', groupcolors{1}, 'edgecolor', 'none', ...
                       'facealpha', 0.45, 'mc', 'k', 'medc', [.3 .3 .3], 'pointsize', 8, 'plotlegend', 0);
            violinplot(Y2, 'x', xg2, 'facecolor', groupcolors{2}, 'edgecolor', 'none', ...
                       'facealpha', 0.45, 'mc', 'k', 'medc', [.3 .3 .3], 'pointsize', 8, 'plotlegend', 0);

            % Pad the x-axis. Without this the first violin sits exactly on the
            % axis line at x = 1 and is clipped by it, which reads as a missing
            % group rather than as a clipped one.
            xlim([min(xg1) - 1.5, max(xg2) + 1.5]);

            set(gca, 'XTick', xg1 + 0.5, 'XTickLabel', labs, 'XTickLabelRotation', 45);
            if isadj
                ylabel(sprintf('cosine similarity, adjusted for %s', strjoin(adjust_for_covs, ', ')));
            else
                ylabel(sprintf('cosine similarity (%s)', scalenames{1}));
            end
            title(sprintf('%s: %s vs %s (%s)', contrasts_done{cc}, groupnames{1}, groupnames{2}, variants{vv}), 'FontSize', 14);
            yline(0, ':', 'Color', [.5 .5 .5]);

            % Significance marks come from the matching column, so each figure
            % is annotated with its own test rather than the other one's.
            yl = ylim; ytxt = yl(2) - 0.03*range(yl);
            Tk = T(keep,:);
            for k = 1:nk
                % ** marks whichever correction is the headline one for this
                % run: Storey if requested, else BH, else the first requested
                % method. Hard-coding q_Storey here would error whenever the
                % user switches it off in corrections_i_want.
                if isadj, sfx = '_adj'; praw = Tk.p_adj(k); else, sfx = '_unadj'; praw = Tk.p_unadj(k); end
                star = '';
                for cand = {'q_Storey','q_BH','q_adaptiveFDR','q_BKY','p_holmSidak'}
                    cn = [cand{1} sfx];
                    if ismember(cn, Tk.Properties.VariableNames)
                        if Tk.(cn)(k) < .05, star = '**'; end
                        break
                    end
                end
                if isempty(star) && praw < .05, star = '*'; end
                if ~isempty(star)
                    text(xg1(k)+0.5, ytxt, star, 'HorizontalAlignment','center', 'FontSize', 16);
                end
            end
            xlabel(sprintf('* p < .05   ** q_{Storey} < .05   (%s)', variants{vv}), 'Interpreter', 'tex');

            plugin_set_figure_size('fig', fh_v);
            drawnow, snapnow

        end

    end

    save(fullfile(resultsdir, sprintf('signature_group_diff_%s_%s%s.mat', simnames{1}, scalenames{1}, sig_results_tag)), ...
         'SIGSUM', 'sig_fdr');

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

%% NPS SUBREGIONS: GROUP DIFFERENCES, TABLE AND PLOT
% -------------------------------------------------------------------------
% Third level of the decomposition: the 8 positive and 7 negative NPS regions.
%
% These come from DAT.NPSsubregions, written by prep_4 via apply_nps using
% whatever similarity metric prep_4 ran with. The field names do NOT record the
% metric, so this script must be run in the same session as (or straight after)
% the COSINE prep_4 variant, or the numbers here will silently describe a
% different metric from the signature results above.
%
% Correction is applied across the positive set and the negative set SEPARATELY,
% since each is a family in its own right, and both q_BH and q_Storey are shown
% for the same reason as above: at 7-8 tests pi0 is barely identifiable.

% The NPS subregion pass is still GROUP-ONLY. Its ten tests are two-sample
% comparisons and its tables report group means, none of which is defined for a
% regression on a continuous predictor. Skipping loudly is the honest option:
% running it with a 91-level "group" would either error deep inside
% ttest2_printout or, worse, produce a table that looks valid.
%
% TODO: give the subregion pass the same branch the signature loop now has.
if isequal(sig_test_type,'continuous') && isfield(DAT, 'NPSsubregions')
    fprintf('\n\n');
    printhdr('NPS SUBREGION ANALYSIS SKIPPED');
    fprintf('\n');
    fprintf(['The subregion pass performs two-sample tests and reports group means,\n' ...
             'which are not defined for sig_test_type = ''continuous''. The signature\n' ...
             'family above IS analysed continuously; only this section is skipped.\n']);
end

if ~isequal(sig_test_type,'continuous') && isfield(DAT, 'NPSsubregions')

    fprintf('\n\n');
    printhdr('NPS SUBREGIONS: GROUP DIFFERENCES');
    fprintf('\n\n');

    subr_sets = {'npspos_by_region_contrasts', 'posnames', 'POSITIVE'
                 'npsneg_by_region_contrasts', 'negnames', 'NEGATIVE'};

    NPSSUBSUM = table();

    for ss = 1:size(subr_sets,1)

        fld = subr_sets{ss,1}; nmfld = subr_sets{ss,2}; lbl = subr_sets{ss,3};
        if ~isfield(DAT.NPSsubregions, fld)
            fprintf('  %s: %s not present, skipped\n', lbl, fld); continue
        end

        regnames = DAT.NPSsubregions.(nmfld);

        for i = 1:kc

            [grp, gnames, gcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i);
            if isempty(grp), continue, end
            if size(grp,2) > 1, grp = grp(:,1); end
            if numel(unique(grp)) > 2, grp = mediansplit(grp); end
            if numel(gcolors) < 2, gcolors = seaborn_colors(2); end
            if numel(gnames)  < 2, gnames  = {'High' 'Low'}; end

            mydat = DAT.NPSsubregions.(fld){i};
            nreg  = size(mydat, 2);

            nm = cell(nreg,1); pv = nan(nreg,1); tv = nan(nreg,1); dv = nan(nreg,1);
            m1 = nan(nreg,1); m2 = nan(nreg,1); Y1 = cell(1,nreg); Y2 = cell(1,nreg);

            for rgn = 1:nreg
                y1 = mydat(grp > 0, rgn); y2 = mydat(grp < 0, rgn);
                [~, pp, ~, st] = ttest2(y1, y2);
                nm{rgn} = regnames{rgn};
                pv(rgn) = pp; tv(rgn) = st.tstat;
                dv(rgn) = st.tstat * sqrt(1/numel(y1) + 1/numel(y2));
                m1(rgn) = mean(y1,'omitnan'); m2(rgn) = mean(y2,'omitnan');
                Y1{rgn} = y1; Y2{rgn} = y2;
            end

            % Subregions use corrections_i_want_subregions, NOT the list used
            % for the two-signature table above: 7-8 regions is a family where
            % an FDR is meaningful, two signatures is not.
            wh_sub = find(ismember(corr_defs(:,1), corrections_i_want_subregions));
            [~, pi0_s, info_s] = LaBGAScore_Storey_FDR(pv);
            Tsub = table(nm, m1, m2, m1-m2, tv, dv, pv, ...
                'VariableNames', {'region','mean_g1','mean_g2','diff','t','cohens_d','p'});
            sub_cols = cell(numel(wh_sub),1);
            for z = 1:numel(wh_sub)
                sub_cols{z} = corr_defs{wh_sub(z),3};
                Tsub.(sub_cols{z}) = LaBGAScore_Storey_FDR(pv, 'method', corr_defs{wh_sub(z),2}, 'verbose', false);
            end
            Tsub = sortrows(Tsub, 'p');
            Tsub.set = repmat({lbl}, height(Tsub), 1);
            Tsub.contrast = repmat(DAT.contrastnames(i), height(Tsub), 1);

            fprintf('\n--- NPS %s subregions, %s (%s n=%d vs %s n=%d), pi0 = %.3f ---\n', ...
                lbl, DAT.contrastnames{i}, gnames{1}, sum(grp>0), gnames{2}, sum(grp<0), pi0_s);
            disp(Tsub(:, [{'region','diff','t','cohens_d','p'}, sub_cols(:)']));
            fprintf('  %d of %d at p < .05', sum(Tsub.p < .05), height(Tsub));
            for z = 1:numel(sub_cols)
                fprintf(', %d at %s < .05', sum(Tsub.(sub_cols{z}) < .05), sub_cols{z});
            end
            fprintf('\n');

            NPSSUBSUM = [NPSSUBSUM; Tsub]; %#ok<AGROW>

            % violin panel, ordered as the table
            fh_s = create_figure(sprintf('NPS %s subregions %s', lbl, DAT.contrastnames{i}));
            [~, ord] = ismember(Tsub.region, nm);
            xg1 = (1:3:3*nreg); xg2 = xg1 + 1;
            violinplot(Y1(ord), 'x', xg1, 'facecolor', gcolors{1}, 'edgecolor', 'none', ...
                'facealpha', 0.45, 'mc', 'k', 'medc', [.3 .3 .3], 'pointsize', 8, 'plotlegend', 0);
            violinplot(Y2(ord), 'x', xg2, 'facecolor', gcolors{2}, 'edgecolor', 'none', ...
                'facealpha', 0.45, 'mc', 'k', 'medc', [.3 .3 .3], 'pointsize', 8, 'plotlegend', 0);
            xlim([min(xg1) - 1.5, max(xg2) + 1.5]);
            set(gca, 'XTick', xg1 + 0.5, 'XTickLabel', Tsub.region, 'XTickLabelRotation', 45);
            ylabel(sprintf('NPS %s subregion response', lower(lbl)));
            title(sprintf('NPS %s subregions: %s vs %s', lbl, gnames{1}, gnames{2}), 'FontSize', 14);
            yline(0, ':', 'Color', [.5 .5 .5]);
            yl = ylim; ytxt = yl(2) - 0.03*range(yl);
            for rgn = 1:nreg
                star = ''; if Tsub.q_Storey(rgn) < .05, star = '**'; elseif Tsub.p(rgn) < .05, star = '*'; end
                if ~isempty(star), text(xg1(rgn)+0.5, ytxt, star, 'HorizontalAlignment','center','FontSize',16); end
            end
            xlabel('* p < .05   ** q_{Storey} < .05', 'Interpreter', 'tex');
            plugin_set_figure_size('fig', fh_s);
            drawnow, snapnow

        end
    end

    if ~isempty(NPSSUBSUM)
        csvsub = fullfile(resultsdir, sprintf('NPS_subregion_group_diff_%s_%s.csv', simnames{1}, scalenames{1}));
        writetable(NPSSUBSUM, csvsub);
        fprintf('\n  saved: %s\n', csvsub);
        save(fullfile(resultsdir, sprintf('NPS_subregion_group_diff_%s_%s.mat', simnames{1}, scalenames{1})), 'NPSSUBSUM');
    end

else
    fprintf('\n\nDAT.NPSsubregions not present - run the cosine prep_4 variant first.\n\n');
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


%% LOCAL FUNCTIONS: covariate-adjusted group comparison
% -------------------------------------------------------------------------


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


function s = ternary_str(c, a, b)
% tiny helper: MATLAB has no inline conditional expression
if c, s = a; else, s = b; end
end


function st = h_sig_continuous_test(yv, xv, covmat, covnames, siglabel, contrastname, kc, i)
% h_sig_continuous_test  Regress a signature response on a continuous predictor.
%
% The continuous counterpart of the two-sample branch: same place in the loop,
% same returned fields, so the correction family and the summary table
% downstream do not care which test produced them.
%
% Returns beta and PARTIAL R rather than a mean difference and Cohen's d.
% Partial r is computed from t and df, so it is the effect size of the
% predictor with any nuisance covariates already partialled out, and it is on a
% bounded, comparable scale across signatures.

ok = ~isnan(yv(:)) & ~isnan(xv(:));
if ~isempty(covmat), ok = ok & all(~isnan(covmat), 2); end

Y = yv(ok); X = xv(ok);
if isempty(covmat)
    mdl = fitlm(X, Y);
else
    mdl = fitlm([X, covmat(ok,:)], Y);
end

% Row 2 is the predictor: row 1 is the intercept and any nuisance covariates
% follow it, because they were appended AFTER X above.
b  = mdl.Coefficients.Estimate(2);
t  = mdl.Coefficients.tStat(2);
pv = mdl.Coefficients.pValue(2);
df = mdl.DFE;
r_partial = sign(t) * sqrt(t^2 / (t^2 + df));

subplot(1, kc, i)
scatter(X, Y, 36, 'filled', 'MarkerFaceAlpha', 0.6); hold on
xl = [min(X) max(X)];
if diff(xl) > 0
    plot(xl, mdl.Coefficients.Estimate(1) + b*xl, '-', 'LineWidth', 2);
end
hold off
title(contrastname)
xlabel('predictor');
ylabel(sprintf('%s Response', siglabel));

% fprintf, not printstr: printstr is an anonymous function created in
% a_set_up_paths_always_run_first, i.e. a SCRIPT-SCOPE VARIABLE, and local
% functions cannot see script variables. Calling it here fails with
% "Undefined function 'printstr'" even though the script body uses it freely.
fprintf('Regression on continuous predictor:\n');
fprintf('  n = %d, beta = %.4f, t(%d) = %.3f, p = %.6f, partial r = %.3f\n', ...
        sum(ok), b, df, t, pv, r_partial);
if ~isempty(covnames)
    fprintf('  adjusted for: %s\n', strjoin(covnames, ', '));
end

st = struct('p', pv, 't', t, 'df', df, 'beta', b, 'r_partial', r_partial, ...
            'n', sum(ok), 'x', X, 'y', Y);
end


function v = h_sig_field(st, f)
% h_sig_field  Value of a field of the continuous-test struct, or NaN.
% The sig_fdr record is built UNCONDITIONALLY for both test types, so the
% continuous-only fields must resolve to something in the group branch too.
if isempty(st) || ~isstruct(st) || ~isfield(st, f)
    v = NaN;
else
    v = st.(f);
end
end


function xv = h_sig_get_covariate(DAT, mygroupnamefield, i, covname)
% h_sig_get_covariate  Pull the continuous predictor for contrast/condition i.
T = DAT.BETWEENPERSON.(mygroupnamefield){i};
if ~istable(T)
    error('DAT.BETWEENPERSON.%s{%d} is not a table, so ''%s'' cannot be read.', ...
          mygroupnamefield, i, covname);
end
if ~ismember(covname, T.Properties.VariableNames)
    error(['sig_covariate_name ''%s'' is not a column of DAT.BETWEENPERSON.%s{%d}.\n' ...
           'available columns are: %s'], covname, mygroupnamefield, i, ...
           strjoin(T.Properties.VariableNames, ', '));
end
xv = double(T.(covname));
xv = xv(:);
end


function [cm, cn] = h_sig_get_covmat(DAT, mygroupnamefield, i, covnames)
% h_sig_get_covmat  Nuisance covariate matrix for the continuous regression.
% Mirrors what h_adjusted_group_test does for the group branch, so
% adjust_for_covs means the same thing under both test types.
cm = []; cn = {};
if isempty(covnames), return, end
T = DAT.BETWEENPERSON.(mygroupnamefield){i};
missing = covnames(~ismember(covnames, T.Properties.VariableNames));
if ~isempty(missing)
    error('adjust_for_covs names column(s) not in the design: %s', strjoin(missing, ', '));
end
cm = zeros(height(T), numel(covnames));
for k = 1:numel(covnames)
    cm(:,k) = double(T.(covnames{k}));
end
cn = covnames;
end
