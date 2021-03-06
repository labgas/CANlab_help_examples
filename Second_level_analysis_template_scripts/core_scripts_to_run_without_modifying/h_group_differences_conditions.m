% Runs a two-sample t-test for a set of signature responses, for each
% conditions.  Also runs NPS subregions. 
% Signatures, scaling, etc. are defined below.
% created by @lukasvo76 from h_group_differences.m, March 2021


%% USER OPTIONS
% -------------------------------------------------------------------------

mysignature =   {'NPS', 'NPSpos', 'NPSneg', 'SIIPS'};   % 'NPS' 'NPSpos' 'NPSneg' 'SIIPS' etc.  See load_image_set('npsplus')
scalenames =    {'raw'};                                % or scaled
simnames =      {'dotproduct'};                         % or 'cosine_sim' 'dotproduct'
mygroupnamefield = 'conditions';  % 'conditions' or 'contrasts'


%% DEFINE GROUPS IN PREP_1b_PREP_BEHAVIORAL_DATA
% -------------------------------------------------------------------------

% There are two ways to define groups. The other is condition- and
% contrast-specific.  These are entered in DAT.BETWEENPERSON.conditions and
% DAT.BETWEENPERSON.contrasts, in cells.  1 -1 codes work best. 
% These are set up in prep_1b_prep_behavioral_data.m
% If they are missing, using DAT.BETWEENPERSON.group will be used as a
% generic option.
% -------------------------------------------------------------------------


%% LOOP THROUGH SIGNATURES, TEST GROUP DIFFERENCE, CREATE ONE PLOT PER CONTRAST
% -------------------------------------------------------------------------
for s = 1:length(mysignature)
    
    % Get data
    % -------------------------------------------------------------------------
    conditiondata = table2array(DAT.SIG_conditions.(scalenames{1}).(simnames{1}).(mysignature{s}));
    
    kc = size(conditiondata, 2);
    
    % Plot
    % -------------------------------------------------------------------------
    printhdr(sprintf('%s responses: Scale = %s Metric = %s', mysignature{s}, scalenames{1}, simnames{1}));
    
    figtitle = sprintf('%s group diffs %s %s', mysignature{s}, scalenames{1}, simnames{1});
    create_figure(figtitle, 1, kc);
    
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
        
        %y = {DAT.(myfield){i}(group > 0) DAT.(myfield){i}(group < 0)};
        y = {conditiondata(group > 0, i) conditiondata(group < 0, i)};
        
        subplot(1, kc, i)
        
        printstr(' ');
        printstr(sprintf('Group differences: %s, %s', mysignature{s}, DAT.conditions{i}));
        printstr(dashes)
        
        barplot_columns(y, 'nofig', 'colors', groupcolors, 'names', groupnames);
        
        title(DAT.conditions{i})
        xlabel('Group');
        ylabel(sprintf('%s Response', mysignature{s}));
        
        printstr('Between-groups test:');
        
        [H,p,ci,stats] = ttest2_printout(y{1}, y{2});
        
        printstr(dashes)
        
    end % panels
    
    drawnow, snapnow
    
end % signature


%% NPS SUBREGIONS, TEST GROUP DIFFERENCE, CREATE ONE PLOT PER CONTRAST
% -------------------------------------------------------------------------

% POSITIVE
% --------

% which variables to use
mysubrfield = 'npspos_by_region'; % 'npspos_by_region_cosinesim';     %'npspos_by_regionsc';
mysubrfieldneg = 'npsneg_by_region'; % 'npsneg_by_region_cosinesim';  % 'npsneg_by_regionsc';

posnames = DAT.NPSsubregions.posnames;
negnames = DAT.NPSsubregions.negnames;

clear means p T

for i = 1:kc  % for each contrast
    
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
        
    mydat = DAT.NPSsubregions.(mysubrfield){i};
    k = size(mydat, 2);
    
    create_figure(sprintf('NPS subregions by group %s', DAT.conditions{i}), 1, k);
    pos = get(gcf, 'Position');
    pos(4) = pos(4) .* 2.5;
    set(gcf, 'Position', pos);
    
    clear means p T
    
    for j = 1:k  % for each subregion
        
        subplot(1, k, j);
        
        y = {mydat(group == 1, j) mydat(group == -1, j)};
        
        printhdr(posnames{j});
        
        barplot_columns(y, 'nofig', 'colors', groupcolors, 'noviolin', 'noind', 'names', groupnames );
        
        title(posnames{j})
        xlabel('Group');
        if j == 1, ylabel('NPS Response'); end
        
        printstr('Between-groups test:');
        [H,p(j, 1),ci,stats] = ttest2_printout(y{1}, y{2});
        
        means(j, :) = stats.means;
        T(j, 1) = stats.tstat;
        
    end
    
    drawnow, snapnow
    
    % Print between-subject Table
    printhdr('Between-group tests');
    Region = posnames';
    regionmeans = table(Region, means, T, p);
    
    disp(regionmeans);
    
end % panels

% NEGATIVE
% --------

for i = 1:kc
    
    [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i);
    
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
    
        
    if isempty(group), continue, end % skip this condition/contrast - no groups
    
    mydat = DAT.NPSsubregions.(mysubrfieldneg){i};
    k = size(mydat, 2);
    
    create_figure(sprintf('NPS neg subregions by group %s', DAT.conditions{i}), 1, k);
    pos = get(gcf, 'Position');
    pos(4) = pos(4) .* 2.5;
    set(gcf, 'Position', pos);
    
    clear means p T
    
    for j = 1:k
        
        subplot(1, k, j);
        
        y = {mydat(group == 1, j) mydat(group == -1, j)};
        
        printhdr(negnames{j});
        
        barplot_columns(y, 'nofig', 'colors', groupcolors, 'noviolin', 'noind', 'names', groupnames );
        
        title(negnames{j})
        xlabel('Group');
        if j == 1, ylabel('NPS Response'); end
        
        printstr('Between-groups test:');
        [H,p(j, 1),ci,stats] = ttest2_printout(y{1}, y{2});
        
        means(j, :) = stats.means;
        T(j, 1) = stats.tstat;
        
    end
    
    drawnow, snapnow
    
    % Print between-subject Table
    printhdr('Between-group tests');
    Region = negnames';
    regionmeans = table(Region, means, T, p);
    
    disp(regionmeans);
    
end