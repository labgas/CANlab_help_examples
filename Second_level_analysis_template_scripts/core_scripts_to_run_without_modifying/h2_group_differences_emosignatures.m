% Runs a two-sample t-test for a set of signature responses, for each
% contrast.
% Signatures, scaling, etc. are defined below.


%% USER OPTIONS
% -------------------------------------------------------------------------

mysignature =   {'PINES', 'Rejection', 'VPS', 'VPS_nooccip'};
scalenames =    {'raw'};                                % or scaled
simnames =      {'dotproduct'};                         % or 'cosine_sim' 'dotproduct'
mygroupnamefield = 'contrasts';  % 'conditions' or 'contrasts'


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
    contrastdata = table2array(DAT.SIG_contrasts.(scalenames{1}).(simnames{1}).(mysignature{s}));
    
    kc = size(contrastdata, 2);
    
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
        
    end % panels
    
    drawnow, snapnow
    
end % signature