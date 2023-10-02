% Check for required DAT fields. Skip analysis and print warnings if missing.
% ---------------------------------------------------------------------
% List required fields in DAT, in cell array:
required_fields = {'conditions', 'colors', 'SIG_conditions'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end

% Controlling for group admin order covariate (mean-centered by default)

group = [];
if isfield(DAT, 'BETWEENPERSON') && isfield(DAT.BETWEENPERSON, 'group')
    group = DAT.BETWEENPERSON.group; % empty for no variable to control/remove
    
    printstr('Controlling for data in DAT.BETWEENPERSON.group');
end

% Format: The prep_4_apply_signatures_and_save script extracts signature responses and saves them.
%
% These fields contain data tables whose columns are conditions or contrasts, 
% with variable names based on DAT.conditions or DAT.contrastnames, 
% but with spaces replaced with underscores:
% DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
% DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
%
% Convert these to numerical arrays using table2array:
% table2array(DAT.SIG_contrasts.scaled.dotproduct.NPSneg)
%
% DAT.SIG_conditions.raw.dotproduct = apply_all_signatures(DATA_OBJ, 'conditionnames', DAT.conditions);
% DAT.SIG_contrasts.raw.dotproduct = apply_all_signatures(DATA_OBJ_CON, 'conditionnames', DAT.conditions);

k = length(DAT.conditions);

for sig = 1:size(signatures_to_plot,2)
    
    mysignames{sig} = strcat(signatures_to_plot{sig}{:});
    nplots{sig} = size(signatures_to_plot{sig},2);
    mypointsize{sig} = get_point_size(nplots{sig}, k);
    
end

myfontsize = get_font_size(k);
myaxislabels = format_strings_for_legend(DAT.conditions);

%% Signature Response - conditions
% ------------------------------------------------------------------------

for sig = 1:size(keyword_sigs,2)
    
    [~,signame] = fileparts(char(keyword_sigs{sig}));

    figtitle = sprintf('%s CONDITIONS %s %s %s', mysignames{sig}, upper(myscaling_sigs), upper(similarity_metric_sigs), upper(signame));

    fprintf('\n\n');
    printhdr(figtitle);
    fprintf('\n\n');

    fighan = create_figure(figtitle, 1, nplots{sig});
    adjust_fig_position_for_long_xlabel_names(fighan);

    clear axh

    for n = 1:nplots{sig}

        printhdr(['SIGNATURE: ', signatures_to_plot{sig}{n}]);

        axh(n) = subplot(1, nplots{sig}, n);

        mysignature = signatures_to_plot{sig}{n};
        
        mydata = table2array(DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(signame).(mysignature));

        % get condition-specific group covariate data if any
        % this will not work yet, because barplot_columns cannot handle
        % condition-specific covariates. See "group diffs" script.
        %     mygroupnamefield = 'conditions';  % 'conditions' or 'contrasts'
        %     [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, n);

        handles = barplot_columns(mydata, 'title', format_strings_for_legend(mysignature), 'colors', DAT.colors, 'MarkerSize', mypointsize{sig}, 'dolines', 'nofig', 'names', myaxislabels, 'covs', group, 'wh_reg', 0);

        set(axh(n), 'FontSize', myfontsize);
        if n == 1
            ylabel(format_strings_for_legend(similarity_metric_sigs));
        else
            ylabel(' ');
        end
        xlabel('');

        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

    end

    
    % kludgy_fix_for_y_axis(axh); % lukasvo76: commented out because causing
    % errors in some situations

    % plugin_save_figure;
    % close

end


%% Signature Response - contrasts
% ------------------------------------------------------------------------

% Check for required DAT fields. Skip analysis and print warnings if missing.
% ---------------------------------------------------------------------
% List required fields in DAT, in cell array:
required_fields = {'contrasts', 'contrastnames', 'contrastcolors', 'SIG_contrasts'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end
% ------------------------------------------------------------------------

kc = size(DAT.contrasts, 1);
myfontsize = get_font_size(kc);
myaxislabels = format_strings_for_legend(DAT.contrastnames);

for sig = 1:size(keyword_sigs,2)
    
    mypointsize{sig} = get_point_size(nplots{sig}, kc);
    
    [~,signame] = fileparts(char(keyword_sigs{sig}));

    figtitle = sprintf('%s CONTRASTS %s %s %s', mysignames{sig}, upper(myscaling_sigs), upper(similarity_metric_sigs), upper(signame));

    fprintf('\n\n');
    printhdr(figtitle);
    fprintf('\n\n');

    fighan = create_figure(figtitle, 1, nplots{sig});
    adjust_fig_position_for_long_xlabel_names(fighan);

    clear axh

    for n = 1:nplots{sig}

        printhdr(['SIGNATURE: ', signatures_to_plot{sig}{n}]);

        axh(n) = subplot(1, nplots{sig}, n);

        mysignature = signatures_to_plot{sig}{n};
        mydata = table2array(DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(signame).(mysignature));

        % get contrast-specific group covariate data if any
        % this will not work yet, because barplot_columns cannot handle
        % condition-specific covariates. See "group diffs" script.
        %     mygroupnamefield = 'contrasts';  % 'conditions' or 'contrasts'
        %     [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, n);

        handles = barplot_columns(mydata, 'title', format_strings_for_legend(mysignature), 'colors', DAT.contrastcolors, 'MarkerSize', mypointsize{sig}, 'nofig', 'names', myaxislabels, 'covs', group, 'wh_reg', 0);

        set(axh(n), 'FontSize', myfontsize);
        if n == 1
            ylabel(format_strings_for_legend(similarity_metric_sigs)); 
        else
            ylabel(' '); 
        end
        xlabel('');

        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

    end

    % kludgy_fix_for_y_axis(axh); % lukasvo76: commented out because causing
    % errors in some situations

    % plugin_save_figure;
    %close

end

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Sub-functions
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%% Get font size
% -------------------------------------------------------------------------

function myfontsize = get_font_size(k)
% get font size
if k > 10
    myfontsize = 8;
elseif k > 8
    myfontsize = 9;
elseif k > 5
    myfontsize = 12;
elseif k > 2
    myfontsize = 14;
else
    myfontsize = 18;
end
end

%% Adjust figure position
% -------------------------------------------------------------------------

function adjust_fig_position_for_long_xlabel_names(fighan)

figpos = get(fighan, 'Position');
figpos(3) = figpos(3) - 100;
figpos(4) = figpos(4) + 100;
set(fighan, 'Position', figpos);

end

%% Fix y-axis values
% -------------------------------------------------------------------------

function kludgy_fix_for_y_axis(axh)
% Matlab is having some trouble with axes for unknown reasons

axis2 = get(axh(2), 'Position');

% re-set axis 1
mypos = get(axh(1), 'Position');
mypos([2 4]) = axis2([2 4]);  % re-set y start and height
set(axh(1), 'Position', mypos);

end


%% set point size
% -------------------------------------------------------------------------
function ptsize = get_point_size(n, k)

ptsize = 12 ./ (.5*n*log(1 + k));

end

% Not used - post hoc setting
function set_point_size(handles, n, k)

myhandles = handles.point_han(:);
myhandles(cellfun(@isempty, myhandles)) = [];
ptsize = get_point_size(n, k);
ptfun = @(x) set(x, 'MarkerSize', ptsize);
cellfun(ptfun, myhandles);

end


%% Get group between-person vectors of covariates and other info
% -------------------------------------------------------------------------

function [group, groupnames, groupcolors] = plugin_get_group_names_colors(DAT, mygroupnamefield, i)

group = []; groupnames = []; groupcolors = [];

if isfield(DAT, 'BETWEENPERSON') && ...
        isfield(DAT.BETWEENPERSON, mygroupnamefield) && ...
        iscell(DAT.BETWEENPERSON.(mygroupnamefield)) && ...
        length(DAT.BETWEENPERSON.(mygroupnamefield)) >= i && ...
        ~isempty(DAT.BETWEENPERSON.(mygroupnamefield){i})
    
    group = DAT.BETWEENPERSON.(mygroupnamefield){i};
    
elseif isfield(DAT, 'BETWEENPERSON') && ...
        isfield(DAT.BETWEENPERSON, 'group') && ...
        ~isempty(DAT.BETWEENPERSON.group)
    
    group = DAT.BETWEENPERSON.group;

end


if isfield(DAT, 'BETWEENPERSON') && isfield(DAT.BETWEENPERSON, 'groupnames')
    groupnames = DAT.BETWEENPERSON.groupnames;
elseif istable(group)
    groupnames = group.Properties.VariableNames(1);
else
    groupnames = {'Group-Pos' 'Group-neg'};
end

if isfield(DAT, 'BETWEENPERSON') && isfield(DAT.BETWEENPERSON, 'groupcolors')
    groupcolors = DAT.BETWEENPERSON.groupcolors;
else
    groupcolors = seaborn_colors(2);
end

if istable(group), group = table2array(group); end

if ~isempty(group), disp('Controlling for between-image group variable.'); end

end
