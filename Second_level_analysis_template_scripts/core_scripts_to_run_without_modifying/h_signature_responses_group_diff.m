%% h_signature_responses_group_diff

% Runs a two-sample t-test for a set of signature responses, for each
% contrast.  Also runs NPS subregions. 
% Signatures, scaling, etc. are defined below.


%% USER OPTIONS
% -------------------------------------------------------------------------

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
                printstr(sprintf('Group differences: %s, %s', signature, DAT.contrastnames{i}));
                printstr(dashes)

                barplot_columns(y, 'nofig', 'colors', groupcolors, 'names', groupnames);

                title(DAT.contrastnames{i})
                xlabel('Group');
                ylabel(sprintf('%s Response', signature));

                printstr('Between-groups test:');

                [H,p,ci,stats] = ttest2_printout(y{1}, y{2});

                printstr(dashes)

            end % panels

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
    
    end % if loop group or individual signature
    
end % signature


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