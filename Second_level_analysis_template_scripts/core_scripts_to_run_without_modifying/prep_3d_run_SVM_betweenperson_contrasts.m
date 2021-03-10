% THIS SCRIPT RUNS BETWEEN-PERSON CONTRASTS
% Assuming that groups are concatenated into contrast image lists.
% Requires DAT.BETWEENPERSON.group field specifying group membership for
% each image.
% --------------------------------------------------------------------


%% Now set in a2_set_default_options
%--------------------------------------------------------------------------
options_needed = {'dosavesvmstats', 'dobootstrap', 'boot_n','myscaling_svm_between'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {true false 1000 'raw'};          % defaults if we cannot find info in a2_set_default_options at all 

plugin_get_options_for_analysis_script


%% Check for required DAT fields. Skip analysis and print warnings if missing.
% ---------------------------------------------------------------------

% List required fields in DAT, in cell array:
required_fields = {'BETWEENPERSON', 'contrastnames', 'contrasts' 'contrastcolors'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end

spath = which('use_spider.m');
if isempty(spath)
    disp('Warning: spider toolbox not found on path; prediction may break')
end

if dobootstrap, svmtime = tic; end


%% Get mask
% -------------------------------------------------------------------------

if exist('maskname_svm', 'var') && ~isempty(maskname_svm)
    
    disp('Masking data')
    svmmask = fmri_data(maskname_svm, 'noverbose');
    
else
    
    disp('No mask found; using full original image data');

end


%% Run between-person SVM for each contrast
% --------------------------------------------------------------------

printhdr('Cross-validated SVM to discriminate between-person contrasts');

kc = size(DAT.contrasts, 1);

svm_stats_results = cell(1, kc);

for c = 1:kc
    
    mygroupnamefield = 'contrasts';  % 'conditions' or 'contrasts'
    [group, ~, ~] = plugin_get_group_names_colors(DAT, mygroupnamefield, c);
    group = group(:,1); % @lukasvo76: added to correctly handle the situation where you added other covariates in addition to group in DAT.BETWEENPERSON.mygroupnamefield{c}
    outcome_value = group;
    
    if isempty(group)
        fprintf('Group not defined for contrast %s. Skipping.\n', DAT.contrastnames{c}); 
        continue
    end
    
    % Define holdout sets: 5-folds, balanced over groups
    % @lukasvo76: changed from the original LOOCV code
    % Assume that subjects are in same position in each input file!
    % --------------------------------------------------------------------
    
    nfolds = 5; % @lukasvo76: hardcoded for now, but can easily be adapted
    cvpart = cvpartition(group,'KFOLD',nfolds);
    [I,J] = find([cvpart.test(1),cvpart.test(2), cvpart.test(3), cvpart.test(4), cvpart.test(5)]);
    holdout_set = sortrows([I,J]);
    holdout_set = holdout_set(:,2);
    
    printstr(DAT.contrastnames{c});
    printstr(dashes)
    
    mycontrast = DAT.contrasts(c, :);
    wh = find(mycontrast);
    
    % Select data for this contrast
    % --------------------------------------------------------------------
    
    switch myscaling_svm_between
        case 'raw'
            printstr('Raw (unscaled) images used in between-person SVM');
            scaling_string = 'no_scaling';
            cat_obj = DATA_OBJ_CON{c};
            
        case 'scaled'
            printstr('Z-scored images used in between-person SVM');
            scaling_string = 'scaling_z_score_conditions';
            cat_obj = DATA_OBJ_CONsc{c};
            
        case 'scaled_contrasts'
            printstr('l2norm scaled contrast images used in between-person GLM');
            scaling_string = 'scaling_l2norm_contrasts';
            cat_obj = DATA_OBJ_CONscc{c};
            
        otherwise
            error('myscaling must be ''raw'' or ''scaled'' or ''scaled_contrasts''');
    end
    
    % Apply mask if specified in a2_ script
    %----------------------------------------------------------------------
    if exist('svmmask', 'var')
        
        disp('Masking data')
        cat_obj = apply_mask(cat_obj, svmmask);
        
    else
        
        disp('No mask found; using full existing image data');
        
    end
    
    % Format and attach outcomes: 1, -1 for pos/neg contrast values
    % --------------------------------------------------------------------
 
    cat_obj.Y = outcome_value;
    
    % Skip if necessary
    % --------------------------------------------------------------------
    
    if all(cat_obj.Y > 0) || all(cat_obj.Y < 0)
        % Only positive or negative weights - nothing to compare
        
        printhdr(' Only positive or negative weights - nothing to compare');
        
        continue
    end
    
    % Run prediction model
    % --------------------------------------------------------------------
    if dobootstrap
        [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, 'bootsamples', boot_n, 'error_type', 'mcr', parallelstr);
        % Threshold, if possible - can re-threshold later with threshold() method
        stats.weight_obj = threshold(stats.weight_obj, .05, 'unc'); 
        
    else
        [cverr, stats, optout] = predict(cat_obj, 'algorithm_name', 'cv_svm', 'nfolds', holdout_set, 'error_type', 'mcr', parallelstr);
    end
    
    % Save stats objects for results later
    % --------------------------------------------------------------------
    
%     stats.weight_obj = enforce_variable_types(stats.weight_obj);
    svm_stats_results{c} = stats;
    
    if exist('svmmask', 'var')
    
        svm_stats_results{c}.mask = svmmask;
        svm_stats_results{c}.maskname = maskname_svm;
    
    end
        
    if dobootstrap, disp('Cumulative run time:'), toc(svmtime); end  
    
end  % between-person contrast


%% Save
% --------------------------------------------------------------------
if dosavesvmstats
    
    savefilenamedata = fullfile(resultsdir,['svm_stats_results_betweenperson_contrasts_',scaling_string,'.mat']);

    save(savefilenamedata, 'svm_stats_results', '-v7.3');
    printhdr('Saved svm_stats_results for contrasts');
    
end



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

end
