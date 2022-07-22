% get holdout set for SVM contrasts, balanced across groups, 2nd-level batch scripts
% adapted from plugin_get_holdout_sets by @lukasvo76, March 2021

% hard-coded options
% 'fivefold_leave_whole_subject_out' : 5-fold cross-validation, leaving
% whole subject out (all images from subject across all conditions)
% 'leave_one_subject_out' : leave one whole subject out

% holdout_set_type = 'fivefold_leave_whole_subject_out';
% nfolds = 5;

% @lukasvo76, July 2022
% removed hard coded options, these are now set in
% a2_set_default_options.m or prep_3c_run_SVMs_on_contrasts_masked.m for
% more flexibility

outcome_value = zeros(size(condition_codes));
holdout_set = {};

switch holdout_set_type
    
% @lukasvo76: did not adapt this option since LOOCV is generally not
% recommended anymore, hence untested
    case 'leave_one_subject_out'
        
        printhdr('Holdout set: leave one whole subject out');
        
        for i = 1:length(wh)
            
            n = sum(condition_codes == i);
            holdout_set{i} = [1:n]';
            
            outcome_value(condition_codes == i) = sign(mycontrast(wh(i)));
            
        end
        
    case 'kfold'
        
        clear n
        
        printhdr('Holdout set:  5-fold cross-validation, leaving whole subject out');
        
        for i = 1:size(wh,2) % wh is which conditions have non-zero contrast weights
            
            n(i) = sum(condition_codes == i);
            
            outcome_value(condition_codes == i) = sign(mycontrast(wh(i)));
            
        end
        
        % define group identifiers for balancing
        if ~isempty(DAT.BETWEENPERSON.group)
            
            group = DAT.BETWEENPERSON.group;
        
        elseif ~isempty(DAT.BETWEENPERSON.contrasts{c}.group)
            
            group = DAT.BETWEENPERSON.contrasts{c}.group;
        
        else
            
            error('\ngroup identifier not defined in DAT.BETWEENPERSON.group or DAT.BETWEENPERSON.contrasts{c}.group, incompatible with "group" method defined in holdout_set_method option\n');
            
        end
       
        % create cvpartition object
        % @lukasvo76: 
        % 1) cvpartition2 is @bogpetre's version of the standard cvpartition
        %    object, with improved stratification options, which we do not
        %    need here - see
        %    LaBGAScore_secondlevel_run_MVPA_regression_single_trial for an
        %    example of its usage
        % 2) group balances holdout sets over groups, stratify option
        %    can be used with subject id can be use to ensure whole subjects
        %    are left out, but the original code below takes care of that
        %    too, so we don't need that option here
        
        cvpart = cvpartition(group,'KFOLD',nfolds);
        
        for i = 1:nfolds
            
            mytest = cvpart.test(i);
            
            for j = 1:length(wh)
                
                holdout_set{j}(mytest(1:n(j))) = i;  % works if some conditions have fewer images, though we expect them to match
                holdout_set{j} = holdout_set{j}';
            end
            
        end
            
    otherwise
        
        error('\ninvalid option "%s" defined in holdout_set_type variable, choose between "kfold" and "leave_one_subject_out"\n',holdout_set_type);

end


holdout_set = cat(1, holdout_set{:});