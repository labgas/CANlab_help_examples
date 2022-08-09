% get holdout set for SVM contrasts, 2nd-level batch scripts

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

switch holdout_set_type_svm
    
    case 'leave_one_subject_out'
        
        printhdr('Holdout set: leave one whole subject out');
        
        for i = 1:length(wh)
            
            n = sum(condition_codes == i);
            holdout_set{i} = [1:n]';
            
            outcome_value(condition_codes == i) = sign(mycontrast(wh(i)));
            
        end
        
    case 'kfold'
        
        clear n
        
        printhdr('Holdout set:  k-fold cross-validation, leaving whole subject out');
        fprintf('\nk = %d\n',nfolds_svm);
        
        for i = 1:length(wh) % wh is which conditions have non-zero contrast weights
            
            n(i) = sum(condition_codes == i);
            holdout_set{i} = zeros(n(i), 1);
            
            outcome_value(condition_codes == i) = sign(mycontrast(wh(i)));
            
        end
        
        % Take largest set and stratify into k folds
        % Or, if equal size sets, pick first
        [~, wh_max_imgs] = max(n);
        
        cvpart = cvpartition(n(wh_max_imgs),'k',nfolds_svm);
        
        % Assign test set for each fold to all images in paired conditions
        % If participants are crossed with conditions, this leaves one whole participant out (all images across conditions)
        
        for i = 1:nfolds_svm
            
            mytest = cvpart.test(i);
            
            for j = 1:length(wh)
                
                holdout_set{j}(mytest(1:n(j))) = i;  % works if some conditions have fewer images, though we expect them to match
                
            end
            
        end
        
    otherwise
        
        error('\ninvalid option "%s" defined in holdout_set_type variable, choose between "kfold" and "leave_one_subject_out"\n',holdout_set_type);

end


holdout_set = cat(1, holdout_set{:});