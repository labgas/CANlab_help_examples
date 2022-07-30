%%% c2f_run_MVPA_regression_single_trial.m
%
% USAGE
%
% This script runs MVPA regression analysis on a continuous outcome Y
% (lasso-pcr, but can easily be adapted to support vector regression or
% other methods available in CANlab's predict function) on an fmri_data_st
% object created using prep_3c_run_SVMs_on_contrasts_masked.m
% That script should be run first, or the present script will load the data
% object if it is saved by the previous script.
%
% Options for this script are set in a2_set_default_options.m, see that
% script for more info. Many of these options get passed into CANlab's
% predict function - help predict in Matlab command window for more info
%
% TUTORIALS AND DOCUMENTATION
%
% A. Canlab predict function
%
% This is the classic CANlab method of running ML models, and can be chosen
% by setting the ml_method_mvpa_reg_st option in a2_set_default_options.m
% to 'predict'
%
% This script is based on the extremely helpful tutorials on 
% single trial MVPA analysis by @bogpetre @CANlab.
%
% Here are Bogdan's walkthroughs:
% https://canlab.github.io/_pages/canlab_single_trials_demo/demo_norming_comparison.html
% https://canlab.github.io/_pages/mlpcr_demo/mlpcr_demo.html (WiP)
%
% Here are two scripts @lukasvo76 adapted from these walkthroughs
% https://www.dropbox.com/sh/e17nl3ew1db1twk/AACO9QAEt6Sy3TejH-n-tbdEa?dl=0
% https://www.dropbox.com/sh/bm0at2dr81isk70/AABD67D_bF8A0NFa4gtt2dHNa?dl=0
% 
% Another highly helpful resource in this context is this Nature Methods
% paper by Tor and Wani Woo
% https://www.nature.com/articles/s41596-019-0289-5
%
% @lukasvo76's version of the script for this paper can be found here
% https://www.dropbox.com/sh/v2nsgoqmbi0cqnk/AAD6I1Gn5KUM6aViom4TLeVJa?dl=0
%
% B. Bogdan's object-oriented ML method for fmri_data objects (and beyond)
%
% This is a newer method inspired by Python's scikit-learn, including more
% flexible options for algorithm and feature selection, 
% hyperparameter optimization, nested cross-validation, etc. However, it
% does require more advanced programming skills and understanding the logic
% of the method, with only example models provided as part of this script
%
% Dependency: https://github.com/canlab/ooFmriDataObjML
%
% Tutorial: https://canlab.github.io/_pages/canlab_pipelines_walkthrough/estimateBestRegionPerformance.html
% Example script: https://github.com/labgas/LaBGAScore/blob/main/secondlevel/LaBGAScore_secondlevel_ooFmriDataObjML_example.m
% 
% NOTE: This script is work in progress, particularly the permutation,
% multilevel, and oofmridataobj options are still undergoing improvement
% and full testing
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   April, 2021
%__________________________________________________________________________
% @(#)% c2f_run_MVPA_regression_single_trial     v3.1        
% last modified: 2022/07/30


%% LOAD FMRI_DATA_ST OBJECT AND OTHER NECESSARY VARIABLES IF NEEDED
%--------------------------------------------------------------------------

if ~exist('resultsdir','var')
    
    a_set_up_paths_always_run_first
    
end

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
    if ~isfield(DAT,'BEHAVIOR')
        error('\n Behavioral data not yet added to DAT structure - run prep_1b script first')
    end
    
end

if ~exist('fmri_dat','var')
    
    load(fullfile(resultsdir,['single_trial_fmri_data_st_object_' DSGN.modelingfilesdir '.mat']));
    
end

% specify which montage to add title to

whmontage = 5; % see region.montage docs


%% DEFINE SUBJECT IDENTIFIERS
%--------------------------------------------------------------------------

subject_id = fmri_dat.metadata_table.(subj_identifier);
[uniq_subject_id, ~, subject_id] = unique(subject_id,'stable');
fmri_dat.metadata_table.subject_id = subject_id;
n_subj = size(uniq_subject_id,1);


%% SCALE AND/OR MASK IMAGES AND BEHAVIORAL OUTCOME ACCORDING TO OPTIONS
%--------------------------------------------------------------------------

% SCALING IMAGES
%---------------

switch myscaling_mvpa_reg_st

    case 'raw'
        
        fprintf('\nNo scaling of input images\n');
        
    case 'centerimages'

        fmri_dat = fmri_dat.rescale('centerimages'); 
        fprintf('\nCentering input images\n');
        
    case 'l2norm_images'

        fmri_dat = fmri_dat.rescale('l2norm_images');
        fprintf('\nNormalizing input images by l2norm\n');
        
    case 'zscore_images'

        fmri_dat_z = fmri_dat.rescale('zscoreimages');
        fprintf('\nZ-scoring input images\n');
        
    otherwise 

        error('\ninvalid scaling option %s specified in myscaling_mvpa_reg_st variable defined in a2_set_default_options CANlabhelpexamples script, please correct\n', myscaling_mvpa_reg_st)

end % switch scaling
   
% MASK
%-----

if exist('maskname_mvpa_reg_st','var') && ~isempty(maskname_mvpa_reg_st) && exist(maskname_mvpa_reg_st, 'file')
    
    [~,maskname_short] = fileparts(maskname_mvpa_reg_st);
    fprintf('\nMasking data with %s\n',maskname_short);
    mask_string = sprintf('within mask %s', maskname_short);
    
    mvpamask = fmri_mask_image(maskname_mvpa_reg_st);
    fmri_dat = fmri_dat.apply_mask(mvpamask);
    fmri_dat.mask_descrip = maskname_mvpa_reg_st;
    
else
    
    fprintf('\nNo mask found; using full original image data\n');

end % if loop mask

% ZSCORE BEHAVIORAL OUTCOME
%--------------------------

% NOTE: useful for more interpretable values of prediction MSE

if zscore_outcome
   
    fmri_dat.Y = zscore(fmri_dat.Y);
    fprintf('\nZ-scoring outcome fmri_dat.Y across subjects\n');
    
end


%% DATA VISUALISATION PRIOR TO MODEL BUILDING
%--------------------------------------------------------------------------

% BETA IMAGES

h1=figure;

    for sub = 1:n_subj
        subj_idx = sub == subject_id;
        this_subj_dat = fmri_dat.dat(:,subj_idx);
        q(sub,:) = quantile(this_subj_dat(:),[0.025,0.5,0.975]);
        mu = mean(mean(this_subj_dat(:)));
        sd = std(this_subj_dat(:));
        h1 = plot([mu-sd, mu+sd],[sub,sub],'-');
        hold on;
        h2 = plot(mu,sub,'o');
        h2.Color = h1.Color;
    end

box off
title('Distribution of beta weights');
xlabel('\beta');
ylabel('Subject');
hold off

p = get(gcf,'Position');
set(gcf,'Position',[p(1:2),1024,2048],'WindowState','Maximized');
drawnow, snapnow;

clear sub

% BEHAVIORAL OUTCOME

% over subjects

b1=figure;
hold off;
b1=histogram(fmri_dat.Y);
box off
title(['Histogram of single trial ' behav_outcome]);
xlabel(behav_outcome);
ylabel('n(observations)');
set(gcf,'WindowState','Maximized');
drawnow, snapnow;

% per subject

b2=figure;

    for sub = 1:n_subj
        this_idx_Y = find(sub == subject_id);
        this_Y = fmri_dat.Y(this_idx_Y);

        subplot(ceil(sqrt(n_subj)), ceil(n_subj/ceil(sqrt(n_subj))), sub);
        hold off
        b2 = histogram(this_Y);
        box off
        title(uniq_subject_id{sub});
        xlabel(behav_outcome);
        ylabel('n(obs)');
    end

set(gcf,'WindowState','Maximized');
drawnow, snapnow;

clear sub


%% CROSS-VALIDATION FOLD SELECTION
%--------------------------------------------------------------------------

% NOTE: balancing over groups, stratifying over subjects (i.e. leave whole
% subject out)

switch holdout_set_method_mvpa_reg_st

    case 'group'
        
        group = fmri_dat.metadata_table.(group_identifier);
        cv = cvpartition2(group, 'Group',subject_id, 'GroupKFold', nfolds_mvpa_reg_st);
            fold_labels = zeros(size(fmri_dat.dat,2),1);
            for sub = 1:cv.NumTestSets
                fold_labels(cv.test(sub)) = sub;
            end
        
    case 'onesample'
        
        cv = cvpartition2(size(fmri_dat.dat,2),'Group',subject_id, 'GroupKFold', nfolds_mvpa_reg_st);
            fold_labels = zeros(size(fmri_dat.dat,2),1);
            for sub = 1:cv.NumTestSets
                fold_labels(cv.test(sub)) = sub;
            end
            
    otherwise
        
        error('\ninvalid option "%s" defined in holdout_set_method_mvpa_reg_st variable, choose between "group" and "onesample"\n',holdout_set_method_mvpa_reg_st);
        
end
    
    
%% FIT SINGLE-LEVEL MVPA MODELS
%--------------------------------------------------------------------------

% NOTE: algorithm choice set in a2_default_options.m
% "help predict" in Matlab command window for more info

% RUN PREDICTIVE REGRESSION MODEL
%--------------------------------
switch ml_method_mvpa_reg_st
    
    case 'predict'
        
        t0 = tic;

        [cverr, stats, optout] = predict(fmri_dat, 'algorithm_name', algorithm_mvpa_reg_st, ...
                    'nfolds', fold_labels, 'error_type', 'mse', parallelstr_mvpa_reg_st, 'verbose', 0);

            % obtain bootstrapped weights if requested

            if dobootstrap_mvpa_reg_st

                [~ , bs_stats] = predict(fmri_dat, 'algorithm_name', algorithm_mvpa_reg_st,...
                    'bootsamples', boot_n_mvpa_reg_st, 'nfolds', 1 , 'error_type', 'mse', ...
                    parallelstr_mvpa_reg_st, 'verbose', 0);
                
            end % if loop bootstrap

            % obtain permutation weights if requested

            % NOTE: code from @pkragel s1_predict_anxiety_ratings_lassopcr.m -
            % NEEDS TESTING AND CHECK BY PHIL

            if doperm_mvpa_reg_st

                null_beta = zeros(perm_n_mvpa_reg_st,size(fmri_dat.dat,1)); % number of permutations, number of voxels

                for it = 1:perm_n_mvpa_reg_st

                    random_inds=randperm(size(fmri_dat.Y,1)); % number of images/ratings over subjects
                    temp_dat=fmri_dat;
                    temp_dat.Y=temp_dat.Y(random_inds);

                    [~, stats_null] = predict(temp_dat, 'algorithm_name', 'cv_lassopcr', 'nfolds', fold_labels, 'error_type', 'mse', parallelstr_mvpa_reg_st, 'verbose', 0);

                    for k = 1:max(fold_labels) % number of cv folds
                        regress_data = fmri_dat;
                        regress_data.X = stats_null.yfit(fold_labels==k);
                        regress_data.dat = regress_data.dat(:,fold_labels==k);
                        regress_stats(k) = regress(regress_data,'nodisplay');
%                         tv = replace_empty(regress_stats(k).b);
                        betas(k,:) = regress_stats(k).dat(:,1);
                    end

                    null_beta(it,:) = mean(betas);
                    null_weights(it,:) = stats_null.weight_obj.dat;

                    it/perm_n_mvpa_reg_st

                end

                clear temp_dat k regress_data regress_stats tv;

                for k = 1:max(fold_labels) % number of cv folds
                        regress_data = fmri_data;
                        regress_data.X = stats.yfit(fold_labels==k);
                        regress_data.dat = regress_data.dat(:,kinds=k);
                        regress_stats(k) = regress(regress_data,'nodisplay');
                        tv = replace_empty(regress_stats(k).b);
                        betas(k,:) = tv.dat(:,1);
                end

                mean_beta = mean(betas);

                sidedness = 'both';

                for i = 1:size(null_beta,2)
                    phat(i) = 1-sum(mean_beta(i) > null_beta(:,i)) / (1+size(null_beta,1));

                    % getting probability of finding observed difference from
                    % random permutations
                    switch sidedness
                        case 'both'
                            phat(i) = (length(find(abs(null_beta(:,i)) > abs(mean_beta(i))))+1) / (it+1); 
                        case 'smaller'
                            phat(i) = (length(find(null_beta(:,i) < mean_beta(i)))+1) / (it+1);
                        case 'larger'
                            phat(i) = (length(find(null_beta(:,i) > mean_beta(i)))+1) / (it+1);
                    end

                end

                perm_stats = stats.weight_obj;
                perm_stats.dat = mean_beta';
                perm_stats.p = phat';

            end % if loop permutation

        d_t = toc(t0);


        % PLOT OBSERVED VERSUS PREDICTED
        %-------------------------------

        fprintf('PCR r = %0.3f\n', corr(stats.yfit, fmri_dat.Y));

        figure
        
        line_plot_multisubject(fmri_dat.Y, stats.yfit, 'subjid', subject_id);
        xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['SVR Estimated ' behav_outcome],'(cross validated)'})
        
        set(gcf,'WindowState','Maximized');
        drawnow, snapnow;


        % PLOT MONTAGE OF UNTHRESHOLDED WEIGHTS
        %--------------------------------------

        fprintf ('\nShowing unthresholded %s results, %s\nScaling: %s\n\n', algorithm_mvpa_reg_st, mask_string, myscaling_mvpa_reg_st);

        figure

        o2 = canlab_results_fmridisplay([], 'compact', 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

        w = region(stats.weight_obj);

        o2 = addblobs(o2, w);
        o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' unthresholded ' mask_string]);

        figtitle = sprintf('%s_unthresholded_montage_%s_%s', algorithm_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;

        clear w, clear o2, clear figtitle


        % PLOT MONTAGE OF THRESHOLDED WEIGHTS AFTER BOOTSTRAPPING OR PERMUTATION TESTING
        %-------------------------------------------------------------------------------

        if dobootstrap_mvpa_reg_st

            fprintf ('\nShowing bootstrapped %s results, %s\nScaling: %s\n\n', algorithm_mvpa_reg_st, mask_string, myscaling_mvpa_reg_st);

            figure

            o2 = canlab_results_fmridisplay([], 'compact', 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

            t = bs_stats.weight_obj;
            t = threshold(t, q_threshold_mvpa_reg_st, 'fdr', 'k', k_threshold_mvpa_reg_st); 
            r = region(bs_stats.weight_obj);

            o2 = addblobs(o2, r);
            o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' bootstrapped ' mask_string]);

            figtitle = sprintf('%s_%1.4f_FDR_montage_%s_%s', algorithm_mvpa_reg_st, q_threshold_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
            set(gcf, 'Tag', figtitle, 'WindowState','maximized');
            drawnow, snapnow;

            clear w, clear o2, clear figtitle
            
        end % if loop bootstrap

        if doperm_mvpa_reg_st

            fprintf ('\nShowing permutation %s results, %s\nScaling: %s\n\n', algorithm_mvpa_reg_st, mask_string, myscaling_mvpa_reg_st);

            figure

            o2 = canlab_results_fmridisplay([], 'compact', 'outline', 'linewidth', 0.5, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]}, 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');

            t = bs_stats.weight_obj;
            t = threshold(t, q_threshold_mvpa_reg_st, 'fdr', 'k', k_threshold_mvpa_reg_st); 
            r = region(bs_stats.weight_obj);

            o2 = addblobs(o2, r);
            o2 = title_montage(o2, whmontage, [algorithm_mvpa_reg_st ' permutation ' mask_string]);

            figtitle = sprintf('%s_%1.4f_FDR_montage_%s_%s', algorithm_mvpa_reg_st, q_threshold_mvpa_reg_st, myscaling_mvpa_reg_st, mask_string);
            set(gcf, 'Tag', figtitle, 'WindowState','maximized');
            drawnow, snapnow;

            clear w, clear o2, clear figtitle

        end % if loop permutation
    
        
    case 'oofmridataobj'
        
        switch opt_method_mvpa_reg_st
            
            case 'gridsearch'
        
                % DEFINE ALGORITHM AND FEATURE EXTRACTOR
                %---------------------------------------
                alg = plsRegressor('numcomponents',5); % intiate alg as a plsRegressor estimator object, other estimators in Github repo/estimators; numcomponents is arbitrary, but typically low for pls - we will optimize this hyperparm later
                alg.fit(fmri_dat.dat', fmri_dat.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
                extractVxl = fmri2VxlFeatTransformer; % initiate extractVxl as an (empty?) fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
                extractVxl.fit(fmri_dat); % transformer takes fmri_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)

                % DEFINE PIPELINE
                %----------------
                fmri_pls = pipeline({{'featExt',extractVxl},{'pls',alg}}); % define fmri_pcr as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
                fmri_pls.fit(fmri_dat,fmri_dat.Y) % fit pipeline - JUST AS A TEST? NOT NEEDED FOR THE BELOW?

                % INNER CROSS-VALIDATION FUNCTION
                %--------------------------------
                switch holdout_set_method_mvpa_reg_st

                    case 'group'
                        innercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id);

                    case 'onesample'
                        innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

                end

                innercv(fmri_dat,fmri_dat.Y); % get cross-validation folds - JUST AS A TEST?

                % DEFINE OPTIMIZATION GRID
                %-------------------------
                gridPoints = table([2 4 5 9 20]','VariableNames',{'pls__numcomponents'}); % note the double underscore linking name of algorithm in pipeline with its parameter we want to optimize

                go = gridSearchCV(fmri_pls,gridPoints,innercv,@get_mse,'verbose',true); % see help gridSearchCV for inputs; scorers in Github repo/scorer; other optimizers in Github repo/estimators
                go.fit(fmri_dat,fmri_dat.Y);

                mdl = go.estimator.transformers{1}.brainModel; % empty .dat at this stage
                mdl.dat = go.estimator.estimator.B(:); % fills mdl.dat with betas
                figure
                mdl.montage;
                set(gcf, 'WindowState','maximized');
                drawnow, snapnow;

                % OUTER CROSS-VALIDATION FUNCTION
                %--------------------------------
                switch holdout_set_method_mvpa_reg_st

                    case 'group'
                        outercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id);

                    case 'onesample'
                        outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

                end
                
                % ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE
                %-------------------------------------------
                cvGS = crossValScore(fmri_pls, outercv, @get_mse, 'verbose', true); % see help/methods crossValScore; more crossvalidators in Github_repo/crossValidators - check out crossValPredict to get cross-validated prediction estimates at this stage
                cvGS.do(fmri_dat, fmri_dat.Y); % HOW DOES THIS TAKE GO INTO ACCOUNT SINCE IT IS NOT PART OF THAT PIPELINE? NOT NESTED? SOMETHING WRONG WITH ORDER OF THINGS?
                cvGS.do_null(); % fits null model - intercept only
                f1 = cvGS.plot; % plots predicted versus observed
                % average MSE over folds = cv model performance
                
            case 'bayes'
                
                % DEFINE ALGORITHM
                %-----------------
                alg = plsRegressor();
%                 alg.fit(fmri_dat.dat', fmri_dat.Y);
                
                % INNER CROSS-VALIDATION FUNCTION
                %--------------------------------
                switch holdout_set_method_mvpa_reg_st

                    case 'group'
                        innercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata.(subj_identifier)); % WHY NOT METADATA_TABLE HERE? DEPENDS ON INPUT (not features here), SO MAYBE BECAUSE OF METADATA CONSTRUCTOR?

                    case 'onesample'
                        innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata.(subj_identifier)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

                end
                
                % DEFINE BAYESIAN OPTIMIZATION AND CROSS-VALIDATE
                %------------------------------------------------
                dims = optimizableVariable('numcomponents',[1,30],'Type','integer','Transform','log');
                bayesOptParams = {dims, 'AcquisitionFunctionName','expected-improvement-plus',...
                    'MaxObjectiveEvaluations',30, 'UseParallel', true, 'verbose',1, 'PlotFcn', {}};
                bo = bayesOptCV(alg,innercv,@get_mse,bayesOptParams);
                
                % DEFINE OTHER COMPONENTS OF PIPELINE
                %------------------------------------
                featConstructor_han = @(X) table(X.metadata_table.(subj_identifier),'VariableNames',{'participant_id'}); % WHY DO WE NEED THIS?
                extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han);
                zTransVxl = zscoreVxlTransformer(@(X) X.metadata_table.participant_id);
                zTransImg = functionTransformer(@(x1) rescale(x1,'zscoreimages'));
                
                % DEFINE PIPELINE
                %------------------------------------
%                 fmri_pls = pipeline({{'zscorevxl', zTransVxl},{'zscoreimg', zTransImg},...
%                     {'featExt',extractVxl},{'bo_pls',bo}});
                fmri_pls = pipeline({{'featExt',extractVxl},{'bo_pls',bo}});
%                 data = features(fmri_dat.dat', fmri_dat.metadata_table.subject_id); % this is an "extended double" that is just a double with metadata in the dat.metadata field
                fmri_pls.fit(fmri_dat,fmri_dat.Y);
                
                % OUTER CROSS-VALIDATION FUNCTION
                %--------------------------------
                switch holdout_set_method_mvpa_reg_st

                    case 'group'
                        outercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id);

                    case 'onesample'
                        outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

                end
                
                % ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE
                %-------------------------------------------
                cvGS = crossValScore(fmri_pls, outercv, @get_mse, 'verbose', true);
                cvGS.do(fmri_dat, fmri_dat.Y);
                cvGS.do_null();
                f1 = cvGS.plot();
                
                % PLOT MODEL
                %-----------
                mdl = fmri_pls.estimator.transformers{1}.brainModel; % empty .dat at this stage
                mdl.dat = fmri_pls.estimator.estimator.B(:); % fills mdl.dat with betas
                mdl.montage;
                
                % FURTHER REPORTING OF RESULTS?
                % boostrapping and thresholding?
                
            otherwise
                
                error('\ninvalid option "%s" defined in opt_method_mvpa_reg_st variable, choose between "gridsearch" and "bayes"\n',opt_method_mvpa_reg_st);

        end % switch optimization method for oofmridataobj
        
    otherwise
        
        error('\ninvalid option "%s" defined in ml_method_mvpa_reg_st variable, choose between "oofmridataobj" and "predict"\n',ml_method_mvpa_reg_st);

end % switch machine learning method


%% SAVE STATS FOR SINGLE LEVEL MODELS
%--------------------------------------------------------------------------

if dosavemvparegstats

    savefilename = fullfile(resultsdir, 'single_trial_multilevel_MVPA_results.mat');
    save(savefilename, 'stats', '-v7.3');

        if dobootstrap_mvpa_reg_st
            save(savefilename,'bs_stats', '-v7.3', '-append');
        end

        if doperm_mvpa_reg_st
            save(savefilename,'perm_stats', '-v7.3', '-append');
        end
        
end
        

%% FIT MULTILEVEL MVPA MODEL
%--------------------------------------------------------------------------

if domultilevel_mvpa_reg_st

    % DETERMINE MAXIMUM NUMBER OF COMPONENTS
    %---------------------------------------

    max_comp = floor(size(fmri_dat.dat,2).*0.75 - n_subj);

    % NOTE: we specify the maximum number of components as < the number of columns in
    % fmri_dat.dat (n_subjects*n_conditions(in every fold)) to avoid overfitting in multilevel models, 
    % where we need to leave df for the random intercepts (upper bound 1 df per random intercept hence subject)

    % FIT SINGLE LEVEL MODEL
    %-----------------------

    [sl_cverr, sl_stats,sl_optout] = fmri_dat.predict('algorithm_name','cv_pcr',...
        'nfolds',fold_labels, 'numcomponents',max_comp);
    fprintf('PCR r = %0.3f\n', corr(sl_stats.yfit,fmri_dat.Y));

    figure

    line_plot_multisubject(fmri_dat.Y, sl_stats.yfit, 'subjid', subject_id);
    xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['PCR Estimated ' behav_outcome],'(cross validated)'})

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    figure

    sl_stats.weight_obj.montage;

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;


    % FIT MULTILEVEL MODEL W/FIXED PARAMETER ESTIMATION
    %--------------------------------------------------

    % split maximum amount of components in between and within

    n_bt_comp = floor(0.75*n_subj);
    n_wi_comp = max_comp - n_bt_comp;

    % NOTE: max between = n_subj IN EVERY FOLD (hence n_subj - 20% in 5-fold CV), 
    % and you want to put more money on within since this typically explains
    % more variance

    % overall model prediction

    [ml_cverr, ml_stats, ml_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr',...
        'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id);
    fprintf('multilevel PCR r = %0.3f\n',corr(ml_stats.yfit, fmri_dat.Y));

    % NOTE: algorithm option 'cv_mlpcr' requires subject identifier,
    % which makes sense since this is a multilevel/mixed model
    % note that fold labels are the same, since they respect subject membership

    figure

    line_plot_multisubject(fmri_dat.Y, ml_stats.yfit, 'subjid', subject_id);
    xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['MLPCR Estimated ' behav_outcome],'(cross validated)'})

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    figure

    ml_stats.weight_obj.montage;

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    figure

    subplot(1,2,1)
    line_plot_multisubject(sl_stats.yfit, ml_stats.yfit, 'subjid', subject_id);
    xlabel({'PCR model prediction'}); ylabel('Multilevel PCR model prediction');
    axis square
    subplot(1,2,2);
    plot(sl_optout{1}(:),ml_optout{1}(:),'.');
    lsline;
    xlabel('PCR model weights'); ylabel('Multilevel PCR model weights');
    axis square

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    % NOTE: contrary to @bogpetre's walkthrough, the
    % pcr and the multilevel pcr models are not exactly equivalent anymore
    % since I have been specifying the number of components

    % get the variance explained by the between and within component

    % NOTE: These functions call the same thing under the hood, 
    % but simply perform cross validation using ONLY between or within
    % subject models.

    [ml_bt_cverr, ml_bt_stats] = fmri_dat.predict('algorithm_name','cv_mlpcr_bt',...
        'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
    pred_bt = ml_bt_stats.yfit;

    [ml_wi_cverr, ml_wi_stats] = fmri_dat.predict('algorithm_name','cv_mlpcr_wi',...
        'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
    pred_wi = ml_wi_stats.yfit;

    % NOTE: algorithm options are created by bogpetre

    fprintf('Between subject PCR components r = %0.3f\n', corr(ml_bt_stats.yfit, fmri_dat.Y));
    fprintf('Within subject PCR components r = %0.3f\n', corr(ml_wi_stats.yfit, fmri_dat.Y));

    figure

    subplot(1,2,1)
    line_plot_multisubject(fmri_dat.Y, pred_bt, 'subjid', subject_id);
    xlabel({['Observed ' behav_outcome]}); ylabel('Between subject components'' prediction');
    axis square
    subplot(1,2,2)
    line_plot_multisubject(fmri_dat.Y, pred_wi, 'subjid', subject_id);
    xlabel({['Observed ' behav_outcome]}); ylabel('Within subject components'' prediction');
    axis square

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;


    % FIT MULTiLEVEL MODEL W/ MIXED PARAMETER ESTIMATION
    %---------------------------------------------------

    % NOTE: this is not part of the walkthrough yet, but @bogpetre pushed
    % mlpcr3 function to CanlabCore

    % NOTE bogpetre:
    % main function is mlpcr3, so help mlpcr3 for usage options.
    % basically it's the same as cv_mlpcr, except there's a randInt, randSlope and fitlmeOpts now
    % fitlmeOpts get passed on to fitlme. It picks some sensible defaults
    % randSlope makes things much slower. randInt is roughly the sae order of magnitude as running the fixed effects version
    % the function defaults to fixed effects by default, so it's a drop in replacement
    % for mlpcr2.m (aka cv_mlpcr)

    % overall model prediction including random intercept only

    [ml3_cverr, ml3_stats, ml3_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr3',...
        'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1);

    fprintf('multilevel PCR r = %0.3f\n',corr(ml3_stats.yfit, fmri_dat.Y));

    % NOTE: compare with code in previous section and note change of algorithm_name option to cv_mlpcr3 and addition of randInt
    % option - see help mlpcr3 for more details

    figure

    line_plot_multisubject(fmri_dat.Y, ml3_stats.yfit, 'subjid', subject_id);
    xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['MLPCR3 Estimated ' behav_outcome],'(cross validated)'})

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    figure

    ml3_stats.weight_obj.montage;

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    figure

    subplot(1,2,1)
    line_plot_multisubject(ml_stats.yfit, ml3_stats.yfit, 'subjid', subject_id);
    xlabel({'Multilevel PCR model prediction'}); ylabel('Multilevel PCR model random int model prediction');
    axis square
    subplot(1,2,2);
    plot(ml_optout{1}(:),ml3_optout{1}(:),'.');
    lsline;
    xlabel('Multilevel PCR model weights'); ylabel('Multilevel PCR model random int model weights');
    axis square

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    % NOTE: models are very similar but not exactly equivalent

    % overall model prediction including random intercept and random slope

    [ml3rs_cverr, ml3rs_stats, ml3rs_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr3',...
        'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1, 'randSlope', 1);

    fprintf('multilevel PCR r = %0.3f\n',corr(ml3rs_stats.yfit, fmri_dat.Y));

    figure

    line_plot_multisubject(fmri_dat.Y, ml3rs_stats.yfit, 'subjid', subject_id);
    xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['MLPCR3 Estimated ' behav_outcome],'(cross validated)'})

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;

    figure

    ml3rs_stats.weight_obj.montage;

    set(gcf, 'WindowState','maximized');
    drawnow, snapnow;


    %% SAVE STATS FOR MULTILEVEL MODELS
    %--------------------------------------------------------------------------
    
    if dosavemvparegstats
    
        savefilename = fullfile(resultsdir, 'single_trial_multilevel_MVPA_results.mat');
        save(savefilename, 'pcr_stats','mlpcr_stats','mlpcr3_stats','mlpcr3rs_stats', '-v7.3');
        
    end
    
end % if loop multilevel

