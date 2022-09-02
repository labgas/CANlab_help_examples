%% c2g_run_multivariate_mediation_single_trial.m
%
%
% USAGE
% -----
%
% This script runs multivariate mediation analysis on a continuous outcome Y
% on an fmri_data_st object created using prep_3f_create_fmri_data_single_trial_object.m
% That script should be run first, or the present script will load the data
% object if it is saved by the previous script.
%
% Options for this script are set in a2_set_default_options.m, see that
% script and below for more info. Many of these options get passed into CANlab's
% multivariateMediation function.
%
% Run this script with Matlab's publish function to generate html report of
% results:
% publish('c2g_run_multivariate_mediation_single_trial','outputDir',htmlsavedir)
%
% TUTORIALS AND DOCUMENTATION
% ---------------------------
%
% Here is the PDM toolbox within CANlab's MediationToolbox on Github
% https://github.com/canlab/MediationToolbox/tree/master/PDM_toolbox
%
% It includes a helpful example scripts and a README referring to papers
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/Multivariate_Mediation_ExampleScript.m
% https://github.com/canlab/MediationToolbox/blob/master/PDM_toolbox/README.md
%
% Here is an older script performing single-level rather than multilevel
% mediation, among other types of mediation analyses
% https://github.com/labgas/proj-emosymp/blob/main/secondlevel/model_1_CANlab_classic_GLM/emosymp_m1_s6_mediation_NPS.m
% 
% OPTIONS
% -------
%
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options, 
% you can make a copy of this script with a letter index (e.g. _s6a_) 
% and change the default options below
%
% GENERAL OPTIONS
%
% zscore_outcome_pdm: default false; true zscores behavioral outcome variable (fmri_dat.Y) prior to fitting models
% maskname_pdm: default which('gray_matter_mask_sparse.img');
%       - default use of sparse gray matter mask
%       - maskdir now defined in a_set_up_paths_always_run_first script
%       - if you do not want to mask, change to []
%       - if you want to use a custom mask, put it in maskdir and change name here.
% myscaling_pdm: default 'raw'; options are 'raw', 'centerimages', 'zscoreimages', 'l2normimages', 'zscorevoxels'
%
% STATISTICS AND RESULTS VISUALIZATION OPTIONS
% --------------------------------------------
% nPDM = 10: default 10; number of PDMs to retain, chances are very low that meaningful variance is explained by PDM # > 10
% dobootstrap_pdm: default false; bootstrapping, does not take an awful lot of time in this case
%     boot_n_pdm: default 5000; number of bootstrap samples, reduce number for quick results, increase to 10k for publication
% dosourcerecon_pdm: default false; source reconstruction/"structure coefficients", i.e. regressing each voxel's activity onto yhat - see Haufe et al NeuroImage 2014
% dosavepdmstats: default true; saves all results as .mat files
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   August, 2022
%__________________________________________________________________________
% @(#)% c2g_run_multivariate_mediation_single_trial     v1.3        
% last modified: 2022/09/02


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% SET/COPY MANDATORY OPTIONS FROM CORRESPONDING PREP_3f_ SCRIPT

results_suffix = ''; % suffix of your choice added to .mat file with saved results
cons2exclude = {}; % cell array of condition names to exclude, separated by commas (or blanks); NOTE: you can expand this if you want to exclude more conditions at this stage
behav_outcome = 'rating'; % name of outcome variable in DAT.BEHAVIOR.behavioral_data_table_st
subj_identifier = 'participant_id'; % name of subject identifier variable in same table
cond_identifier = 'trial_type'; % name of condition identifier variable in same table
% group_identifier = 'group'; % name of group identifier variable in same table; leave commented out if you don't have groups

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run a second version of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m
% 
% See documentation above and a2_set_default_options.m for list of options


%% LOAD FMRI_DATA_ST OBJECT AND OTHER NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('LOADING DATA');
fprintf('\n\n');

if ~exist('DSGN','var') || ~exist('DAT','var')

    load(fullfile(resultsdir,'image_names_and_setup.mat'));

    if ~isfield(DAT,'BEHAVIOR')
        fprintf('\n');
        error('Behavioral data not yet added to DAT structure - run prep_1b script first')
    end

end

if ~exist('fmri_dat','var')
    
    if ~isempty(cons2exclude)
        load(fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome, '_exclude_cond_', char([cons2exclude{:}]), '_', results_suffix, '.mat']));
        
    else
        load(fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome, '_', results_suffix, '.mat']));

    end
end


%% DEFINE SUBJECT AND CONDITION IDENTIFIERS
% -------------------------------------------------------------------------

subject_id = fmri_dat.metadata_table.(subj_identifier);
[uniq_subject_id, ~, subject_id] = unique(subject_id,'stable');
fmri_dat.metadata_table.subject_id = subject_id;
n_subj = size(uniq_subject_id,1);

cond_id = fmri_dat.metadata_table.(cond_identifier);
[uniq_cond_id, ~, cond_id] = unique(cond_id,'stable');
fmri_dat.metadata_table.cond_id = cond_id;

if ~isempty(cons2exclude)

    for cont = 1:size(DAT.contrastnames,2)
        for con = 1:size(cons2exclude,2)
            cont2include(cont,con) = ~contains(DAT.contrastnames{1,cont},cons2exclude{con});
        end
        idx_cont2include = all(cont2include,2);
    end

    clear cont con

    contrastnames2include = DAT.contrastnames(idx_cont2include');
    contrastweights2include = DAT.contrasts(idx_cont2include,:);

else
    
    contrastnames2include = DAT.contrastnames;
    contrastweights2include = DAT.contrasts;
    
end

for cont = 1:size(contrastnames2include,2)
    contweights2include(cont,:) = abs(contrastweights2include(cont,:));
    idx_contweights2include = logical(contweights2include);
end

clear cont con


%% CREATE DIRECTORY STRUCTURE TO WRITE RESULTS
% -------------------------------------------------------------------------

mediationresultsdir = fullfile(resultsdir,'mediation_analysis');
pdmmediationresultsdir = fullfile(mediationresultsdir,'pdm');

if ~exist(mediationresultsdir,'dir')
    mkdir(mediationresultsdir);
    mkdir(pdmmediationresultsdir);
end

%% SCALE AND/OR MASK IMAGES AND BEHAVIORAL OUTCOME ACCORDING TO OPTIONS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('MASKING AND SCALING IMAGES IF REQUESTED IN OPTIONS');
fprintf('\n\n');

% MASKING IMAGES
%---------------

if exist('maskname_pdm','var') && ~isempty(maskname_pdm) && exist(maskname_pdm, 'file')

    [~,maskname_short] = fileparts(maskname_pdm);
    fprintf('\nMasking data with %s\n\n',maskname_short);
    mask_string = sprintf('masked with %s', maskname_short);

    pdmmask = fmri_mask_image(maskname_pdm);
    fmri_dat = fmri_dat.apply_mask(pdmmask);
    fmri_dat.mask_descrip = maskname_pdm;

else

    fprintf('\nNo mask found; using full original image data\n\n');
    mask_string = sprintf('unmasked');
    

end % if loop mask


% SCALING IMAGES
%---------------

switch myscaling_pdm

    case 'raw'

        fprintf('\nNo scaling of input images\n\n');

    case 'centerimages'

        fmri_dat = fmri_dat.rescale('centerimages'); 
        fprintf('\nCentering input images\n\n');

    case 'l2normimages'

        fmri_dat = fmri_dat.rescale('l2norm_images');
        fprintf('\nNormalizing input images by l2norm\n\n');

    case 'zscoreimages'

        fmri_dat = fmri_dat.rescale('zscoreimages');
        fprintf('\nZ-scoring input images\n\n');

    case 'zscorevoxels'

        fmri_dat = fmri_dat.rescale('zscorevoxels');
        fprintf('\nZ-scoring voxels across input images\n\n');

    otherwise 

        error('\ninvalid scaling option %s specified in myscaling_pdm variable defined in a2_set_default_options script, please correct\n\n', myscaling_pdm)

end % switch scaling
   

% ZSCORE BEHAVIORAL OUTCOME
%--------------------------

% NOTE: useful for more interpretable values of prediction MSE

if zscore_outcome_pdm

    fmri_dat.Y = zscore(fmri_dat.Y);
    fprintf('\nZ-scoring outcome fmri_dat.Y across subjects\n\n');

end


%% DATA VISUALISATION PRIOR TO MODEL BUILDING
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('PLOTTING DATA');
fprintf('\n\n');

% BETA IMAGES
%------------

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
%-------------------

% OVER SUBJECTS

b1=figure;
hold off;
b1=histogram(fmri_dat.Y);
box off
title(['Histogram of single trial ' behav_outcome]);
xlabel(behav_outcome);
ylabel('n(observations)');
set(gcf,'WindowState','Maximized');
drawnow, snapnow;

% PER SUBJECT

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

        
%% SET UP AND RUN MULTIVARIATE MEDIATION
% -------------------------------------------------------------------------

for cont = 1:size(contrastnames2include,2)
    
    fprintf('\n\n');
    printhdr(['CONTRAST #', num2str(cont), ': ',upper(contrastnames2include{cont})]);
    fprintf('\n\n');
    
    % SET UP
    %-------
    
    condsincontrast = DAT.conditions(1,idx_contweights2include(cont,:));
    
    contrastpdmmediationresultsdir = fullfile(pdmmediationresultsdir,contrastnames2include{1,cont});
    
    if ~exist(contrastpdmmediationresultsdir,'dir')
        mkdir(contrastpdmmediationresultsdir);
    end
    
    behavcontrastpdmmediationresultsdir = fullfile(contrastpdmmediationresultsdir,behav_outcome);
    
    if ~exist(behavcontrastpdmmediationresultsdir,'dir')
        mkdir(behavcontrastpdmmediationresultsdir);
    end
    
    if isempty(dir(fullfile(behavcontrastpdmmediationresultsdir,'PDM_*'))) || isempty(dir(fullfile(behavcontrastpdmmediationresultsdir,'data_objects_*')))
        
        fprintf('\n\n');
        printhdr('Creating data objects');
        fprintf('\n\n');
        
        for sub = 1:n_subj
            
            idx_sub = fmri_dat.metadata_table.subject_id == sub;
            
            fmri_dat_sub = get_wh_image(fmri_dat,idx_sub);
            
            for con = 1:size(condsincontrast,2)
                idx_con = contains(fmri_dat_sub.metadata_table.(cond_identifier),condsincontrast{1,con});
                fmri_dat_sub_con = get_wh_image(fmri_dat_sub,idx_con);
                Y_temp{con} = fmri_dat_sub_con.metadata_table.(behav_outcome);
                X_temp{con} = fmri_dat_sub_con.metadata_table.cond_id;
                X_temp{con}(X_temp{con}==max(X_temp{con})) = 1;
                X_temp{con}(X_temp{con}==min(X_temp{con})) = -1;
                M_temp{con} = fmri_dat_sub_con.dat;
            end
            
            clear idx_con
            
            Y{sub} = cat(1,Y_temp{:});
            X{sub} = cat(1,X_temp{:});
            M{sub} = cat(2,M_temp{:});
            
        end % for loop subjects
        
        if dosavepdmstats
            
            fprintf('\n\n');
            printhdr('Saving data objects');
            fprintf('\n\n');
        
            if exist('maskname_short', 'var')
                savefilename_data = fullfile(behavcontrastpdmmediationresultsdir, ['data_objects_', myscaling_pdm, '_', maskname_short, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

            else
                savefilename_data = fullfile(behavcontrastpdmmediationresultsdir, ['data_objects_', myscaling_pdm, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

            end
                
            save(savefilename_data,'Y','X','M','-v7.3');
            
        end % if dosavepdmstats
        
    % RUN MULTIVARIATE MEDIATION
    % --------------------------
    
        fprintf('\n\n');
        printhdr('Calculating full PDM path coefficients');
        fprintf('\n\n');

        pdm = multivariateMediation(X,Y,M,'nPDM',nPDM,'plots'); 
        set(gcf,'WindowState','maximized');
        drawnow,snapnow;
        
        pdmfull = pdm;
        
        pathcoeffs = [pdmfull.Theta{:}];
        absabpaths = abs(pathcoeffs(5,:));
        idx_localmax = islocalmax(diff(absabpaths));
        first_localmax = find(idx_localmax);
        first_localmax = first_localmax(1,1);
        nPDM2retain = first_localmax;
        
        fprintf('\n\n');
        printhdr('Path coefficients for PDMs to retain');
        fprintf('\n\n');
        
        pdm = multivariateMediation(pdm,'nPDM',nPDM2retain);
        
        if dosavepdmstats
            
            fprintf('\n\n');
            printhdr('Saving PDM results');
            fprintf('\n\n');
        
            if exist('maskname_short', 'var')
                savefilename_pdm = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_results_', myscaling_pdm, '_', maskname_short, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);
                
            else
                savefilename_pdm = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_results_', myscaling_pdm, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

            end
            
            if dobootstrap_pdm
                
                fprintf('\n\n');
                printhdr('Bootstrapping path coefficients for PDMs to retain');
                fprintf('\n\n');
                
                pdm = multivariateMediation(pdm, 'noPDMestimation','bootPDM',1:nPDM2retain,'bootjPDM','Bsamp',boot_n_pdm,'save2file',savefilename_pdm); 
                
            else
                pdm = multivariateMediation(pdm, 'noPDMestimation','save2file',savefilename_pdm); 
                    
            end
                
            save(savefilename_pdm,'pdmfull','pdm','-v7.3');
            
        else
            
            if dobootstrap_pdm
                
                fprintf('\n\n');
                printhdr('Bootstrapping path coefficients for PDMs to retain');
                fprintf('\n\n');
                
                pdm = multivariateMediation(pdm,'noPDMestimation','bootPDM',1:nPDM2retain,'bootjPDM','Bsamp',boot_n_pdm); 
                
            end
            
        end % if dosavepdmstats
        
    else
        
        fprintf('\n\n');
        printhdr('Loading data objects & PDM results');
        fprintf('\n\n');
        
        if exist('maskname_short', 'var')
            loadfilename_data = fullfile(behavcontrastpdmmediationresultsdir, ['data_objects_', myscaling_pdm, '_', maskname_short, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);
            loadfilename_pdm = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_results_', myscaling_pdm, '_', maskname_short, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

        else
            loadfilename_data = fullfile(behavcontrastpdmmediationresultsdir, ['data_objects_', myscaling_pdm, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);
            loadfilename_pdm = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_results_', myscaling_pdm, '_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);
                            
        end
        
        load(loadfilename_data);
        load(loadfilename_pdm);
        
%         pdm = out;
%         clear out;
          
    end % if results files not yet available
    
    
    % PLOT MONTAGE OF UNTHRESHOLDED PDM RESULTS
    % ------------------------------------------
    
    fprintf('\n\n');
    printhdr('Visualizing unthresholded PDM results');
    fprintf('\n\n');
    
    whmontage = 5;
    
    for comp = 1:size(pdm.Wfull,2)

        fprintf ('\nSHOWING UNTHRESHOLDED PDM RESULTS, CONTRAST: %s, PDM #%s, %s, SCALING: %s\n\n', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);

        dat_unt{comp} = fmri_dat;
        dat_unt{comp}.dat = pdm.Wfull{comp};
        
        reg_unt{comp} = region(dat_unt{comp});

        o2 = montage(reg_unt{comp}, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [contrastnames2include{cont} ' unthresholded PDM #' num2str(comp), ' ', mask_string, ' ', myscaling_pdm]);

        figtitle = sprintf('%s_unthresholded_montage_PDM#%s_%s_%s', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;
            if save_figures_pdm
                plugin_save_figure;
            end
        clear o2, clear figtitle
        
    end
      
    % SOURCE RECONSTRUCTION IF REQUESTED
    % ----------------------------------
    
    if dosourcerecon_pdm
        
        % NOTE: code from emosymp_m1_s6_mediation_NPS.m, checked by Bogdan & Martin
        % "source reconstruction is thus no more than a covariance map which
        % shows how much each voxel covaries with model prediction across images" - Bogdan
        % ref Haufe et al, NeuroImage 2014
        
        fprintf('\n\n');
        printhdr('Calculating source reconstruction maps');
        fprintf('\n\n');
        
        if ~exist(fullfile(behavcontrastpdmmediationresultsdir, ['PDM_source_recon_', behav_outcome, '_', contrastnames2include{cont}, '_', maskname_short, '_', results_suffix, '.mat']),'file') || ...
                ~exist(fullfile(behavcontrastpdmmediationresultsdir, ['PDM_source_recon_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']),'file')
    
            for comp = 1:size(pdm.Wfull,2)

                fmri_dat_cont = fmri_dat;

                for con = 1:size(condsincontrast,2)
                        idx_con(:,con) = contains(fmri_dat_cont.metadata_table.(cond_identifier),condsincontrast{1,con});
                end
                
                idx_con = any(idx_con,2);
                fmri_dat_cont = get_wh_image(fmri_dat_cont,idx_con); % may need to be transposed

                X_source{cont} = fmri_dat_cont;
                Xz_source{cont} = rescale(X_source{cont},'centervoxels');
                M_source{cont,comp} = pdm.Wfull{1,comp};
                P_source{cont,comp} = M_source{cont,comp}'*Xz_source{cont}.dat;
                source{cont,comp} = (P_source{cont,comp}*Xz_source{cont}.dat');
                source{cont,comp} = source{cont,comp}'./(size(P_source{cont,comp},2)-1);
                source_obj{cont,comp} = fmri_dat;
                source_obj{cont,comp}.dat = source{cont,comp};
                source_obj{cont,comp}.image_names = 'source reconstruction images';
                source_obj{cont,comp}.fullpath = '';
                source_obj{cont,comp}.history = {['source reconstructed from single trial beta images and pdm_',num2str(comp)]};
                source_obj_cont{comp} = source_obj{cont,comp};

            end % for loop over pdm components

            if dosavepdmstats
                
                fprintf('\n\n');
                printhdr('Saving source reconstruction results');
                fprintf('\n\n');

                if exist('maskname_short', 'var')
                    savefilename_source_recon = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_source_recon_', behav_outcome, '_', contrastnames2include{cont}, '_', maskname_short, '_', results_suffix, '.mat']);

                else
                    savefilename_source_recon = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_source_recon_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

                end

                save(savefilename_source_recon,'source_obj_cont','-v7.3');

            end % if loop savepdmstats
            
        else
            
            if exist('maskname_short', 'var')
                loadfilename_source_recon = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_source_recon_', behav_outcome, '_', contrastnames2include{cont}, '_', maskname_short, '_', results_suffix, '.mat']);

            else
                loadfilename_source_recon = fullfile(behavcontrastpdmmediationresultsdir, ['PDM_source_recon_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

            end

            load(loadfilename_source_recon);
            
        end % if results files already exist
        
    end % if loop source reconstruction
    
    
    % PLOT MONTAGE OF UNTHRESHOLDED SOURCE RECONSTRUCTION RESULTS
    % -----------------------------------------------------------
    
    fprintf('\n\n');
    printhdr('Visualizing unthresholded source reconstruction results');
    fprintf('\n\n');
    
    for comp = 1:size(pdm.Wfull,2)

        fprintf ('\nSHOWING UNTHRESHOLDED SOURCE RECONSTRUCTION RESULTS, CONTRAST: %s, PDM #%s, %s, SCALING: %s\n\n', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
        
        reg_sc{comp} = region(source_obj_cont{comp});

        o2 = montage(reg_sc{comp}, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
        o2 = title_montage(o2, whmontage, [contrastnames2include{cont} ' source reconstruction unthresholded PDM #' num2str(comp), ' ', mask_string, ' ', myscaling_pdm]);

        figtitle = sprintf('%s_unthresholded_montage_source_reconstruction_PDM#%s_%s_%s', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
        set(gcf, 'Tag', figtitle, 'WindowState','maximized');
        drawnow, snapnow;
            if save_figures_pdm
                plugin_save_figure;
            end
        clear o2, clear figtitle
        
    end
    
    
    % PLOT MONTAGES OF THRESHOLDED RESULTS IF BOOTSTRAPPED
    % ----------------------------------------------------
    
    if dobootstrap_pdm
        
        fprintf('\n\n');
        printhdr('Visualizing thresholded PDM results after bootstrapping');
        fprintf('\n\n');
    
        % FDR-CORRECTED

        fprintf ('\nSHOWING FDR-CORRECTED PDM RESULTS, CONTRAST: %s, PDM #%s, %s, SCALING: %s\n\n', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);

        for comp = 1:size(pdm.Wfull,2)

            dat_fdr{comp} = dat_unt{comp};
            dat_fdr{comp}.dat = pdm.Wfull{comp}.*(pdm.boot.p{comp}<pdm.pThreshold(comp));
            reg_fdr{comp} = region(dat_fdr{comp});

            % Results table

            fprintf ('\nTable bootstrapped PDM results, contrast: %s, PDM #%s, %s, SCALING: %s\n\n', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
            fprintf('FDR q < .05 = p < %3.8f\n', pdm.pThreshold(comp));

            [reg_fdr{comp}, ~, ~ ] = autolabel_regions_using_atlas(reg_fdr{comp});
            [reg_pos_fdr{comp},reg_neg_fdr{comp},results_table_fdr{comp}] = table(reg_fdr{comp},'k',k_threshold_pdm);
            reg_all_fdr{comp} = [reg_pos_fdr{comp},reg_neg_fdr{comp}];

            fprintf('\n\n');

            % Montage

            fprintf ('\nMontage bootstrapped PDM results, contrast: %s, PDM #%s, %s, SCALING: %s\n\n', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
            fprintf('FDR q < .05 = p < %3.8f\n', pdm.pThreshold(comp));

            o2 = montage(reg_all_fdr{comp}, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
            o2 = title_montage(o2, whmontage, [contrastnames2include{cont} ' FDR-corrected PDM # ' num2str(comp), ' ', mask_string, ' ', myscaling_pdm]);

            figtitle = sprintf('%s_FDR_montage_PDM#%s_%s_%s', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
            set(gcf, 'Tag', figtitle, 'WindowState','maximized');
            drawnow, snapnow;
                if save_figures_pdm
                    plugin_save_figure;
                end
            clear o2, clear figtitle

            % Regioncenters

            fprintf ('\nRegioncenters bootstrapped PDM results, contrast: %s, PDM #%s, %s, SCALING: %s\n\n', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
            fprintf('FDR q < .05 = p < %3.8f\n', pdm.pThreshold(comp));

                % Positive weights

                o2 = montage(reg_pos_fdr{comp}, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'regioncenters');

                figtitle = sprintf('%s_FDR_montage_positive_PDM#%s_%s_%s', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_pdm
                        plugin_save_figure;
                    end
                clear o2, clear figtitle

                % Negative weights

                o2 = montage(reg_neg_fdr{comp}, 'colormap', 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]},'regioncenters');

                figtitle = sprintf('%s_FDR_montage_negative_PDM#%s_%s_%s', contrastnames2include{cont}, num2str(comp), mask_string, myscaling_pdm);
                set(gcf, 'Tag', figtitle, 'WindowState','maximized');
                drawnow, snapnow;
                    if save_figures_pdm
                        plugin_save_figure;
                    end
                clear o2, clear figtitle

        end % for loop components visualization FDR thresholded results
        
        % UNCORRECTED
        
        % not yet implemented, neither is visualization of joint PDM - both of less interest
        
        % SAVE REGIONS AND TABLES IF REQUESTED
        
        if dosavepdmstats
        
            if exist('maskname_short', 'var')
                savefilename_regions = fullfile(behavcontrastpdmmediationresultsdir, ['region_objects_tables_', behav_outcome, '_', contrastnames2include{cont}, '_', maskname_short, '_', results_suffix, '.mat']);

            else
                savefilename_regions = fullfile(behavcontrastpdmmediationresultsdir, ['region_objects_tables_', behav_outcome, '_', contrastnames2include{cont}, '_', results_suffix, '.mat']);

            end
            
            fprintf('\n\n');
            printhdr('Saving regions objects and tables');
            fprintf('\n\n');
            
            disp('dat_fdr, etc.  : cell array of thresholded data objects at FDR threshold');
            disp('dat_unc, etc.  : cell array of thresholded data objects at uncorrected threshold');
            disp('In cell arrays of data objects, dat_*{i}.dat contains thresholded data for the i-th pdm');
            disp('reg_fdr, etc.  : cell array of labeled region objects at FDR threshold');
            disp('reg_unc, etc.  : cell array labeled region objects at uncorrected threshold');
            disp('reg_all_fdr, etc.  : cell array of labeled region objects at FDR threshold and extent threshold');
            disp('reg_all_unc, etc.  : cell array labeled region objects at uncorrected threshold and extent threshold');
            disp('In cell arrays of region objects, region_obj{i}(j).dat contains extracted data for the i-th pdm and j-th significant region, averaged over voxels');
            disp('reg_all_x is used to plot montages of all regions');
            disp('Use these data in plots, secondary analyses, or to run mediation.m within individual regions');
            disp('reg_pos/neg_fdr, etc.  : cell array of region objects at FDR threshold, split for positive and negative, and thresholded at extent threshold');
            disp('reg_pos/neg_unc, etc.  : cell array of region objects at FDR threshold, split for positive and negative, and thresholded at extent threshold');
            disp('These data are used to plot montages of individual regions');
            disp('results_table_fdr, etc.  : cell array of tables at FDR threshold');
            disp('results_table_unc, etc.  : cell array of tables at uncorrected threshold');
            disp('Use these data to write tables to excel etc');
            fprintf('\n\n');
                
            save(savefilename_regions, 'dat_fdr', 'dat_unt', 'reg_*', 'results_*'); 
            
        end % if dosavepdmstats
        
    end % if loop bootstrap
    
    
end % for loop contrasts
            