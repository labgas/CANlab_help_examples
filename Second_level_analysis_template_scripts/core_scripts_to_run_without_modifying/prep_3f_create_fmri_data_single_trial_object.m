%% prep_3f_create_fmri_data_single_trial_object
%
%
% USAGE
%
% This script creates and saves an fmri_data_st object from the single trial 
% con images written in rootdir/firstlevel/model_x_yyy/sub-zz by
% LaBGAScore_firstlevel_s2_fit_model.m and 
% adds a convenient metadata_table field containing
% 1. single trial ratings from rootdir/BIDS/phenotype/<phenotype_trial>.tsv
% 2. single trial vifs
%
% OPTIONS
%
% NOTE: defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options (for example
% voxel- and parcelwise regression), you can make a copy of this script with
% a letter index (e.g. _s6a_) and change the default option here
%
% cons2exclude_dat_st: cell array of condition names to exclude, separated by commas (or blanks)
% behav_outcome_dat_st: name of outcome variable in DAT.BEHAVIOR.behavioral_data_table_st
% subj_identifier_dat_st: name of subject identifier variable in same table
% cond_identifier_dat_st: name of condition identifier variable in same table
% group_identifier_dat_st: name of group identifier variable in same table; leave commented out if you don't have groups
% vif_threshold_dat_st: variance inflation threshold to exclude trials
%
% MANDATORY OPTIONS TO BE SPECIFIED IN THIS SCRIPT
%
% results_suffix: name to add to results file to specify model in case of
% multiple models, e.g. 'water_excluded'
%
% IMPORTANT NOTE: this script has not been extensively tested on messy data
% yet (i.e. missing trials, NaNs for outcome, etc)!
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   Dartmouth, March, 2021
%__________________________________________________________________________
% @(#)% prep_3f_create_fmri_data_single_trial_object.m     v3.3       
% last modified: 2024/12/18


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% SET MANDATORY OPTIONS

results_suffix = ''; % adds a suffix of your choice to .mat file with results that will be saved
% NOTE: do NOT delete the latter option, leave empty if not needed
% NOTE: do NOT use to add a suffix specifying the behavioral outcome nor excluded conditions, this will be added automatically

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;
% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT!
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% cons2exclude_dat_st = {'varname1'}; % cell array of condition names to exclude, separated by commas (or blanks)
% behav_outcome_dat_st = 'varname2'; % name of outcome variable in DAT.BEHAVIOR.behavioral_data_table_st
% subj_identifier_dat_st = 'varname3'; % name of subject identifier variable in same table
% cond_identifier_dat_st = 'varname4'; % name of condition identifier variable in same table
% group_identifier_dat_st = 'varname5'; % name of group identifier variable in same table; leave commented out if you don't have groups
% vif_threshold_dat_st = x; % variance inflation threshold to exclude trials


%% LOAD NECESSARY VARIABLES IF NEEDED AND DO PREP WORK
% -------------------------------------------------------------------------
    
if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
    if ~isfield(DAT,'BEHAVIOR')
        error('\n Behavioral data not yet added to DAT structure - run prep_1b script first\n')
    end
    
end

[~,subjs] = fileparts(DSGN.subjects); % assuming use of all subjects for whom firstlevel was run in this particular model
% subjs = firstsubjs';                % use if subjects have been excluded at firstlevel stage; CAUTION: you will also need to exclude them from DAT.BEHAVIOR.behavioral_data_table to avoid mismatch, use lines below (UNTESTED)! 
% idx_subjs = contains(DAT.BEHAVIOR.behavioral_data_table.(subj_identifier_dat_st),subjs);
% DAT.BEHAVIOR.behavioral_data_table = DAT.BEHAVIOR.behavioral_data_table(idx_subjs,:);
% idx_subjs_st = contains(DAT.BEHAVIOR.behavioral_st_data_table.(subj_identifier_dat_st),subjs);

if ~isempty(cons2exclude_dat_st)

    for con2ex = 1:size(cons2exclude_dat_st,2)
        outcome_vars_between_idx(con2ex,:) = contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,behav_outcome_dat_st) & ~contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,cons2exclude_dat_st{con2ex});
    end
    
    for out = 1:size(outcome_vars_between_idx,2)
        outcome_vars_between_idx2(out) = logical(sum(outcome_vars_between_idx(:,out)));
    end
    
    outcome_vars_between = DAT.BEHAVIOR.behavioral_data_table(:,outcome_vars_between_idx2);

    for var = 1:sum(outcome_vars_between_idx2)
        missings{var} = isnan(table2array(outcome_vars_between(:,var))); % index for subjects with missing behavioral data, which we want to exclude
    end
    
    missings = cell2mat(missings);
    
else
    
    outcome_vars_between_idx = contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,behav_outcome_dat_st);
    outcome_vars_between = DAT.BEHAVIOR.behavioral_data_table(:,outcome_vars_between_idx);
    
    for var = 1:sum(outcome_vars_between_idx)
        missings{var} = isnan(table2array(outcome_vars_between(:,var))); % index for subjects with missing behavioral data, which we want to exclude
    end
    
    missings = cell2mat(missings);
    
end

for sub = size(missings,1)
    idx_behav(sub) = sum(missings(sub,:));
end
    
clear sub var out con2ex

idx_behav = ~idx_behav;
subjs2use = subjs(:,idx_behav); % exclude subjects if ratings for any of the conditions of interest are missing
idx_exclude = ~ismember(subjs,subjs2use);
subjs2use = subjs2use';
subjs2exclude = subjs(idx_exclude,:); 

if ~isempty(subjs2exclude)
    warning('subjects %s are missing behavioral data for at least one entire condition - excluding them from analysis',char(subjs2exclude))
else
    warning('no subjects with missing behavioral data for entire condition - including all subjects in analysis')
end

    
%% READ CON IMAGES AND VIFS AND ADD TO FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

conimgs = struct('conpath',{''},'conname',{''},'subjname',{''},'vifname',{''},'vifvalue',[]);

nr_cons = size(DSGN.contrasts,2); % number of contrasts estimated over single trials
nr_cons_no_st = sum(cell2mat(DSGN.singletrials{1}) == 0); % number of conditions for which no single trial analysis is performed

    for sub = 1:size(subjs2use,1)
        
        subjdir = fullfile(DSGN.modeldir,subjs2use(sub,:));
        load(fullfile(char(subjdir),'SPM.mat'));
        load(fullfile(char(subjdir),'diagnostics','vifs.mat'));
        subjconimgs = SPM.xCon;
        
        if ~isempty(cons2exclude_dat_st)
            subjconimgs = subjconimgs(nr_cons+1:end); % we only want to pick the con images corresponding to single trials
            subjvifnames = vifs.name(1,1:end-nr_cons_no_st); % same for vifs
            subjvifvalues = vifs.allvifs(1,1:end-nr_cons_no_st);
            
            for conimg = 1:size(subjconimgs,2)
                for con2ex = 1:size(cons2exclude_dat_st,2)
                    idx_cons(conimg,con2ex) = ~contains(subjconimgs(conimg).name, char(cons2exclude_dat_st{con2ex}));
                end
                idx_cons2(conimg) = logical(sum(idx_cons(conimg,:)) == size(cons2exclude_dat_st,2));
            end
            
            subjconimgs = subjconimgs(idx_cons2'); % excluding cons2exclude_dat_st
            subjvifnames = subjvifnames(idx_cons2);
            subjvifvalues = subjvifvalues(idx_cons2);
            
            clear conimg idx_cons idx_cons2
            
        else
            subjconimgs = subjconimgs(nr_cons+1:end); % we only want to pick the con images corresponding to single trials, but not exclude trials for any conditions
            subjvifnames = vifs.name(1,1:end-nr_cons_no_st);
            subjvifvalues = vifs.allvifs(1,1:end-nr_cons_no_st);
        end
        
        for conimg = 1:size(subjconimgs,2)
            subjconimgs(conimg).path = fullfile(char(subjdir),subjconimgs(conimg).Vcon.fname);
            subjconimgs(conimg).subjname = subjs2use{sub};
        end
        
        conimgs.conpath = [conimgs.conpath;{subjconimgs(:).path}'];
        conimgs.conname = [conimgs.conname;{subjconimgs(:).name}'];
        conimgs.subjname = [conimgs.subjname;{subjconimgs(:).subjname}'];
        conimgs.vifname = [conimgs.vifname;subjvifnames'];
        conimgs.vifvalue = [conimgs.vifvalue;subjvifvalues'];
        conimgs_table = struct2table(conimgs);
        conimgs_table = sortrows(conimgs_table,{'subjname','conname'}); % sort conditions alphabetically!
    end

fmri_dat = fmri_data_st(conimgs_table.conpath);
fmri_dat = remove_empty(fmri_dat);
fmri_dat.source_notes = fprintf('\n loaded single trial con images for %s in %s in fmri_data_st object\n',DSGN.metadata,DSGN.modeldir);
fmri_dat.mask_descrip = fmri_dat.mask.dat_descrip;
fmri_dat.metadata_table.vifvalue = conimgs_table.vifvalue;
fmri_dat.metadata_table.vifname = conimgs_table.vifname;
fmri_dat.metadata_table.conname = conimgs_table.conname;
fmri_dat.metadata_table.subjname = conimgs_table.subjname;

% SANITY CHECK #1

    if ~isequal(height(fmri_dat.metadata_table),size(fmri_dat.dat,2))
        error('\nnumber of vifs (%d rows in fmri_dat.metadata_table) and con images (%d columns in fmri_dat.dat) do not match, please check before proceeding\n',height(fmri_dat.metadata_table),size(fmri_dat.dat,2))
    else
        fprintf('\n');
        warning('passed sanity check #1: number of vifs (%d rows in fmri_dat.metadata_table) and con images (%d columns in fmri_dat.dat) match, proceeding',height(fmri_dat.metadata_table),size(fmri_dat.dat,2))
        fprintf('\n');
    end
    
% SANITY CHECK #2

    for row = height(fmri_dat.metadata_table)
        if ~isequal(fmri_dat.metadata_table.vifname{row}, fmri_dat.metadata_table.conname{row})
            error('\nvifname and conname do not match in row #%d of fmri_data.metadata_table, please check before proceeding\n',row)
        end
    end

fprintf('\n');
warning('passed sanity check #2: vifname and conname match in all rows of fmri_data.metadata_table, proceeding')
fprintf('\n');


%% READ BEHAVIORAL RATINGS AND ADD TO METADATA FIELD OF FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

behav_dat = DAT.BEHAVIOR.behavioral_st_data_table;
    if exist('idx_subjs_st','var')
        behav_dat = behav_dat(idx_subj_st,:);
    end
behav_dat = behav_dat(ismember(behav_dat.(subj_identifier_dat_st),subjs2use),:); % exclude entire subjects to exclude because of lack of rating for all trials within a condition

    if ~isempty(cons2exclude_dat_st)
        behav_dat = behav_dat(~ismember(behav_dat.(cond_identifier_dat_st),cons2exclude_dat_st),:); % exclude conditions if specified
    end
    
behav_dat = sortrows(behav_dat,{subj_identifier_dat_st,cond_identifier_dat_st});

% SANITY CHECK #3

    if ~isequal(height(behav_dat),size(fmri_dat.dat,2))
        error('\nnumber of trials in behavioral data table (%d rows in behav_dat) and con images (%d columns in fmri_dat.dat) do not match, please check before proceeding\n',height(behav_dat),size(fmri_dat.dat,2))
    else
        fprintf('\n');
        warning('passed sanity check #3; number of trials in behavioral data table (%d rows in behav_dat) and con images (%d columns in fmri_dat.dat) match, proceeding',height(behav_dat),size(fmri_dat.dat,2))
        fprintf('\n');
    end    

fmri_dat.metadata_table = [fmri_dat.metadata_table behav_dat];

% SANITY CHECK #4

for row = height(fmri_dat.metadata_table)
    if ~contains(fmri_dat.metadata_table.vifname{row}, fmri_dat.metadata_table.(cond_identifier_dat_st){row})
        error('\nvifname and trial_type do not match in row #%d of fmri_dat.metadata_table, please check before proceeding\n',row)
    end
end

fprintf('\n');
warning('passed sanity check #4: vifname and trial_type do match in all rows of fmri_data.metadata_table, proceeding')
fprintf('\n');

% SANITY CHECK #5

for row = height(fmri_dat.metadata_table)
    if ~isequal(fmri_dat.metadata_table.subjname{row}, fmri_dat.metadata_table.(subj_identifier_dat_st){row})
        error('\nsubject identifiers from fmri_dat and behav_dat do not match in row #%d of fmri_dat.metadata_table, please check before proceeding\n',row)
    end
end

fprintf('\n');
warning('passed sanity check #5: subject identifiers from fmri_dat and behav_dat do match in all rows of fmri_data.metadata_table, proceeding')
fprintf('\n');


%% ADD BEHAVIORAL OUTCOME TO Y FIELD OF FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

fmri_dat.Y = fmri_dat.metadata_table.(behav_outcome_dat_st);
fmri_dat.Y_descrip = behav_outcome_dat_st;
idx_Ynan = ~isnan(fmri_dat.Y);
fmri_dat = get_wh_image(fmri_dat,idx_Ynan); 
% NOTE: @bogpetre's fmri_data_st object nicely excludes the right row in all fields - River Roost IPA earned ;)


%% EXCLUDE TRIALS EXCEEDING VIF THRESHOLD AND MASK
% -------------------------------------------------------------------------

% DEFINE SUBJECT IDENTIFIERS PRIOR TO REMOVING BAD TRIALS

subject_id_vifs = fmri_dat.metadata_table.(subj_identifier_dat_st);
[uniq_subject_id_vifs, ~, subject_id_vifs] = unique(subject_id_vifs,'stable');
n_subj_vifs = size(uniq_subject_id_vifs,1);

% PLOT VIFS AND CALCULATE PERCENTAGE EXCEEDING THRESHOLD DEFINED ABOVE

% over subjects

v1=figure;
hold off
v1 = plot(fmri_dat.metadata_table.vifvalue);
yline(5);
title('vif values over all trials');
xlabel('trial');
ylabel('variance inflation factor');

good_trials_idx = fmri_dat.metadata_table.vifvalue < vif_threshold_dat_st;
bad_trials_perc = sum(~good_trials_idx)./size(fmri_dat.metadata_table.vifvalue,1).*100;
sprintf('%4.2f percent of trials exceeds a vif threshold of %d, indicating multicollinearity with noise regressors; script will remove them',bad_trials_perc,vif_threshold_dat_st);

set(gcf,'WindowState','Maximized');
drawnow, snapnow;

% per subject

v2=figure;

    for sub = 1:n_subj_vifs
        
        this_idx_vifs = find(sub == subject_id_vifs);
        this_vifs = fmri_dat.metadata_table.vifvalue(this_idx_vifs);

        this_good_trials_idx = this_vifs < vif_threshold_dat_st;
        this_bad_trials_perc = sum(~this_good_trials_idx)./size(this_vifs,1).*100;
        sprintf('%4.2f percent of trials for subject %s exceeds a vif threshold of %d, indicating multicollinearity with noise regressors; script will remove them',this_bad_trials_perc,uniq_subject_id_vifs{sub},vif_threshold_dat_st)

        subplot(ceil(sqrt(n_subj_vifs)), ceil(n_subj_vifs/ceil(sqrt(n_subj_vifs))), sub);
        hold off
        v2 = plot(this_vifs);
        yline(5);
        box off
        title(uniq_subject_id_vifs{sub});
        xlabel('trial');
        ylabel('vif');
        
    end
    
set(gcf,'WindowState','Maximized');
drawnow, snapnow;

% REMOVE CON IMAGES CORRESPONDING TO TRIALS EXCEEDING VIF THRESHOLDS FROM
% DATA OBJECT

fmri_dat = fmri_dat.get_wh_image(good_trials_idx);

% CLEAR SUBJECT IDENTIFIERS

clear sub subject_id_vifs uniq_subject_id_vifs n_subj_vifs this_idx_vifs this_vifs;


%% SAVE FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

if ~isempty(cons2exclude_dat_st)
    savefilename = fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome_dat_st, '_exclude_cond_', char([cons2exclude_dat_st{:}]), '_', results_suffix, '.mat']);

else
    savefilename = fullfile(resultsdir, ['single_trial_fmri_data_st_object_', behav_outcome_dat_st, '_', results_suffix, '.mat']);

end

save(savefilename, 'fmri_dat', 'cons2exclude_dat_st', 'behav_outcome_dat_st', 'subj_identifier_dat_st', 'cond_identifier_dat_st','-v7.3');

