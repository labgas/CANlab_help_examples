%% prep_3g_create_fmri_data_runwise_contrast_object
%
%
% USAGE
%
% This script creates and saves an fmri_data_st object with runwise contrast 
% images calculated from condition beta images
% rootdir/firstlevel/model_x_yyy/sub-zz,
% and adds a convenient metadata_table field containing runwise ratings
% from the corresponding contrast from a data file BIDS/phenotype.tsv
%
% OPTIONS
%
% NOTE: defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options (for example
% voxel- and parcelwise regression), you can make a copy of this script with
% a letter index (e.g. _s6a_) and change the default option here
%
% phenofile_dat_rw:         name of phenotype file in BIDS subdataset
% cons2include_dat_rw:      cell array of (maximum 2) condition names to include (as they appear in DSGN/DAT.conditions as well as SPM.Vbeta.descrip), separated by commas (or blanks)
% behav_outcome_dat_rw:     name of outcome variable in phenotype file
% subj_identifier_dat_rw:   name of subject identifier variable in phenotype file
% run_included_dat_rw;      name of index variable in phenotype file identifying runs for which imaging data have been excluded during QC
% group_identifier_dat_rw:  name of group identifier variable in phenotype file; leave commented out if you don't have groups
%
% MANDATORY OPTIONS TO BE SPECIFIED IN THIS SCRIPT
%
% results_suffix: name to add to results file to specify model in case of
% multiple models
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   KU Leuven, December, 2024
%__________________________________________________________________________
% @(#)% prep_3g_create_fmri_data_runwise_contrast_object.m     v1.0       
% last modified: 2025/01/07


%% GET AND SET OPTIONS
% -------------------------------------------------------------------------

% SET MANDATORY OPTIONS

results_suffix = ''; % adds a suffix of your choice to .mat file with results that will be saved
% NOTE: do NOT delete, leave empty if not needed
% NOTE: do NOT use to add a suffix specifying the behavioral outcome nor included conditions, this will be added automatically

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;
% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT!
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m

% phenofile_dat_rw = 'phenotype.tsv';                   % name of phenotype file in BIDS subdataset
% cons2include_dat_rw = {'condname1','condname2'};      % cell array of (maximum 2) condition names to include (as they appear in DSGN/DAT.conditions as well as SPM.Vbeta.descrip), separated by commas (or blanks)
% behav_outcome_dat_rw = 'varnamex';                    % name of outcome variable in phenotype file
% subj_identifier_dat_rw = 'varnamey';                  % name of subject identifier variable in phenotype file
% run_included_dat_rw = 'varnamez';                     % name of index variable in phenotype file identifying runs for which imaging data have been excluded during QC
% group_identifier_dat_rw = 'varnamezz';                % name of group identifier variable in phenotype file; leave commented out if you don't have groups


%% LOAD NECESSARY VARIABLES IF NEEDED AND DO PREP WORK
% -------------------------------------------------------------------------
    
if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));

end

nr_runs = size(rundirnames,1);
contrast_name = [cons2include_dat_rw{1} ' vs ' cons2include_dat_rw{2}];

% [~,subjs] = fileparts(DSGN.subjects); % assuming use of all subjects for whom firstlevel was run in this particular model
subjs = firstsubjs';                % use if subjects have been excluded at firstlevel stage; CAUTION: you will also need to exclude them from DAT.BEHAVIOR.behavioral_data_table to avoid mismatch, use lines below (UNTESTED)! 
behavioral_data_table = readtable(fullfile(BIDSdir,phenofile_dat_rw));
idx_subjs = contains(behavioral_data_table.(subj_identifier_dat_rw),subjs);
behavioral_data_table = behavioral_data_table(idx_subjs,:);
behavioral_data_table = behavioral_data_table((behavioral_data_table.Time > 1),:); % exclude baseline timepoint

subjs_behav = unique(behavioral_data_table.PPID);
subjs_2use = intersect(subjs,subjs_behav);
firstsubjdirs_2use = firstsubjdirs(contains(firstsubjdirs,subjs_2use));

behavioral_data_tables_subjs = cell(size(subjs_2use,1),1);
idx_subjs_behav = cell(size(subjs_2use,1),1);
idx_subjs_mri_qc = cell(size(subjs_2use,1),1);
idx_subjs_combined = cell(size(subjs_2use,1),1);
idx_subjs_betas = cell(size(subjs_2use,1),1);
run_numbers_subjs = cell(size(subjs_2use,1),1);
behav_ratings_subjs = cell(size(subjs_2use,1),1);
betas_subjs = cell(size(subjs_2use,1),1); 
betas_objs_subjs = cell(size(subjs_2use,1),1); 
con_objs_subjs = cell(size(subjs_2use,1),1); 
image_names_subjs = cell(size(subjs_2use,1),1); 
cons_objs_subjs = cell(size(subjs_2use,1),1);
images_names_subjs = cell(size(subjs_2use,1),1);

subj_names_lengths = cellfun(@length, subjs_2use);
max_length_subj_name = max(subj_names_lengths);
max_length_image_name = (size(contrast_name,2) + max_length_subj_name + 6);


%% CREATE FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

% LOOP OVER SUBJECTS

for sub = 1:size(subjs_2use,1)
    
    behavioral_data_tables_subjs{sub} = behavioral_data_table(contains(behavioral_data_table.PPID,subjs_2use{sub}),:); 
    idx_subjs_behav{sub} = ~isnan(behavioral_data_tables_subjs{sub}.(behav_outcome_dat_rw));
    idx_subjs_mri_qc{sub} = (behavioral_data_tables_subjs{sub}.(run_included_dat_rw) == 1);
    idx_subjs_combined{sub} = idx_subjs_behav{sub} & idx_subjs_mri_qc{sub};
    
    behavioral_data_tables_subjs{sub} = behavioral_data_tables_subjs{sub}(idx_subjs_combined{sub},:);
    run_numbers_subjs{sub} = behavioral_data_tables_subjs{sub}.Time - 1;

       if sum(idx_subjs_behav{sub}) == 0 % all behavioral outcome values are NaN

            fprintf('\n');
            warning('subject %s is missing behavioral data for all runs - excluding subject from analysis',subjs_2use{sub})
            fprintf('\n');

            continue

       else

            behav_ratings_subjs{sub} = behavioral_data_tables_subjs{sub}.(behav_outcome_dat_rw); 

            if sum(idx_subjs_behav{sub}) < nr_runs
                
                runs_excluded_behav = true;

                fprintf('\n');
                warning('subject %s is missing behavioral data for %d run(s) - including subject in analysis',subjs_2use{sub},(nr_runs - sum(idx_subjs_behav{sub})))
                fprintf('\n');
                
            else
                
                runs_excluded_behav = false;

            end

       end
            
        if sum(idx_subjs_mri_qc{sub}) < nr_runs

            runs_excluded_mri_qc = true;
            
            fprintf('\n');
            warning('subject %s had %d run(s) excluded during QC - including subject in analysis',subjs_2use{sub},(nr_runs - sum(idx_subjs_mri_qc{sub})))
            fprintf('\n');
            
        else
            
            runs_excluded_mri_qc = false;

        end

   load(fullfile(firstsubjdirs_2use{sub},'SPM.mat'));
   betas_subjs{sub} = SPM.Vbeta(contains({SPM.Vbeta.descrip}',cons2include_dat_rw'));

   if ~runs_excluded_mri_qc
       
        idx_subjs_betas{sub} = repelem(idx_subjs_combined{sub},2);   
        betas_subjs{sub} = betas_subjs{sub}(idx_subjs_betas{sub});
        
   else
       
       if runs_excluded_behav % this (rather unlikely) scenario where a subject has both runs with NaN behavioral outcome values and runs excluded during QC has not been tested
       
            idx_subjs_betas{sub} = repelem(idx_subjs_combined{sub}(run_numbers_subjs{sub},:),2);  
            betas_subjs{sub} = betas_subjs{sub}(idx_subjs_betas{sub}); 
       end
        
   end
 
   betas_objs_subjs{sub} = cell(size(betas_subjs{sub},2),1);
   con_objs_subjs{sub} = cell(size(betas_subjs{sub},2)/2,1);
   image_names_subjs{sub} = cell(size(betas_subjs{sub},2)/2,1);

       for beta = 1:size(betas_subjs{sub},2)
           betas_objs_subjs{sub}{beta} = fmri_data_st(fullfile(firstsubjdirs_2use{sub},betas_subjs{sub}(beta).fname),which('brain_mask_fmriprep20_template_1000.nii'));
       end

       for con = 1:size(betas_subjs{sub},2)/2
           con_objs_subjs{sub}{con} = image_math(betas_objs_subjs{sub}{(con*2)-1},betas_objs_subjs{sub}{con*2},'minus');
           image_names_subjs{sub}{con} = [subjs_2use{sub} ' ' contrast_name ' run ' num2str(run_numbers_subjs{sub}(con,:))];
       end

   cons_objs_subjs{sub} = cat(con_objs_subjs{sub}{:});
   
   matlabdir = matlabroot;
   cd(fullfile(matlabdir,'toolbox','matlab','strfun'));
   pad_strfun = str2func('pad'); % LVO: hacky solution because Matlab built-in function pad.m is superseded by the Osprey function with the same name in the Matlab path on LaBGAS server, which is needed to have Osprey function properly for unknown reasons
   cd(rootdir);
   
   for i = 1:size(image_names_subjs{sub},1)
       if size(image_names_subjs{sub}{i},2) < max_length_image_name
            image_names_subjs{sub}{i} = pad_strfun(image_names_subjs{sub}{i},max_length_image_name,'right');
       end
   end
   
   images_names_subjs{sub} = vertcat(image_names_subjs{sub}{:});
   
end

% GET RID OF EMPTY CELLS TO MAKE CONCATENATION POSSIBLE

nonEmptyCells_cons = ~cellfun('isempty', cons_objs_subjs);
cons_objs_subjs = cons_objs_subjs(nonEmptyCells_cons);

% CREATE OBJECT

fmri_dat = [cons_objs_subjs{:}];
fmri_dat.dat_descrip = cons2include_dat_rw;
fmri_dat.removed_images = logical([]);
fmri_dat.files_exist = logical([]);
fmri_dat.metadata_table = vertcat(behavioral_data_tables_subjs{:});
fmri_dat.Y = vertcat(behav_ratings_subjs{:});
fmri_dat.Y_descrip = behav_outcome_dat_rw;
fmri_dat.image_names = vertcat(images_names_subjs{:});

    
%% SANITY CHECKS ON FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

% SANITY CHECK #1

    if ~isequal(height(fmri_dat.metadata_table),size(fmri_dat.dat,2))
        error('\ntotal number of runs in behavioral data table (%d rows in fmri_dat.metadata_table) and con images (%d columns in fmri_dat.dat) do not match, please check before proceeding\n',height(fmri_dat.metadata_table),size(fmri_dat.dat,2))
    else
        fprintf('\n');
        warning('passed sanity check #1: total number of runs in behavioral data table (%d rows in fmri_dat.metadata_table) and con images (%d columns in fmri_dat.dat) match, proceeding',height(fmri_dat.metadata_table),size(fmri_dat.dat,2))
        fprintf('\n');
    end
    
% SANITY CHECK #2

    if ~isequal(size(fmri_dat.Y,1),size(fmri_dat.dat,2))
        error('\nnumber of run-wise behavioral ratings in outcome variable (%d rows in fmri_dat.Y) and con images (%d columns in fmri_dat.dat) do not match, please check before proceeding\n',size(fmri_dat.Y,1),size(fmri_dat.dat,2))
    else
        fprintf('\n');
        warning('passed sanity check #2: number of run-wise behavioral ratings in outcome variable (%d rows in fmri_dat.Y) and con images (%d columns in fmri_dat.dat) match, proceeding',size(fmri_dat.Y,1),size(fmri_dat.dat,2))
        fprintf('\n');
    end
    
% SANITY CHECK #3

    if ~isequal(size(fmri_dat.image_names,1),size(fmri_dat.dat,2))
        error('\nnumber of run-wise contrast image names (%d rows in fmri_dat.image_names) and con images (%d columns in fmri_dat.dat) do not match, please check before proceeding\n',size(fmri_dat.image_names,1),size(fmri_dat.dat,2))
    else
        fprintf('\n');
        warning('passed sanity check #3: number of run-wise behavioral ratings in outcome variable (%d rows in fmri_dat.Y) and con images (%d columns in fmri_dat.dat) match, proceeding',size(fmri_dat.image_names,1),size(fmri_dat.dat,2))
        fprintf('\n');
    end
    
% SANITY CHECK #4

    for row = height(fmri_dat.metadata_table)
        if ~contains(fmri_dat.image_names(row,:), fmri_dat.metadata_table.(subj_identifier_dat_rw){row})
            error('\nsubject ID in fmri_data.metadata_table and fmri_data.image_names do not match in row #%d, please check before proceeding\n',row)
        end
    end

fprintf('\n');
warning('passed sanity check #4: subject ID in fmri_data.metadata_table and fmri_data.image_names match in all rows of fmri_data.metadata_table, proceeding')
fprintf('\n');


%% SAVE FMRI_DATA_ST OBJECT
% -------------------------------------------------------------------------

savefilename = fullfile(resultsdir, ['runwise_fmri_data_st_object_', behav_outcome_dat_rw, '_', contrast_name, '_', results_suffix, '.mat']);


save(savefilename, 'fmri_dat', 'phenofile_dat_rw', 'cons2include_dat_rw', 'behav_outcome_dat_rw', 'subj_identifier_dat_rw','-v7.3');

