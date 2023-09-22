%% prep_1b_prep_behavioral_data.m
%
% CANLAB NOTES:
%
% Optional: Run these load and attach behavioral data from files (e.g., from Excel)            
%
% This script is an example script only.  You should modify it to fit your
% needs, which will depend on which types of behavioral/non-imaging data
% you have and what variables you want to store and analyze. The basics are
% desribed here:
%
% - Store behavioral data tables in any ad hoc format in DAT.BEHAVIOR.
% This can be a useful reference if you want to change/add custom analyses
% relating brain to behavior. You can create custom scripts that pull data 
% from .BEHAVIOR and use it in analyses.
%
% - Store a between-person grouping variable (e.g., patient vs. control,
% etc.) in DAT.BETWEENPERSON.group. This should be coded with values of 1
% and -1. Also add fields (see below) for names and colors associated with
% each group, and a description of what the 1/-1 codes mean.
% Some analyses consider this variable and run between-group contrasts
% and/or control for them in analyses of the entire sample (e.g.,
% "signature response" analyses).  
% SVM analyses can also be run that use the .group variable. See:
% prep_3d_run_SVM_betweenperson_contrasts and 
% c2b_SVM_betweenperson_contrasts  
%
% - If you have no binary group variable,  it is OK to leave the .group
% field empty. 
%
% - Instead of a single .group variable to be tested with all
% conditions/contrasts, you can also enter different variables for each
% condition and/or contrast.  This is useful if you want to correlate each
% contrast with a different behavioral variable (maybe, e.g., for each contrast,
% reaction time differences for the same contrast). 
% If so, enter DAT.BETWEENPERSON.conditions and/or
% DAT.BETWEENPERSON.contrasts.  These should be cell arrays with one cell
% per condition or contrast.  Each cell contains a matlab "table" object
% with the data and possibly subject IDs (see below).
% 
% - You can run this script as part of a workflow (prep_1...
% prep_2...prep_3 etc)
% You can also run the script AFTER you've prepped all the imaging data, just
% to add behavioral data to the existing DAT structure.  If so, make sure
% you RELOAD the existing DAT structure with b_reload_saved_matfiles.m
% before you run this script.  Otherwise, if you create a new DAT
% structure, important information saved during the data prep (prep_2...,
% prep_3...) process will be missing, and you will need to re-run the whole
% prep sequence.
%
% LABGAS NOTES:
%
% - Always make a study-specific copy of this script in your code subdataset, do NOT edit in the repo!
% - This script is highly study-specific, but since we use a fixed
%       BIDS-compatible directory structure and format for the phenotype
%       files, there should be common elements across studies
% - This is an example from a design with four conditions and 
%       three contrasts of (particular) interest, 
%       see preregistration https://osf.io/z6gy2 for design info etc
%       see datalad dataset as GIN repo https://gin.g-node.org/labgas/proj_erythritol_4a
% - For an example of a design including group comparisons, see https://github.com/labgas/proj-emosymp/blob/main/secondlevel/model_1_CANlab_classic_GLM/prep_1b_emosymp_m1_s4_prep_behavioral_data.m
%    
%
%__________________________________________________________________________
%
% modified by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% prep_1b_prep_behavioral_data.m         v1.2
% last modified: 2023/09/22
%
%
%% RUN SCRIPT A_SET_UP_PATHS_ALWAYS_RUN_FIRST
% -------------------------------------------------------------------------

bit_rew_secondlevel_m1_s0_a_set_up_paths_always_run_first;


%% READ BEHAVIORAL DATA FROM TSV FILES IN BIDS/PHENOTYPE DIR

phenodir = fullfile(BIDSdir,'phenotype');

% MEAN RATINGS PER RUN, AND OFFLINE RATINGS

behavioral_data_filename = 'ratings_means.tsv';
behavioral_fname_path = fullfile(phenodir, behavioral_data_filename);

if ~exist(behavioral_fname_path, 'file') 
    fprintf(1, 'CANNOT FIND FILE: %s\n',behavioral_fname_path); 
end

behavioral_data_table = readtable(behavioral_fname_path,'TreatAsEmpty','n/a','FileType','text','Delimiter','tab'); % read .tsv file into Matlab table variable

% calculate contrasts between conditions for ratings, and zscore them, for
% use as second-level covariates in contrast analyses
% NOTE: respect the order of raw variables in table
behavioral_data_table.intensity_sucro_sucra = zscore((behavioral_data_table.intensity_sucrose - behavioral_data_table.intensity_sucralose),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.intensity_sucro_ery = zscore((behavioral_data_table.intensity_sucrose - behavioral_data_table.intensity_erythritol),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.intensity_ery_sucra = zscore((behavioral_data_table.intensity_erythritol - behavioral_data_table.intensity_sucralose),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.liking_sucro_sucra = zscore((behavioral_data_table.liking_sucrose - behavioral_data_table.liking_sucralose),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.liking_sucro_ery = zscore((behavioral_data_table.liking_sucrose - behavioral_data_table.liking_erythritol),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.liking_ery_sucra = zscore((behavioral_data_table.liking_erythritol - behavioral_data_table.liking_sucralose),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.rating_sucro_sucra = zscore((behavioral_data_table.rating_sucrose - behavioral_data_table.rating_sucralose),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.rating_sucro_ery = zscore((behavioral_data_table.rating_sucrose - behavioral_data_table.rating_erythritol),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1
behavioral_data_table.rating_ery_sucra = zscore((behavioral_data_table.rating_erythritol - behavioral_data_table.rating_sucralose),0,'omitnan'); % respect the order of DAT.contrastnames defined in prep_1

% zscore ratings per condition too, for use as second-level covariates in condition
% analyses
behavioral_data_table.intensity_sucrose = zscore(behavioral_data_table.intensity_sucrose,0,'omitnan');
behavioral_data_table.intensity_erythritol = zscore(behavioral_data_table.intensity_erythritol,0,'omitnan');
behavioral_data_table.intensity_sucralose = zscore(behavioral_data_table.intensity_sucralose,0,'omitnan');
behavioral_data_table.liking_sucrose = zscore(behavioral_data_table.liking_sucrose,0,'omitnan');
behavioral_data_table.liking_erythritol = zscore(behavioral_data_table.liking_erythritol,0,'omitnan');
behavioral_data_table.liking_sucralose = zscore(behavioral_data_table.liking_sucralose,0,'omitnan');
behavioral_data_table.rating_sucrose = zscore(behavioral_data_table.rating_sucrose,0,'omitnan');
behavioral_data_table.rating_erythritol = zscore(behavioral_data_table.rating_erythritol,0,'omitnan');
behavioral_data_table.rating_sucralose = zscore(behavioral_data_table.rating_sucralose,0,'omitnan');
behavioral_data_table.rating_water = zscore(behavioral_data_table.rating_water,0,'omitnan');


% ONLINE RATINGS PER TRIAL
% can be used in single trial analyses including 'signature development'
behavioral_st_data_filename = 'ratings_online.tsv';
behavioral_st_fname_path = fullfile(phenodir, behavioral_st_data_filename);

if ~exist(behavioral_st_fname_path, 'file') 
    fprintf(1, 'CANNOT FIND FILE: %s\n',behavioral_st_fname_path); 
end

behavioral_st_data_table = readtable(behavioral_st_fname_path,'TreatAsEmpty','n/a','FileType','text','Delimiter','tab'); % read .tsv file into Matlab table variable

% Add both to DAT for record, and flexible use later
DAT.BEHAVIOR.behavioral_data_table = behavioral_data_table;
DAT.BEHAVIOR.behavioral_st_data_table = behavioral_st_data_table;


%% INITIALIZE GROUP VARIABLE

% Initialize BETWEENPERSON field in DAT structure
DAT.BETWEENPERSON = [];

% Single group variable, optional, for convenience
% These fields are mandatory, but they can be empty
% If group variable and between-person variables vary by
% condition/contrast, leave these empty
% -------------------------------------------------------------------------
DAT.BETWEENPERSON.group = [];
DAT.BETWEENPERSON.groupnames = {};
DAT.BETWEENPERSON.groupcolors = {};


%% INITIALIZE CONDITION/CONTRAST-SPECIFIC BETWEEN-PERSON DATA TABLES

% Cell array with table of group (between-person) variables for each
% condition, and for each contrast.
% If variables are entered:
% 1. they will be controlled for in analyses of the
% overall condition/contrast effects.
% 2. Analyses relating conditions/contrasts to these variables will be
% performed.
% If no variables are entered (empty elements), only
% within-person/whole-group effects will be analyzed.

% Between-person variables influencing each condition
% Table of [n images in condition x q variables]
% names in table can be any valid name.

DAT.BETWEENPERSON.conditions = cell(1, length(DAT.conditions));
[DAT.BETWEENPERSON.conditions{:}] = deal(table());  % empty tables

% Between-person variables influencing each condition
% Table of [n images in contrast x q variables]
% names in table can be any valid name.

DAT.BETWEENPERSON.contrasts = cell(1, length(DAT.contrastnames));
[DAT.BETWEENPERSON.contrasts{:}] = deal(table());  % empty tables


%% CUSTOM CODE: TRANSFORM INTO between_design_table
%
% Create a table for each condition/contrast with the between-person design, or leave empty.
%
% a variable called 'id' contains subject identfiers.  Other variables will
% be used as regressors.  Variables with only two levels should be effects
% coded, with [1 -1] values.
%
id = DAT.BEHAVIOR.behavioral_data_table.participant_id;
covs = DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames(contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'intensity') | contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'rating'));

for cond = 1:size(DAT.conditions,2)
    if cond < size(DAT.conditions,2)
    DAT.BETWEENPERSON.conditions{cond}.intensity = DAT.BEHAVIOR.behavioral_data_table.(covs{cond}); % we include the intensity ratings for the three non-water conditions here, to be able to include them as covariates in analyses on conditions later;
    end
    DAT.BETWEENPERSON.conditions{cond}.rating = DAT.BEHAVIOR.behavioral_data_table.(covs{cond+4}); % same for ratings
end

for cont = 1:size(DAT.contrastnames,1)
    DAT.BETWEENPERSON.contrasts{cont}.delta_intensity = DAT.BEHAVIOR.behavioral_data_table.(covs{(size(DAT.conditions,2)*2)-1+cont});
    DAT.BETWEENPERSON.contrasts{cont}.delta_rating = DAT.BEHAVIOR.behavioral_data_table.(covs{(size(DAT.conditions,2)*2)-1+cont+3});
end

%% CANLAB EXAMPLE #1
%
% % e.g., Vars of interest
% %
% %   12×1 cell array
% % 
% %     'Subject'
% %     'SALINESEDATION'
% %     'REMISEDATION'
% %     'ORDER'
% %     'MSALINEINTENSITY'
% %     'MSALINEUNPL'
% %     'MREMIINTENSITY'
% %     'MREMIUNPL'
% %     'SSALINEINTENSITY'
% %     'SSALINEUNPL'
% %     'SREMIINTENSITY'
% %     'SREMIUNPL'
% % 
% %   6×1 cell array
% % 
% %     'Sal vs Remi ResistStrong'
% %     'Sal vs Remi AntStrong'
% %     'AntLinear Sal'
% %     'AntLinear Remi'
% %     'ResistStrong vs Weak Sal'
% %     'ResistStrong vs Weak Remi'
% 
% % Enter the behavioral grouping codes (1, -1) for individual differences in each 
% % condition or contrast here.  You can enter them separately for
% % conditions/contrasts, or skip this code block and just enter a single
% % variable in DAT.BETWEENPERSON.group below, which will be used for all
% % conditions and contrasts.
% %
% % Contrast numbers below refer to the contrasts you entered and named in the prep_1_
% % script.  The numbers refer to the rows of the contrast matrix you entered.
% %
% % Often, you will want to use the same behavioral variable for multiple
% % contrasts.  If  you have different between-person covariates for different
% % contrasts, e.g., I-C RT for I-C images and ratings or diagnosis for
% % emotion contrasts, you can enter different behavioral variables in the
% % different cells.  Ditto for conditions.
% 
% id = behavioral_data_table.Subject;
% 
% Cov = behavioral_data_table.SSALINEINTENSITY - behavioral_data_table.SREMIINTENSITY;
% covname = 'Intensity_Saline_vs_Remi';
% mytable = table(id);
% mytable.(covname) = Cov;
% 
% DAT.BETWEENPERSON.contrasts{1} = mytable;
% 
% Cov = behavioral_data_table.SSALINEINTENSITY - behavioral_data_table.SREMIINTENSITY;
% covname = 'Intensity_Saline_vs_Remi';
% mytable = table(id);
% mytable.(covname) = Cov;
% 
% DAT.BETWEENPERSON.contrasts{2} = mytable;
% 
% Cov = behavioral_data_table.SSALINEINTENSITY - behavioral_data_table.MSALINEINTENSITY;
% covname = 'Intensity_Str_vs_Mild_Saline';
% mytable = table(id);
% mytable.(covname) = Cov;
% 
% DAT.BETWEENPERSON.contrasts{3} = mytable;
% 
% Cov = behavioral_data_table.SREMIINTENSITY - behavioral_data_table.MREMIINTENSITY;
% covname = 'Intensity_Str_vs_Mild_Remi';
% mytable = table(id);
% mytable.(covname) = Cov;
% 
% DAT.BETWEENPERSON.contrasts{4} = mytable;
% 
% Cov = behavioral_data_table.SSALINEINTENSITY - behavioral_data_table.MSALINEINTENSITY;
% covname = 'Intensity_Str_vs_Mild_Saline';
% mytable = table(id);
% mytable.(covname) = Cov;
% 
% DAT.BETWEENPERSON.contrasts{5} = mytable;
% 
% Cov = behavioral_data_table.SREMIINTENSITY - behavioral_data_table.MREMIINTENSITY;
% covname = 'Intensity_Str_vs_Mild_Remi';
% mytable = table(id);
% mytable.(covname) = Cov;
% 
% DAT.BETWEENPERSON.contrasts{6} = mytable;
% 
% % Single group variable, optional, for convenience
% % These fields are mandatory, but they can be empty
% %
% % Make sure:
% % - DAT.BETWEENPERSON.group is a numeric vector,  coded 1, -1
% % - If you have string inputs, the function string2indicator can help
% % transform them to numbers
% % - If you have numeric codes that are not 1, -1 then contrast_code can
% % help recode them.
% %
% % Group can be empty.
% % -------------------------------------------------------------------------
% group = behavioral_data_table.ORDER;
% [~, nms, group] = string2indicator(group);
% 
% DAT.BETWEENPERSON.group = contrast_code(group  - 1.5);
% DAT.BETWEENPERSON.group_descrip = '-1 is first group name, 1 is 2nd';
% DAT.BETWEENPERSON.groupnames = nms;
% DAT.BETWEENPERSON.groupcolors = {[.7 .3 .5] [.3 .5 .7]};
%
%
%% CANLAB EXAMPLE #2
%
% see "mastergithubrepo"/CANlab_help_examples/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/prep_1b_prep_behavioral_data_example2.m
%
%
%% LABGAS EXAMPLES
%
% see github.com/labgas

%% CHECK DAT, PRINT WARNINGS, SAVE DAT STRUCTURE

if ~isfield(DAT, 'conditions') 
    printhdr('Incomplete DAT structure');
    disp('The DAT field is incomplete. Run prep_1_set_conditions_contrasts_colors before running prep_1b...')
end
    
if isfield(DAT, 'SIG_conditions') && isfield(DAT, 'gray_white_csf')
    % Looks complete, we already have data, no warnings 
else
    printhdr('DAT structure ready for data prep');
    disp('DAT field does not have info from prep_2, prep_3, or prep_4 sequences');
    disp('prep_2/3/4 scripts should be run before generating results.');
end

printhdr('Save DSGN & DAT structures and directory names in image_names_and_setup.mat');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'dashes','printstr','printhdr','DSGN', 'DAT', 'basedir', 'datadir', 'maskdir', 'resultsdir', 'scriptsdir', 'figsavedir', 'htmlsavedir', '-v7.3');
