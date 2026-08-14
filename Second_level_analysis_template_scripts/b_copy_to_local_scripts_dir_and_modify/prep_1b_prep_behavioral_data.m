%% prep_1b_prep_behavioral_data.m
%
%
% *USAGE*
%
% This optional script attaches behavioral/non-imaging data to the DAT
% structure defined by prep_1_set_conditions_contrasts_colors.m. It
%
% # calls a_set_up_paths_always_run_first, then loads DAT from
%   image_names_and_setup.mat if not already in the workspace (falling back
%   to running prep_1_set_conditions_contrasts_colors if that file doesn't
%   exist yet)
% # reads behavioral data from .tsv files in BIDS/phenotype, computes
%   z-scored condition- and contrast-level rating variables, and stores the
%   resulting tables in DAT.BEHAVIOR
% # initializes DAT.BETWEENPERSON.group/.groupnames/.groupcolors and the
%   per-condition/per-contrast DAT.BETWEENPERSON.conditions/.contrasts table
%   cell arrays
% # (CUSTOM CODE section) populates those tables from DAT.BEHAVIOR for this
%   study's specific design
% # resaves image_names_and_setup.mat, appending the new DAT fields (runs
%   git annex unannex first, since this .mat file is typically already
%   tracked by datalad/git-annex by the time this script runs)
%
%
% *WORKED EXAMPLE, NOT A TEMPLATE*
%
% Like prep_1_set_conditions_contrasts_colors.m, this is a real example from
% one LaBGAS study's design, not a generic template - study-specific
% modifications will typically be extensive. Two more worked examples are
% embedded further down in this script itself ("CANLAB EXAMPLE #1",
% commented out, and "CANLAB EXAMPLE #2"/"LABGAS EXAMPLES", pointers to
% further example scripts/repos). See also README.md and LaBGAS's project
% folders on the KU Leuven server and GIN/Github repos.
%
%
% *CANLAB NOTES*
%
% * store behavioral data in any ad hoc format in DAT.BEHAVIOR - useful if
%   you want to write custom scripts relating brain to behavior
% * store a between-person grouping variable (patient vs. control, etc.) in
%   DAT.BETWEENPERSON.group, coded 1/-1; leave empty if you have no such
%   variable. Some analyses use it to run between-group contrasts and/or
%   control for it across the whole sample
% * instead of (or in addition to) a single .group variable, you can enter a
%   different behavioral variable per condition/contrast in
%   DAT.BETWEENPERSON.conditions/.contrasts - cell arrays with one table per
%   condition/contrast
%
%
% *LABGAS NOTES*
%
% * always make a study-specific copy of this script in your code
%   subdataset - do NOT edit in the repo!
% * this script is highly study-specific, but since LaBGAS uses a fixed
%   BIDS-compatible directory structure and phenotype-file format, there
%   should be common elements across studies
% * this example is from a design with different sessions
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   Dartmouth, May, 2022
%
% -------------------------------------------------------------------------
%
% prep_1b_prep_behavioral_data.m         v1.5
%
% last modified: 2026/08/14
%
%
%% RUN SCRIPT A_SET_UP_PATHS_ALWAYS_RUN_FIRST AND LOAD/CREATE DAT IF NEEDED
% -------------------------------------------------------------------------

a_set_up_paths_always_run_first;

if ~exist('DAT','var')
    
    try
    
        load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
    catch
        
        prep_1_set_conditions_contrasts_colors;
        
    end
    
end


%% READ BEHAVIORAL DATA FROM TSV FILES IN BIDS/PHENOTYPE DIR
% -------------------------------------------------------------------------

phenodir = fullfile(BIDSdir,'phenotype');


% LOAD MEAN RATINGS PER CONDITION/SUBJECT

behavioral_data_filename = 'ratings_means.tsv';
behavioral_fname_path = fullfile(phenodir, behavioral_data_filename);

if ~exist(behavioral_fname_path, 'file') 
    fprintf(1, 'CANNOT FIND FILE: %s\n',behavioral_fname_path); 
end

behavioral_data_table = readtable(behavioral_fname_path,'TreatAsEmpty','n/a','FileType','text','Delimiter','tab'); % read .tsv file into Matlab table variable


% CONDITIONS

% NOTES:
%   respect the order of raw variables in table
%   respect the order of DAT.conditions defined in prep_1
%   no z-scoring needed for conditions with all NaN ratings

% z-score simple conditions already present in table
behavioral_data_table.like_bit_high_calorie = zscore(behavioral_data_table.like_bit_high_calorie,0,'omitnan');
behavioral_data_table.like_bit_low_calorie = zscore(behavioral_data_table.like_bit_low_calorie,0,'omitnan');
% behavioral_data_table.like_bit_neutral = zscore(behavioral_data_table.like_bit_neutral,0,'omitnan');
behavioral_data_table.like_pla_high_calorie = zscore(behavioral_data_table.like_pla_high_calorie,0,'omitnan');
behavioral_data_table.like_pla_low_calorie = zscore(behavioral_data_table.like_pla_low_calorie,0,'omitnan');
% behavioral_data_table.like_pla_neutral = zscore(behavioral_data_table.like_pla_neutral,0,'omitnan');
behavioral_data_table.want_bit_high_calorie = zscore(behavioral_data_table.want_bit_high_calorie,0,'omitnan');
behavioral_data_table.want_bit_low_calorie = zscore(behavioral_data_table.want_bit_low_calorie,0,'omitnan');
behavioral_data_table.want_bit_neutral = zscore(behavioral_data_table.want_bit_neutral,0,'omitnan');
behavioral_data_table.want_pla_high_calorie = zscore(behavioral_data_table.want_pla_high_calorie,0,'omitnan');
behavioral_data_table.want_pla_low_calorie = zscore(behavioral_data_table.want_pla_low_calorie,0,'omitnan');
behavioral_data_table.want_pla_neutral = zscore(behavioral_data_table.want_pla_neutral,0,'omitnan');

% calculate and z-score "complex" conditions
behavioral_data_table.like_bit_high_neu = zscore((behavioral_data_table.like_bit_high_calorie - behavioral_data_table.like_bit_neutral),0);
% no omitnan when conditions with NaN ratings are included
behavioral_data_table.like_bit_low_neu = zscore((behavioral_data_table.like_bit_low_calorie - behavioral_data_table.like_bit_neutral),0);
behavioral_data_table.like_bit_high_low = zscore((behavioral_data_table.like_bit_high_calorie - behavioral_data_table.like_bit_low_calorie),0,'omitnan');
behavioral_data_table.like_pla_high_neu = zscore((behavioral_data_table.like_pla_high_calorie - behavioral_data_table.like_pla_neutral),0);
behavioral_data_table.like_pla_low_neu = zscore((behavioral_data_table.like_pla_low_calorie - behavioral_data_table.like_pla_neutral),0);
behavioral_data_table.like_pla_high_low = zscore((behavioral_data_table.like_pla_high_calorie - behavioral_data_table.like_pla_low_calorie),0,'omitnan');
behavioral_data_table.want_bit_high_neu = zscore((behavioral_data_table.want_bit_high_calorie - behavioral_data_table.want_bit_neutral),0,'omitnan');
behavioral_data_table.want_bit_low_neu = zscore((behavioral_data_table.want_bit_low_calorie - behavioral_data_table.want_bit_neutral),0,'omitnan');
behavioral_data_table.want_bit_high_low = zscore((behavioral_data_table.want_bit_high_calorie - behavioral_data_table.want_pla_high_calorie),0,'omitnan');
behavioral_data_table.want_pla_high_neu = zscore((behavioral_data_table.want_pla_high_calorie - behavioral_data_table.want_pla_neutral),0,'omitnan');
behavioral_data_table.want_pla_low_neu = zscore((behavioral_data_table.want_pla_low_calorie - behavioral_data_table.want_pla_neutral),0,'omitnan');
behavioral_data_table.want_pla_high_low = zscore((behavioral_data_table.want_pla_high_calorie - behavioral_data_table.want_pla_low_calorie),0,'omitnan');

% CALCULATE CONTRASTS AND Z-SCORE

% calculate contrasts between conditions for ratings, and zscore them, for
% use as second-level covariates in contrast analyses
% NOTES: 
%   respect the order of raw variables in table
%   respect the order of DAT.contrastnames defined in prep_1
%   include contrasts with conditions for which ratings do not exist (i.e.
%       have NaN), this makes life easier to convert to between table in
%       section below

behavioral_data_table.like_bit_pla_high = zscore((behavioral_data_table.like_bit_high_calorie - behavioral_data_table.like_pla_high_calorie),0,'omitnan');
behavioral_data_table.like_bit_pla_low = zscore((behavioral_data_table.like_bit_low_calorie - behavioral_data_table.like_pla_low_calorie),0,'omitnan');
behavioral_data_table.like_bit_pla_neu = zscore((behavioral_data_table.like_bit_neutral - behavioral_data_table.like_pla_neutral),0);
behavioral_data_table.like_bit_pla_high_neu = zscore((behavioral_data_table.like_bit_high_calorie - behavioral_data_table.like_bit_neutral - behavioral_data_table.like_pla_high_calorie + behavioral_data_table.like_pla_neutral),0);
behavioral_data_table.like_bit_pla_low_neu = zscore((behavioral_data_table.like_bit_low_calorie - behavioral_data_table.like_bit_neutral - behavioral_data_table.like_pla_low_calorie + behavioral_data_table.like_pla_neutral),0);
behavioral_data_table.like_bit_pla_high_low = zscore((behavioral_data_table.like_bit_high_calorie - behavioral_data_table.like_bit_low_calorie - behavioral_data_table.like_pla_high_calorie + behavioral_data_table.like_pla_low_calorie),0,'omitnan');
behavioral_data_table.want_bit_pla_high = zscore((behavioral_data_table.want_bit_low_calorie - behavioral_data_table.want_pla_low_calorie),0,'omitnan');
behavioral_data_table.want_bit_pla_low = zscore((behavioral_data_table.want_bit_low_calorie - behavioral_data_table.want_pla_low_calorie),0,'omitnan');
behavioral_data_table.want_bit_pla_neu = zscore((behavioral_data_table.want_bit_neutral - behavioral_data_table.want_pla_neutral),0,'omitnan');
behavioral_data_table.want_bit_pla_high_neu = zscore((behavioral_data_table.want_bit_high_calorie - behavioral_data_table.want_bit_neutral - behavioral_data_table.want_pla_high_calorie + behavioral_data_table.want_pla_neutral),0,'omitnan');
behavioral_data_table.want_bit_pla_low_neu = zscore((behavioral_data_table.want_bit_low_calorie - behavioral_data_table.want_bit_neutral - behavioral_data_table.want_pla_low_calorie + behavioral_data_table.want_pla_neutral),0,'omitnan');
behavioral_data_table.want_bit_pla_high_low = zscore((behavioral_data_table.want_bit_high_calorie - behavioral_data_table.want_bit_low_calorie - behavioral_data_table.want_pla_high_calorie + behavioral_data_table.want_pla_low_calorie),0,'omitnan');


% TRIAL-BY-TRIAL RATINGS

% NOTE: can be used in single trial analyses including 'signature development'

behavioral_st_data_filename = 'food_images_wanting_liking.tsv';
behavioral_st_fname_path = fullfile(phenodir, behavioral_st_data_filename);

if ~exist(behavioral_st_fname_path, 'file') 
    fprintf(1, 'CANNOT FIND FILE: %s\n',behavioral_st_fname_path); 
end

behavioral_st_data_table = readtable(behavioral_st_fname_path,'TreatAsEmpty','n/a','FileType','text','Delimiter','tab'); % read .tsv file into Matlab table variable


% ADD BOTH RATING TABLES TO DAT

DAT.BEHAVIOR.behavioral_data_table = behavioral_data_table;
DAT.BEHAVIOR.behavioral_st_data_table = behavioral_st_data_table;


%% INITIALIZE GROUP VARIABLE
% -------------------------------------------------------------------------

% INITIALIZE BETWEENPERSON FIELD

DAT.BETWEENPERSON = [];


% INITIALIZE GROUP VARIABLE

% These fields are mandatory, but they can be empty
% If group variable and between-person variables vary by
% condition/contrast, leave these empty

DAT.BETWEENPERSON.group = [];
DAT.BETWEENPERSON.groupnames = {};
DAT.BETWEENPERSON.groupcolors = {};


%% INITIALIZE CONDITION/CONTRAST-SPECIFIC BETWEEN-PERSON DATA TABLES
% -------------------------------------------------------------------------

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

DAT.BETWEENPERSON.contrasts = cell(1, size(DAT.contrastnames,2));
[DAT.BETWEENPERSON.contrasts{:}] = deal(table());  % empty tables


%% CUSTOM CODE: TRANSFORM INTO BETWEEN DESIGN TABLE
% -------------------------------------------------------------------------

% Create a table for each condition/contrast with the between-person design, or leave empty.
%
% a variable called 'id' contains subject identfiers.  Other variables will
% be used as regressors.  Variables with only two levels should be effects
% coded, with [1 -1] values.

id = DAT.BEHAVIOR.behavioral_data_table.participant_id;
covs = DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames(contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'like') | contains(DAT.BEHAVIOR.behavioral_data_table.Properties.VariableNames,'want'));

for cond = 1:size(DAT.conditions,2)/2
    if ~contains(covs{cond},'neutral') && ~contains(covs{cond},'neu') % neutral conditions have NaN for liking
    DAT.BETWEENPERSON.conditions{cond}.liking = DAT.BEHAVIOR.behavioral_data_table.(covs{cond}); 
    end
    DAT.BETWEENPERSON.conditions{cond}.wanting = DAT.BEHAVIOR.behavioral_data_table.(covs{cond+6});
end

for cond = (size(DAT.conditions,2)/2)+1:size(DAT.conditions,2)
    if ~contains(covs{cond+6},'neutral') && ~contains(covs{cond+6},'neu') % neutral conditions have NaN for liking
    DAT.BETWEENPERSON.conditions{cond}.liking = DAT.BEHAVIOR.behavioral_data_table.(covs{cond+6}); 
    end
    DAT.BETWEENPERSON.conditions{cond}.wanting = DAT.BEHAVIOR.behavioral_data_table.(covs{cond+12});
end

for cont = 1:size(DAT.contrasts,1)
    if ~contains(covs{(size(DAT.conditions,2)*2)+cont},'neu') % neutral contrasts have NaN for liking
        DAT.BETWEENPERSON.contrasts{cont}.delta_liking = DAT.BEHAVIOR.behavioral_data_table.(covs{((size(DAT.conditions,2))*2)+cont});
    end
    DAT.BETWEENPERSON.contrasts{cont}.delta_wanting = DAT.BEHAVIOR.behavioral_data_table.(covs{((size(DAT.conditions,2))*2)+cont+6});
end


%% CANLAB EXAMPLE #1
% -------------------------------------------------------------------------

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


%% CANLAB EXAMPLE #2
% -------------------------------------------------------------------------

% "mastergithubrepo"/CANlab_help_examples/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/prep_1b_prep_behavioral_data_example2.m


%% LABGAS EXAMPLES
% -------------------------------------------------------------------------

% github.com/labgas
%
% gin.g-node.org/labgas


%% CHECK DAT, PRINT WARNINGS, SAVE DAT STRUCTURE
% -------------------------------------------------------------------------

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

cd(resultsdir); % unannex image_names_and_setup.mat file if already datalad saved to prevent write permission problems
! git annex unannex image_names_and_setup.mat
cd(rootdir);

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'dashes','printstr','printhdr','DSGN', 'DAT', 'basedir', 'datadir', 'maskdir', 'resultsdir', 'scriptsdir', 'figsavedir', 'htmlsavedir', '-append');
