%% prep_1_set_conditions_contrasts_colors.m
%
% CANLAB NOTES:
% - Modify to specify image file subdirectories, wildcards to locate images, condition names
%
% LABGAS NOTES:
% - Always make a study-specific copy of this script in your code subdataset, do NOT edit in the repo!
% - Study-specific modifications should in principle be limited to design-related rather
%       than directory structure-related issues since we use a fixed organisation
% - This is an example from a design with four conditions and 
%       three contrasts of (particular) interest, 
%       see preregistration https://osf.io/z6gy2 for design info etc
%       see datalad dataset as GIN repo https://gin.g-node.org/labgas/proj_erythritol_4a
% - For another example, see https://github.com/labgas/proj-emosymp/blob/main/secondlevel/model_1_CANlab_classic_GLM/prep_1_emosymp_m1_s3_set_conditions_contrasts_colors.m
%    
%
%__________________________________________________________________________
%
% modified by: Lukas Van Oudenhove
% date:   Dartmouth, May, 2022
%
%__________________________________________________________________________
% @(#)% prep_1_set_conditions_contrasts_colors.m         v1.0
% last modified: 2022/05/16


%% SET UP CONDITIONS
% ------------------------------------------------------------------------

% EXAMPLE

% conditions = {'C1' 'C2' 'C3' 'etc'};
% structural_wildcard = {'c1*nii' 'c2*nii' 'c3*nii' 'etc*nii'};
% functional_wildcard = {'fc1*nii' 'fc2*nii' 'fc3*nii' 'etc*nii'};
% colors = {'color1' 'color2' 'color3' etc}  One per condition

fprintf('first level beta/con image data should be in /firstleveldir/model_x/sub-xxx organisation\n'); %lukasvo76 adapted to LaBGAS default directory structure

DAT = struct();

% NAMES OF CONDITIONS

% A "condition" usually corresponds to a regressor (beta/COPE image) in a first-level 
% (single-subject) analysis. You will usually have one image per condition per subject.
% You will have the opportunity below to specify contrasts across these
% conditions, which are estimated as part of the prep_ scripts.  
% For example, if you have a task with three types of stimuli, "pain"
% "nausea" and "itch", there would be 3 conditions. 
%
% Within-subjects design:  If this is a within-person design
% (each person experiences all three conditions), you would have 3 images
% per subject in your "data" subfolder for this study. When the images are
% loaded (later, by the prep_2_... script) the subject images should be in
% the same order for each of the three conditions. 
%
% Between-subjects design:  If this is a between-person design, you would have 3 groups 
% of subjects with one image each. You may have different numbers of
% subjects in each condition.
% 
% Designs with multiple events per task condition: Sometimes, you may have
% multiple images corresponding to the same task/trial type. For example,
% you might have 4 levels of "pain", or "early pain" and "late pain".  It
% is recommended to enter each level as a separate "condition". You can
% later create contrasts, e.g., linear contrast across 4 pain levels, or
% "early vs. late". 

% Enter a cell array { } with one cell per condition.  Each cell should
% contain a string specifying the condition name to be used in plots and
% tables. This does not have to correspond to an image/directory name.
% 
% @lukasvo76: it is strongly recommended to have these condition names
% correspond with your first-level condition names defined in
% DSGN.conditions!

DAT.conditions = {'sucrose' 'erythritol' 'sucralose' 'water'};
% DAT.conditions = DSGN.conditions{1}; 
% @lukasvo76: only use if 
% 1) first-level conditions are the same in every run
% 2) all first-level conditions are of interest at second level

DAT.conditions = format_strings_for_legend(DAT.conditions);

% SPECIFYING IMAGE FILE LOCATIONS FOR EACH CONDITION

% The next lines of code specify how to locate the image files associated
% with each condition.  The image load script (i.e., prep_2...m) will attempt to use
% the file system (mac osx/windows/linux) to list the files and store them
% in a cell array of strings, with each string containing the full path name 
% for a valid image.  Images will be listed in the order in which the file
% system lists them.  
% For each condition, the prep_2...m script will construct a string containing
% wildcards (*, ?, [1-9], etc.) that can help you flexibly specify the location of 
% image.  The string is fullfile(datadir, DAT.subfolders{i}, DAT.functional_wildcard{i}),
% where i is the condition number and datadir is the path to the main
% "data" subfolder.
%
% EXAMPLES:
% Here is an example.  Say you have a subfolder in [basedir]/data 
% called 'SAC_NRFX/contrast_2', which contains one image per person, 
% named con_2_subj01.img, con_2_subj02.img, etc.  
% You could enter 'SAC_NRFX/contrast_2' in the .subfolders field below, 
% and 'con*img' in the functional_wildcard field below, and the result
% would be that we'd list files in
% [basedir]/data/SAC_NRFX/contrast_2/con*img.
%
% In a second example, say you have a subfolder in [basedir]/data for each
% subject. Inside that folder are images called "pain.nii" and "itch.nii"
% You could enter '*' in the .subfolders field below, 
% and 'pain.nii' in the functional_wildcard field for condition 1, and "itch.nii"
% for condition 2.  The result
% would be that we'd list files in
% [basedir]/data/*/pain.nii for condition 1, and [basedir]/data/*/itch.nii
% lukasvo76: this is the typical BIDS-compatible situation

% Names of subfolders in [basedir]/data
% Enter a cell array { } with one cell per condition.  Each cell should
% contain a string specifying the subfolder to look in for images of this
% condition. 
% If you do not have subfolders, it is OK to leave this empty, i.e., DAT.subfolders = {};

DAT.subfolders = {'*' '*' '*' '*'}; % @lukasvo76: default option for Linux OS, one wildcard per condition
% DAT.subfolders = {}; % @lukasvo76 fallback option for Windows OS, uses recursive spm_select in prep_2 script to select right con images

% Names of wildcard (expression with *, [1-9], 
% Enter a cell array { } with one cell per condition.  Each cell should
% contain a string specifying the subfolder to look in for images of this
% condition. 

DAT.structural_wildcard = {};
DAT.functional_wildcard = {'con_0001.nii' 'con_0002.nii' 'con_0003.nii' 'con_0004.nii'}; %lukavo76: default option for Linux OS
% DAT.functional_wildcard = {'^con_0001.*\nii$' '^con_0002.*\nii$' '^con_0003.*\nii$' '^con_0004.*\nii$'}; %lukavo76: fallback option for Windows OS, spm_select uses regular expression to filter (like in the GUI)


%% SET UP CONTRASTS
% ------------------------------------------------------------------------

% There are three ways to set up contrasts, which will be displayed as
% maps, run in SVM analyses (if contrast weights are 1 and -1), and used in
% signature and network analyses.
%
% 1. For within-person contrasts, where each individual has an
% image for each condition being compared, use DAT.contrasts, here.
% Important: You must have the same number of images in each condition
% being compared, and the images must be in the SAME SUBJECT ORDER.
% Contrasts are paired tests across these conditions.
% These contrasts should be used if condition is crossed with participant
% (i.e., within-subject design).
% These will be used in c2_SVM_contrasts.m
%
% 2. If your lists of images for each condition include participants from
% different groups, set up prep1b_...behavioral script, which creates
% DAT.BETWEENPERSON.group and group vectors for each condition and
% contrast. These will be used in c2b_SVM_betweenperson_contrasts.m
%
% 3. If conditions being compared include images for different subjects
% i.e., condition{1} and condition{2} include different individuals, 
% use DAT.between_condition_cons below. These contrasts should be used if 
% subjects are nested within conditions (i.e., between-subject design).
% These will be used in c2c_SVM_between_condition_contrasts.

% WITHIN-PERSON CONTRASTS: Vectors across conditions
% There should be one column in DAT.contrasts per condition, and one row
% per contrast. Enter contrasts for within-person comparisons only (matched
% sets of images, where the ith image is from the ith subject for all
% conditions).

DAT.contrasts = [1 0 -1 0; 1 -1 0 0; 0 1 -1 0];
    
% Descriptive names for contrasts to be used in plots and tables. Avoid
% special characters.
DAT.contrastnames = {'sucrose vs sucralose' 'sucrose vs erythritol' 'erythritol vs sucralose'};

DAT.contrastnames = format_strings_for_legend(DAT.contrastnames);


%% SET UP COLORS
% ------------------------------------------------------------------------

% DAT.colors should be a cell array { } with one 3-element rgb vector
% (e.g., [1 0 0] for red) per condition.

% If you do not edit the code below it is ok; sensible default choices are already coded up. 

% There are several options for defining colors for conditions and
% contrasts, or enter your own in a cell array of length(conditions) for
% DAT.colors, and size(contrasts, 1) for DAT.contrastcolors
% It is better if contrasts have distinct colors from conditions

% Some options: scn_standard_colors, custom_colors, colorcube_colors, seaborn_colors, bucknerlab_colors

% DAT.colors = scn_standard_colors(length(DAT.conditions));
% DAT.colors = custom_colors(cm(1, :), cm(end, :), length(DAT.conditions));
% DAT.contrastcolors = custom_colors([.2 .2 .8], [.2 .8 .2], length(DAT.contrasts));

mycolors = colorcube_colors(length(DAT.conditions) + size(DAT.contrasts, 1));

DAT.colors = mycolors(1:length(DAT.conditions));

% DAT.contrastcolors should be a cell array { } with one 3-element rgb vector
% (e.g., [1 0 0] for red) per contrast.

DAT.contrastcolors = mycolors(length(DAT.conditions) + 1:length(mycolors));


disp('SET up conditions, colors, contrasts in DAT structure.');


%% SET UP BETWEEN-CONDITION CONTRASTS, NAMES, AND COLORS
% ------------------------------------------------------------------------

%    If conditions being compared include images for different subjects
%    i.e., condition{1} and condition{2} include different individuals, 
%    enter contrasts in DAT.between_condition_cons below.
%    These will be used in c2c_SVM_between_condition_contrasts.
%    You do not need to have the same number of images in each condition
%    being compared.
%    Contrasts are unpaired tests across these conditions.

% Matrix of [n contrasts x k conditions]

% DAT.between_condition_cons = [1 -1 0;
%                               1 0 -1];
% 
% DAT.between_condition_contrastnames = {'Pain vs Nausea' 'Pain vs Itch'};
%           
% DAT.between_condition_contrastcolors = custom_colors ([.2 .2 .8], [.2 .8 .2], size(DAT.between_condition_cons, 1));
