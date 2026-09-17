
% Edit these two, and possibly the 3rd, before running
% -------------------------------------------------------
% Remember where the study's own setup put the results, so the call below can be
% checked against it (see the guard immediately after).
resultsdir_before_setup = '';
if exist('resultsdir','var'), resultsdir_before_setup = resultsdir; end

a_set_up_paths_always_run_first


% GUARD: did the path setup just move the output directory?
%
% This line is meant to be replaced, in a study's copy, by that study's own
% s0 (e.g. mystudy_secondlevel_m2a_s0_a_set_up_paths_always_run_first). Left
% as the generic call, it RE-DERIVES resultsdir - typically from the
% FIRST-LEVEL model name - and silently overwrites whatever the study's setup
% had already set. Every result then lands in a different model's directory
% while the published report still goes to the right one, so the split is easy
% to miss. This has happened three times: proj_discoverie's SVM wrote into
% secondlevel/model_2_basic, proj_moodbugs wrote into secondlevel/model_3_basic,
% and all seven core scripts of a new discoverie model were about to do the same.
if ~isempty(resultsdir_before_setup) && ~strcmp(resultsdir_before_setup, resultsdir)
    error(['\nPATH SETUP MOVED THE RESULTS DIRECTORY.\n\n' ...
           '  before: %s\n  after : %s\n\n' ...
           'The generic a_set_up_paths_always_run_first re-derived resultsdir and\n' ...
           'discarded the one your study setup had set. In your copy of this script,\n' ...
           'replace that call with your study''s own s0 path script.\n'], ...
           resultsdir_before_setup, resultsdir);
end

prep_1_set_conditions_contrasts_colors

a2_set_default_options

printhdr('BEHAVIORAL DATA - BETWEEN-PERSON DESIGN')

try
    prep_1b_prep_behavioral_data
catch
    printhdr('Behavioral data not included.');
    disp('prep1b_prep_behavioral_data.m did not run correctly. Either configure and test this or omit this script.');
end

% These should not need editing, 
% can run template scripts directly
% -------------------------------------------------------
%% INDIVIDUAL CONDITION PLOTS

printhdr('INDIVIDUAL CONDITION PLOTS')

prep_2_load_image_data_and_save

%% CONTRAST PLOTS

printhdr('CONTRAST PLOTS')
prep_3_calc_univariate_contrast_maps_and_save

% You can also run z_batch_publish_image_prep_and_qc to run these and
% create an .html file with the output.

%% REGRESSIONS

printhdr('REGRESSIONS')
prep_3a_run_second_level_regression_and_save

prep_3b_run_second_level_regression_on_conditions_and_save

%% SVMs with optional bootstrapping

printhdr('CONTRAST SUPPORT VECTOR MACHINES')
prep_3b_run_SVMs_on_contrasts_and_save

prep_3c_run_SVMs_on_contrasts_masked

prep_3d_run_SVM_betweenperson_contrasts

prep_3e_run_SVM_betweenperson_contrasts_on_conditions

% You can also run z_batch_publish_image_prep_and_qc to run these and
% create an .html file with the output.


%% SIGNATURE PREP

printhdr('SIGNATURE EXTRACTION')
prep_4_apply_signatures_and_save

%% PARCELLATION PREP

try 
    printhdr('PARCELLATIONS')
    prep_5_apply_shen_parcellation_and_save
    prep_5b_apply_spmanatomy_parcellation_and_save
catch 
    warning('Parcellation image atlas not on path. Images stored on Canlab Drive. Contact Canlab if you want to include this step.')
end

%% EMOTION MAPS

printhdr('EMOTION MAPS')
prep_6_apply_kragel_emotion_signatures_and_save

