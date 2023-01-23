%% d_signature_responses_generic
%
%
% USAGE
%
% This script plots selected signature responses calculated in
% prep_4_apply_signatures_and_save, and test their significance for
% conditions and contrasts defined in DAT, by calling
% plugin_signature_condition_contrast_plot
%
%
% OPTIONS
% 
% NOTE: 
% defaults are specified in a2_set_default_options for any given model,
% but if you want to run the same model with different options, 
% you can make a copy of this script with a letter index (e.g. _s6a_) 
% and change the default option below
% 
% signatures_to_plot = {'signame1','signame2',...};   
%
%
%__________________________________________________________________________
%
% adapted by: Lukas Van Oudenhove
% date:   Leuven, January, 2023
%
%__________________________________________________________________________
% @(#)% d_signature_responses_generic.m         v1.0
% last modified: 2023/01/23


%% GET PATHS AND OPTIONS AND CHECK OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% GET DEFAULT OPTIONS IF NOT SET IN A2_SET_DEFAULT_OPTIONS

options_needed = {'signatures_to_plot'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {{}};          % defaults if we cannot find info in a2_set_default_options at all 

plugin_get_options_for_analysis_script

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m
% 
signatures_to_plot = {'NCS','NCSfood','NCSdrugs'};


%% LOAD NECESSARY VARIABLES IF NEEDED
% -------------------------------------------------------------------------

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

if ~exist('DATA_OBJ','var') || ~exist('DATA_OBJsc','var')
    
    load(fullfile(resultsdir,'data_objects.mat'));
    load(fullfile(resultsdir,'data_objects_scaled.mat'));
    
end

if ~exist('DATA_OBJ_CON','var') || ~exist('DATA_OBJ_CONsc','var') || ~exist('DATA_OBJ_CONscc','var')
    
    load(fullfile(resultsdir,'contrast_data_objects.mat'));
    
end


%% CHECK REQUIRED DAT FIELDS
% -------------------------------------------------------------------------

% List required fields in DAT, in cell array
    
required_fields = {'SIG_conditions','SIG_contrasts'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end


%% PLOT AND TEST SIGNIFICANCE
%--------------------------------------------------------------------------

if isempty(signatures_to_plot)
    
    signatures_to_plot = DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signaturenames;
    
end

plugin_signature_condition_contrast_plot
