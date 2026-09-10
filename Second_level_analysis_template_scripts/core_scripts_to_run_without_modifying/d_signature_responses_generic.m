%% d_signature_responses_generic.m
%
%
% *USAGE*
%
% This script plots selected signature responses calculated by
% prep_4_apply_signatures_and_save.m, and tests their significance, for
% conditions and contrasts defined in DAT. It
%
% # calls a_set_up_paths_always_run_first, then loads DAT/DATA_OBJ*/DATA_OBJ_CON*
%   if not already in the workspace
% # checks that DAT.SIG_conditions/DAT.SIG_contrasts exist (i.e. that
%   prep_4_apply_signatures_and_save.m has been run)
% # if signatures_to_plot is empty, defaults to all signatures in
%   keyword_sigs (as computed by prep_4)
% # plots and tests significance for each selected signature by calling
%   plugin_signature_condition_contrast_plot
%
% Run this script headless from the Linux command line (the default), which
% publishes the html report and fails loudly if the script errors:
%
%   labgascore_run_headless.sh -d /data/proj_xxx \
%       -s <proj>_secondlevel_m<M>_s0_a_set_up_paths_always_run_first \
%       <proj>_secondlevel_m<M>_s<N>_d_signature_responses_generic
%
% Or, interactively from the Matlab terminal (use this when you want
% higher-resolution figures, or are debugging):
%
%   LaBGAScore_prov_publish('d_signature_responses_generic', htmlsavedir)
%
% NOTE: publish() catches a script error into the html and returns normally, so
% a crashed run looks exactly like a successful one. Prefer the routes above,
% which read the report back and check for a caught error, over a bare
% publish('d_signature_responses_generic','outputDir',htmlsavedir).
%
%
% *OPTIONS*
%
% NOTE:
%       defaults are specified in a2_set_default_options for any given model,
%       but if you want to run the same model with different options, you can
%       make a copy of this script with a letter index (e.g. _s6a_) and
%       change the default option below
%
% * signatures_to_plot     default empty (plot all signatures in keyword_sigs); cell array of selected signature names, e.g. {'signame1','signame2',...}
%
% -------------------------------------------------------------------------
%
% adapted by: Lukas Van Oudenhove
%
% date:   Leuven, January, 2023
%
% -------------------------------------------------------------------------
%
% d_signature_responses_generic.m         v1.2
%
% last modified: 2026/08/14


%% GET PATHS AND OPTIONS AND CHECK OPTIONS
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT
% NOTE: THIS WILL ALSO AUTOMATICALLY CALL A2_SET_DEFAULT_OPTIONS

% SET CUSTOM OPTIONS

% NOTE: only specify if you want to run multiple versions of your model with different options
% than the defaults you set in your model-specific version of a2_set_default_options.m
% 
% signatures_to_plot = {'varname1','varname2',...};


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
    
    signatures_to_plot = cell(1,size(keyword_sigs,2));
    
    for sig = 1:size(keyword_sigs,2)

        [~,signame] = fileparts(char(keyword_sigs{sig}));
    
            if contains(keyword_sigs{sig},filesep) % path to image rather than keyword

                signatures_to_plot{sig} = DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(signame).signaturenames;
        
            else
        
                signatures_to_plot{sig} = DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs{sig}).signaturenames;
        
            end
        
    end
    
end

plugin_signature_condition_contrast_plot

