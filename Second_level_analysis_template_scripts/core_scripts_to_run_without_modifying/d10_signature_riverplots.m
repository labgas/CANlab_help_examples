%% d10_signature_riverplots.m
%
%
% *USAGE*
%
% This script generates riverplots (cosine similarity) for selected
% signature responses against condition and contrast images in DAT. It
%
% # calls a_set_up_paths_always_run_first, then loads DAT/DATA_OBJ*/DATA_OBJ_CON*
%   if not already in the workspace
% # loads signature image sets via load_image_set() - either keyword_sigs
%   (if signatures_to_plot is empty) or each element of signatures_to_plot
% # plots a condition riverplot and a contrast riverplot per signature set,
%   showing statistically significant associations only
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('d10_signature_riverplots','outputDir',htmlsavedir)
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
% * signatures_to_plot     default empty (uses keyword_sigs); cell array of signature GROUP keywords passed individually into load_image_set, e.g. {'signame1','signame2',...}
%
%       NOTE: unlike d_signature_responses_generic.m, this script only works
%       on groups of signatures as defined in load_image_set, not on
%       individual signatures
%
% -------------------------------------------------------------------------
%
% adapted by: Lukas Van Oudenhove
%
% date:   Leuven, January, 2023
%
% -------------------------------------------------------------------------
%
% d10_signature_riverplots.m         v2.1
%
% last modified: 2026/08/14


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


%% LOAD SIGNATURE MAPS
% -------------------------------------------------------------------------
if isempty(signatures_to_plot)
    
    [signatures_obj{1}, signames, ~] = load_image_set(keyword_sigs);
    signatures_obj{1}.image_names = signames;
    
else
    signatures_obj = cell(1,size(signatures_to_plot,2));
    for sig = 1:size(signatures_to_plot,2)
        [signatures_obj{sig}, signames, ~] = load_image_set(signatures_to_plot{sig});
        signatures_obj{sig}.image_names = signames;
    end
    
end


%% CONDITION RIVERPLOTS
% -------------------------------------------------------------------------

printhdr('Cosine Similarity : All conditions');

% k = size(DAT.conditions, 2); % for old code below

% NAME THE CONDITIONS

    for i = 1:length(DATA_OBJ)
        
        DATA_OBJ{i}.image_names = DAT.conditions{i}; 
        
    end

% PLOT SIGNIFICANT ASSOCIATIONS ONLY

    for sig = 1:size(signatures_obj,2)
        
        riverplot(DATA_OBJ, 'layer2', signatures_obj{sig}, 'pos', 'significant_only', 'layer1colors', DAT.colors, 'layer2colors', seaborn_colors(size(signatures_obj{sig}.image_names,2)));
        hh=figure(sig);
        figtitle = 'CANlab signatures riverplot of conditions';
        set(hh, 'Tag', figtitle); plugin_set_figure_size('fig', hh);
        
    end

% Old way: not statistically thresholded, but works:
% Get mean data across subjects
%     m = mean(DATA_OBJ{1});
%     m.image_names = DAT.conditions{1};
%
%     for i = 2:k
%
%         m = cat(m, mean(DATA_OBJ{i}));
%         m.image_names = strvcat(m.image_names, DAT.conditions{i});
%
%     end
%
%
%     riverplot(m, 'layer2', npsplus, 'pos', 'layer1colors', DAT.colors, 'layer2colors', seaborn_colors(length(netnames)), 'thin');
%     pause(2)

% plugin_save_figure;

clear hh figtitle


%% CONTRAST RIVERPLOTS
% -------------------------------------------------------------------------

printhdr('Cosine Similarity : All contrasts');

% k = size(DAT.contrasts, 1); % for old code below

% NAME THE CONTRASTS

    for i = 1:length(DATA_OBJ_CON)

        DATA_OBJ_CON{i}.image_names = DAT.contrastnames{i};

    end

% PLOT SIGNIFICANT ASSOCIATIONS ONLY

    for sig = 1:size(signatures_obj,2)
        
        riverplot(DATA_OBJ_CON, 'layer2', signatures_obj{sig}, 'pos', 'significant_only', 'layer1colors', DAT.contrastcolors, 'layer2colors', seaborn_colors(size(signatures_obj{sig}.image_names,2)));
        hh=figure(sig+size(signatures_obj,2));
        figtitle = 'CANlab signatures riverplot of conditions';
        set(hh, 'Tag', figtitle); plugin_set_figure_size('fig', hh);
        
    end

% Old way: not statistically thresholded, but works:
% Get mean data across subjects
%     m = mean(DATA_OBJ_CON{1});
%     m.image_names = DAT.contrastnames{1};
%
%     for i = 2:k
%
%         m = cat(m, mean(DATA_OBJ_CON{i}));
%         m.image_names = strvcat(m.image_names, DAT.contrastnames{i});
%
%     end
%
%     riverplot(m, 'layer2', npsplus, 'pos', 'layer1colors', DAT.contrastcolors, 'layer2colors', seaborn_colors(length(netnames)), 'thin');

% plugin_save_figure;
