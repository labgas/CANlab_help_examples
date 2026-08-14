%% prep_4_apply_signatures_and_save.m
%
%
% *USAGE*
%
% This script
%
% # calls a_set_up_paths_always_run_first, then loads DAT/DATA_OBJ*/DATA_OBJ_CON*
%   from image_names_and_setup.mat/data_objects*.mat/contrast_data_objects.mat
%   if not already in the workspace
% # calculates selected CANlab signature responses for all conditions and
%   contrasts in DAT, calling apply_all_signatures(), and saves them to
%   DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(signature_name)
%   and the equivalent DAT.SIG_contrasts field
% # if 'nps' is among keyword_sigs (or keyword_sigs = 'all'), additionally
%   computes NPS subregion responses for conditions and contrasts via
%   apply_nps(), saved to DAT.npsresponse/DAT.npscontrasts/DAT.NPSsubregions
% # appends the updated DAT to image_names_and_setup.mat
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('prep_4_apply_signatures_and_save','outputDir',htmlsavedir)
%
%
% *OUTPUT*
%
% DAT.SIG_conditions/DAT.SIG_contrasts contain data tables whose columns are
% conditions or contrasts, with variable names based on DAT.conditions or
% DAT.contrastnames (spaces replaced with underscores):
%
% DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
% DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
%
%
% *OPTIONS*
%
% NOTE:
%       defaults are specified in a2_set_default_options for any given model,
%       but if you want to add signature responses calculated with different
%       options (e.g. scaling), you can make a copy of this script with a
%       letter index (e.g. _s6a_) and change the default option here
%
% * myscaling_sigs             'raw' or 'scaled'
%
% * similarity_metric_sigs     'dotproduct', 'cosine_similarity', or 'correlation'
%
% * keyword_sigs                cell array of signature images and/or keywords passed into load_image_set; 'all' includes every available signature
%
%
% *NOTES*
%
% * NPS subregions are only computed if 'nps' is included in keyword_sigs
%   (or keyword_sigs = 'all') - this could/should later be expanded to
%   other signatures with subregion definitions
%
% -------------------------------------------------------------------------
%
% revamped by: Lukas Van Oudenhove
%
% date:   Leuven, January, 2023
%
% -------------------------------------------------------------------------
%
% prep_4_apply_signatures_and_save.m         v2.2
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
% myscaling_sigs = 'raw'/'scaled';
% similarity_metric_sigs = 'keyword';
% keyword_sigs = 'keyword';


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
    
required_fields = {'contrastnames', 'contrasts' 'contrastcolors', 'conditions', 'colors'};

ok_to_run = plugin_check_required_fields(DAT, required_fields); % Checks and prints warnings
if ~ok_to_run
    return
end


%% SELECTED SIGNATURES
% -------------------------------------------------------------------------

switch myscaling_sigs
    
    case 'raw'
        
        data_object_conds = DATA_OBJ;
        data_object_conts = DATA_OBJ_CON;
        
    case 'scaled'
        
        data_object_conds = DATA_OBJsc;
        data_object_conts = DATA_OBJ_CONsc;
        
end

for sig = 1:size(keyword_sigs,2)

    fprintf('\n\n');

    [~,signame] = fileparts(char(keyword_sigs{sig}));
    printhdr(['APPLYING SIGNATURE ', upper(signame), ' ON ', upper(myscaling_sigs), ' CONDITIONS AND CONTRASTS, SIMILARITY METRIC ', similarity_metric_sigs]);

    fprintf('\n\n');
    
    if contains(keyword_sigs{sig},filesep) % path to image rather than keyword

        DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(signame) = apply_all_signatures(data_object_conds, 'conditionnames', DAT.conditions, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs(sig));
        DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(signame) = apply_all_signatures(data_object_conts, 'conditionnames', DAT.contrastnames, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs(sig));
        
    else
        
        DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(signame) = apply_all_signatures(data_object_conds, 'conditionnames', DAT.conditions, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs{sig});
        DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(signame) = apply_all_signatures(data_object_conts, 'conditionnames', DAT.contrastnames, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs{sig});
        
    end

end


%% NPS SUBREGIONS
% -------------------------------------------------------------------------

if sum(contains(keyword_sigs,'nps')) > 0 || isequal(keyword_sigs{1},'all')

    nr_conds = size(DAT.conditions,2);

    % subregion names
    posnames = {'vermis' 'rIns' 'rV1' 'rThal' 'lIns' 'rdpIns' 'rS2_Op' 'dACC'};
    negnames = {'rLOC' 'lLOC' 'rpLOC' 'pgACC' 'lSTS' 'rIPL' 'PCC'};

    DAT.NPSsubregions.posnames = posnames;
    DAT.NPSsubregions.negnames = negnames;

    printhdr('Extracting NPS, adding to DAT')

    % CONDITIONS
    % ----------

    % NPS
    
    for i = 1:nr_conds

        switch similarity_metric_sigs
            
            case 'dotproduct'

                [DAT.npsresponse(i), ~, ~, DAT.NPSsubregions.npspos_by_region(i), DAT.NPSsubregions.npsneg_by_region(i)] = apply_nps(data_object_conds{i}, 'noverbose', 'notables');
                
            case 'cosine_similarity'

                [DAT.npsresponse(i), ~, ~, DAT.NPSsubregions.npspos_by_region(i), DAT.NPSsubregions.npsneg_by_region(i)] = apply_nps(data_object_conds{i}, 'noverbose', 'notables', similarity_metric_sigs);

        end

    end

    % NPS subregions
    
    printhdr('Extracting NPS Subregions, adding to DAT.NPSsubregions')

    clear posdat negdat spos sneg xx
    
    for i = 1:nr_conds

        % Get averages
        DAT.NPSsubregions.posdat{i} = nanmean(DAT.NPSsubregions.npspos_by_region{i})'; % mean across subjects
        DAT.NPSsubregions.stepos{i} = ste(DAT.NPSsubregions.npspos_by_region{i})'; % ste

        DAT.NPSsubregions.negdat{i} = nanmean(DAT.NPSsubregions.npsneg_by_region{i})'; % mean across subjects
        DAT.NPSsubregions.steneg{i} = ste(DAT.NPSsubregions.npsneg_by_region{i})'; % ste

    end

    
    % CONTRASTS
    % ---------

    printhdr('Defining NPS contrasts, adding to DAT')

    nr_conts = size(DAT.contrasts, 1);

    DAT.npscontrasts = {};

    for c = 1:nr_conts
        
        mycontrast = DAT.contrasts(c, :);
        wh = find(mycontrast);

        DAT.npscontrasts{c} = cat(2, DAT.npsresponse{wh}) * mycontrast(wh)';

        % subregions
        DAT.NPSsubregions.npspos_by_region_contrasts{c} = zeros(size(DAT.NPSsubregions.npspos_by_region{wh(1)}));
        DAT.NPSsubregions.npsneg_by_region_contrasts{c} = zeros(size(DAT.NPSsubregions.npsneg_by_region{wh(1)}));

        for j = 1:length(wh)

            DAT.NPSsubregions.npspos_by_region_contrasts{c} = DAT.NPSsubregions.npspos_by_region_contrasts{c} + DAT.NPSsubregions.npspos_by_region{wh(j)} * mycontrast(wh(j));
            DAT.NPSsubregions.npsneg_by_region_contrasts{c} = DAT.NPSsubregions.npsneg_by_region_contrasts{c} + DAT.NPSsubregions.npsneg_by_region{wh(j)} * mycontrast(wh(j));

        end
    end

    for i = 1:nr_conts

        % Get averages
        DAT.NPSsubregions.posdat_contrasts{i} = nanmean(DAT.NPSsubregions.npspos_by_region_contrasts{i})'; % mean across subjects
        DAT.NPSsubregions.stepos_contrasts{i} = ste(DAT.NPSsubregions.npspos_by_region_contrasts{i})'; % ste

        DAT.NPSsubregions.negdat_contrasts{i} = nanmean(DAT.NPSsubregions.npsneg_by_region_contrasts{i})'; % mean across subjects
        DAT.NPSsubregions.steneg_contrasts{i} = ste(DAT.NPSsubregions.npsneg_by_region_contrasts{i})'; % ste

    end

end


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING SIGNATURE RESPONSES TO DAT');
fprintf('\n\n');

cd(resultsdir); % unannex image_names_and_setup.mat file if already datalad saved to prevent write permission problems
! git annex unannex image_names_and_setup.mat
cd(rootdir);

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'DAT', '-append');

