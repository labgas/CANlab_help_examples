% Format: The prep_4 script extracts signature responses and saves them.
% These fields contain data tables:
% DAT.SIGNATURES.(signaturename).(data scaling).(similarity metric).by_condition
% DAT.SIGNATURES.(signaturename).(data scaling).(similarity metric).contrasts
%
% signaturenames is any of those from load_image_set('npsplus')
% (data scaling) is 'raw' or 'scaled', using DATA_OBJ or DATA_OBJsc
% (similarity metric) is 'dotproduct' or 'cosine_sim'
% 
% Each of by_condition and contrasts contains a data table whose columns
% are conditions or contrasts, with variable names based on DAT.conditions
% or DAT.contrastnames, but with spaces replaced with underscores.


%% USER OPTIONS
% -------------------------------------------------------------------------

% This is a standard block of code that can be used in multiple scripts.
% Each script will have its own options needed and default values for
% these.
% The code: 
% (1) Checks whether the option variables exist
% (2) Runs a2_set_default_options if any are missing
% (3) Checks again and uses the default options if they are still missing
% (e.g., not specified in an older/incomplete copy of a2_set_default_options)

% Now set in a2_set_default_options
options_needed = {'use_scaled_images'};  % Options we are looking for. Set in a2_set_default_options
options_exist = cellfun(@exist, options_needed); 

option_default_values = {false};          % defaults if we cannot find info in a2_set_default_options at all 

plugin_get_options_for_analysis_script


%% NPS RESPONSE AND SUBREGIONS
% -------------------------------------------------------------------------

k = length(DAT.conditions);

% subregion names
posnames = {'vermis'    'rIns'    'rV1'    'rThal'    'lIns'    'rdpIns'    'rS2_Op'    'dACC'};
negnames = {'rLOC'    'lLOC'    'rpLOC'    'pgACC'    'lSTS'    'rIPL'    'PCC'};

DAT.NPSsubregions.posnames = posnames;
DAT.NPSsubregions.negnames = negnames;

printhdr('Extracting NPS, adding to DAT')

% CONDITIONS
% ----------

% NPS, dot product and cosine sim
for i = 1:k
    
    if use_scaled_images
        
        [DAT.npsresponse(i), ~, ~, DAT.NPSsubregions.npspos_by_region(i), DAT.NPSsubregions.npsneg_by_region(i)] = apply_nps(DATA_OBJsc{i}, 'noverbose', 'notables');
        
        [DAT.npsresponse_cosinesim(i), ~, ~, DAT.NPSsubregions.npspos_by_region_cosinesim(i), DAT.NPSsubregions.npsneg_by_region_cosinesim(i)] = apply_nps(DATA_OBJsc{i}, 'noverbose', 'notables', 'cosine_similarity');

    else
        
        [DAT.npsresponse(i), ~, ~, DAT.NPSsubregions.npspos_by_region(i), DAT.NPSsubregions.npsneg_by_region(i)] = apply_nps(DATA_OBJ{i}, 'noverbose', 'notables');
        
        [DAT.npsresponse_cosinesim(i), ~, ~, DAT.NPSsubregions.npspos_by_region_cosinesim(i), DAT.NPSsubregions.npsneg_by_region_cosinesim(i)] = apply_nps(DATA_OBJ{i}, 'noverbose', 'notables', 'cosine_similarity');
        
    end
    
end

% , ~, ~, DAT.NPSsubregions.npspos_by_region_cosinesim(i), DAT.NPSsubregions.npsneg_by_region_cosinesim(i)

% NPS subregions, dot product and cosine sim
printhdr('Extracting NPS Subregions, adding to DAT.NPSsubregions')

clear posdat negdat spos sneg xx
for i = 1:k
    
    % Get averages
    DAT.NPSsubregions.posdat{i} = nanmean(DAT.NPSsubregions.npspos_by_region{i})'; % mean across subjects
    DAT.NPSsubregions.stepos{i} = ste(DAT.NPSsubregions.npspos_by_region{i})'; % ste
    
    DAT.NPSsubregions.negdat{i} = nanmean(DAT.NPSsubregions.npsneg_by_region{i})'; % mean across subjects
    DAT.NPSsubregions.steneg{i} = ste(DAT.NPSsubregions.npsneg_by_region{i})'; % ste
        
end

% CONTRASTS
% ---------

printhdr('Defining NPS contrasts, adding to DAT')

kc = size(DAT.contrasts, 1);

DAT.npscontrasts = {};

for c = 1:kc
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

for i = 1:kc
    
    % Get averages
    DAT.NPSsubregions.posdat_contrasts{i} = nanmean(DAT.NPSsubregions.npspos_by_region_contrasts{i})'; % mean across subjects
    DAT.NPSsubregions.stepos_contrasts{i} = ste(DAT.NPSsubregions.npspos_by_region_contrasts{i})'; % ste
    
    DAT.NPSsubregions.negdat_contrasts{i} = nanmean(DAT.NPSsubregions.npsneg_by_region_contrasts{i})'; % mean across subjects
    DAT.NPSsubregions.steneg_contrasts{i} = ste(DAT.NPSsubregions.npsneg_by_region_contrasts{i})'; % ste
    
end


%% ALL SIGNATURES
% -------------------------------------------------------------------------

printhdr('Extracting all signatures');

% CONDITIONS
% ----------

% RAW CONDITION IMAGES

% Dot product metric
DAT.SIG_conditions.raw.dotproduct = apply_all_signatures(DATA_OBJ, 'conditionnames', DAT.conditions);

% Cosine similarity
DAT.SIG_conditions.raw.cosine_sim = apply_all_signatures(DATA_OBJ, 'conditionnames', DAT.conditions, 'similarity_metric', 'cosine_similarity');

% SCALED CONDITION IMAGES  
% apply_all_signatures will do scaling as well, but we did this in image
% loading, so use those here

% Dot product metric
DAT.SIG_conditions.scaled.dotproduct = apply_all_signatures(DATA_OBJsc, 'conditionnames', DAT.conditions);

% Cosine similarity
DAT.SIG_conditions.scaled.cosine_sim = apply_all_signatures(DATA_OBJsc, 'conditionnames', DAT.conditions, 'similarity_metric', 'cosine_similarity');

% CONTRASTS
% ---------

if exist('DATA_OBJ_CON', 'var') && iscell(DATA_OBJ_CON) && ~isempty(DATA_OBJ_CON{1})

% RAW CONTRAST IMAGES
    DAT.SIG_contrasts.raw.dotproduct = apply_all_signatures(DATA_OBJ_CON, 'conditionnames', DAT.contrastnames);
    DAT.SIG_contrasts.raw.cosine_sim = apply_all_signatures(DATA_OBJ_CON, 'conditionnames', DAT.contrastnames, 'similarity_metric', 'cosine_similarity');

% SCALED CONTRAST IMAGES   
    DAT.SIG_contrasts.scaled.dotproduct = apply_all_signatures(DATA_OBJ_CONsc, 'conditionnames', DAT.contrastnames);
    DAT.SIG_contrasts.scaled.cosine_sim = apply_all_signatures(DATA_OBJ_CONsc, 'conditionnames', DAT.contrastnames, 'similarity_metric', 'cosine_similarity');
    
end


%% SAVE
% -------------------------------------------------------------------------

printhdr('Save results');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'DAT', '-append');