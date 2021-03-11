%% LOAD SIGNATURE MAPS
% -------------------------------------------------------------------------
[npsplus, netnames, imgnames] = load_image_set('npsplus');
npsplus.image_names = netnames;


%% CONDITIONS
% -------------------------------------------------------------------------

% RIVERPLOT
% ---------

printhdr('Cosine Similarity : All conditions');

% Name the conditions
for i = 1:length(DATA_OBJ), DATA_OBJ{i}.image_names = DAT.conditions{i}; end

% plot significant associations only
riverplot(DATA_OBJ, 'layer2', npsplus, 'pos', 'significant_only', 'layer1colors', DAT.colors, 'layer2colors', seaborn_colors(length(netnames)));
hh=figure(1);
figtitle = 'CANlab signatures riverplot of conditions';
set(hh, 'Tag', figtitle);

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

plugin_save_figure; % @lukasvo76: figure does not save correctly, 


%% CONTRASTS
% -------------------------------------------------------------------------

if ~isfield(DAT, 'contrasts') || isempty(DAT.contrasts)
    % skip
    return
end

printhdr('Cosine Similarity : All contrasts');

k = size(DAT.contrasts, 1);

% RIVERPLOT
% ---------

for i = 1:length(DATA_OBJ_CON), DATA_OBJ_CON{i}.image_names = DAT.contrastnames{i}; end

% plot significant associations only
riverplot(DATA_OBJ_CON, 'layer2', npsplus, 'pos', 'significant_only', 'layer1colors', DAT.contrastcolors, 'layer2colors', seaborn_colors(length(netnames)));
hh=figure(2);
figtitle = 'CANlab signatures riverplot of contrasts';
set(hh, 'Tag', figtitle);

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

plugin_save_figure;