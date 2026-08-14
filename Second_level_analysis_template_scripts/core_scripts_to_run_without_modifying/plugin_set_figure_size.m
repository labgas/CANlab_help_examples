function fh = plugin_set_figure_size(varargin)
% Plugin helper function for Second_level_analysis_template_scripts
% See those scripts for usage.
%
% Sets a fixed, screen-resolution- AND screen-DPI-independent Position for
% the current figure (gcf), for consistent publish() HTML output regardless
% of which client (e.g. X2go) machine/screen was connected to the server
% when the script was run. Replaces set(gcf,'WindowState','maximized'),
% which ties figure pixel dimensions to the active client's screen and
% therefore makes fixed-point-size text (the default for all figures in
% these scripts, since none set FontSize/FontUnits explicitly) look
% inconsistently too large or too small depending on who ran the script.
%
% Figure size is set in INCHES, not pixels. This matters: MATLAB's default
% font sizes are specified in points (1 point = 1/72 inch, a fixed physical
% unit, independent of screen DPI). If the figure canvas were instead fixed
% in pixels, the font-to-canvas ratio would still drift across sessions
% whose ScreenPixelsPerInch differs (which X2go sessions can, depending on
% the connecting client's display), because point-based text would render
% at a different pixel size on each session while the pixel-fixed canvas
% would not. Anchoring the canvas in inches - another physical unit, related
% to points by an exact, DPI-independent conversion - keeps the font-to-
% canvas ratio constant no matter what DPI the session reports. Do not
% change this back to pixels.
%
% Call this immediately before drawnow/snapnow (and before plugin_save_figure,
% if used), after all plotting/display calls for the figure are complete.
%
% USAGE:
% plugin_set_figure_size()                        % default fixed size
% plugin_set_figure_size('nrows', nrows)           % multi-row montage (canlab_results_fmridisplay 'multirow')
% plugin_set_figure_size('width', w, 'height', h)  % fully explicit override
%
% OPTIONAL NAME-VALUE ARGUMENTS:
% 'width'   figure width in inches (default 16)
% 'height'  figure height in inches (default 10 if 'nrows' not given)
% 'nrows'   number of montage rows passed to canlab_results_fmridisplay's
%           'multirow' option (e.g. num_effects). Does NOT scale height:
%           canlab_results_fmridisplay's own 'multirow' code (see its
%           'multirow' case) allocates each row a FIXED normalized-
%           coordinate band and sizes the whole figure the same way
%           regardless of how many rows it holds (1-4 per figure, extra
%           rows spill into new figures) - so a canvas shrunk for fewer
%           rows starves every row of physical space rather than saving
%           any. Confirmed via real-project testing: an earlier version
%           of this function that shrank height for low nrows produced
%           montages with the title clipped and slices squeezed into a
%           sliver at the top. 'nrows' is accepted (so multirow call
%           sites can still document their row count) but currently only
%           uses the same default height as any other figure.
%
% OUTPUT:
% fh   handle of the resized figure
%
% Example:
% figure; plot(1:10);
% plugin_set_figure_size();
% drawnow, snapnow;

p = inputParser;
addParameter(p, 'width', 16, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'height', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'nrows', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
parse(p, varargin{:});

width = p.Results.width;
height = p.Results.height;

if isempty(height)
    height = 10;                       % default figure height; not varied by 'nrows' - see note above
end

fh = gcf;

% Defensive: WindowState 'maximized' silently overrides any Position set
% while it is active, so it must be cleared first.
set(fh, 'WindowState', 'normal');
set(fh, 'Units', 'inches');

current_pos = get(fh, 'Position');
set(fh, 'Position', [current_pos(1), current_pos(2), width, height]);

end % function
