function [fh, actual_size] = plugin_set_figure_size(varargin)
% Plugin helper function for Second_level_analysis_template_scripts
% See those scripts for usage.
%
% Sizes the current figure (gcf) for capture into a publish() HTML report,
% choosing the largest canvas that the current session's display can
% actually capture. Replaces set(gcf,'WindowState','maximized').
%
% WHY NOT 'maximized': maximizing ties figure pixel dimensions to whichever
% client screen (e.g. an X2go session) happened to be connected when the
% script ran. All figures in these scripts use MATLAB's default font sizes,
% which are specified in POINTS (1 pt = 1/72 inch, a fixed physical unit).
% A canvas whose physical size changes per session therefore renders text
% inconsistently too large or too small depending on who ran the script.
%
% WHY INCHES, NOT PIXELS: anchoring the canvas in inches - a physical unit
% related to points by an exact, DPI-independent conversion - keeps the
% font-to-canvas ratio constant across sessions whose ScreenPixelsPerInch
% differs. A pixel-fixed canvas would not: point-based text would render at
% a different pixel size on each session while the canvas stayed put. Do
% not change the primary specification back to pixels.
%
% WHY IT NOW FITS TO THE SCREEN: publish() captures what is on screen. A
% figure larger than the display - or positioned partly off it - is
% captured at display size instead, AND at a different aspect ratio than
% requested, silently. get(fh,'Position') still reports the size you asked
% for, so nothing warns you. Measured on the LaBGAS server (1718x1360 at
% 133 DPI): the previous fixed 16x10 inch default needs 2128x1330 px, does
% not fit, and captured 1718x1254 - aspect 1.37 instead of 1.60, i.e. the
% same result as 'maximized', which is what this function exists to avoid.
% The failure appears on high-DPI sessions, which is precisely the case it
% was written for.
%
% WHY THE DEFAULT IS 12 x 7.5 INCHES, NOT 16 x 10: 16 x 10 in is not
% reachable on any lab laptop. It needs 16*DPI x 10*DPI pixels of window,
% so a 1366x768 client would have to run at 72 DPI and a 1600x900 client at
% 84 DPI - both impractically small to work in. 12 x 7.5 in (same 16:10
% aspect) is reachable on every screen in the lab at 96 DPI, and leaves
% headroom up to ~135 DPI on a 1920x1080 client and ~140 DPI on half of a
% 3440x1440 ultrawide. A default nobody can actually achieve guarantees the
% inconsistency this function exists to prevent. Run
% LaBGAScore_check_display (LaBGAScore/clean) to see what your own session
% can do, and pass 'width'/'height' explicitly if you want something else.
%
% So the requested size is now treated as an upper bound. If it does not
% fit the display, both dimensions are scaled down by the same factor, so
% the ASPECT RATIO IS ALWAYS HONOURED and only the absolute size gives way.
% The font-to-canvas ratio is then constant across every session whose
% display can hold the requested size, and degrades gracefully (text
% relatively larger, which keeps a smaller capture legible) below that.
% The figure is also repositioned fully on-screen, since a window hanging
% off the edge is clamped at capture no matter how it was sized.
%
% Call this immediately before drawnow/snapnow (and before
% plugin_save_figure, if used), after all plotting/display calls for the
% figure are complete.
%
% USAGE:
% plugin_set_figure_size()                         % default, fitted to screen
% plugin_set_figure_size('nrows', nrows)           % multi-row montage (canlab_results_fmridisplay 'multirow')
% plugin_set_figure_size('width', w, 'height', h)  % explicit upper bound
%
% OPTIONAL NAME-VALUE ARGUMENTS:
% 'width'   maximum figure width in inches (default 12)
% 'height'  maximum figure height in inches (default 7.5 if 'nrows' not
%           given; keeps the 16:10 aspect of the previous 16x10 default)
% 'nrows'   number of montage rows passed to canlab_results_fmridisplay's
%           'multirow' option (e.g. num_effects). Does NOT scale height:
%           canlab_results_fmridisplay's own 'multirow' code allocates each
%           row a FIXED normalized-coordinate band and sizes the whole
%           figure the same way regardless of how many rows it holds (1-4
%           per figure, extra rows spill into new figures) - so a canvas
%           shrunk for fewer rows starves every row of physical space
%           rather than saving any. Confirmed via real-project testing: an
%           earlier version that shrank height for low nrows produced
%           montages with the title clipped and slices squeezed into a
%           sliver at the top. 'nrows' is accepted (so multirow call sites
%           can still document their row count) but currently only uses the
%           same default height as any other figure.
% 'margin'  [horizontal vertical] fraction of the screen to leave free for
%           window decorations and panels (default [0.02 0.06])
% 'minsize' warn if fitting forces the canvas below this width in inches
%           (default 7). A capture much smaller than this makes montage
%           text hard to read in the report.
% 'verbose' print a line when the requested size had to be reduced
%           (default true). Reported only ONCE per session per distinct
%           request/display combination, since a published report calls
%           this once per figure. Silence it entirely with false.
%
% OUTPUT:
% fh            handle of the resized figure
% actual_size   [width height] in inches actually applied
%
% Example:
% figure; plot(1:10);
% plugin_set_figure_size();
% drawnow, snapnow;
%
% See also: LaBGAScore_prov_publish (LaBGAScore/clean), which records the
% session's screen size and DPI in the published report, so figure
% differences between machines are diagnosable after the fact.

p = inputParser;
addParameter(p, 'width', 12, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'height', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'nrows', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
addParameter(p, 'margin', [0.02 0.06], @(x) isnumeric(x) && numel(x) == 2 && all(x >= 0 & x < 0.5));
addParameter(p, 'minsize', 7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, varargin{:});

width = p.Results.width;
height = p.Results.height;
margin = p.Results.margin;

if isempty(height)
    height = 7.5;                      % keeps the 16:10 aspect of the old 16x10 default
end

fh = gcf;

% WindowState 'maximized' silently overrides any Position set while it is
% active, so it must be cleared first.
set(fh, 'WindowState', 'normal');


%% FIT THE REQUEST TO WHAT THIS DISPLAY CAN CAPTURE
% -------------------------------------------------------------------------

screen_px = get(0, 'ScreenSize');           % [1 1 width height], pixels
dpi = get(0, 'ScreenPixelsPerInch');

usable_w_in = screen_px(3) * (1 - margin(1)) / dpi;
usable_h_in = screen_px(4) * (1 - margin(2)) / dpi;

% one scale factor for both dimensions, so the aspect ratio survives
scale = min([1, usable_w_in / width, usable_h_in / height]);

actual_size = [width height] * scale;

% Report at most once per session per distinct situation. A published
% report calls this once per figure - often eight or more times - and the
% same notice repeated down the page is noise, not information.
persistent announced
if isempty(announced), announced = {}; end

situation = sprintf('%g_%g_%d_%d_%g', width, height, screen_px(3), screen_px(4), dpi);
firsttime = ~ismember(situation, announced);

if firsttime
    announced{end+1} = situation;
end

if p.Results.verbose && scale < 1 && firsttime
    fprintf(['plugin_set_figure_size: %.3g x %.3g in does not fit this display ' ...
             '(%d x %d px at %g DPI); using %.3g x %.3g in instead, aspect ratio ' ...
             'preserved. Reported once per session.\n'], width, height, ...
             screen_px(3), screen_px(4), dpi, actual_size(1), actual_size(2));
end

if p.Results.verbose && actual_size(1) < p.Results.minsize && firsttime
    warning('plugin_set_figure_size:smallCanvas', ...
        ['this display only allows a %.3g in wide figure, below the %.3g in ' ...
         'guideline; montage text may be hard to read in the published report. ' ...
         'Consider running from a session with a larger or lower-DPI display.'], ...
        actual_size(1), p.Results.minsize);
end


%% APPLY, KEEPING THE WHOLE WINDOW ON SCREEN
% -------------------------------------------------------------------------
% A window that extends past the screen edge is clamped at capture just as a
% too-large one is, so position matters as much as size. Anchor near the
% bottom-left, which is safe for any size that fits.

set(fh, 'Units', 'inches');

left = margin(1) * screen_px(3) / dpi / 2;
bottom = margin(2) * screen_px(4) / dpi / 2;

set(fh, 'Position', [left, bottom, actual_size(1), actual_size(2)]);

end % function
