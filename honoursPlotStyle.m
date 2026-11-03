% HONOURSPLOTSTYLE
function col = honoursPlotStyle()

% Have one consistent plotting style across the entire honours project.

% Global figure defaults

set(groot, 'defaultFigureColor', 'w');

set(groot, 'defaultAxesFontName', 'Helvetica');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultAxesLineWidth', 1.1);

set(groot, 'defaultLineLineWidth', 2.7);

set(groot, 'defaultAxesBox', 'off');

set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultAxesGridAlpha', 0.15);
set(groot, 'defaultAxesMinorGridAlpha', 0.08);

set(groot, 'defaultLegendBox', 'off');

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% CHOSEN COLOUR PALETTE: NAVY / TEAL / CORAL

col.data  = [0.10 0.11 0.13];

% Tumour variables
col.total = [0.20 0.27 0.43];     % navy
col.A     = [0.18 0.53 0.63];     % teal-blue
col.B     = [0.85 0.39 0.36];     % coral

% Immune variables
col.T     = [0.12 0.58 0.52];     % teal
col.D     = [0.91 0.61 0.25];     % amber

% Treatments
col.no    = [0.47 0.49 0.52];     % cool grey
col.ivm   = [0.91 0.55 0.24];     % amber-orange
col.pd1   = [0.18 0.43 0.66];     % navy blue
col.combo = [0.58 0.31 0.55];     % muted magenta

end