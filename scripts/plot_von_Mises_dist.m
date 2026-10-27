%% plot_von_Mises_dist.m
% Plots von Mises distributions
%   f(theta) = exp(beta*cos(theta)) / (2*pi*I0(beta))
% for several concentration parameters beta, in two versions:
%   1. Color             (von_mises_distributions.svg/png)
%   2. Black and white   (von_mises_distributions_bw.svg/png)
% Formatted to match the advection test plots.
%
% Brzenski, Jared
% 06/26/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clearvars;

svgfilename    = 'von_mises_distributions.svg';
filename       = 'von_mises_distributions.png';
svgfilename_bw = 'von_mises_distributions_bw.svg';
filename_bw    = 'von_mises_distributions_bw.png';

% Single Column figure size (pixels): half the width of the two-column
% advection figures, so the fonts print at the same size
imagesizevector = [0 0 821 616];

%% Settings
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter',        'latex');
set(groot, ...
    'defaultAxesFontSize', 16, ...
    'defaultTextFontSize', 16, ...
    'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial', ...
    'defaultLineLineWidth', 2.0);

theta  = linspace(-pi, pi, 500);
kappas = [0, 0.5, 1, 2, 5, 10];   % concentration parameter beta

% Black and white line styles, one per beta
bw_lstyles = {'-', '--', ':', '-.', '-', '--'};

titlefontsize  = 14;
labelfontsize  = 16;
legendfontsize = 12;

legendStrings = arrayfun(@(k) sprintf('$\\beta = %g$', k), ...
                         kappas, 'UniformOutput', false);

%% Von Mises densities, one row per beta
F = zeros(numel(kappas), numel(theta));
for i = 1:numel(kappas)
    k = kappas(i);
    F(i,:) = exp(k*cos(theta)) / (2*pi*besseli(0,k));
end

%% =========================================================
%  Figure 1: Color
% ==========================================================
fig1 = figure('Name', 'Von Mises Distributions', ...
              'Units', 'pixels', 'Position', imagesizevector);
ax = axes(fig1);
hold(ax, 'on');

h = gobjects(numel(kappas), 1);   % store handles
for i = 1:numel(kappas)
    h(i) = plot(ax, theta, F(i,:));
end
hold(ax, 'off');

set(ax, 'FontSize', labelfontsize, 'TickLabelInterpreter', 'latex');
xticks(ax, [-pi -pi/2 0 pi/2 pi]);
xticklabels(ax, {'$-\pi$', '$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$', '$\pi$'});
xlim(ax, [-pi pi]);

xlabel(ax, '$\theta$',    'Interpreter', 'latex', 'FontSize', labelfontsize);
ylabel(ax, '$f(\theta)$', 'Interpreter', 'latex', 'FontSize', labelfontsize);
%title(ax, 'Von Mises Distributions', 'Interpreter', 'latex', 'FontSize', titlefontsize);

legend(ax, h, legendStrings, 'Interpreter', 'latex', ...
       'Location', 'northwest', 'FontSize', legendfontsize);
grid(ax, 'on');
box(ax, 'on');

print(fig1, svgfilename, '-dsvg');
exportgraphics(fig1, filename, 'Resolution', 300);

%% =========================================================
%  Figure 2: Black and white
% ==========================================================
fig2 = figure('Name', 'Von Mises Distributions (B&W)', ...
              'Units', 'pixels', 'Position', imagesizevector, 'Color', 'w');
ax = axes(fig2);
hold(ax, 'on');

h = gobjects(numel(kappas), 1);
for i = 1:numel(kappas)
    h(i) = plot(ax, theta, F(i,:), 'Color', 'k', 'LineStyle', bw_lstyles{i});
end
hold(ax, 'off');

set(ax, 'FontSize', labelfontsize, 'TickLabelInterpreter', 'latex');
xticks(ax, [-pi -pi/2 0 pi/2 pi]);
xticklabels(ax, {'$-\pi$', '$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$', '$\pi$'});
xlim(ax, [-pi pi]);

xlabel(ax, '$\theta$',    'Interpreter', 'latex', 'FontSize', labelfontsize);
ylabel(ax, '$f(\theta)$', 'Interpreter', 'latex', 'FontSize', labelfontsize);
%title(ax, 'Von Mises Distributions', 'Interpreter', 'latex', 'FontSize', titlefontsize);

legend(ax, h, legendStrings, 'Interpreter', 'latex', ...
       'Location', 'northwest', 'FontSize', legendfontsize, 'Box', 'off');
grid(ax, 'on');
box(ax, 'on');

print(fig2, svgfilename_bw, '-dsvg');
exportgraphics(fig2, filename_bw, 'Resolution', 300);
