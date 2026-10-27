function run_all(what)
%% run_all.m
% Regenerates the .mat data files and the figures for the semi-Lagrangian
% fish advection tests.
%
%   run_all            generate every .mat file, then every figure
%   run_all('data')    only the data generation scripts
%   run_all('figures') only the figure scripts (needs the .mat files)
%
% Each script runs from the scripts/ directory, in its own workspace, so
% the .mat files and figures are written to scripts/. Before every script
% the random number generator is reset to its start-up state, and after
% every script the groot defaults are restored, so each script gives the
% same result as when it is run on its own in a fresh MATLAB session.
%
% Brzenski, Jared
% 10/27/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1, what = 'all'; end
what = validatestring(what, {'all', 'data', 'figures'});

%% Scripts, in run order
data_scripts = { ...
    'test1_Linear_A', 'test1_Linear_RA', 'test1_Linear_A_D', ...                    % horizontal, one file per CFL
    'test2_Diag_A',   'test2_Diag_RA',   'test2_Diag_AD', ...                       % diagonal,   one file per CFL
    'test3_CFL_A_error',      'test3_CFL_RA_error',      'test3_CFL_AD_error', ...  % horizontal speed x dt sweeps
    'test5_diag_CFL_A_error', 'test5_diag_CFL_RA_error', 'test5_diag_CFL_AD_error', ... % diagonal speed x dt sweeps
    'test4_kappa_AD'};                                                              % kappa sweep at CFL = 1

figure_scripts = { ...
    'horizontal_advect_plot', 'diagonal_advect_plot', ...
    'analyze_test1_plot2',    'analyze_test2_plot2', ...
    'error_bar_plot',         'cfl_contour_plot', ...
    'analyze_mass_cons',      'cfl_shift_plot', ...
    'plot_von_Mises_dist'};

switch what
    case 'all',     scripts = [data_scripts, figure_scripts];
    case 'data',    scripts = data_scripts;
    case 'figures', scripts = figure_scripts;
end

%% Run
root     = fileparts(mfilename('fullpath'));
out_dir  = fullfile(root, 'scripts');
old_dir  = cd(out_dir);
defaults = get(groot, 'default');
cleanup  = onCleanup(@() restore_state(old_dir, defaults));

n       = numel(scripts);
elapsed = zeros(1, n);
for k = 1:n
    fprintf('\n[%2d/%d] %s\n', k, n, scripts{k});
    t = tic;
    rng('default');   % start-up random state, so RA results do not depend on run order
    run_script(scripts{k});
    elapsed(k) = toc(t);
    reset_defaults(defaults);
end

% The scripts clear the command window, so print the log once at the end
fprintf('\n%-26s %10s\n', 'Script', 'Time [s]');
for k = 1:n
    fprintf('%-26s %10.1f\n', scripts{k}, elapsed(k));
end
fprintf('Finished %d scripts in %.1f min. Output is in %s\n', ...
        n, sum(elapsed) / 60, out_dir);
end

%% =========================================================
%  Local functions
%% =========================================================
function run_script(name)
% Run a script in this function's workspace, so its clearvars cannot
% clear run_all's variables.
    eval(name);
end

function reset_defaults(defaults)
% Remove every groot default a script set, then put back the originals.
    names = fieldnames(get(groot, 'default'));
    for i = 1:numel(names)
        set(groot, names{i}, 'remove');
    end
    names = fieldnames(defaults);
    for i = 1:numel(names)
        set(groot, names{i}, defaults.(names{i}));
    end
end

function restore_state(old_dir, defaults)
    cd(old_dir);
    reset_defaults(defaults);
end
