clear;
close all;
clc;

%% ========================================================================
% General space setup and variable declarations
auburn_color_map = Generate_Auburn_Colormap_v0();

% --- constitutive parameters ---
eps_0 = 8.854e-12; 
eps_r1 = 2; eps_r2 = 3; % relative permittivities of layers
mu_0 = (4e-7)*pi;
% -------------------------------

c = 1 / sqrt(mu_0 * eps_0); % speed of light

f0 = 20e6;   % design frequency: 20 MHz
lambda0 = c / f0; % free-space wavelength (~14.99 m)

% --- make the simulation space ---
dx = 0.025; % dx = 2.5cm
x_max = 150; % 150m
x = 0:dx:x_max;
Nx = length(x);

% Courant stability factor (starting with magic time step)
S = 1.0; 
dt = S * dx / c; % ~83.39 ps
t_max = 450e-9;  % 450 ns
t = 0:dt:t_max;
Nt = length(t);

E = zeros(Nx, 3); % pre-allocate E array
eps_r = ones(1, Nx);  % pre-allocate permittivity array
% ---------------------------------


% --- build anti-reflective coating ---
% quarter-wavelength layer thicknesses
d1 = lambda0 / (4 * sqrt(eps_r1)); % ~2.6498 m
d2 = lambda0 / (4 * sqrt(eps_r2)); % ~2.1636 m
d3 = d1;                           % ~2.6498 m

% positions in meters (roughly centered in the space to prevent PEC 
% boundary reflections before t_max)
x_src   = 50.0;  % source excitation position
x_obs   = 70.0;  % reflection observation point
x_start = 75.0;  % front face of coating

% convert meters to indices
idx_src = round(x_src / dx) + 1;
idx_obs = round(x_obs / dx) + 1;

% map permittivity layers onto spatial grid
% this method allows me to change the thicknesses and have the grid update
idx_L1 = (x >= x_start) & (x < x_start + d1);
idx_L2 = (x >= x_start + d1) & (x < x_start + d1 + d2);
idx_L3 = (x >= x_start + d1 + d2) & (x < x_start + d1 + d2 + d3);

eps_r(idx_L1) = eps_r1;
eps_r(idx_L2) = eps_r2;
eps_r(idx_L3) = eps_r1;
% -------------------------------------

%% ========================================================================
% Figure setup

% set up figure to animate during runtime
fig = figure('Color', 'w', 'Theme','Light');

% RGB colors directly from the Auburn colormap function
navy_bg   = auburn_color_map(1, :);     
white_eps = auburn_color_map(50, :);    
orange_wf = auburn_color_map(end, :);
gray_eps  = [0.65, 0.65, 0.65]; % gray for eps_r = 3

ax = axes('Parent', fig, 'Color', navy_bg);
hold(ax, 'on'); y_lim = [-1.5, 1.5];

% draw dielectric layer patches behind the waveform
p1 = patch(ax, [x_start, x_start+d1, x_start+d1, x_start], ...
    [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], white_eps, ...
    'EdgeColor', [0.65, 0.65, 0.65], 'LineWidth', 0.8, 'FaceAlpha', 0.9);
p2 = patch(ax, [x_start+d1, x_start+d1+d2, x_start+d1+d2, x_start+d1], ...
    [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], gray_eps, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.9);
p3 = patch(ax, [x_start+d1+d2, x_start+d1+d2+d3, x_start+d1+d2+d3, x_start+d1+d2], ...
    [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], white_eps, ...
    'EdgeColor', [0.65, 0.65, 0.65], 'FaceAlpha', 0.9);

% dummy patch for legend entry for the free-space background
p0 = patch(ax, [NaN NaN NaN NaN], [NaN NaN NaN NaN], navy_bg, ...
    'EdgeColor', 'k');

% initialize waveform line handle to update YData in loop
h_wave = plot(ax, x, E(:, 2), 'Color', orange_wf, 'LineWidth', 2.2);

% mark observation location (may or may not use this)
xline(ax, x_obs, ':', 'Obs', 'Color', [0.85 0.85 0.85], ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% formatting axes and labeling figure
ylim(ax, y_lim);
xlim(ax, [0, x_max]);
grid(ax, 'on');
ax.GridColor = [0.4, 0.45, 0.55];
ax.GridAlpha = 0.35;
xlabel(ax, 'Position x (m)', 'FontSize', 12);
ylabel(ax, 'Electric Field Intensity E_z (V/m)', 'FontSize', 12);

% title (says 20MHz initially)
title(ax, ...
    sprintf('1D FDTD Wave Propagation at %.1f MHz (Step 0 / %d, t = 0.0 ns)', ...
    20/1e6, Nt), 'FontSize', 13);

% legend
legend(ax, [h_wave, p0, p1, p2], ...
    {'E_z Waveform', 'Free Space (\epsilon_r = 1)', ...
     'ARC Layer 1 & 3 (\epsilon_r = 2)', 'ARC Layer 2 (\epsilon_r = 3)'}, ...
    'Location', 'southwest', 'FontSize', 10);

%% ========================================================================
% 1D FDTD LOOP

% factor multiplying the first three terms in the update equation
eps_r_col = eps_r(:); % make column vector to use in element-wise
factor = (1 ./ (mu_0 * eps_0 * eps_r_col(2:end-1))) * ((dt / dx)^2);

% set frequency for this run
f = f0;

tic % the whole loop takes ~2.3s to run
for ii = 1:Nt
    % source excitation 
    % (simply turning a sine wave on for 2 periods)
    if t(ii) < 2/f
        E(1,2) = sin(2*pi*f*t(ii));
    else
        E(1,2) = 0;
    end

    % forward stepping update equation [ E(i, n+1) ]
    E(2:end-1, 3) = factor.*(E(3:end, 2) - 2*E(2:end-1, 2) + E(1:end-2, 2)) ...
        + 2*E(2:end-1, 2) - E(2:end-1, 1);

    % enforce PEC boundary conditions at edges (E = 0)
    E(1, 3)   = 0;
    E(end, 3) = 0;

    % scoot array forward in time (n)
    E(:, 1) = E(:, 2); % E(i, n-1) = E(i, n)
    E(:, 2) = E(:, 3); % E(i, n)   = E(i, n+1)

    % update animation plot (every 50th timestep)
    if mod(ii, 10) == 0 || ii == Nt
        % update wave amplitude handle
        set(h_wave, 'YData', E(:, 2));

        % update title with accurate timestep and frequency
        title(ax, ...
            sprintf('1D FDTD Wave Propagation at %.1f MHz (Step %d / %d, t = %.1f ns)', ...
            f/1e6, ii, Nt, t(ii)*1e9), 'FontSize', 13);
        drawnow;
    end
end
toc
% -------------------------------

%% ========================================================================
% Comparison of Analytical and Numerical Solutions

% calculate percent difference relative to analytical baseline
% percent_diff = abs(V_out_numerical - ...
%                         V_out_analytical) / abs(V_out_analytical) * 100;

% display comparison results in command window
% fprintf('\n=======================================================\n');
% fprintf('   ELECTRIC FIELD MILL SOLUTION COMPARISON RESULTS\n');
% fprintf('=======================================================\n');
% fprintf('Analytical Output Voltage  : %10.4f uV\n', V_out_analytical * 1e6);
% fprintf('Numerical Output Voltage   : %10.4f uV (at steady-state)\n', V_out_numerical * 1e6);
% fprintf('Absolute Difference        : %10.4f uV\n', abs(V_out_numerical - V_out_analytical) * 1e6);
% fprintf('Percent Difference         : %10.2f %%\n', percent_diff);
% fprintf('=======================================================\n\n');

%% ========================================================================
% Plotting of results


%% ========================================================================
% Functions

% --- auburn colormap generation ---
% function from TRACE repo:

function auburn_color_map = Generate_Auburn_Colormap_v0()
% This function generates the Auburn colors colormap where
% blue is low value and orange is high value. RGB values
% according to 
% https://ocm.auburn.edu/brand-center/_assets/pdf/au-colorpalette-primary-supporting.pdf

% INPUTS: None.

% OUTPUTS: 100 by 3 array of Auburn blue to orange

% DEPENDENCIES: none

% Author: Clint Snider
% last edited: 02/05/2025

% DESIRED UPDATES: Maybe make length 256.

% RGB values from [0:255]
orange = [232,97,0];
blue   = [11,35,65];
white  = [255,255,255];

%keep as an even number
length_colormap = 100;

% make blue to white
auburn_color_map = zeros(length_colormap,3);
auburn_color_map(1:length_colormap/2,1) = (interp1([1,length_colormap/2], ...
    [blue(1) white(1)], 1:length_colormap/2))';
auburn_color_map(1:length_colormap/2,2) = (interp1([1,length_colormap/2], ...
    [blue(2) white(2)], 1:length_colormap/2))';
auburn_color_map(1:length_colormap/2,3) = (interp1([1,length_colormap/2], ...
    [blue(3) white(3)], 1:length_colormap/2))';

% make white to orange
auburn_color_map(length_colormap/2+1:length_colormap,1) = ...
    (interp1([length_colormap/2+1,length_colormap], [white(1) orange(1)], ...
    length_colormap/2+1:length_colormap))';
auburn_color_map(length_colormap/2+1:length_colormap,2) = ...
    (interp1([length_colormap/2+1,length_colormap], [white(2) orange(2)], ...
    length_colormap/2+1:length_colormap))';
auburn_color_map(length_colormap/2+1:length_colormap,3) = ...
    (interp1([length_colormap/2+1,length_colormap], [white(3) orange(3)], ...
    length_colormap/2+1:length_colormap))';

% RGB values normalized from [0:1]
auburn_color_map = auburn_color_map/255;
end

% ----------------------------------