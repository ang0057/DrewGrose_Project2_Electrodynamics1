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
x_max = 200; % 200m
x = 0:dx:x_max;
Nx = length(x);

% Courant stability factor (starting with magic time step)
S = 1.0; 
dt = S * dx / c; % ~83.39 ps
t_max = 600e-9;  % 600 ns
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
x_obs   = 70.0;  % reflection observation point
x_start = 75.0;  % front face of coating

% convert meters to indices
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
% Animated figure setup

% set up figure to animate during runtime
fig1 = figure('Color', 'w', 'Theme','Light');

% RGB colors directly from the Auburn colormap function
navy_au   = auburn_color_map(1, :);     
white_eps = auburn_color_map(50, :);    
orange_au = auburn_color_map(end, :);
gray_eps  = [0.65, 0.65, 0.65]; % gray for eps_r = 3

ax1 = axes('Parent', fig1, 'Color', navy_au);
hold(ax1, 'on'); y_lim = [-1.5, 1.5];

% draw dielectric layer patches behind the waveform
p1 = patch(ax1, [x_start, x_start+d1, x_start+d1, x_start], ...
    [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], white_eps, ...
    'EdgeColor', [0.65, 0.65, 0.65], 'LineWidth', 0.8, 'FaceAlpha', 0.9);
p2 = patch(ax1, [x_start+d1, x_start+d1+d2, x_start+d1+d2, x_start+d1], ...
    [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], gray_eps, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.9);
p3 = patch(ax1, [x_start+d1+d2, x_start+d1+d2+d3, x_start+d1+d2+d3, x_start+d1+d2], ...
    [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], white_eps, ...
    'EdgeColor', [0.65, 0.65, 0.65], 'FaceAlpha', 0.9);

% dummy patch for legend entry for the free-space background
p0 = patch(ax1, [NaN NaN NaN NaN], [NaN NaN NaN NaN], navy_au, ...
    'EdgeColor', 'k');

% initialize waveform line handle to update YData in loop
h_wave = plot(ax1, x, E(:, 2), 'Color', orange_au, 'LineWidth', 2.2);

% mark observation location (may or may not use this)
xline(ax1, x_obs, ':', 'Obs', 'Color', [0.85 0.85 0.85], ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

% formatting axes and labeling figure
ylim(ax1, y_lim);
xlim(ax1, [0, x_max]);
grid(ax1, 'on');
ax1.GridColor = [0.4, 0.45, 0.55];
ax1.GridAlpha = 0.35;
xlabel(ax1, 'Position x (m)', 'FontSize', 12);
ylabel(ax1, 'Electric Field Intensity E_z (V/m)', 'FontSize', 12);

% title (says 20MHz initially)
title(ax1, ...
    sprintf('1D FDTD Wave Propagation at %.1f MHz (Step 0 / %d, t = 0.0 ns)', ...
    f0/1e6, Nt), 'FontSize', 13);

% legend
legend(ax1, [h_wave, p0, p1, p2], ...
    {'E_z Waveform', 'Free Space (\epsilon_r = 1)', ...
     'ARC Layer 1 & 3 (\epsilon_r = 2)', 'ARC Layer 2 (\epsilon_r = 3)'}, ...
    'Location', 'southwest', 'FontSize', 10);

%% ========================================================================
% 1D FDTD LOOP

% factor multiplying the first three terms in the update equation
eps_r_col = eps_r(:); % make column vector to use in element-wise
factor_coat = (1 ./ (mu_0 * eps_0 * eps_r_col(2:end-1))) * ((dt / dx)^2);
factor_free = (1 ./ (mu_0 * eps_0 * ones(Nx-2, 1)))     * ((dt / dx)^2);

% array for frequency sweep
rel_freqs_fdtd = 0.1:0.1:2.0; % 20 points to build a curve
freqs = rel_freqs_fdtd * f0;
Nf = length(freqs);

% pre-allocate vector for the computed reflection coefficients
Gamma_fdtd = zeros(1, Nf);

% last 200ns as the observation window for steady-state reflection analysis
idx_win = t >= (t_max - 200e-9);
t_win   = t(idx_win);

% quick anonymous function for 4 term blackman harris window
bh_win = @(N) (0.35875 ...
    - 0.48829 * cos(2*pi*(0:N-1)' / max(N-1, 1)) ...
    + 0.14128 * cos(4*pi*(0:N-1)' / max(N-1, 1)) ...
    - 0.01168 * cos(6*pi*(0:N-1)' / max(N-1, 1))).';

% info message for starting the big loop
fprintf('Running 1D FDTD frequency sweep (%d frequencies)...\n', Nf);

tic % the whole loop takes ~30s to run
for ff = 1:Nf
% set frequency for this run
f = freqs(ff);

% flag to animate only at 20 MHz
animate_this_run = (abs(f - f0) < 1);

% --- source excitation: blackman harris windowed sine wave ---
% [src_wave, ~] = Blackman_Harris_Windowed_Sin_Wave_v0(f, dt, t, 0, 1);
% I was having issues with the function from the TRACE repo because of the
% time resolution in this FDTD, so I built the waveform from scratch
T_period = 1 / f;
% using 6 cycles for mid/high frequencies, or 250 ns for VLF
T_burst = min(6 * T_period, 250e-9);
n_burst = t <= T_burst;
N_b     = sum(n_burst);

src_wave = zeros(1, Nt);
src_wave(n_burst) = bh_win(N_b) .* sin(2 * pi * f * t(n_burst));
% -------------------------------------------------------------

% --- free-space reference simulation (to isolate E_inc at x_obs) ---
% pre-allocate E arrays
E_free = zeros(Nx, 3);
E_inc_obs = zeros(1, Nt);

% run loop
for ii = 1:Nt
    % source excitation
    E_free(1, 2) = src_wave(ii);

    % forward stepping update equation [ E(i, n+1) ]
    E_free(2:end-1, 3) = factor_free .* ...
        (E_free(3:end, 2) - 2*E_free(2:end-1, 2) + E_free(1:end-2, 2)) ...
        + 2*E_free(2:end-1, 2) - E_free(2:end-1, 1);
    % enforce boundary conditions at edges (E = 0)
    E_free(1, 3)   = 0;
    E_free(end, 3) = 0;

    % scoot array forward in time (n)
    E_free(:, 1)   = E_free(:, 2); % E(i, n-1) = E(i, n)
    E_free(:, 2)   = E_free(:, 3); % E(i, n)   = E(i, n+1)

    % store field data
    E_inc_obs(ii)  = E_free(idx_obs, 2);
end
% -------------------------------------------------------------------

% --- FDTD with coating (E_tot = E_inc + E_ref) ---
% pre-allocate E arrays
E = zeros(Nx, 3);
E_tot_obs = zeros(1, Nt);

% run loop
for ii = 1:Nt
    % source excitation
    E(1, 2) = src_wave(ii);

    % forward stepping update equation [ E(i, n+1) ]
    E(2:end-1, 3) = factor_coat .* (E(3:end, 2) - 2*E(2:end-1, 2) + E(1:end-2, 2)) ...
        + 2*E(2:end-1, 2) - E(2:end-1, 1);

    % enforce boundary conditions at edges (E = 0)
    E(1, 3)   = 0;
    E(end, 3) = 0;

    % scoot array forward in time (n)
    E(:, 1) = E(:, 2); % E(i, n-1) = E(i, n)
    E(:, 2) = E(:, 3); % E(i, n)   = E(i, n+1)
    
    % store field data
    E_tot_obs(ii) = E(idx_obs, 2);

    % update animation plot (every 10th timestep)
    if animate_this_run && (mod(ii, 10) == 0 || ii == Nt)
        % update wave amplitude handle
        set(h_wave, 'YData', E(:, 2));

        % update title with accurate timestep and frequency
        title(ax1, sprintf('1D FDTD Wave Propagation at %.1f MHz (Step %d / %d, t = %.1f ns)', ...
            f/1e6, ii, Nt, t(ii)*1e9), 'FontSize', 13);
        drawnow;
    end
end
% -------------------------------------------------

% pure reflected field at observation plane (since |E_inc| = 1 V/m)
E_ref_obs = E_tot_obs - E_inc_obs;

% summed complex phasor amplitudes for both reflected and incident waves
S_ref = sum(E_ref_obs .* exp(-1j * 2 * pi * f * t));
S_inc = sum(E_inc_obs .* exp(-1j * 2 * pi * f * t));

% calulate gamma magnitude and print it for this frequency
Gamma_fdtd(ff) = abs(S_ref) / max(abs(S_inc), eps);
fprintf('  f/f0 = %4.2f (%5.1f MHz) -> |Gamma_FDTD| = %.4f\n', ...
    rel_freqs_fdtd(ff), f/1e6, Gamma_fdtd(ff));
end
toc
% -------------------------------

%% ========================================================================
% Analytical solution
% build a big vector for relative frequencies to plot a smooth curve later
rel_freqs_ana = linspace(0.05, 2.0, 400);
f_ana        = rel_freqs_ana * f0;

% define intrinsic impedances for each medium
eta0 = sqrt(mu_0 / eps_0);
eta1 = eta0 / sqrt(eps_r1);
eta2 = eta0 / sqrt(eps_r2);
eta3 = eta1;
etaL = eta0;

% exact analytical reflection coefficient for symmetric N = 3 layers
Gamma_exact = zeros(size(f_ana)); % pre-allocate
for k = 1:length(f_ana)
    % omega and beta calculations
    w  = 2 * pi * f_ana(k);
    b1 = w * sqrt(mu_0 * eps_0 * eps_r1);
    b2 = w * sqrt(mu_0 * eps_0 * eps_r2);
    b3 = b1;

    % calculate input impedances for each medium
    Z3 = eta3 * (etaL + 1j * eta3 * tan(b3 * d3)) / (eta3 + 1j * etaL * tan(b3 * d3));
    Z2 = eta2 * (Z3   + 1j * eta2 * tan(b2 * d2)) / (eta2 + 1j * Z3   * tan(b2 * d2));
    Z1 = eta1 * (Z2   + 1j * eta1 * tan(b1 * d1)) / (eta1 + 1j * Z2   * tan(b1 * d1));

    % calculate gamma magnitudes
    Gamma_exact(k) = abs((Z1 - eta0) / (Z1 + eta0));
end

%% ========================================================================
% Comparison of solutions plot
fig2 = figure('Color', 'w', 'Theme', 'Light');
ax2 = axes('Parent', fig2);
hold(ax2, 'on');

plot(ax2, rel_freqs_ana, Gamma_exact, '-', 'Color', navy_au, 'LineWidth', 2.2);
plot(ax2, rel_freqs_fdtd, Gamma_fdtd, 'o', 'Color', orange_au, ...
    'MarkerFaceColor', orange_au, 'MarkerSize', 6.5, 'LineWidth', 1.5);

grid(ax2, 'on');
xlim(ax2, [0, 2.0]);
ylim(ax2, [0, 0.5]);
xticks(ax2, 0:0.5:2.0);
xlabel(ax2, 'Relative Frequency (f / f_0)', 'FontSize', 12);
ylabel(ax2, 'Input Reflection Coefficient |\Gamma_{in}|', 'FontSize', 12);
title(ax2, 'Frequency Response of 3-Layer Quarter-Wave ARC (\epsilon_r = [2, 3, 2], f_0 = 20 MHz)', ...
    'FontSize', 13);
legend(ax2, {'Analytical Solution', ...
             'Numerical Solution'}, ...
    'Location', 'north', 'FontSize', 10.5);

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