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

c     = 1 / sqrt(mu_0 * eps_0); % speed of light

f0      = 20e6;   % design frequency: 20 MHz
lambda0 = c / f0; % free-space wavelength (~14.99 m)

% --- make the simulation space ---
dx = 0.025; % dx = 2.5cm
x_max = 150; % 150m
x = 0:dx:x_max;
Nx = length(x);

% Courant stability factor (starting with magic time step)
S = 1.0; 
dt = S * dx / c; % ~83.39 ps
t_max = 400e-9;  % 400 ns (8 full cycles at 20 MHz)
t = 0:dt:t_max;
Nt = length(t);

E = zeros(3, length(x)); % pre-allocate E array
eps_r = ones(1, Nx);  % pre-allocate permittivity array
% ---------------------------------


% --- build anti-reflective coating ---
% quarter-wavelength layer thicknesses
d1 = lambda0 / (4 * sqrt(eps_r1)); % ~2.6498 m
d2 = lambda0 / (4 * sqrt(eps_r2)); % ~2.1636 m
d3 = d1;                           % ~2.6498 m

% positions in meters (roughly centered in the space to prevent PEC 
% boundary reflections before t_max)
x_src   = 65.0;  % source excitation position
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


% plot the space with initial conditions (example purposes)
% fig1 = figure;
% colormap(auburn_color_map);
% c1 = colorbar; c1.Label.String = '(V/m)'; c1.Label.FontSize = 14;
% axis xy;
% title(sprintf('Voltage Distribution of the Space at t=0s'), 'FontSize', 14);
% xlabel('Length (m)', 'FontSize', 12); ylabel('Height (m)', 'FontSize', 12);

%% ========================================================================
% 1D FDTD LOOP

% set up Figure 2 to animate the motion during runtime
fig2 = figure;
colormap(auburn_color_map); 

% set frequency for this run
f = f0;

tic % the whole loop takes ~XXs to run
for ii = 1:Nt
    % source excitation 
    % (simply turning a sine wave on for 4 periods)
    if t(ii) < 4/f
        E(2,1) = sin(2*pi*f*t(ii));
    else
        E(2,1) = 0;
    end

    % factor multiplying the first three terms
    factor = (1/(mu_0*eps_0*eps_r(ii)))*((dt/dx)^2);

    % forward stepping update equation [ E(i, n+1) ]
    E(2:end-1, 3) = factor*(E(3:end, 2) - 2*E(2:end-1, 2) + E(1:end-2, 2)) ...
        + 2*E(2:end-1, 2) - E(1:end-2, 1);

    % scoot array in n (time)
    E(2:end-1, 0) = E(2:end-1, 1); % E(i, n-2) = E(i, n-1)
    E(2:end-1, 1) = E(2:end-1, 2); % E(i, n-1) = E(i, n)

    % update animation plot
    colormap(auburn_color_map);
    plot(z, E(2,:), 'LineWidth', 3);
    % c2 = colorbar; c2.Label.String = 'Electric Field Intensity (V/m)';
    % c2.Label.FontSize = 14;
    % axis xy;
    title(sprintf('Propagation (Step %d / %d)', ...
                                      ii, Nt), 'FontSize', 14);
    xlabel('Length (m)', 'FontSize', 12); 
    % ylabel('Height (m)', 'FontSize', 12);
    drawnow;
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
% Plotting


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