clear;
close all;
clc;

%% ========================================================================
% General space setup and variable declarations
% auburn_color_map = Generate_Auburn_Colormap_v0();

% --- make the simulation space ---
% dx = 0.001; % dx = 1mm
% x_max = 20; % 20m
% x = 0:dx:x_max;

% dt = 0;
% ---------------------------------


% constitutive parameters
% eps_0 = 8.854e-12; 



% implement boundary conditions


% --- build anti-reflective coating ---

% define plate locations/dimensions as index ranges
% z_shutter = 120:124;
% z_sense   = 109:113;
% z_ground  = 101:105;
% 
% x_sense1  = 150:245; % sense plate 1 (left)
% x_sense2  = 255:350; % sense plate 2 (right)
% x_ground  = 150:350; % bottom ground plate

% -------------------------------


% plot the space with initial conditions (example purposes)
% fig1 = figure;
% imagesc(x, z, V');
% colormap(auburn_color_map);
% c1 = colorbar; c1.Label.String = '(V/m)'; c1.Label.FontSize = 14;
% axis xy;
% title(sprintf('Voltage Distribution of the Space at t=0s'), 'FontSize', 14);
% xlabel('Length (m)', 'FontSize', 12); ylabel('Height (m)', 'FontSize', 12);

%% ========================================================================
% 1D FDTD LOOP

% set up Figure 2 to animate the motion during runtime
% fig2 = figure;
% colormap(auburn_color_map); 

% tic % the whole loop takes ~XXs to run
% for ii = 1:num_timesteps
% 
%     % update animation plot
%     colormap(auburn_color_map);
%     c2 = colorbar; c2.Label.String = 'Electric Field Intensity (V/m)';
%     c2.Label.FontSize = 14;
%     axis xy;
%     title(sprintf('Voltage Distribution (Step %d / %d)', ...
%                                       ii, num_timesteps), 'FontSize', 14);
%     xlabel('Length (m)', 'FontSize', 12); 
%     ylabel('Height (m)', 'FontSize', 12);
%     drawnow;
% end
% toc
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
% according to https://ocm.auburn.edu/brand-center/_assets/pdf/au-colorpalette-primary-supporting.pdf

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
auburn_color_map(1:length_colormap/2,1) = (interp1([1,length_colormap/2], [blue(1) white(1)], 1:length_colormap/2))';
auburn_color_map(1:length_colormap/2,2) = (interp1([1,length_colormap/2], [blue(2) white(2)], 1:length_colormap/2))';
auburn_color_map(1:length_colormap/2,3) = (interp1([1,length_colormap/2], [blue(3) white(3)], 1:length_colormap/2))';

% make white to orange
auburn_color_map(length_colormap/2+1:length_colormap,1) = (interp1([length_colormap/2+1,length_colormap], [white(1) orange(1)], length_colormap/2+1:length_colormap))';
auburn_color_map(length_colormap/2+1:length_colormap,2) = (interp1([length_colormap/2+1,length_colormap], [white(2) orange(2)], length_colormap/2+1:length_colormap))';
auburn_color_map(length_colormap/2+1:length_colormap,3) = (interp1([length_colormap/2+1,length_colormap], [white(3) orange(3)], length_colormap/2+1:length_colormap))';

% RGB values normalized from [0:1]
auburn_color_map = auburn_color_map/255;
end

% ----------------------------------