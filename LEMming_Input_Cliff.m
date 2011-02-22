% LEMming - a MATLAB-based landscape evolution model
% Copyright (C) 2011 Dylan Ward
% 
% 
% Developer can be contacted at djward@unm.edu 
% 
% and 
% 
% University of New Mexico, 
% Dept. of Earth and Planetary Sciences, 
% Albuquerque, NM 87131
% 
% 
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 2 of 
% the License, or (at your option) any later version. See License.txt 
% and gpl-2.0.txt for the full terms of this license.
% -------------------------------------------------------------------

% LEMming_Input.m - Input script to set parameters for a LEMming (V015 and up) run.

clear all
close all
clc

Setup_Name = 'Cliff';  % As in, LEMming_Input_<Setup_Name>.m
run_name = 'Same, debris rdot 0.005';

% Load a topo file
inFile = 'landscapes/500mCliff200mBack_nscales.mat';
%% Comment this line out to generate 
%% a synthetic DEM according to the parameters below:

% grid size input

x = 80;     % pixels x
y = 120;     % pixels y
dx = 10;     % meters per pixel
dy = dx;

% grid z parameters

z_init = 500;              % Initial grid elevation (m)
z_roughness = 5;            % Height scale factor (relative to initial slope) of roughness in initial grid (m)
z_bound = 0;                % z_init;       % Fixed-elevation boundary elevation (m)
tilt = 0.4;                   % Apply a slope to the intial grid (m/m)
DO_RIDGE = 1;               % Nonzero to make a ridge axial at y=ymax/2

borderwidth = 5;        % Width of edge over which boundary condition is applied, cells

% timesteps and model duration

dt_init = .001;                   % year
dt_slope_threshold = 1;         % As a ratio. This is the most slope change allowed during a given time step, as a fraction of the current slope
dt_min = 1;                   % Minimum timestep (years)
dt_max = 1000;                 % Maximum timestep

z_max_plot = z_init;
cinterval = 50;                % contour interval for plotting, m
tmax = 1500;                   % years
plottime = 500;                % years between plots
tracktime = 0.1 * plottime;                % years between recording the current erosion rates

VE = 1;                         % (times) Vertical exaggeration of plot

% Toggle saving mode. 
% 0: Only saves beginning and end states.
% 1: Saves each figure, and beginning and end states.
% 2: Saves entire workspace each time a figure is plotted.
SAVEMODE = 2;          

% physical parameters
rock_uplift = 0;        % Rock uplift rate relative to boundary condition (m/yr)
stream_width_coeff = .01; % Multiplies the width = sqrt(area) function. Zero for constant-width channels, set by w0
w0 = 10;   % Base width if stream_width_coeff == 0. 

% rocktype 0 - default substrate
k0 = 5e-6;              % fluvial erodibility constant
m0 = 1;
n0 = 1;
rdot0 = .01;         % bare-rock regolith production rate, m/yr 
rstar0 = .01;         % e-folding depth for falloff of regolith production
rfslope0 = inf;         % Positive slope (L/L) above which qualifies a rockfall source

% additional rocktypes, beginning with 1
%        [index k    m n rdot rstar rfslope]
RTproplist = [1 1e-7 1 1 1e-6 .003  3; ...
              2 1e-7 1 1 .001   .1  inf]';
%              2 1e-7 1 1 1e-6 .1    0.6]';
%              2 5e-7 1 1 5e-3 .1    inf]';
              %3 1e-1 0.7 1.5  140   1]';
              
% Regolith properties
kappa = .1;      % Regolith transport coefficient (m/year)
sc = inf;        % Critical slope (L/L) above which diffusion is nonlinear with slope (Roering) or e-folding scale for slope (exponential models)
Ac = 1000;       % Approximate critical drainage area for channelized bedload transport
mt = 1;         % area exponent on fluvial bedload transport equation
nt = 1;         % slope exponent on fluvial bedload transport equation
Max_Mobile_H = 0.05; %0.1; % meters - maximum thickness of mobile layer of regolith

kt = (kappa * mean(dx,dy)) / Ac.^mt;        % fluvial bedload transport efficiency
% UregSuspRef = .1;                               % m/yr - when bedload is moving this fast...
% SuspRef = 10;                               % ... the suspended fraction is this many times the bedload

ArDL = 70000;     % Reference area at which suspension slope is defined
SrDL = 0.1;   % Slope at which all sediment moves in suspension given the reference drainage area



% Rockfall parameters
DO_ROCKFALL = 1;

RF_Debris_Rtype = 2;    % Which of the above defined rocktypes represents rockfall debris?
RF_Debris_PlotCutoff = 0; % m - thicknesses of debris less than this are plotted as the underlying substrate

% Rockfall source
RFSource_curv = 0; % Curvature below which qualifies a rockfall source (universal - threshold slope is set by rock type.)

DepoAngleCutoff = 30;   % Degrees - angle from surface normal beyond which deposition probability is zero
distStar = 100; % e-folding lengthscale (meters) for falloff of deposition with distance
distMax = inf; % No deposition allowed more than this far from a source
curvStar = 0.25; % Approximate curvature where probability asymptotes to zero (when -) and one (when +)
slopeExp = 1;  % Deposition falls off as 1/slope^slopeExp
sDepCrit = 1.1; % Above this slope deposition probability is zero.

Backwearing_Rate = inf;    % m/yr - Backwearing erosion rate cap on rockfall
RF_Ht = 50; % m - Typical height of a rockfall event, usually the thickness of the hardcap is a good default

% Convert assigned backwearing rate to a volumetric erosion rate
RFerodeRate = Backwearing_Rate * (x * dx) * RF_Ht;

% additional stratigraphy
%stratlist = [1 10 10 -500 40 40 -600; 1 30 30 -800 40 40 -1800; 1 10 20 -1500 60 60 -1600];
%stratlist = [1 1 1 800 80 80 750; 1 1 1 -50 80 80 -100; 1 1 1 -300 40 80 -1000];
%bookcliffs160 = [1 1 1 1700 160 160 1600; 1 1 1 1100 160 160 900; 1 1 1 700 160 160 650; 1 1 1 550 160 160 400; 1 1 1 300 160 160 200];
%bookcliffssimple100 = [1 1 1 1200 100 100 1100; 1 1 1 0 100 100 -100];
%multicliffs = [1 1 1 1700 x y 1600; 1 1 1 1100 x y 900; 1 1 1 700 x y 650; 1 1 1 550 x y 400; 1 1 1 300 x y 200];
%twocliffs = [1 1 1 750 x y 650; 1 1 1 300 x y 250];
smallcliff = [1 1 1 z_init+10 x y z_init-50];
%nocliff = [1 1 1 z_init+1000 x y z_init+500];
%oneblob = [1 20 20 1000 60 60 800];
%clams = ones(1000,7);clams(:,2) = ceil(rand(1000,1).*x-1);clams(:,3) = ceil(rand(1000,1).*x-1);clams(:,4) = ceil(rand(1000,1).*4000);clams(:,4) = clams(:,4) - 3000;clams(:,5) = clams(:,2)+1;clams(:,6) = clams(:,3)+1;clams(:,7) = clams(:,4)-10;clams(clams == 0) = 1;
%HI_50 = [1 1 1 1000 20 y -50000; 1 40 1 1000 60 y -3000];
%tiltstrat = ones(y,7);tiltstrat(:,3) = 1:1:y;tiltstrat(:,4) = .015 * dy .* (tiltstrat(:,3)-round(4*y/5)) + z_init;tiltstrat(:,6) = tiltstrat(:,3);tiltstrat(:,5) = tiltstrat(:,5).*x;tiltstrat(:,7) = tiltstrat(:,4)-50;tiltstrat(tiltstrat == 0) = 1;
%lefttiltstrat = ones(x,7);lefttiltstrat(:,2) = 1:1:x;lefttiltstrat(:,4) = .07 * dx .* (lefttiltstrat(:,2)-round(4*y/5)) + z_init;lefttiltstrat(:,5) = lefttiltstrat(:,2);lefttiltstrat(:,6) = lefttiltstrat(:,6).*y;lefttiltstrat(:,7) = lefttiltstrat(:,4)-50;lefttiltstrat(lefttiltstrat == 0) = 1;
%backtiltstrat = ones(y,7);backtiltstrat(:,3) = 1:1:y;backtiltstrat(:,4) = -.015 * dy .* (backtiltstrat(:,3)-round(y/5)) + z_init;backtiltstrat(:,6) = backtiltstrat(:,3);backtiltstrat(:,5) = backtiltstrat(:,5).*x;backtiltstrat(:,7) = backtiltstrat(:,4)-50;backtiltstrat(backtiltstrat == 0) = 1;

stratlist = smallcliff';


%% CALL LEMMING %%
%LEMming015_silent
LEMming015


%%%%%% End of LEMming_Input.m %%%%%%----------------------------------------