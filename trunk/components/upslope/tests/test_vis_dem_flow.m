function suite = test_vis_dem_flow

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function test_data = setup
s = load('milford_ma_dem');
test_data.E = s.Z;
[test_data.R, test_data.S] = dem_flow(test_data.E);

function teardown(test_data)
close all

function test_smallInput(test_data);
% Make sure vis_dem_flow executes without error for DEM smaller than 50x50.
E = test_data.E(1:40, 1:40);
R = test_data.R(1:40, 1:40);
S = test_data.S(1:40, 1:40);
vis_dem_flow(E, R, S);

function test_largeInput(test_data)
% Make sure vis_dem_flow executes without error for DEM larger than 50x50.
vis_dem_flow(test_data.E, test_data.R, test_data.S);

