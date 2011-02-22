function suite = test_vis_map

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function setup

function teardown
close all

function test_example
% Make sure vis_map example executes without error.
E = peaks;
R = dem_flow(E);
T = flow_matrix(E, R);
D = dependence_map(E, T, 12, 24);
vis_map(D, E, 12, 24)
title('Dependence map')
