function suite = test_dependence_map

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function E = setup

% Use an E matrix containing a two-element regional minimum, some plateau
% pixels, and some border NaNs.

E = [ ...
    2 2 2 2 2 3 NaN
    2 1 2 2 2 3 NaN
    2 1 2 2 2 3 NaN
    2 2 2 2 2 3 NaN];



function test_normalCase(E)

R = dem_flow(E);
T = flow_matrix(E, R);

D = dependence_map(E, T, 2, 2);

% Expected dependence map for (2,2) element was hand-computed.
test_data.exp_D = [ ...
    1, 1, 1, 1,   (1/2 + 1/3), (1/2 + 1/3), 0
    1, 1, 1, 2/3, 2/3,         2/3,         0
    0, 0, 0, 1/3, 1/3,         1/3,         0
    0, 0, 0, 0,   1/6,         1/6,         0];

assertElementsAlmostEqual(D, test_data.exp_D);

function test_binaryImageInput(E)

bw = false(size(E));
bw(3, 2) = true;
R = dem_flow(E);
T = flow_matrix(E, R);
D1 = dependence_map(E, T, 3, 2);
D2 = dependence_map(E, T, bw);
assertEqual(D1, D2);