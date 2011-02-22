function suite = test_influence_map

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

I = influence_map(E, T, 3, 5);

% Expected influence map for (3, 5) element was hand-computed.
exp_I = [ ...
    0, 0,   1/9,  0,   0, 0, 0
    0, 1/3, 2/9,  1/3, 0, 0, 0
    0, 2/3, 7/18, 1/3, 1, 0, 0
    0, 0,   5/18, 1/3, 0, 0, 0];

assertElementsAlmostEqual(I, exp_I);

function test_binaryImageInput(E)

bw = false(size(E));
bw(3, 5) = true;
bw(1, 7) = true;
R = dem_flow(E);
T = flow_matrix(E, R);
I1 = influence_map(E, T, [3 1], [5 7]);
I2 = influence_map(E, T, bw);
assertEqual(I1, I2);