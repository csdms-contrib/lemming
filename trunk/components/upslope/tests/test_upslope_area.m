function suite = test_upslope_area

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function test_normalCase

E = [ ...
    2 2 2 2 2 3 NaN
    2 1 2 2 2 3 NaN
    2 1 2 2 2 3 NaN
    2 2 2 2 2 3 NaN];

% Expected upslope area computed by hand for this E.
exp_A = [ ...
    1, 1,  31/9, 8/3,  2, 1, 0
    1, 12, 41/9, 10/3, 2, 1, 0
    1, 12, 41/9, 10/3, 2, 1, 0
    1, 1,  31/9, 8/3,  2, 1, 0];

R = dem_flow(E);
T = flow_matrix(E, R);
A = upslope_area(E, T);

assertElementsAlmostEqual(A, exp_A);