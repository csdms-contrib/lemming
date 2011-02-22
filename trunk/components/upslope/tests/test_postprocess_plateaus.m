function suite = test_postprocess_plateaus

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

Ap = postprocess_plateaus(A, E);

% The 4th and 5th columns of E form a plateau (pixels with no downhill
% neighbors).  So we expect the 4th and 5th columns of Ap to contain the mean of
% the 4th and 5th columns of A.  The elements (2,2) and (2,3) also form a
% plateau, but A(2,2) and A(2,3) already equal the mean of those two values.
exp_Ap = Ap;
exp_Ap(:, [4 5]) = mean(mean(A(:, [4 5])));

assertElementsAlmostEqual(Ap, exp_Ap);