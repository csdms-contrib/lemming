function suite = test_facet_flow

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

%===============================================================================
function test_flow_within_facet
e0 = 10;    e1 = 8;    e2 = 6;
d1 = 2;     d2 = 3;

[r, s] = facet_flow(e0, e1, e2, d1, d2);

% Expected values calculated by hand using equations 1-5 in Tarboton paper.
expected_r = 0.588002603547568;
expected_s = 1.201850425154663;

assertElementsAlmostEqual(r, expected_r);
assertElementsAlmostEqual(s, expected_s);

%===============================================================================
function test_flow_on_east_facet_edge

% Flow for the following values should be directly to the east.
e0 = 10;    e1 = 9;    e2 = 10;
d1 = 2;     d2 = 3;

% Expected values calculated by hand using equations 1-5 in Tarboton paper.
[r, s] = facet_flow(e0, e1, e2, d1, d2);

expected_r = 0;
expected_s = 0.5;

assertElementsAlmostEqual(r, expected_r);
assertElementsAlmostEqual(s, expected_s);

%===============================================================================
function test_flow_on_northeast_facet_edge

% Flow for the following values should be directly towards the northeast
% neighbor.
e0 = 10;     e1 = 6;     e2 = -3;
d1 = 2;      d2 = 3;

% Expected values calculated by hand using equations 1-5 in Tarboton paper.
[r, s] = facet_flow(e0, e1, e2, d1, d2);
expected_r = atan2(d2, d1);
expected_s = sqrt(13);

assertElementsAlmostEqual(r, expected_r);
assertElementsAlmostEqual(s, expected_s);

%===============================================================================
function test_flow_outside_of_facet

% In Tarboton's method, the direction of facet flow is constrained to lie within
% the facet. If the equations indicate an angle outside the range 0 to atan2(d2,
% d1), the angle is clipped to that range, and the slope (s) is modified to
% match.

% For the following values, equations 1-3 in Tarboton indicate a flow angle less
% than 0.  It should be clipped to 0.
e0 = 10;     e1 = 9;     e2 = 10;
d1 = 2;      d2 = 3;

[r, s] = facet_flow(e0, e1, e2, d1, d2);
expected_r = 0;
expected_s = 0.5;
assertElementsAlmostEqual(r, 0);
assertElementsAlmostEqual(s, 0.5);

% For the following values, equations 1-3 in Tarboton indicate a flow angle
% greater than atan2(d2, d1).  It should be clipped to atan2(d2, d1).
e0 = 10;     e1 = 6;     e2 = -5;
d1 = 2;      d2 = 3;

[r, s] = facet_flow(e0, e1, e2, d1, d2);
expected_r = atan2(d2, d1);
expected_s = sqrt(13);

%===============================================================================
function test_flat_facet

% For a flat facet the angle is theoretically arbitrary, but we expect it to be
% returned as 0 because that's what atan2(0,0) conventionally returns.  The
% slope should be 0.

e0 = 0;     e1 = 0;     e2 = 0;
d1 = 2;     d2 = 3;

[r, s] = facet_flow(e0, e1, e2, d1, d2);
assertElementsAlmostEqual(r, 0);
assertElementsAlmostEqual(s, 0);

%===============================================================================
function test_uphill_facet

% For a facet going uphill, the angle would be the negative of the angle for the
% negated facet, but then the angle will be clipped as usual to the range 0 to
% atan2(d2, d1).  The slope will be negative.

e0 = -10;    e1 = -8;    e2 = -6;
d1 = 2;     d2 = 3;

[r, s] = facet_flow(e0, e1, e2, d1, d2);

% Expected values calculated by hand using equations 1-5 in Tarboton paper.
expected_r = 0;
expected_s = -1;

assertElementsAlmostEqual(r, expected_r);
assertElementsAlmostEqual(s, expected_s);

%===============================================================================
function test_matrix_input

% facet_flow is an elementwise function with respect to e0, e1, and e2. The
% outputs are expected to be the same size as e0, e1, and e2.

e0 = [1 1 1; 1 1 1];
e1 = [0 0 0; 0 0 0];
e2 = [1 .9 .8; .7 .6 .5];

d1 = 2; d2 = 3;
[r, s] = facet_flow(e0, e1, e2, d1, d2);

[M, N] = size(e0);
expected_r = zeros(M, N);
expected_s = zeros(M, N);
for p = 1:M
    for q = 1:N
        [expected_r(p, q), expected_s(p, q)] = facet_flow(e0(p, q), e1(p, q), ...
            e2(p, q), d1, d2);
    end
end

assertElementsAlmostEqual(r, expected_r);
assertElementsAlmostEqual(s, expected_s);

%===============================================================================
function test_default_input

% If d1 and d2 are not passed in, they should be assumed to be 1.
e0 = 1;    e1 = 0.8;    e2 = 0.7;

[r, s] = facet_flow(e0, e1, e2);
[expected_r, expected_s] = facet_flow(e0, e1, e2, 1, 1);

assertElementsAlmostEqual(r, expected_r);
assertElementsAlmostEqual(s, expected_s);



