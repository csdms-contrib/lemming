function suite = test_flow_matrix

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function test_flows
% Use an E matrix containing a two-element regional minimum, some plateau
% pixels, and some border NaNs.

E = [ ...
    2 2 2 2 2 3 NaN
    2 1 2 2 2 3 NaN
    2 1 2 2 2 3 NaN
    2 2 2 2 2 3 NaN];

% Expected flows hand-computed, using rules for regional minima, border NaNs,
% and plateau handling.

expected_flows = [ ...
    6 1 -1
    6 2 -1
    7 3 -1
    7 4 -1
    6 5 -1
    7 8 -1
    6 9 -1
    6 10 -1
    7 11 -1
    7 12 -1
    9 13 -1/2
    10 13 -1/2
    9 14 -1/3
    10 14 -1/3
    11 14 -1/3
    10 15 -1/3
    11 15 -1/3
    12 15 -1/3
    11 16 -1/2
    12 16 -1/2
    13 17 -1/2
    14 17 -1/2
    13 18 -1/3
    14 18 -1/3
    15 18 -1/3
    14 19 -1/3
    15 19 -1/3
    16 19 -1/3
    15 20 -1/2
    16 20 -1/2
    17 21 -1
    18 22 -1
    19 23 -1
    20 24 -1];

expected_flows = [expected_flows;
    (1:28)', (1:28)', ones(28, 1)];

expected_T = sparse(expected_flows(:, 1), expected_flows(:, 2), ...
    expected_flows(:, 3), 28, 28);

R = dem_flow(E);
T = flow_matrix(E, R);

assertEqual(T, expected_T);

function test_defaultDs
% d1 and d2 should be assumed to be 1 by default.

E = magic(5);
R = dem_flow(E);
T1 = flow_matrix(E, R);
T2 = flow_matrix(E, R, 1, 1);

assertEqual(T1, T2);


