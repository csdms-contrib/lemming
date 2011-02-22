function suite = test_border_nans

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

%===============================================================================
function test_typical_case

E = [NaN  0   0   0   0
     0    NaN 0   NaN 0
     NaN  0   0   0   0
     0    0   0   0   NaN];

result = border_nans(E);
expected_result = logical([1 0 0 0 0
                           0 1 0 0 0
                           1 0 0 0 0
                           0 0 0 0 1]);
                       
assertEqual(result, expected_result);

%===============================================================================
function test_no_border_nans

E = [1 2 3 4 5
     6 7 8 9 10];
 
result = border_nans(E);
expected_result = false(size(E));

assertEqual(result, expected_result);

%===============================================================================
function test_empty_input

% If input is empty, output should be a logical empty with the same size.
E = zeros(1, 0, 3);

result = border_nans(E);
expected_result = false(size(E));

assertEqual(result, expected_result);
