function test_suite = test_fill_sinks

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function testDiagToBorderCase
% Diagonal lines that run to border should not be filled.

E = [1 2 2 2 1
    2 1 2 1 2
    2 2 1 2 2
    2 1 2 1 2
    1 2 2 2 1];

assertEqual(E, fill_sinks(E));
