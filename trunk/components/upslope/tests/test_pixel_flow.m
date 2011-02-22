function suite = test_pixel_flow

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function test_angles
% Move the downhill facet around to all eight positions and make sure the angle
% moves appropriately.

% Use E - NE as the reference facet.
E = [1.0 1.0 0.4
     1.0 1.0 0.5
     1.0 1.0 1.0];
 
[r, s] = pixel_flow(E, 2, 2, 1, 1);

% Test NE - N facet.
E2 = [1.0 0.5 0.4
      1.0 1.0 1.0
      1.0 1.0 1.0];
[r2, s2] = pixel_flow(E2, 2, 2, 1, 1);
assertElementsAlmostEqual(r2, pi/2 - r);
assertElementsAlmostEqual(s2, s);

% Test N - NW facet.
E3 = [0.4 0.5 1.0
      1.0 1.0 1.0
      1.0 1.0 1.0];
[r3, s3] = pixel_flow(E3, 2, 2, 1, 1);
assertElementsAlmostEqual(r3, pi/2 + r);
assertElementsAlmostEqual(s3, s);

% Test NW - W facet.
E4 = [0.4 1.0 1.0
      0.5 1.0 1.0
      1.0 1.0 1.0];
[r4, s4] = pixel_flow(E4, 2, 2, 1, 1);
assertElementsAlmostEqual(r4, pi - r);
assertElementsAlmostEqual(s4, s);

% Test W - SW facet.
E5 = [1.0 1.0 1.0
      0.5 1.0 1.0
      0.4 1.0 1.0];
[r5, s5] = pixel_flow(E5, 2, 2, 1, 1);
assertElementsAlmostEqual(r5, pi + r);
assertElementsAlmostEqual(s5, s);

% Test SW - S facet.
E6 = [1.0 1.0 1.0
      1.0 1.0 1.0
      0.4 0.5 1.0];
[r6, s6] = pixel_flow(E6, 2, 2, 1, 1);
assertElementsAlmostEqual(r6, 1.5*pi - r);
assertElementsAlmostEqual(s6, s);

% Test S - SE facet.
E7 = [1.0 1.0 1.0
      1.0 1.0 1.0
      1.0 0.5 0.4];
[r7, s7] = pixel_flow(E7, 2, 2, 1, 1);
assertElementsAlmostEqual(r7, 1.5*pi + r);
assertElementsAlmostEqual(s7, s);

% Test SE - E facet.
E8 = [1.0 1.0 1.0
      1.0 1.0 0.5
      1.0 1.0 0.4];
[r8, s8] = pixel_flow(E8, 2, 2, 1, 1);
assertElementsAlmostEqual(r8, 2*pi - r);
assertElementsAlmostEqual(s8, s);

function test_allNeighborsUphill
E = [2 1 2; 1 0 1; 2 1 2];
[r, s] = pixel_flow(E, 2, 2, 1, 1);
assertEqual(r, NaN);
assertEqual(s, -1);
     
function test_plateau
[r, s] = pixel_flow(zeros(3, 3), 2, 2, 1, 1);
assertEqual(r, NaN);
assertEqual(s, 0);

function test_borderNaNs
% NaNs on border should act like pixels just a bit higher than the maximum of
% the other pixels.
E = [NaN NaN NaN
     NaN 0   0
     NaN 0   0];
[r, s] = pixel_flow(E, 2, 2, 1, 1);
assertEqual(r, NaN);
assertEqual(s, 0);

function test_vectorizedComputation
E = [0.4 0.5 1.0 1.0 0.5 0.4
     1.0 1.0 1.0 1.0 1.0 1.0
     1.0 1.0 1.0 1.0 1.0 1.0
     1.0 1.0 1.0 1.0 1.0 1.0
     1.0 1.0 1.0 1.0 1.0 1.0
     0.4 0.5 1.0 1.0 0.5 0.4];
 
 i = [2 2; 5 5];
 j = [2 5; 5 2];
 
 [r, s] = pixel_flow(E, i, j, 1, 1);
 assertEqual(size(r), [2 2]);
 assertEqual(size(s), [2 2]);
 
 assertElementsAlmostEqual(r(1, 2), pi - r(1, 1));
 assertElementsAlmostEqual(s(1, 2), s(1, 1));
 
 assertElementsAlmostEqual(r(2, 1), pi + r(1, 1));
 assertElementsAlmostEqual(s(2, 1), s(1, 1));
 
 assertElementsAlmostEqual(r(2, 2), 2*pi - r(1, 1));
 assertElementsAlmostEqual(s(2, 2), s(1, 1));
 
function test_defaultDs
E = magic(3);
[r, s] = pixel_flow(E, 2, 2);
[r2, s2] = pixel_flow(E, 2, 2, 1, 1);
assertEqual(r, r2);
assertEqual(s, s2);