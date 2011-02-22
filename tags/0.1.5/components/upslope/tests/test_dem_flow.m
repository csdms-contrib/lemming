function suite = test_dem_flow

%   Steven L. Eddins
%   Copyright 2007-2009 The MathWorks, Inc.

initTestSuite;

function E = setup
% Matrix E extracted from values (50:55, 60:65) of
% Z inside milford_ma_dem.mat.  Then we pad by replication around all the
% borders.

E = [ ...
   153   154   154   153   151   151
   150   151   152   151   148   146
   149   149   150   149   145   141
   148   148   148   146   142   138
   146   145   144   142   138   134
   143   142   141   139   134   138];

function test_normalCase(E)

% Pad E and then compute the pixel flow for all the interior of elements
% of the padded matrix.  That should give the same result as calling 
% dem_flow.

expected_R = zeros(size(E));
expected_S = zeros(size(E));
[M, N] = size(expected_R);

% Could use padarray here, but let's choose a computation that's more
% independent of what dem_flow done.

Ep = [E(1, :); E];
Ep = [Ep; Ep(end, :)];
Ep = [Ep(:, 1), Ep];
Ep = [Ep, Ep(:, end)];

for p = 1:M
    for q = 1:N
        [expected_R(p,q), expected_S(p,q)] = pixel_flow(Ep, p+1, q+1);
    end
end

[R, S] = dem_flow(E);

assertEqual(R, expected_R);
assertEqual(S, expected_S);

function test_defaultDs(E)

% d1 and d2 should be assumed to be 1 if not provided.
[R1, S1] = dem_flow(E);

[R2, S2] = dem_flow(E, 1, 1);

assertEqual(R1, R2);
assertEqual(S1, S2);
