% Topoman.m - Mighty maker of arbitrary topography for LEMming initial
% conditions

% Could load a dem here...
clc
clear all

filename = '500mCliff200mBack_nscales.mat';

TopoPoints(1,:) = [0 0 100 100 0 100 0 100];       % X coordinates, Percent
TopoPoints(2,:) = [0 2 0   2  10  10  100 100];       % Y coordinates, Percent
TopoPoints(3,:) = [0 0  0   0   500 500 490 490];    % Z coordinates, Meters

% Grid size to create. Remember that LEMming automatically subsets large
% grids, so it's ok to make it bigger than will be used in a typical run...
x = 1000;           % pixels
y = 1000;           % 
cellsize = 10;      % This cellsize will be used by LEMming unless overridden
z_roughness = 20;            % Height scale of roughness in initial grid (m)
%tilt = 1;           % steepest slope in initial grid (only affects noise)

topo = zeros(y,x);

[XI YI] = meshgrid(1:x,1:y);

% Convert percent coordinates into grid coordinates
Xs = ceil((TopoPoints(1,:) ./ 100) .* x);
Ys = ceil((TopoPoints(2,:) ./ 100) .* y);

topoResult = griddata(Xs,Ys,TopoPoints(3,:),XI,YI,'linear');

topoBase = isnan(topoResult);
topo(topoBase) = min(min(topoResult));
topo(~topoBase) = topoResult(~topoBase);

% Make random roughness at all scales of the DEM

filtpad = 4 * z_roughness;  % pad boundaries of noise so filter border effects don't show in the topo

 %   z_roughness = z_roughness * tilt * GridDelta;
    noise = rand(y+2*filtpad,x+2*filtpad);
    
    noisei = noise;
for nscale = 3:2:z_roughness
    % Filter noise to reduce sharply closed depressions 
    noisei = noisei + filter2(ones(nscale)/nscale^2,noise);
    % Add it to the topography
end
    
    noiset = noisei(filtpad+1:end-filtpad,filtpad+1:end-filtpad);
    noiset = noiset - mean((mean(noiset)));
    noiset = noiset .* z_roughness/(2 * max(max(noiset)));
    topo = topo + noiset;

figure(13)
clf
imagesc(topo) 
% figure(666)
% clf
% surfc(topo)
%shading interp
%view(-55,55)

save(filename,'topo','cellsize');