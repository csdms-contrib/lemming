% Dinf.m = A function that gathers all of the scripts to extract a
% D-infinity flow grid based on the Steven Eddins Mathowrks scripts

function [R slopes AccGrid] = Dinf( filename, dx, dy, NDVAL )

        load(filename)
        
    try
        if ~exist('dx','var')
            dx = cellsize;
        end
        
        if ~exist('dy','var')
            dy = cellsize;
        end
        
    catch
        input(dx,'Please specify a grid spacing in meters: ')
        dy = dx;
    end
    
    if ~exist('NDVAL','var')
        NDVAL = 0; %min(min(topo));
    end
    
    % First, make nodata values into NaNs
    topo(topo <= NDVAL) = nan;
    
    
        % Fill DEM sinks
            disp 'Filling DEM sinks...'
            
            % fill_sinks doesn't like nans
            topof = topo;
            topof(isnan(topof)) = -1;
    
            % do filling
            topof = fill_sinks(topof);   
            
            
        % Slopes and flow directions
            disp 'Computing slopes and flow directions...'
            [R, slopes] = dem_flow(topof,dx,dy);

            

            
        % Accumulation areas
            disp 'Creating flow matrix...'
            Ts = flow_matrix(topof,R);

            % restore nodata nans
            topo(~isnan(topo)) = topof(~isnan(topo));           
            
            disp 'Calculating upstream areas...'
            AccGrid = (dx * dy) * upslope_area(topo,Ts);
            
        % Save output
            disp(['Saving output as ' filename '_dinf.mat...'])
            save([filename '_dinf.mat']);
            
        % Display the accumulation grid
            [y x] = size(topo);
            Xs = 1:x;
            Ys = 1:y;
            Xs = dx .* Xs;
            Ys = dy .* Ys;
            
        
            figure
            imagesc(Xs, Ys, AccGrid); axis xy; axis image; colorbar; title('Accumulation area (m^2)','fontsize',14)
            xlabel('Easting (m)'); ylabel('Northing (m)')
            
            disp('Done!')
end