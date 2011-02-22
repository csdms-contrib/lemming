% LEMming - a MATLAB-based landscape evolution model
% Copyright (C) 2011 Dylan Ward
% 
% 
% Developer can be contacted at djward@unm.edu 
% 
% and 
% 
% University of New Mexico, 
% Dept. of Earth and Planetary Sciences, 
% Albuquerque, NM 87131
% 
% 
% This program is free software; you can redistribute it and/or 
% modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; either version 2 of 
% the License, or (at your option) any later version. See License.txt 
% and gpl-2.0.txt for the full terms of this license.
% -------------------------------------------------------------------

%% A quick advection-diffusion LEM based on Taylor Perron's branching
%% ideas, with the ability to add layers of variable lithology

%% D. Ward, 12/16/2007: Version 002
%% Version 002 adds dynamic timestepping
%% Version 003 adds file input for DEMs
%% Version 004 changes the accumulation area algorithm, hopefully to a faster
%% one. Also adds movie support, improvements to slope calculations.
%% Version 005 adds support for multiple rock types
%% Version 006 reorders the calculation so diffusion is applied to the
%% post-advection topography, in an attempt to improve stability at larger
%% timesteps. Various improvements to plotting.
%% Version 007 fixes stratigraphy bugs, adds options for boundary
%% conditions, and fixes several typos in the introductory comments.
%% Version 009 fixes more stratigraphy bugs, slope calculation bug, other
%% bugs, adds rockfall, file saving and MakeLEMmingMov.m support
%% Version 010 changes the slope calculation to use the pixel flow slopes
%% output by dem_flow.m, and improves rockfall rate control
%% Version 011 reduces the nonlinearities in rockfall distribution in an
%% effort to improve stability at longer timesteps, adds Roering nonlinear
%% transport, and adds the ability to adjust the critical slope for
%% rockfall on a rocktype by rocktype basis
%% Version 012 adds regolith tracking, exponential regolith production, as
%% well as exponentially increasing nonlinear transport, which should be
%% more stable than Roering transport. Also upgrades the timestep
%% calculation.
%% Version 013 was a failed experiment.
%% Version 014 adds a discharge threshold to allow the upper extent of
%% all-bedrock channels to be set manually, sidestepping suspended sediment
%% transport. Also improves the regolith flux calculations. Also adjusts
%% the way rockfall debris is distributed as a function of slope,
%% curvature, and regolith thickness.
%% Version 015 changes the channel behavior (again), allowing a
%% suspended-load fraction to be set, based on stream power. Adds
%% directional rockfall distribution. 

% clear all
% close all
% clc

% Variables should be specified in LEMming_Input_<setupname>.m
% or otherwise set by script before running LEMming


%%%%%%%%%%%%%%%%%%%% EXECUTION %%%%%%%%%%%%%%%%%%%%

% Set path so components can be found
addpath('components','components/upslope','tools','landscapes');

% Prompt for a model run name
run_name = input('What would you like to call this model run? > ','s');
disp('Creating folder for this model run...')
run_filename = ['RUNS/' Setup_Name '/' run_name ' ' datestr(now,'dd-mmmm-yyyy HH.MMh')];
mkdir(run_filename);        % Create a subfolder for this run

% Copy current codes into the new subfolder so the run can be replicated
copyfile('LEMming015.m',['./' run_filename]);
copyfile('tools/MakeLEMmingMov.m',['./' run_filename]);
copyfile('components/RF_Spread_lemLink.m',['./' run_filename]);
copyfile(['LEMming_Input_' Setup_Name '.m'],['./' run_filename]);

disp('Executing...')


% Make XY absolute position grids
    [Xs, Ys] = meshgrid(1:x,1:y);
    Xs = Xs .* dx;
    Ys = Ys .* dy;

% Grid reference areas
    CellArea = (dx * dy);
    GridArea = (x * y) * CellArea;
    GridDelta = mean(dx,dy);  % Use as "universal" grid spacing when dx ~= dy
    
    

if exist('inFile') %#ok<EXIST>
load(inFile);

topo = topo(1:y,1:x); % Subset to X,Y.
%dx = cellsize; dy = dx; % Use cellsize from file (comment out to override to dx specified above)
z_bound = min(min(topo));
topo = fill_sinks(topo);    % Fill DEM sinks

else


% Make the grid
topo = ones(y,x) .* z_init;

% Apply the slope
init_slope = zeros( size(topo) ); 

for row = 1:1:length( topo(:,1) )
    init_slope(row,:) = tilt * row * dy;
end

topo = topo + init_slope;

    if DO_RIDGE
        % Symmetrical about y = ymax/2
        topo1 = min(topo, flipud(topo));
        topo2 = min(topo', fliplr(topo')); %#ok<UDIM>
        topo = min(topo1,topo2);
    end 
    
    
    
        % Make random roughness at all scales of the DEM

        filtpad = 4 * z_roughness;  % pad boundaries of noise so filter border effects don't show in the topo

         %   z_roughness = z_roughness * tilt * GridDelta;
            noise = rand(y+2*filtpad,x+2*filtpad);
            
            noisei = noise;
        for nscale = 3:2:z_radius
            % Filter noise to reduce sharply closed depressions 
            noisei = noisei + filter2(ones(nscale)/nscale^2,noise);
            % Add it to the topography
        end

            noiset = noisei(filtpad+1:end-filtpad,filtpad+1:end-filtpad);
            noiset = noiset - mean((mean(noiset)));
            noiset = noiset .* z_roughness/(2 * max(max(noiset)));
            topo = topo + noiset;

            % Correct for boundary area
            topo = topo - tilt*(borderwidth+1)*GridDelta; %min(topo(borderwidth+1,:));
    

    
    
end % if inFile

topo_init = topo;       % For double-checking erosion volumes later
z_max_plot = max(z_max_plot,max(max(topo))) + rock_uplift * tmax;
colordef none

% Miscellaneous constants
    e = exp(1);
    sp_susp = ArDL.^mt .* SrDL.^nt;
    
% Prepare the stratigraphy grids
    stratlist = [stratlist; 1:length(stratlist(1,:))];   %add unit ids to strat list row 8
    numStypes = length(RTproplist(1,:));

% Initialize property grids
    k = ones(size(topo)) .* k0;          
    m = ones(size(topo)) .* m0;           
    n = ones(size(topo)) .* n0;        
    rdot = ones(size(topo)) .* rdot0;          
    rstar = ones(size(topo)) .* rstar0;          
    rfslope = ones(size(topo)) .* rfslope0; 

    RTGrid = zeros(size(topo));
    rtype_here = zeros(size(topo));

% Initialize the surface stratigraphy layers
    RF_Debris_H = zeros(size(topo));
    Regolith_H = zeros(size(topo));
    Reg_overshoot = zeros(size(topo));
    Flux_overshoot = zeros(size(topo));
    
    dzCum = zeros(size(topo));
    RFdzdt = zeros(size(topo));
    dzSuspSed = zeros(size(topo));
    
    Adep = zeros(size(topo));
    
% Initialize flow-related grids
    StreamWs = ones(size(topo)) * w0;


% Initialize the time loop and counters
    t = 0;
    dt = dt_init;
    t_plot = 0;
    t_track = 0;
    erodeVol = 0;
    trackstep = 0;
    NoDistTargets = 0;
    
% Set up erosion rate tracking
    numRecSteps = ceil(tmax/tracktime) + 1;
    FLUVtrackVol = zeros(1,numRecSteps);
    DIFFtrackVol = zeros(2,numRecSteps);
    RFtrackVol = zeros(2,numRecSteps);
    TIMEtrackVol = zeros(1,numRecSteps);

    clear numRecSteps;

    if SAVEMODE == 1 || SAVEMODE == 2
        stateNo = 0;
    end

% Set up slope filter
    window = 3;
    filtpatch = ones(window) / window^2;

% Set up flow direction filter
    Rfiltpatch = ones(3) / (3^2 - 1);
    Rfiltpatch(2,2) = 0;


while t <= tmax       % MAIN TIME LOOP ---------------------------------%%

        % Apply elevation boundary condition
% Flat border
   %      topo(:,1:borderwidth) = z_bound;  % WEST
         topo(end-(borderwidth-1):end,:) = z_bound; % min(topo(end - borderwidth,:)); % NORTH %
   %      topo(:,end-(borderwidth-1):end) = z_bound; % EAST
    
        % Base Level
        topo(1:borderwidth,:) = z_bound;  % SOUTH

%         % Wraparound edges
        topo(:,1) = mean([topo(:,2) topo(:,end-1)],2);
        topo(:,end) = topo(:,1);
        
        

        % Don't let anything erode past base level. A bandaid for instability at big timesteps.
        topo(topo < z_bound) = z_bound;     
        
        % Figure out local slopes and flow directions
            [R, slopes] = dem_flow(topo,dx,dy);
            RNANS = isnan(R);   % Find closed depressions for special treatment

        % Apply slope boundary condition (zero slope; to be replaced with constant-slope in the future)
            slopes(1,:) = 0;
            slopes(:,1) = 0;
            slopes(end,:) = 0;
            slopes(:,end) = 0;
            slopes(slopes < 0) = 0;

%         % Wraparound edges
            slopes(:,1) = mean([slopes(:,2) slopes(:,end-1)],2);
            slopes(:,end) = slopes(:,1);
        
        
        % Constant regolith thickness boundaries
            Regolith_H(end-(borderwidth-1):end,:) = 0;  % NORTH
            Regolith_H(1:borderwidth,:) = 0;    % SOUTH
%             Regolith_H(:,end-(borderwidth-1):end) = 0;  % EAST
%             Regolith_H(:,1:borderwidth) = 0;    % WEST
        
        % Wraparound edges
            Regolith_H(:,1) = mean([Regolith_H(:,2) Regolith_H(:,end-1)],2);  % EAST
            Regolith_H(:,end) = Regolith_H(:,1);        % WEST

        % Accumulation areas
            Ts = flow_matrix(topo,R);
            AccGrid = (dx * dy) * upslope_area(topo,Ts);
        
        % Stream widths
            if stream_width_coeff ~= 0
                StreamWs = stream_width_coeff .* realsqrt(AccGrid);
            end
  
        
    %%%% STRATIGRAPHY HANDLER %%%%
    
        % For each defined unit, Read top left and lower right x,y,z-values
        
            RTGrid(:,:) = 0;  % Reset rocktype grid - Default rocktype is zero

            % Make an array of all unit IDs with a top bound above z_bound and
            % a bottom bound below max(topo), then use it, so the surface
            % intersection calculation doesn't need to go through the whole list
            % when there's lots of stratigraphy

            nearsurf = stratlist(7,:) <= max(max(topo)) & stratlist(4,:) >= z_bound;
            unitlist = stratlist(8,nearsurf);
      
        for unit = unitlist
        
        % RTGrid where topo(Xtl:Xlr,Ytl:Ylr) is less than top value and more than bottom value is assigned
        % the rock type index of that stratigraphic unit. If units are
        % defined to overlap, later units supercede earlier.
        rtype_here(:,:) = 0;
        topoabove = topo(stratlist(3,unit):stratlist(6,unit),stratlist(2,unit):stratlist(5,unit)) <= stratlist(4,unit);
        topobelow = topo(stratlist(3,unit):stratlist(6,unit),stratlist(2,unit):stratlist(5,unit)) >= stratlist(7,unit);
        rtype_here(stratlist(3,unit):stratlist(6,unit),stratlist(2,unit):stratlist(5,unit)) = topoabove .* topobelow;
        rtype_here = logical(rtype_here);
        RTGrid(rtype_here) = stratlist(1,unit);
        
        % Update the stratigraphy list z-values given the rock uplift rate
        stratlist(4,unit) = stratlist(4,unit) + (rock_uplift * dt);
        stratlist(7,unit) = stratlist(7,unit) + (rock_uplift * dt);
        
        end % for unit
        
        % Now add in any surface layers tracked
        if DO_ROCKFALL
            RTGrid(RF_Debris_H > RF_Debris_PlotCutoff) = RF_Debris_Rtype;
        end
        
    %%%% PROPERTY HANDLER %%%%
    
        % Establish grids for each property. Loop through units present in RTGrid 
        % and assign property values where RTGrid is each.       
       
        k(:,:) = k0;              % fluvial erodibility constant
        m(:,:) = m0;                % area exponent
        n(:,:) = n0;                % slope exponent
        rdot(:,:) = rdot0;          % regolith production rate, m/yr
        rstar(:,:) = rstar0;        % e-folding depth for falloff of regolith production rate
        rfslope(:,:) = rfslope0;        % rockfall threshold slope
        
        for rtype = 1:numStypes;
        
        prtype_here = (RTGrid == rtype);
            
        k(prtype_here) = RTproplist(2,rtype);              % fluvial erodibility constant
        m(prtype_here) = RTproplist(3,rtype);              % area exponent
        n(prtype_here) = RTproplist(4,rtype);              % slope exponent
        rdot(prtype_here) = RTproplist(5,rtype);           % regolith production rate, m/yr 
        rstar(prtype_here) = RTproplist(6,rtype);           % e-folding depth for falloff of regolith production rate
        rfslope(prtype_here) = RTproplist(7,rtype);        % rockfall threshold slope
        
        end % for prtype
    
    %%%% PROCESS HANDLER %%%%

        % Raw stream power
            stream_power = AccGrid.^m .* (abs(slopes)).^n;
            stream_power(RNANS) = 0;    % Where no downhill neighbors, no stream power!
            
        % Apply fluvial erosion
            Edot_fluv = -k .* stream_power;
            
        % Correct for stream width
            Edot_fluv = Edot_fluv .* (StreamWs./GridDelta);
        
       % Find bedrock channels
            BRChan = stream_power > (sp_susp);

        % Figure out curvatures
            curves = 4*del2(topo,dx,dy);

        % Generate regolith
            rdot_H = rdot .* exp(-Regolith_H ./ rstar);     % Heimsathian
            rdot_H(BRChan) = 0;
            Regolith_H = Regolith_H + rdot_H .* dt;
            
        % Regolith transport
        %    UregDiff = (kappa .* abs(slopes))./((1-abs(slopes)./sc).^2);    % Regolith transport rate, Roering 2001
            UregDiff = (kappa .* abs(slopes)) .* exp(sqrt(abs(slopes))./sc);    % Regolith transport rate, exponential model
            UregFluv =  kt .* AccGrid.^mt .* (abs(slopes)).^nt;    % Regolith transport rate by stream transport
        
        % Correct for stream width
            UregFluv = UregFluv .* (StreamWs / GridDelta);
            
        % Mobile regolith thickness
            Mobile_H = min(Regolith_H,Max_Mobile_H);
            
        % Total regolith flux
            Qreg = abs(UregDiff + UregFluv) .* Mobile_H;          
            
        % Regolith flux boundaries
            Qreg(end-(borderwidth-1):end,:) = 2 * max(Qreg(end - borderwidth-2,:));
            Qreg(1:borderwidth,:) = 2 * max(Qreg(borderwidth+1,:));

        % Wraparound edges
            Qreg(:,1) = mean([Qreg(:,2) Qreg(:,end-1)],2);
            Qreg(:,end) = Qreg(:,1);
            
        % Decompose into x and y components
            
            % First, assign values to NANs in R
              Rfilt = R;
              RRand = rand(y,x) * 2 * pi;
              Rfilt(RNANS) = RRand(RNANS);
              Rfilt = filter2(Rfiltpatch,Rfilt);
              R(RNANS) = Rfilt(RNANS);
        
            % Get topographic gradient for its sign
            [gradX gradY] = gradient(topo,dx,dy);
            %   grads = sqrt(gradX.^2 + gradY.^2);  
            
            QregX = Qreg .* abs(cos(R)) .* -sign(gradX) .* GridDelta; % x component of regolith flux
            QregY = Qreg .* abs(sin(R)) .* -sign(gradY) .* GridDelta; % y component of regolith flux
        
        % Convert to erosion rates
            Edot_reg = -divergence(Xs,Ys,QregX,QregY) / CellArea;
            Edot_reg(isnan(Edot_reg)) = 0;

            
        % Elevation change due to loss in supended transport
            %Mobile_H = Mobile_H + Edot_reg.*dt;
            dzSuspSed = -(stream_power / sp_susp) .* Mobile_H;
            dzSuspSed(dzSuspSed < -Regolith_H) = -Mobile_H(dzSuspSed < -Regolith_H);    % No removal of more than Regolith_H worth of regolith
            
            
        % Old topo diffusion code :)
            %             % Diffuse
            %             Edot_diff = kappa .* curves;
            
        % Subdivide timestep so regolith erosion rates are not applied to
        % bedrock
            dt_prime = (Regolith_H + dzSuspSed) ./ -Edot_reg;  % time needed to strip regolith
            dt_prime(dt_prime > dt) = dt;       % if not stripped, use normal dt
            dt_prime(dt_prime < 0) = dt;         % if negative, then thickness increases, use normal dt
            dt_prime(Edot_reg == 0 & Regolith_H > 0) = dt;        % if no reg erosion, no bedrock erosion
            dt_prime(BRChan) = 0;         % if a bedrock channel, dt_prime is 0
            
        % Rockfall
        if DO_ROCKFALL
            erodeMaxRate = RFerodeRate;
            RF_Spread_lemLink;
        else
            RFdzdt(:,:) = 0;
        end
            
            
        % Integrated erosion
            
            % Elevation change due to regolith motion
            dzRegFrac = dt_prime .* Edot_reg;   

            % Total elevation change for the timestep
            dz = dzRegFrac + dzSuspSed + RFdzdt; 
           
            % Keep track of total elevation change - should equal
            % init topo-final topo at the end of the run...
            dzCum = dzCum + dz;     

            % Regolith thickness change due to fluvial bedrock erosion
            dzFluvFrac = (dt-dt_prime) .* Edot_fluv;
            
        % Update topography
            topo = topo + dz + (dt * rock_uplift);
            topo = fill_sinks(topo);    % Fill DEM sinks
            
        % Update surface layer thicknesses
            Regolith_H = Regolith_H + dz;    % Update regolith thickness for erosion
            Regolith_H = Regolith_H - dzFluvFrac;    % Add fluvial sediment to regolith
            Regolith_H(Regolith_H < 0) = 0;       % No negative thicknesses
            
        % Time of last deposition
            eroded = dz < 0;
            deposited = dz >= 0;
            %Adep(eroded) = 0;
            Adep(deposited) = t;

        % Rockfall layer magic
            if DO_ROCKFALL
                newDebris = RFdzdt >= Mobile_H;
                RF_Debris_H(newDebris) = RF_Debris_H(newDebris) + RFdzdt(newDebris);     % Add new rockfall debris into the debris thickness layer
                RF_Debris_H(eroded) = RF_Debris_H(eroded) + dz(eroded);   % Subtract eroded material
                RF_Debris_H(RF_Debris_H < 0) = 0;       % No negative thicknesses
                RF_Debris_H(newDebris) = RF_Debris_H(newDebris) + Regolith_H(newDebris); % Incorporate existing regolith into rockfall debris (where debris thicker than regolith)
                Regolith_H(newDebris) = 0;          % New rockfall debris wipes out existing regolith (no depositional stratigraphy)
                Regolith_H(~newDebris & RFdzdt > 0) = Regolith_H(~newDebris & RFdzdt > 0) + RFdzdt(~newDebris & RFdzdt > 0);  % Thin RF debris is incorporated into regolith
            end
        
        % Track volumetric erosion rates by each process every 'tracktime' years
            if t_track >= tracktime || t == 0 || t >= tmax

                t_track = 0 + rem(t,tracktime);   % Reset tracking timer such that additional time due to 
                                                % variable dt doesn't accumulate

                trackstep = trackstep + 1;

                TIMEtrackVol(trackstep) = t;
                FLUVtrackVol(trackstep) = (CellArea / dt) * sum(sum(dzFluvFrac));
                DIFFtrackVol(1,trackstep) = (CellArea / dt) * sum(sum(dzRegFrac(dzRegFrac < 0))); % erosion
                DIFFtrackVol(2,trackstep) = (CellArea / dt) * sum(sum(dzRegFrac(dzRegFrac > 0))); % deposition
                RFtrackVol(1,trackstep) = (CellArea / dt) * sum(sum(RFdzdt(RFdzdt < 0)));
                RFtrackVol(2,trackstep) = (CellArea / dt) * sum(sum(RFdzdt(RFdzdt > 0)));
                
                % Check conservation of mass
                Delta_Z = topo_init - (topo - rock_uplift .* t);   % Total elevation change
                Delta_V = sum(sum(Delta_Z)) * CellArea;
                
                Boundary_Loss = (x*dx)*(2*borderwidth*dy) * rock_uplift;
                Mass_Loss = sum(sum(dzCum)) * CellArea;
                
                t
                Mass_Conservation = (-Mass_Loss + Boundary_Loss) / Delta_V;
                
                
            end % if tracktime
        
        % Plot every "plottime" years
            if t_plot >= plottime || t == 0 || t >= tmax

                t_plot = 0 + rem(t,plottime);   % Reset plot timer such that additional time due to 
                                                % variable dt doesn't accumulate

                stype = (RTGrid .* (z_max_plot / numStypes)) + (topo/numStypes);
                figure(1); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; surfc(Xs,Ys,topo,stype); set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp;
                title([run_name ' Topography and Lithology, Year ' int2str(t)]);

                regmap = Regolith_H; 
                RegHCutoff = min( max(max(Regolith_H)), 20 );
                %regmap(AccGrid > 10*dx*dy) = 256;
                figure(2); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([min(min(Regolith_H)) RegHCutoff]); surfc(Xs,Ys,topo,regmap); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp;
                title([run_name ' Regolith Thickness (m), Year ' int2str(t)]);

    %             qregmap = Qreg;
    %             qregmap(isinf(qregmap)) = max(max(qregmap(~isinf(qregmap))));
    %             figure(3); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([min(min(qregmap)) max(max(qregmap))]); surfc(Xs,Ys,topo,qregmap); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp; 
    %             title([run_name ' Regolith Flux (m^3/yr), Year ' int2str(t)]);

%                 depoage = t - Adep;
%                 depoage(depoage < 1) = 1;
%                 depoage = log10(depoage);
%                 daCut = min( max(max(depoage)), inf);
%                 figure(4); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([min(min(depoage)) daCut]); surfc(Xs,Ys,topo,depoage); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp; 
%                 title([run_name ' Log_1_0 time since last deposition (years), Year ' int2str(t)]);

                instedot = dz / dt;
                figure(5); clf; colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; caxis([min(min(instedot)) max(max(instedot))]); surfc(Xs,Ys,topo,instedot); colorbar; set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp; 
                title([run_name ' Instantaneous Erosion Rate (m/yr), Year ' int2str(t)]);
                
                drawnow

                if SAVEMODE == 1      % Save only the figure
                    disp('Saving current figure... ')

                    hgsave(['./' run_filename '/' run_name '_Plot' int2str(stateNo) '.fig'])

                    % Save the entire workspace the first time
                    if t == 0
                        save(['./' run_filename '/' run_name '_State' int2str(stateNo) '.mat']);
                    end

                    stateNo = stateNo + 1;

                    disp('Saved.')

                elseif SAVEMODE == 2    % Save the entire workspace
                    disp('Saving current state... ')

                    % Save the workspace
                    save(['./' run_filename '/' run_name '_State' int2str(stateNo) '.mat']);
                    stateNo = stateNo + 1;

                    disp('Saved.')

                end   % if SAVEMODE

                disp('Max erosion over last timestep: ')
                min(min(dz))
                disp('Max deposition over last timestep: ')
                max(max(dz))
                disp('Timestep: ')
                dt %#ok<NOPTS>
                disp('Time: ')
                t %#ok<NOPTS>
                
                if NoDistTargets ~= 0
                    disp([int2str(NoDistTargets) ' rockfall events had no distribution targets.'])
                    NoDistTargets = 0;
                end
                
            end % if t_plot
        
        % Figure out new timestep

            max_e = abs(dz / dt);   % This timestep's erosion
            tauE = min(min(rstar ./ max_e));    % Don't erode more than one rstar worth of material!
            tauReg = min(min(rstar ./ rdot_H));  % Don't make more than one rstar worth of regolith!

            dt = 0.1 * min(tauE, tauReg);
            
        % Double-check against diffusion stability
            dtDM = (0.25 * GridDelta^2) ./ (kappa .* slopes.^2 .* exp(slopes/sc));
            dtDMm = 0.1 * min(min(dtDM));
            
            dt = min(dt, dtDMm);
            
        % Double-check against user timestep constraints    
            if dt < dt_min; dt = dt_min; end
            if dt > dt_max; dt = dt_max; end
        
        t = t + dt;             % Update time
        t_plot = t_plot + dt;   % Update plot time counter
        t_track = t_track + dt; % Update track time counter
        
end                % END MAIN TIME LOOP ---------------------------------%%

% Save landscape matrix and specified variables to the run's subfolder
% -------------------------------------------------------------------------

disp('Saving final state... ')

% Save the landscape matrix
save(['./' run_filename '/' run_name '_EndState.mat']);

disp('Saved.')

if SAVEMODE == 1 || SAVEMODE == 2

    % Prompt to run MakeLEMmingMov.m
    makemovie = input('Do you want to make a movie of this run? Y/N [Y]: ', 's');
    if isempty(makemovie) || makemovie == 'y';
        makemovie = 'Y';
    end

else
    disp('Done!')
end

if makemovie == 'Y'
    MakeLEMmingMov
else
    disp('Done!')
end

%%%%%%%% END OF LEMming.m %%%%%%%%