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

% Stochastic rockfall runout simulator. This version gets run inline with
% LEMming v015 and above.

RFdzdt(:,:) = 0;

sourceMask = curves < RFSource_curv & slopes > rfslope;

[sourceYs,sourceXs] = find(sourceMask);
sourceZs = topo(sourceMask);

n_sources = length(sourceZs);
n_events_do = n_sources;

eventList = ceil( rand(1,n_events_do) .* (n_sources-.1) );
[eventListFilt uniqueEventIdx] = unique(eventList);  % Screen for dupes

erodeBucket = erodeMaxRate * dt;    % Figure out max erosion allowed in this timestep

% Now simulate the collection of events
if ~isempty(eventList);
    
    % Steps to compute aspect from surface normals
        [Xnorm Ynorm Znorm] = surfnorm(topo);
        AspectDeg=atand(abs(Xnorm)./abs(Ynorm));
        AspectDeg(Xnorm < 0 & Ynorm > 0) = 360-AspectDeg(Xnorm < 0 & Ynorm > 0);
        AspectDeg(Xnorm < 0 & Ynorm < 0) = 180+AspectDeg(Xnorm < 0 & Ynorm < 0);
        AspectDeg(Xnorm > 0 & Ynorm < 0) = 180-AspectDeg(Xnorm > 0 & Ynorm < 0);
    
    % Debris distribution probability functions not dependent on event
    % location
%         curvtemp = curves * curvStar;
%         curvtemp(curvtemp > 500) = 500;     % On DW's system, e^709 is finite but e^710 and up returns Inf. This threshold prevents NaNts later in the calculation.
        cfilt = filter2(filtpatch,curves);
        
        curvProb = ( erf(cfilt ./ (curvStar/e)) + 1 ).^2;
        curvProb = curvProb ./ max(max(curvProb));

        sfilt = slopes; %filter2(filtpatch,slopes);
        
        slopeProb = sDepCrit^slopeExp - sfilt.^slopeExp;
        slopeProb(slopeProb < 0) = 0;
        slopeProb = slopeProb ./ max(max(slopeProb));
    
        
    for event = uniqueEventIdx
        if erodeVol < erodeBucket;
            
        sourceX = sourceXs(event);
        sourceY = sourceYs(event);
        sourceZ = sourceZs(event);

        % Erosion. Maintain relief; subtract source cell to height of highest cell
        % below it.
        event_rows = sourceY-1:sourceY+1;
            event_rows(event_rows > y) = event_rows(event_rows > y) - y; % Wraparound boundary
            event_rows(event_rows < 1) = event_rows(event_rows < 1) + y; % Wraparound boundary
        event_cols = sourceX-1:sourceX+1;
            event_cols(event_cols > x) = event_cols(event_cols > x) - x; % Wraparound boundary
            event_cols(event_cols < 1) = event_cols(event_cols < 1) + x; % Wraparound boundary
        
        localDrops = sourceZ - topo(event_rows,event_cols);
        erodeHt = min(min( localDrops(localDrops > 0) )); % could use mean of these.
        if isempty(erodeHt); erodeHt = 0; end
        
        
        % Deposition - convolve probability distributions to determine
        % deposition patterns
        IsBelow = topo < (sourceZ-(1*erodeHt)); % Don't allow deposition above the base of the failure
        
        dists = sqrt((abs(Xs-sourceX*dx).^2) + (abs(Ys-sourceY*dy)).^2);

        distProb = exp(-dists ./ distStar);
        distProb(dists > distMax) = 0;
        distProb = distProb ./ max(max(distProb));

        % Only deposit in the viewshed (needs Mapping Toolbox)        
        % Could use a small sphere radius to filter outlying paths...
%         CPD = 111.11e3 / GridDelta;    % cells/degree
%         llLat = 39;     % Lower left cell north latitude
%         llLong = 108;  % Lower left cell west longitude
%         [visZ,visrefvec] = viewshed(filter2(filtpatch,topo),[CPD llLat llLong],llLat - (y - sourceY)/CPD,llLong+sourceX/CPD, (0.5 * GridDelta * slopes(sourceY,sourceX)));
%         viewProb = filter2(filtpatch,visZ);  % Filter viewshed, 'cause things roll
         XDiffs = Xs-sourceX*dx;
         YDiffs = Ys-sourceY*dy;
        % Angular distances from source
         AngDists = atand(XDiffs ./ YDiffs);
         AngDists(XDiffs < 0 & YDiffs >= 0) = 360+AngDists(XDiffs < 0 & YDiffs >= 0);
         AngDists(XDiffs < 0 & YDiffs < 0) = 180+AngDists(XDiffs < 0 & YDiffs < 0);
         AngDists(XDiffs >= 0 & YDiffs < 0) = 180+AngDists(XDiffs >= 0 & YDiffs < 0);
%          
     % Angles between source azimuth and each pixel
        PixDepoAngs = AspectDeg(sourceY,sourceX) - AngDists;
        PixDepoAngs = min(abs(PixDepoAngs),abs(360 - AngDists + AspectDeg(sourceY,sourceX)));
     % Compute angular probabilities
        angProb = PixDepoAngs <= DepoAngleCutoff;
        
        
         cumEventProb = angProb .* distProb .* IsBelow .* slopeProb .* curvProb;       
%        cumEventProb = viewProb .* distProb .* IsBelow .* slopeProb .* curvProb;
%        cumEventProb = IsBelow .* distProb .* (slopeProb + curvProb);
%        cumEventProb(distProb <= 0 | slopeProb <= 0 | curvProb <= 0) = 0;
        
        cumEventTotal = sum(sum(cumEventProb));

        if cumEventTotal > 0;
            cumEventProb = cumEventProb ./ cumEventTotal;
        else
            NoDistTargets = NoDistTargets + 1;
        end
             
        % Update erosion/deposition amount grid
        RFdzdt = RFdzdt + (cumEventProb .* erodeHt); % no volume necessary, because cells all the same size
        RFdzdt(sourceY, sourceX) = RFdzdt(sourceY, sourceX) - erodeHt;

        % Volume tracking for regulation of erosion rate
        eventVol = erodeHt * CellArea;
        erodeVol = erodeVol + eventVol;
        
        else break % the event loop
        end  % if erodebucket
        
    end  % for event

end % if ~isempty(eventList)

erodeVol = erodeVol - erodeBucket; % get excess erosion
erodeVol = max(erodeVol,0); % reset to zero if no leftovers