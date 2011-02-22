% MakeLEMmingMov.m - restores each saved step in the LEMming model run, plots it,
% captures the frame, and makes an .avi movie from it.

clc
close all

FigWidth = 800; % Specify the desired movie size
FigHeight = 600;    % in pixels
nth_frame = 1;      % integer, plot every nth_frame
colordef none

if exist('M','var'); clear M; end  % Reinitialize the movie matrix

if ~exist('run_name','var')
    run_name = ' ';      % Copy the run name here or open any statefile from the run before running MakeLEMmingMov.m
    run_filename = ' ';  % Copy the run filename (Folder name) here or "
end

% Load the initial state 
try     % in case something isn't right in the run folder
    
    % Load the final workspace and read the SAVEMODE and final
    % stateNo
    load(['./' run_filename '/' run_name '_EndState.mat'],'SAVEMODE','stateNo');

    stateNo = stateNo - 1; % Because it would have been incremented past the final state
    
    frame = 1;  % initial frame
    
    % Loop through states and load each file
    for state = 0:nth_frame:stateNo
                                
            if SAVEMODE == 1      % Saved only the figures

                close all
                open(['./' run_filename '/' run_name '_Plot' int2str(state) '.fig'])
                set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1])
                           
            elseif SAVEMODE == 2    % Saved the entire workspaces


                % Load the workspace
                load(['./' run_filename '/' run_name '_State' int2str(state) '.mat']);

                % Remake the plot
                stype = (RTGrid .* (z_max_plot / numStypes)) + (topo/numStypes);
                figure(1); clf; set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]); colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; surfc(Xs,Ys,topo,stype); set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp;
                title([run_name ' Topography and Lithology, Year ' int2str(t)]);

%                 regmap = (Regolith_H ./ max(max(Regolith_H)) * max(max(topo)));
%                 %regmap(AccGrid > 10*dx*dy) = 256;
%                 figure(1); clf; set(gcf,'Position',[1 1 FigWidth+1 FigHeight+1]); colormap jet; contour3(Xs,Ys,topo,cinterval:cinterval:z_max_plot,'k'); hold on; surfc(Xs,Ys,topo,regmap); set(gca,'DataAspectRatio',[1 1 1/VE]); view(-50,50); zlim([z_bound z_max_plot]); shading interp;
%                 title([run_name ' Regolith Thickness (m), Year ' int2str(t)]);

                drawnow
                % ADDITIONAL PLOTS we may wish to make
                %             figure(2); imagesc(dz); view(2); shading flat; colorbar;
                %             figure(3); imagesc(slopes); view(2); shading flat; colorbar;
                %             figure(4); imagesc(curves); view(2); shading flat; colorbar;
                %             figure(5); imagesc(AccGrid); view(2); shading flat; colorbar;
                %             drawnow
                
            end   % if SAVEMODE

        % Get the movie frame from the state. If multiple plots, this code
        % needs to be put inside the IF SAVEMODE bit after each plot
        % command and separate M variables assigned.
        pause(1); % Short pause to make sure everything draws before the fram is captured
        M(frame) = getframe(gcf); %#ok<AGROW>
        frame = frame + 1;

    end % for state

    disp 'Saving movie...'
    MOVERROR = 0;
    try
        movie2avi(M,['./' run_filename '/' run_name '.avi'],'fps',10)
    catch
        MOVERROR = 1;
    end

    if ~MOVERROR
        % Prompt to toggle deletion of intermediate files to free disk space. Does not
        % remove Year0 or EndState.
        DO_CLEANUP = input('Remove all but the initial and final state savefiles? Y/N [Y]: ', 's');
        if isempty(DO_CLEANUP) || DO_CLEANUP == 'y';
            DO_CLEANUP = 'Y';
        end

        if DO_CLEANUP == 'Y'
           for state = 1:stateNo 
               if SAVEMODE == 1
                    delete(['./' run_filename '/' run_name '_Plot' int2str(state) '.fig'])
               elseif SAVEMODE ==2
                    delete(['./' run_filename '/' run_name '_State' int2str(state) '.mat'])
               end
           end

        end
    end
    disp 'Done!'
catch
        disp('One or more required files missing or run name not specified. Movie not made.')
end

if MOVERROR
    disp 'Error converting movie to AVI. Movie not saved.'
end

% END of MakeLEMmingMovie.m