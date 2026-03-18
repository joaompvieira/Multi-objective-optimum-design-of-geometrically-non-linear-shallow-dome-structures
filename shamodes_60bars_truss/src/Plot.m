%% Plot Class
%
% This class defines a plotting object in the NUMA-TF program.
% A plotting object is responsible for providing the analysis results in
% graphical format in Matlab figure windows.
%
classdef Plot < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        sim       = [];  % object of the simulation class
        modelName = [];  % model name for figure titles
        
        % Matlab figures identification names
        figName_model     = 'Model';
        figName_mode      = 'Mode';
        figName_eqPath    = 'Equilibrium Path';
        figName_timeAccel = 'Time_x_Acceleration';
        figName_timeVeloc = 'Time_x_Velocity';
        figName_timeDispl = 'Time_x_Displacement';
        figName_freqAccel = 'Freq_x_Acceleration';
        figName_freqVeloc = 'Freq_x_Velocity';
        figName_freqDispl = 'Freq_x_Displacement';
        figName_phasePort = 'Phase Portrait';
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function plot = Plot(sim,fin)
            % Set simulation object
            plot.sim = sim;
            
            % Plot results
            if (sim.opt.plot)
                % Get model name
                [~,filename,~] = fileparts(fopen(fin));
                plot.modelName = filename;
                
                % Execute plotting
                plot.execute();
            end
        end
    end
    
    %% Public methods
    methods
        %% Main function
        %------------------------------------------------------------------
        % Execute plotting process according to analysis type.
        function execute(plot)
            c = Constants();            
            if (plot.sim.anl.type == c.LINEAR_ELASTIC)
                % Initial configuration / Deformed configuration / View limits
                plot.drawModel(plot.sim.mdl);
                plot.drawDeformed(plot.sim.mdl,plot.sim.anl,plot.sim.res);
                plot.adjustLimits(plot.sim.anl,plot.sim.res)
            elseif (plot.sim.anl.type == c.NONLINEAR_GEOM)
                % Equilibrium path / Initial configuration / Deformed configuration / View limits
                plot.plotEquilibriumPath(plot.sim.mdl,plot.sim.res);
                plot.drawModel(plot.sim.mdl);
                plot.drawDeformed(plot.sim.mdl,plot.sim.anl,plot.sim.res);
                plot.adjustLimits(plot.sim.anl,plot.sim.res)
            elseif (plot.sim.anl.type == c.VIBRATION_FREE)
                % Phase portrait / Time response / Frequency response / initial configuration / View limits
                plot.plotPhasePortrait(plot.sim.mdl,plot.sim.res);
                plot.plotTimeResponse(plot.sim.mdl,plot.sim.res);
                plot.plotFreqResponse(plot.sim.mdl,plot.sim.anl,plot.sim.res);
                plot.drawModel(plot.sim.mdl);
                plot.adjustLimits(plot.sim.anl,plot.sim.res)
            elseif (plot.sim.anl.type == c.MODAL_BUCKLING || plot.sim.anl.type == c.MODAL_VIBRATION)
                % Modes / View limits
                plot.drawModes(plot.sim.mdl,plot.sim.res);
                plot.adjustLimits(plot.sim.anl,plot.sim.res)
            end
        end
        
        %% Functions for plotting results
        %------------------------------------------------------------------
        % Draw original configuration of the structural model.
        function drawModel(self,mdl)
            c = Constants();
            if (self.sim.opt.overlay)
                figname = self.figName_model;
            else
                figname = [self.figName_model,' ',self.modelName];
            end
            if (isempty(findall(0,'type','figure','name',figname)))
                figure('Name',figname);
            else
                set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
            end
            hold on
            grid off
            axis equal
            for i = 1:mdl.nel
                x1 = mdl.elems(i).nodes(1).coord(1);
                y1 = mdl.elems(i).nodes(1).coord(2);
                z1 = mdl.elems(i).nodes(1).coord(3);
                x2 = mdl.elems(i).nodes(2).coord(1);
                y2 = mdl.elems(i).nodes(2).coord(2);
                z2 = mdl.elems(i).nodes(2).coord(3);
                x = [x1 x2];
                y = [y1 y2];
                z = [z1 z2];
                scatter3(x,y,z,10,'filled','black');
                plot3(x,y,z,'black');
            end
            if (mdl.anm.type == c.TRUSS3D || mdl.anm.type == c.FRAME3D)
                view(45,45)
                grid on
            end
        end
        
        %------------------------------------------------------------------
        % Draw deformed configuration.
        function drawDeformed(self,mdl,anl,res)
            c = Constants();
            if (self.sim.opt.overlay)
                figname = self.figName_model;
            else
                figname = [self.figName_model,' ',self.modelName];
            end
            if (isempty(findall(0,'type','figure','name',figname)))
                figure('Name',figname);
            else
                set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
            end
            hold on
            axis equal
            
            % Scale factor for linear analysis
            if (anl.type == c.LINEAR_ELASTIC)
                maxDisp = 0;
                maxSize = 0;
                % Compute maximum nodal displacement
                for i = 1:mdl.neqf
                    localDOFnum = find(mdl.ID == i);
                    if ( (mdl.anm.type == c.FRAME2D && localDOFnum == 3) ||...
                         (mdl.anm.type == c.FRAME3D && localDOFnum >  3))
                        continue;
                    else
                        if (abs(res.U(i)) > abs(maxDisp))
                            maxDisp = abs(res.U(i));
                        end  
                    end
                end
                % Compute maximum element length
                for i = 1:mdl.nel
                    if (mdl.elems(i).len_ini > maxSize)
                        maxSize = mdl.elems(i).len_ini;
                    end
                end
                fct = maxSize / (10*maxDisp);
            else
                fct = 1;
            end
            
            size = res.steps+1;
            
            for i = 1:mdl.nel
                if (mdl.anm.type == c.TRUSS2D)
                    x1 = mdl.elems(i).nodes(1).coord(1) + fct*res.U(mdl.elems(i).gle(1),size);
                    y1 = mdl.elems(i).nodes(1).coord(2) + fct*res.U(mdl.elems(i).gle(2),size);
                    z1 = 0.0;
                    x2 = mdl.elems(i).nodes(2).coord(1) + fct*res.U(mdl.elems(i).gle(3),size);
                    y2 = mdl.elems(i).nodes(2).coord(2) + fct*res.U(mdl.elems(i).gle(4),size);
                    z2 = 0.0;
                elseif (mdl.anm.type == c.FRAME2D)
                    x1 = mdl.elems(i).nodes(1).coord(1) + fct*res.U(mdl.elems(i).gle(1),size);
                    y1 = mdl.elems(i).nodes(1).coord(2) + fct*res.U(mdl.elems(i).gle(2),size);
                    z1 = 0.0;
                    x2 = mdl.elems(i).nodes(2).coord(1) + fct*res.U(mdl.elems(i).gle(4),size);
                    y2 = mdl.elems(i).nodes(2).coord(2) + fct*res.U(mdl.elems(i).gle(5),size);
                    z2 = 0.0;
                elseif (mdl.anm.type == c.TRUSS3D)
                    x1 = mdl.elems(i).nodes(1).coord(1) + fct*res.U(mdl.elems(i).gle(1),size);
                    y1 = mdl.elems(i).nodes(1).coord(2) + fct*res.U(mdl.elems(i).gle(2),size);
                    z1 = mdl.elems(i).nodes(1).coord(3) + fct*res.U(mdl.elems(i).gle(3),size);
                    x2 = mdl.elems(i).nodes(2).coord(1) + fct*res.U(mdl.elems(i).gle(4),size);
                    y2 = mdl.elems(i).nodes(2).coord(2) + fct*res.U(mdl.elems(i).gle(5),size);
                    z2 = mdl.elems(i).nodes(2).coord(3) + fct*res.U(mdl.elems(i).gle(6),size);
                elseif (mdl.anm.type == c.FRAME3D)
                    x1 = mdl.elems(i).nodes(1).coord(1) + fct*res.U(mdl.elems(i).gle(1),size);
                    y1 = mdl.elems(i).nodes(1).coord(2) + fct*res.U(mdl.elems(i).gle(2),size);
                    z1 = mdl.elems(i).nodes(1).coord(3) + fct*res.U(mdl.elems(i).gle(3),size);
                    x2 = mdl.elems(i).nodes(2).coord(1) + fct*res.U(mdl.elems(i).gle(7),size);
                    y2 = mdl.elems(i).nodes(2).coord(2) + fct*res.U(mdl.elems(i).gle(8),size);
                    z2 = mdl.elems(i).nodes(2).coord(3) + fct*res.U(mdl.elems(i).gle(9),size);
                end
                x = [x1 x2];
                y = [y1 y2];
                z = [z1 z2];
                plot3(x,y,z,'Color','red','LineWidth',0.5);
            end
            if (anl.type == c.LINEAR_ELASTIC)
                title(['Deformed configuration: Scale fator ',num2str(fct)]);
            else
                title('Deformed configuration');
            end
        end
        
        %------------------------------------------------------------------
        % Draw buckling or vibration modes with given scale factor.
        function drawModes(self,mdl,res)
            c = Constants();
            num_pos = 50; % Number of positions to discretize each element draw
            for i = res.nmodes:-1:1
                if (self.sim.opt.overlay)
                    figname = [self.figName_mode,' ',int2str(i)];
                else
                    figname = [self.figName_mode,' ',int2str(i),' - ',self.modelName];
                end
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                end
                hold on
                grid off
                axis equal
                title([self.figName_mode,int2str(res.mode(i))]);
                
                for j = 1:mdl.nel
                    % Node coordinates
                    x1 = mdl.elems(j).nodes(1).coord(1);
                    y1 = mdl.elems(j).nodes(1).coord(2);
                    z1 = mdl.elems(j).nodes(1).coord(3);
                    x2 = mdl.elems(j).nodes(2).coord(1);
                    y2 = mdl.elems(j).nodes(2).coord(2);
                    z2 = mdl.elems(j).nodes(2).coord(3);
                    
                    % Draw initial configuration
                    scatter3([x1,x2],[y1,y2],[z1,z2],10,'filled','black');
                    plot3([x1,x2],[y1,y2],[z1,z2],'black','LineWidth',0.9);
                    
                    % Orientation cosines
                    L = mdl.elems(j).len_cur;
                    cx = (x2 - x1)/L;
                    cy = (y2 - y1)/L;
                    cz = (z2 - z1)/L;
                    
                    % Basis rotation matrix
                    if (mdl.anm.type == c.TRUSS2D || mdl.anm.type == c.FRAME2D)
                        T = mdl.elems(j).T(1:2,1:2);
                    elseif (mdl.anm.type == c.TRUSS3D || mdl.anm.type == c.FRAME3D)
                        T = mdl.elems(j).T;
                        view(45,45)
                        grid on
                    end
                    
                    % Nodal displacements at element ends in local systems
                    rot = mdl.elems(j).rot;
                    ul  = rot * res.U(mdl.elems(j).gle,res.mode(i));
                    
                    % Compute element displacements in each internal position and draw mode
                    stride = L / (num_pos-1);
                    for x = 0:stride:L
                        % Compute internal displacements of current position in global system
                        N  = mdl.anm.shapeFcnMtx(mdl.elems(j),x);
                        ui = T' * N * ul;
                        
                        % Global coordinates of current position on deformed configuration
                        xd = (x1 + x * cx) + res.scale(i) * ui(1);
                        yd = (y1 + x * cy) + res.scale(i) * ui(2);
                        if (mdl.anm.type == c.TRUSS2D || mdl.anm.type == c.FRAME2D)
                            zd = 0;
                        elseif (mdl.anm.type == c.TRUSS3D || mdl.anm.type == c.FRAME3D)
                            zd = (z1 + x * cz) + res.scale(i) * ui(3);
                        end
                        
                        % Connect previous position to current
                        if (x == 0)
                            xi = x1 + res.scale(i) * ui(1);
                            yi = y1 + res.scale(i) * ui(2);
                            if (mdl.anm.type == c.TRUSS2D || mdl.anm.type == c.FRAME2D)
                                zi = 0;
                            elseif (mdl.anm.type == c.TRUSS3D || mdl.anm.type == c.FRAME3D)
                                zi = z1 + res.scale(i) * ui(3);
                            end
                        end
                        line([xi,xd],[yi,yd],[zi,zd],'Color','red','LineWidth',0.4);
                        
                        % Update previous cross-section coordinates
                        xi = xd;
                        yi = yd;
                        zi = zd;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Adjust model view limits.
        function adjustLimits(self,anl,res)
            c = Constants();
            if (anl.type == c.LINEAR_ELASTIC || anl.type == c.NONLINEAR_GEOM || anl.type == c.VIBRATION_FREE)
                num_figs = 1;
            elseif (anl.type == c.MODAL_BUCKLING || anl.type == c.MODAL_VIBRATION)
                num_figs = res.nmodes;
            end
            for i = 1:num_figs
                if (anl.type == c.LINEAR_ELASTIC || anl.type == c.NONLINEAR_GEOM || anl.type == c.VIBRATION_FREE)
                    if (self.sim.opt.overlay)
                        figname = self.figName_model;
                    else
                        figname = [self.figName_model,' ',self.modelName];
                    end
                elseif (anl.type == c.MODAL_BUCKLING || anl.type == c.MODAL_VIBRATION)
                    if (self.sim.opt.overlay)
                        figname = [self.figName_mode,' ',int2str(i)];
                    else
                        figname = [self.figName_mode,' ',int2str(i),' - ',self.modelName];
                    end
                end
                if (isempty(findall(0,'type','figure','name',figname)))
                    continue;
                else
                    set(0,'CurrentFigure', findall(0,'type','figure','name',figname));
                end
                axis equal
                axis auto
                axis manual
                xLimits = get(gca,'XLim');
                yLimits = get(gca,'YLim');
                zLimits = get(gca,'ZLim');
                gapX    = (xLimits(2)-xLimits(1))/10;
                gapY    = (yLimits(2)-yLimits(1))/10;
                gapZ    = (zLimits(2)-zLimits(1))/10;
                minX    = xLimits(1) - gapX;
                maxX    = xLimits(2) + gapX;
                minY    = yLimits(1) - gapY;
                maxY    = yLimits(2) + gapY;
                minZ    = zLimits(1) - gapZ;
                maxZ    = zLimits(2) + gapZ;
                xlim([minX maxX]);
                ylim([minY maxY]);
                zlim([minZ maxZ]);
            end
        end
        
        %------------------------------------------------------------------
        % Plot equilibrium paths (load ratio vs. displacements).
        function plotEquilibriumPath(self,mdl,res)
            if (res.ncurves > 0)
                if (self.sim.opt.overlay)
                    figname = self.figName_eqPath;
                else
                    figname = [self.figName_eqPath,' ',self.modelName];
                end
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca,'children'),'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    x = -res.U(d,1:res.steps+1);
                    y = (res.lbd(1:res.steps+1))*44500;
                    plot(x,y);
                end
                xlabel('Displacement')
                ylabel('Load Ratio')
                legend(leg);
            end
        end
        
        %------------------------------------------------------------------
        % Plot time response (displacement, velocity, acceleration).
        function plotTimeResponse(self,mdl,res)
            c = Constants();
            if (res.ncurves > 0 && res.domain == c.TIME)
                % Plot acceleration
                if (self.sim.opt.overlay)
                    figname = self.figName_timeAccel;
                else
                    figname = [self.figName_timeAccel,' ',self.modelName];
                end
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    x = res.times;
                    y = res.A(d,:);
                    plot(x,y);
                end
                title('Time vs. Acceleration')
                xlabel('Time')
                ylabel('Acceleration')
                legend(leg);
                
                % Plot velocity
                figname = [self.figName_timeVeloc,' ',self.modelName];
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    x = res.times;
                    y = res.V(d,:);
                    plot(x,y);
                end
                title('Time vs. Velocity')
                xlabel('Time')
                ylabel('Velocity')
                legend(leg);
                
                % Plot displacement
                figname = [self.figName_timeDispl,' ',self.modelName];
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    x = res.times;
                    y = res.U(d,:);
                    plot(x,y);
                end
                title('Time vs. Displacement')
                xlabel('Time')
                ylabel('Displacement')
                legend(leg);
            end
        end
        
        %------------------------------------------------------------------
        % Plot frequency response (displacement, velocity, acceleration).
        function plotFreqResponse(self,mdl,anl,res)
            c = Constants();
            if (res.ncurves > 0 && res.domain == c.FREQUENCY)
                % Plot Acceleration
                if (self.sim.opt.overlay)
                    figname = self.figName_freqAccel;
                else
                    figname = [self.figName_freqAccel,' ',self.modelName];
                end
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d   = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    pts = 2^nextpow2(length(res.times));
                    D2  = abs(res.A(d,:)/length(res.times));
                    D1  = 2 * D2(1:pts/2+1);
                    fS  = 1/anl.increment;
                    fD  = 0 : fS/pts : fS*(1/2-1/pts);
                    maxD = max(max(D1));
                    minD = min(min(D1));
                    if (minD == maxD)
                        maxD = maxD + 1;
                        minD = minD - 1;
                    elseif (isempty(minD) || isempty(maxD))
                        maxD =  1;
                        minD = -1;
                    else
                        medD = (maxD + minD) * 0.5;
                        intD = maxD - minD;
                        maxD = medD + 0.6*intD;
                        minD = medD - 0.6*intD;
                    end
                    plot(fD,D1(1:pts/2));
                end
                ylim([minD, maxD])
                title('Frequency of Acceleration')
                xlabel('f(Hz)')
                ylabel(' ')
                legend(leg);
                
                % Plot Velocity
                figname = [self.figName_freqVeloc,' ',self.modelName];
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d   = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    pts = 2^nextpow2(length(res.times));
                    D2  = abs(res.V(d,:)/length(res.times));
                    D1  = 2 * D2(1:pts/2+1);
                    fS  = 1/anl.increment;
                    fD  = 0 : fS/pts : fS*(1/2-1/pts);
                    maxD = max(max(D1));
                    minD = min(min(D1));
                    if (minD == maxD)
                        maxD = maxD + 1;
                        minD = minD - 1;
                    elseif (isempty(minD) || isempty(maxD))
                        maxD =  1;
                        minD = -1;
                    else
                        medD = (maxD + minD) * 0.5;
                        intD = maxD - minD;
                        maxD = medD + 0.6*intD;
                        minD = medD - 0.6*intD;
                    end
                    plot(fD,D1(1:pts/2));
                end
                ylim([minD, maxD])
                title('Frequency of Velocity')
                xlabel('f(Hz)')
                ylabel(' ')
                legend(leg);
                
                % Plot displacement
                figname = [self.figName_freqDispl,' ',self.modelName];
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d   = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    pts = 2^nextpow2(length(res.times));
                    D2  = abs(res.U(d,:)/length(res.times));
                    D1  = 2 * D2(1:pts/2+1);
                    fS  = 1/anl.increment;
                    fD  = 0 : fS/pts : fS*(1/2-1/pts);
                    maxD = max(max(D1));
                    minD = min(min(D1));
                    if (minD == maxD)
                        maxD = maxD + 1;
                        minD = minD - 1;
                    elseif (isempty(minD) || isempty(maxD))
                        maxD =  1;
                        minD = -1;
                    else
                        medD = (maxD + minD) * 0.5;
                        intD = maxD - minD;
                        maxD = medD + 0.6*intD;
                        minD = medD - 0.6*intD;
                    end
                    plot(fD,D1(1:pts/2));
                end
                ylim([minD, maxD])
                title('Frequency of Displacement')
                xlabel('f(Hz)')
                ylabel(' ')
                legend(leg);
            end
        end
        
        %------------------------------------------------------------------
        % Plot time response (displacement, velocity, acceleration).
        function plotPhasePortrait(self,mdl,res)
            c = Constants();
            if (res.ncurves > 0 && res.domain == c.TIME)
                if (self.sim.opt.overlay)
                    figname = self.figName_phasePort;
                else
                    figname = [self.figName_phasePort,' ',self.modelName];
                end
                if (isempty(findall(0,'type','figure','name',figname)))
                    figure('Name',figname);
                    leg = res.name;
                else
                    set(0, 'CurrentFigure', findall(0,'type','figure','name',figname));
                    strings = flipud(get(get(gca, 'children'), 'displayname'));
                    leg = [strings;res.name'];
                end
                hold on
                grid on
                for i = 1:res.ncurves
                    d = mdl.ID(mdl.anm.gla==res.dof(i),res.node(i));
                    x = res.U(d,:);
                    y = res.V(d,:);
                    z = res.A(d,:);
                    plot3(x,y,z);
                end
                title('Phase Portrait')
                xlabel('Displacement')
                ylabel('Velocity')
                zlabel('Acceleration')
                legend(leg);
                view(45,45)
            end
        end
    end
end