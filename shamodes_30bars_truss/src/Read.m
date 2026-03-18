%% Read Class
%
% This class defines a reader object in the NUMA-TF program.
% A reader object is responsible for reading the input file with model/analysis
% information and store it in the program data structure.
%
classdef Read < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        status = 1;   % flag for read success
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function read = Read(sim,fin)
            if (nargin > 0)
                read.inputFile(sim,fin);
            end
        end
    end
    
    %% Public methods
    methods
        %% Main function
        %------------------------------------------------------------------
        % Main function to read input file tags.
        function inputFile(read,sim,fin)
            % Create Model and Result objects
            sim.mdl = Model();
            sim.res = Result();
            
            while read.status == 1 && ~feof(fin)
                % Get file line without blank spaces
                string = deblank(fgetl(fin));
                
                % Look for tag strings
                switch string
                    % Model tags
                    case '%MODEL.ANALYSIS_MODEL'
                        read.analysisModel(fin,sim.mdl);
                    case '%MODEL.NODE_COORD'
                        read.nodeCoord(fin,sim.mdl);
                    case '%MODEL.NODE_SUPPORT'
                        read.nodeSupport(fin,sim.mdl);
                    case '%MODEL.NODE_PRESC_DISPL'
                        read.nodePrescDispl(fin,sim.mdl);
                    case '%MODEL.INIT_CONDITION'
                        read.nodeInitCondition(fin,sim.mdl);
                    case '%MODEL.MATERIAL_PROPERTY'
                        read.materialProperty(fin,sim.mdl);
                    case '%MODEL.SECTION_PROPERTY'
                        read.sectionProperty(fin,sim.mdl);
                    case '%MODEL.ELEMENT'
                        read.element(fin,sim.mdl);
                    case '%MODEL.LOAD.NODE_LOAD'
                        read.nodeLoad(fin,sim.mdl);
                    case '%MODEL.LOAD.ELEM_UNIFORM'
                        read.elemUnifLoad(fin,sim.mdl);
                    case '%MODEL.LOAD.ELEM_LINEAR'
                        read.elemLinearLoad(fin,sim.mdl);
                    case '%MODEL.LOAD.ELEM_THERMAL'
                        read.elemThermalLoad(fin,sim.mdl);
                    
                    % Analysis tags
                    case '%ANALYSIS.TYPE'
                        read.analysisType(fin,sim);
                    case '%ANALYSIS.TANGENT_STIFFNESS'
                        read.tangStiff(fin,sim.anl);
                    case '%ANALYSIS.MASS_MATRIX'
                        read.massMtx(fin,sim.anl);
                    case '%ANALYSIS.DAMPING'
                        read.damping(fin,sim.anl);
                    case '%ANALYSIS.METHOD_NONLINEAR'
                        read.methodNonlin(fin,sim.anl);
                    case '%ANALYSIS.METHOD_VIBRATION'
                        read.methodVibration(fin,sim.anl);
                    case '%ANALYSIS.INCREMENT_TYPE'
                        read.incrementType(fin,sim.anl);
                    case '%ANALYSIS.ITERATION_TYPE'
                        read.iterationType(fin,sim.anl);
                    case '%ANALYSIS.PARAMETERS_NONLINEAR'
                        read.analysisParamNonlin(fin,sim.anl);
                    case '%ANALYSIS.PARAMETERS_VIBRATION'
                        read.analysisParamVib(fin,sim.anl);
                        
                    % Result tags
                    case '%ANALYSIS.RESULT_DOMAIN'
                        read.resultDomain(fin,sim.res);
                    case '%ANALYSIS.RESULT_CURVES'
                        read.resultCurves(fin,sim.res,sim.mdl);
                    case '%ANALYSIS.RESULT_MODES'
                        read.resultModes(fin,sim.res,sim.mdl);
                end
                if (strcmp(string,'%END'))
                    break;
                end
            end
            if (read.status == 1)
                read.checkInput(sim);
            end
        end
    
        %% Functions for reading model TAGS
        %------------------------------------------------------------------
        function analysisModel(read,fin,mdl)
            string = deblank(fgetl(fin));
            if (strcmp(string,'''TRUSS2D'''))
                mdl.anm = Anm_Truss2D();
            elseif (strcmp(string,'''FRAME2D'''))
                mdl.anm = Anm_Frame2D();
            elseif (strcmp(string,'''TRUSS3D'''))
                mdl.anm = Anm_Truss3D();
            elseif (strcmp(string,'''FRAME3D'''))
                mdl.anm = Anm_Frame3D();
            else
                fprintf('Invalid input data: analysis model!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function nodeCoord(read,fin,mdl)
            % Initial check
            if (isempty(mdl.anm))
                fprintf('Analysis model must be provided before node coordinates!\n');
                read.status = 0;
                return;
            end
            
            % Total number of nodes
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of nodes');
            if (read.status == 0)
                return;
            end
            
            % Create vector of Node objects and compute total number of equations
            mdl.nnp = n;
            nodes(n,1) = Node();
            mdl.nodes = nodes;
            mdl.neq = n * mdl.anm.ndof;
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,n,'node ID for coordinate specification');
                if (read.status == 0)
                    return;
                end
                
                % Coordinates (X,Y,Z)
                [coord,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid coordinates of node %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).id = id;
                mdl.nodes(id).coord = coord;
            end
        end
        
        %------------------------------------------------------------------
        function nodeSupport(read,fin,mdl)
            c = Constants();
            
            % Initial check
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before node supports!\n');
                read.status = 0;
                return;
            end
            
            % Total number of nodes with supports
            n = fscanf(fin,'%d',1);
            read.chkInt(n,mdl.nnp,'number of nodes with support conditions');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nnp,'node ID for support condition specification');
                if (read.status == 0)
                    return;
                end
                
                % Support condition flags
                [supp,count] = fscanf(fin,'%d',6);
                if (count~= 6 || ~all(ismember(supp,c.FREE:c.FIXED)))
                    fprintf('Invalid support conditions of node %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).supp = supp;
            end
        end
        
        %------------------------------------------------------------------
        function nodePrescDispl(read,fin,mdl)
            % Initial check
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before prescribed displacements!\n');
                read.status = 0;
                return;
            end
            
            % Total number of nodes with prescribed displacements
            n = fscanf(fin,'%d',1);
            read.chkInt(n,mdl.nnp,'number of nodes with prescribed displacements');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nnp,'node ID for prescribed displacement specification');
                if (read.status == 0)
                    return;
                end
                
                % Prescribed displacements values
                [disp,count] = fscanf(fin,'%f',6);
                if (count ~= 6)
                    fprintf('Invalid prescribed displacements of node %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).prescDispl = disp;
            end
        end
        
        %------------------------------------------------------------------
        function nodeInitCondition(read,fin,mdl)
            % Initial check
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before initial conditions!\n');
                read.status = 0;
                return;
            end
            
            % Total number of nodes with initial conditions
            n = fscanf(fin,'%d',1);
            read.chkInt(n,mdl.nnp,'number of nodes with initial conditions');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nnp,'node ID for initial condition specification');
                if (read.status == 0)
                    return;
                end
                
                % Initial condition values
                [ic,count] = fscanf(fin,'%f',18);
                if (count ~= 18)
                    fprintf('Invalid initial conditions of node %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).U0 = ic(1:6);
                mdl.nodes(id).V0 = ic(7:12);
                mdl.nodes(id).A0 = ic(13:18);
            end
        end
        
        %------------------------------------------------------------------
        function materialProperty(read,fin,mdl)
            % Total number of materials
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of materials');
            if (read.status == 0)
                return;
            end
            
            % Create vector of Material objects
            mdl.nmat = n;
            materials(n,1) = Material();
            mdl.materials  = materials;
            
            for i = 1:n
                % Material ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,n,'material ID');
                if (read.status == 0)
                    return;
                end
                
                % Material properties: E,v,rho,alpha
                [p,count] = fscanf(fin,'%f',4);
                read.chkMatProp(p,count,id);
                if (read.status == 0)
                    return;
                end
                
                % Store data
                mdl.materials(id) = Material(id,p(1),p(2),p(3),p(4));
            end
        end
        
        %------------------------------------------------------------------
        function sectionProperty(read,fin,mdl)
            % Total number of cross-sections
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of cross-sections');
            if (read.status == 0)
                return;
            end
            
            load('areas.mat', 'Solucao_candidata');
            
            % Create vector of Section objects
            mdl.nsec = n;
            sections(n,1) = Section();
            mdl.sections  = sections;
            
            for i = 1:n
                % Cross-section ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,n,'cross-section ID');
                if (read.status == 0)
                    return;
                end
                
                % Cross-section properties: Ax,Ay,Az,Ix,Iy,Iz
                [p,count] = fscanf(fin,'%f',8);
                for ind = 1:4
                    p(1) = Solucao_candidata(1,i);
                    p(2) = p(1);
                    p(3) = p(1);
                    r = sqrt(p(1)/pi);
                    p(4) = 0.5*pi*(r^4);
                    p(5) = 0.25*pi*(r^4);
                    p(6) = p(5);
                    p(7) = 2*r;
                    p(8) = p(7);
                    
                    if (count ~= 8 ||...
                        p(1) <= 0  || p(2) <= 0 || p(3) <= 0 ||...
                        p(4) <= 0  || p(5) <= 0 || p(6) <= 0 ||...
                        p(7) <= 0  || p(8) <= 0)
                        fprintf('Invalid properties for cross-section %d\n',id);
                        fprintf('Cross-section properties must be positive values!\n');
                        read.status = 0;
                    return;
                    end
                end 
                
                % Store data
                mdl.sections(id) = Section(id,p(1),p(2),p(3),...
                                              p(4),p(5),p(6),...
                                              p(7),p(8));
            end
        end
        
        %------------------------------------------------------------------
        function element(read,fin,mdl)
            % Initial check
            if (isempty(mdl.anm)      ||...
                isempty(mdl.nodes)    ||...
                isempty(mdl.materials)||...
                isempty(mdl.sections))
                fprintf('Analysis model, node coordinates, materials, and cross-sections must be provided before elements!\n');
                read.status = 0;
                return;
            end
            
            % Number of elements
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of elements');
            if (read.status == 0)
                return;
            end
            
            % Create vector of Element Objects
            mdl.nel = n;
            elems(n,1) = Element();
            mdl.elems  = elems;
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,n,'element ID');
                if (read.status == 0)
                    return;
                end
                
                % Element properties: type, material, section, nodes, hinges, vz
                [p,count] = fscanf(fin,'%f',10);
                read.chkElemProp(mdl,id,p,count,10);
                if (read.status == 0)
                    return;
                end

                % Store data
                mdl.elems(id) = Element(id,p(1),mdl.anm,...
                                        mdl.materials(p(2)),mdl.sections(p(3)),...
                                       [mdl.nodes(p(4));mdl.nodes(p(5))],...
                                       [p(6);p(7)],...
                                       [p(8);p(9);p(10)]);
            end
        end
        
        %------------------------------------------------------------------
        function nodeLoad(read,fin,mdl)
            % Initial check
            if (isempty(mdl.nodes))
                fprintf('Node coordinates must be provided before nodal loads!\n');
                read.status = 0;
                return;
            end
            
            % Total number of nodes with load
            n = fscanf(fin,'%d',1);
            read.chkInt(n,mdl.nnp,'number of nodes with load');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Node ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nnp,'node ID for load specification');
                if (read.status == 0)
                    return;
                end
                
                % Load values
                [load,count] = fscanf(fin,'%f',6);
                if (count ~= 6)
                    fprintf('Invalid load of node %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.nodes(id).load = load;
            end
        end
        
        %------------------------------------------------------------------
        function elemUnifLoad(read,fin,mdl)
            c = Constants();
            
            % Initial check
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before distributed loads!\n');
                read.status = 0;
                return;
            end
            
            % Total number of elements with uniform distributed load
            n = fscanf(fin,'%d',1);
            read.chkInt(n,mdl.nel,'number of elements with uniform distributed load');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nel,'element ID for uniform distributed load specification');
                if (read.status == 0)
                    return;
                end
                
                % Load direction
                dir = fscanf(fin,'%d',1);
                if (dir ~= c.GLOBAL && dir ~= c.LOCAL)
                    fprintf('Invalid uniform distributed load direction of element %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Load information: Qx,Qy,Qz
                [load,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid uniform distributed load specification of element %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.elems(id).unifdir  = dir;
                mdl.elems(id).unifload = load;
            end
        end
        
        %------------------------------------------------------------------
        function elemLinearLoad(read,fin,mdl)
            c = Constants();
            
            % Initial check
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before distributed loads!\n');
                read.status = 0;
                return;
            end
            
            % Total number of elements with linear distributed load
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of elements with linear distributed load');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nel,'element ID for linear distributed load specification');
                if (read.status == 0)
                    return;
                end
                
                % Load direction
                dir = fscanf(fin,'%d',1);
                if (dir ~= c.GLOBAL && dir ~= c.LOCAL)
                    fprintf('Invalid linear distributed load direction of element %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Load information: Qx1,Qy1,Qz2,Qx2,Qy2,Qz2
                [load,count] = fscanf(fin,'%f',6);
                if (count ~= 6)
                    fprintf('Invalid linear distributed load specification of element %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.elems(id).lindir  = dir;
                mdl.elems(id).linload = load;
            end
        end
        
        %------------------------------------------------------------------
        function elemThermalLoad(read,fin,mdl)
            % Initial check
            if (isempty(mdl.elems))
                fprintf('Elements must be provided before thermal loads!\n');
                read.status = 0;
                return;
            end
            
            % Total number of elements with thermal load
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of elements with thermal load');
            if (read.status == 0)
                return;
            end
            
            for i = 1:n
                % Element ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.nel,'element ID for thermal load specification');
                if (read.status == 0)
                    return;
                end
                
                % Load information: DTX, DTY, DTZ
                [load,count] = fscanf(fin,'%f',3);
                if (count ~= 3)
                    fprintf('Invalid thermal load specification of element %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                mdl.elems(id).thermload = load;
            end
        end
        
        %% Functions for reading analysis TAGS
        %------------------------------------------------------------------
        function analysisType(read,fin,sim)
            string = deblank(fgetl(fin));
            if (strcmp(string,'''LINEAR_ELASTIC'''))
                sim.anl = Anl_Linear();
            elseif (strcmp(string,'''NONLINEAR_GEOM'''))
                sim.anl = Anl_Nonlinear();
            elseif (strcmp(string,'''VIBRATION_FREE'''))
                sim.anl = Anl_Vibration();
            elseif (strcmp(string,'''MODAL_BUCKLING'''))
                sim.anl = Anl_ModalBuck();
            elseif (strcmp(string,'''MODAL_VIBRATION'''))
                sim.anl = Anl_ModalVib();
            else
                fprintf('Invalid input data: analysis type!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function tangStiff(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.NONLINEAR_GEOM)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''CR'''))
                anl.tang_mtx = c.COROTATIONAL;
            elseif (strcmp(string,'''UL'''))
                anl.tang_mtx = c.UPDATED_LAGRANGIAN;
            else
                fprintf('Invalid input data: tangent stiffness matrix!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function massMtx(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.VIBRATION_FREE && anl.type ~= c.MODAL_VIBRATION)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''CONSISTENT'''))
                anl.mass_mtx = c.CONSISTENT;
            elseif (strcmp(string,'''LUMPED'''))
                anl.mass_mtx = c.LUMPED;
            else
                fprintf('Invalid input data: mass matrix!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function damping(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.VIBRATION_FREE)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''UNDAMPED'''))
                anl.damping = c.UNDAMPED;
            elseif (strcmp(string,'''PROPORTIONAL'''))
                anl.damping = c.PROPORTIONAL;
            else
                fprintf('Invalid input data: damping type!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function methodNonlin(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.NONLINEAR_GEOM)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''EULER'''))
                anl.method = c.EULER;
            elseif (strcmp(string,'''LCM'''))
                anl.method = c.LOAD_CONTROL;
            elseif (strcmp(string,'''WCM'''))
                anl.method = c.WORK_CONTROL;
            elseif (strcmp(string,'''ALCM_FNP'''))
                anl.method = c.ARC_LENGTH_FNP;
            elseif (strcmp(string,'''ALCM_UNP'''))
                anl.method = c.ARC_LENGTH_UNP;
            elseif (strcmp(string,'''ALCM_CYL'''))
                anl.method = c.ARC_LENGTH_CYL;
            elseif (strcmp(string,'''ALCM_SPH'''))
                anl.method = c.ARC_LENGTH_SHP;
            elseif (strcmp(string,'''MNCM'''))
                anl.method = c.MINIMUM_NORM;
            elseif (strcmp(string,'''ORCM'''))
                anl.method = c.ORTHOGONAL_RES;
            elseif (strcmp(string,'''GDCM'''))
                anl.method = c.GENERAL_DISPL;
            else
                fprintf('Invalid input data: nonlinear solution method!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function methodVibration(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.VIBRATION_FREE)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''RK4'''))
                anl.method = c.RUNGE_KUTTA4;
            elseif (strcmp(string,'''NEWMARK'''))
                anl.method = c.NEWMARK;
            elseif (strcmp(string,'''WILSONT'''))
                anl.method = c.WILSON_THETA;
            elseif (strcmp(string,'''AM'''))
                anl.method = c.ADAMS_MOULTON3;
            else
                fprintf('Invalid input data: vibration solution method!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function incrementType(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.NONLINEAR_GEOM)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''CONSTANT'''))
                anl.incr_type = c.CONSTANT;
            elseif (strcmp(string,'''ADJUSTED'''))
                anl.incr_type = c.ADJUSTED;
            else
                fprintf('Invalid input data: increment type!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function iterationType(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.NONLINEAR_GEOM)
                return;
            end
            
            string = deblank(fgetl(fin));
            if (strcmp(string,'''STANDARD'''))
                anl.iter_type = c.STANDARD;
            elseif (strcmp(string,'''MODIFIED'''))
                anl.iter_type = c.MODIFIED;
            else
                fprintf('Invalid input data: iteration type!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function analysisParamNonlin(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.NONLINEAR_GEOM)
                return;
            end
            
            [p,count] = fscanf(fin,'%f %f %d %d %d %f',6);
            
            if (count ~= 6 || p(1) <= 0 || ~isnumeric(p(2)) || p(3) <= 0 ||...
                              p(4) <= 0 || p(5) <= 0        || p(6) < 0)
                fprintf('Invalid analysis numerical parameter!\n');
                read.status = 0;
                return;
            end
            
            anl.increment  = p(1);
            anl.max_lratio = p(2);
            anl.max_step   = p(3);
            anl.max_iter   = p(4);
            anl.trg_iter   = p(5);
            anl.tol        = p(6);
        end
        
        %------------------------------------------------------------------
        function analysisParamVib(read,fin,anl)
            c = Constants();
            if (anl.type ~= c.VIBRATION_FREE)
                return;
            end
            
            [p,count] = fscanf(fin,'%f %d %f %f',4);
            
            if (count ~= 4 || p(1) <= 0 || p(2) <= 0 || ~isnumeric(p(3)) || ~isnumeric(p(4)))
                fprintf('Invalid analysis numerical parameter!\n');
                read.status = 0;
                return;
            end
            
            anl.increment = p(1);
            anl.max_step  = p(2);
            anl.rrA       = p(3);
            anl.rrB       = p(4);
        end
        
        %% Functions for reading result TAGS
        %------------------------------------------------------------------
        function resultDomain(read,fin,res)
            c = Constants();
            string = deblank(fgetl(fin));
            if (strcmp(string,'''TIME'''))
                res.domain = c.TIME;
            elseif (strcmp(string,'''FREQUENCY'''))
                res.domain = c.FREQUENCY;
            else
                fprintf('Invalid input data: result domain type!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function resultCurves(read,fin,res,mdl)
            % Total number of curves
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of result curves');
            if (read.status == 0)
                return;
            end
            
            res.ncurves = n;
            
            for i = 1:n
                % Curve name
                name = fscanf(fin,'%s',1);
                
                % Node ID
                node = fscanf(fin,'%d',1);
                read.chkInt(node,mdl.nnp,'node ID for result curve specification');
                if (read.status == 0)
                    return;
                end
                
                % Degree of freedom
                dof = fscanf(fin,'%d',1);
                if (~isnumeric(dof) || ~ismember(dof,mdl.anm.gla))
                    fprintf('Invalid degree of freedom for result curve specification!\n');
                    read.status = 0;
                    return;
                end
                
                % Store data
                res.name{i} = name;
                res.node(i) = node;
                res.dof(i)  = dof;
            end
        end
        
        %------------------------------------------------------------------
        function resultModes(read,fin,res,mdl)
            % Total number of modes
            n = fscanf(fin,'%d',1);
            read.chkInt(n,inf,'number of result modes');
            if (read.status == 0)
                return;
            end
            
            res.nmodes = n;
            
            for i = 1:n
                % Mode ID
                id = fscanf(fin,'%d',1);
                read.chkInt(id,mdl.neq,'result mode ID!');
                if (read.status == 0)
                    return;
                end
                
                % Plot scale
                scale = fscanf(fin,'%f',1);
                if (~isnumeric(scale))
                    fprintf('Invalid scale for result mode %d\n',id);
                    read.status = 0;
                    return;
                end
                
                % Store data
                res.mode(i)  = id;
                res.scale(i) = scale;
            end
        end
        
        %% Functions for checking input data
        %------------------------------------------------------------------
        function checkInput(read,sim)
            if isempty(sim.mdl.anm)
                fprintf('Analysis model type not provided!\n');
                read.status = 0;
            elseif isempty(sim.mdl.nodes)
                fprintf('Nodes not provided!\n');
                read.status = 0;
            elseif isempty(sim.mdl.elems)
                fprintf('Elements not provided!\n');
                read.status = 0;
            elseif isempty(sim.mdl.materials)
                fprintf('Materials not provided!\n');
                read.status = 0;
            elseif isempty(sim.mdl.sections)
                fprintf('Cross-section not provided!\n');
                read.status = 0;
            elseif isempty(sim.anl)
                fprintf('Analysis type not provided!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function chkInt(read,val,max,string)
            if (~isnumeric(val))
                string = ['Invalid',' ',string,'!\n'];
                fprintf(string);
                if (max == inf)
                    fprintf('It must be a positive integer!\n');
                else
                    fprintf('It must be a positive integer less/equal to %d!\n',max);
                end
                read.status = 0;
            end
            if (floor(val) ~= ceil(val) || val <= 0 || val > max)
                string = ['Invalid',' ',string,'!\n'];
                fprintf(string);
                if (max == inf)
                    fprintf('It must be a positive integer!\n');
                else
                    fprintf('It must be a positive integer less/equal to %d!\n',max);
                end
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function chkMatProp(read,p,count,id)
            if (count ~= 4)
                fprintf('Invalid properties for material %d\n',id);
                read.status = 0;
            elseif (p(1) <= 0)
                fprintf('Invalid properties for material %d\n',id);
                fprintf('Elasticity modulus must be positive!\n');
                read.status = 0;
            elseif (p(2) <= -1 || p(2) > 0.5)
                fprintf('Invalid properties for material %d\n',id);
                fprintf('Poisson ratio must be between -1.0 and +0.5!\n');
                read.status = 0;
            elseif (p(3) <= 0)
                fprintf('Invalid properties for material %d\n',id);
                fprintf('Density must be a positive value!\n');
                read.status = 0;
            elseif (p(4) <= 0)
                fprintf('Invalid properties for material %d\n',id);
                fprintf('Thermal expansion coeff. must be positive!\n');
                read.status = 0;
            end
        end
        
        %------------------------------------------------------------------
        function chkElemProp(read,mdl,id,p,count,num)
            c = Constants();
            
            read.chkInt(p(2),mdl.nmat,'material ID for element definition')
            read.chkInt(p(3),mdl.nsec,'cross-section ID for element definition')
            read.chkInt(p(4),mdl.nnp,'node ID for element definition')
            read.chkInt(p(5),mdl.nnp,'node ID for element definition')
            
            if (count ~= num)
                fprintf('Invalid properties of element %d\n',id);
                read.status = 0;
            elseif (p(1) ~= c.EULER_BERNOULLI && p(1) ~= c.TIMOSHENKO)
                fprintf('Invalid type flag of element %d\n',id);
                read.status = 0;
            elseif (p(4) == p(5))
                fprintf('Invalid node IDs of element %d\n',id);
                read.status = 0;
            elseif (p(6) ~= c.CONTINUOUS && p(6) ~= c.HINGED)
                fprintf('Invalid end fixity flag of element %d\n',id);
                read.status = 0;
            elseif (p(7) ~= c.CONTINUOUS && p(7) ~= c.HINGED)
                fprintf('Invalid end fixity flag of element %d\n',id);
                read.status = 0;
            elseif (mdl.anm.type == c.TRUSS2D || mdl.anm.type == c.FRAME2D)
                if (p(8) ~= 0 || p(9) ~= 0)
                    fprintf('Elements direction vector "vz" must be (0,0,1) in 2D models!\n');
                    read.status = 0;
                end
            elseif (mdl.anm.type == c.TRUSS3D || mdl.anm.type == c.FRAME3D)
                dx = mdl.nodes(p(5)).coord(1) - mdl.nodes(p(4)).coord(1);
                dy = mdl.nodes(p(5)).coord(2) - mdl.nodes(p(4)).coord(2);
                dz = mdl.nodes(p(5)).coord(3) - mdl.nodes(p(4)).coord(3);
                x  = [dx;dy;dz];
                vz = [p(8);p(9);p(10)];
                if (norm(cross(x,vz)) == 0)
                    fprintf('Elements direction vector "vz" cannot be aligned with its longitudinal axis!\n');
                    read.status = 0;
                end
            end
        end
    end
end