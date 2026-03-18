%% Model Class
%
% This class defines a model object in the NUMA-TF program.
% A model object is responsible for storing the global variables of a
% structural model.
% The methods of a model object are functions that are not dependent on the 
% analysis model or analysis type.
%
classdef Model < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % General:
        anm = [];          % object of the Anm (Analysis Model) class
        
        % Model properties:
        nnp       = 0;     % number of nodes
        nodes     = [];    % vector of objects of the Node class
        nel       = 0;     % number of elements
        elems     = [];    % vector of objects of the Element class
        nmat      = 0;     % number of materials
        materials = [];    % vector of objects of the Material class
        nsec      = 0;     % number of cross-sections
        sections  = [];    % vector of objects of the Section class
        
        % Degree-of-freedom numbering:
        ID   = [];         % global d.o.f. numbering matrix
        neq  = 0;          % number of d.o.f.'s
        neqf = 0;          % number of free d.o.f.'s
        neqc = 0;          % number of constrained (fixed) d.o.f.'s
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function mdl = Model()
            return;
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Pre-process model data to assemble global ID matrix and local
        % element gather vectors.
        function preProcess(mdl)
            % Initialize global d.o.f. numbering matrix
            mdl.setupDOFNum();
            
            % Assemble global d.o.f. numbering matrix
            mdl.assembleDOFNum();
            
            % Assemble element gather vectors
            mdl.assembleGle();
        end
        
        %------------------------------------------------------------------
        % Initialize global d.o.f numbering matrix with ones and zeros,
        % and count total number of equations of free and fixed  d.o.f.'s.
        %  ID matrix initialization:
        %  if ID(k,n) = 0, d.o.f. k of node n is free.
        %  if ID(k,n) = 1, d.o.f. k of node n is fixed.
        function setupDOFNum(mdl)
            c = Constants();
            
            % Dimension global d.o.f. numbering matrix
            mdl.ID = zeros(mdl.anm.ndof,mdl.nnp);
            
            % Initialize number of fixed d.o.f.
            mdl.neqc = 0;
            
            % Count number of fixed d.o.f. and setup ID matrix
            for i = 1:mdl.nnp
                for j = 1:mdl.anm.ndof
                    % Get d.o.f number
                    dof = mdl.anm.gla(j);
                    
                    % Check for fixed d.o.f
                    if (mdl.nodes(i).supp(dof) == c.FIXED)
                        mdl.neqc = mdl.neqc + 1;
                        mdl.ID(j,i) = 1;
                    end
                end
            end
            
            % Compute total number of free d.o.f.
            mdl.neqf = mdl.neq - mdl.neqc;
        end
        
        %------------------------------------------------------------------
        % Assemble global d.o.f numbering matrix:
        %  ID(k,n) = equation number of d.o.f. k of node n
        %  Free d.o.f.'s have the initial numbering.
        %  countF --> counts free d.o.f.'s. (numbered first)
        %  countC --> counts fixed d.o.f.'s. (numbered later)
        function assembleDOFNum(mdl)
            % Initialize equation numbers for free and fixed d.o.f.'s.
            countF = 0;
            countC = mdl.neqf;
            
            % Check if each d.o.f. is free or fixed to increment equation
            % number and store it in the ID matrix.
            for i = 1:mdl.nnp
                for j = 1:mdl.anm.ndof
                    if (mdl.ID(j,i) == 0)
                        countF = countF + 1;
                        mdl.ID(j,i) = countF;
                    else
                        countC = countC + 1;
                        mdl.ID(j,i) = countC;
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assemble element gather vector (gle) that stores element global
        % d.o.f. numbers.
        function assembleGle(mdl)
            for i = 1:mdl.nel
                % Initialize element gather vector
                mdl.elems(i).gle = zeros(2*mdl.anm.ndof,1);
                k = 0;
                
                % Assemble d.o.f.'s of initial node to element gather vector
                for j = 1:mdl.anm.ndof
                    k = k + 1;
                    mdl.elems(i).gle(k) = mdl.ID(j,mdl.elems(i).nodes(1).id);
                end
                
                % Assemble d.o.f.'s of final node to element gather vector
                for j = 1:mdl.anm.ndof
                    k = k + 1;
                    mdl.elems(i).gle(k) = mdl.ID(j,mdl.elems(i).nodes(2).id);
                end
            end
        end
        
        %------------------------------------------------------------------
        % Set initial model cofiguration (length, angles, and rotation
        % matrix of each element).
        function setInitConfig(mdl)
            for i = 1:mdl.nel
                mdl.anm.setInitGeom(mdl.elems(i));
            end
        end
        
        %------------------------------------------------------------------
        % Update current model cofiguration (length, angles, and rotation
        % matrix of each element).
        function updCurConfig(mdl,d_U,U)
            for i = 1:mdl.nel
                mdl.anm.updCurGeom(mdl.elems(i),d_U,U);
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global initial conditions matrix.
        function IC = gblInitCondMtx(mdl)
            % Initialize initial conditions matrix
            IC = zeros(mdl.neqf,3);
            
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).U0) ||...
                    ~isempty(mdl.nodes(i).V0) ||...
                    ~isempty(mdl.nodes(i).A0))
                    for j = 1:mdl.anm.ndof
                        % Get d.o.f numbers
                        id  = mdl.ID(j,i);
                        dof = mdl.anm.gla(j);
                        
                        % Add nodal initial conditions to global matrix
                        if (id <= mdl.neqf)
                            IC(id,1) = mdl.nodes(i).U0(dof);
                            IC(id,2) = mdl.nodes(i).V0(dof);
                            IC(id,3) = mdl.nodes(i).A0(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add contribution of nodal loads to reference load vector.
        function Pref = addNodalLoad(mdl,Pref)
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).load))
                    for j = 1:mdl.anm.ndof
                        % Get d.o.f numbers
                        id  = mdl.ID(j,i);
                        dof = mdl.anm.gla(j);
                        
                        % Add load to reference load vector
                        Pref(id) = Pref(id) + mdl.nodes(i).load(dof);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add contribution of equivalent nodal loads from element distributed
        % and thermal loads to reference load vector.
        function Pref = addEquivLoad(mdl,Pref)
            for i = 1:mdl.nel
                % Check if element has distributed loads
                if (~isempty(mdl.elems(i).unifload) || ~isempty(mdl.elems(i).linload))
                    % Get element equivalent nodal load vector in global system
                    feg = mdl.anm.equivLoadVctDistrib(mdl.elems(i));
                    
                    % Assemble element vector to global vector
                    gle = mdl.elems(i).gle;
                    Pref(gle) = Pref(gle) + feg;
                end
                
                % Check if element has thermal loads
                if (~isempty(mdl.elems(i).thermload))
                    % Get element equivalent nodal load vector in global system
                    feg = mdl.anm.equivLoadVctThermal(mdl.elems(i));
                    
                    % Assemble element vector to global vector
                    gle = mdl.elems(i).gle;
                    Pref(gle) = Pref(gle) + feg;
                end
            end
        end
        
        %------------------------------------------------------------------
        % Add contribution of prescribed displacements to reference load vector.
        function Pref = addPrescDispl(mdl,Pref)
            % Get vector of prescribed displacement values
            Uc = mdl.gblPrescDisplVct();
            
            if (any(Uc))
                % Get global elastic stiffness matrix
                Ke = mdl.gblElastStiffMtx();
                
                % Get free and fixed d.o.f's numbers
                f = 1:mdl.neqf;
                c = mdl.neqf+1:mdl.neq;
                
                % Add contribution of prescribed displacements
                %   [ Kff Kfc ] * [ Uf ] = [ Pf ] -> Kff*Uf = Pf - Kfc*Uc
                %   [ Kcf Kcc ]   [ Uc ] = [ Pc ]
                Pref(f) = Pref(f) - Ke(f,c) * Uc;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble vector of prescribed displacements values, considering
        % only fixed d.o.f's.
        function Uc = gblPrescDisplVct(mdl)
            % Initialize prescribed displacements vector
            Uc = zeros(mdl.neqc,1);
            
            for i = 1:mdl.nnp
                if (~isempty(mdl.nodes(i).prescDispl))
                    for j = 1:mdl.anm.ndof
                        % Get d.o.f numbers
                        id  = mdl.ID(j,i);
                        dof = mdl.anm.gla(j);
                        
                        % Add prescribed displacement
                        if (id > mdl.neqf && mdl.nodes(i).prescDispl(dof) ~= 0)
                            Uc(id-mdl.neqf) = mdl.nodes(i).prescDispl(dof);
                        end
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global elastic stiffness matrix.
        function Ke = gblElastStiffMtx(mdl)
            % Initialize global elastic matrix
            Ke = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element elastic matrix in global system
                keg = mdl.anm.elastStiffMtx(mdl.elems(i));
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                Ke(gle,gle) = Ke(gle,gle) + keg;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global geometric stiffness matrix.
        function Kg = gblGeomStiffMtx(mdl)
            % Initialize global geometric matrix
            Kg = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element geometric matrix in global system
                kgg = mdl.anm.geomStiffMtx(mdl.elems(i));
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                Kg(gle,gle) = Kg(gle,gle) + kgg;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global tangent stiffness matrix.
        function Kt = gblTangStiffMtx(mdl,tang_mtx)
            % Initialize global tangent matrix
            Kt = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element tangent stiffness matrix in global system
                ktg = mdl.anm.tangStiffMtx(mdl.elems(i),tang_mtx);
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                Kt(gle,gle) = Kt(gle,gle) + ktg;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global mass matrix.
        function M = gblMassMtx(mdl,mass_mtx)
            % Initialize global mass matrix
            M = zeros(mdl.neq,mdl.neq);
            
            for i = 1:mdl.nel
                % Get element mass matrix in global system
                mg = mdl.anm.massMtx(mdl.elems(i),mass_mtx);
                
                % Assemble element matrix to global matrix
                gle = mdl.elems(i).gle;
                M(gle,gle) = M(gle,gle) + mg;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global internal force vector due to nodal displacements
        % considering linear effects (small displacements).
        function Fi = gblIntForceVctLinear(mdl,U)
            % Initialize global internal force vector
            Fi = zeros(mdl.neq,1);
            
            for i = 1:mdl.nel
                % Get element internal force vector in global system
                fig = mdl.anm.intForceVctLinear(mdl.elems(i),U);
                
                % Assemble element vector to global vector
                gle = mdl.elems(i).gle;
                Fi(gle) = Fi(gle) + fig;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble global internal force vector due to nodal displacements
        % considering nonlinear effects (large displacements).
        function Fi = gblIntForceVctNonLin(mdl,U)
            % Initialize global internal force vector
            Fi = zeros(mdl.neq,1);
            
            for i = 1:mdl.nel
                % Get element internal force vector in global system
                fig = mdl.anm.intForceVctNonLin(mdl.elems(i),U);
                
                % Assemble element vector to global vector
                gle = mdl.elems(i).gle;
                Fi(gle) = Fi(gle) + fig;
            end
        end
        
        %------------------------------------------------------------------
        % Compute element end forces by adding the contributions of
        % global analysis (nodal displ.) and local analysis (element loads).
        function elemEndForce(mdl)
            for i = 1:mdl.nel
                % Recover end forces of global analysis (corner forces)
                % assuming that these forces are already computed
                fc = mdl.elems(i).Fc;
                
                % Compute end forces of local analysis (fixed-end-forces)
                % in global system
                feg = zeros(length(fc),1);
                
                if (~isempty(mdl.elems(i).unifload) || ~isempty(mdl.elems(i).linload))
                    feg = -mdl.anm.equivLoadVctDistrib(mdl.elems(i));
                end
                if (~isempty(mdl.elems(i).thermload))
                    feg = feg - mdl.anm.equivLoadVctThermal(mdl.elems(i));
                end
                
                % Transform end forces of local analysis to local system
                fel = mdl.elems(i).rot * feg;
                
                % Add contributions
                mdl.elems(i).Fi = fc + fel;
            end
        end
    end
end