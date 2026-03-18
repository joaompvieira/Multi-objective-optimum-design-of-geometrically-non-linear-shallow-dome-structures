%% Anm (Analysis Model) Class
%
% This is an abstract super-class that generically specifies an analysis 
% model in the NUMA-TF program.
%
% Essentially, this super-class declares abstract methods that define
% the particular behavior of an analysis model. These abstract methods
% are the functions that should be implemented in a derived sub-class
% that deals with specific types of analysis models.
%
classdef Anm < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = 0;    % flag for type of analysis model
        ndof = 0;    % number of degrees-of-freedom per node
        gla  = [];   % gather vector (stores local d.o.f. numbers of a node)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm(type,ndof,gla)
            anm.type = type;
            anm.ndof = ndof;
            anm.gla  = gla;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Set element initial geometry (length, angles, rotation matrix).
        setInitGeom(anm,elem);
        
        %------------------------------------------------------------------
        % Update element geometry (length, angles, rotation matrix) with
        % given displacements (increment and total).
        updCurGeom(anm,elem,d_U,U);
        
        %------------------------------------------------------------------
        % Assemble element equivalent nodal load vector in global system
        % for a generic linear distributed load.
        feg = equivLoadVctDistrib(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble element equivalent nodal load vector in global system
        % for a thermal load.
        feg = equivLoadVctThermal(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble element internal force vector due to nodal displacements
        % considering linear effects (small displacements).
        fig = intForceVctLinear(anm,elem,U);
        
        %------------------------------------------------------------------
        % Assemble element internal force vector due to nodal displacements
        % considering nonlinear effects (large displacements).
        fig = intForceVctNonLin(anm,elem,U);
        
        %------------------------------------------------------------------
        % Assemble element tangent stiffness matrix in global system.
        ktg = tangStiffMtx(anm,elem,tang_mtx);
        
        %------------------------------------------------------------------
        % Assemble element elastic stiffness matrix in global system.
        keg = elastStiffMtx(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble element geometric stiffness matrix in global system.
        kgg = geomStiffMtx(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble element corotational stiffness matrix in global system.
        kcr = corotStiffMtx(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble element mass matrix in global system.
        mg = massMtx(anm,elem,mass_mtx);
        
        %------------------------------------------------------------------
        % Assemble element consistent mass matrix in global system.
        mcg = consistMassMtx(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble element lumped mass matrix in global system.
        mlg = lumpedMassMtx(anm,elem);
        
        %------------------------------------------------------------------
        % Assemble displacement shape function matrix evaluated in a given
        % element position.
        N = shapeFcnMtx(anm,elem,x);
    end
end