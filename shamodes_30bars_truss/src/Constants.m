%% Constants Class
%
% This class defines an object with constant flag variables as properties
% to be used in the NUMA-TF program.
% An object of this class must be created in every function that uses a
% constant variable.
%
classdef Constants < handle
    %% Public properties (Access Only)
    properties (SetAccess = private, GetAccess = public)
        % Types of model
        TRUSS2D = 1;
        FRAME2D = 2;
        TRUSS3D = 3;
        FRAME3D = 4;
        
        % Degrees of freedom
        DX = 1;
        DY = 2;
        DZ = 3;
        RX = 4;
        RY = 5;
        RZ = 6;
        
        % Types of element
        EULER_BERNOULLI = 1;
        TIMOSHENKO      = 2;
        
        % Support conditions
        FREE  = 0;
        FIXED = 1;
        
        % Element end conditions
        HINGED     = 0;
        CONTINUOUS = 1;
        
        % Load directions
        GLOBAL = 0;
        LOCAL  = 1;
        
        % Types of analysis
        LINEAR_ELASTIC  = 1;
        NONLINEAR_GEOM  = 2;
        VIBRATION_FREE  = 3;
        MODAL_BUCKLING  = 4;
        MODAL_VIBRATION = 5;
        
        % Types of tangent stiffness matrix
        COROTATIONAL       = 1;
        UPDATED_LAGRANGIAN = 2;
        
        % Types of mass matrix
        CONSISTENT = 1;
        LUMPED     = 2;
        
        % Types of damping
        UNDAMPED     = 1;
        PROPORTIONAL = 2;
        
        % Solution methods for nonlinear analysis
        EULER          = 1;
        LOAD_CONTROL   = 2;
        WORK_CONTROL   = 3;
        ARC_LENGTH_FNP = 4;
        ARC_LENGTH_UNP = 5;
        ARC_LENGTH_CYL = 6;
        ARC_LENGTH_SHP = 7;
        MINIMUM_NORM   = 8;
        ORTHOGONAL_RES = 9;
        GENERAL_DISPL  = 10;
        
        % Solution methods for vibration analysis
        RUNGE_KUTTA4   = 1;
        NEWMARK        = 2;
        WILSON_THETA   = 3;
        ADAMS_MOULTON3 = 4;
        
        % Types of increment size
        CONSTANT = 1;
        ADJUSTED = 2;
        
        % Types of iteration strategy
        STANDARD = 1;
        MODIFIED = 2;
        
        % Types of result domain
        TIME      = 1;
        FREQUENCY = 2;
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function c = Constants()
            return;
        end
    end
end