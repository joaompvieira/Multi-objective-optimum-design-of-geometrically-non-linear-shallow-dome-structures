%% Node Class
%
% This class defines a node object in the NUMA-TF program.
% A node is a joint between two or more elements, or any element end,
% used to discretize the model. It is always considered as a three-dimensional
% entity, with six degrees-of-freedom, which may be free or fixed by supports.
% It may have applied loads, prescribed displacements towards fixed
% degrees-of-freedom, and initial conditions towards free degrees-of-freedom.
%
classdef Node < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id         = 0;   % identification number
        coord      = [];  % vector of coordinates in global system [X Y Z]
        supp       = [];  % vector of support flags [DX DY DZ RX RY RZ]
        load       = [];  % vector of applied load components [FX FY FZ MX MY MZ]
        prescDispl = [];  % vector of prescribed displacement values [DX DY DZ RX RY RZ]
        U0         = [];  % vector of initial displacement components [dDX dDY dDZ dRX dRY dRZ]
        V0         = [];  % vector of initial velocities components [vDX vDY vDZ vRX vRY vRZ]
        A0         = [];  % vector of initial acceleration components [aDX aDY aDZ aRX aRY aRZ]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function node = Node()
            node.supp = zeros(6,1);
        end
    end
end