%% Element Class
%
% This class defines an element object in the NUMA-TF program.
%
% In 2D models, the local axes of an element are defined uniquely in the
% following manner:
% * The local x-axis of an element is its longitudinal axis, from its
%   initial node to its final node.
% * The local z-axis of an element is always in the direction of the
%   global Z-axis, which is perpendicular to the model plane and its
%   positive direction points out of the screen.
% * The local y-axis of an element lays on the global XY-plane and is
%   perpendicular to the element x-axis in such a way that the
%   cross-product (x-axis * y-axis) results in a vector in the
%   global Z direction.
%
% In 3D models, the local y-axis and z-axis are defined by an auxiliary
% vector vz = (vzx,vzy,vzz), which is an element property and should be
% specified as an input data of each element:
% * The local x-axis of an element is its longitudinal axis, from its
%   initial node to its final node.
% * The auxiliary vector lays in the local xz-plane of an element, and the
%   cross-product (vz * x-axis) defines the the local y-axis vector.
% * The direction of the local z-axis is then calculated with the
%   cross-product (x-axis * y-axis).
%
% In 2D models, the auxiliary vector vz must be (0,0,1).
% In 3D models, it is important that the auxiliary vector is not parallel
% to the local x-axis; otherwise, the cross-product (vz * x-axis) is zero.
%
classdef Element < handle
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        id        = 0;   % identification number
        type      = 0;   % flag for type of element (Euler-Bernoulli or Timoshenko)
        anm       = [];  % object of the Anm (Analysis Model) class
        gle       = [];  % gather vector (stores element global d.o.f. numbers)
        material  = [];  % object of the Material class
        section   = [];  % object of the Section class
        nodes     = [];  % vector of objects of the Node class [initial_node, final_node]
        fixend    = [];  % vector of flags of rotation liberation at both ends [initial_end, final_end]
        vz        = [];  % auxiliary direction vector in local xz-plane [vz_X, vz_Y, vz_Z]
        len_ini   = 0;   % initial length
        len_cur   = 0;   % current length
        ang_ini   = 0;   % initial angle with horizontal axis (for 2D models)
        ang_cur   = 0;   % current angle with horizontal axis (for 2D models)
        T_ini     = [];  % initial basis rotation transformation matrix of element (for 3D frame models)
        T         = [];  % current basis rotation transformation matrix of element
        TN1       = [];  % basis rotation transformation matrix of node 1 (for 3D frame models)
        TN2       = [];  % basis rotation transformation matrix of node 2 (for 3D frame models)
        rot       = [];  % element rotation transformation matrix from global system to local system
        Rm        = [];  % element mean rotation matrix for large rotations (for 3D Corotational tangent matrix)
        unifdir   = 0;   % flag for direction of uniform distributed loads
        unifload  = [];  % vector of uniform distributed load components [QX QY QZ]
        lindir    = 0;   % flag for direction of linear distributed loads
        linload   = [];  % vector of linear distributed load components [QX1 QY1 QZ1 QX2 QY2 QZ2]
        thermload = [];  % vector of thermal load components [DTX, DTY, DTZ]
        Fc        = [];  % vector of corner forces (internal forces due to nodal displacements only)
        Fi        = [];  % vector of internal forces (from nodal displacements and element loads)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function elem = Element(id,type,anm,mat,sec,nodes,fixend,vz)
            if (nargin > 0)
                % Set properties
                elem.id       = id;
                elem.type     = type;
                elem.anm      = anm;
                elem.material = mat;
                elem.section  = sec;
                elem.nodes    = nodes;
                elem.fixend   = fixend;
                elem.vz       = vz;
                elem.Fc       = zeros(2*elem.anm.ndof,1);
                elem.Fi       = zeros(2*elem.anm.ndof,1);
                
                % Set initial configuration
                elem.anm.setInitGeom(elem);
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Assemble a vector of distributed load values at end nodes in
        % element local system: [Qx1 Qy1 Qz1 Qx2 Qy2 Qz2].
        function Q = distribLoads(elem)
            c = Constants();
            
            % Initialize vector of load values
            Q = zeros(6,1);
            
            % Add uniform load contribution in local system
            if ~isempty(elem.unifload)
                % Rotate uniform load to local system
                if (elem.unifdir == c.GLOBAL)
                    lcl_unif = elem.T * elem.unifload;
                elseif (elem.unifdir == c.LOCAL)
                    lcl_unif = elem.unifload;
                end
                
                % Add loads in local system
                Q(1) = lcl_unif(1);
                Q(2) = lcl_unif(2);
                Q(3) = lcl_unif(3);
                Q(4) = lcl_unif(1);
                Q(5) = lcl_unif(2);
                Q(6) = lcl_unif(3);
            end  
            
            % Add linear load contribution in local system
            if ~isempty(elem.linload)
                % Rotate linear load to local system
                if (elem.lindir == c.GLOBAL)
                    Trot = blkdiag(elem.T,elem.T);
                    lcl_linear = Trot * elem.linload;
                elseif (elem.lindir == c.LOCAL)
                    lcl_linear = elem.linload;
                end
                
                % Add loads in local system
                Q(1) = Q(1) + lcl_linear(1);
                Q(2) = Q(2) + lcl_linear(2);
                Q(3) = Q(3) + lcl_linear(3);
                Q(4) = Q(4) + lcl_linear(4);
                Q(5) = Q(5) + lcl_linear(5);
                Q(6) = Q(6) + lcl_linear(6);
            end
        end
        
        %------------------------------------------------------------------
        % Static condensation of element matrix to consider end liberations.
        function K = condenseEndLib(elem,K)
            c = Constants();
            ndof = elem.anm.ndof;
            for i = 1:2*ndof
                if ((i >= ndof/2+1 && i <= ndof && elem.fixend(1) == c.HINGED) ||...
                    (i >= 2*ndof-(2*(ndof/3-1)) && elem.fixend(2) == c.HINGED))
                    row = K(i,:);
                    for j = 1:2*ndof
                        if (K(j,i) == 0)
                            continue;
                        end
                        aux = K(j,i)/row(i);
                        K(j,:) = K(j,:) - aux*row;
                    end
                end
            end
        end
    end
end