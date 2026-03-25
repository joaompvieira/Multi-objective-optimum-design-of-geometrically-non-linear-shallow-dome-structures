%% Material Class
%
% This class defines a material object in the NUMA-TF program.
% All materials are considered to have a linear-elastic bahavior.
% In adition, homogeneous and isotropic properties are also considered,
% that is, all materials have the same properties at every point and in all
% directions.
%
classdef Material < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id    = 0;  % identification number
        E     = 0;  % elasticity modulus
        v     = 0;  % poisson ratio
        G     = 0;  % shear modulus
        rho   = 0;  % density (specific weight divided by gravity acceleration)
        alpha = 0;  % thermal expansion coefficient
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function mat = Material(id,E,v,rho,alpha)
            if (nargin > 0)
                mat.id    = id;
                mat.E     = E;
                mat.v     = v;
                mat.G     = E/(2*(1+v));
                mat.rho   = rho;
                mat.alpha = alpha;
            end
        end
    end
end