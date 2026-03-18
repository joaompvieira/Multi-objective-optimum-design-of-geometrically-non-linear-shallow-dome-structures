%% Section Class
%
% This class defines a cross-section object in the NUMA-TF program.
% All cross-sections are considered to be of a generic type, which means
% that particular shapes are not specified, only their geometric properties
% are provided, such as area and moment of inertia.
%
classdef Section < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        id = 0;  % identification number
        Ax = 0;  % area relative to local x-axis (full area)
        Ay = 0;  % area relative to local y-axis (effective shear area)
        Az = 0;  % area relative to local z-axis (effective shear area)
        Ix = 0;  % moment of inertia relative to local x-axis (torsion inertia)
        Iy = 0;  % moment of inertia relative to local y-axis (bending inertia)
        Iz = 0;  % moment of inertia relative to local z-axis (bending inertia)
        Hy = 0;  % height relative to local y-axis
        Hz = 0;  % height relative to local z-axis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function sec = Section(id,Ax,Ay,Az,Ix,Iy,Iz,Hy,Hz)
            if (nargin > 0)
                sec.id = id;
                sec.Ax = Ax;
                sec.Ay = Ay;
                sec.Az = Az;
                sec.Ix = Ix;
                sec.Iy = Iy;
                sec.Iz = Iz;
                sec.Hy = Hy;
                sec.Hz = Hz;
            end
        end
    end
end