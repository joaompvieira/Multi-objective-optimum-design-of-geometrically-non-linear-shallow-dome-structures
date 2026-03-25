%% Result Class
%
% This class defines a result object in the NUMA-TF program.
% A result object is responsible for storing the analysis results.
%
classdef Result < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Result options
        domain  = 0;   % flag for type of domain: time or frequency (for vibration analysis)
        ncurves = 0;   % number of curves to plot
        name    = [];  % vector of curves names
        node    = [];  % vector of curves nodes
        dof     = [];  % vector of curves d.o.f's
        nmodes  = 0;   % number of modes to plot
        mode    = [];  % vector of modes ID numbers
        scale   = [];  % vector of modes plotting scales
        
        % Result storage
        steps = 0;     % number of performed steps
        lbd   = [];    % vector of load ratios/eigenvalues of all steps/modes
        times = [];    % vector of times of all steps
        F     = [];    % vector of nodal forces and reactions
        U     = [];    % matrix of nodal displacement vectors of all steps/modes
        V     = [];    % matrix of nodal velocity vectors of all steps
        A     = [];    % matrix of nodal acceleration vectors of all steps
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function res = Result()
            % Default result options
            c = Constants();
            res.domain = c.TIME;
        end
    end
end