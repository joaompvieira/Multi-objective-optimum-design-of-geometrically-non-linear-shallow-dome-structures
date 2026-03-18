%% Simulation Class
%
% This class defines a simulation object in the NUMA-TF program.
% A simulation object is the one with the highest hierarchical level.
% Its properties are the three fundamental objects to perform the analysis:
% Model, Analysis, and Result.
%
classdef Simulation < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        mdl = [];   % object of the Model class
        anl = [];   % object of the Anl (Analysis) class
        res = [];   % object of the Result class
        opt = [];   % additional options for running analysis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function sim = Simulation(opt)
            sim.opt = opt;
            sim.runAll();
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Open files and execute each simulation.
        function runAll(sim)
            % Close all Matlab figure windows
            if (~sim.opt.overlay)
                close all;
            end
            
            % Get input file names
            file = 'truss3d_t60_novo.txt';
            path = 'models\nonlinear\';
            if (isequal(file,0))
                return;
            end
            if (~iscell(file))
                file = {file};
            end
            
            % Open input files
            nfiles = length(file);
            fin = zeros(nfiles,1);
            for i=1:nfiles
                namein = strcat(path,string(file(i)));
                fin(i) = fopen(namein,'rt');
                if (fin(i) < 0)
                    fprintf('Error opening input file: %s\n',namein);
                    return;
                end
            end
            
            % Execute simulations
            for i=1:nfiles
%                 fprintf('Start of analysis: %s\n',string(file(i)));
                sim.run(fin(i));
                fclose(fin(i));
%                 fprintf('\n\n');
            end
        end
        
        %------------------------------------------------------------------
        % Run one simulation: read file, pre-process/process/pos-process data,
        % and plot/print results.
        function run(sim,fin)
            % Read input file
            if (Read(sim,fin).status == 0)
                return;
            end
            
            % Pre-process
            sim.mdl.preProcess();
            
            % Analysis process and output results
            if (sim.anl.process(sim))
%                 Plot(sim,fin);
%                 Print(sim,fin);
            end
        end
    end
end