function nlinear = main()
%% NUMA-TF - Numerical Analysis of Trusses and Frames
% This is the main script file of NUMA-TF.
% Run this script to select the input file with model and analysis data for
% performing the simulation.
% Multiple files can be selected, by holding down the Shift or Ctrl key and
% clicking file names, to run the simulations sequentially.
% For more information, see the README and model instructions files in the
% documents folder.

%% Options
% Set options with flags 0 (NO) or 1 (YES):
%  opt.feedback -> Show analysis feedback on command window
%  opt.report   -> Write results in a report file
%  opt.plot     -> Plot results (deformation, equilibrium paths, transient response)
%  opt.overlay  -> Keep graphs from different analyses and overlay results

opt.feedback = 1;
opt.report   = 0;
opt.plot     = 1;
opt.overlay  = 0;

%% Run Analysis
% clc; clearvars -except opt
addpath('src');
Simulation(opt);
etp = ans.res.steps + 1;
nlinear = ans.res.U(:,etp);
end