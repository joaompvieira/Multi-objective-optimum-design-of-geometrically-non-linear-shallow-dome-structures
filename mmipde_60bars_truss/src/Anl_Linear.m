%% Anl_Nonlinear Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with linear-elastic analysis.
%
classdef Anl_Linear < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl_Linear()
            c = Constants();
            anl = anl@Anl(c.LINEAR_ELASTIC);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function status = process(anl,sim)
            status = 1;
            mdl = sim.mdl;
            res = sim.res;
            
            % Initialize vector of nodal loads and add contributions of
            % nodal forces, element loads, and prescribed displacements
            P = zeros(mdl.neq,1);
            P = mdl.addNodalLoad(P);
            P = mdl.addEquivLoad(P);
            P = mdl.addPrescDispl(P);
            
            % Assemble elastic stiffness matrix and check model stability
            Ke = mdl.gblElastStiffMtx();
            if (anl.singularMtx(mdl,Ke))
			    status = 0;
                fprintf('Status: Unstable model!\n');
                return;
            end
            
            % Solve linear-elastic analysis
            [res.U,res.F] = anl.solveSystem(mdl,Ke,P,true,true);
            
            % Compute and store element internal forces from nodal displacements
            % (corner forces)
            mdl.gblIntForceVctLinear(res.U);
            
            % Compute element end forces with stored corner forces
            mdl.elemEndForce();
            
            % Print results on command window
            anl.printResults(sim);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Print results on command window.
        function printResults(sim)
            fprintf('Status: Analysis completed!\n\n');
            if (sim.opt.feedback)
                fprintf('Free displacements:\n');
                for i = 1:sim.mdl.neqf
                    [dof,node] = find(sim.mdl.ID == i);
                    fprintf('Node %d | DOF %d: %.6f\n',node,dof,sim.res.U(i));
                end
                fprintf('\nFixed reactions:\n');
                for i = sim.mdl.neqf+1:sim.mdl.neq
                    [dof,node] = find(sim.mdl.ID == i);
                    fprintf('Node %d | DOF %d: %.6f\n',node,dof,sim.res.F(i));
                end
                fprintf('\nEnd forces:\n');
                for i = 1:sim.mdl.nel
                    fprintf('Element %d: {',sim.mdl.elems(i).id);
                    for j = 1:length(sim.mdl.elems(i).Fi)
                        fprintf(" %10.6f ",sim.mdl.elems(i).Fi(j));
                    end
                    fprintf("}\n");
                end
            end
        end
    end
end