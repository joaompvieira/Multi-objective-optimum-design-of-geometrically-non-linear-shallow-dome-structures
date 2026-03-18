%% Anl_ModalBuck Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with stability analyses
% (buckling modes and critical loads).
%
classdef Anl_ModalBuck < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        tang_mtx = 0;   % flag for type of tangent stiffness matrix
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl_ModalBuck()
            c = Constants();
            anl = anl@Anl(c.MODAL_BUCKLING);
            
            % Default analysis options
            anl.tang_mtx = c.UPDATED_LAGRANGIAN;
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
            tic
            
            % Check for any free d.o.f
            if (mdl.neqf == 0)
                status = 0;
                fprintf('Model with no free degree-of-freedom!\n');
                return;
            end
            
            % Assemble elastic stiffness matrix and check model stability
            Ke = mdl.gblElastStiffMtx();
            if (anl.singularMtx(mdl,Ke))
                status = 0;
                fprintf('Status: Unstable model!\n');
                return;
            end
            
            % Initialize vector of nodal loads and add contributions of
            % nodal forces and distributed loads
            P = zeros(mdl.neq,1);
            P = mdl.addNodalLoad(P);
            P = mdl.addEquivLoad(P);
            
            % Solve linear-elastic analysis
            U = anl.solveSystem(mdl,Ke,P,false,false);
            
            % Compute and store element internal forces from nodal displacements
            mdl.gblIntForceVctLinear(U);
            
            % Assemble geometric stiffness matrix
            Kg = mdl.gblGeomStiffMtx();
            
            % Extract free d.o.f.'s terms and change sign of geometric matrix
            Ke =  Ke(1:mdl.neqf,1:mdl.neqf);
            Kg = -Kg(1:mdl.neqf,1:mdl.neqf);
            
            % Solve generalized eigenproblem
            [eigvec,eigval] = eig(Ke,Kg,'vector');
            
            % Sort order by absolute value and normalize eigenvectors
            [eigvec,eigval] = anl.sortNorm(mdl,eigvec,eigval);
            
            % Store results including all d.o.f.'s
            res.lbd = eigval;
            res.U = zeros(mdl.neq,length(eigval));
            res.U(1:mdl.neqf,:) = eigvec;
            
            % Print results on command window
            anl.printResults(sim,toc);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Print results on command window.
        function printResults(sim,toc)
            fprintf('Status: Analysis completed!\n');
            fprintf('Time: %.6fs\n\n',toc);
            if (sim.opt.feedback)
                fprintf('Critical load ratios:\n');
                for i = 1:sim.mdl.neqf
                    fprintf('%d -> %f\n',i,sim.res.lbd(i));
                end
                fprintf('\nBuckling modes:\n');
                disp(sim.res.U(1:sim.mdl.neqf,:));
            end
        end
    end
end