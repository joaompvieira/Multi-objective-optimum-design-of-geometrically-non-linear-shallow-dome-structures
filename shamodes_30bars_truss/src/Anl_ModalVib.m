%% Anl_ModalVib Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with natural vibration analyses
% (vibration modes and natural frequencies).
%
classdef Anl_ModalVib < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        mass_mtx = 0   % flag for type of mass matrix
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl_ModalVib()
            c = Constants();
            anl = anl@Anl(c.MODAL_VIBRATION);
            
            % Default analysis options
            anl.mass_mtx = c.CONSISTENT;
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
            
            % Assemble mass matrix
            M = mdl.gblMassMtx(anl.mass_mtx);
            
            % Extract free d.o.f.'s terms
            Ke = Ke(1:mdl.neqf,1:mdl.neqf);
            M  =  M(1:mdl.neqf,1:mdl.neqf);
            
        %------------------------------------------------------------------
        % Massas não estruturais T25
            
%             for cont = 1:18
%                 M(cont,cont) = M(cont,cont) + 45;
%             end
        %------------------------------------------------------------------
            
            % Solve generalized eigenproblem
            [eigvec,eigval] = eig(Ke,M,'vector');
            
            % Compute natural frequencies
            eigval = sqrt(eigval)/(2*pi);
            
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
                fprintf('Natural Frequencies:\n');
                for i = 1:sim.mdl.neqf
                    fprintf('%d -> %f\n',i,sim.res.lbd(i));
                end
                fprintf('\nVibration modes:\n');
                disp(sim.res.U(1:sim.mdl.neqf,:));
            end
        end
    end
end