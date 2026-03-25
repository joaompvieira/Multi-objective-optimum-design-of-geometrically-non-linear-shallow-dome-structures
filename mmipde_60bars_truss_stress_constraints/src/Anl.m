%% Anl (Analysis) Class
%
% This is an abstract super-class that generically specifies an analysis 
% type in the NUMA-TF program.
%
% Essentially, this super-class declares abstract methods that are the
% functions that should be implemented in a derived sub-class that deals
% with a specific type of analysis.
%
classdef Anl < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        type = 0;   % flag for type of analysis
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl(type)
            anl.type = type;
        end
    end
    
    %% Abstract methods
    % Declaration of abstract methods implemented in derived sub-classes.
    methods (Abstract)
        %------------------------------------------------------------------
        % Process model data to compute results.
        status = process(anl,sim);
    end
    
    %% Public methods
    % Generic routines used by multiple analysis options.
    methods
        %------------------------------------------------------------------
        % Check matrix singularity by checking its reciprocal condition number.
        % A very low reciprocal condition number indicates that the matrix
        % is badly conditioned and may be singular.
        function singular = singularMtx(~,mdl,K)
            if (rcond(K(1:mdl.neqf,1:mdl.neqf)) < 10e-15)
                singular = 1;
            else
                singular = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Partition and solve a linear system of equations.
        %  f --> free d.o.f. (natural B.C. - unknown) 
        %  c --> constrained d.o.f. (essential B.C. - known) 
        %
        % [ Kff Kfc ] * [ Uf ] = [ Pf ]
        % [ Kcf Kcc ]   [ Uc ] = [ Pc ]
        %
        function [U,P] = solveSystem(~,mdl,K,P,addUc,compPc)
            % Partition system of equations
            Kff = K(1:mdl.neqf,1:mdl.neqf);
            Kcf = K(mdl.neqf+1:mdl.neq,1:mdl.neqf);
            Kcc = K(mdl.neqf+1:mdl.neq,mdl.neqf+1:mdl.neq);
            Pf  = P(1:mdl.neqf);
            Pc  = P(mdl.neqf+1:mdl.neq);
            
            % Solve for free d.o.f. displacements
            Uf = Kff\Pf;
            
            % Reconstruct global displacements vector
            if (addUc)
                Uc = mdl.gblPrescDisplVct();
            else
                Uc = zeros(mdl.neqc,1);
            end
            U = [Uf;Uc];
            
            % Recover forcing unknown values (reactions) at essential B.C.
            % It is assumed that the Pc vector currently stores nodal loads
            % applied directly to fixed d.o.f's.
            % Superimpose computed reaction values to nodal loads, with
            % inversed direction, that were applied directly to fixed d.o.f.'s.
            if (compPc)
                Pc = -Pc + Kcf * Uf + Kcc * Uc;
                P = [Pf;Pc];
            end
        end
        
        %------------------------------------------------------------------
        % Sort order of eigenvalues, and corresponding eigenvectors, by
        % absolute value, and normalize eigenvectors w.r.t. largest value.
        function [v,d] = sortNorm(~,mdl,v,d)
            [~,ind] = sort(abs(d));
            d = d(ind);
            v = v(:,ind);
            
            % Normalize eigenvectors
            for i = 1:mdl.neqf
                % Get maximum absolute value of current eigenvector
                [~,dof] = max(abs(v(:,i)));
                m = v(dof,i);
                
                % normalize eigenvector w.r.t. maximum value
                for j = 1:mdl.neqf
                    v(j,i) = v(j,i)/m;
                end
            end
        end
    end
end