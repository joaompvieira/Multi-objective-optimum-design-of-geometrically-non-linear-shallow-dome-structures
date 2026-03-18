%% Print Class
%
% This class defines a printing object in the NUMA-TF program.
% A printing object is responsible for providing the analysis results in
% text format in an output report file.
%
classdef Print < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        sim   = [];  % object of the simulation class
        fout  = [];  % file ID of output report
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function print = Print(sim,fin)
            % Set simulation object
            print.sim = sim;
            
            % Print results
            if (sim.opt.report)
                % Open output file
                [pathname,filename,~] = fileparts(fopen(fin));
                nameout = [pathname,'\',filename,'.pos'];
                print.fout = fopen(nameout,'w');
                
                % Execute printing
                print.execute();
                
                % Close output file
                fclose(print.fout);
            end
        end
    end
    
    %% Public methods
    methods
        %% Main function
        %------------------------------------------------------------------
        % Execute printing process according to analysis type.
        function execute(print)
            c = Constants();
            if (print.sim.anl.type == c.LINEAR_ELASTIC)
                print.reportLinear(print.sim.mdl,print.sim.res);
            elseif (print.sim.anl.type == c.NONLINEAR_GEOM)
                print.reportNonlinear(print.sim.mdl,print.sim.anl,print.sim.res);
            elseif (print.sim.anl.type == c.VIBRATION_FREE)
                print.reportVibration(print.sim.mdl,print.sim.res);
            elseif (print.sim.anl.type == c.MODAL_BUCKLING)
                print.reportModalBuck(print.sim.res);
            elseif (print.sim.anl.type == c.MODAL_VIBRATION)
                print.reportModalVib(print.sim.res);
            end
        end
        
        %% Functions for printing results
        %------------------------------------------------------------------
        % Print results of linear-elastic analysis.
        function reportLinear(self,mdl,res)
            c = Constants();
            fprintf(self.fout,"LINEAR-ELASTIC ANALYSIS REPORT\n\n");
            if (mdl.anm.type == c.TRUSS2D)
                fprintf(self.fout,"Analysis model: TRUSS2D\n\n");
            elseif (mdl.anm.type == c.FRAME2D)
                fprintf(self.fout,"Analysis model: FRAME2D\n\n");
            elseif (mdl.anm.type == c.TRUSS3D)
                fprintf(self.fout,"Analysis model: TRUSS3D\n\n");
            elseif (mdl.anm.type == c.FRAME3D)
                fprintf(self.fout,"Analysis model: FRAME3D\n\n");
            end
            fprintf(self.fout,'Free displacements:\n');
            for i = 1:mdl.neqf
                [d,n] = find(mdl.ID == i);
                fprintf(self.fout,'Node %d | DOF %d: %.6f\n',n,d,res.U(i));
            end
            fprintf(self.fout,'\nFixed reactions:\n');
            for i = mdl.neqf+1:mdl.neq
                [d,n] = find(mdl.ID == i);
                fprintf(self.fout,'Node %d | DOF %d: %.6f\n',n,d,res.F(i));
            end
            fprintf(self.fout,'\nEnd forces:\n');
            for i = 1:mdl.nel
                fprintf(self.fout,'Element %d:\n',mdl.elems(i).id);
                if (mdl.anm.type == c.TRUSS2D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(3));
                elseif (mdl.anm.type == c.FRAME2D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fy1: %.6f\n',mdl.elems(i).Fi(2));
                    fprintf(self.fout,'   Mz1: %.6f\n',mdl.elems(i).Fi(3));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(4));
                    fprintf(self.fout,'   Fy2: %.6f\n',mdl.elems(i).Fi(5));
                    fprintf(self.fout,'   Mz2: %.6f\n',mdl.elems(i).Fi(6));
                elseif (mdl.anm.type == c.TRUSS3D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(4));
                elseif (mdl.anm.type == c.FRAME3D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fy1: %.6f\n',mdl.elems(i).Fi(2));
                    fprintf(self.fout,'   Fz1: %.6f\n',mdl.elems(i).Fi(3));
                    fprintf(self.fout,'   Mx1: %.6f\n',mdl.elems(i).Fi(4));
                    fprintf(self.fout,'   My1: %.6f\n',mdl.elems(i).Fi(5));
                    fprintf(self.fout,'   Mz1: %.6f\n',mdl.elems(i).Fi(6));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(7));
                    fprintf(self.fout,'   Fy2: %.6f\n',mdl.elems(i).Fi(8));
                    fprintf(self.fout,'   Fz2: %.6f\n',mdl.elems(i).Fi(9));
                    fprintf(self.fout,'   Mx2: %.6f\n',mdl.elems(i).Fi(10));
                    fprintf(self.fout,'   My2: %.6f\n',mdl.elems(i).Fi(11));
                    fprintf(self.fout,'   Mz2: %.6f\n',mdl.elems(i).Fi(12));
                end
            end
        end
        
        %------------------------------------------------------------------
        % Print iteration results of nonlinear analysis.
        function reportNonlinear(self,mdl,anl,res)
            c = Constants();
            fprintf(self.fout,"NONLINEAR ANALYSIS REPORT\n\n");
            fprintf(self.fout,"Analysis Information:\n");
            fprintf(self.fout,"Analysis model:...%d\n",  mdl.anm.type);
            fprintf(self.fout,"Stiffness matrix:.%d\n",  anl.tang_mtx);
            fprintf(self.fout,"Solution method...%d\n",  anl.method);
            fprintf(self.fout,"Increment type....%d\n",  anl.incr_type);
            fprintf(self.fout,"Iteration type....%d\n",  anl.iter_type);
            fprintf(self.fout,"Increment value:..%f\n",  anl.increment);
            fprintf(self.fout,"Max. load ratio:..%f\n",  anl.max_lratio);
            fprintf(self.fout,"Max. steps:.......%d\n",  anl.max_step);
            fprintf(self.fout,"Desired iter.:....%d\n",  anl.trg_iter);
            fprintf(self.fout,"Max. iter.:.......%d\n",  anl.max_iter);
            fprintf(self.fout,"tolerance:........%f\n\n",anl.tol);
            for i = 1:res.steps
                fprintf(self.fout,"Step %d\n",i);
                fprintf(self.fout,"Load ratio: %.15f\n",res.lbd(i+1));                
                fprintf(self.fout,"Displ.: {");
                for j = 1:mdl.neq
                    fprintf(self.fout, " %15.8f ",res.U(j,i+1));
                end
                fprintf(self.fout,"}\n\n");
            end
            fprintf(self.fout,'\nEnd forces:\n');
            for i = 1:mdl.nel
                fprintf(self.fout,'Element %d:\n',mdl.elems(i).id);
                if (mdl.anm.type == c.TRUSS2D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(3));
                elseif (mdl.anm.type == c.FRAME2D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fy1: %.6f\n',mdl.elems(i).Fi(2));
                    fprintf(self.fout,'   Mz1: %.6f\n',mdl.elems(i).Fi(3));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(4));
                    fprintf(self.fout,'   Fy2: %.6f\n',mdl.elems(i).Fi(5));
                    fprintf(self.fout,'   Mz2: %.6f\n',mdl.elems(i).Fi(6));
                elseif (mdl.anm.type == c.TRUSS3D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(4));
                elseif (mdl.anm.type == c.FRAME3D)
                    fprintf(self.fout,'   Fx1: %.6f\n',mdl.elems(i).Fi(1));
                    fprintf(self.fout,'   Fy1: %.6f\n',mdl.elems(i).Fi(2));
                    fprintf(self.fout,'   Fz1: %.6f\n',mdl.elems(i).Fi(3));
                    fprintf(self.fout,'   Mx1: %.6f\n',mdl.elems(i).Fi(4));
                    fprintf(self.fout,'   My1: %.6f\n',mdl.elems(i).Fi(5));
                    fprintf(self.fout,'   Mz1: %.6f\n',mdl.elems(i).Fi(6));
                    fprintf(self.fout,'   Fx2: %.6f\n',mdl.elems(i).Fi(7));
                    fprintf(self.fout,'   Fy2: %.6f\n',mdl.elems(i).Fi(8));
                    fprintf(self.fout,'   Fz2: %.6f\n',mdl.elems(i).Fi(9));
                    fprintf(self.fout,'   Mx2: %.6f\n',mdl.elems(i).Fi(10));
                    fprintf(self.fout,'   My2: %.6f\n',mdl.elems(i).Fi(11));
                    fprintf(self.fout,'   Mz2: %.6f\n',mdl.elems(i).Fi(12));
                end
            end
        end
        
        %------------------------------------------------------------------
        % Print step results of vibration analysis.
        function reportVibration(self,mdl,res)
            fprintf(self.fout,"FREE VIBRATION ANALYSIS REPORT\n\n");
            fprintf(self.fout,"Initial Conditions:\n");
            fprintf(self.fout,"Displ.: {");
            for j = 1:mdl.neq
                fprintf(self.fout, " %15.8f ",res.U(j,1));
            end
            fprintf(self.fout,"}\n");
            fprintf(self.fout,"Veloc.: {");
            for j = 1:mdl.neq
                fprintf(self.fout," %15.8f ",res.V(j,1));
            end
            fprintf(self.fout,"}\n");
            fprintf(self.fout,"Accel.: {");
            for j = 1:mdl.neq
                fprintf(self.fout," %15.8f ",res.A(j,1));
            end
            fprintf(self.fout,"}\n\n");
            for i = 1:res.steps
                fprintf(self.fout,"---------------------------------------------- STEP %d ----------------------------------------------\n",i);
                fprintf(self.fout,"Time: %.8f\n",res.times(i+1));
                fprintf(self.fout,"Displ.: {");
                for j = 1:mdl.neq
                    fprintf(self.fout, " %15.8f ",res.U(j,i+1));
                end
                fprintf(self.fout,"}\n");
                fprintf(self.fout,"Veloc.: {");
                for j = 1:mdl.neq
                    fprintf(self.fout," %15.8f ",res.V(j,i+1));
                end
                fprintf(self.fout,"}\n");
                fprintf(self.fout,"Accel.: {");
                for j = 1:mdl.neq
                    fprintf(self.fout," %15.8f ",res.A(j,i+1));
                end
                fprintf(self.fout,"}\n\n");
            end
        end
        
        %------------------------------------------------------------------
        % Print eigenproblem results of buckling modes analysis.
        function reportModalBuck(self,res)
            fprintf(self.fout,"STABILITY ANALYSIS REPORT\n\n");
            fprintf(self.fout,"Critical Load Ratios:\n");
            for i = 1:length(res.lbd)
                fprintf(self.fout,"%d -> %.15f\n", i, res.lbd(i));
            end
            fprintf(self.fout,"\n------------------------------------------------------\n");
            fprintf(self.fout,"Buckling Modes:\n\n");
            for i = 1:length(res.lbd)
                fprintf(self.fout,"MODE %d\n", i);
                fprintf(self.fout,"Critical load ratio: %.15f\n", res.lbd(i));
                fprintf(self.fout,"Normalized buckling mode:\n");
                for j = 1:length(res.U(:,i))
                    fprintf(self.fout,"%d%25.15f\n", j, res.U(j,i));
                end
                fprintf(self.fout,"\n");
            end
        end
        
        %------------------------------------------------------------------
        % Print eigenproblem results of vibration modes analysis.
        function reportModalVib(self,res)
            fprintf(self.fout,"NATURAL VIBRATION ANALYSIS REPORT\n\n");
            fprintf(self.fout,"Natural Frequencies:\n");
            for i = 1:length(res.lbd)
                fprintf(self.fout,"%d -> %.15f\n", i, res.lbd(i));
            end
            fprintf(self.fout,"\n------------------------------------------------------\n");
            fprintf(self.fout,"Vibration Modes:\n\n");
            for i = 1:length(res.lbd)
                fprintf(self.fout,"MODE %d\n", i);
                fprintf(self.fout,"Natural frequency: %.15f\n", res.lbd(i));
                fprintf(self.fout,"Normalized vibration mode:\n");
                for j = 1:length(res.U(:,i))
                    fprintf(self.fout,"%d%25.15f\n", j, res.U(j,i));
                end
                fprintf(self.fout,"\n");
            end
        end
    end
end