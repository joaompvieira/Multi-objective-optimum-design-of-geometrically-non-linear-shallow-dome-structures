%% Anl_Vibration Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anl to deal with vibration analysis.
%
classdef Anl_Vibration < Anl
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        mass_mtx  = 0    % flag for type of mass matrix
        damping   = 0;   % flag for type of damping
        method    = 0;   % flag for vibration solution method
        increment = 0;   % increment of time in each step
        max_step  = 0;   % maximum number of steps
        rrA       = 0;   % Rayleigh coefficient A: [C] = A*[M] + B*[K]
        rrB       = 0;   % Rayleigh coefficient B: [C] = A*[M] + B*[K]
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anl = Anl_Vibration()
            c = Constants();
            anl = anl@Anl(c.VIBRATION_FREE);
            
            % Default analysis options
            anl.mass_mtx  = c.CONSISTENT;
            anl.damping   = c.UNDAMPED;
            anl.method    = c.NEWMARK;
            anl.increment = 0.001;
            anl.max_step  = 5000;
            anl.rrA       = 0.01;
            anl.rrB       = 0.01;
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anl
    methods
        %------------------------------------------------------------------
        % Process model data to compute results.
        function status = process(anl,sim)
            c = Constants();
            status = 1;
            mdl = sim.mdl;
            res = sim.res;
            
            % Check for any free d.o.f
            if (mdl.neqf == 0)
                status = 0;
                fprintf('Model with no free degree-of-freedom!\n');
                return;
            end
            
            % Assemble elastic stiffness matrix and check model stability
            K = mdl.gblElastStiffMtx();
            if (anl.singularMtx(mdl,K))
                status = 0;
                fprintf('Unstable model!\n');
                return;
            end
            
            % Assemble Mass matrix
            M = mdl.gblMassMtx(anl.mass_mtx);
            
            % Assemble Damping matrix
            if (anl.damping == c.UNDAMPED)
                C = zeros(mdl.neq);
            elseif (anl.damping == c.PROPORTIONAL)
                C = anl.rrA * M + anl.rrB * K;
            end
            
            % Extract free-free terms of global matrices
            K = K(1:mdl.neqf,1:mdl.neqf);
            M = M(1:mdl.neqf,1:mdl.neqf);
            C = C(1:mdl.neqf,1:mdl.neqf);
            
            % Assemble global initial conditions matrix
            IC = mdl.gblInitCondMtx();
            
            % Initialize results
            res.steps = anl.max_step;
            res.times = zeros(1,anl.max_step+1);
            res.U     = zeros(mdl.neq,anl.max_step+1);
            res.V     = zeros(mdl.neq,anl.max_step+1);
            res.A     = zeros(mdl.neq,anl.max_step+1);
            
            % Solve transient problem
            fprintf('Solving system of differential equations...\n');
            tic
            if (anl.method == c.RUNGE_KUTTA4)
                step = anl.rungeKutta4(anl,mdl,res,K,M,C,IC);
            elseif (anl.method == c.NEWMARK)
                step = anl.newmark(anl,mdl,res,K,M,C,IC);
            elseif (anl.method == c.WILSON_THETA)
                step = anl.wilsonTheta(anl,mdl,res,K,M,C,IC);
            elseif (anl.method == c.ADAMS_MOULTON3)
                step = anl.adamsMoulton3(anl,mdl,res,K,M,C,IC);
            end
            
            % Clean unused steps
            if (step < anl.max_step)
                res.steps = step;
                res.times = res.times(1:step+1);
                res.U = res.U(:,1:step+1);
                res.V = res.V(:,1:step+1);
                res.A = res.A(:,1:step+1);
            end
            
            % Compute results in frequency domain
            if (res.domain == c.FREQUENCY)
                res.U = fft(res.U')';
                res.V = fft(res.V')';
                res.A = fft(res.A')';  
            end
            
            % Print results on command window
            fprintf('Analysis completed!\n');
            fprintf('Time: %.6fs\n',toc);
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Solve transient problem using 4th order Runge-Kutta method.
        function step = rungeKutta4(anl,mdl,res,K,M,C,IC)
            nf   = mdl.neqf;
            inc  = anl.increment;
            endt = anl.max_step * inc;
            
            % Initialize displacement, velocity, and acceleration matrices
            U = zeros(nf, anl.max_step+1);
            V = zeros(nf, anl.max_step+1);
            A = zeros(nf, anl.max_step+1);
            
            % Set initial conditions
            U(:,1) = IC(:,1);
            V(:,1) = IC(:,2);
            A(:,1) = IC(:,3);
            
            % Loop through all steps
            step = 0;
            res.times(1) = 0;
            for i = inc:inc:endt
                step = step + 1;
                Us   = U(:,step);
                Vs   = V(:,step);
                
                % Compute derivatives
                f1v = Vs;
                f1a = M \ (-K*Us - C*Vs);
                
                x   = Us + 0.5*inc*f1v;
                y   = Vs + 0.5*inc*f1a;
                f2v = y;
                f2a = M \ (-K*x - C*y);
                
                x   = Us + 0.5*inc*f2v;
                y   = Vs + 0.5*inc*f2a;
                f3v = y;
                f3a = M \ (-K*x - C*y);
                
                x   = Us + inc*f3v;
                y   = Vs + inc*f3a;
                f4v = y;
                f4a = M \ (-K*x - C*y);
                
                % Compute displacement, velocity and acceleration on the next step
                u = Us + (1/6)*inc*(f1v + 2*f2v + 2*f3v + f4v);
                v = Vs + (1/6)*inc*(f1a + 2*f2a + 2*f3a + f4a);
                a =      (1/6)*    (f1a + 2*f2a + 2*f3a + f4a);
                
                % Store step results
                res.times(step+1) = i;
                U(:,step+1) = u;
                V(:,step+1) = v;
                A(:,step+1) = a;
            end
            
            % Store results
            res.U(1:nf,:) = U;
            res.V(1:nf,:) = V;
            res.A(1:nf,:) = A;
        end
        
        %------------------------------------------------------------------
        % Solve transient problem using Newmark method.
        function step = newmark(anl,mdl,res,K,M,C,IC)
            nf   = mdl.neqf;
            inc  = anl.increment;
            endt = anl.max_step * inc;
            
            % Notation simplification
            inc_2  = 2/inc;
            inc_4  = 4/inc;
            inc2_4 = 4/inc^2;
            
            % Initialize displacement, velocity, and acceleration matrices
            U = zeros(nf, anl.max_step+1);
            V = zeros(nf, anl.max_step+1);
            A = zeros(nf, anl.max_step+1);
            
            % Set initial conditions
            U(:,1) = IC(:,1);
            V(:,1) = IC(:,2);
            A(:,1) = IC(:,3);
            
            % Auxiliar matrices
            Aux_B = inc_2 * C + inc2_4 * M;
            Aux_A = K + Aux_B;
            Aux_C = C + inc_4 * M;
            
            % Loop through all steps
            step = 0;
            res.times(1) = 0;
            for i = inc:inc:endt
                step = step + 1;
                Us   = U(:,step);
                Vs   = V(:,step);
                As   = A(:,step);
                
                % Compute displacement, velocity and acceleration on the next step
                u = Aux_A \ (Aux_B*Us + Aux_C*Vs + M*As);
                v = (u-Us)*inc_2  - Vs;
                a = (u-Us)*inc2_4 - Vs*inc_4 - As;
                
                % Store step results
                res.times(step+1) = i;
                U(:,step+1) = u;
                V(:,step+1) = v;
                A(:,step+1) = a;
            end
            
            % Store results
            res.U(1:nf,:) = U;
            res.V(1:nf,:) = V;
            res.A(1:nf,:) = A;
        end
        
        %------------------------------------------------------------------
        % Solve transient problem using Wilson-Theta method.
        function step = wilsonTheta(anl,mdl,res,K,M,C,IC)
            nf    = mdl.neqf;
            inc   = anl.increment;
            endt  = anl.max_step * inc;
            theta = 1.4;
            
            % Notation simplification
            inc2   = inc^2;
            inc2_6 = inc2/6;
            t_3    = 3/theta;
            ti     = theta * inc;
            ti_2   = ti/2;
            ti_3   = 3/ti;
            ti21   = theta^2*inc;
            ti22   = theta^2*inc^2;
            ti32   = theta^3*inc^2;
            ti32_6 = 6/ti32;
            ti21_6 = 6/ti21;
            
            % Initialize displacement, velocity, and acceleration matrices
            U = zeros(nf, anl.max_step+1);
            V = zeros(nf, anl.max_step+1);
            A = zeros(nf, anl.max_step+1);
            
            % Set initial conditions
            U(:,1) = IC(:,1);
            V(:,1) = IC(:,2);
            A(:,1) = M \ (-K*U(:,1) - C*V(:,1));
            
            % Pseudo stiffness matrix
            Kk = K + (6/ti22)*M + (3/ti)*C;
            
            % Loop through all steps
            step = 0;
            res.times(1) = 0;
            for i = inc:inc:endt
                step = step + 1;
                
                % Pseudo force vector
                Ff  = 6*(A(:,step)/3      + V(:,step)/ti + U(:,step)/ti22)'   * M +...
                        (A(:,step)*(ti_2) + V(:,step)*2  + U(:,step)*(ti_3))' * C;
                
                % Pseudo displacement vector
                Uu  = Kk \ Ff';
                
                % Compute acceleration, velocity and displacement on the next step
                a = (ti32_6)*(Uu-U(:,step)) - (ti21_6)*V(:,step) + (1-t_3)*A(:,step);
                v =  V(:,step) + 0.5*inc*(a+A(:,step));
                u =  U(:,step) + inc*V(:,step) + (inc2_6)*(a+2*A(:,step));
                
                % Store step results
                res.times(step+1) = i;
                U(:,step+1) = u;
                V(:,step+1) = v;
                A(:,step+1) = a;
            end
            
            % Store results
            res.U(1:nf,:) = U;
            res.V(1:nf,:) = V;
            res.A(1:nf,:) = A;
        end
        
        %------------------------------------------------------------------
        % Solve transient problem using 3rd order Adams-Moulton method.
        function step = adamsMoulton3(anl,mdl,res,K,M,C,IC)
            nf   = mdl.neqf;
            inc  = anl.increment;
            endt = anl.max_step * inc;
            
            % Initialize displacement, velocity, and acceleration matrices
            U = zeros(nf, anl.max_step+1);
            V = zeros(nf, anl.max_step+1);
            A = zeros(nf, anl.max_step+1);
            
            % Set initial conditions
            U(:,1) = IC(:,1);
            V(:,1) = IC(:,2);
            A(:,1) = IC(:,3);
            
            % Auxiliar matrices
            Aux_B = 2/inc * C + 4/inc^2 * M;
            Aux_A = K + Aux_B;
            Aux_C = C + 4/inc * M;
            
            % Inverse or pseudo-inverse of mass matrix
            if (rcond(M) > 10e-15)
                invM = inv(M);
            else
                invM = pinv(M);
            end
            
            % First step -> Newmark method
            u = Aux_A \ (Aux_B * U(:,1) + Aux_C * V(:,1) + M * A(:,1));
            v = (u-U(:,1)) * (2/inc)   - V(:,1);
            a = (u-U(:,1)) * (2/inc)^2 - V(:,1) * (4/inc) - A(:,1);
            
            res.times(1) = 0;
            res.times(2) = inc;
            U(:,2) = u;
            V(:,2) = v;
            A(:,2) = a;
            
            % Second step -> second order Adams-Moulton method
            if (2*inc <= endt)
                % Compute two previous derivatives
                f1v = V(:,1);
                f1a = invM*(-K*U(:,1)-C*V(:,1));
                f1  = [f1v;f1a];
                
                f2v = V(:,2);
                f2a = invM*(-K*U(:,2)-C*V(:,2));
                f2  = [f2v;f2a];
                
                D = [ eye(nf)            -(5/12)*inc*eye(nf);
                     (5/12)*inc*(invM*K)  (5/12)*inc*(invM*C) + eye(nf) ];

                % Compute displacement, velocity and acceleration on the next step
                x = [U(:,2);V(:,2)] + (1/12)*inc*(8*f2-f1);
                y = D \ x;
                
                u = y(1:nf,:);
                v = y(nf+1:end,:);
                a = invM*(-K*u - C*v);
                
                % Store step results
                res.times(3) = 2*inc;
                U(:,3) = u;
                V(:,3) = v;
                A(:,3) = a;
            end
            
            % Compute derivative matrix for entering loop
            D = [ eye(nf)            -(9/24)*inc*eye(nf);
                 (9/24)*inc*(invM*K)  (9/24)*inc*(invM*C) + eye(nf) ];
            
            % Loop through all remaining steps
            step = 2;
            for i = 3*inc:inc:endt
                step = step + 1;
                
                % Compute three previous derivatives
                f1v = V(:,step-2);
                f1a = invM*(-K*U(:,step-2)-C*V(:,step-2));
                f1  = [f1v;f1a];
                
                f2v = V(:,step-1);
                f2a = invM*(-K*U(:,step-1)-C*V(:,step-1));
                f2  = [f2v;f2a];
                
                f3v = V(:,step);
                f3a = invM*(-K*U(:,step)-C*V(:,step));
                f3  = [f3v;f3a];
                
                % Compute displacement and speed on the next step
                x = [U(:,step);V(:,step)] + (1/24)*inc*(19*f3-5*f2+f1);
                y  = D \ x;
                
                u = y(1:nf,:);
                v = y(nf+1:end,:);
                a = invM*(-K*u - C*v);
                
                % Store step results
                res.times(step+1) = i;
                U(:,step+1) = u;
                V(:,step+1) = v;
                A(:,step+1) = a;
            end
            
            % Store results
            res.U(1:nf,:) = U;
            res.V(1:nf,:) = V;
            res.A(1:nf,:) = A;
        end
    end
end