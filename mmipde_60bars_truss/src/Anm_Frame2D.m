%% Anm_Frame2D Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anm to deal with 2D frame models.
%
classdef Anm_Frame2D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Frame2D()
            c = Constants();
            anm = anm@Anm(c.FRAME2D,3,[1,2,6]);
        end
    end
    
    %% Public methods
    % Implementation of the abstract methods declared in super-class Anm
    methods
        %------------------------------------------------------------------
        % Set element initial geometry (length, angles, rotation matrix).
        function setInitGeom(~,elem)
            % Get nodal coordinates
            x1 = elem.nodes(1).coord(1);
            y1 = elem.nodes(1).coord(2);
            x2 = elem.nodes(2).coord(1);
            y2 = elem.nodes(2).coord(2);
            
            % Calculate element length
            dx = x2 - x1;
            dy = y2 - y1;
            L  = sqrt(dx^2 + dy^2);
            
            % Calculate element angle with global X-axis
            ang = atan2(dy,dx);
            
            % Calculate unity element local axes
            x = [dx,dy,0];
            x = x/norm(x);
            y = cross(elem.vz,x);
            y = y/norm(y);
            z = cross(x,y);
            
            % Assemble basis transformation matrix
            T = [ x(1) x(2) x(3);
                  y(1) y(2) y(3);
                  z(1) z(2) z(3) ];
              
            % Set element properties
            elem.len_ini = L;
            elem.len_cur = L;
            elem.ang_ini = ang;
            elem.ang_cur = ang;
            elem.T       = T;
            elem.rot     = blkdiag(T,T);
        end
        
        %------------------------------------------------------------------
        % Update element geometry (length, angles, rotation matrix) with
        % given displacements (increment and total).
        function updCurGeom(~,elem,d_U,U)
            % Get total element displacements in global system
            Ug = U(elem.gle);
            
            % Evaluate new nodal coordinates
            x1 = elem.nodes(1).coord(1) + Ug(1);
            y1 = elem.nodes(1).coord(2) + Ug(2);
            x2 = elem.nodes(2).coord(1) + Ug(4);
            y2 = elem.nodes(2).coord(2) + Ug(5);
            
            % Calculate element length
            dx = x2 - x1;
            dy = y2 - y1;
            L  = norm([dx,dy]);

            % Get increment of element displacements in local system
            d_Ul = elem.rot * d_U(elem.gle);
            
            % Compute increment of angle with global X-axis
            dx    = d_Ul(4) - d_Ul(1);
            dy    = d_Ul(5) - d_Ul(2);
            d_ang = atan2(dy,(elem.len_cur+dx));
            
            % Assemble element basis rotation transformation matrix
            c = cos(elem.ang_cur+d_ang);
            s = sin(elem.ang_cur+d_ang);
            
            T = [ c  s  0 ; 
                 -s  c  0 ;
                  0  0  1 ];
            
            % Update element properties
            elem.len_cur = L;
            elem.ang_cur = elem.ang_cur + d_ang;
            elem.T       = T;
            elem.rot     = blkdiag(elem.T,elem.T);
        end
        
        %------------------------------------------------------------------
        % Assemble element equivalent nodal load vector in global system
        % for a generic linear distributed load.
        function feg = equivLoadVctDistrib(anm,elem)
            c = Constants();
            
            % Get distributed load values at end nodes in element local system
            Q = elem.distribLoads();
            
            % Separate uniform load portion from linear load portion
            qx0 = Q(1);
            qy0 = Q(2);
            qx1 = Q(4) - Q(1);
            qy1 = Q(5) - Q(2);
            
            % Compute element fixed-end-forces vector in local system
            if (elem.type == c.EULER_BERNOULLI)
                fel = anm.fixEndForceVctDistrib_EB(elem,qx0,qy0,qx1,qy1);
            elseif (elem.type == c.TIMOSHENKO)
                fel = anm.fixEndForceVctDistrib_Tim(elem,qx0,qy0,qx1,qy1);
            end
            
            % Transform fixed-end-forces in local system to equivalent
            % nodal loads in global system
            feg = -elem.rot' * fel;
        end
        
        %------------------------------------------------------------------
        % Assemble element equivalent nodal load vector in global system
        % for a thermal load.
        function feg = equivLoadVctThermal(anm,elem)
            c = Constants();
            
            % Compute element fixed-end-forces vector in local system
            if (elem.type == c.EULER_BERNOULLI)
                fel = anm.fixEndForceVctThermal_EB(elem);
            elseif (elem.type == c.TIMOSHENKO)
                fel = anm.fixEndForceVctThermal_Tim(elem);
            end
            
            % Transform fixed-end-forces in local system to equivalent
            % nodal loads in global system
            feg = -elem.rot' * fel;
        end
        
        %------------------------------------------------------------------
        % Assemble element internal force vector due to nodal displacements
        % considering linear effects (small displacements).
        function fig = intForceVctLinear(anm,elem,U)
            % Element elastic stiffness matrix in global system
            keg = anm.elastStiffMtx(elem);
            
            % Element nodal displacements in global system
            ug = U(elem.gle);
            
            % Internal forces in global system
            fig = keg*ug;
            
            % Update element corner forces
            elem.Fc = elem.rot * fig;
        end
        
        %------------------------------------------------------------------
        % Assemble element internal force vector due to nodal displacements
        % considering nonlinear effects (large displacements).
        function fig = intForceVctNonLin(anm,elem,U)
            c = Constants();
            
            % Element properties
            E    = elem.material.E;
            A    = elem.section.Ax;
            L0   = elem.len_ini;
            L1   = elem.len_cur;
            ang0 = elem.ang_ini;
            ang1 = elem.ang_cur;
            
            % Length increment since beginning of analysis
            dL = L1 - L0;
            
            % Rigid body rotation since beginning of analysis
            rbr = ang1 - ang0;
            
            % Relative rotations
            r1 = U(elem.gle(3)) - rbr;
            r2 = U(elem.gle(6)) - rbr;
            
            % Internal forces in local system
            N1 = -E*A*dL/L0;
            N2 =  E*A*dL/L0;
            
            if (elem.type == c.EULER_BERNOULLI)
                [M1,M2] = anm.endMoment_EB(elem,r1,r2);
            elseif (elem.type == c.TIMOSHENKO)
                [M1,M2] = anm.endMoment_Tim(elem,r1,r2);
            end
            
            Q1 = (M1+M2)/L0;
            Q2 = -Q1;
            
            % Assemble element internal force vector in local system
            fil = [N1; Q1; M1; N2; Q2; M2];
            
            % Update element corner forces
            elem.Fc = fil;
            
            % Transform element internal force from local to global system
            fig = elem.rot' * fil;
        end
        
        %------------------------------------------------------------------
        % Assemble element tangent stiffness matrix in global system.
        function ktg = tangStiffMtx(anm,elem,tang_mtx)
            c = Constants();
            if (tang_mtx == c.COROTATIONAL)
                ktg = anm.corotStiffMtx(elem);
            elseif (tang_mtx == c.UPDATED_LAGRANGIAN)
                keg = anm.elastStiffMtx(elem);
                kgg = anm.geomStiffMtx(elem);
                ktg = keg + kgg;
            end
        end
        
        %------------------------------------------------------------------
        % Assemble element elastic stiffness matrix in global system.
        function keg = elastStiffMtx(anm,elem)
            c = Constants();
            
            % Compute element elastic matrix in local system
            if (elem.type == c.EULER_BERNOULLI)
                kel = anm.elastStiffMtx_EB(elem);
            elseif (elem.type == c.TIMOSHENKO)
                kel = anm.elastStiffMtx_Tim(elem);
            end
            
            % Transform element matrix from local to global system
            keg = elem.rot' * kel * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble element geometric stiffness matrix in global system.
        function kgg = geomStiffMtx(anm,elem)
            c = Constants();
            
            % Compute element geometric matrix in local system
            if (elem.type == c.EULER_BERNOULLI)
                kgl = anm.geomStiffMtx_EB(elem);
            elseif (elem.type == c.TIMOSHENKO)
                kgl = anm.geomStiffMtx_Tim(elem);
            end
            
            % Transform element matrix from local to global system
            kgg = elem.rot' * kgl * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble element corotational stiffness matrix in global system.
        function kcr = corotStiffMtx(anm,elem)
            c = Constants();
            if (elem.type == c.EULER_BERNOULLI)
                kcr = anm.corotStiffMtx_EB(elem);
            elseif (elem.type == c.TIMOSHENKO)
                kcr = anm.corotStiffMtx_Tim(elem);
            end
        end
        
        %------------------------------------------------------------------
        % Assemble element mass matrix in global system.
        function mg = massMtx(anm,elem,mass_mtx)
            c = Constants();
            if (mass_mtx == c.CONSISTENT)
                mg = anm.consistMassMtx(elem);
            elseif (mass_mtx == c.LUMPED)
                mg = anm.lumpedMassMtx(elem);
            end
        end
        
        %------------------------------------------------------------------
        % Assemble element consistent mass matrix in global system.
        function mcg = consistMassMtx(anm,elem)
            c = Constants();
            
            % Compute element consistent mass matrix in local system
            if (elem.type == c.EULER_BERNOULLI)
                mcl = anm.consistMassMtx_EB(elem);
            elseif (elem.type == c.TIMOSHENKO)
                mcl = anm.consistMassMtx_Tim(elem);
            end
            
            % Transform element matrix from local to global system
            mcg = elem.rot' * mcl * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble element lumped mass matrix in global system.
        function mlg = lumpedMassMtx(~,elem)
            % Element properties
            rho = elem.material.rho;
            A   = elem.section.Ax;
            L   = elem.len_cur;
            M   = rho * A * L;
            
            % Assemble element lumped mass matrix in local system
            mll = M/2 *...
                  [ 1  0  0  0  0  0;
                    0  1  0  0  0  0;
                    0  0  0  0  0  0;
                    0  0  0  1  0  0;
                    0  0  0  0  1  0;
                    0  0  0  0  0  0 ];
            
            % Transform element matrix from local to global system
            mlg = elem.rot' * mll * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble displacement shape function matrix evaluated in a given
        % element position.
        function N = shapeFcnMtx(anm,elem,x)
            c = Constants();
            if (elem.type == c.EULER_BERNOULLI)
                N = anm.shapeFcnMtx_EB(elem,x);
            elseif (elem.type == c.TIMOSHENKO)
                N = anm.shapeFcnMtx_Tim(elem,x);
            end
        end
    end
    
    %% Static methods
    methods (Static)
        %------------------------------------------------------------------
        % Compute element fixed-end-forces vector in local system for a
        % distributed load considering Euler-Bernoulli beam theory.
        function fel = fixEndForceVctDistrib_EB(elem,qx0,qy0,qx1,qy1)
            c = Constants();
            
            % Element length
            L  = elem.len_cur;
            L2 = L^2;
            
            % Compute axial fixed-end-forces
            N1 = -(qx0*L/2 + qx1*L/6);
            N2 = -(qx0*L/2 + qx1*L/3);
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Q1 = -(qy0*L/2   + qy1*3*L/20);
                M1 = -(qy0*L2/12 + qy1*L2/30);
                Q2 = -(qy0*L/2   + qy1*7*L/20);
                M2 =   qy0*L2/12 + qy1*L2/20;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Q1 = -(qy0*3*L/8 + qy1*L/10);
                M1 =   0;
                Q2 = -(qy0*5*L/8 + qy1*2*L/5);
                M2 =   qy0*L2/8  + qy1*L2/15;
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Q1 = -(qy0*5*L/8 + qy1*9*L/40);
                M1 = -(qy0*L2/8  + qy1*7*L2/120);
                Q2 = -(qy0*3*L/8 + qy1*11*L/40);
                M2 =   0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Q1 = -(qy0*L/2 + qy1*L/6);
                M1 =   0;
                Q2 = -(qy0*L/2 + qy1*L/3);
                M2 =   0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Q1; M1; N2; Q2; M2];
        end
        
        %------------------------------------------------------------------
        % Compute element fixed-end-forces vector in local system for a
        % distributed load considering Timoshenko beam theory.
        function fel = fixEndForceVctDistrib_Tim(elem,qx0,qy0,qx1,qy1)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ay = elem.section.Ay;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            L2 = L^2;
            
            % Timoshenko parameters
            om       = E * Iz / (G * Ay * L2);
            mu       = 1 + 12 * om;
            lambda   = 1 + 3  * om;
            zeta     = 1 + 40 * om/3;
            xi       = 1 + 5  * om;
            eta      = 1 + 15 * om;
            vartheta = 1 + 4  * om;
            psi      = 1 + 12 * om/5;
            varpi    = 1 + 20 * om/9;
            epsilon  = 1 + 80 * om/7;
            varrho   = 1 + 10 * om;
            upsilon  = 1 + 5  * om/2;
            varsigma = 1 + 40 * om/11;
            
            % Compute axial fixed-end-forces
            N1 = -(qx0*L/2 + qx1*L/6);
            N2 = -(qx0*L/2 + qx1*L/3);
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Q1 = -(qy0*L/2   + qy1*3*L*zeta/(20*mu));
                M1 = -(qy0*L2/12 + qy1*L2*eta/(30*mu));
                Q2 = -(qy0*L/2   + qy1*7*L*epsilon/(20*mu));
                M2 =   qy0*L2/12 + qy1*L2*varrho/(20*mu);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Q1 = -(qy0*3*L*vartheta/(8*lambda) + qy1*L*xi/(10*lambda));
                M1 =   0;
                Q2 = -(qy0*5*L*psi/(8*lambda)      + qy1*2*L*upsilon/(5*lambda));
                M2 =   qy0*L2/(8*lambda)           + qy1*L2/(15*lambda);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Q1 = -(qy0*5*L*psi/(8*lambda)      + qy1*9*L*varpi/(40*lambda));
                M1 = -(qy0*L2/(8*lambda)           + qy1*7*L2/(120*lambda));
                Q2 = -(qy0*3*L*vartheta/(8*lambda) + qy1*11*L*varsigma/(40*lambda));
                M2 =   0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Q1 = -(qy0*L/2 + qy1*L/6);
                M1 =   0;
                Q2 = -(qy0*L/2 + qy1*L/3);
                M2 =   0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Q1; M1; N2; Q2; M2];
        end
        
        %------------------------------------------------------------------
        % Compute element fixed-end-forces vector in local system for a
        % thermal load considering Euler-Bernoulli beam theory.
        function fel = fixEndForceVctThermal_EB(elem)
            c = Constants();
            
            % Element properties
            E     = elem.material.E;
            alpha = elem.material.alpha;
            A     = elem.section.Ax;
            I     = elem.section.Iz;
            H     = elem.section.Hy;
            L     = elem.len_cur;
            EI    = E*I;
            
            % Temperature variations
            dtx = elem.thermload(1);
            dty = elem.thermload(2);
            
            % Unitary dimensionless temperature gradient
            tg = -(alpha * dty) / H;
            
            % Compute axial fixed-end-forces
            N1 =  E * A * alpha * dtx;
            N2 = -E * A * alpha * dtx;
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Q1 =  0;
                M1 =  tg * EI;
                Q2 =  0;
                M2 = -tg * EI;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Q1 = -tg * 3*EI/(2*L);
                M1 =  0;
                Q2 =  tg * 3*EI/(2*L);
                M2 = -tg * 3*EI/2;
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Q1 =  tg * 3*EI/(2*L);
                M1 =  tg * 3*EI/2;
                Q2 = -tg * 3*EI/(2*L);
                M2 =  0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Q1 = 0;
                M1 = 0;
                Q2 = 0;
                M2 = 0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Q1; M1; N2; Q2; M2];
        end
        
        %------------------------------------------------------------------
        % Compute element fixed-end-forces vector in local system for a
        % thermal load considering Timoshenko beam theory.
        function fel = fixEndForceVctThermal_Tim(elem)
            c = Constants();
            
            % Element properties
            E     = elem.material.E;
            G     = elem.material.G;
            alpha = elem.material.alpha;
            Ax    = elem.section.Ax;
            Ay    = elem.section.Ay;
            I     = elem.section.Iz;
            H     = elem.section.Hy;
            L     = elem.len_cur;
            L2    = L*L;
            EI    = E*I;
            GAy   = G*Ay;
            
            % Timoshenko parameters
            om = EI / (GAy * L2);
            la = 1 + 3  * om;
            
            % Temperature variations
            dtx = elem.thermload(1);
            dty = elem.thermload(2);
            
            % Unitary dimensionless temperature gradient
            tg = -(alpha * dty) / H;
            
            % Compute axial fixed-end-forces
            N1 =  E * Ax * alpha * dtx;
            N2 = -E * Ax * alpha * dtx;
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Q1 =  0;
                M1 =  tg * EI;
                Q2 =  0;
                M2 = -tg * EI;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Q1 = -tg * 3*EI/(2*L*la);
                M1 =  0;
                Q2 =  tg * (3*EI/(2*L*la));
                M2 = -tg * (3*EI/(2*la));
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Q1 =  tg * 3*EI/(2*L*la);
                M1 =  tg * 3*EI/(2*la);
                Q2 = -tg * 3*EI/(2*L*la);
                M2 =  0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Q1 = 0;
                M1 = 0;
                Q2 = 0;
                M2 = 0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Q1; M1; N2; Q2; M2];
        end
        
        %------------------------------------------------------------------
        % Compute element end moments considering Euler-Bernoulli beam theory.
        function [M1,M2] = endMoment_EB(elem,r1,r2)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            Iz = elem.section.Iz;
            L  = elem.len_ini;
            
            % Compute end moments
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS)
                M1 = 2*E*Iz*(2*r1+r2)/L;
                M2 = 2*E*Iz*(r1+2*r2)/L;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                M1 = 0;
                M2 = 3*E*Iz*r2/L;
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                M1 = 3*E*Iz*r1/L;
                M2 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                M1 = 0;
                M2 = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Compute element end moments considering Timoshenko beam theory.
        function [M1,M2] = endMoment_Tim(elem,r1,r2)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ay = elem.section.Ay;
            Iz = elem.section.Iz;
            L  = elem.len_ini;
            
            % Timoshenko parameters
            om = E * Iz / (G * Ay * L * L);
            ni = 1 + 12 * om;
            la = 1 + 3  * om;
            ga = 1 - 6  * om;
            
            % Compute end moments
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS)
                M1 = 2*E*Iz*(2*la*r1+ga*r2)/(ni*L);
                M2 = 2*E*Iz*(ga*r1+2*la*r2)/(ni*L);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                M1 = 0;
                M2 = 3*E*Iz*la*r2/(ni*L);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                M1 = 3*E*Iz*la*r1/(ni*L);
                M2 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                M1 = 0;
                M2 = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Compute element elastic stiffness matrix in local system
        % considering Euler-Bernoulli beam theory.
        function kel = elastStiffMtx_EB(elem)
            % Element properties
            E  = elem.material.E;
            Ax = elem.section.Ax;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            
            % Notation simplifications
            L2 = L^2;
            L3 = L2 * L;
            EA = E * Ax;
            EI = E * Iz;
            
            % Assemble element elastic stiffness matrix in local system
            kel = [ EA/L     0          0        -EA/L    0          0;
                    0        12*EI/L3   6*EI/L2   0      -12*EI/L3   6*EI/L2;
                    0        6*EI/L2    4*EI/L    0      -6*EI/L2    2*EI/L;
                   -EA/L     0          0         EA/L    0          0;
                    0       -12*EI/L3  -6*EI/L2   0       12*EI/L3  -6*EI/L2;
                    0        6*EI/L2    2*EI/L    0      -6*EI/L2    4*EI/L ];
            
            % Static condensation to consider beam end liberations
            kel = elem.condenseEndLib(kel);
        end
        
        %------------------------------------------------------------------
        % Compute element elastic stiffness matrix in local system
        % considering Timoshenko beam theory.
        function kel = elastStiffMtx_Tim(elem)
            % Element properties
            E   = elem.material.E;
            G   = elem.material.G;
            Ax  = elem.section.Ax;
            Ay  = elem.section.Ay;
            Iz  = elem.section.Iz;
            L   = elem.len_cur;
            
            % Timoshenko parameters
            om = E * Iz / (G * Ay * L * L);
            mu = 1 + 12 * om;
            la = 1 + 3  * om;
            ga = 1 - 6  * om;
            
            % Notation simplifications
            EA   = E    * Ax;
            EI   = E    * Iz;
            muL  = mu   * L;
            muL2 = muL  * L;
            muL3 = muL2 * L;
            
            % Assemble element elastic stiffness matrix in local system
            kel = [ EA/L     0            0           -EA/L    0            0;
                    0        12*EI/muL3   6*EI/muL2    0      -12*EI/muL3   6*EI/muL2;
                    0        6*EI/muL2    4*la*EI/muL  0      -6*EI/muL2    2*ga*EI/muL;
                   -EA/L     0            0            EA/L    0            0;
                    0       -12*EI/muL3  -6*EI/muL2    0       12*EI/muL3  -6*EI/muL2;
                    0        6*EI/muL2    2*ga*EI/muL  0      -6*EI/muL2    4*la*EI/muL ];
            
            % Static condensation to consider beam end liberations
            kel = elem.condenseEndLib(kel);
        end
        
        %------------------------------------------------------------------
        % Compute element geometric stiffness matrix in local system
        % considering Euler-Bernoulli beam theory.
        function kgl = geomStiffMtx_EB(elem)
            % Element properties
            N  = elem.Fc(4);
            L  = elem.len_cur;
            L2 = L*L;
            
            % Assemble element geometric stiffness matrix in local system
            kgl = N/L *...
                  [ 1        0        0       -1        0        0;
                    0        6/5      L/10     0       -6/5      L/10;
                    0        L/10     2*L2/15  0       -L/10    -L2/30;
                   -1        0        0        1        0        0;
                    0       -6/5     -L/10     0        6/5     -L/10;
                    0        L/10    -L2/30    0       -L/10    2*L2/15 ];
            
            % Static condensation to consider beam end liberations
            kgl = elem.condenseEndLib(kgl);
        end
        
        %------------------------------------------------------------------
        % Compute element geometric stiffness matrix in local system
        % considering Timoshenko beam theory.
        function kgl = geomStiffMtx_Tim(elem)
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ay = elem.section.Ay;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            N  = elem.Fc(4);
            
            % Timoshenko parameters
            om = E * Iz / (G * Ay * L * L);
            mu = 1 + 12 * om;
            a  = 1 + 20 * om + 120 * om * om;
            b  = 1 + 15 * om + 90  * om * om;
            c  = 1 + 60 * om + 360 * om * om;
            
            % Notation simplifications
            L2  = L*L;
            mu2 = mu * mu;
            
            % Assemble element geometric stiffness matrix in local system
            kgl = N/L *...
                  [ 1    0             0               -1    0             0;
                    0    6*a/(5*mu2)   L/(10*mu2)       0   -6*a/(5*mu2)   L/(10*mu2);
                    0    L/(10*mu2)    2*L2*b/(15*mu2)  0   -L/(10*mu2)   -L2*c/(30*mu2);
                   -1    0             0                1    0             0;
                    0   -6*a/(5*mu2)  -L/(10*mu2)       0    6*a/(5*mu2)  -L/(10*mu2);
                    0    L/(10*mu2)   -L2*c/(30*mu2)    0   -L/(10*mu2)    2*L2*b/(15*mu2) ];
            
            % Static condensation to consider beam end liberations
            kgl = elem.condenseEndLib(kgl);
        end
        
        %------------------------------------------------------------------
        % Compute element corotational tangent stiffness matrix in global system
        % considering Euler-Bernoulli beam theory.
        function kcr = corotStiffMtx_EB(elem)
            % Element properties
            E  = elem.material.E;
            Ax = elem.section.Ax;
            Iz = elem.section.Iz;
            L0 = elem.len_ini;
            L  = elem.len_cur;
            N  = elem.Fc(4);
            M1 = elem.Fc(3);
            M2 = elem.Fc(6);
            c  = cos(elem.ang_cur);
            s  = sin(elem.ang_cur);
            
            % Notation simplifications
            EA   = E*Ax;
            EI   = E*Iz;
            L2   = L*L;
            LL0  = L*L0;
            L2L0 = L2*L0;
            c2   = c*c;
            s2   = s*s;
            cs   = c*s;
            M12  = M1+M2;
            
            % Compute stiffness coefficients
            k1 = (EA*L2*c2 + N*L*L0*s2 - 2*M12*L0*cs + 12*EI*s2)/L2L0;
            k2 = (EA*L2*s2 + N*L*L0*c2 + 2*M12*L0*cs + 12*EI*c2)/L2L0;
            k3 = (cs*(12*EI - EA*L2 + N*LL0) - M12*(c2-s2)*L0)/L2L0;
            k4 = (6*EI*c)/LL0;
            k5 = (6*EI*s)/LL0;
            k6 = (2*EI)/L0;
            k7 = (4*EI)/L0;
            
            % Assemble tangent stiffness matrix in global system
            kcr = [ k1  -k3  -k5  -k1   k3  -k5 
                   -k3   k2   k4   k3  -k2   k4 
                   -k5   k4   k7   k5  -k4   k6 
                   -k1   k3   k5   k1  -k3   k5 
                    k3  -k2  -k4  -k3   k2  -k4 
                   -k5   k4   k6   k5  -k4   k7 ];
            
            % Static condensation to consider beam end liberations
            kcr = elem.condenseEndLib(kcr);
        end
        
        %------------------------------------------------------------------
        % Compute element corotational tangent stiffness matrix in global system
        % considering Timoshenko beam theory.
        function kcr = corotStiffMtx_Tim(elem)
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ax = elem.section.Ax;
            Ay = elem.section.Ay;
            Iz = elem.section.Iz;
            L0 = elem.len_ini;
            L  = elem.len_cur;
            N  = elem.Fc(4);
            M1 = elem.Fc(3);
            M2 = elem.Fc(6);
            c  = cos(elem.ang_cur);
            s  = sin(elem.ang_cur);
            
            % Timoshenko parameters
            om = E * Iz / (G * Ay * L0 * L0);
            mu = 1 + 12 * om;
            la = 1 + 3  * om;
            ga = 1 - 6  * om;
            
            % Notation simplifications
            c2     = c*c;
            s2     = s*s;
            cs     = c*s;
            EA     = E*Ax;
            EI     = E*Iz;
            L2     = L*L;
            LL0    = L*L0;
            L0mu   = L0*mu;
            LL0mu  = L*L0mu;
            L2L0mu = L2*L0mu;
            M12    = M1 + M2;
            g2l    = ga+2*la;
            
            % Compute stiffness coefficients
            k1 = (EA*L2*c2*mu + N*L*L0*s2*mu - 2*M12*L0*cs*mu + 4*EI*s2*g2l)/L2L0mu;
            k2 = (EA*L2*s2*mu + N*L*L0*c2*mu + 2*M12*L0*cs*mu + 4*EI*c2*g2l)/L2L0mu;
            k3 = (cs*(4*EI*g2l - EA*L2*mu + N*LL0*mu) - M12*(c2-s2)*L0*mu)/L2L0mu;
            k4 = 2*EI*c*g2l/LL0mu;
            k5 = 2*EI*s*g2l/LL0mu;
            k6 = (2*EI*ga)/L0mu;
            k7 = (4*EI*la)/L0mu;
            
            % Assemble tangent stiffness matrix in global system
            kcr = [ k1  -k3  -k5  -k1   k3  -k5
                   -k3   k2   k4   k3  -k2   k4
                   -k5   k4   k7   k5  -k4   k6
                   -k1   k3   k5   k1  -k3   k5
                    k3  -k2  -k4  -k3   k2  -k4
                   -k5   k4   k6   k5  -k4   k7 ];
            
            % Static condensation to consider beam end liberations
            kcr = elem.condenseEndLib(kcr);
        end
        
        %------------------------------------------------------------------
        % Compute element consistent mass matrix in local system
        % considering Euler-Bernoulli beam theory.
        function mcl = consistMassMtx_EB(elem)
            % Element properties
            rho = elem.material.rho;
            A   = elem.section.Ax;
            L   = elem.len_cur;
            L2  = L^2;
            M   = rho * A * L;
            
            % Assemble element consistent mass matrix in local system
            mcl = M/420*...
                  [ 140    0      0      70     0      0;
                    0      156    22*L   0      54    -13*L;
                    0      22*L   4*L2   0      13*L  -3*L2;
                    70     0      0      140    0      0;
                    0      54     13*L   0      156   -22*L;
                    0     -13*L  -3*L2   0     -22*L   4*L2 ];
            
            % Static condensation to consider beam end liberations
            mcl = elem.condenseEndLib(mcl);
        end
        
        %------------------------------------------------------------------
        % Compute element consistent mass matrix in local system
        % considering Timoshenko beam theory.
        function mcl = consistMassMtx_Tim(elem)
            % Element properties
            E   = elem.material.E;
            G   = elem.material.G;
            rho = elem.material.rho;
            Ax  = elem.section.Ax;
            Ay  = elem.section.Ay;
            Iz  = elem.section.Iz;
            L   = elem.len_cur;
            M   = rho * Ax * L;
            
            % Timoshenko parameters
            om = E * Iz / (G * Ay * L * L);
            mu = 1 + 12 * om;
            th = 1 + 4  * om;
            
            % Notation simplifications
            L2  = L  * L;
            mu2 = mu * mu;
            
            % Assemble element consistent mass matrix in local system
            m1 =  264 + 48  * th / mu - 48 * om / mu2;
            m2 = (44  - 108 * om / mu - 24 * om / mu2) * L;
            m3 =  108 + 360 * om / mu + 72 * om * th / mu2;
            m4 = (26  + 108 * om / mu + 24 * om / mu2) * L;
            m5 = (7 + 1 / mu2) * L2;
            m6 = (7 - 1 / mu2) * L2;
            
            mcl = M/840 *...
                  [ 280  0    0    140  0    0 ;
                    0    m1   m2   0    m3  -m4
                    0    m2   m5   0    m4  -m6 ;
                    140  0    0    280  0    0
                    0    m3   m4   0    m1  -m2 ;
                    0   -m4  -m6   0   -m2   m5 ];
            
            % Static condensation to consider beam end liberations
            mcl = elem.condenseEndLib(mcl);
        end
        
        %------------------------------------------------------------------
        % Assemble displacement shape function matrix evaluated in a given
        % element position, considering Euler-Bernoulli beam theory.
        function N = shapeFcnMtx_EB(elem,x)
            c  = Constants();
            L  = elem.len_cur;
            L2 = L  * L;
            L3 = L2 * L;
            x2 = x  * x;
            x3 = x2 * x;
            
            % Compute axial displacements
            Nu1 = 1 - x/L;
            Nu2 = x/L;
            
            % Compute transversal displacements
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 3*x2/L2 + 2*x3/L3;
                Nv2 = x - 2*x2/L + x3/L2;
                Nv3 = 3*x2/L2 - 2*x3/L3;
                Nv4 = -x2/L + x3/L2;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nv2 = 0;
                Nv3 = 3*x/(2*L) - x3/(2*L3);
                Nv4 = -x/2 + x3/(2*L2);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nv2 = x - 3*x2/(2*L) + x3/(2*L2);
                Nv3 = 3*x2/(2*L2) - x3/(2*L3);
                Nv4 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - x/L;
                Nv2 = 0;
                Nv3 = x/L;
                Nv4 = 0;
            end
            
            % Assemble displacement shape function matrix
            N = [ Nu1  0    0    Nu2  0    0;
                  0    Nv1  Nv2  0    Nv3  Nv4 ];
        end
        
        %------------------------------------------------------------------
        % Assemble displacement shape function matrix evaluated in a given
        % element position, considering Timoshenko beam theory.
        function N = shapeFcnMtx_Tim(elem,x)
            c  = Constants();
            
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ay = elem.section.Ay;
            I  = elem.section.Iz;
            L  = elem.len_cur;  
            
            % Timoshenko parameters
            om = E * I / (G * Ay * L * L);
            mu = 1 + 12 * om;
            la = 1 + 3  * om;
            ga = 1 - 6  * om;
            
            % Notation simplifications
            L2 = L^2;
            L3 = L2 * L;
            x2 = x^2;
            x3 = x * x2;
            
            % Compute axial displacements
            Nu1 = 1 - x/L;
            Nu2 = x/L;
            
            % Compute transversal displacements
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 12*om*x/(L*mu) - 3*x2/(L2*mu) + 2*x3/(L3*mu);
                Nv2 = x - 6*om*x/mu - 2*la*x2/(L*mu) + x3/(L2*mu);
                Nv3 = 12*om*x/(L*mu) + 3*x2/(L2*mu) - 2*x3/(L3*mu);
                Nv4 = -6*om*x/mu - ga*x2/(L*mu) + x3/(L2*mu);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 3*x/(2*L*la) - 3*om*x/(L*la) + x3/(2*L3*la);
                Nv2 = 0;
                Nv3 = 3*x/(2*L*la) + 3*om*x/(L*la) - x3/(2*L3*la);
                Nv4 = -ga*x/(2*la) - 3*om*x/la + x3/(2*L2*la);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - 3*om*x/(L*la) - 3*x2/(2*L2*la) + x3/(2*L3*la);
                Nv2 = x - 3*om*x/la - 3*x2/(2*L*la) + x3/(2*L2*la);
                Nv3 = 3*om*x/(L*la) + 3*x2/(2*L2*la) - x3/(2*L3*la);
                Nv4 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - x/L;
                Nv2 = 0;
                Nv3 = x/L;
                Nv4 = 0;
            end
            
            % Assemble displacement shape function matrix
            N = [ Nu1  0    0    Nu2  0    0;
                  0    Nv1  Nv2  0    Nv3  Nv4 ];
        end
    end
end