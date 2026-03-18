%% Anm_Frame3D Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anm to deal with 3D frame models.
%
classdef Anm_Frame3D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Frame3D()
            c = Constants();
            anm = anm@Anm(c.FRAME3D,6,[1,2,3,4,5,6]);
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
            z1 = elem.nodes(1).coord(3);
            x2 = elem.nodes(2).coord(1);
            y2 = elem.nodes(2).coord(2);
            z2 = elem.nodes(2).coord(3);
            
            % Calculate element length
            dx = x2 - x1;
            dy = y2 - y1;
            dz = z2 - z1;
            L  = sqrt(dx^2 + dy^2 + dz^2);
            
            % Calculate unity element local axes
            x = [dx,dy,dz];
            x = x/norm(x);
            y = cross(elem.vz,x);
            y = y/norm(y);
            z = cross(x,y);
            
            % Assemble basis rotation transformation matrix
            T = [ x(1) x(2) x(3);
                  y(1) y(2) y(3);
                  z(1) z(2) z(3) ];
            
            % Set element properties
            elem.len_ini = L;
            elem.len_cur = L;
            elem.T_ini   = T;
            elem.T       = T;
            elem.TN1     = T';
            elem.TN2     = T';
            elem.rot     = blkdiag(T,T,T,T);
            elem.Rm      = eye(3);
        end
        
        %------------------------------------------------------------------
        % Update element geometry (length, angles, rotation matrix) with
        % given displacements (increment and total).
        function updCurGeom(anm,elem,d_U,U)
            % Get increment and total element displacements in global system
            d_Ug = d_U(elem.gle);
            Ug   = U(elem.gle);
            
            % Evaluate new nodal coordinates
            x1 = elem.nodes(1).coord(1) + Ug(1);
            y1 = elem.nodes(1).coord(2) + Ug(2);
            z1 = elem.nodes(1).coord(3) + Ug(3);
            x2 = elem.nodes(2).coord(1) + Ug(7);
            y2 = elem.nodes(2).coord(2) + Ug(8);
            z2 = elem.nodes(2).coord(3) + Ug(9);
            
            % Assemble element local x-axis
            dx = x2 - x1;
            dy = y2 - y1;
            dz = z2 - z1;
            x = [dx;dy;dz];
            
            % Calculate element length
            L = norm(x);
            
            % Calculate mean rotation from last configuration
            rm = 0.5*[d_Ug(4)+d_Ug(10); d_Ug(5)+d_Ug(11); d_Ug(6)+d_Ug(12)];
            Rm = anm.largeRotMtx(anm,rm);
            Rm = Rm*elem.T';
            
            % Get columns of mean rotation matrix
            r1 = Rm(:,1);
            r2 = Rm(:,2);
            r3 = Rm(:,3);
            
            % Compute local axes vectors
            x = x/norm(x);
            y = r2 - (r2'*x)/(1+r1'*x) * (r1+x);
            z = r3 - (r3'*x)/(1+r1'*x) * (r1+x);
            
            % Ensure vectors have unit length
            y = y/norm(y);
            z = z/norm(z);
            
             % Assemble element basis rotation transformation matrix
            T = [ x(1) x(2) x(3);
                  y(1) y(2) y(3);
                  z(1) z(2) z(3) ];
            
            % Compute large rotation matrix of element nodes
            RN1 = anm.largeRotMtx(anm,d_Ug(4:6));
            RN2 = anm.largeRotMtx(anm,d_Ug(10:12));
            
            % Set element properties
            elem.len_cur = L;
            elem.T       = T;
            elem.TN1     = RN1 * elem.TN1;
            elem.TN2     = RN2 * elem.TN2;
            elem.rot     = blkdiag(T,T,T,T);
            elem.Rm      = Rm;
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
            qz0 = Q(3);
            qx1 = Q(4) - Q(1);
            qy1 = Q(5) - Q(2);
            qz1 = Q(6) - Q(3);
            
            % Compute element fixed-end-forces vector in local system
            if (elem.type == c.EULER_BERNOULLI)
                fel = anm.fixEndForceVctDistrib_EB(elem,qx0,qy0,qz0,qx1,qy1,qz1);
            elseif (elem.type == c.TIMOSHENKO)
                fel = anm.fixEndForceVctDistrib_Tim(elem,qx0,qy0,qz0,qx1,qy1,qz1);
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
        function fig = intForceVctNonLin(anm,elem,~)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            A  = elem.section.Ax;
            J  = elem.section.Ix;
            L0 = elem.len_ini;
            L1 = elem.len_cur;
            
            % Element local basis
            x = elem.T(1,:);
            y = elem.T(2,:);
            z = elem.T(3,:);
            
            % Nodes local basis
            nx1 = elem.TN1(:,1);
            ny1 = elem.TN1(:,2);
            nz1 = elem.TN1(:,3);
            nx2 = elem.TN2(:,1);
            ny2 = elem.TN2(:,2);
            nz2 = elem.TN2(:,3);
            
            % Length increment since beginning of analysis
            dL = L1 - L0;
            
            % Relative rotations
            rx  = asin(0.5*(z*ny2-y*nz2))-asin(0.5*(z*ny1-y*nz1));
            ry1 = asin(0.5*(x*nz1-z*nx1));
            rz1 = asin(0.5*(y*nx1-x*ny1));
            ry2 = asin(0.5*(x*nz2-z*nx2));
            rz2 = asin(0.5*(y*nx2-x*ny2));
            
            % Internal forces in local system
            N1 = -E*A*dL/L0;
            N2 =  E*A*dL/L0;
            
            T1 = -G*J*rx/L0;
            T2 =  G*J*rx/L0;
            
            if (elem.type == c.EULER_BERNOULLI)
                [My1,Mz1,My2,Mz2] = anm.endMoment_EB(elem,ry1,rz1,ry2,rz2);
            elseif (elem.type == c.TIMOSHENKO)
                [My1,Mz1,My2,Mz2] = anm.endMoment_Tim(elem,ry1,rz1,ry2,rz2);
            end
            
            Qy1 =  (Mz1+Mz2)/L0;
            Qz1 = -(My1+My2)/L0;
            Qy2 = -Qy1;
            Qz2 = -Qz1;
            
            % Assemble element internal force vector in local system
            fil = [N1; Qy1; Qz1; T1; My1; Mz1; N2; Qy2; Qz2; T2; My2; Mz2];
            
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
            
            % Get elastic stiffness matrix in natural system
            if (elem.type == c.EULER_BERNOULLI)
                ken = anm.natElastStiffMtx_EB(elem);
            elseif (elem.type == c.TIMOSHENKO)
                ken = anm.natElastStiffMtx_Tim(elem);
            end
            
            % Element length and internal forces
            L   = elem.len_cur;
            N   = elem.Fc(7);
            Mx1 = elem.Fc(4);
            My1 = elem.Fc(5);
            Mz1 = elem.Fc(6);
            Mx2 = elem.Fc(10);
            My2 = elem.Fc(11);
            Mz2 = elem.Fc(12);            
            
            % Element local basis
            x = elem.T(1,:);
            y = elem.T(2,:);
            z = elem.T(3,:);
            
            % Nodes local basis
            nx1 = elem.TN1(:,1);
            ny1 = elem.TN1(:,2);
            nz1 = elem.TN1(:,3);
            nx2 = elem.TN2(:,1);
            ny2 = elem.TN2(:,2);
            nz2 = elem.TN2(:,3);
            
            % Natural rotations
            rx1 = asin(0.5*(z*ny1-y*nz1));
            ry1 = asin(0.5*(x*nz1-z*nx1));
            rz1 = asin(0.5*(y*nx1-x*ny1));
            rx2 = asin(0.5*(z*ny2-y*nz2));
            ry2 = asin(0.5*(x*nz2-z*nx2));
            rz2 = asin(0.5*(y*nx2-x*ny2));
            
            % Columns of mean rotation matrix
            r1 = elem.Rm(:,1);
            r2 = elem.Rm(:,2);
            r3 = elem.Rm(:,3);
            
            % Skew matrices
            Sx   = anm.skewMtx(x);
            Sy   = anm.skewMtx(y);
            Sz   = anm.skewMtx(z);
            Sr1  = anm.skewMtx(r1);
            Sr2  = anm.skewMtx(r2);
            Sr3  = anm.skewMtx(r3);
            Snx1 = anm.skewMtx(nx1);
            Sny1 = anm.skewMtx(ny1);
            Snz1 = anm.skewMtx(nz1);
            Snx2 = anm.skewMtx(nx2);
            Sny2 = anm.skewMtx(ny2);
            Snz2 = anm.skewMtx(nz2);
            
            % Auxiliary matrix A
            A = (eye(3)-x'*x)/L;
            
            % Axiliary matrices Lr
            L21 = 0.50 * (r2'*x'*A + A*r2*(x+r1'));
            L31 = 0.50 * (r3'*x'*A + A*r3*(x+r1'));
            L22 = 0.25 * (2*Sr2 - r2'*x'*Sr1 - Sr2*x'*(x+r1'));
            L32 = 0.25 * (2*Sr3 - r3'*x'*Sr1 - Sr3*x'*(x+r1'));
            
            Lr2 = [L21; L22; -L21; L22];
            Lr3 = [L31; L32; -L31; L32];
            
            % Auxiliary matrices H
            Hx1 = [ 0;0;0;  Sny1*z'-Snz1*y';  0;0;0;  0;0;0 ];
            Hy1 = [ A*nz1;  Snx1*z'-Snz1*x'; -A*nz1;  0;0;0 ];
            Hz1 = [ A*ny1;  Snx1*y'-Sny1*x'; -A*ny1;  0;0;0 ];
            Hx2 = [ 0;0;0;  0;0;0;            0;0;0;  Sny2*z'-Snz2*y' ];
            Hy2 = [ A*nz2;  0;0;0;           -A*nz2;  Snx2*z'-Snz2*x' ];
            Hz2 = [ A*ny2;  0;0;0;           -A*ny2;  Snx2*y'-Sny2*x' ];
            
            % Lines of transformation matrix from natural system to global system
            t1  =  [-x 0 0 0 x 0 0 0];
            t2  =  (Lr2*nx1 + Hz1)' / (2*cos(rz1));
            t3  =  (Lr2*nx2 + Hz2)' / (2*cos(rz2));
            t4  = -(Lr3*nx1 + Hy1)' / (2*cos(ry1));
            t5  = -(Lr3*nx2 + Hy2)' / (2*cos(ry2));
            t6a =  (Lr3*ny1 - Lr2*nz1 + Hx1)' / (2*cos(rx1));
            t6b =  (Lr3*ny2 - Lr2*nz2 + Hx2)' / (2*cos(rx2));
            t6  =  t6b - t6a;
            
            % Assemble transformation matrix from natural system to global system
            T = [t1;t2;t3;t4;t5;t6];
            
            % Scaled basic forces
            m2  = Mz1/(2*cos(rz1));
            m3  = Mz2/(2*cos(rz2));
            m4  = My1/(2*cos(ry1));
            m5  = My2/(2*cos(ry2));
            m6a = Mx2/(2*cos(rx1));
            m6b = Mx2/(2*cos(rx2));
            
            % Auxiliary matrices M
            Mny1 = anm.mMtx(ny1,L,A,x);
            Mnz1 = anm.mMtx(nz1,L,A,x);
            Mny2 = anm.mMtx(ny2,L,A,x);
            Mnz2 = anm.mMtx(nz2,L,A,x);
            
            % Axiliary matrices G
            G1 = anm.gMtx(anm,r2,nx1,L,A,x,r1);
            G2 = anm.gMtx(anm,r2,nz1,L,A,x,r1);
            G3 = anm.gMtx(anm,r2,nx2,L,A,x,r1);
            G4 = anm.gMtx(anm,r2,nz2,L,A,x,r1);
            G5 = anm.gMtx(anm,r3,nx1,L,A,x,r1);
            G6 = anm.gMtx(anm,r3,ny1,L,A,x,r1);
            G7 = anm.gMtx(anm,r3,nx2,L,A,x,r1);
            G8 = anm.gMtx(anm,r3,ny2,L,A,x,r1);
            
            % Compute geometric stiffness submatrices
            kg1 = N * [ A         zeros(3) -A         zeros(3);
                        zeros(3)  zeros(3)  zeros(3)  zeros(3);
                       -A         zeros(3)  A         zeros(3);
                        zeros(3)  zeros(3)  zeros(3)  zeros(3) ];
            
            kg2 = Mz1*(t2'*t2)*tan(rz1)  +...
                  Mz2*(t3'*t3)*tan(rz2)  +...
                  My1*(t4'*t4)*tan(ry1)  +...
                  My2*(t5'*t5)*tan(ry2)  +...
                  Mx1*(t6b'*t6b*tan(rx2) - t6a'*t6a*tan(rx1));
            
            kg3 = m2*G1 + m3*G3 - m4*G5 - m5*G7 + m6a*(G8-G4) - m6b*(G6-G2);
            
            kg4a = -Lr2*(m2*Snx1 + m6a*Snz1) + Lr3*(m4*Snx1 + m6a*Sny1);
            kg4b = -Lr2*(m3*Snx2 - m6b*Snz2) + Lr3*(m5*Snx2 - m6b*Sny2);
            kg4  = [zeros(12,3)  kg4a  zeros(12,3)  kg4b];
            
            kg5 = kg4';
            
            kg6_11 = -m2*Mny1 - m3*Mny2 + m4*Mnz1 + m5*Mnz2;
            kg6_12 = -m2*A*Sny1 + m4*A*Snz1;
            kg6_14 = -m3*A*Sny2 + m5*A*Snz2;
            kg6_22 =  m2*(Sy*Snx1-Sx*Sny1) - m4*(Sz*Snx1-Sx*Snz1) - m6a*(Sz*Sny1-Sy*Snz1);
            kg6_44 =  m3*(Sy*Snx2-Sx*Sny2) - m5*(Sz*Snx2-Sx*Snz2) + m6b*(Sz*Sny2-Sy*Snz2);
            
            kg6 = [ kg6_11    kg6_12   -kg6_11    kg6_14;
                    kg6_12'   kg6_22   -kg6_12'   zeros(3);
                   -kg6_11   -kg6_12    kg6_11   -kg6_14;
                    kg6_14'   zeros(3) -kg6_14'   kg6_44 ];
            
            % Assemble geometric stiffness matrix in global system
            kgg = kg1 + kg2 + kg3 + kg4 + kg5 + kg6;
            
            % Assemble tangent stiffness matrix in global system
            kcr = T'*ken*T + kgg;
            
            % Static condensation to consider beam end liberations
            kcr = elem.condenseEndLib(kcr);
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
                  [ 1  0  0  0  0  0  0  0  0  0  0  0;
                    0  1  0  0  0  0  0  0  0  0  0  0;
                    0  0  1  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0;    
                    0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  1  0  0  0  0  0;
                    0  0  0  0  0  0  0  1  0  0  0  0;
                    0  0  0  0  0  0  0  0  1  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0 ];
            
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
        function fel = fixEndForceVctDistrib_EB(elem,qx0,qy0,qz0,qx1,qy1,qz1)
            c = Constants();
            
            % Element length
            L  = elem.len_cur;
            L2 = L^2;
            
            % Compute axial fixed-end-forces
            N1 = -(qx0*L/2 + qx1*L/6);
            N2 = -(qx0*L/2 + qx1*L/3);
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Qy1 = -(qy0*L/2   + qy1*3*L/20);
                Qz1 = -(qz0*L/2   + qz1*3*L/20);
                My1 =   qz0*L2/12 + qz1*L2/30;
                Mz1 = -(qy0*L2/12 + qy1*L2/30);
                Qy2 = -(qy0*L/2   + qy1*7*L/20);
                Qz2 = -(qz0*L/2   + qz1*7*L/20);
                My2 = -(qz0*L2/12 + qz1*L2/20);
                Mz2 =   qy0*L2/12 + qy1*L2/20;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Qy1 = -(qy0*3*L/8 + qy1*L/10);
                Qz1 = -(qz0*3*L/8 + qz1*L/10);
                My1 =   0;
                Mz1 =   0;
                Qy2 = -(qy0*5*L/8 + qy1*2*L/5);
                Qz2 = -(qz0*5*L/8 + qz1*2*L/5);
                My2 = -(qz0*L2/8  + qz1*L2/15);
                Mz2 =   qy0*L2/8  + qy1*L2/15;
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Qy1 = -(qy0*5*L/8 + qy1*9*L/40);
                Qz1 = -(qz0*5*L/8 + qz1*9*L/40);
                My1 =   qz0*L2/8  + qz1*7*L2/120;
                Mz1 = -(qy0*L2/8  + qy1*7*L2/120);
                Qy2 = -(qy0*3*L/8 + qy1*11*L/40);
                Qz2 = -(qz0*3*L/8 + qz1*11*L/40);
                My2 =   0;
                Mz2 =   0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Qy1 = -(qy0*L/2 + qy1*L/6);
                Qz1 = -(qz0*L/2 + qz1*L/6);
                My1 =   0;
                Mz1 =   0;
                Qy2 = -(qy0*L/2 + qy1*L/3);
                Qz2 = -(qz0*L/2 + qz1*L/3);
                My2 =   0;
                Mz2 =   0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Qy1; Qz1; 0; My1; Mz1; N2; Qy2; Qz2; 0; My2; Mz2];
        end
        
        %------------------------------------------------------------------
        % Compute element fixed-end-forces vector in local system for a
        % distributed load considering Timoshenko beam theory.
        function fel = fixEndForceVctDistrib_Tim(elem,qx0,qy0,qz0,qx1,qy1,qz1)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ay = elem.section.Ay;
            Az = elem.section.Az;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            L2 = L^2;
            
            % Compute axial fixed-end-forces
            N1 = -(qx0*L/2 + qx1*L/6);
            N2 = -(qx0*L/2 + qx1*L/3);
            
            % Timoshenko parameters for bending about z direction
            omz      = E * Iz / (G * Ay * L2);
            mu       = 1 + 12 * omz;
            lambda   = 1 + 3  * omz;
            zeta     = 1 + 40 * omz/3;
            xi       = 1 + 5  * omz;
            eta      = 1 + 15 * omz;
            vartheta = 1 + 4  * omz;
            psi      = 1 + 12 * omz/5;
            varpi    = 1 + 20 * omz/9;
            epsilon  = 1 + 80 * omz/7;
            varrho   = 1 + 10 * omz;
            upsilon  = 1 + 5  * omz/2;
            varsigma = 1 + 40 * omz/11;
            
            % Compute flexural fixed-end-forces for bending about z direction
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Qy1 = -(qy0*L/2   + qy1*3*L*zeta/(20*mu));
                Mz1 = -(qy0*L2/12 + qy1*L2*eta/(30*mu));
                Qy2 = -(qy0*L/2   + qy1*7*L*epsilon/(20*mu));
                Mz2 =   qy0*L2/12 + qy1*L2*varrho/(20*mu);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Qy1 = -(qy0*3*L*vartheta/(8*lambda) + qy1*L*xi/(10*lambda));
                Mz1 =   0;
                Qy2 = -(qy0*5*L*psi/(8*lambda)      + qy1*2*L*upsilon/(5*lambda));
                Mz2 =   qy0*L2/(8*lambda)           + qy1*L2/(15*lambda);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Qy1 = -(qy0*5*L*psi/(8*lambda)      + qy1*9*L*varpi/(40*lambda));
                Mz1 = -(qy0*L2/(8*lambda)           + qy1*7*L2/(120*lambda));
                Qy2 = -(qy0*3*L*vartheta/(8*lambda) + qy1*11*L*varsigma/(40*lambda));
                Mz2 =   0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Qy1 = -(qy0*L/2 + qy1*L/6);
                Mz1 =   0;
                Qy2 = -(qy0*L/2 + qy1*L/3);
                Mz2 =   0;
            end
            
            % Timoshenko parameters for bending about y direction
            omy      = E * Iy / (G * Az * L2);
            mu       = 1 + 12 * omy;
            lambda   = 1 + 3  * omy;
            zeta     = 1 + 40 * omy/3;
            xi       = 1 + 5  * omy;
            eta      = 1 + 15 * omy;
            vartheta = 1 + 4  * omy;
            psi      = 1 + 12 * omy/5;
            varpi    = 1 + 20 * omy/9;
            epsilon  = 1 + 80 * omy/7;
            varrho   = 1 + 10 * omy;
            upsilon  = 1 + 5  * omy/2;
            varsigma = 1 + 40 * omy/11;
            
            % Compute flexural fixed-end-forces for bending about y direction
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Qz1 = -(qz0*L/2   + qz1*3*L*zeta/(20*mu));
                My1 =   qz0*L2/12 + qz1*L2*eta/(30*mu);
                Qz2 = -(qz0*L/2   + qz1*7*L*epsilon/(20*mu));
                My2 = -(qz0*L2/12 + qz1*L2*varrho/(20*mu));
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Qz1 = -(qz0*3*L*vartheta/(8*lambda) + qz1*L*xi/(10*lambda));
                My1 =   0;
                Qz2 = -(qz0*5*L*psi/(8*lambda)      + qz1*2*L*upsilon/(5*lambda));
                My2 = -(qz0*L2/(8*lambda)           + qz1*L2/(15*lambda));
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Qz1 = -(qz0*5*L*psi/(8*lambda)      + qz1*9*L*varpi/(40*lambda));
                My1 =   qz0*L2/(8*lambda)           + qz1*7*L2/(120*lambda);
                Qz2 = -(qz0*3*L*vartheta/(8*lambda) + qz1*11*L*varsigma/(40*lambda));
                My2 =   0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Qz1 = -(qz0*L/2 + qz1*L/6);
                My1 =   0;
                Qz2 = -(qz0*L/2 + qz1*L/3);
                My2 =   0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Qy1; Qz1; 0; My1; Mz1; N2; Qy2; Qz2; 0; My2; Mz2];
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
            Iy    = elem.section.Iy;
            Iz    = elem.section.Iz;
            Hy    = elem.section.Hy;
            Hz    = elem.section.Hz;
            L     = elem.len_cur;
            EIy   = E*Iy;
            EIz   = E*Iz;
            
            % Temperature variations
            dtx = elem.thermload(1);
            dty = elem.thermload(2);
            dtz = elem.thermload(3);
            
            % Unitary dimensionless temperature gradient
            tgy = -(alpha * dty) / Hy;
            tgz = -(alpha * dtz) / Hz;
            
            % Compute axial fixed-end-forces
            N1 =  E * A * alpha * dtx;
            N2 = -E * A * alpha * dtx;
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Qy1 =  0;
                Qz1 =  0;
                My1 = -tgz * EIy;
                Mz1 =  tgy * EIz;
                Qy2 =  0;
                Qz2 =  0;
                My2 =  tgz * EIy;
                Mz2 = -tgy * EIz;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Qy1 = -tgy * 3*EIz/(2*L);
                Qz1 = -tgz * 3*EIy/(2*L);
                My1 =  0;
                Mz1 =  0;
                Qy2 =  tgy * 3*EIz/(2*L);
                Qz2 =  tgz * 3*EIy/(2*L);
                My2 =  tgz * 3*EIy/2;
                Mz2 = -tgy * 3*EIz/2;
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Qy1 =  tgy * 3*EIz/(2*L);
                Qz1 =  tgz * 3*EIy/(2*L);
                My1 = -tgz * 3*EIy/2;
                Mz1 =  tgy * 3*EIz/2;
                Qy2 = -tgy * 3*EIz/(2*L);
                Qz2 = -tgz * 3*EIy/(2*L);
                My2 =  0;
                Mz2 =  0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Qy1 = 0;
                Qz1 = 0;
                My1 = 0;
                Mz1 = 0;
                Qy2 = 0;
                Qz2 = 0;
                My2 = 0;
                Mz2 = 0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Qy1; Qz1; 0; My1; Mz1; N2; Qy2; Qz2; 0; My2; Mz2];
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
            Az    = elem.section.Az;
            Iy    = elem.section.Iy;
            Iz    = elem.section.Iz;
            Hy    = elem.section.Hy;
            Hz    = elem.section.Hz;
            L     = elem.len_cur;
            L2    = L*L;
            EIy   = E*Iy;
            EIz   = E*Iz;
            GAy   = G*Ay;
            GAz   = G*Az;
            
            % Timoshenko parameters
            omy = EIy / (GAz * L2);
            omz = EIz / (GAy * L2);
            lay = 1 + 3  * omy;
            laz = 1 + 3  * omz;
            
            % Temperature variations
            dtx = elem.thermload(1);
            dty = elem.thermload(2);
            dtz = elem.thermload(3);
            
            % Unitary dimensionless temperature gradient
            tgy = -(alpha * dty) / Hy;
            tgz = -(alpha * dtz) / Hz;
            
            % Compute axial fixed-end-forces
            N1 =  E * Ax * alpha * dtx;
            N2 = -E * Ax * alpha * dtx;
            
            % Compute flexural fixed-end-forces
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS) 
                Qy1 =  0;
                Qz1 =  0;
                My1 = -tgz * EIy;
                Mz1 =  tgy * EIz;
                Qy2 =  0;
                Qz2 =  0;
                My2 =  tgz * EIy;
                Mz2 = -tgy * EIz;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Qy1 = -tgy * 3*EIz/(2*L*laz);
                Qz1 = -tgz * 3*EIy/(2*L*lay);
                My1 =  0;
                Mz1 =  0;
                Qy2 =  tgy * (3*EIz/(2*L*laz));
                Qz2 =  tgz * (3*EIy/(2*L*lay));
                My2 =  tgz * (3*EIy/(2*lay));
                Mz2 = -tgy * (3*EIz/(2*laz));
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Qy1 =  tgy * 3*EIz/(2*L*laz);
                Qz1 =  tgz * 3*EIy/(2*L*lay);
                My1 = -tgz * 3*EIy/(2*lay);
                Mz1 =  tgy * 3*EIz/(2*laz);
                Qy2 = -tgy * 3*EIz/(2*L*laz);
                Qz2 = -tgz * 3*EIy/(2*L*lay);
                My2 =  0;
                Mz2 =  0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Qy1 =  0;
                Qz1 =  0;
                My1 =  0;
                Mz1 =  0;
                Qy2 =  0;
                Qz2 =  0;
                My2 =  0;
                Mz2 =  0;
            end
            
            % Assemble element fixed-end-force vector in local system
            fel = [N1; Qy1; Qz1; 0; My1; Mz1; N2; Qy2; Qz2; 0; My2; Mz2];
        end
        
        %------------------------------------------------------------------
        % Compute element end moments considering Euler-Bernoulli beam theory.
        function [My1,Mz1,My2,Mz2] = endMoment_EB(elem,ry1,rz1,ry2,rz2)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_ini;
            
            % Compute end moments
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS)
                My1 = 2*E*Iy*(2*ry1+ry2)/L;
                Mz1 = 2*E*Iz*(2*rz1+rz2)/L;
                My2 = 2*E*Iy*(ry1+2*ry2)/L;
                Mz2 = 2*E*Iz*(rz1+2*rz2)/L;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                My1 = 0;
                Mz1 = 0;
                My2 = 3*E*Iy*ry2/L;
                Mz2 = 3*E*Iz*rz2/L;
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                My1 = 3*E*Iy*ry1/L;
                Mz1 = 3*E*Iz*rz1/L;
                My2 = 0;
                Mz2 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                My1 = 0;
                Mz1 = 0;
                My2 = 0;
                Mz2 = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Compute element end moments considering Timoshenko beam theory.
        function [My1,Mz1,My2,Mz2] = endMoment_Tim(elem,ry1,rz1,ry2,rz2)
            c = Constants();
            
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ay = elem.section.Ay;
            Az = elem.section.Az;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_ini;
            
            % Timoshenko parameters
            omz = E * Iz / (G * Ay * L * L);
            omy = E * Iy / (G * Az * L * L);
            niz = 1 + 12 * omz;
            niy = 1 + 12 * omy;
            laz = 1 + 3  * omz;
            lay = 1 + 3  * omy;
            gaz = 1 - 6  * omz;
            gay = 1 - 6  * omy;
            
            % Compute end moments
            if (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.CONTINUOUS)
                My1 = 2*E*Iy*(2*lay*ry1+gay*ry2)/(niy*L);
                Mz1 = 2*E*Iz*(2*laz*rz1+gaz*rz2)/(niz*L);
                My2 = 2*E*Iy*(gay*ry1+2*lay*ry2)/(niy*L);
                Mz2 = 2*E*Iz*(gaz*rz1+2*laz*rz2)/(niz*L);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                My1 = 0;
                Mz1 = 0;
                My2 = 3*E*Iy*lay*ry2/(niy*L);
                Mz2 = 3*E*Iz*laz*rz2/(niz*L);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                My1 = 3*E*Iy*lay*ry1/(niy*L);
                Mz1 = 3*E*Iz*laz*rz1/(niz*L);
                My2 = 0;
                Mz2 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                My1 = 0;
                Mz1 = 0;
                My2 = 0;
                Mz2 = 0;
            end
        end
        
        %------------------------------------------------------------------
        % Compute element elastic stiffness matrix in local system
        % considering Euler-Bernoulli beam theory.
        function kel = elastStiffMtx_EB(elem)
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ax = elem.section.Ax;
            Ix = elem.section.Ix;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            
            % Notation simplifications
            L2  = L^2;
            L3  = L2 * L;
            EA  = E * Ax;
            GIx = G * Ix;
            EIy = E * Iy;
            EIz = E * Iz;
            
            % Assemble element elastic stiffness matrix in local system
            kel = [ EA/L    0          0         0       0          0        -EA/L    0          0         0       0         0;
                    0       12*EIz/L3  0         0       0          6*EIz/L2  0      -12*EIz/L3  0         0       0         6*EIz/L2;
                    0       0          12*EIy/L3 0      -6*EIy/L2   0         0       0         -12*EIy/L3 0      -6*EIy/L2  0;
                    0       0          0         GIx/L   0          0         0       0          0        -GIx/L   0         0;
                    0       0         -6*EIy/L2  0       4*EIy/L    0         0       0          6*EIy/L2  0       2*EIy/L   0;
                    0       6*EIz/L2   0         0       0          4*EIz/L   0      -6*EIz/L2   0         0       0         2*EIz/L;
                   -EA/L    0          0         0       0          0         EA/L    0          0         0       0         0;
                    0      -12*EIz/L3  0         0       0         -6*EIz/L2  0       12*EIz/L3  0         0       0        -6*EIz/L2;
                    0       0         -12*EIy/L3 0       6*EIy/L2   0         0       0          12*EIy/L3 0       6*EIy/L2  0;
                    0       0          0        -GIx/L    0         0         0       0          0         GIx/L   0         0;
                    0       0         -6*EIy/L2  0       2*EIy/L    0         0       0          6*EIy/L2  0       4*EIy/L   0;
                    0       6*EIz/L2   0         0       0          2*EIz/L   0      -6*EIz/L2   0         0       0         4*EIz/L ];
            
            % Static condensation to consider beam end liberations
            kel = elem.condenseEndLib(kel);
        end
        
        %------------------------------------------------------------------
        % Compute element elastic stiffness matrix in local system
        % considering Timoshenko beam theory.
        function kel = elastStiffMtx_Tim(elem)
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ax = elem.section.Ax;
            Ay = elem.section.Ay;
            Az = elem.section.Az;
            Ix = elem.section.Ix;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            
            % Timoshenko parameters
            omz = E * Iz / (G * Ay * L * L);
            omy = E * Iy / (G * Az * L * L);
            muz = 1 + 12 * omz;
            muy = 1 + 12 * omy;
            laz = 1 + 3  * omz;
            lay = 1 + 3  * omy;
            gaz = 1 - 6  * omz;
            gay = 1 - 6  * omy;
            
            % Notation simplifications
            EA    = E * Ax;
            GIx   = G * Ix;
            EIz   = E * Iz;
            EIy   = E * Iy;
            muzL  = muz   * L;
            muzL2 = muzL  * L;
            muzL3 = muzL2 * L;
            muyL  = muy   * L;
            muyL2 = muyL  * L;
            muyL3 = muyL2 * L;
            
            % Assemble element elastic stiffness matrix in local system
            kel = [ EA/L  0             0             0      0               0              -EA/L  0             0             0      0               0;
                    0     12*EIz/muzL3  0             0      0               6*EIz/muzL2     0    -12*EIz/muzL3  0             0      0               6*EIz/muzL2;
                    0     0             12*EIy/muyL3  0     -6*EIy/muyL2     0               0     0            -12*EIy/muyL3  0     -6*EIy/muyL2     0;
                    0     0             0             GIx/L  0               0               0     0             0            -GIx/L  0               0;
                    0     0            -6*EIy/muyL2   0      4*lay*EIy/muyL  0               0     0             6*EIy/muyL2   0      2*gay*EIy/muyL  0;
                    0     6*EIz/muzL2   0             0      0               4*laz*EIz/muzL  0    -6*EIz/muzL2   0             0      0               2*gaz*EIz/muzL;
                   -EA/L  0             0             0      0               0               EA/L  0             0             0      0               0;
                    0    -12*EIz/muzL3  0             0      0              -6*EIz/muzL2     0     12*EIz/muzL3  0             0      0              -6*EIz/muzL2;
                    0     0            -12*EIy/muyL3  0      6*EIy/muyL2     0               0     0             12*EIy/muyL3  0      6*EIy/muyL2     0;
                    0     0             0            -GIx/L  0               0               0     0             0             GIx/L  0               0;
                    0     0            -6*EIy/muyL2   0      2*gay*EIy/muyL  0               0     0             6*EIy/muyL2   0      4*lay*EIy/muyL  0;
                    0     6*EIz/muzL2   0             0      0               2*gaz*EIz/muzL  0    -6*EIz/muzL2   0             0      0               4*laz*EIz/muzL ];
            
            % Static condensation to consider beam end liberations
            kel = elem.condenseEndLib(kel);
        end
        
        %------------------------------------------------------------------
        % Compute element geometric stiffness matrix in local system
        % considering Euler-Bernoulli beam theory.
        function kgl = geomStiffMtx_EB(elem)
            % Element properties
            A  = elem.section.Ax;
            Ix = elem.section.Ix;
            N  = elem.Fc(7);
            L  = elem.len_cur;
            L2 = L*L;
            
            % Assemble element geometric stiffness matrix in local system
            kgl = N/L *...
                  [ 1        0        0        0        0        0       -1        0        0        0        0        0;
                    0        6/5      0        0        0        L/10     0       -6/5      0        0        0        L/10;
                    0        0        6/5      0       -L/10     0        0        0       -6/5      0       -L/10     0;
                    0        0        0        Ix/A     0        0        0        0        0       -Ix/A     0        0;
                    0        0       -L/10     0        2*L2/15  0        0        0        L/10     0       -L2/30    0;    
                    0        L/10     0        0        0        2*L2/15  0       -L/10     0        0        0       -L2/30;
                   -1        0        0        0        0        0        1        0        0        0        0        0;
                    0       -6/5      0        0        0       -L/10     0        6/5      0        0        0       -L/10;
                    0        0       -6/5      0        L/10     0        0        0        6/5      0        L/10     0;
                    0        0        0       -Ix/A     0        0        0        0        0        Ix/A     0       0;
                    0        0       -L/10     0       -L2/30    0        0        0        L/10     0        2*L2/15  0
                    0        L/10     0        0        0       -L2/30    0       -L/10     0        0        0        2*L2/15];
            
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
            Ax = elem.section.Ax;
            Ay = elem.section.Ay;
            Az = elem.section.Az;
            Ix = elem.section.Ix;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_cur;
            N  = elem.Fc(7);
            
            % Timoshenko parameters
            omz = E * Iz / (G * Ay * L * L);
            omy = E * Iy / (G * Az * L * L);
            muz = 1 + 12 * omz;
            muy = 1 + 12 * omy;
            az  = 1 + 20 * omz + 120 * omz * omz;
            ay  = 1 + 20 * omy + 120 * omy * omy;
            bz  = 1 + 15 * omz + 90  * omz * omz;
            by  = 1 + 15 * omy + 90  * omy * omy;
            cz  = 1 + 60 * omz + 360 * omz * omz;
            cy  = 1 + 60 * omy + 360 * omy * omy;
            
            % Notation simplifications
            L2   = L*L;
            muz2 = muz * muz;
            muy2 = muy * muy;
            
            % Assemble element geometric stiffness matrix in local system
            kgl = N/L *...
                  [ 1       0              0              0       0                  0                -1       0              0              0        0                 0;
                    0       6*az/(5*muz2)  0              0       0                  L/(10*muz2)       0      -6*az/(5*muz2)  0              0        0                 L/(10*muz2);
                    0       0              6*ay/(5*muy2)  0      -L/(10*muy2)        0                 0       0             -6*ay/(5*muy2)  0       -L/(10*muy2)       0;
                    0       0              0              Ix/Ax   0                  0                 0       0              0             -Ix/Ax    0                 0;
                    0       0             -L/(10*muy2)    0       2*L2*by/(15*muy2)  0                 0       0              L/(10*muy2)    0       -L2*cy/(30*muy2)   0;    
                    0       L/(10*muz2)    0              0       0                  2*L2*bz/(15*muz2) 0      -L/(10*muz2)    0              0        0                -L2*cz/(30*muz2);
                   -1       0              0              0       0                  0                 1       0              0              0        0                 0;
                    0      -6*az/(5*muz2)  0              0       0                 -L/(10*muz2)       0       6*az/(5*muz2)  0              0        0                -L/(10*muz2);
                    0       0             -6*ay/(5*muy2)  0       L/(10*muy2)        0                 0       0              6*ay/(5*muy2)  0        L/(10*muy2)       0;
                    0       0              0             -Ix/Ax   0                  0                 0       0              0              Ix/Ax    0                 0;
                    0       0             -L/(10*muy2)    0      -L2*cy/(30*muy2)    0                 0       0              L/(10*muy2)    0        2*L2*by/(15*muy2) 0
                    0       L/(10*muz2)    0              0       0                 -L2*cz/(30*muz2)   0      -L/(10*muz2)    0              0        0                 2*L2*bz/(15*muz2) ];
            
            % Static condensation to consider beam end liberations
            kgl = elem.condenseEndLib(kgl);
        end
        
        %------------------------------------------------------------------
        % Compute element elastic stiffness matrix in natural system
        % considering Euler-Bernoulli beam theory.
        function ken = natElastStiffMtx_EB(elem)
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ax = elem.section.Ax;
            Ix = elem.section.Ix;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L0 = elem.len_ini;
            
            % Notation simplifications
            EA  = E*Ax;
            GIx = G*Ix;
            EIy = E*Iy;
            EIz = E*Iz;
            
            % Assemble elastic stiffness matrix in natural system
            ken = [ EA/L0     0         0         0         0         0;
                    0         4*EIz/L0  2*EIz/L0  0         0         0;
                    0         2*EIz/L0  4*EIz/L0  0         0         0;
                    0         0         0         4*EIy/L0  2*EIy/L0  0;
                    0         0         0         2*EIy/L0  4*EIy/L0  0;
                    0         0         0         0         0         GIx/L0 ];
        end
        
        %------------------------------------------------------------------
        % Compute element elastic stiffness matrix in natural system
        % considering Timoshenko beam theory.
        function ken = natElastStiffMtx_Tim(elem)
            % Element properties
            E  = elem.material.E;
            G  = elem.material.G;
            Ax = elem.section.Ax;
            Ay = elem.section.Ay;
            Az = elem.section.Az;
            Ix = elem.section.Ix;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L0 = elem.len_ini;
            
            % Timoshenko parameters
            omz = E * Iz / (G * Ay * L0 * L0);
            omy = E * Iy / (G * Az * L0 * L0);
            muz = 1 + 12 * omz;
            muy = 1 + 12 * omy;
            laz = 1 + 3  * omz;
            lay = 1 + 3  * omy;
            gaz = 1 - 6  * omz;
            gay = 1 - 6  * omy;
            
            % Notation simplifications
            EA    = E   * Ax;
            GIx   = G   * Ix;
            EIz   = E   * Iz;
            EIy   = E   * Iy;
            muzL0 = muz * L0;
            muyL0 = muy * L0;
            
            % Assemble elastic stiffness matrix in natural system
            ken = [ EA/L0     0                0                0                0                0;
                    0         4*laz*EIz/muzL0  2*gaz*EIz/muzL0  0                0                0;
                    0         2*gaz*EIz/muzL0  4*laz*EIz/muzL0  0                0                0;
                    0         0                0                4*lay*EIy/muyL0  2*gay*EIy/muyL0  0;
                    0         0                0                2*gay*EIy/muyL0  4*lay*EIy/muyL0  0;
                    0         0                0                0                0                GIx/L0 ];
        end
        
        %------------------------------------------------------------------
        % Compute element consistent mass matrix in local system
        % considering Euler-Bernoulli beam theory.
        function mcl = consistMassMtx_EB(elem)
            % Element properties
            rho = elem.material.rho;
            Ax  = elem.section.Ax;
            Ix  = elem.section.Ix;
            L   = elem.len_cur;
            L2  = L^2;
            M   = rho * Ax * L;
            At  = Ix / Ax;
            
            % Assemble element elastic stiffness matrix in local system
            mcl = M/420*...
                  [ 140    0      0      0      0      0      70     0      0      0      0      0;
                    0      156    0      0      0      22*L   0      54     0      0      0     -13*L;
                    0      0      156    0     -22*L   0      0      0      54     0      13*L   0;
                    0      0      0      140*At 0      0      0      0      0      70*At  0      0;
                    0      0     -22*L   0      4*L2   0      0      0     -13*L   0     -3*L2   0;    
                    0      22*L   0      0      0      4*L2   0      13*L   0      0      0     -3*L2;
                    70     0      0      0      0      0      140    0      0      0      0      0;
                    0      54     0      0      0      13*L   0      156    0      0      0     -22*L;
                    0      0      54     0     -13*L   0      0      0      156    0      22*L   0;
                    0      0      0      70*At  0      0      0      0      0      140*At 0      0;
                    0      0      13*L   0     -3*L2   0      0      0      22*L   0      4*L2   0;
                    0     -13*L   0      0      0     -3*L2   0     -22*L   0      0      0      4*L2 ];
            
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
            Az  = elem.section.Az;
            Ix  = elem.section.Ix;
            Iy  = elem.section.Iy;
            Iz  = elem.section.Iz;
            L   = elem.len_cur;
            M   = rho * Ax * L;
            At  = Ix / Ax;
            
            % Timoshenko parameters for bending about z direction
            omz = E * Iz / (G * Ay * L * L);
            muz = 1 + 12 * omz;
            thz = 1 + 4  * omz;
            
            % Timoshenko parameters for bending about y direction
            omy = E * Iy / (G * Az * L * L);
            muy = 1 + 12 * omy;
            thy = 1 + 4  * omy;
            
            % Notation simplifications
            L2   = L  * L;
            muz2 = muz * muz;
            muy2 = muy * muy;
            
            % Assemble element consistent mass matrix in local system
            m1z =  264 + 48  * thz / muz - 48 * omz / muz2;
            m1y =  264 + 48  * thy / muy - 48 * omy / muy2;
            m2z = (44  - 108 * omz / muz - 24 * omz / muz2) * L;
            m2y = (44  - 108 * omy / muy - 24 * omy / muy2) * L;
            m3z =  108 + 360 * omz / muz + 72 * omz * thz / muz2;
            m3y =  108 + 360 * omy / muy + 72 * omy * thz / muy2;
            m4z = (26  + 108 * omz / muz + 24 * omz / muz2) * L;
            m4y = (26  + 108 * omy / muy + 24 * omy / muy2) * L;
            m5z = (7 + 1 / muz2) * L2;
            m5y = (7 + 1 / muy2) * L2;
            m6z = (7 - 1 / muz2) * L2;
            m6y = (7 - 1 / muy2) * L2;
            
            mcl = M/840 *...
                  [ 280    0      0      0      0      0      140    0      0      0      0      0;
                    0      m1z    0      0      0      m2z    0      m3z    0      0      0      -m4z;
                    0      0      m1y    0     -m2y    0      0      0      m3y    0      m4y    0;
                    0      0      0      280*At 0      0      0      0      0      140*At 0      0;
                    0      0     -m2y    0      m5y    0      0      0     -m4y    0     -m6y    0;    
                    0      m2z    0      0      0      m5z    0      m4z    0      0      0      -m6z;
                    140    0      0      0      0      0      280    0      0      0      0      0;
                    0      m3z    0      0      0      m4z    0      m1z    0      0      0      -m2z;
                    0      0      m3y    0     -m4y    0      0      0      m1y    0      m2y    0;
                    0      0      0      140*At 0      0      0      0      0      280*At 0      0;
                    0      0      m4y    0     -m6y    0      0      0      m2y    0      m5y    0;
                    0     -m4z    0      0      0     -m6z    0     -m2z    0      0      0      m5z ];
            
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
                Nw1 = 1 - 3*x2/L2 + 2*x3/L3;
                Nv2 = x - 2*x2/L + x3/L2;
                Nw2 = -x + 2*x2/L - x3/L2;
                Nv3 = 3*x2/L2 - 2*x3/L3;
                Nw3 = 3*x2/L2 - 2*x3/L3;
                Nv4 = -x2/L + x3/L2;
                Nw4 = x2/L - x3/L2;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nw1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nv2 = 0;
                Nw2 = 0;
                Nv3 = 3*x/(2*L) - x3/(2*L3);
                Nw3 = 3*x/(2*L) - x3/(2*L3);
                Nv4 = -x/2 + x3/(2*L2);
                Nw4 = x/2 - x3/(2*L2);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nw1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nv2 = x - 3*x2/(2*L) + x3/(2*L2);
                Nw2 = -x + 3*x2/(2*L) - x3/(2*L2);
                Nv3 = 3*x2/(2*L2) - x3/(2*L3);
                Nw3 = 3*x2/(2*L2) - x3/(2*L3);
                Nv4 = 0;
                Nw4 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - x/L;
                Nw1 = 1 - x/L;
                Nv2 = 0;
                Nw2 = 0;
                Nv3 = x/L;
                Nw3 = x/L;
                Nv4 = 0;
                Nw4 = 0;
            end
            
            % Assemble displacement shape function matrix
            N = [ Nu1  0    0    0    0    0    Nu2  0    0    0    0    0;
                  0    Nv1  0    0    0    Nv2  0    Nv3  0    0    0    Nv4;
                  0    0    Nw1  0    Nw2  0    0    0    Nw3  0    Nw4  0 ];
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
            Az = elem.section.Az;
            Iy = elem.section.Iy;
            Iz = elem.section.Iz;
            L  = elem.len_cur;  
            
            % Timoshenko parameters
            omz = E * Iz / (G * Ay * L * L);
            omy = E * Iy / (G * Az * L * L);
            muz = 1 + 12 * omz;
            muy = 1 + 12 * omy;
            laz = 1 + 3  * omz;
            lay = 1 + 3  * omy;
            gaz = 1 - 6  * omz;
            gay = 1 - 6  * omy;
            
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
                Nv1 = 1 - 12*omz*x/(L*muz) - 3*x2/(L2*muz) + 2*x3/(L3*muz);
                Nw1 = 1 - 12*omy*x/(L*muy) - 3*x2/(L2*muy) + 2*x3/(L3*muy);
                Nv2 = x - 6*omz*x/muz - 2*laz*x2/(L*muz) + x3/(L2*muz);
                Nw2 = -x + 6*omy*x/muy + 2*lay*x2/(L*muy) - x3/(L2*muy);
                Nv3 = 12*omz*x/(L*muz) + 3*x2/(L2*muz) - 2*x3/(L3*muz);
                Nw3 = 12*omy*x/(L*muy) + 3*x2/(L2*muy) - 2*x3/(L3*muy);
                Nv4 = -6*omz*x/muz - gaz*x2/(L*muz) + x3/(L2*muz);
                Nw4 = 6*omy*x/muy + gay*x2/(L*muy) - x3/(L2*muy);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 3*x/(2*L*laz) - 3*omz*x/(L*laz) + x3/(2*L3*laz);
                Nw1 = 1 - 3*x/(2*L*lay) - 3*omy*x/(L*lay) + x3/(2*L3*lay);
                Nv2 = 0;
                Nw2 = 0;
                Nv3 = 3*x/(2*L*laz) + 3*omz*x/(L*laz) - x3/(2*L3*laz);
                Nw3 = 3*x/(2*L*lay) + 3*omy*x/(L*lay) - x3/(2*L3*lay);
                Nv4 = -gaz*x/(2*laz) - 3*omz*x/laz + x3/(2*L2*laz);
                Nw4 = gay*x/(2*lay) + 3*omy*x/lay - x3/(2*L2*lay);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - 3*omz*x/(L*laz) - 3*x2/(2*L2*laz) + x3/(2*L3*laz);
                Nw1 = 1 - 3*omy*x/(L*lay) - 3*x2/(2*L2*lay) + x3/(2*L3*lay);
                Nv2 = x - 3*omz*x/laz - 3*x2/(2*L*laz) + x3/(2*L2*laz);
                Nw2 = -x + 3*omy*x/lay + 3*x2/(2*L*lay) - x3/(2*L2*lay);
                Nv3 = 3*omz*x/(L*laz) + 3*x2/(2*L2*laz) - x3/(2*L3*laz);
                Nw3 = 3*omy*x/(L*lay) + 3*x2/(2*L2*lay) - x3/(2*L3*lay);
                Nv4 = 0;
                Nw4 = 0;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - x/L;
                Nw1 = 1 - x/L;
                Nv2 = 0;
                Nw2 = 0;
                Nv3 = x/L;
                Nw3 = x/L;
                Nv4 = 0;
                Nw4 = 0;
            end
            
            % Assemble displacement shape function matrix
            N = [ Nu1  0    0    0    0    0    Nu2  0    0    0    0    0;
                  0    Nv1  0    0    0    Nv2  0    Nv3  0    0    0    Nv4;
                  0    0    Nw1  0    Nw2  0    0    0    Nw3  0    Nw4  0 ];
        end
        
        %------------------------------------------------------------------
        % Assemble skew matrix from a given vector.
        function S = skewMtx(v)
            S = [ 0    -v(3)  v(2);
                  v(3)  0    -v(1);
                 -v(2)  v(1)  0 ];
        end
        
        %------------------------------------------------------------------
        % Compute large rotation matrix about a rotational pseudo-vector
        % using Rodrigues formula.
        function R = largeRotMtx(anm,d_r)
            r = norm(d_r);
            if (r ~= 0)
                S = anm.skewMtx(d_r);
                R = eye(3) + sin(r)/r * S + (1-cos(r))/r^2 * S^2;
            else
                R = eye(3);
            end
        end
        
        %------------------------------------------------------------------
        % Compute auxiliary matrix G for corotational tangent stiffness matrix.
        function G = gMtx(anm,r,z,L,A,x,r1)
            Mr = anm.mMtx(r,L,A,x);
            Mz = anm.mMtx(z,L,A,x);
            Sr = anm.skewMtx(r);
            Sz = anm.skewMtx(z);
            Sx = anm.skewMtx(x);
            S1 = anm.skewMtx(r1);
            
            g11 = -0.500*(A*z*r'*A + A*r*z'*A + r'*x'*Mz + (x+r1')*z*Mr);
            g12 = -0.250*(A*z*x*Sr + (x+r1')*z*A*Sr + A*r*z'*S1);
            g22 =  0.125*(-(r'*x')*Sz*S1 + Sr*x'*z'*S1 + 2*Sz*Sr + S1*z*x*Sr - (x+r1')*z*Sx*Sr);
            
            G = [ g11   g12  -g11   g12;
                  g12'  g22  -g12'  g22;
                 -g11  -g12   g11  -g12;
                  g12'  g22  -g12'  g22 ];
        end
        
        %------------------------------------------------------------------
        % Compute auxiliary matrix M for corotational tangent stiffness matrix.
        function M = mMtx(z,L,A,x)
            M = -1/L * (A*z*x + (A*z*x)' + A*(x*z));
        end
    end
end