%% Anm_Truss3D Class
%
% This is a sub-class in the NUMA-TF program that implements abstract 
% methods declared in super-class Anm to deal with 3D truss models.
%
classdef Anm_Truss3D < Anm
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function anm = Anm_Truss3D()
            c = Constants();
            anm = anm@Anm(c.TRUSS3D,3,[1,2,3]);
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
            elem.T       = T;
            elem.rot     = blkdiag(T,T);
        end
        
        %------------------------------------------------------------------
        % Update element geometry (length, angles, rotation matrix) with
        % given displacements (increment and total).
        function updCurGeom(~,elem,~,U)
            % Get total element displacements in global system
            Ug = U(elem.gle);
            
            % Evaluate new nodal coordinates
            x1 = elem.nodes(1).coord(1) + Ug(1);
            y1 = elem.nodes(1).coord(2) + Ug(2);
            z1 = elem.nodes(1).coord(3) + Ug(3);
            x2 = elem.nodes(2).coord(1) + Ug(4);
            y2 = elem.nodes(2).coord(2) + Ug(5);
            z2 = elem.nodes(2).coord(3) + Ug(6);
            
            % Assemble element local x-axis
            dx = x2 - x1;
            dy = y2 - y1;
            dz = z2 - z1;
            x  = [dx,dy,dz];
            
            % Calculate element length
            L = norm(x);
            
            % Initialize element local y-axis
            y = elem.T(2,:);
            
            % Calculate z-axis and recompute y-axis
            z = cross(x,y);
            y = cross(z,x);
            
            % Normalize local axes vectors
            x = x/norm(x);
            y = y/norm(y);
            z = z/norm(z);
            
            % Assemble basis rotation transformation matrix
            T = [ x(1) x(2) x(3);
                  y(1) y(2) y(3);
                  z(1) z(2) z(3) ];
            
            % Set element properties
            elem.len_cur = L;
            elem.T       = T;
            elem.rot     = blkdiag(T,T);
        end
        
        %------------------------------------------------------------------
        % Assemble element equivalent nodal load vector in global system
        % for a generic linear distributed load.
        function feg = equivLoadVctDistrib(~,~)
            % Element distributed loads are not considered on truss models
            feg = [0;0;0;0;0;0];
        end
        
        %------------------------------------------------------------------
        % Assemble element equivalent nodal load vector in global system
        % for a thermal load.
        function feg = equivLoadVctThermal(~,~)
            % Element thermal loads are not considered on truss models
            feg = [0;0;0;0;0;0];
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
        function fig = intForceVctNonLin(~,elem,~)
            % Element properties
            E  = elem.material.E;
            A  = elem.section.Ax;
            L0 = elem.len_ini;
            L1 = elem.len_cur;
            
            % Length increment since beginning of analysis
            dL = L1 - L0;
            
            % Internal forces in local system
            N1 = -E*A*dL/L0;
            N2 =  E*A*dL/L0;
            
            % Assemble element internal force vector in local system
            fil = [ N1; 0; 0; N2; 0; 0 ];
            
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
        function keg = elastStiffMtx(~,elem)
            % Element properties
            E = elem.material.E;
            A = elem.section.Ax;
            L = elem.len_cur;
            
            % Assemble element elastic stiffness matrix in local system
            kel = [ E*A/L   0       0      -E*A/L   0      0;
                    0       0       0       0       0      0;
                    0       0       0       0       0      0;
                   -E*A/L   0       0       E*A/L   0      0;
                    0       0       0       0       0      0;
                    0       0       0       0       0      0 ];
            
            % Transform element matrix from local to global system
            keg = elem.rot' * kel * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble element geometric stiffness matrix in global system.
        function kgg = geomStiffMtx(~,elem)
            % Element properties
            L = elem.len_cur;
            N = elem.Fc(4);
            
            % Assemble element geometric stiffness matrix in local system
            kgl = [ 0       0       0       0       0       0;
                    0       N/L     0       0      -N/L     0;
                    0       0       N/L     0       0      -N/L;
                    0       0       0       0       0       0;
                    0      -N/L     0       0       N/L     0;
                    0       0      -N/L     0       0       N/L];
            
            % Transform element matrix from local to global system
            kgg = elem.rot' * kgl * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble element corotational stiffness matrix in global system.
        function kcr = corotStiffMtx(~,elem)
            % Element properties
            E  = elem.material.E;
            A  = elem.section.Ax;
            L0 = elem.len_ini;
            L1 = elem.len_cur;
            N  = elem.Fc(4);
            EA = E*A;
            
            % Assemble element elastic stiffness matrix in local system
            kel = [ EA/L0  0       0      -EA/L0  0       0;
                    0      0       0       0      0       0;
                    0      0       0       0      0       0;
                   -EA/L0  0       0       EA/L0  0       0;
                    0      0       0       0      0       0;
                    0      0       0       0      0       0 ];
            
            % Assemble element geometric stiffness matrix in local system
            kgl = [ 0      0       0       0      0       0;
                    0      N/L1    0       0     -N/L1    0;
                    0      0       N/L1    0      0      -N/L1;
                    0      0       0       0      0       0;
                    0     -N/L1    0       0      N/L1    0;
                    0      0      -N/L1    0      0       N/L1];
            
            % Compute element tangent stiffness matrix in local system
            ktl = kel + kgl;
            
            % Transform element matrix from local to global system
            kcr = elem.rot' * ktl * elem.rot;
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
        function mcg = consistMassMtx(~,elem)
            % Element properties
            rho = elem.material.rho;
            A   = elem.section.Ax;
            L   = elem.len_cur;
            M   = rho * A * L;
            
            % Assemble element consistent mass matrix in local system
            mcl = M/6 *...
                  [ 2  0  0  1  0  0;
                    0  0  0  0  0  0;
                    0  0  0  0  0  0;
                    1  0  0  2  0  0;
                    0  0  0  0  0  0;
                    0  0  0  0  0  0 ];
            
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
                    0  0  0  0  0  0;
                    0  0  0  0  0  0;
                    0  0  0  1  0  0;
                    0  0  0  0  0  0;
                    0  0  0  0  0  0 ];
            
            % Transform element matrix from local to global system
            mlg = elem.rot' * mll * elem.rot;
        end
        
        %------------------------------------------------------------------
        % Assemble displacement shape function matrix evaluated in a given
        % element position.
        function N = shapeFcnMtx(~,elem,x)
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
                Nv2 = 3*x2/L2 - 2*x3/L3;
                Nw2 = 3*x2/L2 - 2*x3/L3;
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.CONTINUOUS)
                Nv1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nw1 = 1 - 3*x/(2*L) + x3/(2*L3);
                Nv2 = 3*x/(2*L) - x3/(2*L3);
                Nw2 = 3*x/(2*L) - x3/(2*L3);
            elseif (elem.fixend(1) == c.CONTINUOUS && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nw1 = 1 - 3*x2/(2*L2) + x3/(2*L3);
                Nv2 = 3*x2/(2*L2) - x3/(2*L3);
                Nw2 = 3*x2/(2*L2) - x3/(2*L3);
            elseif (elem.fixend(1) == c.HINGED && elem.fixend(2) == c.HINGED)
                Nv1 = 1 - x/L;
                Nw1 = 1 - x/L;
                Nv2 = x/L;
                Nw2 = x/L;
            end
            
            % Assemble displacement shape function matrix
            N = [ Nu1  0    0    Nu2  0    0;
                  0    Nv1  0    0    Nv2  0;
                  0    0    Nw1  0    0    Nw2 ];
        end
    end
end