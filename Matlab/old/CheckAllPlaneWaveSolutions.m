%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModelPlaneWaveUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a single
% internal wave for all possible wavenumbers.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% November 17th, 2016   Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 8;
Ny = 8;
Nz = 8;

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wavemodel = InternalWaveModelSlow([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = 0.01; % m/s

for k_loop=(-Nx/2 + 1):1:(Nx/2-1)
    for l_loop=(-Ny/2 + 1):1:(Ny/2-1)
        fprintf('(k0,l0)=(%d,%d) ',k_loop,l_loop);
        for j0=1:Nz/2
            for sign=-1:2:1
                wavemodel = InternalWaveModel([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
                wavemodel.InitializeWithPlaneWave(k_loop,l_loop,j0,U,sign);
                
                t = 4*86400;
                [u,v] = wavemodel.VelocityFieldAtTime(t);
                [w,zeta] = wavemodel.VerticalFieldsAtTime(t);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Create a single plane-wave with the known analytical solution
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                f0 = wavemodel.f0;
                k = wavemodel.k;
                l = wavemodel.l;
                h = wavemodel.h;
                x = wavemodel.x;
                y = wavemodel.y;
                X = wavemodel.X;
                Y = wavemodel.Y;
                Z = wavemodel.Z;
                
                if (k_loop < 0)
                    k0 = Nx + k_loop;
                else
                    k0 = k_loop;
                end
                if (l_loop < 0)
                    l0 = Ny + l_loop;
                else
                    l0 = l_loop;
                end
                omega = sign*2*pi/wavemodel.period;
                m = j0*pi/Lz;
                % u_unit = U*cos( k(k0+1)*X + omega*t ).*cos(m*Z);
                % v_unit = -(f0/omega)*U*sin( k(k0+1)*X + omega*t ) .* cos(m*Z);
                % w_unit = (U*k(k0+1)/m) * sin(k(k0+1)*X + omega*t) .* sin(m*Z);
                % zeta_unit = -(U*k(k0+1)/m/omega) * cos(k(k0+1)*X + omega*t) .* sin(m*Z);
                
                alpha=atan2(l(l0+1),k(k0+1));
                K = sqrt( k(k0+1)^2 + l(l0+1)^2);
                u_unit = U*(cos(alpha)*cos( k(k0+1)*X + l(l0+1)*Y + omega*t ) + (f0/omega)*sin(alpha)*sin( k(k0+1)*X + l(l0+1)*Y + omega*t )).*cos(m*Z);
                v_unit = U*(sin(alpha)*cos( k(k0+1)*X + l(l0+1)*Y + omega*t ) - (f0/omega)*cos(alpha)*sin( k(k0+1)*X + l(l0+1)*Y + omega*t )).*cos(m*Z);
                w_unit = (U*K/m) * sin(k(k0+1)*X + l(l0+1)*Y + omega*t) .* sin(m*Z);
                zeta_unit = -(U*K/m/omega) * cos(k(k0+1)*X + l(l0+1)*Y + omega*t) .* sin(m*Z);
                
                % Compute the relative error
                max_speed = max(max(max( sqrt(u.*u + v.*v) )));
                max_u = max( [max(max(max( u ))), 1e-15] );
                u_error = max(max(max(abs(u-u_unit)/max_u)));
                max_v = max( [max(max(max( v ))), 1e-15] );
                v_error = max(max(max(abs(v-v_unit)/max_v)));
                max_w = max( [max(max(max( abs(w) ))), 1e-15] );
                w_error = max( [max(max(max(abs(w-w_unit)/max_w))), 1e-15] );
                max_zeta = max( [max(max(max( zeta ))), 1e-15] );
                zeta_error = max( [max(max(max(abs(zeta-zeta_unit)/max_zeta))), 1e-15] );
                
                max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(zeta_error)))]);
                
                if max_error > -3
                    if sign > 0
                        fprintf('\nFound at large error at +(k,l,j)=(%d,%d,%d):\n',k_loop,l_loop,j0);
                    else
                        fprintf('\nFound at large error at -(k,l,j)=(%d,%d,%d):\n',k_loop,l_loop,j0);
                    end
                    fprintf('The model solution for (u,v) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
                    fprintf('The model solution for (w,zeta) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(w_error))), round((log10(zeta_error))),t);
                end
            end
        end
    end
    fprintf('\n');
end
