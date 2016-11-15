%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModel
%
% This implements a simple internal wave model for constant stratification.
%
% The usage is simple. First call,
%   wavemodel = InternalWaveModel(dims, n, latitude, N0);
% to initialize the model with,
%   dims        a vector containing the length scales of x,y,z
%   n           a vector containing the number of grid points of x,y,z
%   latitude    the latitude of the model (e.g., 45)
%   N0          the buoyancy frequency of the stratification
%
% You must now intialize the model by calling either,
%   wavemodel.InitializeWithPlaneWave(k0, l0, j0, UAmp, sign);
% or
%   wavemodel.InitializeWithGMSpectrum(Amp);
% where Amp sets the relative GM amplitude.
%
% Finally, you can compute u,v,w,zeta at time t by calling,
%   [u,v] = wavemodel.VelocityFieldAtTime(t);
%   [w,zeta] = wavemodel.VerticalFieldsAtTime(t);
%
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% March 25th, 2016      Version 1.0
% March 30th, 2016      Version 1.1
%

classdef InternalWaveModel < handle
    properties
        Lx, Ly, Lz % Domain size
        Nx, Ny, Nz % Number of points in each direction
        latitude
        N0
        
        x, y, z
        k, l, j
        X,Y,Z
        K,L,J
        
        K2, Kh, F, G, M, h, Omega, f0
        u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, zeta_plus, zeta_minus
        period
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = InternalWaveModel(dims, n, latitude, N0)
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            obj.Lx = dims(1);
            obj.Ly = dims(2);
            obj.Lz = dims(3);
            
            obj.Nx = n(1);
            obj.Ny = n(2);
            obj.Nz = n(3);
            
            obj.latitude = latitude;
            obj.N0 = N0;
            
            dx = obj.Lx/obj.Nx;
            dy = obj.Ly/obj.Ny;
            dz = obj.Lz/obj.Nz;
            
            obj.x = dx*(0:obj.Nx-1)'; % periodic basis
            obj.y = dy*(0:obj.Ny-1)'; % periodic basis
            obj.z = dz*(0:obj.Nz-1)'; % cosine basis (not your usual dct basis, however)
            
            % Spectral domain, in radians
            dk = 1/obj.Lx;          % fourier frequency
            obj.k = 2*pi*([0:ceil(obj.Nx/2)-1 -floor(obj.Nx/2):-1]*dk)';
            
            dl = 1/obj.Ly;          % fourier frequency
            obj.l = 2*pi*([0:ceil(obj.Ny/2)-1 -floor(obj.Ny/2):-1]*dl)';
            
            nModes = obj.Nz;
            obj.j = (1:nModes)';
            
            % Using meshgrid means that,
            % K varies across dim 1
            % L varies across dim 2
            % J varies across dim 3
            [L,K,J] = meshgrid(obj.l,obj.k,obj.j);
            [Y,X,Z] = meshgrid(obj.y,obj.x,obj.z);
            
            obj.L = L; obj.K = K; obj.J = J;
            obj.X = X; obj.Y = Y; obj.Z = Z;
            
            obj.InitializeWaveProperties();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Begin initializing the wave field (internal use only)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWaveProperties(obj)
            obj.f0 = 2 * 7.2921E-5 * sin( obj.latitude*pi/180 );
            g = 9.81;
            
            obj.M = obj.J*pi/obj.Lz;        % Vertical wavenumber
            obj.K2 = obj.K.*obj.K + obj.L.*obj.L;   % Square of the horizontal wavenumber
            obj.Kh = sqrt(obj.K2);
            
            C2 = (obj.N0*obj.N0-obj.f0*obj.f0)./(obj.M.*obj.M+obj.K2);
            C = sqrt( C2 );                         % Mode speed
            obj.h = C2/g;                               % Mode depth
            obj.Omega = sqrt(C.*C.*obj.K2 + obj.f0*obj.f0);         % Mode frequency
            
            % F contains the coefficients for the U-V modes
            obj.F = (obj.h.*obj.M)*sqrt(2*g/(obj.Lz*(obj.N0*obj.N0-obj.f0*obj.f0)));
            
            % G contains the coefficients for the W-modes
            obj.G = sqrt(2*g/(obj.Lz*(obj.N0*obj.N0-obj.f0*obj.f0)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithPlaneWave(obj, k0, l0, j0, UAmp, sign)
            imagAmp = 1;
            if (k0 < -obj.Nx/2 || k0 > obj.Nx/2)
                error('Invalid choice for k0. Must be an integer %d <= k0 <= %d',-obj.Nx/2,obj.Nx/2);
            elseif (k0 < 0)
                k0 = abs(k0); % Because 'MakeHermitian' only grabs from the lower half, otherwise k0 = Nx + k0;
                imagAmp = -sqrt(-1);
            end
       
            if (l0 < -obj.Ny/2 || l0 > obj.Ny/2)
                error('Invalid choice for l0. Must be an integer %d <= l0 <= %d',-obj.Ny/2,obj.Ny/2);
            elseif (l0 < 0)
                l0 = obj.Ny + l0;
            end
            
            if (j0 < 1 || j0 >= obj.Nz)
                error('Invalid choice for j0. Must be an integer 0 < j < %d', obj.Nz);
            end
            
            myH = obj.h(k0+1,l0+1,j0);
            m = j0*pi/obj.Lz;
            g = 9.81;
            F_coefficient = myH * m * sqrt(2*g/obj.Lz)/sqrt(obj.N0^2 - obj.f0^2);
            
            U = zeros(size(obj.K));
            U(k0+1,l0+1,j0) = imagAmp*UAmp*sqrt(myH)/F_coefficient/2;
            if sign > 0
                A_plus = -MakeHermitian(U); % Careful with this, it doubles energy for some k,l
                A_plus(1,1,:) = 2*A_plus(1,1,:); % Double the zero frequency
                A_plus(obj.Nx/2+1,1,:) = -2*A_plus(obj.Nx/2+1,1,:); % Double the Nyquist frequency
                A_plus(1,obj.Ny/2+1,:) = -2*A_plus(1,obj.Ny/2+1,:); % Double the Nyquist frequency
                A_minus = zeros(size(U));
            else
                A_plus = zeros(size(U)); % Careful with this, it doubles energy for some k,l
                A_minus = MakeHermitian(U);
                A_minus(1,1,:) = 2*A_minus(1,1,:); % Double the zero frequency
                A_minus(obj.Nx/2+1,1,:) = -2*A_minus(obj.Nx/2+1,1,:); % Double the Nyquist frequency
                A_minus(1,obj.Ny/2+1,:) = -2*A_minus(1,obj.Ny/2+1,:); % Double the Nyquist frequency
                A_minus(1,1,:) = 0; % Intertial motions go only one direction!
            end
            
            obj.GenerateWavePhases(A_plus,A_minus);

            % This is the empirical method for setting the maximum U speed
            % to the desired value. For now we use the analytical method.
%             [u,~] = obj.VelocityFieldAtTime(0);
%             max_u = max(max(max( u )));
%    
%             obj.u_plus = obj.u_plus*(UAmp/max_u);
%             obj.u_minus = obj.u_minus*(UAmp/max_u);
%             obj.v_plus = obj.v_plus*(UAmp/max_u);
%             obj.v_minus = obj.v_minus*(UAmp/max_u);
%             obj.w_plus = obj.w_plus*(UAmp/max_u);
%             obj.w_minus = obj.w_minus*(UAmp/max_u);
%             obj.zeta_plus = obj.zeta_plus*(UAmp/max_u);
%             obj.zeta_minus = obj.zeta_minus*(UAmp/max_u);
            
            omega = obj.Omega(k0+1,l0+1,j0);
            obj.period = 2*pi/omega;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithGMSpectrum(obj, Amp)
            j_star = 3;
           
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*Amp;
                      
            % Compute the proper normalization with lots of modes
            j2 = (1:512);
            H = 2*(j_star.^(3/2))./(j2+j_star).^(5/2);
            H_norm = sum(H);
            
            % Analytical function for the GM spectrum in one horizontal
            % direction.
            D = obj.Lz;
%             GM2D_function = @(k,j) E*( (2/H_norm)*(j_star.^(3/2))./(j+j_star).^(5/2) ) .* (2/pi)*obj.f0*(j*pi/D).*(j*pi/D).*sqrt( (obj.N0*obj.N0-obj.f0*obj.f0)./(k.*k+(j*pi/D).*(j*pi/D)))./(obj.N0*obj.N0*k.*k+obj.f0*obj.f0*(j*pi/D).*(j*pi/D));
            alpha2 = obj.N0*obj.N0/(obj.f0*obj.f0);
            GM2D_analytic_int = @(k0,k1,j) E*( (2/H_norm)*(j_star.^(3/2))./(j+j_star).^(5/2) ) .* (2/pi)*(atan(k1/(j*pi/D)*sqrt( (alpha2-1)/(k0*k0/((j*pi/D)*(j*pi/D)) + 1))) - atan(k0/(j*pi/D)*sqrt( (alpha2-1)/(k0*k0/((j*pi/D)*(j*pi/D)) + 1)))) ;
            
            % Now create a wavenumber basis that uses the smallest
            % increment, and only goes to the smallest max wavenumber
            dk = (obj.k(2)-obj.k(1));
            dl = (obj.l(2)-obj.l(1));
            if (dk < dl)
                dm = dk;
            else
                dm = dl;
            end
            if ( max(obj.k) < max(obj.l) )
                m_max = max(obj.k);
            else
                m_max = max(obj.l);
            end
            m = (0:dm:m_max)';
            
            % Isotropically spread out the energy, but ignore the stuff
            % beyond m_max for now (we could come back and fix this).
            GM3D = zeros(size(obj.Kh));
            fprintf('Manually distributing the GM spectrum isotropically in horizontal wavenumber space. This may take a minute or two...\n');
            for mode=1:max(obj.j)
                fullIndices = false(size(GM3D));
%                 func = @(k) GM2D_function(k,mode);
%                 GM3D(1,1,mode) = integral(func,m(1),m(1)+dm/2);
                GM3D(1,1,mode) = GM2D_analytic_int(m(1),m(1)+dm/2,mode);
                for i=2:length(m)
                    m_lower = m(i)-dm/2;
                    m_upper = m(i)+dm/2;
                    indices = obj.Kh(:,:,mode) >= m_lower & obj.Kh(:,:,mode) < m_upper;
                    n = sum(sum(sum(indices)));
%                     total = integral(func,m_lower,m_upper)/n;
                    total = GM2D_analytic_int(m_lower,m_upper,mode)/n;
                    fullIndices(:,:,mode) = indices;
                    GM3D(fullIndices) = total;
                end
            end
            
            
            
            A = sqrt(GM3D/2); % Now split this into even and odd.
            
            A_plus = A.*GenerateHermitianRandomMatrix( size(obj.K) );
            A_minus = A.*GenerateHermitianRandomMatrix( size(obj.K) );
            A_minus(1,1,:) = 0; % Intertial motions go only one direction!
            
            GM_sum = sum(sum(sum(GM3D)))/E;
            GM_random_sum = sum(sum(sum(A_plus.*conj(A_plus) + A_minus.*conj(A_minus)  )))/E;
            fprintf('The coefficients sum to %.2f GM given the scales, and the randomized field sums to %f GM\n', GM_sum, GM_random_sum);
            
            obj.GenerateWavePhases(A_plus,A_minus);
        end
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GenerateWavePhases(obj, U_plus, U_minus)
            
%             denominator = obj.Kh.*obj.Omega.*sqrt(obj.h);
%             obj.u_plus = U_plus .* MakeHermitian((sqrt(-1)*obj.L*obj.f0 - obj.K.*obj.Omega)./denominator) .* obj.F;
%             obj.u_minus = U_minus .* MakeHermitian((sqrt(-1)*obj.L*obj.f0 + obj.K.*obj.Omega)./denominator) .* obj.F;
%             
%             obj.v_plus = U_plus .* MakeHermitian((-sqrt(-1)*obj.K*obj.f0 - obj.L.*obj.Omega)./denominator) .* obj.F;
%             obj.v_minus = U_minus .* MakeHermitian( (-sqrt(-1)*obj.K*obj.f0 + obj.L.*obj.Omega)./denominator) .* obj.F;
            
            theta = atan2(obj.L,obj.K);
            denominator = obj.Omega.*sqrt(obj.h);
            obj.u_plus = U_plus .* MakeHermitian( (sqrt(-1)*obj.f0 .* sin(theta) - obj.Omega .* cos(theta))./denominator ) .* obj.F;
            obj.u_minus = U_minus .* MakeHermitian( (sqrt(-1)*obj.f0 .* sin(theta) + obj.Omega .* cos(theta))./denominator ) .* obj.F;
            
            obj.v_plus = U_plus .* MakeHermitianNoInertial( ( -sqrt(-1)*obj.f0 .* cos(theta) - obj.Omega .* sin(theta) )./denominator ) .* obj.F;
            obj.v_minus = U_minus .* MakeHermitianNoInertial( ( -sqrt(-1)*obj.f0 .* cos(theta) + obj.Omega .* sin(theta) )./denominator ) .* obj.F;
            
            obj.w_plus = U_plus .* MakeHermitianNoInertial(sqrt(-1) *  obj.Kh .* sqrt(obj.h) ) * obj.G;
            obj.w_minus = U_minus .* MakeHermitianNoInertial( -sqrt(-1) * obj.Kh .* sqrt(obj.h) ) * obj.G;
            
            obj.zeta_plus = U_plus .* MakeHermitian( obj.Kh .* sqrt(obj.h) ./ obj.Omega ) * obj.G;
            obj.zeta_minus = U_minus .* MakeHermitian( obj.Kh .* sqrt(obj.h) ./ obj.Omega ) * obj.G;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u,v] = VelocityFieldAtTime(obj, t)
            phase_plus = exp(sqrt(-1)*obj.Omega*t);
            phase_minus = exp(-sqrt(-1)*obj.Omega*t);
            u_bar = MakeHermitianOnlyInertial( obj.u_plus.*phase_plus + obj.u_minus.*phase_minus );
            v_bar = MakeHermitianOnlyInertial( obj.v_plus.*phase_plus + obj.v_minus.*phase_minus );
            
            % Re-order to convert to an fast cosine transform
            u = TransformToSpatialDomainWithF(u_bar, obj.Nx, obj.Ny, obj.Nz);
            v = TransformToSpatialDomainWithF(v_bar, obj.Nx, obj.Ny, obj.Nz);
        end
        
        function [w,zeta] = VerticalFieldsAtTime(obj, t)
            phase_plus = MakeHermitian( exp(sqrt(-1)*obj.Omega*t) );
            phase_minus = MakeHermitian( exp(-sqrt(-1)*obj.Omega*t) );
            w_bar = obj.w_plus.*phase_plus + obj.w_minus.*phase_minus;
            zeta_bar = obj.zeta_plus.*phase_plus + obj.zeta_minus.*phase_minus;
            
            % Re-order to convert to an fast cosine transform
            w = TransformToSpatialDomainWithG(w_bar, obj.Nx, obj.Ny, obj.Nz);
            zeta = TransformToSpatialDomainWithG(zeta_bar, obj.Nx, obj.Ny, obj.Nz);
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forces a 3D matrix to be Hermitian, ready for an FFT (public)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = MakeHermitian(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

for k=1:K
   for i=M:-1:1
       for j=N:-1:1
           ii = mod(M-i+1, M) + 1;
           jj = mod(N-j+1, N) + 1;
           if i == ii && j == jj
               A(i,j,k) = real(A(i,j,k)); % self-conjugate term
           else
               A(i,j,k) = conj(A(ii,jj,k));
           end
       end
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forces the *non-inertial* terms of a 3D matrix to be Hermitian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = MakeHermitianNoInertial(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

for k=1:K
   for i=M:-1:1
       for j=N:-1:1
           ii = mod(M-i+1, M) + 1;
           jj = mod(N-j+1, N) + 1;
           if i == 1 && j == 1
               continue;
           elseif i == ii && j == jj
               A(i,j,k) = real(A(i,j,k)); % self-conjugate term
           else
               A(i,j,k) = conj(A(ii,jj,k));
           end
       end
   end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forces the *inertial* terms of a 3D matrix to be Hermitian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = MakeHermitianOnlyInertial(A)
A(1,1,:) = real(A(1,1,:));
end



function A = GenerateHermitianRandomMatrix( size )

nX = size(1); nY = size(2); nZ = size(3);
A = MakeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
A(1,1,:) = 2*real(A(1,1,:)); % Double the zero frequency
A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
A(:,:,nZ) = zeros(nX,nY); % Because we can't resolve the last mode.

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the phase information given the amplitudes (public)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = TransformToSpatialDomainWithF( u_bar, Nx, Ny, Nz )
% Here we use what I call the 'Fourier series' definition of the ifft, so
% that the coefficients in frequency space have the same units in time.
    u = Nx*Ny*ifft(ifft(u_bar,Nx,1,'symmetric'),Ny,2,'symmetric');
    u = fft(cat(3, zeros(Nx,Ny), 0.5*u(:,:,1:Nz-1), u(:,:,Nz), 0.5*u(:,:,Nz-1:-1:1)),2*Nz,3);
    u = u(:,:,1:Nz);
end

function w = TransformToSpatialDomainWithG( w_bar, Nx, Ny, Nz )
% Here we use what I call the 'Fourier series' definition of the ifft, so
% that the coefficients in frequency space have the same units in time.
    w = Nx*Ny*ifft(ifft(w_bar,Nx,1,'symmetric'),Ny,2,'symmetric');
    w = fft( sqrt(-1)*cat(3, zeros(Nx,Ny), 0.5*w(:,:,1:Nz-1), w(:,:,Nz), -0.5*w(:,:,Nz-1:-1:1)),2*Nz,3);
    w = w(:,:,1:Nz);
end