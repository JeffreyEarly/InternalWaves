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
            obj.Kh(1,1,:) = 1;      % prevent divide by zero
            
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
            U = zeros(size(obj.K));
            U(k0+1,l0+1,j0) = UAmp;
            if sign > 0
                A_plus = -MakeHermitian(U); % Careful with this, it doubles energy for some k,l
                A_minus = zeros(size(U));
            else
                A_plus = zeros(size(U)); % Careful with this, it doubles energy for some k,l
                A_minus = MakeHermitian(U);
            end
            
            obj.GenerateWavePhases(A_plus,A_minus);
            
            [u,~] = obj.VelocityFieldAtTime(0);
            max_u = max(max(max( u )));
            obj.u_plus = obj.u_plus*(UAmp/max_u);
            obj.u_minus = obj.u_minus*(UAmp/max_u);
            obj.v_plus = obj.v_plus*(UAmp/max_u);
            obj.v_minus = obj.v_minus*(UAmp/max_u);
            obj.w_plus = obj.w_plus*(UAmp/max_u);
            obj.w_minus = obj.w_minus*(UAmp/max_u);
            obj.zeta_plus = obj.zeta_plus*(UAmp/max_u);
            obj.zeta_minus = obj.zeta_minus*(UAmp/max_u);
            
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
            GM2D_function = @(k,j) E*( (2/H_norm)*(j_star.^(3/2))./(j+j_star).^(5/2) ) .* (2/pi)*obj.f0*(j*pi/D).*(j*pi/D).*sqrt( (obj.N0*obj.N0-obj.f0*obj.f0)./(k.*k+(j*pi/D).*(j*pi/D)))./(obj.N0*obj.N0*k.*k+obj.f0*obj.f0*(j*pi/D).*(j*pi/D));
            
            
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
                func = @(k) GM2D_function(k,mode);
                GM3D(1,1,mode) = integral(func,m(1),m(1)+dm/2);
                for i=2:length(m)
                    m_lower = m(i)-dm/2;
                    m_upper = m(i)+dm/2;
                    indices = obj.Kh(:,:,mode) >= m_lower & obj.Kh(:,:,mode) < m_upper;
                    n = sum(sum(sum(indices)));
                    total = integral(func,m_lower,m_upper)/n;
                    fullIndices(:,:,mode) = indices;
                    GM3D(fullIndices) = total;
                end
            end
            
            
            
            A = sqrt(GM3D/2); % Now split this into even and odd.
            
            A_plus = A.*MakeHermitian(randn(size(obj.K)) + sqrt(-1)*randn(size(obj.K)) )/sqrt(2);
            A_minus = A.*MakeHermitian(randn(size(obj.K)) + sqrt(-1)*randn(size(obj.K)) )/sqrt(2);
            
            GM_sum = sum(sum(sum(GM3D)))/E;
            GM_random_sum = sum(sum(sum(A_plus.*conj(A_plus) + A_minus.*conj(A_minus)  )))/E;
            fprintf('The coefficients sum to %.2f GM given the scales, and the randomized field gives %f\n', GM_sum, GM_random_sum);
            
            obj.GenerateWavePhases(A_plus,A_minus);
            
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GenerateWavePhases(obj, U_plus, U_minus)
            
            denominator = obj.Kh.*obj.Omega.*sqrt(obj.h);
            obj.u_plus = U_plus.*MakeHermitian(obj.F.*(sqrt(-1)*obj.L*obj.f0 - obj.K.*obj.Omega)./denominator);
            obj.u_minus = U_minus.*MakeHermitian(obj.F.*(sqrt(-1)*obj.L*obj.f0 + obj.K.*obj.Omega)./denominator);
            
            obj.v_plus = U_plus.*MakeHermitian(obj.F.*(-sqrt(-1)*obj.K*obj.f0 - obj.L.*obj.Omega)./denominator);
            obj.v_minus = U_minus.*MakeHermitian(obj.F.*(-sqrt(-1)*obj.K*obj.f0 + obj.L.*obj.Omega)./denominator);
            
            obj.w_plus = U_plus .* MakeHermitian(sqrt(-1) *  obj.Kh .* sqrt(obj.h) ) * obj.G;
            obj.w_minus = U_minus .* MakeHermitian( -sqrt(-1) * obj.Kh .* sqrt(obj.h) ) * obj.G;
            
            obj.zeta_plus = U_plus .* MakeHermitian( obj.Kh .* sqrt(obj.h) ./ obj.Omega ) * obj.G;
            obj.zeta_minus = U_minus .* MakeHermitian( obj.Kh .* sqrt(obj.h) ./ obj.Omega ) * obj.G;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u,v] = VelocityFieldAtTime(obj, t)
            phase_plus = MakeHermitian( exp(sqrt(-1)*obj.Omega*t) );
            phase_minus = MakeHermitian( exp(-sqrt(-1)*obj.Omega*t) );
            u_bar = obj.u_plus.*phase_plus + obj.u_minus.*phase_minus;
            v_bar = obj.v_plus.*phase_plus + obj.v_minus.*phase_minus;
            
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
% Computes the phase information given the amplitudes (public)
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
%            if k==1 && i == ii && j == jj
%               fprintf('i,j = (%d,%d)\n',i,j); 
%            end
           A(i,j,k) = conj(A(ii,jj,k));
       end
   end
end

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
    u = fft(cat(3, zeros(Nx,Ny), 0.5*u(:,:,1:Nz-2), u(:,:,Nz-1), zeros(Nx,Ny), u(:,:,Nz-1), 0.5*u(:,:,Nz-2:-1:1)),2*Nz,3);
    u = u(:,:,1:Nz);
end

function w = TransformToSpatialDomainWithG( w_bar, Nx, Ny, Nz )
% Here we use what I call the 'Fourier series' definition of the ifft, so
% that the coefficients in frequency space have the same units in time.
    w = Nx*Ny*ifft(ifft(w_bar,Nx,1,'symmetric'),Ny,2,'symmetric');
    w = fft( sqrt(-1)*cat(3, zeros(Nx,Ny), 0.5*w(:,:,1:Nz-2), w(:,:,Nz-1), zeros(Nx,Ny), -w(:,:,Nz-1), -0.5*w(:,:,Nz-2:-1:1)),2*Nz,3);
    w = w(:,:,1:Nz);
end