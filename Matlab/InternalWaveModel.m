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
% November 17th, 2016   Version 1.2
% December 9th, 2016    Version 1.3

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
        
        K2, Kh, F, G, M, h, Omega, Omega_plus, Omega_minus, f0, C
        u_plus, u_minus, v_plus, v_minus, w_plus, w_minus, zeta_plus, zeta_minus
        period
        version = 1.3
        dctScratch, dstScratch;
        performSanityChecks = 0
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
            
            % Preallocate this array for a faster dct
            obj.dctScratch = zeros(obj.Nx,obj.Ny,2*obj.Nz);
            obj.dstScratch = complex(zeros(obj.Nx,obj.Ny,2*obj.Nz));
            
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
            obj.C = sqrt( C2 );                         % Mode speed
            obj.h = C2/g;                               % Mode depth
            obj.Omega = sqrt(obj.C.*obj.C.*obj.K2 + obj.f0*obj.f0);         % Mode frequency
            
            % F contains the coefficients for the U-V modes
            obj.F = (obj.h.*obj.M)*sqrt(2*g/(obj.Lz*(obj.N0*obj.N0-obj.f0*obj.f0)));
            
            % G contains the coefficients for the W-modes
            obj.G = sqrt(2*g/(obj.Lz*(obj.N0*obj.N0-obj.f0*obj.f0)));
            
            % Create the hermitian conjugates of the phase vectors;
            obj.Omega(:,(obj.Ny/2+1):end,:) = -obj.Omega(:,(obj.Ny/2+1):end,:);
            obj.Omega((obj.Nx/2+1):end,1,:) = -obj.Omega((obj.Nx/2+1):end,1,:);      
        end
        
        function ShowDiagnostics(obj)
            omega = abs(obj.Omega);
            fprintf('Model resolution is %.2f x %.2f x %.2f meters.\n', obj.x(2)-obj.x(1), obj.y(2)-obj.y(1), obj.z(2)-obj.z(1));
            fprintf('The ratio N0/f0 is %.1f.\n', obj.N0/obj.f0);
            fprintf('Discretization effects will become apparent after %.1f hours in the frequency domain as the fastest modes traverse the domain.\n', max([obj.Lx obj.Ly])/max(max(max(obj.C)))/3600);
            sortedOmega = sort(unique(reshape(omega(:,:,1),1,[])));
            fprintf('j=1 mode has discrete frequencies (%.4f f0, %.4f f0, ..., %.4f N0, %.4f N0)\n', sortedOmega(1)/obj.f0, sortedOmega(2)/obj.f0, sortedOmega(end-1)/obj.N0, sortedOmega(end)/obj.N0);
            sortedOmega = sort(unique(reshape(omega(:,:,end),1,[])));
            fprintf('j=%d mode has discrete frequencies (%.4f f0, %.4f f0, ..., %.4f N0, %.4f N0)\n', obj.Nz, sortedOmega(1)/obj.f0, sortedOmega(2)/obj.f0, sortedOmega(end-1)/obj.N0, sortedOmega(end)/obj.N0);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithPlaneWave(obj, k0, l0, j0, UAmp, sign)
            % User input sanity checks. We don't deal with the Nyquist.
            if (k0 <= -obj.Nx/2 || k0 >= obj.Nx/2)
                error('Invalid choice for k0 (%d). Must be an integer %d < k0 < %d',k0,-obj.Nx/2+1,obj.Nx/2-1);
            end
            if (l0 <= -obj.Ny/2 || l0 >= obj.Ny/2)
                error('Invalid choice for l0 (%d). Must be an integer %d < l0 < %d',l0,-obj.Ny/2+1,obj.Ny/2+1);
            end
            if (j0 < 1 || j0 >= obj.Nz)
                error('Invalid choice for j0 (%d). Must be an integer 0 < j < %d',j0, obj.Nz);
            end
            
            % Deal with the negative wavenumber cases (and inertial)
            if l0 == 0 && k0 == 0 % inertial
                sign=1;
            elseif l0 == 0 && k0 < 0
                k0 = -k0;
                sign = -1*sign;
                UAmp = -1*UAmp;
            elseif l0 < 0
                l0 = -l0;
                k0 = -k0;
                sign = -1*sign;
                UAmp = -1*UAmp;
            end
            
            % Rewrap (k0,l0) to follow standard FFT wrapping. l0 should
            % already be correct.
            if (k0 < 0)
                k0 = obj.Nx + k0;
            end
                
            myH = obj.h(k0+1,l0+1,j0);
            m = j0*pi/obj.Lz;
            g = 9.81;
            F_coefficient = myH * m * sqrt(2*g/obj.Lz)/sqrt(obj.N0^2 - obj.f0^2);
            
            U = zeros(size(obj.K));
            U(k0+1,l0+1,j0) = UAmp*sqrt(myH)/F_coefficient/2;
            if sign > 0
                A_plus = MakeHermitian(U);
                A_minus = zeros(size(U));
                A_minus(1,1,:) = A_plus(1,1,:); % Inertial oscillations are created using this trick.                
            else
                A_plus = zeros(size(U));
                A_minus = MakeHermitian(U);
            end
            
            obj.GenerateWavePhases(A_plus,A_minus);
            
            omega = obj.Omega(k0+1,l0+1,j0);
            obj.period = 2*pi/omega;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a full Garrett-Munk spectrum (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function InitializeWithGMSpectrum(obj, Amp)
            % GM Parameters
            j_star = 3;
            L_gm = 1.3e3; % thermocline exponential scale, meters
            invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
            E_gm = 6.3e-5; % non-dimensional energy parameter
            E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*Amp;
                      
            % Compute the proper vertical function normalization
            H = (j_star+(1:1024)).^(-5/2);
            H_norm = sum(H);
            
            % Do the same for the frequency function.
            B_norm = 1/atan(sqrt(obj.N0*obj.N0/(obj.f0*obj.f0)-1));
            
            % This function tells you how much energy you need between two
            % frequencies for a given vertical mode.
            GM2D_int = @(omega0,omega1,j) (E/H_norm)*((j+j_star).^(-5/2))*(atan(obj.f0/sqrt(omega0*omega0-obj.f0*obj.f0)) - atan(obj.f0/sqrt(omega1*omega1-obj.f0*obj.f0)))*B_norm;
            
            GM2D_uv_int = @(omega0,omega1,j) B_norm*(E/H_norm)*((j+j_star).^(-5/2))*( obj.f0*sqrt(omega1*omega1-obj.f0*obj.f0)/(2*omega1*omega1) - (3/2)*atan(obj.f0/sqrt(omega1*omega1-obj.f0*obj.f0)) - obj.f0*sqrt(omega0*omega0-obj.f0*obj.f0)/(2*omega0*omega0) + (3/2)*atan(obj.f0/sqrt(omega0*omega0-obj.f0*obj.f0)));
            GM2D_w_int = @(omega0,omega1,j) B_norm*(E/H_norm)*((j+j_star).^(-5/2))*( obj.f0*sqrt(omega1*omega1-obj.f0*obj.f0) - obj.f0*obj.f0*atan(obj.f0/sqrt(omega1*omega1-obj.f0*obj.f0)) - obj.f0*sqrt(omega0*omega0-obj.f0*obj.f0) + obj.f0*obj.f0*atan(obj.f0/sqrt(omega0*omega0-obj.f0*obj.f0)));
            GM2D_zeta_int = @(omega0,omega1,j) B_norm*(E/H_norm)*((j+j_star).^(-5/2))*( ((omega1*omega1-obj.f0*obj.f0)^(3/2))/(2*obj.f0*omega1*omega1) - (1/2)*atan(obj.f0/sqrt(omega1*omega1-obj.f0*obj.f0)) - sqrt(omega1*omega1-obj.f0*obj.f0)/(2*obj.f0) - ((omega0*omega0-obj.f0*obj.f0)^(3/2))/(2*obj.f0*omega0*omega0) + (1/2)*atan(obj.f0/sqrt(omega0*omega0-obj.f0*obj.f0)) + sqrt(omega0*omega0-obj.f0*obj.f0)/(2*obj.f0) );
            
            % Do a quick check to see how much energy is lost due to
            % limited vertical resolution.
            totalEnergy = 0;
            for mode=1:(max(obj.j)/1)
                totalEnergy = totalEnergy + GM2D_int(obj.f0,obj.N0,mode);
            end
            fprintf('You are missing %.2f%% of the energy due to limited vertical modes.\n',100-100*totalEnergy/E);
            
            % Find the *second* lowest frequency
            [sortedOmegas, indices] = sort(reshape(abs(obj.Omega(:,:,max(obj.j)/2)),1,[]));
            omegaStar = sortedOmegas(2);
            
            wVariancePerMode = [];
            for mode=1:(max(obj.j)/2)
                wVariancePerMode(mode) = GM2D_w_int(obj.f0+(min(min(obj.Omega(2:end,2:end,mode)))-obj.f0)/2,max(max(obj.Omega(:,:,mode))),1);
                wVariancePerModeStar(mode) = GM2D_w_int(obj.f0+(omegaStar-obj.f0)/2,max(max(obj.Omega(:,:,mode))),1);
            end
            
            shouldUseOmegaStar = 0;
            
            % Sort the frequencies (for each mode) and distribute energy.
            GM3D = zeros(size(obj.Kh));
            for iMode = 1:(max(obj.j)/1)
                % Stride to the linear index for the full 3D matrix
                modeStride = (iMode-1)*size(obj.Omega,1)*size(obj.Omega,2);
                
                % Sort the linearized frequencies for this mode.
                [sortedOmegas, indices] = sort(reshape(abs(obj.Omega(:,:,iMode)),1,[]));
                
                % Then find where the omegas differ.
                omegaDiffIndices = find(diff(sortedOmegas) > 0);
            
                lastIdx = 1;
                omega0 = sortedOmegas(lastIdx);
                leftDeltaOmega = 0;
                for idx=omegaDiffIndices
                    currentIdx = idx+1;
                    nOmegas = currentIdx-lastIdx;
                    
%                     if omega0 ~= obj.f0
%                         continue;
%                     end
                    
                    if shouldUseOmegaStar && omega0 == obj.f0
                        omega1 = omegaStar;
                    else
                        omega1 = sortedOmegas(idx + 1);
                    end
                    rightDeltaOmega = (omega1-omega0)/2;
%                      energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
%                     energyPerFrequency = (GM2D_uv_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(omega0*omega0/(omega0*omega0 + obj.f0*obj.f0));
                    
%                     if omega0 == obj.f0
%                         energyPerFrequency = GM2D_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas;
%                     else 
%                         energyPerFrequency = (GM2D_zeta_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(omega0*omega0/(omega0*omega0 - obj.f0*obj.f0));
%                         energyPerFrequency = (GM2D_w_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(1/(omega0*omega0 - obj.f0*obj.f0));
%                     end
                    
                    energyPerFrequency = (GM2D_uv_int(omega0-leftDeltaOmega,omega0+rightDeltaOmega,iMode)/nOmegas)*(omega0*omega0/(omega0*omega0 + obj.f0*obj.f0));
                    
                    GM3D(indices(lastIdx:(currentIdx-1))+modeStride) = energyPerFrequency;
                    
                    if shouldUseOmegaStar && omega0 == obj.f0
                        leftDeltaOmega = sortedOmegas(idx + 1) - (omega0+rightDeltaOmega);
                        omega0 = sortedOmegas(idx + 1);
                    else
                        omega0 = omega1;
                        leftDeltaOmega = rightDeltaOmega;
                    end
                    
%                     omega0 = sortedOmegas(idx + 1); % same as omega0 = omega1, except when omega0 == f0
                     lastIdx = currentIdx;
%                     leftDeltaOmega = rightDeltaOmega;
                end
                % Still have to deal with the last point.
            end
            fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*sum(sum(sum(GM3D)))/E);
            fprintf('Due to restricted domain size, the j=1,k=l=0 mode contains %.2f%% the total energy.\n',100*GM3D(1,1,1)/sum(sum(sum(GM3D))));
            
            A = sqrt(GM3D/2); % Now split this into even and odd.
            
            % Randomize phases, but keep unit length
            A_plus = GenerateHermitianRandomMatrix( size(obj.K) );
            A_minus = GenerateHermitianRandomMatrix( size(obj.K) );
            
            goodIndices = abs(A_plus) > 0;
            A_plus(goodIndices) = A_plus(goodIndices)./abs(A_plus(goodIndices));
            A_plus = A.*A_plus;
            goodIndices = abs(A_minus) > 0;
            A_minus(goodIndices) = A_minus(goodIndices)./abs(A_minus(goodIndices));
            A_minus = A.*A_minus;
            
%              A_plus = A;
%              A_minus = A;        
            
%             A_plus = A.*GenerateHermitianRandomMatrix( size(obj.K) );
%             A_minus = A.*GenerateHermitianRandomMatrix( size(obj.K) );
            A_minus(1,1,:) = conj(A_plus(1,1,:)); % Intertial motions go only one direction!
            
            GM_sum = sum(sum(sum(GM3D)))/E;
            GM_random_sum = sum(sum(sum(A_plus.*conj(A_plus) + A_minus.*conj(A_minus)  )))/E;
            fprintf('The coefficients sum to %.2f%% GM given the scales, and the randomized field sums to %.2f%% GM\n', 100*GM_sum, 100*GM_random_sum);
            
            obj.GenerateWavePhases(A_plus,A_minus);
        end
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function GenerateWavePhases(obj, U_plus, U_minus)
            alpha = atan2(obj.L,obj.K);
            omega = abs(obj.Omega); % The following definitions assume omega > 0.
            denominator = omega.*sqrt(obj.h);
            
            % Without calling MakeHermitian, this doesn't deal with l=0.
            obj.u_plus = U_plus .* MakeHermitian( ( -sqrt(-1)*obj.f0 .* sin(alpha) + omega .* cos(alpha) )./denominator ) .* obj.F;
            obj.u_minus = U_minus .* MakeHermitian( (sqrt(-1)*obj.f0 .* sin(alpha) + omega .* cos(alpha) )./denominator ) .* obj.F;
            
            obj.v_plus = U_plus .* MakeHermitian( ( sqrt(-1)*obj.f0 .* cos(alpha) + omega .* sin(alpha) )./denominator ) .* obj.F;
            obj.v_minus = U_minus .* MakeHermitian( ( -sqrt(-1)*obj.f0 .* cos(alpha) + omega .* sin(alpha) )./denominator ) .* obj.F;
            
            obj.w_plus = U_plus .* MakeHermitian(-sqrt(-1) *  obj.Kh .* sqrt(obj.h) ) * obj.G;
            obj.w_minus = U_minus .* MakeHermitian( -sqrt(-1) * obj.Kh .* sqrt(obj.h) ) * obj.G;
            
            obj.zeta_plus = U_plus .* MakeHermitian( -obj.Kh .* sqrt(obj.h) ./ omega ) * obj.G;
            obj.zeta_minus = U_minus .* MakeHermitian( obj.Kh .* sqrt(obj.h) ./ omega ) * obj.G;  
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Compute the dynamical fields at a given time (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [u,v] = VelocityFieldAtTime(obj, t)
            phase_plus = exp(sqrt(-1)*obj.Omega*t);
            phase_minus = exp(-sqrt(-1)*obj.Omega*t);
            u_bar = obj.u_plus.*phase_plus + obj.u_minus.*phase_minus;
            v_bar = obj.v_plus.*phase_plus + obj.v_minus.*phase_minus;
            
            if obj.performSanityChecks == 1
                CheckHermitian(u_bar);CheckHermitian(v_bar);
            end
            
            u = obj.TransformToSpatialDomainWithF(u_bar, obj.Nx, obj.Ny, obj.Nz);
            v = obj.TransformToSpatialDomainWithF(v_bar, obj.Nx, obj.Ny, obj.Nz);
        end
        
        function [w,zeta] = VerticalFieldsAtTime(obj, t)
            phase_plus = exp(sqrt(-1)*obj.Omega*t);
            phase_minus = exp(-sqrt(-1)*obj.Omega*t);
            w_bar = obj.w_plus.*phase_plus + obj.w_minus.*phase_minus;
            zeta_bar = obj.zeta_plus.*phase_plus + obj.zeta_minus.*phase_minus;
            
            if obj.performSanityChecks == 1
                CheckHermitian(w_bar);CheckHermitian(zeta_bar);
            end
            
            w = obj.TransformToSpatialDomainWithG(w_bar, obj.Nx, obj.Ny, obj.Nz);
            zeta = obj.TransformToSpatialDomainWithG(zeta_bar, obj.Nx, obj.Ny, obj.Nz);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = TransformToSpatialDomainWithF(obj, u_bar, Nx, Ny, Nz )
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            u = Nx*Ny*ifft(ifft(u_bar,Nx,1),Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            obj.dctScratch = cat(3, zeros(Nx,Ny), 0.5*u(:,:,1:Nz-1), u(:,:,Nz), 0.5*u(:,:,Nz-1:-1:1));
  
            u = fft(obj.dctScratch,2*Nz,3);
            if obj.performSanityChecks == 1
                ratio = max(max(max(abs(imag(u)))))/max(max(max(abs(real(u)))));
                if ratio > 1e-6
                    fprintf('WARNING: The inverse cosine transform reports an unreasonably large imaginary part, %.2g.\n',ratio);
                end
            end
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            u = u(:,:,1:Nz);
        end
        
        function w = TransformToSpatialDomainWithG(obj, w_bar, Nx, Ny, Nz )
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            w = Nx*Ny*ifft(ifft(w_bar,Nx,1),Ny,2,'symmetric');
            
            % Re-order to convert to an fast cosine transform
            obj.dstScratch = sqrt(-1)*cat(3, zeros(Nx,Ny), 0.5*w(:,:,1:Nz-1), w(:,:,Nz), -0.5*w(:,:,Nz-1:-1:1));
            
            w = fft( obj.dstScratch,2*Nz,3);
            if obj.performSanityChecks == 1
                ratio = max(max(max(abs(imag(w)))))/max(max(max(abs(real(w)))));
                if ratio > 1e-6
                    fprintf('WARNING: The inverse sine transform reports an unreasonably large imaginary part, %.2g.\n',ratio);
                end
            end
            % should not have to call real, but for some reason, with enough
            % points, it starts generating some small imaginary component.
            w = w(:,:,1:Nz);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Forces a 3D matrix to be Hermitian, ready for an FFT (internal)
%
% The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
% numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
% conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
% numbers are set to zero to avoid complications.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = MakeHermitian(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

% The order of the for-loop is chosen carefully here.
for k=1:K
    for j=1:(N/2+1)
        for i=1:M
            ii = mod(M-i+1, M) + 1;
            jj = mod(N-j+1, N) + 1;
            if i == ii && j == jj
                % A(i,j,k) = real(A(i,j,k)); % self-conjugate term
                % This is normally what you'd do, but we're being tricky
                if i == 1 && j == 1
                    continue;
                else
                    A(i,j,k) = 0;
                end
            elseif j == N/2+1 % Kill the Nyquist, rather than fix it.
                A(i,j,k) = 0;
            else % we are letting l=0, k=Nx/2+1 terms set themselves again, but that's okay 
                A(ii,jj,k) = conj(A(i,j,k));
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Checks that the matrix is Hermitian.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = CheckHermitian(A)
M = size(A,1);
N = size(A,2);
K = size(A,3);

for k=1:K
   for i=M:-1:1
       for j=N:-1:1
           ii = mod(M-i+1, M) + 1;
           jj = mod(N-j+1, N) + 1;
           if A(i,j,k) ~= conj(A(ii,jj,k))
               error('Not hermitian conjugate')
           end
       end
   end
end

end

function A = GenerateHermitianRandomMatrix( size )

nX = size(1); nY = size(2); nZ = size(3);
A = MakeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
% A(1,1,:) = 2*real(A(1,1,:)); % Double the zero frequency
A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
A(:,:,nZ) = zeros(nX,nY); % Because we can't resolve the last mode.

end


