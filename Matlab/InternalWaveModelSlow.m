classdef InternalWaveModelSlow < InternalWaveModel
    
    properties
       S
       Sprime
       nModes
    end
    
    methods
        function self = InternalWaveModelSlow(dims, n, latitude, N0)
            self@InternalWaveModel(dims,n,latitude,N0);
            
            
        end
        
        function InitializeWaveProperties(self)
            InitializeWaveProperties@InternalWaveModel(self)
            
            self.nModes = max(self.j);
            
            % We no longer want to assume the normalization
            self.F = ones(size(self.F));
            self.G = ones(size(self.G));
            
            self.S = zeros(self.Nz, self.nModes, self.Nx, self.Ny);
            self.Sprime = zeros(self.Nz, self.nModes, self.Nx, self.Ny);

            g = 9.81;
            
            rho0 = 1025;
            rho = @(z) -(self.N0*self.N0*rho0/g)*z + rho0;
            im = InternalModes(rho,[-self.Lz 0],self.z,self.latitude,'nEVP', 64);
            im.nModes = self.nModes;
            
            K2 = self.K2(:,:,1);
            [K2_unique,~,iK2_unique] = unique(K2);
            fprintf('Solving the EVP for %d unique wavenumbers.\n',length(K2_unique));
            for iUnique=1:length(K2_unique)
                kk = K2_unique(iUnique);  
                
%                 [F,G,h] = self.ModesForConstantStratificationAtWavenumber(sqrt(kk));
                [F,G,h] = im.ModesAtWavenumber(sqrt(kk));
                h = reshape(h,[1 1 self.nModes]);
                
                indices = find(iK2_unique==iUnique);
                for iIndex=1:length(indices)
                    currentIndex = indices(iIndex);
                    [i,j] = ind2sub([self.Nx self.Ny], currentIndex);
                    self.h(i,j,:) = h;
                    self.S(:,:,i,j) = G;
                    self.Sprime(:,:,i,j) = F;
                end
            end

            % Slower algoritm
%             for i=1:self.Nx
%                 for j=1:self.Ny
%                     kk = self.K2(i,j,1);
%                                          
%                     [F,G,h] = im.ModesAtWavenumber(sqrt(kk));
%                     self.h(i,j,:) = reshape(h,[1 1 self.nModes]);
%                     self.S(:,:,i,j) = G;
%                     self.Sprime(:,:,i,j) = F;
%                 end
%             end
            
            % Need to set omega

        end
        
        function [F,G,h] = ModesForConstantStratificationAtWavenumber(self, k)
            kk = k*k;
            z = reshape(self.z,[self.Nz 1 1 1]);
            mode = reshape(1:self.nModes,[1 self.nModes 1 1]);
            kz = (pi/self.Lz)*mode;
            g = 9.81;
            
            h = (self.N0*self.N0-self.f0*self.f0)./(g*(kz.*kz + kk));
            G = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)))*sin(kz.*z);
            F = sqrt(2*g/(self.Lz*(self.N0*self.N0-self.f0*self.f0)))* h.*kz.*cos(kz.*z);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Computes the phase information given the amplitudes (internal)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = TransformToSpatialDomainWithF(self, u_bar)
            u_bar = permute(u_bar,[3 1 2]); % Speed optimization: keep matrices adjacent in memory
            
            u_temp = zeros(size(u_bar));
            for i=1:self.Nx
                for j=1:self.Ny
                    u_temp(i,j,:) = self.Sprime(:,:,i,j)*u_bar(:,i,j);
                end
            end
            
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            u = self.Nx*self.Ny*ifft(ifft(u_temp,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function w = TransformToSpatialDomainWithG(self, w_bar )
            w_bar = permute(w_bar,[3 1 2]); % Speed optimization: keep matrices adjacent in memory
            
            w_temp = zeros(size(w_bar));
            for i=1:self.Nx
                for j=1:self.Ny
                    w_temp(i,j,:) = self.S(:,:,i,j)*w_bar(:,i,j);
                end
            end
            
            % Here we use what I call the 'Fourier series' definition of the ifft, so
            % that the coefficients in frequency space have the same units in time.
            w = self.Nx*self.Ny*ifft(ifft(w_temp,self.Nx,1),self.Ny,2,'symmetric');
        end
    end
end