% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';
file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';

WM = WintersModel(file);
wavemodel = WM.WaveModelFromInitialConditionsFile;

% Create a reasonable wavenumber axis
allKs = unique(reshape(abs(wavemodel.Kh),[],1),'sorted');
deltaK = max(diff(allKs));
k = 0:deltaK:max(allKs);

% Create a map from the line wavenumber index to K
nK = length(k)-1;
indicesForK = cell(nK,1);
for i=1:nK
   indicesForK(i) = {find( k(i) <= squeeze(wavemodel.Kh(:,:,1)) & squeeze(wavemodel.Kh(:,:,1)) < k(i+1))};
end
k(end) = [];

% Conversion factor to go from B to P0 (depth integrated energy)
g = wavemodel.g;
K2 = wavemodel.K2;
h = wavemodel.h;
D = wavemodel.Lz;
f = wavemodel.f0;
N = wavemodel.N0;
omega = wavemodel.Omega;

% P0_TE_factor = sqrt( (g/f/f) * omega.*omega .* (N*N - omega.*omega - f*f) / (N*N - f*f) ); % m^(1/2)/s

Ppm_HKE_factor = (1 + f*f./(omega.*omega)) .* (N*N - omega.*omega) / (2 * (N*N - f*f) );
Ppm_VKE_factor = (omega.*omega - f*f) / (2 * (N*N - f*f) );
Ppm_PE_factor = N*N .* (omega.*omega - f*f) ./ (2 * (N*N - f*f) * omega.*omega );
P0_HKE_factor = (g/(f*f)) * (omega.*omega - f*f) .* (N*N - omega.*omega) / (2 * (N*N - f*f) );
P0_PE_factor = g*N*N/(N*N-f*f)/2;


nFiles = WM.NumberOfTimeSteps;
fileIncrements = 1:100:nFiles;
fileIncrements = 1:1;

fprintf('Analyzing %d files\n',length(fileIncrements))


nJ = size(wavemodel.Kh,3);
P_p_KE = zeros([length(fileIncrements) nK nJ]);
P_m_KE = zeros([length(fileIncrements) nK nJ]);
P_0_KE = zeros([length(fileIncrements) nK nJ]);

P_p_PE = zeros([length(fileIncrements) nK nJ]);
P_m_PE = zeros([length(fileIncrements) nK nJ]);
P_0_PE = zeros([length(fileIncrements) nK nJ]);

for iTime = 1:length(fileIncrements)
    iFile = fileIncrements(iTime);
    [t,u,v,w,rho_prime] = WM.VariableFieldsAtTimeIndex(iFile,'t','u','v','w','rho_prime');
%         [u,v,w,rho_prime] = WM.VariableFieldsAtTimeIndex(0,'u','v','w','rho_prime'); t = 0;

    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
    
    P_p_KE_full = wavemodel.Amp_plus.*conj(wavemodel.Amp_plus).*(Ppm_HKE_factor + Ppm_VKE_factor);
    P_m_KE_full = wavemodel.Amp_minus.*conj(wavemodel.Amp_minus).*(Ppm_HKE_factor + Ppm_VKE_factor);
    P_0_KE_full = wavemodel.B.*conj(wavemodel.B).*P0_HKE_factor;
    
    P_p_PE_full = wavemodel.Amp_plus.*conj(wavemodel.Amp_plus).*Ppm_PE_factor;
    P_m_PE_full = wavemodel.Amp_minus.*conj(wavemodel.Amp_minus).*Ppm_PE_factor;
    P_0_PE_full = wavemodel.B.*conj(wavemodel.B).*P0_PE_factor;
    
    offset = size(P_p_PE_full,1)*size(P_p_PE_full,2);
    
    for iJ = 1:nJ
        for iK = 1:nK
            P_p_KE(iTime,iK,iJ) = sum( P_p_KE_full( (iJ-1)*offset + indicesForK{iK}) );
            P_m_KE(iTime,iK,iJ) = sum( P_m_KE_full( (iJ-1)*offset + indicesForK{iK}) );
            P_0_KE(iTime,iK,iJ) = sum( P_0_KE_full( (iJ-1)*offset + indicesForK{iK}) );
            
            P_p_PE(iTime,iK,iJ) = sum( P_p_PE_full( (iJ-1)*offset + indicesForK{iK}) );
            P_m_PE(iTime,iK,iJ) = sum( P_m_PE_full( (iJ-1)*offset + indicesForK{iK}) );
            P_0_PE(iTime,iK,iJ) = sum( P_0_PE_full( (iJ-1)*offset + indicesForK{iK}) );
        end
    end  
    
    eta = g * rho_prime /(wavemodel.rho0 * N * N);
    integratedEnergy = trapz(wavemodel.z,mean(mean( u.^2 + v.^2 + w.^2 + N*N*eta.*eta, 1 ), 2 ) )/2;
    fprintf('total integrated energy: %f m^3/s^2\n', integratedEnergy);
    
    spectralEnergy = sum(sum(sum( P_p_KE(iTime,:,:) + P_m_KE(iTime,:,:) + P_0_KE(iTime,:,:) + P_p_PE(iTime,:,:) + P_m_PE(iTime,:,:) + P_0_PE(iTime,:,:))));
    fprintf('total spectral energy: %f m^3/s^2\n', spectralEnergy);
end


iTime = 1;
figure

plot(k,squeeze(sum(P_p_KE(iTime,:,:),3))), ylog, hold on
plot(k,squeeze(sum(P_m_KE(iTime,:,:),3)))
plot(k,squeeze(sum(P_0_KE(iTime,:,:),3)))
plot(k,squeeze(sum(P_p_PE(iTime,:,:),3)))
plot(k,squeeze(sum(P_m_PE(iTime,:,:),3)))
plot(k,squeeze(sum(P_0_PE(iTime,:,:),3)))
legend('KE_+','KE_-','KE_0','PE_+','PE_-','PE_0')