 file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';

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
f0 = wavemodel.f0;
N0 = wavemodel.N0;
P0_factor = sqrt(g*D*(g*h.*K2 + f0*f0)./(4*f0*f0*h)); % m^(1/2)/s



nFiles = WM.NumberOfTimeSteps;
fileIncrements = 1:100:nFiles;

fprintf('Analyzing %d files\n',length(fileIncrements))


nJ = size(wavemodel.Kh,3);
P_p = zeros([length(fileIncrements) nK nJ]);
P_m = zeros([length(fileIncrements) nK nJ]);
P_0 = zeros([length(fileIncrements) nK nJ]);

for iTime = 1:1
    iFile = fileIncrements(iTime);
    [t,u,v,w,rho_prime] = WM.VariableFieldsAtTimeIndex(iFile,'t','u','v','w','rho_prime');
%         [u,v,w,rho_prime] = WM.VariableFieldsAtTimeIndex(0,'u','v','w','rho_prime'); t = 0;

    wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);
    
    P0 = wavemodel.B .* P0_factor;
    for iJ = 1:nJ
        Amp_plus = wavemodel.Amp_plus(:,:,iJ);
        Amp_minus = wavemodel.Amp_minus(:,:,iJ);
        Amp_zero = P0(:,:,iJ);
        for iK = 1:nK
            P_p(iTime,iK,iJ) = sum(abs(Amp_plus(indicesForK{iK})).^2);
            P_m(iTime,iK,iJ) = sum(abs(Amp_minus(indicesForK{iK})).^2);
            P_0(iTime,iK,iJ) = sum(abs(Amp_zero(indicesForK{iK})).^2);
        end
    end  
    
    eta = g * rho_prime /(wavemodel.rho0 * N0 * N0);
    integratedEnergy = trapz(wavemodel.z,mean(mean( u.^2 + v.^2 + w.^2 + N0*N0*eta.*eta, 1 ), 2 ) )/2;
    fprintf('total integrated energy: %f m^3/s^2\n', integratedEnergy);
    
    spectralEnergy = sum(sum(sum(wavemodel.Amp_plus.*conj(wavemodel.Amp_plus) + wavemodel.Amp_minus.*conj(wavemodel.Amp_minus) + P0.*conj(P0) )));
    fprintf('total spectral energy: %f m^3/s^2\n', spectralEnergy);
end

iTime = 1;
figure
plot(k,squeeze(sum(P_p(iTime,:,:),3))), ylog, hold on
plot(k,squeeze(sum(P_m(iTime,:,:),3)))
plot(k,squeeze(sum(P_0(iTime,:,:),3)))