%% Define the beamline
DefineCLARABeamline
clearvars -except beamline driftlist quadlist quaderrlist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
% DefineCLARABeamline
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);
rng(1)

%% Setup
% Load the functions which will be used in this script.
% Access the functions using the object lt
lt = linTools;
lt.beamline = beamline;
lt.ML_Algorithm = 'cwls';

%% Variables

% Fixed variables
numSeeds   = 100;
dcorrStrength   = 1e-5; % Tesla

% Varying
ML_numobs       = linspace(1000,1000,1);
BPMnoise        = linspace(10,70,7)*1e-6;   % unit
FocusingError   = linspace(1,10,10)*1e-3;   % unit

AlignmentError  = linspace(1,90,11)*1e-6; % m 

%%

% fixedFocusingError = 0.005; % percentage
% 
RespMatC = lt.calcRespMatC(dcorrStrength);
% 
% figure(1)
% subplot(3,1,1)
% hold off
% plot(1e3*lt.track_getBPMreadings(), '-.ok')
% ylabel('BPM readings (mm)')
% title('No errors')
% 
% lt.setQuadFerrors(true, fixedFocusingError);
% 
% subplot(3,1,2)
% hold off
% plot(1e3*lt.track_getBPMreadings(), '-.ok')
% ylabel('BPM readings (mm)')
% title('Focusing error')
% 
% lt.setQuadAerrors(true, kickStrength);
% 
% subplot(3,1,3)
% hold off
% plot(1e3*lt.track_getBPMreadings(), '-.ok')
% ylabel('BPM readings (mm)')
% title('Focusing then alignment error')
% xlabel('BPM index')
% 
% lt.setQuadFerrors(false, fixedFocusingError);
% lt.setQuadAerrors(false, kickStrength);
% 
% figure(2)
% subplot(3,1,1)
% hold off
% plot(1e3*lt.track_getBPMreadings(), '-.ok')
% ylabel('BPM readings (mm)')
% title('No errors')
% 
% lt.setQuadAerrors(true, kickStrength);
% 
% subplot(3,1,2)
% hold off
% plot(1e3*lt.track_getBPMreadings(), '-.ok')
% ylabel('BPM readings (mm)')
% title('Alignment error')
% 
% lt.setQuadFerrors(true, fixedFocusingError);
% 
% subplot(3,1,3)
% hold off
% plot(1e3*lt.track_getBPMreadings(), '-.ok')
% ylabel('BPM readings (mm)')
% xlabel('BPM index')
% title('Alignment then focusing error')

rmsBPMvals        = zeros(numSeeds,1);
mean_rmsBPMvals   = zeros(numel(AlignmentError),1);
err_rmsBPMvals    = zeros(numel(AlignmentError),2);

rmsBPMvals_C      = zeros(numSeeds,1);
mean_rmsBPMvals_C = zeros(numel(AlignmentError),1);
err_rmsBPMvals_C  = zeros(numel(AlignmentError),2);

FocusingError = 0.01;

for i=1:numel(AlignmentError)
   for n=1:numSeeds
        rng(n)
        lt.setQuadFerrors(true, FocusingError);
        lt.setQuadAerrors(true, AlignmentError(i));
        
        bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*70e-6;
        
        rmsBPMvals(n)   = rms(bpmValsY);
        rmsBPMvals_C(n) = rms(lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), 70e-6));
        
        lt.setQuadAerrors(false, AlignmentError(i));
        lt.setQuadFerrors(false, FocusingError);
   end
   mean_rmsBPMvals(i)     = mean(rmsBPMvals(:));
   err_rmsBPMvals(i,:)    = prctile(rmsBPMvals,[5 95]);
   
   mean_rmsBPMvals_C(i)   = mean(rmsBPMvals_C(:));
   err_rmsBPMvals_C(i,:)  = prctile(rmsBPMvals_C,[5 95]);
end

figure(2)
hold off
errorbar(AlignmentError*1e6,mean_rmsBPMvals(:)*1e6,err_rmsBPMvals(:,1)*1e6,err_rmsBPMvals(:,2)*1e6, '-ok')
hold on
errorbar(AlignmentError*1e6,mean_rmsBPMvals_C(:)*1e6,err_rmsBPMvals_C(:,1)*1e6,err_rmsBPMvals_C(:,2)*1e6, '-.xr')

ylabel('RMS vertical trajectory [\mum]')
xlabel('Quadrupole equivalent displacement [\mum]')






