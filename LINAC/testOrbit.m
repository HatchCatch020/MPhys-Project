%% Define the beamline
DefineCLARABeamline
clearvars -except beamline driftlist quadflist quadflist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
% DefineCLARABeamline
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);
rng(200)

%% Setup
% Load the functions which will be used in this script.
% Access the functions using the object lt
lt = linTools;
lt.beamline = beamline;
lt.ML_Algorithm = 'cwls';

%% Variables

% Fixed variables
numIterations   = 5;
dcorrStrength   = 1e-4; % Tesla
kickStrength    = 100;   % Tesla
kick = randn(numel(beamline.quaderrlist),1)*kickStrength;

% Varying
ML_numobs       = linspace(1000,1000,1);
BPMnoise        = linspace(10,70,7)*1e-6;   % unit
FocusingError   = linspace(1,10,10)*1e-3;   % unit

%%

fixedFocusingError = 0.9;

RespMatC = lt.calcRespMatC(dcorrStrength);

figure(1)
subplot(3,1,1)
hold off
plot(1e3*lt.track_getBPMreadings(), '-.ok')
title('No errors')

lt.setQuadFerrors(true, fixedFocusingError);

subplot(3,1,2)
hold off
plot(1e3*lt.track_getBPMreadings(), '-.ok')
title('Focusing error')

lt.setQuadAerrors(true, kick);

subplot(3,1,3)
hold off
plot(1e3*lt.track_getBPMreadings(), '-.ok')
title('Focusing then alignment error')

lt.setQuadFerrors(false, fixedFocusingError);
lt.setQuadAerrors(false, kick);

figure(2)
subplot(3,1,1)
hold off
plot(1e3*lt.track_getBPMreadings(), '-.ok')
title('No errors')

lt.setQuadAerrors(true, kick);

subplot(3,1,2)
hold off
plot(1e3*lt.track_getBPMreadings(), '-.ok')
title('Alignment error')

lt.setQuadFerrors(true, fixedFocusingError);

subplot(3,1,3)
hold off
plot(1e3*lt.track_getBPMreadings(), '-.ok')
title('Alignment then focusing error')


