% C. Monaghan
% Simple ring example 
% ringCorrectionExample
%% Define the beamline
DefineBeamline
clearvars -except beamline driftlist quadflist quadflist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
% DefineCLARABeamline
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Setup
% Load the functions which will be used in this script.
% Access the functions using the object rt
rt = ringTools;
rt.beamline = beamline;
rt.ML_Algorithm = 'mvn';

%% Variables

numSeeds = 10;

FocusingError   = 0.005;
dcorrStrength   = 1e-6;
AlignmentError  = 20e-6;
ML_numobs       = 100;
BPMnoise        = 30e-6;

%%


%%

bpmValsY            = zeros(numel(beamline.bpmlist),numSeeds);
bpmValsCorrected_C  = zeros(numel(beamline.bpmlist),numSeeds);
bpmValsCorrected_ML = zeros(numel(beamline.bpmlist),numSeeds);

RespMatC = rt.calcRespMatC(dcorrStrength);

for n=1:numSeeds
    rng(n)
    rt.setQuadFerrors(true, FocusingError);
    rt.setQuadAerrors(true, AlignmentError);
    
    RespMatML = rt.calcRespMatML(ML_numobs, BPMnoise, dcorrStrength);
    
    bpmValsY(:,n) = rt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise;

    bpmValsCorrected_C(:,n)  = rt.getBPMvalues_corr(pinv(RespMatC,   1e-3)*(-bpmValsY), BPMnoise);
    bpmValsCorrected_ML(:,n) = rt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), BPMnoise);

    rt.setQuadFerrors(false, FocusingError);
    rt.setQuadAerrors(false, AlignmentError);
end

figure(1)
hold off
plot(bpmValsY(:,1)*1e6,'-k')
hold on
plot(bpmValsCorrected_C(:,1)*1e6,'-r')
plot(bpmValsY(:,:)*1e6,'-k')
plot(bpmValsCorrected_C(:,:)*1e6,'-r')

xlabel('BPM index')
ylabel('y [\mum]')
legend('Uncorrected','Conventional correction')

figure(2)
hold off
plot(bpmValsY(:,1)*1e6,'-k')
hold on
plot(bpmValsCorrected_C(:,1)*1e6, '-r')
plot(bpmValsCorrected_ML(:,1)*1e6,'-b')
plot(bpmValsY(:,:)*1e6,'-k')
plot(bpmValsCorrected_C(:,:)*1e6, '-r')
plot(bpmValsCorrected_ML(:,:)*1e6,'-b')

xlabel('BPM index')
ylabel('y [\mum]')
legend('Uncorrected','Conventional correction', 'ML correction')








