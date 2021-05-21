% C. Monaghan
% 27/03/2021
% A systematic review of the CLARA beamline and the use of machine learning
% to determine the response matrix for orbit correction.
%% Define the beamline
DefineCLARABeamline
clearvars -except beamline driftlist quadflist quadflist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
% DefineCLARABeamline
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Setup
% Load the functions which will be used in this script.
% Access the functions using the object lt
lt = linTools;
lt.beamline = beamline;
lt.ML_Algorithm = 'mvn';

%% Variables

% Fixed variables
numSeeds   = 60;
dcorrStrength   = 100e-6; % tesla (a)

% Varying
ML_numobs       = 60;
BPMnoise        = linspace(1,40,11)*1e-6;   % metre
FocusingError   = 0.005;                    % % tesla   (dg)
AlignmentError  = linspace(1,40,11)*1e-6;   % metre     (dy)

% Setup fixed values for tests
fixedBPMnoise = 20e-6;      % metre
fixedFocusingError = 0.01;  % % tesla   (dg)
fixedAlignment     = 50e-6; % metre     (dy)

%% Test BPM noise
% Vary the BPM noise with a fixed focusing error and 'alignment' error, use
% multiple number of observations for ML.
% Generate a plot of the rms values for the conventional method
% and multiple ML_nomobs.

% Setup arrays
rmsBPMvals_C        = zeros(numel(BPMnoise), numSeeds);
rmsBPMvals_ML       = zeros(numel(BPMnoise), numSeeds);
mean_rmsBPMvals_C   = zeros(numel(BPMnoise),1);
mean_rmsBPMvals_ML  = zeros(numel(BPMnoise),1);
err_rmsBPMvals_C    = zeros(numel(BPMnoise),4);
err_rmsBPMvals_ML   = zeros(numel(BPMnoise),4);

for i = 1:numel(BPMnoise)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = lt.calcRespMatC(dcorrStrength);
    
    % Generate corrections for multiple different seeds
    for k = 1:numSeeds
        % set seed and turn on a set of errors
        rng(k)
        lt.setQuadFerrors(true, fixedFocusingError);
        lt.setQuadAerrors(true, fixedAlignment);
        
        bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise(i);
        
        bpmValsCorrected_C  = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), BPMnoise(i));
        rmsBPMvals_C(i,k) = rms(bpmValsCorrected_C);
        
        RespMatML = lt.calcRespMatML(ML_numobs, BPMnoise(i), dcorrStrength);
        bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), BPMnoise(i));
        rmsBPMvals_ML(i,k) = rms(bpmValsCorrected_ML);
            
        % Turn off errors
        lt.setQuadFerrors(false, fixedFocusingError);
        lt.setQuadAerrors(false, fixedAlignment);
    end
    % Calc the mean of rms values over all the seeds, and also the std of
    % the rms values.
    mean_rmsBPMvals_C(i)  =  mean(rmsBPMvals_C(i,:));
    err_rmsBPMvals_C(i,:)   =   prctile(rmsBPMvals_C(i,:),[10 90 20 80]);
    mean_rmsBPMvals_ML(i) = mean(rmsBPMvals_ML(i,:));
    err_rmsBPMvals_ML(i,:)  =  prctile(rmsBPMvals_ML(i,:),[10 90 20 80]);
end

figure(1)
hold off
plot(BPMnoise*1e6,1e3*mean_rmsBPMvals_C(:),'--+r', 'DisplayName', 'Conventional')
hold on
plot(BPMnoise*1e6,1e3*mean_rmsBPMvals_ML(:), ':*', 'DisplayName', 'ML')
shade(BPMnoise*1e6,1e3*err_rmsBPMvals_C(:,3),BPMnoise*1e6,1e3*err_rmsBPMvals_C(:,4),'FillType',[2 1],'LineStyle', 'none','FillColor', 'r')
shade(BPMnoise*1e6,1e3*err_rmsBPMvals_ML(:,3),BPMnoise*1e6,1e3*err_rmsBPMvals_ML(:,4),'FillType',[2 1],'LineStyle', 'none', 'FillColor', 'cyan')
xlabel('BPM noise [\mum]')
ylabel('RMS y [mm]')
%xlim([1 20]); ylim([0 0.07])
pbaspect([1 1 1])
legend('Conventional','ML')
title(sprintf('(Focusing error = %g, BPM noise = %g)', fixedFocusingError, fixedBPMnoise))

% Generate a plot of the correction from the final interation 
figure(2)
hold off
plot(1e3*bpmValsY,'-ok')
hold on
plot(1e3*bpmValsCorrected_ML,'--+r')
plot(1e3*bpmValsCorrected_C,'-.xb')
xlabel('BPM index')
ylabel('Vertical position (mm)')
legend('Original',sprintf('ML corrected numObs = %d', ML_numobs(end)), 'Conventionally corrected')
title(sprintf('Focusing error scale = %d T/m, Alignment error scale = %d T, BPM noise scale = %d',fixedFocusingError, fixedAlignment, BPMnoise(end)))



%% Test Alignment error

% Setup arrays
rmsBPMvals_C        = zeros(numel(AlignmentError), 1, numSeeds);
rmsBPMvals_ML       = zeros(numel(AlignmentError), numel(ML_numobs), numSeeds);
mean_rmsBPMvals_C   = zeros(numel(AlignmentError),1);
mean_rmsBPMvals_ML  = zeros(numel(AlignmentError), numel(ML_numobs));
err_rmsBPMvals_C    = zeros(numel(AlignmentError),6);
err_rmsBPMvals_ML   = zeros(numel(AlignmentError),6, numel(ML_numobs));

for i = 1:numel(AlignmentError)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = lt.calcRespMatC(dcorrStrength);
    
    % Generate corrections for multiple different seeds
    for k = 1:numSeeds
        % set seed and turn on a set of errors
        rng(k*i)
        lt.setQuadFerrors(true, fixedFocusingError);
        lt.setQuadAerrors(true, AlignmentError(i));
        
        bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*fixedBPMnoise;
        
        bpmValsCorrected_C  = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), fixedBPMnoise);
        rmsBPMvals_C(i,1,k) = rms(bpmValsCorrected_C);
        
        for j = 1:numel(ML_numobs)
%            dcorrStrength_alt = [1e-6 1e-5 1e-4];
            RespMatML = lt.calcRespMatML(ML_numobs(j), fixedBPMnoise, dcorrStrength);
%            RespMatML = lt.calcRespMatML(ML_numobs(j), fixedBPMnoise, dcorrStrength);
            bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), fixedBPMnoise);
            rmsBPMvals_ML(i,j,k) = rms(bpmValsCorrected_ML);
        end
        % Turn off errors
        lt.setQuadFerrors(false, fixedFocusingError);
        lt.setQuadAerrors(false, AlignmentError(i));
    end
    % Calc the mean of rms values over all the seeds, and also the std of
    % the rms values.
    mean_rmsBPMvals_C(i)  =  mean(rmsBPMvals_C(i,1,:));
    err_rmsBPMvals_C(i,:)   =   prctile(rmsBPMvals_C(i,1,:),[10 90 20 80 30 70]);
    for j=1:numel(ML_numobs) 
        mean_rmsBPMvals_ML(i,j) = mean(rmsBPMvals_ML(i,j,:));
        err_rmsBPMvals_ML(i,:,j)  =  prctile(rmsBPMvals_ML(i,j,:),[10 90 20 80 30 70]);
    end
end

figure(1)
hold off
plot(AlignmentError*1e6,1e3*mean_rmsBPMvals_C(:),'--+r', 'DisplayName', 'Conventional')
hold on
for j=1:numel(ML_numobs)
    plot(AlignmentError*1e6,1e3*mean_rmsBPMvals_ML(:), ':*', 'DisplayName', sprintf('ML', ML_numobs(j)))
end
shade(AlignmentError*1e6,1e3*err_rmsBPMvals_C(:,5),AlignmentError*1e6,1e3*err_rmsBPMvals_C(:,6),'FillType',[2 1],'LineStyle', 'none','FillColor', 'r')
shade(AlignmentError*1e6,1e3*err_rmsBPMvals_ML(:,5),AlignmentError*1e6,1e3*err_rmsBPMvals_ML(:,6),'FillType',[2 1],'LineStyle', 'none', 'FillColor', 'cyan')
xlabel('\Delta{}y [\mum]')
ylabel('RMS y [mm]')
%xlim([1 20]); ylim([0 0.07])
legend('Conventional','ML')
pbaspect([1 1 1])
title(sprintf('(Focusing error = %g, BPM noise = %g)', fixedFocusingError, fixedBPMnoise))

% generate a plot of the correction from the final interation 
figure(3)
hold off
plot(1e3*bpmValsY,'-ok')
hold on
plot(1e3*bpmValsCorrected_ML,'--+r')
plot(1e3*bpmValsCorrected_C,'-.xb')
xlabel('BPM index')
ylabel('Vertical position (mm)')
legend('Original',sprintf('ML corrected numObs = %d', ML_numobs(end)), 'Conventionally corrected')
title(sprintf('Focusing error scale = %d T/m, Alignment error scale = %d T, BPM noise scale = %d',fixedFocusingError, kickStrength, BPMnoise(end)))

