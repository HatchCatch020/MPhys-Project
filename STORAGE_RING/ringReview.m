% C. Monaghan
% 10/04/2021
% A systematic review of a storage ring and the use of machine learning to 
% determine the response matrix for orbit correction.
%% Define the beamline
DefineBeamline
clearvars -except beamline

% Set the rf frequency
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Setup
rt = ringTools;
rt.beamline = beamline;
rt.ML_Algorithm = 'mvn';

%% Variables

% Fixed variables
numSeeds        = 60;
dcorrStrength   = 1e-6; % tesla (a)

ML_numobs       = 60;
BPMnoise        = linspace(1,40,11)*1e-6;   % metre
FocusingError   = 0.01;                     % % tesla   (dg)
AlignmentError  = linspace(1,40,11)*1e-6;   % metre     (dy)

% Setup fixed values for tests
fixedBPMnoise      = 20e-6; % metre
fixedFocusingError = 0.01;  % % tesla   (dg)
fixedAlignment     = 50e-6; % metre     (dy)

%% Test BPM noise

% Setup arrays
STDresiduals_BPMnoise_C  = zeros(length(BPMnoise), 1, numIterations);
STDresiduals_BPMnoise_ML = zeros(length(BPMnoise), length(ML_numobs), numIterations);
avgSTD_BPMnoise_ML       = zeros(length(BPMnoise), numel(ML_numobs));
avgSTD_BPMnoise_C        = zeros(length(BPMnoise), 1);
RespMatML                = zeros(numel(beamline.bpmlist), numel(beamline.corrlist), numel(ML_numobs));

for i = 1:numel(BPMnoise)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = rt.calcRespMatC(dcorrStrength);
    
    % Turn on errors
    rt.setQuadFerrors(true, fixedFocusingError);
    rt.setQuadAerrors(true, kick);
    
    % Calc resp mat with ML in model with errors
    for j = 1:numel(ML_numobs)
        % Gen response matrix with ML method
        RespMatML(:,:,j) = rt.calcRespMatML(ML_numobs(j), BPMnoise(i), dcorrStrength);
    end
    
    for k = 1:numIterations
        bpmorbitdy = rt.track_getBPMreadings();
        
        corrected_C = rt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmorbitdy), BPMnoise(i));
        STDresiduals_BPMnoise_C(i,1,k)  = std(corrected_C);
        
        for j = 1:numel(ML_numobs)
            corrected_ML = rt.getBPMvalues_corr(pinv(RespMatML(:,:,j),  1e-3)*(-bpmorbitdy), BPMnoise(i));
            STDresiduals_BPMnoise_ML(i,j,k) = std(corrected_ML);
        end
    end
    
    % Turn off errors
    rt.setQuadFerrors(false, fixedFocusingError);
    rt.setQuadAerrors(false, kick);
end

for i=1:numel(BPMnoise)
    avgSTD_BPMnoise_C(i,1) = mean(STDresiduals_BPMnoise_C(i,1,:));
    for j=1:numel(ML_numobs)
        avgSTD_BPMnoise_ML(i,j) = mean(STDresiduals_BPMnoise_ML(i,j,:));
    end
end

figure(1)
hold off
plot(BPMnoise,avgSTD_BPMnoise_C(:,1), '--+r', 'DisplayName', 'Conventional')
hold on
for j=1:numel(ML_numobs)
    plot(BPMnoise,avgSTD_BPMnoise_ML(:,j), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(j)))
end
xlabel('BPM Noise')
ylabel('Standard Deviation')
legend
title(sprintf('Standard deviation of the correction (focusing error = %d)', fixedFocusingError))

figure(2)
subplot(2,1,1)
hold off
plot(1e3*bpmorbitdy,'-ok')
hold on
plot(1e3*corrected_ML,'--+r')
plot(1e3*corrected_C,'-.xb')
xlabel('BPM index')
ylabel('Vertical position (mm)')
legend('Original','ML corrected', 'Conventionally corrected')
title('Change in vertical position from random correctors and corrected positions')

subplot(2,1,2)
hold off
plot(1e3*(corrected_ML - bpmorbitdy),'--+r')
hold on
plot(1e3*(corrected_C - bpmorbitdy),'-.xb')
xlabel('BPM index')
ylabel('Residual (mm)')
legend('ML corrected', 'Conventionally corrected')
title('Difference between corrected and uncorrected positions')

%% Test Alignment error

% Setup arrays
rmsBPMvals          = zeros(numel(AlignmentError), numSeeds);
rmsBPMvals_C        = zeros(numel(AlignmentError), numSeeds);
rmsBPMvals_ML       = zeros(numel(AlignmentError), numSeeds);
mean_rmsBPMvals     = zeros(numel(AlignmentError),1);
mean_rmsBPMvals_C   = zeros(numel(AlignmentError),1);
mean_rmsBPMvals_ML  = zeros(numel(AlignmentError),1);
err_rmsBPMvals      = zeros(numel(AlignmentError),6);
err_rmsBPMvals_C    = zeros(numel(AlignmentError),6);
err_rmsBPMvals_ML   = zeros(numel(AlignmentError),6);

for i = 1:numel(AlignmentError)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = rt.calcRespMatC(dcorrStrength);
    
    % Generate corrections for multiple different seeds
    for k = 1:numSeeds
        % set seed and turn on a set of errors
        rng(k*i)
        rt.setQuadFerrors(true, fixedFocusingError);
        rt.setQuadAerrors(true, AlignmentError(i));
        
        bpmValsY = rt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*fixedBPMnoise;
        rmsBPMvals(i,k) = rms(bpmValsY);
        
        bpmValsCorrected_C  = rt.getBPMvalues_corr(pinv(RespMatC,  1e-2)*(-bpmValsY), fixedBPMnoise);
        rmsBPMvals_C(i,k) = rms(bpmValsCorrected_C);
        
        RespMatML = rt.calcRespMatML(ML_numobs, fixedBPMnoise, dcorrStrength);
        bpmValsCorrected_ML  = rt.getBPMvalues_corr(pinv(RespMatML,  1e-2)*(-bpmValsY), fixedBPMnoise);
        rmsBPMvals_ML(i,k) = rms(bpmValsCorrected_ML);
       
        % Turn off errors
        rt.setQuadFerrors(false, fixedFocusingError);
        rt.setQuadAerrors(false, AlignmentError(i));
    end
    % Calc the mean of rms values over all the seeds, and also the std of
    % the rms values.
    mean_rmsBPMvals(i)      = mean(rmsBPMvals(i,:));
    err_rmsBPMvals(i,:)     = prctile(rmsBPMvals(i,:),[10 90 20 80 30 70]);
    mean_rmsBPMvals_C(i)    = mean(rmsBPMvals_C(i,:));
    err_rmsBPMvals_C(i,:)   = prctile(rmsBPMvals_C(i,:),[10 90 20 80 30 70]);
    mean_rmsBPMvals_ML(i)   = mean(rmsBPMvals_ML(i,:));
    err_rmsBPMvals_ML(i,:)  = prctile(rmsBPMvals_ML(i,:),[10 90 20 80 30 70]);
end

figure(4)
hold off
plot(AlignmentError*1e6,1e3*mean_rmsBPMvals_C(:),'--+r', 'DisplayName', 'Conventional')
hold on
plot(AlignmentError*1e6,1e3*mean_rmsBPMvals_ML(:), '-.x', 'DisplayName', sprintf('ML (Observations = %d)', ML_numobs))
shade(AlignmentError*1e6,1e3*err_rmsBPMvals_C(:,3),AlignmentError*1e6,1e3*err_rmsBPMvals_C(:,4),'FillType',[2 1],'LineStyle', 'none','FillColor', 'r')
shade(AlignmentError*1e6,1e3*err_rmsBPMvals_ML(:,5),AlignmentError*1e6,1e3*err_rmsBPMvals_ML(:,6),'FillType',[2 1],'LineStyle', 'none', 'FillColor', 'cyan')
xlabel('Alignment Error [\mum]')
ylabel('RMS BPM Values [mm]')
%xlim([1 20]); ylim([0 0.07])
legend('Conventional',sprintf('ML (Observations = %d)', ML_numobs))
%title(sprintf('(Focusing error = %g, BPM noise = %g)', fixedFocusingError, fixedBPMnoise))

% generate a plot of the correction from the final interation 
figure(3)
subplot(2,1,1)
hold off
plot(1e3*bpmValsY,'-ok')
hold on
plot(1e3*bpmValsCorrected_ML,'--+r')
plot(1e3*bpmValsCorrected_C,'-.xb')
xlabel('BPM index')
ylabel('Vertical position (mm)')
legend('Original',sprintf('ML corrected numObs = %d', ML_numobs), 'Conventionally corrected')
%title(sprintf('Focusing error scale = %d T/m, Alignment error scale = %d T, BPM noise scale = %d',fixedFocusingError, kickStrength, BPMnoise(end)))

subplot(2,1,2)

