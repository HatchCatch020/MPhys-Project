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
numSeeds   = 10;
dcorrStrength   = 1e-3; % Tesla
kickStrength    = 1e-6;   % m

% Varying
ML_numobs       = 1000;
BPMnoise        = linspace(10,50,5)*1e-6;   % metres
FocusingError   = linspace(10,80,8)*1e-2;   % percentage
AlignmentError  = linspace(10,50,5)*1e-6;   % metres

% Setup fixed values for tests
fixedBPMnoise = 10e-6;
fixedFocusingError = 0.01;

% BPMs to be 

%% Test BPM noise
% % Vary the BPM noise with a fixed focusing error and 'alignment' error, use
% % multiple number of observations for ML.
% % Generate a plot of the rms values for the conventional method
% % and multiple ML_nomobs.
% 
% % Setup arrays
% rmsBPMvals_C        = zeros(numel(BPMnoise), 1, numSeeds);
% rmsBPMvals_ML       = zeros(numel(BPMnoise), numel(ML_numobs), numSeeds);
% mean_rmsBPMvals_C   = zeros(numel(BPMnoise),1);
% mean_rmsBPMvals_ML  = zeros(numel(BPMnoise), numel(ML_numobs));
% err_rmsBPMvals_C    = zeros(numel(BPMnoise),1);
% err_rmsBPMvals_ML   = zeros(numel(BPMnoise), numel(ML_numobs));
% 
% for i = 1:numel(BPMnoise)
%     % Calc resp mat with conventional method in 'perfect' model
%     RespMatC = lt.calcRespMatC(dcorrStrength);
%     
%     % Generate corrections for multiple different seeds
%     for k = 1:numSeeds
%         % set seed and turn on a set of errors
%         rng(k)
%         lt.setQuadFerrors(true, fixedFocusingError);
%         lt.setQuadAerrors(true, kickStrength);
%         
%         bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise(i);
%         
%         bpmValsCorrected_C  = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), BPMnoise(i));
%         rmsBPMvals_C(i,1,k) = rms(bpmValsCorrected_C);
%         
%         for j = 1:numel(ML_numobs)
%             RespMatML = lt.calcRespMatML(ML_numobs(j), BPMnoise(i), dcorrStrength);
%             bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), BPMnoise(i));
%             rmsBPMvals_ML(i,j,k) = rms(bpmValsCorrected_ML);
%         end
%         % Turn off errors
%         lt.setQuadFerrors(false, fixedFocusingError);
%         lt.setQuadAerrors(false, kickStrength);
%     end
%     % Calc the mean of rms values over all the seeds, and also the std of
%     % the rms values.
%     mean_rmsBPMvals_C(i)  =  mean(rmsBPMvals_C(i,1,:));
%     err_rmsBPMvals_C(i,:)   =   prctile(rmsBPMvals_C(i,1,:),[10 90]);
%     for j=1:numel(ML_numobs) 
%         mean_rmsBPMvals_ML(i,j) = mean(rmsBPMvals_ML(i,j,:));
%         err_rmsBPMvals_ML(i,:,j)  =  prctile(rmsBPMvals_ML(i,j,:),[10 90]);
%     end
% end
% 
% figure(1)
% hold off
% errorbar(BPMnoise,1e3*mean_rmsBPMvals_C(:),1e3*err_rmsBPMvals_C(:,1),1e3*err_rmsBPMvals_C(:,2), '--+r', 'DisplayName', 'Conventional')
% hold on
% for j=1:numel(ML_numobs)
%     errorbar(BPMnoise,1e3*mean_rmsBPMvals_ML(:,1,j),1e3*err_rmsBPMvals_ML(:,j),1e3*err_rmsBPMvals_ML(:,2,j), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(j)))
% end
% xlabel('BPM Noise')
% ylabel('RMS BPM Values (mm)')
% legend
% title(sprintf('(Focusing error = %d)', fixedFocusingError))
% 
% % generate a plot of the correction from the final interation 
% figure(2)
% subplot(2,1,1)
% hold off
% plot(1e3*bpmValsY,'-ok')
% hold on
% plot(1e3*bpmValsCorrected_ML,'--+r')
% plot(1e3*bpmValsCorrected_C,'-.xb')
% xlabel('BPM index')
% ylabel('Vertical position (mm)')
% legend('Original',sprintf('ML corrected numObs = %d', ML_numobs(end)), 'Conventionally corrected')
% title(sprintf('Focusing error scale = %d T/m, Alignment error scale = %d T, BPM noise scale = %d',fixedFocusingError, kickStrength, BPMnoise(end)))
% 
% subplot(2,1,2)


%% Test Alignment error

% Setup arrays
rmsBPMvals_ML       = zeros(numel(AlignmentError), numSeeds);
mean_rmsBPMvals_ML  = zeros(numel(AlignmentError),1);
err_rmsBPMvals_ML   = zeros(numel(AlignmentError),6);

for i = 1:numel(AlignmentError)
    
    % Generate corrections for multiple different seeds
    for k = 1:numSeeds
        % set seed and turn on a set of errors
        rng(k*i)
        lt.setQuadFerrors(true, fixedFocusingError);
        lt.setQuadAerrors(true, AlignmentError(i));
        
        bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*fixedBPMnoise;
        
        RespMatML = lt.calcRespMatML(ML_numobs, fixedBPMnoise, dcorrStrength);
        bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), fixedBPMnoise);
        rmsBPMvals_ML(i,k) = rms(bpmValsCorrected_ML);
        % Turn off errors
        lt.setQuadFerrors(false, fixedFocusingError);
        lt.setQuadAerrors(false, AlignmentError(i));
    end
    % Calc the mean of rms values over all the seeds, and also the std of
    % the rms values.
    mean_rmsBPMvals_ML(i) = mean(rmsBPMvals_ML(i,:));
    err_rmsBPMvals_ML(i,:)  =  prctile(rmsBPMvals_ML(i,:),[10 90 20 80 30 70]);
end

figure(3)
hold off
errorbar(AlignmentError*1e6,1e3*mean_rmsBPMvals_ML(:),1e3*err_rmsBPMvals_ML(:,1),1e3*err_rmsBPMvals_ML(:,2), '-.x', 'DisplayName', sprintf('ML (Observations = %d)', ML_numobs))

xlabel('Alignment Error [\mum]')
ylabel('RMS BPM Values [mm]')
legend
title(sprintf('(Focusing error = %g, BPM noise = %g)', fixedFocusingError, fixedBPMnoise))

% % generate a plot of the correction from the final interation 
% figure(3)
% subplot(2,1,1)
% hold off
% plot(1e3*bpmValsY,'-ok')
% hold on
% plot(1e3*bpmValsCorrected_ML,'--+r')
% plot(1e3*bpmValsCorrected_C,'-.xb')
% xlabel('BPM index')
% ylabel('Vertical position (mm)')
% legend('Original',sprintf('ML corrected numObs = %d', ML_numobs(end)), 'Conventionally corrected')
% title(sprintf('Focusing error scale = %d T/m, Alignment error scale = %d T, BPM noise scale = %d',fixedFocusingError, kickStrength, BPMnoise(end)))
% 
% subplot(2,1,2)

% %% Test alignment error
% 
% kickStrength    = linspace(10,90,9);   % Tesla
% kick = randn(numel(beamline.quaderrlist),numel(kickStrength)).*kickStrength;
% 
% % Setup arrays
% STDresiduals_Aerror_C = zeros(length(kickStrength), 1, numIterations);
% STDresiduals_Aerror_ML = zeros(length(kickStrength), length(ML_numobs), numIterations);
% avgSTD_Aerror_ML = zeros(length(kickStrength),numel(ML_numobs));
% avgSTD_Aerror_C  = zeros(length(kickStrength), 1);
% RespMatML = zeros(numel(beamline.bpmlist), numel(beamline.corrlist), numel(ML_numobs));
% 
% for i = 1:numel(kickStrength)
%     % Calc resp mat with conventional method in 'perfect' model
%     RespMatC = lt.calcRespMatC(dcorrStrength);
%     
%     % Turn on errors
%     lt.setQuadAerrors(true, kick(:,i));
%     lt.setQuadFerrors(true, fixedFocusingError);
%     
%     % Calc resp mat with ML in model with errors
%     for j = 1:numel(ML_numobs)
%         % Gen response matrix with ML method
%         RespMatML(:,:,j) = lt.calcRespMatML(ML_numobs(j), fixedBPMnoise, dcorrStrength);
%     end
%     
%     for k = 1:numIterations
%         bpmorbitdy = lt.track_getBPMreadings();
%         
%         corrected_C = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmorbitdy), fixedBPMnoise);
%         STDresiduals_Aerror_C(i,1,k)  = std(corrected_C);
%         
%         for j = 1:numel(ML_numobs)
%             corrected_ML = lt.getBPMvalues_corr(pinv(RespMatML(:,:,j),  1e-3)*(-bpmorbitdy), fixedBPMnoise);
%             STDresiduals_Aerror_ML(i,j,k) = std(corrected_ML);
%         end
%     end
%     
%     % Turn off errors
%     lt.setQuadFerrors(false, fixedFocusingError);
%     lt.setQuadAerrors(false, kick(:,i));
% end
% 
% for i=1:numel(kickStrength)
%     avgSTD_Aerror_C(i,1)     = mean(STDresiduals_Aerror_C(i,1,:));
%     avgSTD_Aerror_C_err(i,1) = std(STDresiduals_Aerror_C(i,1,:))/sqrt(numel(STDresiduals_Aerror_C(i,1,:)));
%     for j=1:numel(ML_numobs)
%         avgSTD_Aerror_ML(i,j)      = mean(STDresiduals_Aerror_ML(i,j,:));
%         avgSTD_Aerror_ML_err(i,j)  = std(STDresiduals_Aerror_ML(i,j,:))/sqrt(numel(STDresiduals_Aerror_ML(i,1,:)));
%     end
% end
% 
% figure(3)
% hold off
% errorbar(kickStrength,avgSTD_Aerror_C(:,1), avgSTD_Aerror_C_err(:,1), '--+r', 'DisplayName', 'Conventional')
% hold on
% for j=1:numel(ML_numobs)
%     errorbar(kickStrength,avgSTD_Aerror_ML(:,j), avgSTD_Aerror_ML_err(:,j), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(j)))
% end
% xlabel('Alignment error')
% ylabel('Standard Deviation')
% legend
% title(sprintf('Standard deviation of the correction (BPM noise = %6.5g, focusing error = %6.5g)', fixedBPMnoise, fixedFocusingError))
% 
% figure(4)
% subplot(2,1,1)
% hold off
% plot(1e3*bpmorbitdy,'-ok')
% hold on
% plot(1e3*corrected_ML,'--+r')
% plot(1e3*corrected_C,'-.xb')
% xlabel('BPM index')
% ylabel('Vertical position (mm)')
% legend('Original','ML corrected', 'Conventionally corrected')
% title('Change in vertical position from random correctors and corrected positions')
% 
% subplot(2,1,2)
% hold off
% plot(1e3*(corrected_ML - bpmorbitdy),'--+r')
% hold on
% plot(1e3*(corrected_C - bpmorbitdy),'-.xb')
% xlabel('BPM index')
% ylabel('Residual (mm)')
% legend('ML corrected', 'Conventionally corrected')
% title('Difference between corrected and uncorrected positions')



