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

rng(2)

%% Variables

% Fixed variables
numIterations   = 100;
dcorrStrength   = 1e-4; % Tesla
kickStrength    = 10;   % Tesla
kick = randn(numel(beamline.quaderrlist),1)*kickStrength;

% Varying
ML_numobs       = linspace(100,1000,2);
BPMnoise        = linspace(10,70,7)*1e-6;   % unit
FocusingError   = linspace(10,90,9)*1e-2;   % unit

% Setup fixed values for tests
fixedBPMnoise = 70e-6;
fixedFocusingError = 0.01;

% %% Test BPM noise
% % Vary the BPM noise with a fixed focusing error and 'alignment' error, use
% % multiple number of observations for ML.
% % Generate a plot of the standard deviations for the conventional method
% % and multiple ML_nomobs.
% 
% % Setup fixed values
% fixedFocusingError = 0.01;
% 
% % Setup arrays
% STDresiduals_BPMnoise_C = zeros(length(BPMnoise), 1, numIterations);
% STDresiduals_BPMnoise_ML = zeros(length(BPMnoise), length(ML_numobs), numIterations);
% avgSTD_BPMnoise_ML = zeros(length(BPMnoise),numel(ML_numobs));
% avgSTD_BPMnoise_C  = zeros(length(BPMnoise), 1);
% RespMatML = zeros(numel(beamline.bpmlist), numel(beamline.corrlist), numel(ML_numobs));
% 
% for i = 1:numel(BPMnoise)
%     % Calc resp mat with conventional method in 'perfect' model
%     RespMatC = lt.calcRespMatC(dcorrStrength);
%     
%     % Turn on errors
%     lt.setQuadFerrors(true, fixedFocusingError);
%     lt.setQuadAerrors(true, kick);
%     
%     % Calc resp mat with ML in model with errors
%     for j = 1:numel(ML_numobs)
%         % Gen response matrix with ML method
%         RespMatML(:,:,j) = lt.calcRespMatML(ML_numobs(j), BPMnoise(i), dcorrStrength);
%     end
%     
%     for k = 1:numIterations
%         bpmorbitdy = lt.track_getBPMreadings();
%         
%         corrected_C = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmorbitdy), BPMnoise(i));
%         STDresiduals_BPMnoise_C(i,1,k)  = std(corrected_C);
%         
%         for j = 1:numel(ML_numobs)
%             corrected_ML = lt.getBPMvalues_corr(pinv(RespMatML(:,:,j),  1e-3)*(-bpmorbitdy), BPMnoise(i));
%             STDresiduals_BPMnoise_ML(i,j,k) = std(corrected_ML);
%         end
%     end
%     
%     % Turn off errors
%     lt.setQuadFerrors(false, fixedFocusingError);
%     lt.setQuadAerrors(false, kickStrength);
% end
% 
% for i=1:numel(BPMnoise)
%     avgSTD_BPMnoise_C(i,1) = mean(STDresiduals_BPMnoise_C(i,1,:));
%     for j=1:numel(ML_numobs)
%         avgSTD_BPMnoise_ML(i,j) = mean(STDresiduals_BPMnoise_ML(i,j,:));
%     end
% end
% 
% figure(1)
% hold off
% plot(BPMnoise,avgSTD_BPMnoise_C(:,1), '--+r', 'DisplayName', 'Conventional')
% hold on
% for j=1:numel(ML_numobs)
%     plot(BPMnoise,avgSTD_BPMnoise_ML(:,j), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(j)))
% end
% xlabel('BPM Noise')
% ylabel('Standard Deviation')
% legend
% title(sprintf('Standard deviation of the correction (focusing error = %d)', fixedFocusingError))
% 
% figure(2)
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
% 

%% Test focusing error

% Setup arrays
STDresiduals_Ferror_C = zeros(length(FocusingError), 1, numIterations);
STDresiduals_Ferror_ML = zeros(length(FocusingError), length(ML_numobs), numIterations);
avgSTD_Ferror_ML = zeros(length(FocusingError),numel(ML_numobs));
avgSTD_Ferror_C  = zeros(length(FocusingError), 1);
RespMatML = zeros(numel(beamline.bpmlist), numel(beamline.corrlist), numel(ML_numobs));

for i = 1:numel(FocusingError)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = lt.calcRespMatC(dcorrStrength);
    
    % Turn on errors
    lt.setQuadFerrors(true, FocusingError(i));
    lt.setQuadAerrors(true, kick);
    
    % Calc resp mat with ML in model with errors
    for j = 1:numel(ML_numobs)
        % Gen response matrix with ML method
        RespMatML(:,:,j) = lt.calcRespMatML(ML_numobs(j), fixedBPMnoise, dcorrStrength);
    end
    
    for k = 1:numIterations
        bpmorbitdy = lt.track_getBPMreadings();
        
        corrected_C = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmorbitdy), fixedBPMnoise);
        STDresiduals_Ferror_C(i,1,k)  = std(corrected_C);
        
        for j = 1:numel(ML_numobs)
            corrected_ML = lt.getBPMvalues_corr(pinv(RespMatML(:,:,j),  1e-3)*(-bpmorbitdy), fixedBPMnoise);
            STDresiduals_Ferror_ML(i,j,k) = std(corrected_ML);
        end
    end
    
    % Turn off errors
    lt.setQuadFerrors(false, FocusingError);
    lt.setQuadAerrors(false, kickStrength);
end

for i=1:numel(FocusingError)
    avgSTD_Ferror_C(i,1)     = mean(STDresiduals_Ferror_C(i,1,:));
    avgSTD_Ferror_C_err(i,1) = std(STDresiduals_Ferror_C(i,1,:))/sqrt(numel(STDresiduals_Ferror_C(i,1,:)));
    for j=1:numel(ML_numobs)
        avgSTD_Ferror_ML(i,j)      = mean(STDresiduals_Ferror_ML(i,j,:));
        avgSTD_Ferror_ML_err(i,j)  = std(STDresiduals_Ferror_ML(i,j,:))/sqrt(numel(STDresiduals_Ferror_ML(i,1,:)));
    end
end

figure(1)
hold off
errorbar(FocusingError,avgSTD_Ferror_C(:,1), avgSTD_Ferror_C_err(:,1), '--+r', 'DisplayName', 'Conventional')
hold on
for j=1:numel(ML_numobs)
    errorbar(FocusingError,avgSTD_Ferror_ML(:,j), avgSTD_Ferror_ML_err(:,j), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(j)))
end
xlabel('Focusing error')
ylabel('Standard Deviation')
legend
title(sprintf('Standard deviation of the correction (BPM noise = %6.5g, alignment error scale = %6.5g)', fixedBPMnoise,kickStrength))

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

%% Test alignment error

kickStrength    = linspace(1,90,9);   % Tesla
kick = randn(numel(beamline.quaderrlist),numel(kickStrength)).*kickStrength;

% Setup arrays
STDresiduals_Aerror_C = zeros(length(kickStrength), 1, numIterations);
STDresiduals_Aerror_ML = zeros(length(kickStrength), length(ML_numobs), numIterations);
avgSTD_Aerror_ML = zeros(length(kickStrength),numel(ML_numobs));
avgSTD_Aerror_C  = zeros(length(kickStrength), 1);
RespMatML = zeros(numel(beamline.bpmlist), numel(beamline.corrlist), numel(ML_numobs));

for i = 1:numel(kickStrength)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = lt.calcRespMatC(dcorrStrength);
    
    % Turn on errors
    lt.setQuadFerrors(true, fixedFocusingError);
    lt.setQuadAerrors(true, kick(:,i));
    
    % Calc resp mat with ML in model with errors
    for j = 1:numel(ML_numobs)
        % Gen response matrix with ML method
        RespMatML(:,:,j) = lt.calcRespMatML(ML_numobs(j), fixedBPMnoise, dcorrStrength);
    end
    
    for k = 1:numIterations
        bpmorbitdy = lt.track_getBPMreadings();
        
        corrected_C = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmorbitdy), fixedBPMnoise);
        STDresiduals_Aerror_C(i,1,k)  = std(corrected_C);
        
        for j = 1:numel(ML_numobs)
            corrected_ML = lt.getBPMvalues_corr(pinv(RespMatML(:,:,j),  1e-3)*(-bpmorbitdy), fixedBPMnoise);
            STDresiduals_Aerror_ML(i,j,k) = std(corrected_ML);
        end
    end
    
    % Turn off errors
    lt.setQuadFerrors(false, fixedFocusingError);
    lt.setQuadAerrors(false, kick(:,i));
end

for i=1:numel(kickStrength)
    avgSTD_Aerror_C(i,1)     = mean(STDresiduals_Aerror_C(i,1,:));
    avgSTD_Aerror_C_err(i,1) = std(STDresiduals_Aerror_C(i,1,:))/sqrt(numel(STDresiduals_Ferror_C(i,1,:)));
    for j=1:numel(ML_numobs)
        avgSTD_Aerror_ML(i,j)      = mean(STDresiduals_Aerror_ML(i,j,:));
        avgSTD_Aerror_ML_err(i,j)  = std(STDresiduals_Aerror_ML(i,j,:))/sqrt(numel(STDresiduals_Ferror_ML(i,1,:)));
    end
end

figure(3)
hold off
errorbar(kickStrength,avgSTD_Aerror_C(:,1), avgSTD_Aerror_C_err(:,1), '--+r', 'DisplayName', 'Conventional')
hold on
for j=1:numel(ML_numobs)
    errorbar(kickStrength,avgSTD_Aerror_ML(:,j), avgSTD_Aerror_ML_err(:,j), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(j)))
end
xlabel('Alignment error')
ylabel('Standard Deviation')
legend
title(sprintf('Standard deviation of the correction (BPM noise = %6.5g, focusing error = %6.5g)', fixedBPMnoise, fixedFocusingError))

figure(4)
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




