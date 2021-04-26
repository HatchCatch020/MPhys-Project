% Investigation of the scale of dcorrector strength used to calculate
% response matrices and the resulting change in the orbit from such
% corrector strengths.

%% Define the beamline
DefineCLARABeamline
clearvars -except beamline driftlist quadlist quaderrlist corrlist bpmlist Lcavitylist Sbendlist beam bl

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

numSeeds = 100;

dcorrStrength   = linspace(1,150,21)*1e-6;%[1 4.5 10 45 100 145]*1e-6;    % Tesla
BPMnoise        = [10 40 90]*1e-6;                  % Metres

FocusingError   = 0.005;                            % Percentage
AlignmentError  = 10e-6;                            % Metres

%%

% lt.setQuadFerrors(true, FocusingError);
% lt.setQuadAerrors(true, AlignmentError);
% 
% originalTraj     = lt.track_getBPMreadings();
% RMSoriginalTraj  = rms(originalTraj);
% 
% RMScorrectorsTraj       = zeros(numSeeds,1);
% RMScorrBPMTraj          = zeros(numSeeds, numel(BPMnoise));
% meanRMScorrectorsTraj   = zeros(numel(dcorrStrength),1);
% meanRMScorrBPMTraj      = zeros(numel(dcorrStrength),numel(BPMnoise));
% ERRmeanRMScorrectorsTraj= zeros(numel(dcorrStrength),2);
% ERRmeanRMScorrBPMTraj   = zeros(numel(dcorrStrength),numel(BPMnoise),2);

% for i=1:numel(dcorrStrength)
%     for n=1:numSeeds
%         dcorr       =  dcorrStrength(i)*randn(numel(beamline.corrlist),1);
% 
%         correctorsTraj  = lt.getBPMvalues_corr(dcorr, 0);
%         RMScorrectorsTraj(n)  = rms(correctorsTraj);
% 
%         for j=1:(numel(BPMnoise))
%             RMScorrBPMTraj(n,j) = rms(correctorsTraj + randn(numel(beamline.bpmlist),1)*BPMnoise(j));
%         end
%     end
%     meanRMScorrectorsTraj(i) = mean(RMScorrectorsTraj);
%     
%     ERRmeanRMScorrectorsTraj(i,1) = prctile(RMScorrectorsTraj, [10]);
%     ERRmeanRMScorrectorsTraj(i,2) = prctile(RMScorrectorsTraj, [90]);
%     for j=1:(numel(BPMnoise))
%         meanRMScorrBPMTraj(i,j) = mean(RMScorrBPMTraj(:,j)); 
%         
%         ERRmeanRMScorrBPMTraj(i,j,1) = prctile(RMScorrBPMTraj(:,j), [10]);
%         ERRmeanRMScorrBPMTraj(i,j,2) = prctile(RMScorrBPMTraj(:,j), [90]);
%     end
% end
% 
% figure(1)
% hold off
% plot(dcorrStrength*1e6,meanRMScorrectorsTraj*1e3, '-+k','DisplayName', 'Trajectory with correctors')
% hold on
% plot(dcorrStrength*1e6,meanRMScorrBPMTraj(:,1)*1e3, '-.r', 'DisplayName', 'BPMnoise = 10e-6')
% plot(dcorrStrength*1e6,meanRMScorrBPMTraj(:,2)*1e3, '--ok','DisplayName', 'BPMnoise = 40e-6')
% plot(dcorrStrength*1e6,meanRMScorrBPMTraj(:,3)*1e3, '-.b', 'DisplayName', 'BPMnoise = 90e-6')
% 
% shade(dcorrStrength*1e6,meanRMScorrBPMTraj(:,1)*1e3,dcorrStrength*1e6,meanRMScorrBPMTraj(:,3)*1e3,'FillType',[2 1],'LineStyle', 'none')
% 
% xlabel('\Deltacorr strength [\muT]')
% ylabel('RMS BPM values [mm]')
% legend('Trajectory with correctors', 'BPMnoise = 10 \mum','BPMnoise = 40 \mum','BPMnoise = 90 \mum', 'Location', 'northwest')
% %plot(meanRMScorrectorsTraj, 'DisplayName', 'Trajectory with correctors')

%%
lt.setQuadFerrors(true, FocusingError);
lt.setQuadAerrors(true, AlignmentError);

originalTraj     = lt.track_getBPMreadings();
RMSoriginalTraj  = rms(originalTraj);

RMScorrectorsTraj       = zeros(numSeeds,1);
RMScorrBPMTraj          = zeros(numSeeds, numel(BPMnoise));
meanRMScorrectorsTraj   = zeros(numel(dcorrStrength),1);
meanRMScorrBPMTraj      = zeros(numel(dcorrStrength),numel(BPMnoise));
ERRmeanRMScorrectorsTraj= zeros(numel(dcorrStrength),2);
ERRmeanRMScorrBPMTraj   = zeros(numel(dcorrStrength),numel(BPMnoise),2);

for i=1:numel(dcorrStrength)
    for n=1:numSeeds
        dcorr       =  dcorrStrength(i)*randn(numel(beamline.corrlist),1);

        correctorsTraj  = lt.getBPMvalues_corr(dcorr, 0);
        RMScorrectorsTraj(n)  = rms(correctorsTraj);

        for j=1:(numel(BPMnoise))
            RMScorrBPMTraj(n,j) = rms((correctorsTraj + randn(numel(beamline.bpmlist),1)*BPMnoise(j)));
        end
    end
    meanRMScorrectorsTraj(i) = mean(RMScorrectorsTraj);
    
    ERRmeanRMScorrectorsTraj(i,1) = prctile(RMScorrectorsTraj, [10]);
    ERRmeanRMScorrectorsTraj(i,2) = prctile(RMScorrectorsTraj, [90]);
    for j=1:(numel(BPMnoise))
        meanRMScorrBPMTraj(i,j) = mean(RMScorrBPMTraj(:,j)); 
        
        ERRmeanRMScorrBPMTraj(i,j,1) = prctile(RMScorrBPMTraj(:,j), [10]);
        ERRmeanRMScorrBPMTraj(i,j,2) = prctile(RMScorrBPMTraj(:,j), [90]);
    end
end

figure(1)
hold off
plot(dcorrStrength*1e6,meanRMScorrectorsTraj*1e3, '-+k','DisplayName', 'Trajectory with correctors')
hold on
plot(dcorrStrength*1e6,meanRMScorrBPMTraj(:,1)*1e3, '-.r', 'DisplayName', 'BPMnoise = 10e-6')
plot(dcorrStrength*1e6,meanRMScorrBPMTraj(:,2)*1e3, '--ok','DisplayName', 'BPMnoise = 40e-6')
plot(dcorrStrength*1e6,meanRMScorrBPMTraj(:,3)*1e3, '-.b', 'DisplayName', 'BPMnoise = 90e-6')

shade(dcorrStrength*1e6,meanRMScorrBPMTraj(:,1)*1e3,dcorrStrength*1e6,meanRMScorrBPMTraj(:,3)*1e3,'FillType',[2 1],'LineStyle', 'none')

xlabel('\Deltacorr strength [\muT]')
ylabel('\Deltay_{rms} (BPM values) [mm]')
legend('Trajectory with correctors', 'BPMnoise = 10 \mum','BPMnoise = 40 \mum','BPMnoise = 90 \mum', 'Location', 'northwest')

f0 = fit((dcorrStrength*1e6)',meanRMScorrectorsTraj*1e3, 'poly3');
f1 = fit((dcorrStrength*1e6)',meanRMScorrBPMTraj(:,1)*1e3, 'poly3');
f2 = fit((dcorrStrength*1e6)',meanRMScorrBPMTraj(:,2)*1e3, 'poly3');
f3 = fit((dcorrStrength*1e6)',meanRMScorrBPMTraj(:,3)*1e3, 'poly3');

figure(2)
hold off
plot(dcorrStrength*1e6,meanRMScorrectorsTraj*1e3, '-+k','DisplayName', 'Trajectory with correctors')
hold on
plot(f1, '-r')
plot(f2, '--k')
plot(f3, '-b')


%plot(meanRMScorrectorsTraj, 'DisplayName', 'Trajectory with correctors')

