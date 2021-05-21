% C. Monaghan
% CLARAdcorrTest

% Investigation of the scale of dcorrector strength used to calculate
% response matrices and the resulting change in the orbit from such
% corrector strengths. 

% Otherwise known as paramter (a) in the report.

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

dcorrStrength   = linspace(1,150,21)*1e-6       % tesla     (a)
BPMnoise        = [10 40 90]*1e-6;              % metre

FocusingError   = 0.01;                         % % tesla   (dg)
AlignmentError  = 10e-6;                        % metre     (dy)

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
            RMScorrBPMTraj(n,j) = rms(correctorsTraj + randn(numel(beamline.bpmlist),1)*BPMnoise(j));
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
ylabel('RMS BPM values [mm]')
legend('Trajectory with correctors', 'BPMnoise = 10 \mum','BPMnoise = 40 \mum','BPMnoise = 90 \mum', 'Location', 'northwest')
%plot(meanRMScorrectorsTraj, 'DisplayName', 'Trajectory with correctors')

%%
rng(142)
lt.setQuadFerrors(true, FocusingError);
lt.setQuadAerrors(true, AlignmentError);

numObs = 60;
dcorrs = [1e-6 1e-5 1e-4 1e-3];

originalTraj     = lt.track_getBPMreadings();
correctorsTraj   = zeros(numel(dcorrs),numel(beamline.bpmlist));
correctorsTraj2  = zeros(numObs,numel(beamline.bpmlist));
correctorsTraj3  = zeros(numObs,numel(beamline.bpmlist));
correctorsTraj4  = zeros(numObs,numel(beamline.bpmlist));
dy               = zeros(numel(dcorrs),numel(beamline.bpmlist));


for i=1:numel(dcorrs)
    rng(1)
    dcorr       =  dcorrs(i)*randn(numel(beamline.corrlist),1);
    
    correctorsTraj(i,:)  = lt.getBPMvalues_corr(dcorr, 0);
    dy(i,:) = lt.getBPMorbitdy(dcorr);
end

for i=1:numObs
   rng(i)
   dcorr      =  dcorrs(1)*randn(numel(beamline.corrlist),1);
   correctorsTraj2(i,:)  = lt.getBPMvalues_corr(dcorr, 0);
   dcorr      =  dcorrs(2)*randn(numel(beamline.corrlist),1);
   correctorsTraj3(i,:)  = lt.getBPMvalues_corr(dcorr, 0);
      dcorr      =  dcorrs(3)*randn(numel(beamline.corrlist),1);
   correctorsTraj4(i,:)  = lt.getBPMvalues_corr(dcorr, 0);
end


figure(1)
subplot(2,1,2)
hold off
plot(dy(3,:)*1e3,'-r')
hold on
plot(dy(4,:)*1e3,'-b')
legend('100\mu{}T','1 mT')
xlabel('BPM index')
ylabel('\Deltay [mm]')

subplot(2,1,1)
hold off
plot(dy(1,:)*1e3, '-k')
hold on
plot(dy(2,:)*1e3,'-m')
plot(dy(3,:)*1e3,'-r')
legend('1 \mu{}T','10 \mu{}T','100 \mu{}T')
xlabel('BPM index')
ylabel('\Deltay [mm]')

figure(2)
subplot(3,1,1)
hold off
plot(correctorsTraj2'*1e3,'-k')
xlabel('BPM index')
ylabel('y [mm]')
title('a = 1 \mu{}T')
subplot(3,1,2)
hold off
plot(correctorsTraj3'*1e3,'-k')
xlabel('BPM index')
ylabel('y [mm]')
title('a = 10 \mu{}T')
subplot(3,1,3)
hold off
plot(correctorsTraj4'*1e3,'-k')
xlabel('BPM index')
ylabel('y [mm]')
title('a = 100 \mu{}T')

%% Old test 
% kept for reference - do not use

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
% 
% for i=1:numel(dcorrStrength)
%     for n=1:numSeeds
%         rng(n*i)
%         dcorrs       =  dcorrStrength(i)*randn(numel(beamline.corrlist),1);
% 
%         correctorsTraj  = lt.getBPMvalues_corr(dcorrs, 0);
%         RMScorrectorsTraj(n)  = rms(correctorsTraj);
% 
%         for j=1:(numel(BPMnoise))
%             RMScorrBPMTraj(n,j) = rms((correctorsTraj + randn(numel(beamline.bpmlist),1)*BPMnoise(j)));
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
% ylabel('\Deltay_{rms} (BPM values) [mm]')
% legend('Trajectory with correctors', 'BPMnoise = 10 \mum','BPMnoise = 40 \mum','BPMnoise = 90 \mum', 'Location', 'northwest')
% 
% f0 = fit((dcorrStrength*1e6)',meanRMScorrectorsTraj*1e3, 'poly3');
% f1 = fit((dcorrStrength*1e6)',meanRMScorrBPMTraj(:,1)*1e3, 'poly3');
% f2 = fit((dcorrStrength*1e6)',meanRMScorrBPMTraj(:,2)*1e3, 'poly3');
% f3 = fit((dcorrStrength*1e6)',meanRMScorrBPMTraj(:,3)*1e3, 'poly3');
% 
% figure(2)
% hold off
% plot(dcorrStrength*1e6,meanRMScorrectorsTraj*1e3, '-+k','DisplayName', 'Trajectory with correctors')
% hold on
% plot(f1, '-r')
% plot(f2, '--k')
% plot(f3, '-b')
% 
% 
% %plot(meanRMScorrectorsTraj, 'DisplayName', 'Trajectory with correctors')
% 
