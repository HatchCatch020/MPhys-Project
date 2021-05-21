% C. Monaghan
% CLARA_RMinvestigation
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
FocusingError   = 0.01;     % %tesla    (dg)
dcorrStrength   = 1e-2;     % tesla     (dy)
kickStrength    = 10e-6;    % metre     (a)
ML_numobs       = 100;      %
BPMnoise        = 10e-6;    % metre

%% Calc RMs
% An older method of comparison, keep for reference later.


% rng(3)
% RespMatC = lt.calcRespMatC(1e-6);
% RespMatC2 = lt.calcRespMatC(1e-2);
% 
% figure(1)
% subplot(2,1,1)
% hold off
% hist(RespMatC(:) - RespMatC2(:),50)
% xlabel('Residual [m/T]')
% ylabel('Frequency')
% title('Residual between RM calculated with conventional method (dcorr = 1e-6,1e-2 T)')
% 
% % Turn on errors
% lt.setQuadFerrors(true, FocusingError);
% lt.setQuadAerrors(true, kickStrength);
% 
% RespMatML = lt.calcRespMatML(ML_numobs, BPMnoise, 1e-6);
% RespMatML2 = lt.calcRespMatML(ML_numobs, BPMnoise, 1e-2);
% 
% subplot(2,1,2)
% hold off
% hist(RespMatML(:) - RespMatML2(:),50')
% xlabel('Residual [m/T]')
% ylabel('Frequency')
% title('Residual between RM calculated with ML method (dcorr = 1e-6,1e-2 T)')
% 
% 
% %% Compare Response Matrices
% figure(2)
% subplot(2,1,1)
% hold off
% plot(RespMatC2(:),RespMatML2(:),'.k')
% hold on
% xlabel('Conventional Resp Mat [m/T]')
% ylabel('ML Resp Mat [m/T]')
% title('Values of response matrix elements (dcorr = 1e-2 T)')
% 
% subplot(2,1,2)
% hold off
% hist(RespMatML2(:) - RespMatC2(:),25)
% xlabel('Residual [m/T]')
% ylabel('Frequency')
% title('Difference between ML method and conventional method (dcorr = 1e-2 T)')
% 
% figure(3)
% subplot(2,1,1)
% hold off
% plot(RespMatC(:),RespMatML(:),'.k')
% hold on
% xlabel('Conventional Resp Mat [m/T]')
% ylabel('ML Resp Mat [m/T]')
% title('Values of response matrix elements (dcorr = 1e-6 T)')
% 
% subplot(2,1,2)
% hold off
% hist(RespMatML(:) - RespMatC(:),25)
% xlabel('Residual [m/T]')
% ylabel('Frequency')
% title('Difference between ML method and conventional method (dcorr = 1e-6 T)')

%%
dcorr = [1 10 100 1000 10000 100000]*1e-6;%[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];

figure(4)
for i=1:numel(dcorr)
    RespMatC    = lt.calcRespMatC(dcorr(i));
    % Turn on errors
    lt.setQuadFerrors(true, FocusingError);
    lt.setQuadAerrors(true, kickStrength);
    RespMatML   = lt.calcRespMatML(ML_numobs, BPMnoise, dcorr(i));
    lt.setQuadFerrors(false, FocusingError);
    lt.setQuadAerrors(false, kickStrength);
    
    % gen plot
    subplot(2,3,i)
    plot(RespMatC(:),RespMatML(:), '.k')
    xlabel('Conventional Resp Mat [m/T]')
    ylabel('ML Resp Mat [m/T]')
    title(sprintf('a = %g \\mu{}T',dcorr(i)*1e6))
    pbaspect([1 1 1])
end

%%
dcorr = linspace(1,10,4)*1e-5;%[1e-5 1e-5 1e-5 1e-5 1e-5 10e-5];

figure(5)
for i=1:numel(dcorr)
    RespMatC    = lt.calcRespMatC(dcorr(i));
    % Turn on errors
    lt.setQuadFerrors(true, FocusingError);
    lt.setQuadAerrors(true, kickStrength);
    RespMatML   = lt.calcRespMatML(ML_numobs, BPMnoise, dcorr(i));
    lt.setQuadFerrors(false, FocusingError);
    lt.setQuadAerrors(false, kickStrength);
    
    % gen plot
    subplot(2,2,i)
    plot(RespMatC(:),RespMatML(:), '.k')
    xlabel('Conventional Resp Mat [m/T]')
    ylabel('ML Resp Mat [m/T]')
    title(sprintf('\\Deltacorr = %g T',dcorr(i)))
    pbaspect([1 1 1])
end

