% C. Monaghan
% 27/03/2021
% Fixed parameter scan for BPM test
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
FocusingError   = 0.1;

% Varying
ML_numobs       = linspace(100,1000,2);
BPMnoise        = linspace(40,90,6)*1e-6;

kickStrength    = linspace(10,90,3);   % Tesla

%% Test

STDresiduals_BPMnoise_C = zeros(length(BPMnoise), 1, numIterations);
STDresiduals_BPMnoise_ML = zeros(length(BPMnoise), length(ML_numobs), numIterations);
avgSTD_BPMnoise_ML = zeros(length(BPMnoise),numel(ML_numobs),numel(kickStrength));
avgSTD_BPMnoise_C  = zeros(length(BPMnoise), 1,numel(kickStrength));
RespMatML = zeros(numel(beamline.bpmlist), numel(beamline.corrlist), numel(ML_numobs));

for k = 1:numel(kickStrength)
    for i = 1:numel(BPMnoise)
        % Calc resp mat with conventional method in 'perfect' model
        RespMatC = lt.calcRespMatC(dcorrStrength);

        kickVals        = randn(numel(beamline.quaderrlist),1)*kickStrength(k);
        % Turn on errors
        lt.setQuadFerrors(true, FocusingError);
        lt.setQuadAerrors(true, kickVals);

        % Calc resp mat with ML in model with errors
        for j = 1:numel(ML_numobs)
            % Gen response matrix with ML method
            RespMatML(:,:,j) = lt.calcRespMatML(ML_numobs(j), BPMnoise(i), dcorrStrength);
        end

        for n = 1:numIterations
            bpmorbitdy = lt.track_getBPMreadings();

            corrected_C = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmorbitdy), BPMnoise(i));
            STDresiduals_BPMnoise_C(i,1,n)  = std(corrected_C);

            for j = 1:numel(ML_numobs)
                corrected_ML = lt.getBPMvalues_corr(pinv(RespMatML(:,:,j),  1e-3)*(-bpmorbitdy), BPMnoise(i));
                STDresiduals_BPMnoise_ML(i,j,n) = std(corrected_ML);
            end
        end

        % Turn off errors
        lt.setQuadFerrors(false, FocusingError);
        lt.setQuadAerrors(false, kickStrength);
    end

    for i=1:numel(BPMnoise)
        avgSTD_BPMnoise_C(i,1,k) = mean(STDresiduals_BPMnoise_C(i,1,:));
        for j=1:numel(ML_numobs)
            avgSTD_BPMnoise_ML(i,j,k) = mean(STDresiduals_BPMnoise_ML(i,j,:));
        end
    end
end

%% Construct plot

x_t = [];
x_b = [BPMnoise BPMnoise BPMnoise];
for i=1:numel(kickStrength), x_t = [x_t ones(numel(BPMnoise))*kickStrength(i)]; end

y_C    = [];
y_ML_1 = [];
y_ML_2 = [];
for k=1:numel(kickStrength)
    y_C    = [y_C avgSTD_BPMnoise_C(:,:,k)]; 
    y_ML_1 = [y_ML_1 avgSTD_BPMnoise_ML(:,j,1)]; 
    y_ML_2 = [y_ML_2 avgSTD_BPMnoise_ML(:,j,2)]; 
end

figure(1)

ax1 = axes();
plot(y_C(:), '--+r', 'DisplayName', 'Conventional')
plot(x_b,y_ML_1(:), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(1)))
plot(x_b,y_ML_2(:), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(2)))
ax2 = axes();
plot(x_t,y_C(:), '--+r', 'DisplayName', 'Conventional')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'left';
ax1.Box = 'off';
ax2.Box = 'off';
