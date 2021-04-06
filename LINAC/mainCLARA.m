%% Define the beamline
DefineCLARABeamline
clearvars -except beamline driftlist quadflist quadflist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
% DefineCLARABeamline
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);


%% Define variables
ML_numobs       = 100;
BPMnoise        = 70.0e-6;
dcorrstrength   = 1e-4;
quaderror       = 1e-3;

global MLalgorithm;
MLalgorithm = 'cwls';

%% Gen response matrices
% Generate the response matrix using the conventional method with an ideal
% beamline that includes no errors
RespMatC = calcRespMatC(dcorrstrength, beamline);

% Set quad errors for the ML respmat 
setQuadErrors(beamline, quaderror);

% Gen response matrix with ML method
RespMatML = calcRespMatML(ML_numobs, BPMnoise, dcorrstrength, beamline);

%% Compare Response Matrices
figure(1)
subplot(2,1,1)
hold off
plot(RespMatC(:),RespMatML(:),'.k')
hold on
% plot(RespMatML(:),'--+r')
xlabel('Conventional Resp Mat [m/T]')
ylabel('ML Resp Mat [m/T]')
% legend('Conventional','Machine learning')
title('Values of response matrix elements')

subplot(2,1,2)
hold off
hist(RespMatML(:) - RespMatC(:),25)
xlabel('Residual [m/T]')
ylabel('Frequency')
title('Difference between ML method and Conventional method')

%% Check the result for ML
% Generate a random set of corrector changes
dcorr      =  dcorrstrength*randn(numel(beamline.corrlist),1);

% ...and find the resulting change in the closed orbit
bpmorbitdy = getBPMorbitdy(beamline, dcorr);

% See what we get from the machine learning analysis...
bpmorbitdyML = RespMatML*dcorr;

% Calculate residual of ML and machine model
residual = (bpmorbitdyML - bpmorbitdy);

% Plot a comparison of the fitted vs actual response matrix elements
figure(2)
subplot(2,1,1)
hold off
plot(1e3*bpmorbitdy,'-ok')
hold on
plot(1e3*bpmorbitdyML,'--+r')
xlabel('BPM index')
ylabel('Vertical position (mm)')
legend('Model','Machine learning')
title('Change in vertical position from random correctors')

subplot(2,1,2)
hold off
plot(1e3*residual,'-ok')
xlabel('BPM index')
ylabel('Residual (mm)')
title('Difference between ML prediction and machine model')

%%

% Generate a random set of corrector changes
dcorr      =  dcorrstrength*randn(numel(beamline.corrlist),1);

% ...and find the resulting change in the closed orbit
bpmorbitdy = getBPMorbitdy(beamline, dcorr);

dcorr_ML = pinv(RespMatML)*(bpmorbitdy);
dcorr_c = pinv(RespMatC)*(bpmorbitdy);

bpmorbitdy_correctedML = getBPMorbitdy(beamline, pinv(RespMatML, 1e-3)*(bpmorbitdy))+ randn(numel(beamline.bpmlist),1)*BPMnoise;
bpmorbitdy_correctedC = getBPMorbitdy(beamline, pinv(RespMatC, 1e-3)*(bpmorbitdy))+ randn(numel(beamline.bpmlist),1)*BPMnoise;

figure(3)
subplot(2,1,1)
hold off
plot(1e3*bpmorbitdy,'-ok')
hold on
plot(1e3*bpmorbitdy_correctedML,'--+r')
plot(1e3*bpmorbitdy_correctedC,'-.xb')
xlabel('BPM index')
ylabel('Vertical position (mm)')
legend('Original','ML corrected', 'Conventionally corrected')
title('Change in vertical position from random correctors and corrected positions')

subplot(2,1,2)
hold off
plot(1e3*(bpmorbitdy_correctedML - bpmorbitdy),'--+r')
hold on
plot(1e3*(bpmorbitdy_correctedC - bpmorbitdy),'-.xb')
xlabel('BPM index')
ylabel('Residual (mm)')
legend('ML corrected', 'Conventionally corrected')
title('Difference between corrected and uncorrected positions')

fprintf('std ML = %6.5g \n', std(bpmorbitdy_correctedML - bpmorbitdy));
fprintf('std C = %6.5g \n', std(bpmorbitdy_correctedC - bpmorbitdy));

figure(4)
subplot(2,1,1)
hold off
plot(dcorr,'-ok')
hold on
plot(dcorr_ML,'--+r')
plot(dcorr_c,'-.xb')
xlabel('Corrector index')
ylabel('Corrector Strenght ()')
legend('Original','ML predication', 'Conventional prediction')
title('Corrector magnet values')

subplot(2,1,2)
hold off
plot((dcorr_ML - dcorr),'--+r')
hold on
plot((dcorr_c - dcorr),'-.xb')
xlabel('Corrector index')
ylabel('Residual ()')
legend('ML predication', 'Conventional prediction')
title('Residual difference of corrector strengths for each RM')

%% Functions

function respmatML = calcRespMatML(numobs_, BPMnoise_, dcorrstrengthScale, beamline_)
    last_dbpm = csvread("last_dbpm.csv");
    last_dcorra = csvread("last_dcorra.csv");
    global MLalgorithm;
    dcorrstrength = dcorrstrengthScale;    % Set the scale of corrector field strengths
    bpmresn       = BPMnoise_; % Set the resolution of the bpms
    
    dbpm   = zeros(numel(beamline_.bpmlist),  numobs_);
    dcorra = zeros(numel(beamline_.corrlist), numobs_);
    
    % Generate input data for ML
    for i = 1:numobs_

        % Generate a random set of changes to corrector strengths...
        dcorr       =  dcorrstrength*randn(numel(beamline_.corrlist),1);

        % ...and find the change in the closed orbit with these strengths
        bpmorbitdy  = getBPMorbitdy(beamline_, dcorr) + randn(numel(beamline_.bpmlist),1)*bpmresn;

%         for n = 1:numel(beamline_.corrlist)
%            dcorr(n) = dcorr(n)*(beamline_.corrlist{n}.length/beamline_.beam.rigidity); 
%         end
        
        % Record the observation
        dbpm(:,i)   = bpmorbitdy;
        dcorra(:,i) = dcorr;

    end
    
    dbpm    = [dbpm last_dbpm];
    dcorra  = [dcorra last_dcorra];
    
    xmat  = dcorra';
    ymat  = dbpm';
    
    xcell = cell(1,length(dbpm));

    for i = 1:length(dbpm)
        xcell{i} = [kron([xmat(i,:)],eye(numel(beamline_.bpmlist)))];
    end

%     xcell = cell(1,numobs_);
% 
%     for i = 1:numobs_
%         xcell{i} = [kron([xmat(i,:)],eye(numel(beamline_.bpmlist)))];
%     end

    % Fit a response matrix to the observed changes in trajectory resulting
    % from given changes in corrector strengths
    [beta,sigma,E,V] = mvregress(xcell,ymat, 'algorithm', MLalgorithm);

    % Calculate the error on the fit
    se = sqrt(diag(V));

    % beta is the response matrix found from the machine learning approach
    respmatML = reshape(beta, [numel(beamline_.bpmlist), numel(beamline_.corrlist)]);
    
    dlmwrite('last_dbpm.csv', dbpm ,'delimiter',',');
    dlmwrite('last_dcorra.csv', dcorra ,'delimiter',',');
    
end

function respmatC = calcRespMatC(dcorrfieldScale, beamline_)
    % Make sure all kickers are set to zero
    for i = 1:numel(beamline_.corrlist)
        beamline_.corrlist{i}.field = [0, 0];
    end
    
    % Get the bpm readings with no correctors applied
    beamline_.beam.momentum = beamline_.momentum0;
    beamline_.beam.particles = [0 0 0 0 0 0]';
    bpmdy_nocorr = track_getBPMreadings(beamline_);
    
    % Prepare the matrix for the response matrix
    respmatC = zeros(numel(beamline_.bpmlist), numel(beamline_.corrlist));
    
    for i = 1:numel(beamline_.corrlist)
        beamline_.corrlist{i}.field = [dcorrfieldScale, 0];
        
        % Get the bpm readings with the new correctors applied
        beamline_.beam.momentum = beamline_.momentum0;
        beamline_.beam.particles = [0 0 0 0 0 0]';
        bpmdy = track_getBPMreadings(beamline_) - bpmdy_nocorr;
        
        for j = 1:numel(beamline_.bpmlist)
            % Calculate the components of the response matrix
            respmatC(j,i) = (bpmdy(j)/dcorrfieldScale);%*(beamline_.beam.rigidity/beamline_.corrlist{i}.length);
        end
        
        % Turn off this corrector
        beamline_.corrlist{i}.field = [0, 0];
    end
end

function bpmorbitdy_ = getBPMorbitdy(beamline_, dcorrStrength_)
    beamline_.beam.momentum = beamline_.momentum0;
    % Create a particle
    beamline_.beam.particles = [0 0 0 0 0 0]';
    
    % Track the particle and retreive the BPM readings 
    bpmorbitdy = track_getBPMreadings(beamline_);
    
    % Apply changes in correctors
    for i = 1:numel(beamline_.corrlist)
       beamline_.corrlist{i}.field = beamline_.corrlist{i}.field + [dcorrStrength_(i), 0];
    end
    
    % Reset particle to the beginning of beamline
    beamline_.beam.momentum = beamline_.momentum0;
    beamline_.beam.particles = [0 0 0 0 0 0]';
    
    % Get the new BPM readings and calculate the change
    bpmorbitdy_ = track_getBPMreadings(beamline_) - bpmorbitdy;
    
    % Return corrector strengths to previous value
    for i = 1:numel(beamline_.corrlist)
       beamline_.corrlist{i}.field = beamline_.corrlist{i}.field - [dcorrStrength_(i), 0];
    end
    
end

function bpmreadings_ = track_getBPMreadings(beamline_)
    bpmorbitdy = zeros(numel(beamline_.bpmlist), 1);
    
    % Reset the BPM buffers
    for j = 1:numel(beamline_.bpmlist)
        beamline_.bpmlist{j}.ResetBuffer(1);
    end

    % Track the particle through the beamline
    for n = 1:length(beamline_.bl.componentlist)
        beamline_.bl.Track([n n],beamline_.beam);
    end
    
    % Retreive the BPM readings
    for j = 1:numel(beamline_.bpmlist)
        bpmorbitdy(j) = beamline_.bpmlist{j}.buffer(2);
    end
    
    bpmreadings_ = bpmorbitdy;
end

function setQuadErrors(beamline_, quadError_)
    for n = 1:numel(beamline_.quadflist)
        fgradient = beamline_.quadflist{n}.gradient;
        beamline_.quadflist{n}.gradient = fgradient * (1 + quadError_*randn);
    end
    for n = 1:numel(beamline_.quaddlist)
        dgradient = beamline_.quaddlist{n}.gradient;
        beamline_.quaddlist{n}.gradient = dgradient * (1 + quadError_*randn);
    end
end


