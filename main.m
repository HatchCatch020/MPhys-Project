   % Tidy the workspace
%clearvars % -except beamline

% Define the beamline
DefineBeamline

% Define a strucuture to store all the variables relating to the beamline.
% This makes it easier to pass all the variables when calling a function.
% Could possibly add more variables to this as and when they are needed.
beamline.bl   = bl;    
beamline.corr = corr;
beamline.bpm  = bpm;
beamline.beam = beam;
beamline.quadsf = quadflist;
beamline.quadsd = quaddlist;

% Apply errors to the quadrupole magnets
quaderror = 0.01;
  
for n = 1:ncells
   
    fgradient = quadflist{n}.gradient;
    quadflist{n}.gradient = fgradient * (1 + quaderror*randn);
    
    dgradient = quaddlist{n}.gradient;
    quaddlist{n}.gradient = dgradient * (1 + quaderror*randn);
    
end

% Set the rf frequency
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

% Start by recording the closed orbit with all correctors off
zerocorrs        = zeros(numel(beamline.corr), 1);
bpmorbity_nocorr = getBPMorbitdy(beamline, zerocorrs);

% Specify the number of obersvations i.e. the number of columns in the
% arrays of corrector settings and BPM readings
numobs = 125;
% Initialize the arrays
dbpm   = zeros(numel(beamline.bpm),  numobs);
dcorra = zeros(numel(beamline.corr), numobs);

%% Main loop to generate an input for mvregress
dcorrstrength = 1e-4;    % Set the scale of corrector field strengths
bpmresn       = 70.0e-6; % Set the resolution of the bpms

for i = 1:numobs
    
    % Generate a random set of changes to corrector strengths...
    dcorr       =  dcorrstrength*randn(numel(beamline.corr),1);
    
    % ...and find the change in the closed orbit with these strengths
    bpmorbitdy  = getBPMorbitdy(beamline, dcorr) + randn(numel(beamline.bpm),1)*bpmresn;
    
    % Record the observation
    dbpm(:,i)   = bpmorbitdy;
    dcorra(:,i) = dcorr;
    
end

%% Analyse the data using mvregress
xmat  = dcorra';
ymat  = dbpm';

xcell = cell(1,numobs);

for i = 1:numobs
    xcell{i} = [kron([xmat(i,:)],eye(numel(beamline.bpm)))];
end

% Fit a response matrix to the observed changes in trajectory resulting
% from given changes in corrector strengths
[beta,sigma,E,V] = mvregress(xcell,ymat);

% Calculate the error on the fit
se = sqrt(diag(V));

% beta is the response matrix found from the machine learning approach
respmatML = reshape(beta, [numel(beamline.bpm), numel(beamline.corr)]);

%% Check the result

% Generate a random set of corrector changes
dcorr      =  dcorrstrength*randn(numel(beamline.corr),1);

% ...and find the resulting change in the closed orbit
bpmorbitdy = getBPMorbitdy(beamline, dcorr);

% See what we get from the machine learning analysis...
bpmorbitdyML = respmatML*dcorr;

% Calculate residual of ML and machine model
residual = (bpmorbitdyML - bpmorbitdy);

% Plot a comparison of the fitted vs actual response matrix elements
figure(1)
subplot(2,1,1)
hold off
plot(1e3*bpmorbitdy,'-ok')
hold on
plot(1e3*bpmorbitdyML,'--+r')
xlabel('BPM index')
ylabel('Vertical orbit (mm)')
legend('Model','Machine learning')
title('Change in closed orbit from random correctors')

subplot(2,1,2)
hold off
plot(1e3*residual,'-ok')
xlabel('BPM index')
ylabel('Residual (mm)')
title('Difference between ML prediction and machine model')

std(residual)

%% Functions
function respmatML=calcRespMatML(numobs_, BPMnoise_, dcorrstrengthScale, beamline_)
    dcorrstrength = dcorrstrengthScale;    % Set the scale of corrector field strengths
    bpmresn       = BPMnoise_; % Set the resolution of the bpms
    
    dbpm   = zeros(numel(beamline_.bpm),  numobs_);
    dcorra = zeros(numel(beamline_.corr), numobs_);

    for i = 1:numobs_

        % Generate a random set of changes to corrector strengths...
        dcorr       =  dcorrstrength*randn(numel(beamline_.corr),1);

        % ...and find the change in the closed orbit with these strengths
        bpmorbitdy  = getBPMorbitdy(beamline_, dcorr) + randn(numel(beamline_.bpm),1)*bpmresn;

        % Record the observation
        dbpm(:,i)   = bpmorbitdy;
        dcorra(:,i) = dcorr;

    end
    
    xmat  = dcorra';
    ymat  = dbpm';

    xcell = cell(1,numobs_);

    for i = 1:numobs_
        xcell{i} = [kron([xmat(i,:)],eye(numel(beamline_.bpm)))];
    end

    % Fit a response matrix to the observed changes in trajectory resulting
    % from given changes in corrector strengths
    [beta,sigma,E,V] = mvregress(xcell,ymat);

    % Calculate the error on the fit
    se = sqrt(diag(V));

    % beta is the response matrix found from the machine learning approach
    respmatML = reshape(beta, [numel(beamline_.bpm), numel(beamline_.corr)]);
end

function respmatC=calcRespMatC(beamline_, dcorrfieldScale)
    for j = 1:numel(beamline_.bpm)
        beamline_.bpm{j}.ResetBuffer(1);
    end

    for i = 1:numel(beamline_.corr)
        beamline.corr{i}.field = [0, 0];
    end
    
    closedorbit_nocorr = ComputeClosedOrbit(beamline_.bl, beamline_.beam);
    
    beamline_.beam.particles = closedorbit_nocorr(:,1);
    beamline_.bl.Track([1, numel(beamline_.bl.componentlist)], beamline_.beam);
    
    bpmorbity_nocorr = zeros(numel(beamline_.bpm), 1);
    
    for j=1:numel(beamline_.bpm)
       bpmorbity_nocorr(j) =  beamline_.bpm{j}.buffer(2);
    end
    
    respmatC = zeros(numel(beamline_.bpm), numel(beamline_.corr));
    
    for i=1:numel( beamline_.corr)

        dcorrfield = dcorrfieldScale;

        beamline_.corr{i}.field = [dcorrfield, 0];

        closedorbit = ComputeClosedOrbit( beamline_.bl,  beamline_.beam);

        for j = 1:numel( beamline_.bpm)
             beamline_.bpm{j}.ResetBuffer(1);
        end

         beamline_.beam.particles = closedorbit(:,1);
         beamline_.bl.Track([1, numel( beamline_.bl.componentlist)], beamline_.beam);

        for j=1:numel( beamline_.bpm)
            dy =  beamline_.bpm{j}.buffer(2) -  bpmorbity_nocorr(j);
            respmatC(j,i) = dy/dcorrfield;
        end
        
        beamline_.corr{i}.field = [0, 0];

    end
    
end

% This function makes for a tidy way to get the change in y BPM readings for
% a set of changes in corrector magnet strengths (corrStrenghth).
% In the future this function can be extended to get the change in x BPM
% readings as well. 
function bpmorbitdy_ = getBPMorbitdy(beamline_, dcorrStrength_)

    bpmorbitdy = zeros(numel(beamline_.bpm),1);

    % Find the closed orbit at the start
    closedorbit = ComputeClosedOrbit(beamline_.bl, beamline_.beam);
    
    for j = 1:numel(beamline_.bpm)
        beamline_.bpm{j}.ResetBuffer(1);
    end
    
    beamline_.beam.particles = closedorbit(:,1);
    beamline_.bl.Track([1, numel(beamline_.bl.componentlist)], beamline_.beam);

    for j = 1:numel(beamline_.bpm)
        bpmorbitdy(j) = beamline_.bpm{j}.buffer(2);
    end
    
    % Now apply the changes in corrector strength
    for i = 1:numel(beamline_.corr)
       beamline_.corr{i}.field = beamline_.corr{i}.field + [dcorrStrength_(i), 0];
    end
    
    % Find the change in closed orbit resulting from the change in
    % corrector strengths
    closedorbit = ComputeClosedOrbit(beamline_.bl,beamline_.beam);
    
    for j = 1:numel(beamline_.bpm)
        beamline_.bpm{j}.ResetBuffer(1);
    end
    
    beamline_.beam.particles = closedorbit(:,1);
    beamline_.bl.Track([1, numel(beamline_.bl.componentlist)], beamline_.beam);

    for j = 1:numel(beamline_.bpm)
        bpmorbitdy(j) = beamline_.bpm{j}.buffer(2) - bpmorbitdy(j);
    end
    
    % Put the corrector strengths back where they were
    for i = 1:numel(beamline_.corr)
       beamline_.corr{i}.field = beamline_.corr{i}.field - [dcorrStrength_(i), 0];
    end
    
    % Finally, return the result
    bpmorbitdy_ = bpmorbitdy;
    
end

