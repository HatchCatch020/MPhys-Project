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

% Specify some prperties 
numobs = 100;            % Specify the number of obersvations
dcorrstrength = 1e-4;   % The strength of the correctors
bpmresn = 1.0e-5;      % The resolution of the BPMs
quaderror = 0.001;      % The error on the quads  

% Set the rf frequency
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Generate Response Matrices
RespMatC = calcRespMatC(dcorrstrength, beamline);

% Now apply errors to the quadrupole magnets for the ML method
for n = 1:ncells
   
    fgradient = quadflist{n}.gradient;
    quadflist{n}.gradient = fgradient * (1 + quaderror*randn);
    
    dgradient = quaddlist{n}.gradient;
    quaddlist{n}.gradient = dgradient * (1 + quaderror*randn);
    
end

RespMatML = calcRespMatML(numobs, bpmresn, dcorrstrength, beamline);


%% Compare Response Matrices

figure(1)
subplot(2,1,1)
hold off
plot(RespMatC(:).*(beam.rigidity/corrlength),RespMatML(:).*(beam.rigidity/corrlength),'.k')
hold on
% plot(RespMatML(:),'--+r')
xlabel('Conventional Resp Mat [m/Rad]')
ylabel('ML Resp Mat [m/Rad]')
% legend('Conventional','Machine learning')
title('Values of response matrix elements')

subplot(2,1,2)
hold off
hist(RespMatML(:).*(beam.rigidity/corrlength) - RespMatC(:).*(beam.rigidity/corrlength),25)
xlabel('Residual [m/Rad]')
ylabel('Frequency')
title('Difference between ML method and Conventional method')
%% Check the result for ML

% Generate a random set of corrector changes
dcorr      =  dcorrstrength*randn(numel(beamline.corr),1);

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

%%

% Generate a random set of corrector changes
dcorr      =  dcorrstrength*randn(numel(beamline.corr),1);

% ...and find the resulting change in the closed orbit
bpmorbitdy = getBPMorbitdy(beamline, dcorr);

dcorr_correct = pinv(RespMatML)*(bpmorbitdy);

bpmorbitdy_correctedML = getBPMorbitdy(beamline, pinv(RespMatML)*(bpmorbitdy));
bpmorbitdy_correctedC = getBPMorbitdy(beamline, pinv(RespMatC)*(bpmorbitdy));

figure(3)
subplot(2,1,1)
hold off
plot(1e3*bpmorbitdy,'-ok')
hold on
plot(1e3*bpmorbitdy_correctedML,'--+r')
plot(1e3*bpmorbitdy_correctedC,'-.xb')
xlabel('BPM index')
ylabel('Vertical orbit (mm)')
legend('Original orbit','ML corrected orbit', 'Conventionally corrected orbit')
title('Change in closed orbit from random correctors and corrected orbit')

subplot(2,1,2)
hold off
plot(1e3*(bpmorbitdy_correctedML - bpmorbitdy),'--+r')
hold on
plot(1e3*(bpmorbitdy_correctedC - bpmorbitdy),'-.xb')
xlabel('BPM index')
ylabel('Residual (mm)')
legend('ML corrected orbit', 'Conventionally corrected orbit')
title('Difference between corrected and uncorrected orbit')

%% Output some information

ML  = std(bpmorbitdy_correctedML - bpmorbitdy) 
C   = std(bpmorbitdy_correctedC - bpmorbitdy)

dlmwrite('output.csv',[numobs dcorrstrength bpmresn quaderror std(bpmorbitdy_correctedML - bpmorbitdy) std(bpmorbitdy_correctedC - bpmorbitdy)],'delimiter',',','-append');

%% Functions
function respmatML=calcRespMatML(numobs_, BPMnoise_, dcorrstrengthScale, beamline_)
    dcorrstrength = dcorrstrengthScale;    % Set the scale of corrector field strengths
    bpmresn       = BPMnoise_; % Set the resolution of the bpms
    
    dbpm   = zeros(numel(beamline_.bpm),  numobs_);
    dcorra = zeros(numel(beamline_.corr), numobs_);
    
    % Generate input data for ML
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

% Generate the response matrix for a conventional method, this method turns
% on each corrector magnet in turn and then turns it off.
function respmatC=calcRespMatC(dcorrfieldScale, beamline_)
    for j = 1:numel(beamline_.bpm)
        beamline_.bpm{j}.ResetBuffer(1);
    end

    for i = 1:numel(beamline_.corr)
        beamline.corr{i}.field = [0, 0];
    end
    
    % compute the closed orbit without changes in corr strength
    closedorbit_nocorr = ComputeClosedOrbit(beamline_.bl, beamline_.beam);
    
    beamline_.beam.particles = closedorbit_nocorr(:,1);
    beamline_.bl.Track([1, numel(beamline_.bl.componentlist)], beamline_.beam);
    
    % Prepare the array BPM values without correctors
    bpmorbity_nocorr = zeros(numel(beamline_.bpm), 1);
    
    % fill the array with the values
    for j=1:numel(beamline_.bpm)
       bpmorbity_nocorr(j) =  beamline_.bpm{j}.buffer(2);
    end
    
    % Prepare the matrix for the response matrix
    respmatC = zeros(numel(beamline_.bpm), numel(beamline_.corr));
    
    % Main loop for applying the changes in corrector strength. Turning on
    % a different corr in each iteration.
    for i=1:numel( beamline_.corr)

        % Set the corr strength as the input of the function.
        dcorrfield = dcorrfieldScale;

        beamline_.corr{i}.field = [dcorrfield, 0];

        % Compute the closed orbit with this change in corrector strength 
        closedorbit = ComputeClosedOrbit( beamline_.bl,  beamline_.beam);

        for j = 1:numel( beamline_.bpm)
             beamline_.bpm{j}.ResetBuffer(1);
        end

         beamline_.beam.particles = closedorbit(:,1);
         beamline_.bl.Track([1, numel( beamline_.bl.componentlist)], beamline_.beam);

        for j=1:numel( beamline_.bpm)
            % Calculate the change in BPM as a result of this change in
            % corr strength.
            dy =  beamline_.bpm{j}.buffer(2) -  bpmorbity_nocorr(j);
            % Calculate this component of the response matrix 
            respmatC(j,i) = dy/dcorrfield;
        end
        
        % Turn off this corrector in preperation for the next itteration
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

