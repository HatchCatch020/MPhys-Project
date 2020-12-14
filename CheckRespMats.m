% Define the beamline
DefineBeamline

% Define a strucuture to store all the variables relating to the beamline.
% This makes it easier to pass all the variables when calling a function.
% Could possibly add more variables to this as and when they are needed.
beamline.bl     = bl;    
beamline.corr   = corr;
beamline.bpm    = bpm;
beamline.beam   = beam;
beamline.quadsf = quadflist;
beamline.quadsd = quaddlist;

% Specify some prperties 
numobs          = 50;       % Specify the number of obersvations
dcorrstrength   = 1e-4;     % The strength of the correctors
bpmresn         = 70.0e-6;  % The resolution of the BPMs
quaderror       = 0.00;     % The error on the quads
  
% Apply errors to the quadrupole magnets
for n = 1:ncells
   
    fgradient = quadflist{n}.gradient;
    quadflist{n}.gradient = fgradient * (1 + quaderror*randn);
    
    dgradient = quaddlist{n}.gradient;
    quaddlist{n}.gradient = dgradient * (1 + quaderror*randn);
    
end

% Set the rf frequency
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Main

% Check ML response matrix

% Get the response matrix from the Ml method, dcorr will also store the
% random set of corr changes that are used here.
[RespMatML,dcorr] = calcRespMatML(numobs, bpmresn, dcorrstrength, beamline);

% Get the bpm values again for this same set of cor changes
bpmorbitdy = getBPMorbitdy(beamline, dcorr);

% Calculate the dcorr values using the new bpm readings and the response
% matrix from the ML model
dcorrCalculated = pinv(RespMatML)*bpmorbitdy;

figure(1)
subplot(2,1,1)  % compare the calulated dcorr to the orignal one used
hold off
plot(dcorr, '-ob');
hold on
plot(dcorrCalculated, '-or');
legend("dcorr", "dcorrCalculated");
ylabel("Corrector Strength"); xlabel("Corrector index");

subplot(2,1,2)  % calculate the residual of the previous plot
hold off
plot(dcorrCalculated - dcorr, '-ok');
ylabel("Residual"); xlabel("Corrector index");

std(dcorrCalculated - dcorr)

%% Functions
function [respmatML,dcorr]=calcRespMatML(numobs_, BPMnoise_, dcorrstrengthScale, beamline_)
    dcorrstrength = dcorrstrengthScale; % Set the scale of corrector field strengths
    bpmresn       = BPMnoise_;          % Set the resolution of the bpms
    
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
function respmatC=calcRespMatC(dcorrfieldarray, beamline_)
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

        % Set the corr strength as the function input.
        dcorrfield = dcorrfieldarray(i);

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

