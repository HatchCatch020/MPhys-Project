DefineBeamline

% Define a strucuture to store all the variables relating to
% DefineBeamline.m
% Could possible add more variables to this as and when they are needed.
beamline.bl = bl;       % This makes it easier to pass all the variables in 
beamline.beam = beam;   % on go when calling a function.
beamline.corr = corr;
beamline.bpm = bpm;

% Tidy the workspace of all the DefineBeamline variables
clearvars -except beamline

f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

% The number of obersvations i.e the number of columns in the array.
numobs = 50;

% % Test the getBPMorbity function. 
% figure(1)
% hold on
% plot(getBPMorbity(dcorr, beamline), '--b.')
% plot(getBPMorbity(zeros(numel(beamline.corr), 1), beamline), '-ro')

bpmorbity_nocorr = getBPMorbity(beamline, zeros(numel(beamline.corr), 1));

dbpm = zeros(numel(beamline.bpm), numobs);
dcorra = zeros(numel(beamline.corr), numobs);

% Main loop to generate an input for mvregress
for i = 1:numobs
    dcorrstrength = 1e-4;
    dcorr   =  dcorrstrength*randn(numel(beamline.corr),1);
    
    bpmorbity = getBPMorbity(beamline, dcorr);
    
    dbpm(:,i) = getBPMorbity(beamline, dcorr);
    dcorra(:,i) = dcorr;
end

%% mvregress testing

xmat    = dcorra';
ymat    = dbpm';

xcell   = cell(1,numobs);
for i = 1:numobs
    xcell{i} = [kron([xmat(i,:)],eye(numel(beamline.bpm)))];
end

% Fit a response matrix to the observed changes in trajectory resulting
% from given changes in corrector strengths
[beta,sigma,E,V] = mvregress(xcell,ymat);

% Calculate the error on the fit
se = sqrt(diag(V));

% beta is the response matrix found from the machine learning approach

% Plot a comparison of the fitted vs actual response matrix elements
errorbar(respmat(:),beta,se,'.k')
xlabel('Real value')
ylabel('Fitted value')
title('Response matrix elements')

%% Functions
% This function makes for a tidy way to get all the BPM readings for a set
% of corrector magnet strengths (corrStrenghth). In the future this single
% fucntion can be extended to get the BPM readings in x also. 
function bpmorbity_ = getBPMorbity(beamline_, corrStrength_)
    for j = 1:numel(beamline_.bpm)
        beamline_.bpm{j}.ResetBuffer(1);
    end
    for i = 1:numel(beamline_.corr)
       beamline_.corr{i}.field = [0, 0];
    end
    
    closedorbit = ComputeClosedOrbit(beamline_.bl,beamline_.beam);
    
    for i = 1:numel(beamline_.corr)
       beamline_.corr{i}.field = [corrStrength_(i), 0];
    end
    
    beamline_.beam.particles = closedorbit(:,1);
    beamline_.bl.Track([1, numel(beamline_.bl.componentlist)], beamline_.beam);
    bpmorbity = zeros(numel(beamline_.bpm), 1);
    
    for j = 1:numel(beamline_.bpm)
        bpmorbity(j) = beamline_.bpm{j}.buffer(2);
    end
    
    bpmorbity_ = bpmorbity;
    
end

