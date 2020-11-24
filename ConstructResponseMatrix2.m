
DefineBeamline

f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

% Define blmorbity_nocorr as the bpm readings with no corector magnets.

for j = 1:numel(bpm)
   bpm{j}.ResetBuffer(1);
end

for i = 1:numel(corr)
   corr{i}.field = [0, 0];
end

closedorbit_nocorr = ComputeClosedOrbit(bl, beam);

beam.particles = closedorbit_nocorr(:, 1);
bl.Track([1, numel(bl.componentlist)], beam);

bpmorbity_nocorr = zeros(numel(bpm), 1);

for j = 1:numel(bpm)
    bpmorbity_nocorr(j) = bpm{j}.buffer(2);
end

respmatrix = zeros(numel(bpm), numel(corr));
figure(1)
hold off
plot(bpmorbity_nocorr);


%dcorrfield = randn*5e-3;

for i=1:numel(corr)
    
    dcorrfield = randn*5e-3;
    
    corr{i}.field = [dcorrfield, 0];
    
    closedorbit = ComputeClosedOrbit(bl, beam);
    
    for j = 1:numel(bpm)
        bpm{j}.ResetBuffer(1);
    end
    
    beam.particles = closedorbit(:,1);
    bl.Track([1, numel(bl.componentlist)],beam);
    
    for j=1:numel(bpm)
        dy = bpm{j}.buffer(2) -  bpmorbity_nocorr(j);
        respmatrix(j,i) = dy/dcorrfield;
    end
    
    corr{i}.field = [0, 0];
    
end

% Test the response matrix %

% Set some random strengths for the correctors
corrfields = dcorrfield*randn(numel(corr),1);

for i=1:numel(corr)
    
    corr{i}.field = [corrfields(i), 0];
    
end

% Find the closed orbit with the random correctors
closedorbit = ComputeClosedOrbit(bl, beam);

for j = 1:numel(bpm)
    bpm{j}.ResetBuffer(1);
end

beam.particles = closedorbit(:,1);
bl.Track([1, numel(bl.componentlist)],beam);

dbpmorbit = zeros(numel(bpm),1);

for j = 1:numel(bpm)
    dbpmorbit(j) = bpm{j}.buffer(2) -  bpmorbity_nocorr(j);
end

predictedorbit = respmatrix * corrfields;

figure(1)
hold on
plot(dbpmorbit,'or-')
plot(bpmorbity_nocorr,'ob--')


% Find the predicted orbit change

%predictedorbit = respmatrix * corrfields;

figure(2)
subplot(2,1,1)
hold off
plot((dbpmorbit-predictedorbit)./abs(dbpmorbit),'ob-')
% hold on
% plot(predictedorbit,'xr-')

% Find the predicted corrector strengths

predictedcorrs = pinv(respmatrix) * dbpmorbit;

subplot(2,1,2)
hold off
plot((corrfields-predictedcorrs)./abs(corrfields),'ob-')

