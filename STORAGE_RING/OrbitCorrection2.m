ConstructResponseMatrix

% Specify the number of observations
numobs  =  50;

% Specify the number of BPMs
d       =   numel(bpm);

% Specify the number of trajectory correctors
ncorr   =   numel(corr);

% Generate a random response matrix
respmat = respmatrix;

% Generate a set of random changes to corrector strengths
dcorrstrength = 1e-5;
dcorr   = dcorrstrength*randn(ncorr/2,numobs);

% Calculate the corresponding changes in BPM readings (with some errors)
relnoise = 0.9;
dbpm    = respmat*dcorr + relnoise*dcorrstrength*randn(d,numobs);

% Set up the data structures for the linear regression
xmat    = dcorr';
ymat    = dbpm';

xcell   = cell(1,numobs);
for i = 1:numobs
    xcell{i} = [kron([xmat(i,:)],eye(d))];
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
