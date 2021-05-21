% C. Monaghan
% CLARA_MLobsCheck
% Investigation of the scale of dcorrector strength used to calculate
% response matrices and the resulting change in the orbit from such
% corrector strengths.

%% Define the beamline
DefineCLARABeamline
clearvars -except beamline driftlist quadlist quaderrlist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Setup
% Load the functions which will be used in this script.
% Access the functions using the object lt
lt = linTools;
lt.beamline = beamline;
lt.ML_Algorithm = 'mvn';

%% Variables

numSeeds = 5;

MLnumobs = [52, 60, 70, 80, 90, 100, 1000];

dcorrStrength   = [1 10 100 1000]*1e-6;   % tesla   (a)
BPMnoise        = 10e-6;                  % metres 

FocusingError   = 0.01;                   % % tesla (dg)
AlignmentError  = 10e-6;                  % metre   (dy)

%%

figure(1)

for j=1:numel(dcorrStrength)
    rmsBPMvals_ML = zeros(numel(MLnumobs), numSeeds);
    mean_rmsBPMvals_ML = zeros(numel(MLnumobs), 1);
    err_rmsBPMvals_ML = zeros(numel(MLnumobs), 6);

    for i=1:numel(MLnumobs)
       for n=1:numSeeds
            rng(n*i)
            % Turn on errors
            lt.setQuadFerrors(true, FocusingError);
            lt.setQuadAerrors(true, AlignmentError);

            bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise;

            RespMatML = lt.calcRespMatML(MLnumobs(i), BPMnoise, dcorrStrength(j));

            bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), BPMnoise);
            rmsBPMvals_ML(i,n) = rms(bpmValsCorrected_ML);

            % Turn off errors
            lt.setQuadFerrors(false, FocusingError);
            lt.setQuadAerrors(false, AlignmentError);
       end
       mean_rmsBPMvals_ML(i) = mean(rmsBPMvals_ML(i,:));
       err_rmsBPMvals_ML(i,:)  =  prctile(rmsBPMvals_ML(i,:),[10 90 20 80 30 70]);
    end

    f0 = fit([1 2 3 4 5 6 7]',mean_rmsBPMvals_ML*1e6, 'poly1')

    subplot(2,2,j)
    hold off
    plot([1 2 3 4 5 6 7]', mean_rmsBPMvals_ML*1e6,'sk')
    hold on
    plot(f0, '--k')
    shade([1 2 3 4 5 6 7]',err_rmsBPMvals_ML(:,3)*1e6,[1 2 3 4 5 6 7]',err_rmsBPMvals_ML(:,4)*1e6,'FillType',[2 1],'LineStyle', 'none')
    
    legend('Linear fit', 'Location', 'NorthWest')
    xlim([0 8])
    xlabel('Number of observations')
    ylabel('RMS BPM values [\mum]')
    title(sprintf('a = %g \\mu{}T', dcorrStrength(j)*1e6))

    set(gca, 'XTickLabel', ['52','60','70','80','90','100','1000'])
    
end





