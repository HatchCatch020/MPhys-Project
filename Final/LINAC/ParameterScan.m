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

%% Variables

% Fixed variables
numSeeds   = 100;

% Varying
ML_numobs       = linspace(100,100,1);
BPMnoise        = linspace(10,60,6)*1e-6;

dcorrStrength   = 1e-5;%[1e-4 1e-5 1e-6];
AlignmentError  = linspace(10,50,3)*1e-6;   % m

FocusingError   = 0.01;%linspace(1,9,5)*1e-1;

%% Test

% Setup arrrays
rmsBPMvals_C        = zeros(numel(BPMnoise), 1, numSeeds);
rmsBPMvals_ML       = zeros(numel(BPMnoise), numel(ML_numobs), numSeeds);
mean_rmsBPMvals_C   = zeros(numel(BPMnoise), 1, numel(AlignmentError));
mean_rmsBPMvals_ML  = zeros(numel(BPMnoise), numel(ML_numobs), numel(AlignmentError));
err_rmsBPMvals_C    = zeros(numel(BPMnoise), 6, numel(AlignmentError));
err_rmsBPMvals_ML   = zeros(numel(BPMnoise), 6, numel(AlignmentError));

for k = 1:numel(AlignmentError)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = lt.calcRespMatC(dcorrStrength);
    for i = 1:numel(BPMnoise)
        for n = 1:numSeeds
            % set seed and turn on a set of errors
            rng(n*i*k)
            lt.setQuadFerrors(true, FocusingError);
            lt.setQuadAerrors(true, AlignmentError(k));
            
            bpmValsY = lt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise(i);

            bpmValsCorrected_C  = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), BPMnoise(i));
            rmsBPMvals_C(i,1,n) = rms(bpmValsCorrected_C);

            for j = 1:numel(ML_numobs)
                RespMatML = lt.calcRespMatML(ML_numobs(j), BPMnoise(i), dcorrStrength);
                bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), BPMnoise(i));
                rmsBPMvals_ML(i,j,n) = rms(bpmValsCorrected_ML);
            end
            % Turn off errors
            lt.setQuadFerrors(false, FocusingError);
            lt.setQuadAerrors(false, AlignmentError(k));
        end % n (numSeeds)
        % Calc the mean of rms values over all the seeds, and also the std of
        % the rms values.
        mean_rmsBPMvals_C(i,k)  =  mean(rmsBPMvals_C(i,1,:));
        err_rmsBPMvals_C(i,:,k)   =   prctile(rmsBPMvals_C(i,1,:), [10 90 20 80 30 70]);
        for j=1:numel(ML_numobs) 
            mean_rmsBPMvals_ML(i,j,k) = mean(rmsBPMvals_ML(i,j,:));
            err_rmsBPMvals_ML(i,:,k)  =  prctile(rmsBPMvals_ML(i,j,:), [10 90 20 80 30 70]);
        end % j 
    end % i (BPMnoise) 
end % k (kickStrength)

%% Construct plot

% x_t = [];
% x_b = [AlignmentError AlignmentError AlignmentError];
% for i=1:numel(dcorrStrength), x_t = [x_t ones(numel(AlignmentError))*dcorrStrength(i)]; end
% 
% y_ML_1 = [];
% y_ML_2 = [];
% y_ML_1err = [];
% %y_ML_2err = [];
% y_C     = mean_rmsBPMvals_C(:);
% y_Cerr  =  errp_rmsBPMvals_C(:);
% for k=1:numel(dcorrStrength)
%     y_ML_1 = [y_ML_1 mean_rmsBPMvals_ML(:,1,k)]; 
%     %y_ML_2 = [y_ML_2 mean_rmsBPMvals_ML(:,2,k)];
%     y_ML_1err  =  [y_ML_1err errp_rmsBPMvals_ML(:,1,k)];
%     %y_ML_2err  =  [y_ML_2err std_rmsBPMvals_ML(:,2,k)];
% end
% y_ML_1 = y_ML_1(:);
% %y_ML_2 = y_ML_2(:);
% y_ML_1err = y_ML_1err(:);
% %y_ML_2err = y_ML_2err(:);
% 
% output = [(x_b')*1e6 x_t(i,:)' y_C*1e3 y_ML_1*1e3  y_Cerr*1e3 y_ML_1err*1e3];

LtopVal = length(AlignmentError);
LbotVal = length(BPMnoise);
yConv = NaN(LtopVal,LbotVal*LtopVal+1);
yML = NaN(LtopVal,LbotVal*LtopVal+1,numel(ML_numobs));

for i=1:LtopVal
    yConv(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = mean_rmsBPMvals_C(:,i)';
    for j = 1:numel(ML_numobs)
        yML(i,(i*LbotVal)-LbotVal+1:i*LbotVal,j) = mean_rmsBPMvals_ML(:,j,i)';
    end
end
yConverrN = NaN(LtopVal,LbotVal*LtopVal+1);
yConverrP = NaN(LtopVal,LbotVal*LtopVal+1);
yMLerrN = NaN(LtopVal,LbotVal*LtopVal+1,numel(ML_numobs));
yMLerrP = NaN(LtopVal,LbotVal*LtopVal+1,numel(ML_numobs));

for i=1:LtopVal
    yConverrN(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = err_rmsBPMvals_C(:,5,i)';
    yConverrP(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = err_rmsBPMvals_C(:,6,i)';
    yMLerrN(i,(i*LbotVal)-LbotVal+1:i*LbotVal,j) = err_rmsBPMvals_ML(:,5,i)';
    yMLerrP(i,(i*LbotVal)-LbotVal+1:i*LbotVal,j) = err_rmsBPMvals_ML(:,6,i)';
end

labels_bot = strings(LtopVal*LbotVal,1);
labels_top = strings(LtopVal*LbotVal,1);
g=1;
for i=1:LtopVal
    for k=1:LbotVal
        labels_bot(g) = sprintf('%1g', BPMnoise(k)*1e6);
        labels_top(g) = sprintf('%1g', AlignmentError(i)*1e6);
        g = g+1;
    end
end
labels_bot(end) = [];

x = 0:LbotVal*LtopVal;

xT = 0:LbotVal*LtopVal-1;
xTL = 0:1:LbotVal*LtopVal;

% fill1   = [8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
% fill2   = [NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN];
% fill3   = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4];

fill1   = [8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
fill2   = [NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN];
fill3   = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4];

figure

ax2 = axes('Position',[0.1 0.15 0.85 0.72]);
    %plot(ax2,x,y_ML_1c*1e3, '-.xb', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(1)))
    hold on
    area(ax2,xTL, fill1*1e3,'FaceAlpha',0.5, 'EdgeColor', 'none')
    area(ax2,xTL, fill2*1e3,'FaceAlpha',0.5, 'EdgeColor', 'none')
    area(ax2,xTL, fill3*1e3,'FaceAlpha',0.5, 'EdgeColor', 'none')
    xlabel('Alignment Error [\mum]','LineWidth',3,'FontWeight','bold','FontSize',18)
    xlim([0 18])
    
ax1 = axes('Position',[0.1 0.15 0.85 0.72]);
    hold on
    eC(1) = errorbar(ax1,x,yConv(1,:)*1e3, yConverrN(i,:)*1e3, yConverrP(i,:)*1e3, '--+r');
    eML(1) = errorbar(ax1,x,yML(1,:)*1e3, yMLerrN(1,:)*1e3, yMLerrP(1,:)*1e3,':*k');
    for i=2:LtopVal
        eC(i) = errorbar(ax1,x,yConv(i,:)*1e3, yConverrN(i,:)*1e3, yConverrP(i,:)*1e3, '--+r');
    end

    for i=2:LtopVal
        eML(i) = errorbar(ax1,x,yML(i,:)*1e3, yMLerrN(i,:)*1e3, yMLerrP(i,:)*1e3, ':*k');
    end
    xlabel('BPM Noise [\mum]','LineWidth',3,'FontWeight','bold','FontSize',18)
    ylabel('\Deltay_{RMS} [mm]','LineWidth',3,'FontWeight','bold','FontSize',18)
    xlim([0 18])
    
box(ax1,'on');
ax2.Box = 'off';

legend('Conventional','ML (NumObs = 100)')

set(ax2, 'FontSize',18,'FontWeight','bold','LineWidth',3,'YColor','none', ...
    'XTick',xTL, 'XTickLabel', labels_top, 'XTickLabelRotation',45, ...
    'XAxisLocation', 'top','YAxisLocation','right')

set(ax1, 'FontSize',18,'FontWeight','bold','LineWidth',3,'Color','none', ...
    'XTick',x, 'XTickLabel', labels_bot)

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 5])
print('-dpng','fig.png','-r600')

%output = [(x_b')*1e6 x_t(i,:)' y_C*1e3 y_ML_1*1e3 y_ML_2*1e3 y_Cerr*1e3 y_ML_1err*1e3 y_ML_2err*1e3];

% figure(1)
% 
% ax1 = axes();
% plot(y_C(:), '--+r', 'DisplayName', 'Conventional')
% plot(x_b,y_ML_1(:), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(1)))
% plot(x_b,y_ML_2(:), '-.x', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(2)))
% ax2 = axes();
% plot(x_t,y_C(:), '--+r', 'DisplayName', 'Conventional')
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'left';
% ax1.Box = 'off';
% ax2.Box = 'off';
