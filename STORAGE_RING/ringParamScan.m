% C. Monaghan
% 27/03/2021
% ringParamScan
%% Define the beamline
DefineBeamline
clearvars -except beamline driftlist quadflist quadflist corrlist bpmlist Lcavitylist Sbendlist beam bl

% Set the master oscillaor, currently this is as previosuly define in
% DefineCLARABeamline
f0 = MasterOscillator.GetFrequency();
MasterOscillator.SetFrequency(f0*1.0);

%% Setup
% Load the functions which will be used in this script.
% Access the functions using the object rt
rt = ringTools;
rt.beamline = beamline;
rt.ML_Algorithm = 'mvn';

%% Variables

% Fixed variables
numSeeds   = 50;

% Varying
ML_numobs       = 60;
BPMnoise        = linspace(10,70,5)*1e-6;       % metre

dcorrStrength   = [1 10 100]*1e-6;              % tesla     (a)
AlignmentError  = 10e-6;                        % metre     (dy)

FocusingError   = 0.01;%linspace(1,9,5)*1e-1;   % % tesla   (dg)

%% Test

% Setup arrrays
rmsBPMvals_C        = zeros(numel(BPMnoise), numSeeds);
rmsBPMvals_ML       = zeros(numel(BPMnoise), numSeeds);
mean_rmsBPMvals_C   = zeros(numel(BPMnoise), numel(dcorrStrength));
mean_rmsBPMvals_ML  = zeros(numel(BPMnoise), numel(dcorrStrength));
err_rmsBPMvals_C    = zeros(numel(BPMnoise), 6, numel(dcorrStrength));
err_rmsBPMvals_ML   = zeros(numel(BPMnoise), 6, numel(dcorrStrength));

for k = 1:numel(dcorrStrength) % top
    RespMatC = rt.calcRespMatC(dcorrStrength(k));
    % Calc resp mat with conventional method in 'perfect' model
    for i = 1:numel(BPMnoise) % bot
        for n = 1:numSeeds
            % set seed and turn on a set of errors
            rng(n*i*k)
            rt.setQuadFerrors(true, FocusingError);
            rt.setQuadAerrors(true, AlignmentError);
            
            bpmValsY = rt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise(i);

            bpmValsCorrected_C  = rt.getBPMvalues_corr(pinv(RespMatC,  10e-3)*(-bpmValsY), BPMnoise(i));
            rmsBPMvals_C(i,n) = rms(bpmValsCorrected_C);

            RespMatML = rt.calcRespMatML(ML_numobs, BPMnoise(i), dcorrStrength(k));
            bpmValsCorrected_ML  = rt.getBPMvalues_corr(pinv(RespMatML,  10e-3)*(-bpmValsY), BPMnoise(i));
            rmsBPMvals_ML(i,n) = rms(bpmValsCorrected_ML);
            % Turn off errors
            rt.setQuadFerrors(false, FocusingError);
            rt.setQuadAerrors(false, AlignmentError);
        end % n (numSeeds)
        % Calc the mean of rms values over all the seeds, and also the std of
        % the rms values.
        mean_rmsBPMvals_C(i,k)  =  mean(rmsBPMvals_C(i,:));
        err_rmsBPMvals_C(i,:,k)   =   prctile(rmsBPMvals_C(i,:), [10 90 20 80 30 70]);
        
        mean_rmsBPMvals_ML(i,k) = mean(rmsBPMvals_ML(i,:));
        err_rmsBPMvals_ML(i,:,k)  =  prctile(rmsBPMvals_ML(i,:), [10 90 20 80 30 70]);
    end % i (bot) 
end % k (top)

%% Test No ML
% Use this to generate data without the ML RM

% 
% % Setup arrrays
% rmsBPMvals_C        = zeros(numel(AlignmentError), numSeeds);
% mean_rmsBPMvals_C_new   = zeros(numel(AlignmentError), numel(FocusingError));
% err_rmsBPMvals_C_new    = zeros(numel(AlignmentError), 6, numel(FocusingError));
% 
% RespMatC = rt.calcRespMatC(dcorrStrength);
% 
% for k = 1:numel(FocusingError) % top
%     % Calc resp mat with conventional method in 'perfect' model
%     for i = 1:numel(AlignmentError) % bot
%         for n = 1:numSeeds
%             % set seed and turn on a set of errors
%             rng('shuffle')
%             rt.setQuadFerrors(true, FocusingError(k));
%             rt.setQuadAerrors(true, AlignmentError(i));
%             
%             bpmValsY = rt.track_getBPMreadings() + randn(numel(beamline.bpmlist),1)*BPMnoise;
% 
%             bpmValsCorrected_C  = rt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), BPMnoise);
%             rmsBPMvals_C(i,n) = rms(bpmValsCorrected_C);
% 
%             % Turn off errors
%             rt.setQuadFerrors(false, FocusingError(k));
%             rt.setQuadAerrors(false, AlignmentError(i));
%         end % n (numSeeds)
%         % Calc the mean of rms values over all the seeds, and also the std of
%         % the rms values.
%         mean_rmsBPMvals_C_new(i,k)  =  mean(rmsBPMvals_C(i,:));
%         err_rmsBPMvals_C_new(i,:,k)   =   prctile(rmsBPMvals_C(i,:), [10 90 20 80 30 70]);
%     end % i (bot) 
% end % k (top)
% 
% % sprintf('done')

%% Construct plot

LtopVal = length(dcorrStrength);
LbotVal = length(BPMnoise);
yConv = NaN(LtopVal,LbotVal*LtopVal+1);
yML = NaN(LtopVal,LbotVal*LtopVal+1);

for i=1:LtopVal
    yConv(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = mean_rmsBPMvals_C(:,i)';
    yML(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = mean_rmsBPMvals_ML(:,i)';
end
yConverrN = NaN(LtopVal,LbotVal*LtopVal+1);
yConverrP = NaN(LtopVal,LbotVal*LtopVal+1);
yMLerrN = NaN(LtopVal,LbotVal*LtopVal+1);
yMLerrP = NaN(LtopVal,LbotVal*LtopVal+1);

for i=1:LtopVal
    yConverrN(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = err_rmsBPMvals_C(:,3,i)';
    yConverrP(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = err_rmsBPMvals_C(:,4,i)';
    yMLerrN(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = err_rmsBPMvals_ML(:,3,i)';
    yMLerrP(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = err_rmsBPMvals_ML(:,4,i)';
end

labels_bot = strings(LtopVal*LbotVal,1);
labels_top = strings(LtopVal*LbotVal,1);
g=1;
for i=1:LtopVal
    for k=1:LbotVal
        labels_bot(g) = sprintf('%1g', BPMnoise(k)*1e6);
        labels_top(g) = sprintf('%1g', dcorrStrength(i)*1e6);
        g = g+1;
    end
end
labels_bot(end) = [];

x = 0:LbotVal*LtopVal;

xT = 0:LbotVal*LtopVal-1;
xTL = 0:1:LbotVal*LtopVal;

fill1   = [8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
fill2   = [NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN];
fill3   = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4];

% fill1   = [8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
% fill2   = [NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN];
% fill3   = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4];

figure1 = figure('PaperType','<custom>','PaperSize',[13 6.8],...
    'Color',[1 1 1])

ax2 = axes('Parent',figure1,'Position',[0.116531165311653 0.169724770642202 0.833468834688347 0.674311926605505]);
    %plot(ax2,x,y_ML_1c*1e3, '-.xb', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(1)))
    hold on
    area(ax2,xTL, fill1*1e3,'FaceAlpha',0.3, 'EdgeColor', 'none', 'FaceColor',[0 0.447058826684952 0.74117648601532])
    area(ax2,xTL, fill2*1e3,'FaceAlpha',0.3, 'EdgeColor', 'none', 'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625])
    area(ax2,xTL, fill3*1e3,'FaceAlpha',0.3, 'EdgeColor', 'none', 'FaceColor',[0.929411768913269 0.694117665290833 0.125490203499794])
    xlabel('a [\mu{}T]','LineWidth',1.25,'FontSize',12)
    xlim([0 15])
    
ax1 = axes('Parent',figure1,'Position',[0.116531165311653 0.165137614678899 0.833468834688347 0.678899082568809]);
    hold on
    eC(1) = plot(ax1,x,yConv(1,:)*1e3, '--+r');
    eML(1) = plot(ax1,x,yML(1,:)*1e3,':*k');
    for i=2:LtopVal
        eC(i) = plot(ax1,x,yConv(i,:)*1e3, '--+r');
    end

    for i=2:LtopVal
        eML(i) = plot(ax1,x,yML(i,:)*1e3, ':*k');
    end
    for i=1:LtopVal
        shade(ax1, (i*LbotVal)-LbotVal:i*LbotVal-1 ,err_rmsBPMvals_ML(:,1,i)*1e3,(i*LbotVal)-LbotVal:i*LbotVal-1,err_rmsBPMvals_ML(:,2,i)*1e3,'FillType',[2 1],'LineStyle', 'none')
        shade(ax1,(i*LbotVal)-LbotVal:i*LbotVal-1,err_rmsBPMvals_C(:,1,i)*1e3,(i*LbotVal)-LbotVal:i*LbotVal-1,err_rmsBPMvals_C(:,2,i)*1e3,'FillType',[2 1],'LineStyle', 'none')
    end
    xlabel('BPM noise [\mum]','LineWidth',1.25,'FontSize',12)
    ylabel('RMS y [mm]','LineWidth',1.25,'FontSize',12)
    xlim([0 15])
    
box(ax1,'on');
ax2.Box = 'off';

legend('Conventional','ML (NumObs = 100)')

set(ax2, 'FontSize',12,'LineWidth',1.25,'YColor','none', ...
    'XTick',xTL, 'XTickLabel', labels_top, 'XTickLabelRotation',45, ...
    'XAxisLocation', 'top','YAxisLocation','right')

set(ax1, 'FontSize',12,'LineWidth',1.25,'Color','none', ...
    'XTick',x, 'XTickLabel', labels_bot)

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 5])
print('-dpng','fig.png','-r600')

%output = [(x_b')*1e6 x_t(i,:)' y_C*1e3 y_ML_1*1e3 y_ML_2*1e3 y_Cerr*1e3 y_ML_1err*1e3 y_ML_2err*1e3];


