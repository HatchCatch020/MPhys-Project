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
numSeeds   = 10;

% Varying
ML_numobs       = linspace(100,100,1);
BPMnoise        = 70e-6;%linspace(40,90,6)*1e-6;

dcorrStrength   = [1e-4 1e-5 1e-6];
AlignmentError    = linspace(10,50,5)*1e-6;   % m

FocusingError   = 0.01;%linspace(1,9,5)*1e-1;

%% Test

% Setup arrrays
rmsBPMvals_C        = zeros(numel(AlignmentError), 1, numSeeds);
rmsBPMvals_ML       = zeros(numel(AlignmentError), numel(ML_numobs), numSeeds);
mean_rmsBPMvals_C   = zeros(numel(AlignmentError), 1, numel(dcorrStrength));
mean_rmsBPMvals_ML  = zeros(numel(AlignmentError), numel(ML_numobs), numel(dcorrStrength));
errp_rmsBPMvals_C   = zeros(numel(AlignmentError), 1, numel(dcorrStrength));
errn_rmsBPMvals_C   = zeros(numel(AlignmentError), 1, numel(dcorrStrength));
errp_rmsBPMvals_ML  = zeros(numel(AlignmentError), numel(ML_numobs), numel(dcorrStrength));
errn_rmsBPMvals_ML  = zeros(numel(AlignmentError), numel(ML_numobs), numel(dcorrStrength));

for k = 1:numel(dcorrStrength)
    % Calc resp mat with conventional method in 'perfect' model
    RespMatC = lt.calcRespMatC(dcorrStrength(k));
    for i = 1:numel(AlignmentError)
        for n = 1:numSeeds
            % set seed and turn on a set of errors
            rng(n)
            lt.setQuadFerrors(true, FocusingError);
            lt.setQuadAerrors(true, AlignmentError(i));
            
            bpmValsY = lt.track_getBPMreadings();

            bpmValsCorrected_C  = lt.getBPMvalues_corr(pinv(RespMatC,  1e-3)*(-bpmValsY), BPMnoise);
            rmsBPMvals_C(i,1,n) = rms(bpmValsCorrected_C);

            for j = 1:numel(ML_numobs)
                RespMatML = lt.calcRespMatML(ML_numobs(j), BPMnoise, dcorrStrength(k));
                bpmValsCorrected_ML  = lt.getBPMvalues_corr(pinv(RespMatML,  1e-3)*(-bpmValsY), BPMnoise);
                rmsBPMvals_ML(i,j,n) = rms(bpmValsCorrected_ML);
            end
            % Turn off errors
            lt.setQuadFerrors(false, FocusingError);
            lt.setQuadAerrors(false, AlignmentError(i));
        end % n (numSeeds)
        % Calc the mean of rms values over all the seeds, and also the std of
        % the rms values.
        mean_rmsBPMvals_C(i,k)  =  mean(rmsBPMvals_C(i,1,:));
        errn_rmsBPMvals_C(i,k)   =   prctile(rmsBPMvals_C(i,1,:), [5]);
        errp_rmsBPMvals_C(i,k)   =   prctile(rmsBPMvals_C(i,1,:), [95]);
        for j=1:numel(ML_numobs) 
            mean_rmsBPMvals_ML(i,j,k) = mean(rmsBPMvals_ML(i,j,:));
            errn_rmsBPMvals_ML(i,j,k)  =  prctile(rmsBPMvals_ML(i,j,:), [5]);
            errp_rmsBPMvals_ML(i,j,k)  =  prctile(rmsBPMvals_ML(i,j,:), [95]);
        end % j 
    end % i (BPMnoise) 
end % k (kickStrength)

%% Construct plot

x_t = [];
x_b = [AlignmentError AlignmentError AlignmentError];
for i=1:numel(dcorrStrength), x_t = [x_t ones(numel(AlignmentError))*dcorrStrength(i)]; end

y_ML_1 = [];
y_ML_2 = [];
y_ML_1err = [];
%y_ML_2err = [];
y_C     = mean_rmsBPMvals_C(:);
y_Cerr  =  errp_rmsBPMvals_C(:);
for k=1:numel(dcorrStrength)
    y_ML_1 = [y_ML_1 mean_rmsBPMvals_ML(:,1,k)]; 
    %y_ML_2 = [y_ML_2 mean_rmsBPMvals_ML(:,2,k)];
    y_ML_1err  =  [y_ML_1err errp_rmsBPMvals_ML(:,1,k)];
    %y_ML_2err  =  [y_ML_2err std_rmsBPMvals_ML(:,2,k)];
end
y_ML_1 = y_ML_1(:);
%y_ML_2 = y_ML_2(:);
y_ML_1err = y_ML_1err(:);
%y_ML_2err = y_ML_2err(:);

output = [(x_b')*1e6 x_t(i,:)' y_C*1e3 y_ML_1*1e3  y_Cerr*1e3 y_ML_1err*1e3];

LtopVal = length(dcorrStrength);
LbotVal = length(AlignmentError);
yConv = NaN(LtopVal,LbotVal*LtopVal);
yML = NaN(LtopVal,LbotVal*LtopVal,numel(ML_numobs));
for i=1:LtopVal
    yConv(i,(i*LbotVal)-LbotVal+1:i*LbotVal) = mean_rmsBPMvals_C(:,i)';
    for j = 1:numel(ML_numobs)
        yML(i,(i*LbotVal)-LbotVal+1:i*LbotVal,j) = mean_rmsBPMvals_ML(:,j,i)';
    end
end
yConverr = NaN(LtopVal,2,LbotVal*LtopVal);
yMLerr = NaN(LtopVal,2,LbotVal*LtopVal,numel(ML_numobs));
for i=1:LtopVal
    yConverr(i,1,(i*LbotVal)-LbotVal+1:i*LbotVal) = errn_rmsBPMvals_C(:,i)';
    yConverr(i,2,(i*LbotVal)-LbotVal+1:i*LbotVal) = errp_rmsBPMvals_C(:,i)';
    for j = 1:numel(ML_numobs)
        yMLerr(i,1,(i*LbotVal)-LbotVal+1:i*LbotVal,j) = errn_rmsBPMvals_ML(:,j,i)';
        yMLerr(i,2,(i*LbotVal)-LbotVal+1:i*LbotVal,j) = errp_rmsBPMvals_ML(:,j,i)';
    end
end

labels_bot = strings(LtopVal*LbotVal,1);
labels_top = strings(LtopVal*LbotVal,1);
g=1;
for i=1:LtopVal
    for k=1:LbotVal
        labels_bot(g) = sprintf('%1g', AlignmentError(k)*1e6);
        labels_top(g) = sprintf('%1g', dcorrStrength(i)*1e6);
        g = g+1;
    end
end

x = 0:LbotVal*LtopVal-1;

xT = 0:LbotVal*LtopVal-1;
xTL = 0:2:LbotVal*LtopVal;

fill1   = [8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN NaN NaN NaN NaN NaN ];
fill2   = [NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4 8e-4 NaN NaN NaN NaN ];
fill3   = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 8e-4 8e-4 8e-4 8e-4 8e-4];

figure

ax2 = axes('Position',[0.1 0.15 0.85 0.72]);
    %plot(ax2,x,y_ML_1c*1e3, '-.xb', 'DisplayName', sprintf('ML (NumObs = %d)', ML_numobs(1)))
    hold on
    area(ax2,xT, fill1*1e3,'FaceAlpha',0.5, 'EdgeColor', 'none')
    area(ax2,xT, fill2*1e3,'FaceAlpha',0.5, 'EdgeColor', 'none')
    area(ax2,xT, fill3*1e3,'FaceAlpha',0.5, 'EdgeColor', 'none')
    xlabel('\Deltacorr strength [\muT]','LineWidth',3,'FontWeight','bold','FontSize',18)
    
ax1 = axes('Position',[0.1 0.15 0.85 0.72]);
    hold on
    for i=1:LtopVal
        eC(i) = errorbar(ax1,x,yConv(i,:)*1e3, reshape(yConverr(i,1,:), [1 numel(yConverr(i,1,:))])*1e3,...
            reshape(yConverr(i,2,:), [1 numel(yConverr(i,2,:))])*1e3, '--+r'); 
    end
    set(eC(1), 'DisplayName', 'Conventional')

    for i=1:LtopVal
        eML(i) = errorbar(ax1,x,yML(i,:)*1e3,reshape(yMLerr(i,1,:),[1 numel(yMLerr(i,1,:))])*1e3,...
            reshape(yMLerr(i,2,:),[1 numel(yMLerr(i,2,:))])*1e3, ':*k')
    end
    set(eML(1), 'DisplayName', 'ML (NumObs = 100)')
    xlabel('Kick Strength [\mum]','LineWidth',3,'FontWeight','bold','FontSize',18)
    ylabel('\Deltay_{RMS} [mm]','LineWidth',3,'FontWeight','bold','FontSize',18)
    
box(ax1,'on');
ax2.Box = 'off';

legend

set(ax2, 'FontSize',18,'FontWeight','bold','LineWidth',3,'YColor','none', ...
    'XTick',xT, 'XTickLabel', labels_top, 'XTickLabelRotation',45, ...
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
