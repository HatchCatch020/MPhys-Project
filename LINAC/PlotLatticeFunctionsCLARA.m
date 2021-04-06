DefineCLARABeamline

for i = 1:(length(bl.componentlist)-length(bpmlist))
    if bl.componentlist{i}.name == "BPM"
        svals(i) = [];
        pvals(i) = [];
        bl.componentlist(i) = [];
    end
end

nmax = length(bl.componentlist)+1;

betax  =  7.007755;
alphax = -0.612920;
gammax = (1 + alphax^2)/betax;

betay  =  6.604553;
alphay = -0.549430;
gammay = (1 + alphay^2)/betay;

sigz   = 1e-3;
sigdp  = 1e-3;

epsz   = sigz*sigdp;
betaz  = sigz^2/epsz;
alphaz = 0;
gammaz = (1 + alphaz^2)/betaz;

beta = zeros(6,6,3,nmax);

beta(1,1,1,1) = betax;
beta(1,2,1,1) =-alphax;
beta(2,1,1,1) =-alphax;
beta(2,2,1,1) = gammax;

beta(3,3,2,1) = betay;
beta(3,4,2,1) =-alphay;
beta(4,3,2,1) =-alphay;
beta(4,4,2,1) = gammay;

beta(5,5,3,1) = betaz;
beta(5,6,3,1) =-alphaz;
beta(6,5,3,1) =-alphaz;
beta(6,6,3,1) = gammaz;

m = ComputeTransferMatrix(bl,[1 nmax-1],beam,zeros(6,1));

for n = 2:nmax
    m1 = m(:,:,n);
    beta(:,:,1,n) = m1*beta(:,:,1,1)*m1';
    beta(:,:,2,n) = m1*beta(:,:,2,1)*m1';
    beta(:,:,3,n) = m1*beta(:,:,3,1)*m1';        
end

beta1 = permute(beta,[4 1 2 3]);

% close all
figure(1)
betavals = [svals; ...
    (pvals.*real(beta1(:,1,1,1)))'/pvals(1); ...
    (pvals.*real(beta1(:,3,3,2)))'/pvals(1)]';

subplot(3,1,1)
hold off
plot(svals, betavals(:,2), '-k');
hold on
plot(svals, betavals(:,3), '-r');
axis([0 svals(end) 0 40])
% xlabel('s [m]');
ylabel('\beta [m]');
legend('\beta_x', '\beta_y', 'Location', 'northwest');
title("SAMM Lattice Functions")

subplot(3,1,2)
plot(svals, pvals*PhysicalConstants.SpeedOfLight/PhysicalUnits.MeV, '-k');
ylabel('P_0 (MeV/c)');

subplot(3,1,3)
plot(svals, real(beta1(:,1,6,3)./beta1(:,6,6,3)), '-k');
axis([0 svals(end) -inf inf])
xlabel('s [m]');
ylabel('\eta [m]');

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 5])
print('-dpng','CLARA Lattice Functions.png','-r600')

figure(2)
hold off
plot(svals, '-r')
hold on
plot(NumericInfo(:,2), '-.b')
xlabel('Component list');
ylabel('s [m]');
legend('SAMM s positions', 'MADX s positions', 'Location', 'northwest');

figure(3)
subplot(3,1,1)
hold off
plot(NumericInfo(:,2), NumericInfo(:,4), '-k');
hold on
plot(NumericInfo(:,2), NumericInfo(:,7), '-r');
axis([0 svals(end) 0 40])
% xlabel('s [m]');
ylabel('\beta [m]');
legend('\beta_x', '\beta_y', 'Location', 'northwest');
title("MAD-X Lattice Functions")

% subplot(3,1,2)
% hold off
% plot(NumericInfo(:,2), NumericInfo(:,4)-betavals(:,2), '-k')
% hold on
% plot(NumericInfo(:,2), NumericInfo(:,7)-betavals(:,3), '-r')
% xlabel('s [m]');
% ylabel('\beta_{residual} [m]');
% legend('\beta_x', '\beta_y', 'Location', 'northwest');
% title("Difference of lattice functions in SAMM and MAD-X")

subplot(3,1,3)
plot(NumericInfo(:,2), NumericInfo(:,6), '-k');
%axis([0 svals(end) -inf inf])
xlabel('s [m]');
ylabel('\eta [m]');

set(gcf,'PaperUnits','inches')
set(gcf,'PaperPosition',[1 1 9 5])
print('-dpng','CLARA Lattice Functions MADX.png','-r600')

