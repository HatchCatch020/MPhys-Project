% This script is a modification of the original script DefineBeamline
% exmaple in SAMM. It has been modified by CM for the application of
% various beam imperfections and the analysis of orbit corrections.

beam        = Beam(Positron);
beam.energy = 2.0 * PhysicalUnits.GeV;

ncells      = 16;

corrlength  = 0.2;                       % corrector length in metres

quaderrorLength = 1e-6;

drift1 = Drift;
  drift1.length = 5.00 - corrlength-quaderrorLength;     % metres 4.8
  
drift2 = Drift;
  drift2.length = 6.00;                  % metres
  
drift3 = Drift;
  drift3.length = 0.10;                  % metres
  
drift4 = Drift;
  drift4.length = 4.75;                  % metres
  
drift5 = Drift;
  drift5.length = 5.75 - corrlength-quaderrorLength;      % metres
  
dipole1 = Dipole;
  dipole1.length = 2.50;                 % metres
  dipole1.angle  = 2*pi/ncells;          % radians
  dipole1.field  = beam.rigidity * dipole1.angle / dipole1.length; % tesla
  
sextF  = Sextupole;
  sextF.length   = 0.15;                 % metres
  sextF.gradient = 0.27 * beam.rigidity; % tesla/metre^2

sextD  = Sextupole;
  sextD.length   = 0.15;                 % metres
  sextD.gradient =-0.37 * beam.rigidity; % tesla/metre^2

rfcav1 = RFCavity;
  rfcav1.length   = 0.2965;              % metres
  rfcav1.voltage  = 600e3*sign(beam.species.charge); % volts
  % Note that we define a separate variable for the nominal rf frequency.
  % The actual frequency of the cavity is set below, to an harmonic of
  % the revolution frequency.
  rffreq          = 500 * PhysicalUnits.megahertz;

drift6 = Drift;
   drift6.length = drift5.length - rfcav1.length; % metres
  
bpm    = cell(1,2*ncells);
corr   = cell(1,  ncells);
quaderrlist  = cell(1,2*ncells);

quadlist = cell(1, 2*ncells);
quaddlist = cell(1, ncells);

bl = Beamline;
  
for n = 1:ncells
    
    quadlist{2*n-1} = Quadrupole; % d
       quadlist{2*n-1}.length   = 0.25;                 % metres
       quadlist{2*n-1}.gradient = 0.45 * beam.rigidity; % tesla/metre
  
    quadlist{2*n} = Quadrupole;   % f
       quadlist{2*n}.length   = 0.25;                 % metres
       quadlist{2*n}.gradient =-0.40 * beam.rigidity; % tesla/metre
       
    quaderrlist{2*n-1} = OrbitCorrector;
       quaderrlist{2*n-1}.length = quaderrorLength;
       quaderrlist{2*n-1}.field  = [0, 0];
    quaderrlist{2*n} = OrbitCorrector;
       quaderrlist{2*n}.length = quaderrorLength;
       quaderrlist{2*n}.field  = [0, 0];    
       
    corr{n} = OrbitCorrector;
       corr{n}.length = corrlength;
       corr{n}.field  = [0, 0];
  
    bl.AppendComponent(drift2);
    bpm{2*n-1}     = BeamPositionMonitor;
    bl.AppendComponent(bpm{2*n-1});
    bl.AppendComponent(quadlist{2*n-1});
    bl.AppendComponent(quaderrlist{2*n-1});
    bl.AppendComponent(drift3);
    bl.AppendComponent(sextF);
    bl.AppendComponent(drift4);
    bl.AppendComponent(dipole1);
    bl.AppendComponent(corr{n});
    bl.AppendComponent(drift1);
    bpm{2*n}       = BeamPositionMonitor;
    bl.AppendComponent(bpm{2*n});
    bl.AppendComponent(quadlist{2*n});
    bl.AppendComponent(quaderrlist{2*n});
    bl.AppendComponent(drift3);
    bl.AppendComponent(sextD);
    bl.AppendComponent(drift5);      
end

bl.componentlist{end} = drift6;
bl.AppendComponent(rfcav1);

svals = bl.ComputePositions();

% Set Master Oscillator frequency before setting RF cavity harmonic and phase
MasterOscillator.SetFrequency(beam.beta*PhysicalConstants.SpeedOfLight/svals(end));
rfcav1.harmonic = floor(rffreq/MasterOscillator.GetFrequency());
rfcav1.phase    = pi - 2*pi*svals(end-1)*rfcav1.frequency/beam.beta/PhysicalConstants.SpeedOfLight;

% Check the stability of the longitudinal motion
m = ComputeTransferMatrix(bl,[1 numel(bl.componentlist)],beam,zeros(6,1));
if max(abs(eig(m(:,:,end))))-1>1e-6 % if true, then we need to reverse the polarity of the rf
    rfcav1.voltage = -1*rfcav1.voltage;
    m = ComputeTransferMatrix(bl,[1 numel(bl.componentlist)],beam,zeros(6,1));
end

beamline.quadlist       = quadlist;
beamline.quaderrlist    = quaderrlist;
beamline.corrlist       = corr;
beamline.bpmlist        = bpm;
beamline.beam           = beam;
beamline.bl             = bl;

% Define variables containing original quad strengths
for n = 1:numel(beamline.quadlist)
    beamline.gradient0(n) = beamline.quadlist{n}.gradient;
end

% m(:,:,end)