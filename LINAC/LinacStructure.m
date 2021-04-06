classdef LinacStructure < handle
    % LinacStructure class
    % 
    % LinacStructure:
    %   name
    %   structuretype
    %   length
    %   ncell
    %   voltage
    %   harmonic
    %   phase
    %   globalclock
    %   aperture
    %
    % Methods:
    %   Track
    
    properties
        name          = ''; % string
        structuretype = 'TravellingWave'; % 'TravellingWave' or 'StandingWave'
        length        = 0;  % in metres
        ncell         = 1;  % number of cells
        voltage       = 0;  % in volts
        harmonic      = 1;  % frequency (in units of master oscillator frequency)
        phase         = pi; % relative to master oscillator, in radians
        globalclock   = 1;  % 1 = synchronise with global clock
        aperture      = []; % 1x2 array of elliptical aperture half-axes, in metres
    end % properties

    properties (Dependent=true)
        frequency;
    end % properties (dependent)

    methods
        
        function set.frequency(linacstructure,f)
            f1 = linacstructure.harmonic * MasterOscillator.GetFrequency();
            MasterOscillator.SetFrequency(MasterOscillator.GetFrequency()*f/f1);
        end
        
        function f = get.frequency(linacstructure)
            f = linacstructure.harmonic * MasterOscillator.GetFrequency();
        end
        
        function beam = Track(linacstructure,beam)
            % beam2 = RFAcceleratingStructure.Track(beam1)
            % Applies the dynamical map for a linac structure to the particles
            % in beam1.  The reference energy is changed by qV*cos(phase).
                        
            [x0, px0, y0, py0, ct0, dp0] = beam.GetParticles();
            
            mofreq = MasterOscillator.GetFrequency();
            
            beam.globaltime = beam.globaltime - floor(beam.globaltime*mofreq)/mofreq;            
            gt  = beam.globaltime * linacstructure.globalclock;
            
            nc     = linacstructure.ncell;
            f      = linacstructure.harmonic*mofreq;
            dL     = linacstructure.length / nc;
            phi    = linacstructure.phase;
            cosphi = cos(phi);

            dE     = beam.species.charge * linacstructure.voltage * cosphi / nc;
            eta    = 1 - cos(2*phi);
            c1     = 2*sqrt(eta/8);
            
            for n = 1:nc
                
                % First, apply a drift map through L/2
                % to the longitudinal coordinate
                beta0  = beam.beta;
                d1     = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
                ct0    = ct0 + dL*(1 - (1 + beta0*dp0)./d1)/beta0/2;

                % Calculate the transverse part of the map
                E0    = beam.energy;
                alpha = c1*log(1 + dE/E0)/cosphi/2;
                c2    = cosphi*sin(alpha)/c1;

                if abs(dE/E0)>1e-9
                    r11 =  cos(alpha) - c2;
                    r12 =  2*c2*E0*dL/dE;
                    r21 = -(c2 + sin(alpha)^2/c2)*dE/dL/E0/2;
                    r22 =  (cos(alpha) + c2)*E0/(E0+dE);
                else
                    r11 = 1;
                    r12 = dL;
                    r21 = 0;
                    r22 = 1;
                end
                
                % Now apply the map to the transverse variables
                x1  = r11*x0 + r12*px0;
                px0 = r21*x0 + r22*px0;
                x0  = x1;

                y1  = r11*y0 + r12*py0;
                py0 = r21*y0 + r22*py0;
                y0  = y1;

                % Calculate the initial energy deviation of each particle
                P0  = beam.momentum;
                Edc = (dp0 + 1/beta0)*P0;

                % Increase reference energy on passing through each cell
                beam.energy = beam.energy + dE;
                
                % Calculate the final momentum deviation of each particle
                beta0 = beam.beta;
                P1    = beam.momentum;
                
                t = gt - ct0/PhysicalConstants.SpeedOfLight;
                Edc = Edc + ...
                    beam.species.charge * linacstructure.voltage * cos(2*pi*f*t + phi) / nc / PhysicalConstants.SpeedOfLight;

                dp0 = Edc/P1 - 1/beta0;

                % Finally, apply a drift map through L/2
                % to the longitudinal coordinate
                d1     = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
                ct0    = ct0 + dL*(1 - (1 + beta0*dp0)./d1)/beta0/2;
            
            end
            
            % Set the new values for the dynamical variables
            beam.SetParticles(x0, px0, y0, py0, ct0, dp0);

        end % function Track
        
        function TrackLibrary(rfstructure,trackingMethod,libroutine)
            
            rfstructure1.length   = rfstructure.length;
            rfstructure1.voltage  = rfstructure.voltage;
            rfstructure1.harmonic = rfstructure.harmonic;
            rfstructure1.phase    = rfstructure.phase;
            
            rfstructure1.masteroscillatorfrequency = MasterOscillator.GetFrequency();

            if(~isempty(rfstructure.aperture))
               rfstructure1.apertureX = rfstructure.aperture(1);
               rfstructure1.apertureY = rfstructure.aperture(2);
            end
            
            calllib(trackingMethod,[libroutine 'RFAcceleratingStructure'],rfstructure1);
            
        end % function TrackLibrary

        
    end % methods
    
end % classdef RFAcceleratingStructure
