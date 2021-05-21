% C. Monaghan
% Define the CLARA beamline based on the input of a csv file containing the
% information for each component.       

clear all

%% Open File
filename = 'CLARA_lattice.csv';
fileID = fopen(filename,'r','n','UTF-8');

CLARA_lattice = textscan(fileID, '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]', 'Delimiter', ',', 'TextType', 'string',  'ReturnOnError', false);

fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(CLARA_lattice{1}),length(CLARA_lattice)-1);
for col=1:length(CLARA_lattice)-1
    raw(1:length(CLARA_lattice{col}),col) = mat2cell(CLARA_lattice{col}, ones(length(CLARA_lattice{col}), 1));
end
numericData = NaN(size(CLARA_lattice{1},1),size(CLARA_lattice,2));

for col=[1,4,5,6,7,8,9,10,11,12,13,14,15,16]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = CLARA_lattice{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [1,4,5,6,7,8,9,10,11,12,13,14,15,16]);
rawStringColumns = string(raw(:, [2,3]));
NumericInfo = cell2mat(rawNumericColumns);

%% DefineBeamline
beam  = Beam(Electron);
momentum0 = 35.0 * PhysicalUnits.MeV / PhysicalConstants.SpeedOfLight;

beam.momentum = momentum0;

driftlist   = cell(0,0);
quadflist   = cell(0,0);
quadlist    = cell(0,0);
quaddlist   = cell(0,0);
quaderrlist = cell(0,0);
corrlist    = cell(0,0);
bpmlist     = cell(0,0);
Lcavitylist = cell(0,0);
Sbendlist   = cell(0,0);

pvals    = zeros(length(NumericInfo)+1,1);
pvals(1) = beam.momentum;

bl = Beamline;
MasterOscillator.SetFrequency(2998500000);

quaderrorLength   = 1e-6;   % metre

c=0;

% For each component in the file execute the corresponding block of code
for n=1:length(NumericInfo)
    switch rawStringColumns(n,2)
        case "Drift"
            driftlist{end+1} = Drift;
            driftlist{end}.length = NumericInfo(n,3);
            bl.AppendComponent(driftlist{end});
        case "Quadrupole"
            quadlist{end+1}         = Quadrupole;
                quadlist{end}.length    = NumericInfo(n,3);
                quadlist{end}.gradient  = NumericInfo(n,11) * beam.rigidity;
                bl.AppendComponent(quadlist{end});
            % Comment out here for PlotLatticeFunctionsCLARA
            quaderrlist{end+1} = OrbitCorrector;
                quaderrlist{end}.length = quaderrorLength;
                quaderrlist{end}.field = [0, 0];
                bl.AppendComponent(quaderrlist{end});
        case "Kicker"
            corrlist{end+1} = OrbitCorrector;
            corrlist{end}.length = NumericInfo(n,3);
            corrlist{end}.field = [0, 0];
            bl.AppendComponent(corrlist{end});
        case "Monitor"
            bpmlist{end+1} = BeamPositionMonitor;
            bpmlist{end}.name = "BPM";
            bl.AppendComponent(bpmlist{end});
            pvals(numel(bl.componentlist)+1) = beam.momentum;
            % Uncomment here for PlotLatticeFunctionsCLARA
%             pvals(c+1) = beam.momentum;
%             c=c+1;
            % Drift added by AW 12/2/2021
            driftlist{end+1} = Drift;
            driftlist{end}.length = NumericInfo(n,3);
            bl.AppendComponent(driftlist{end});
        case "Lcavity"
            % case "Lcavity" modified by AW 12/2/2021
            Lcavitylist{end+1} = LinacStructure;
            Lcavitylist{end}.structuretype = 'TravellingWave';
            Lcavitylist{end}.length = NumericInfo(n,3);
            Lcavitylist{end}.ncell = round(3*NumericInfo(n,3)*NumericInfo(n,14)/PhysicalConstants.SpeedOfLight);
            Lcavitylist{end}.harmonic = round(NumericInfo(n,14)/MasterOscillator.GetFrequency());
            Lcavitylist{end}.voltage = -NumericInfo(n,12) * NumericInfo(n,3);
            Lcavitylist{end}.phase = NumericInfo(n,13)*2*pi/360;
            Lcavitylist{end}.globalclock = false;
            bl.AppendComponent(Lcavitylist{end});
            beam.energy = beam.energy + beam.species.charge * Lcavitylist{end}.voltage * cos(Lcavitylist{end}.phase);
        case "Sbend"
            Sbendlist{end+1} = Dipole;
            Sbendlist{end}.length = NumericInfo(n,3);
            Sbendlist{end}.curvature = 1/NumericInfo(n,10);
            theta = NumericInfo(n,3)/NumericInfo(n,10);
            Sbendlist{end}.e1 = theta/2;
            Sbendlist{end}.e2 = theta/2;
            Sbendlist{end}.field = beam.rigidity * Sbendlist{end}.curvature;
            bl.AppendComponent(Sbendlist{end});
        case {"Ecollimator","Rcollimator","Instrument"}
            driftlist{end+1} = Drift;
            driftlist{end}.length = NumericInfo(n,3);
            bl.AppendComponent(driftlist{end});
        otherwise
            % Uncomment here for PlotLatticeFunctionsCLARA
%             driftlist{end+1} = Drift;
%             driftlist{end}.length = NumericInfo(n,3);
%             bl.AppendComponent(driftlist{end});
%             pvals(c+1) = beam.momentum;
%             c=c+1;
            continue;
    end
    pvals(numel(bl.componentlist)+1) = beam.momentum;
    % Uncomment here for PlotLatticeFunctionsCLARA
%     pvals(c+1) = beam.momentum;
%     c=c+1;
end

fprintf('Initial beam momentum = %6.5g MeV/c\n', momentum0 / PhysicalUnits.MeV * PhysicalConstants.SpeedOfLight);
fprintf('  Final beam momentum = %6.5g MeV/c\n', beam.momentum / PhysicalUnits.MeV * PhysicalConstants.SpeedOfLight);

beam.momentum = momentum0;

svals = bl.ComputePositions();

fprintf('Total beamline length = %6.5g m\n', svals(end) );

% Setup structure for all useful variables in this script
beamline.driftlist      = driftlist;
beamline.quadlist       = quadlist;
beamline.quaderrlist    = quaderrlist;
beamline.corrlist       = corrlist;
beamline.bpmlist        = bpmlist;
beamline.Lcavitylist    = Lcavitylist;
beamline.Sbendlist      = Sbendlist;
beamline.beam           = beam;
beamline.bl             = bl;
beamline.momentum0 = momentum0;

% Define variables containing original quad strengths
for n = 1:numel(beamline.quadlist)
    beamline.gradient0(n) = beamline.quadlist{n}.gradient;
end

